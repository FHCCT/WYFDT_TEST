#include "dt.h"

namespace {
using Vec = Eigen::Vector3d;
using Vec2 = Eigen::Vector2d;
using Polygon = std::vector<Vec2, Eigen::aligned_allocator<Vec2>>;
using Face = std::array<int, 4>;
using Edge = std::pair<int, int>;
using Tree = dt::BinaryTree<Eigen::Matrix3d>;

Edge edgeKey(int a, int b) { return std::minmax(a, b); }
double cross2(const Vec2& a, const Vec2& b) { return a.x()*b.y()-a.y()*b.x(); }

// The index, feature priorities, anchors and pending transactions belong to one
// input-surface cleanup. No volume topology exists while this workspace runs.
class SurfaceImprint {
    dt::Mesh& mesh;
    Tree tree;
    std::vector<std::vector<int>> incident;
    std::map<Edge, std::vector<int>> edges;
    std::map<Edge, std::vector<int>> segments;
    std::set<Edge> imprintEdges;
    std::vector<std::pair<int,int>> priority;
    std::vector<unsigned char> fixed, anchored, queued;
    std::queue<int> pending;
    std::map<int,int> faceCounts;
    std::vector<int> sourceFace;
    // A complete overlay handles an original pair once. Its child diagonals
    // are not new sheets to imprint repeatedly against that same pair.
    std::set<Edge> completedOverlays;
    const size_t originalVertices;
    int indexed = 0;
    const double nearRatio = 0.02;
    const double parallelCos = std::cos(5.0*PI/180.0);
    const double turnCos = std::cos(15.0*PI/180.0);

    Vec point(int n) const { return Vec(mesh.V[n].data()); }
    bool live(int f) const { return f>=0 && f<static_cast<int>(mesh.F.size()) && mesh.F[f][0]>=0; }
    Vec normal(const Face& f) const { return (point(f[1])-point(f[0])).cross(point(f[2])-point(f[0])); }
    double length(const Face& f) const {
        return std::max({(point(f[0])-point(f[1])).norm(),(point(f[1])-point(f[2])).norm(),(point(f[2])-point(f[0])).norm()});
    }
    static Tree::BoundingBox bounds(const dt::Mesh& m) {
        Tree::BoundingBox b;
        for(int k=0;k<3;++k)b[k]=b[k+3]=m.V[0][k];
        for(const auto& p:m.V)for(int k=0;k<3;++k){b[k]=std::min(b[k],p[k]);b[k+3]=std::max(b[k+3],p[k]);}
        double scale=0,coordinate=0;
        for(int k=0;k<3;++k){scale=std::max(scale,b[k+3]-b[k]);coordinate=std::max(coordinate,std::max(std::abs(b[k]),std::abs(b[k+3])));}
        const double pad=1e-8*scale+64*std::numeric_limits<double>::epsilon()*coordinate;
        for(int k=0;k<3;++k){b[k]-=pad;b[k+3]+=pad;}
        return b;
    }
    Tree::BoundingBox box(const Face& f, double pad=0) const {
        Tree::BoundingBox b;
        for(int k=0;k<3;++k){b[k]=std::min({mesh.V[f[0]][k],mesh.V[f[1]][k],mesh.V[f[2]][k]})-pad;b[k+3]=std::max({mesh.V[f[0]][k],mesh.V[f[1]][k],mesh.V[f[2]][k]})+pad;}
        return b;
    }
    void attach(int t) {
        const Face& f=mesh.F[t]; Eigen::Matrix3d tri;
        for(int k=0;k<3;++k){incident[f[k]].push_back(t);edges[edgeKey(f[k],f[(k+1)%3])].push_back(t);tri.row(k)=point(f[k]).transpose();}
        if(t<indexed)tree.reinsert(tri,t);else{tree.insert(tri,t);++indexed;}
        ++faceCounts[f[3]];
    }
    void detach(int t) {
        const Face& f=mesh.F[t];tree.remove(t,t);--faceCounts[f[3]];
        for(int k=0;k<3;++k){
            auto& v=incident[f[k]];v.erase(std::remove(v.begin(),v.end(),t),v.end());
            auto it=edges.find(edgeKey(f[k],f[(k+1)%3]));auto& list=it->second;
            list.erase(std::remove(list.begin(),list.end(),t),list.end());if(list.empty())edges.erase(it);
        }
    }
    void enqueue(int t) {
        if(queued.size()<mesh.F.size())queued.resize(mesh.F.size(),0);
        if(live(t)&&!queued[t]){queued[t]=1;pending.push(t);}
    }
    std::array<int,3> key(const Face& f) const {
        std::array<int,3> k={{f[0],f[1],f[2]}};std::sort(k.begin(),k.end());return k;
    }

    // Positive projected area is required: ordinary adjacent triangles sharing
    // an edge or vertex must never be mistaken for a double surface.
    bool overlapping(int a,int b,double& tolerance) const {
        const Face& fa=mesh.F[a];const Face& fb=mesh.F[b];
        Vec na=normal(fa),nb=normal(fb);double la=na.norm(),lb=nb.norm();
        if(!(la>0&&lb>0)||std::abs(na.dot(nb))<parallelCos*la*lb)return false;
        na/=la;nb/=lb;
        tolerance=nearRatio*std::min(length(fa),length(fb));
        const Vec origin=point(fa[0]),u=(point(fa[1])-origin).normalized(),v=na.cross(u);
        std::array<Vec2,3> clip;
        Polygon polygon;
        for(int i=0;i<3;++i){Vec p=point(fa[i])-origin;clip[i]=Vec2(p.dot(u),p.dot(v));p=point(fb[i])-origin;polygon.push_back(Vec2(p.dot(u),p.dot(v)));}
        for(int e=0;e<3&&!polygon.empty();++e){
            Polygon next;const Vec2 dir=clip[(e+1)%3]-clip[e];
            for(size_t i=0;i<polygon.size();++i){const Vec2 p=polygon[i],q=polygon[(i+1)%polygon.size()];double dp=cross2(dir,p-clip[e]),dq=cross2(dir,q-clip[e]);
                if(dp>=0)next.push_back(p);
                if((dp<0)!=(dq<0))next.push_back(p+(q-p)*(dp/(dp-dq)));
            }
            polygon.swap(next);
        }
        if(polygon.size()<3)return false;
        double area=0;
        for(size_t i=0;i<polygon.size();++i){area+=cross2(polygon[i],polygon[(i+1)%polygon.size()]);
            const Vec p=origin+u*polygon[i].x()+v*polygon[i].y();
            if(std::abs((point(fb[0])-p).dot(nb)/na.dot(nb))>tolerance)return false;
        }
        return std::abs(area)>1e-8*std::min(la,lb);
    }

    bool intersects(const Face& a,const Face& b,int moved,const Vec& old,bool beforeA,bool beforeB) const {
        int ia[3]={a[0],a[1],a[2]},ib[3]={b[0],b[1],b[2]};
        Vec pa[3],pb[3];double* ap[3];double* bp[3];
        for(int k=0;k<3;++k){pa[k]=beforeA&&a[k]==moved?old:point(a[k]);pb[k]=beforeB&&b[k]==moved?old:point(b[k]);ap[k]=pa[k].data();bp[k]=pb[k].data();}
        return dt::GEOM_FUNC::tri_tri_intersect3d_fast(ia,ib,ap,bp)!=0;
    }

    // Either snap an existing vertex onto the master face/edge, or identify the
    // crossing of two nearby projected edges. Split ALL incident faces, also
    // those with another geometric ID, so no hanging edge or crack is created.
    bool transact(int master,int slave,int moved,const Vec& destination,Edge target,Edge source) {
        const Face& reference=mesh.F[master];
        Vec exactPoint=destination;
        // Nearly coplanar is not enough here: a rounding-sized gap in an
        // imprinted constraint can become a positive, almost-zero-volume tet.
        // Keep the input unchanged when double coordinates cannot represent
        // the required incidence exactly.
        const bool retainVertex=moved>=0&&target.first<0&&(point(moved)-destination).squaredNorm()==0;
        if(!retainVertex&&dt::GEOM_FUNC::orient3d(mesh.V[reference[0]].data(),mesh.V[reference[1]].data(),
            mesh.V[reference[2]].data(),exactPoint.data())!=0)return false;
        if(target.first>=0)for(int k=0;k<3;++k){const int j=(k+1)%3;
            double a[2]={mesh.V[target.first][k],mesh.V[target.first][j]};
            double b[2]={mesh.V[target.second][k],mesh.V[target.second][j]};
            double p[2]={destination[k],destination[j]};
            if(dt::GEOM_FUNC::orient2d(a,b,p)!=0)return false;
        }
        const bool added=moved<0;const int n=added?static_cast<int>(mesh.V.size()):moved;
        const Vec old=added?destination:point(n);
        std::set<int> affected;
        if(!added)affected.insert(incident[n].begin(),incident[n].end());
        affected.insert(master);
        for(const Edge& e:{target,source})if(e.first>=0){auto it=edges.find(e);if(it==edges.end())return false;affected.insert(it->second.begin(),it->second.end());}
        struct RestorePoint {
            dt::Mesh& mesh;int n;Vec old;bool added,committed;
            ~RestorePoint(){if(!committed){if(added)mesh.V.pop_back();else for(int k=0;k<3;++k)mesh.V[n][k]=old[k];}}
        };
        if(added)mesh.V.push_back({{destination.x(),destination.y(),destination.z()}});
        else for(int k=0;k<3;++k)mesh.V[n][k]=destination[k];
        RestorePoint restore{mesh,n,old,added,false};
        std::vector<Face> replacement;std::vector<int> parents;
        for(int t:affected){const Face& f=mesh.F[t];int split=-1;
            if(std::find(f.begin(),f.begin()+3,n)==f.begin()+3)
                for(int k=0;k<3;++k)if(edgeKey(f[k],f[(k+1)%3])==target||edgeKey(f[k],f[(k+1)%3])==source)split=k;
            const size_t first=replacement.size();
            if(split>=0){int a=f[split],b=f[(split+1)%3],c=f[(split+2)%3];replacement.push_back({{a,n,c,f[3]}});replacement.push_back({{n,b,c,f[3]}});}
            else if(t==master&&std::find(f.begin(),f.begin()+3,n)==f.begin()+3){for(int k=0;k<3;++k)replacement.push_back({{f[k],f[(k+1)%3],n,f[3]}});}
            else replacement.push_back(f);
            parents.resize(replacement.size(),t);
            Vec p[3];for(int k=0;k<3;++k)p[k]=!added&&f[k]==n?old:point(f[k]);
            const Vec before=(p[1]-p[0]).cross(p[2]-p[0]);const double size=before.norm();
            const double scale=std::max({(p[1]-p[0]).squaredNorm(),(p[2]-p[0]).squaredNorm(),(p[2]-p[1]).squaredNorm()});
            for(size_t j=first;j<replacement.size();++j){const Vec after=normal(replacement[j]);const double norm=after.norm();
                if(!(norm>1e-8*scale)||after.dot(before)<turnCos*norm*size)return false;
            }
        }
        std::map<std::array<int,3>,size_t> unique;
        std::vector<unsigned char> discard(replacement.size(),0);
        int duplicates=0;
        for(size_t i=0;i<replacement.size();++i){auto inserted=unique.emplace(key(replacement[i]),i);if(!inserted.second){
            if(replacement[inserted.first->second][3]!=replacement[i][3])return false;
            discard[i]=1;++duplicates;
        }}
        for(size_t i=0;i<replacement.size();++i)if(!discard[i]){
            std::vector<size_t> nearby;tree.query(box(replacement[i]),nearby);
            for(size_t id:nearby){const int t=static_cast<int>(id);if(affected.count(t)||!live(t))continue;
                if(key(replacement[i])==key(mesh.F[t])){
                    if(replacement[i][3]!=mesh.F[t][3])return false;
                    discard[i]=1;++duplicates;break;
                }
                if(intersects(replacement[i],mesh.F[t],n,old,false,false)&&
                   !intersects(mesh.F[parents[i]],mesh.F[t],n,old,true,true))return false;
            }
            if(discard[i])continue;
            for(size_t j=0;j<i;++j)if(!discard[j]&&intersects(replacement[i],replacement[j],n,old,false,false))return false;
        }
        // A contact alone can pinch a still-separated sheet into a sliver.
        // Commit only a complete local overlap removal, not intermediate
        // scaffolding that depends on later, possibly rejected transactions.
        if(!duplicates)return false;
        std::map<int,int> remaining;
        for(int t:affected){const int id=mesh.F[t][3];if(!remaining.count(id))remaining[id]=faceCounts[id];--remaining[id];}
        for(size_t i=0;i<replacement.size();++i)if(!discard[i])++remaining[replacement[i][3]];
        for(const auto& count:remaining)if(count.second<=0)return false;
        const Face masterFace=mesh.F[master];
        std::vector<int> sources;for(int p:parents)sources.push_back(sourceFace[p]);
        for(int t:affected)detach(t);
        if(added){incident.emplace_back();priority.push_back({1,1});fixed.push_back(0);anchored.push_back(1);
            if(!mesh.pointSize.empty()){
                double size=DBL_MAX;for(int j=0;j<3;++j)for(int node:{masterFace[j],mesh.F[slave][j]}){double s=mesh.pointSize[node];if(s>0)size=std::min(size,s);}mesh.pointSize.push_back(size==DBL_MAX?0:size);
            }
        }
        anchored[n]=1;for(int j=0;j<3;++j)anchored[masterFace[j]]=1;
        auto slot=affected.begin();
        for(size_t i=0;i<replacement.size();++i)if(!discard[i]){
            int t;if(slot!=affected.end()){t=*slot++;mesh.F[t]=replacement[i];}else{t=static_cast<int>(mesh.F.size());mesh.F.push_back(replacement[i]);}
            if(sourceFace.size()<=static_cast<size_t>(t))sourceFace.resize(t+1);
            sourceFace[t]=sources[i];attach(t);enqueue(t);
        }
        while(slot!=affected.end()){mesh.F[*slot][0]=-1;++slot;}
        for(const Edge& e:{target,source})if(imprintEdges.erase(e)){
            imprintEdges.insert(edgeKey(e.first,n));imprintEdges.insert(edgeKey(n,e.second));
        }
        for(const Edge& e:{target,source}){auto it=segments.find(e);if(it==segments.end())continue;
            fixed[n]=1;
            const auto indices=it->second;segments.erase(it);
            for(int index:indices){auto segment=mesh.S[index];mesh.S[index]={{segment[0],n,segment[2]}};
                segments[edgeKey(segment[0],n)].push_back(index);segments[edgeKey(n,segment[1])].push_back(static_cast<int>(mesh.S.size()));mesh.S.push_back({{n,segment[1],segment[2]}});
            }
        }
        restore.committed=true;++patches;changedIDs.insert(masterFace[3]);
        return true;
    }

    bool attempt(int master,int slave,double tolerance) {
        const Face a=mesh.F[master],b=mesh.F[slave];
        const Vec origin=point(a[0]),normalA=normal(a).normalized();
        const Vec u=(point(a[1])-origin).normalized(),v=normalA.cross(u);
        const double tiny=1e-10*std::max(length(a),length(b));
        std::pair<int,int> masterPriority=std::max({priority[a[0]],priority[a[1]],priority[a[2]]});
        std::array<int,3> order={{b[0],b[1],b[2]}};
        std::stable_sort(order.begin(),order.end(),[&](int x,int y){return priority[x]<priority[y];});
        for(int n:order){if(std::find(a.begin(),a.begin()+3,n)!=a.begin()+3)continue;
            const Vec p=point(n);Vec destination=p-normalA*(p-origin).dot(normalA);Edge target(-1,-1);
            bool inside=true;
            for(int k=0;k<3;++k){const Vec e=point(a[(k+1)%3])-point(a[k]);if(e.cross(destination-point(a[k])).dot(normalA)<=tiny*e.norm())inside=false;}
            // First preserve the existing feature/corner position. If the
            // complete fold can be replaced by one nonoverlapping triangle
            // fan, no projection or change to incident feature curves is needed.
            if(inside&&(p-destination).norm()<=tolerance&&
                transact(master,slave,n,p,Edge(-1,-1),Edge(-1,-1)))return true;
            double best=DBL_MAX;
            for(int k=0;k<3;++k){const Vec start=point(a[k]),delta=point(a[(k+1)%3])-start;const double t=(p-start).dot(delta)/delta.squaredNorm();
                if(!(t>1e-8&&t<1-1e-8))continue;const Vec q=start+t*delta;const double distance=(q-p).norm();
                if(distance<=tolerance&&distance<best){best=distance;destination=q;target=edgeKey(a[k],a[(k+1)%3]);}
            }
            if(target.first<0&&!inside)continue;
            const double distance=(destination-p).norm();
            if(distance>tolerance||(!std::isfinite(distance)))continue;
            if(distance>tiny&&(fixed[n]||anchored[n]||priority[n]>masterPriority))continue;
            if(distance<=tiny)destination=p;
            if(transact(master,slave,n,destination,target,Edge(-1,-1)))return true;
        }
        // Partial overlaps can have no vertex inside the other triangle.
        // Split both crossing edges in one transaction, accepting only a
        // complete local overlap removal without moving corner nodes.
        for(int i=0;i<3;++i)for(int j=0;j<3;++j){Edge ea=edgeKey(a[i],a[(i+1)%3]),eb=edgeKey(b[j],b[(j+1)%3]);
            // Artificial diagonals from a previous face split are not new
            // imprint constraints. Crossing them again recursively refines
            // the same overlap instead of aligning the original sheets.
            if(!imprintEdges.count(ea)||!imprintEdges.count(eb))continue;
            if(ea.first==eb.first||ea.first==eb.second||ea.second==eb.first||ea.second==eb.second||segments.count(eb))continue;
            Vec pa=point(ea.first),da=point(ea.second)-pa,pb=point(eb.first),db=point(eb.second)-pb;
            Vec2 e(da.dot(u),da.dot(v)),f(db.dot(u),db.dot(v)),offset((pb-pa).dot(u),(pb-pa).dot(v));
            double den=cross2(e,f);if(std::abs(den)<=1e-12*e.norm()*f.norm())continue;
            double ta=cross2(offset,f)/den,tb=cross2(offset,e)/den;
            if(!(ta>1e-8&&ta<1-1e-8&&tb>1e-8&&tb<1-1e-8))continue;
            const Vec destination=pa+ta*da;if((destination-(pb+tb*db)).norm()>tolerance)continue;
            if(transact(master,slave,-1,destination,ea,eb))return true;
        }
        return false;
    }

    // Replace a nearly flat triangle by its two short boundary segments.
    // Keep all existing vertices: the opposite edge is subdivided through the
    // middle vertex in EVERY incident facet before the sliver is removed.
    bool mergeSliver(int t) {
        const Face old=mesh.F[t];int edge=0;double longest=0;
        for(int k=0;k<3;++k){const double d=(point(old[k])-point(old[(k+1)%3])).squaredNorm();if(d>longest){longest=d;edge=k;}}
        if(!(longest>0)||normal(old).norm()>0.01*longest)return false;
        const double flatness=normal(old).norm()/longest;
        const int a=old[edge],b=old[(edge+1)%3],n=old[(edge+2)%3];
        const Vec delta=point(b)-point(a);const double parameter=(point(n)-point(a)).dot(delta)/longest;
        if(!(parameter>1e-8&&parameter<1-1e-8))return false;
        const Edge e=edgeKey(a,b);if(segments.count(e))return false;
        const auto incidentEdge=edges.find(e);if(incidentEdge==edges.end())return false;
        const std::set<int> affected(incidentEdge->second.begin(),incidentEdge->second.end());
        bool sameID=false;for(int f:affected)if(f!=t&&mesh.F[f][3]==old[3])sameID=true;
        if(!sameID)return false;
        std::vector<Face> replacement;std::vector<int> parents;
        for(int f:affected)if(f!=t){const Face& before=mesh.F[f];
            if(std::find(before.begin(),before.begin()+3,n)!=before.begin()+3)return false;
            for(int k=0;k<3;++k)if(edgeKey(before[k],before[(k+1)%3])==e){
                replacement.push_back({{before[k],n,before[(k+2)%3],before[3]}});
                replacement.push_back({{n,before[(k+1)%3],before[(k+2)%3],before[3]}});
                parents.push_back(f);parents.push_back(f);break;
            }
        }
        std::vector<unsigned char> discard(replacement.size(),0);
        std::map<std::array<int,3>,size_t> unique;
        for(size_t i=0;i<replacement.size();++i){const Face& f=replacement[i];const Vec before=normal(mesh.F[parents[i]]),after=normal(f);
            if(!(after.norm()>1e-8*length(mesh.F[parents[i]])*length(mesh.F[parents[i]]))||after.dot(before)<turnCos*after.norm()*before.norm())return false;
            if(after.norm()<=1.05*flatness*length(f)*length(f))return false;
            auto inserted=unique.emplace(key(f),i);if(!inserted.second){if(replacement[inserted.first->second][3]!=f[3])return false;discard[i]=1;}
        }
        for(size_t i=0;i<replacement.size();++i)if(!discard[i]){
            std::vector<size_t> nearby;tree.query(box(replacement[i]),nearby);
            for(size_t id:nearby){int f=static_cast<int>(id);if(affected.count(f)||!live(f))continue;
                if(key(replacement[i])==key(mesh.F[f])){if(replacement[i][3]!=mesh.F[f][3])return false;discard[i]=1;break;}
                if(intersects(replacement[i],mesh.F[f],-1,Vec::Zero(),false,false)&&
                    !intersects(mesh.F[parents[i]],mesh.F[f],-1,Vec::Zero(),false,false))return false;
            }
            if(discard[i])continue;
            for(size_t j=0;j<i;++j)if(!discard[j]&&intersects(replacement[i],replacement[j],-1,Vec::Zero(),false,false))return false;
        }
        std::map<int,int> remaining;
        for(int f:affected){int id=mesh.F[f][3];if(!remaining.count(id))remaining[id]=faceCounts[id];--remaining[id];}
        for(size_t i=0;i<replacement.size();++i)if(!discard[i])++remaining[replacement[i][3]];
        for(const auto& entry:remaining)if(entry.second<=0)return false;
        std::vector<int> sources;for(int p:parents)sources.push_back(sourceFace[p]);
        for(int f:affected)detach(f);
        auto slot=affected.begin();
        for(size_t i=0;i<replacement.size();++i)if(!discard[i]){int f;
            if(slot!=affected.end()){f=*slot++;mesh.F[f]=replacement[i];}else{f=static_cast<int>(mesh.F.size());mesh.F.push_back(replacement[i]);}
            if(sourceFace.size()<=static_cast<size_t>(f))sourceFace.resize(f+1);
            sourceFace[f]=sources[i];attach(f);enqueue(f);
        }
        while(slot!=affected.end()){mesh.F[*slot][0]=-1;++slot;}
        imprintEdges.erase(e);imprintEdges.insert(edgeKey(a,n));imprintEdges.insert(edgeKey(n,b));
        for(int f:incident[n])enqueue(f);
        ++mergedSlivers;return true;
    }

    // Overlay the WHOLE pair before committing. Shared intersection vertices
    // replace both original sheets, including every incident edge face. This
    // avoids the intermediate contact/sliver left by a one-intersection edit.
    bool overlay(int master,int slave,double tolerance) {
        const Edge sourcePair=edgeKey(sourceFace[master],sourceFace[slave]);
        if(sourcePair.first==sourcePair.second||completedOverlays.count(sourcePair))return false;
        const Face original[2]={mesh.F[master],mesh.F[slave]};
        const Vec origin=point(original[0][0]);
        const Vec u=(point(original[0][1])-origin).normalized();
        const Vec axis=normal(original[0]).normalized(),v=axis.cross(u);
        const double scale=std::max(length(original[0]),length(original[1]));
        const double eps=1e-10*std::max(scale,std::sqrt(std::max(normal(mesh.F[master]).norm(),normal(mesh.F[slave]).norm())));
        const double areaEps=1e-8*std::min(normal(mesh.F[master]).norm(),normal(mesh.F[slave]).norm());
        const size_t oldSize=mesh.V.size();
        struct RestoreVertices {
            dt::Mesh& mesh;size_t size;bool committed;
            ~RestoreVertices(){if(!committed)mesh.V.resize(size);}
        } restore{mesh,oldSize,false};
        std::vector<int> nodes;
        std::vector<std::array<double,2>> xy;
        std::map<int,int> weld;
        auto project=[&](int n){const Vec p=point(n)-origin;return Vec2(p.dot(u),p.dot(v));};
        auto at=[&](int n){return Vec2(xy[n][0],xy[n][1]);};
        int corners[2][3];
        for(int side=0;side<2;++side)for(int k=0;k<3;++k){
            const int n=original[side][k];const Vec2 p=project(n);int found=-1;
            for(int j=0;j<static_cast<int>(nodes.size());++j)
                if(nodes[j]==n||((at(j)-p).norm()<=tolerance&&(point(nodes[j])-point(n)).norm()<=tolerance)){found=j;break;}
            if(found<0){found=static_cast<int>(nodes.size());nodes.push_back(n);xy.push_back({{p.x(),p.y()}});}
            else if(nodes[found]!=n){
                int keep=nodes[found],drop=n;
                if((point(keep)-point(drop)).norm()>tolerance)return false;
                if(fixed[drop]>fixed[keep]||(fixed[drop]==fixed[keep]&&priority[drop]>priority[keep]))std::swap(keep,drop);
                if(fixed[drop])return false;
                weld[drop]=keep;nodes[found]=keep;
                const Vec2 q=project(keep);xy[found]={{q.x(),q.y()}};
            }
            corners[side][k]=found;
        }
        auto mapped=[&](int n){auto it=weld.find(n);while(it!=weld.end()){n=it->second;it=weld.find(n);}return n;};
        for(int& n:nodes)n=mapped(n);
        for(int side=0;side<2;++side)
            if(corners[side][0]==corners[side][1]||corners[side][1]==corners[side][2]||corners[side][2]==corners[side][0])return false;
        // Every polygon intersection is interned once. Even when a rounded
        // coordinate is not exactly on the old line, ALL incident triangles
        // use the same subdivided edge, so no independent micro-gap remains.
        auto intern=[&](const Vec2& p)->int {
            for(int j=0;j<static_cast<int>(nodes.size());++j)if((at(j)-p).norm()<=eps)return j;
            Vec destination;bool found=false;
            for(int side=0;side<2&&!found;++side)for(int k=0;k<3;++k){
                int ia=corners[side][k],ib=corners[side][(k+1)%3];const Vec2 a=at(ia),e=at(ib)-a;
                const double t=(p-a).dot(e)/e.squaredNorm();
                if(t>=-1e-9&&t<=1+1e-9&&std::abs(cross2(e,p-a))<=eps*e.norm()){
                    destination=point(nodes[ia])+std::max(0.0,std::min(1.0,t))*(point(nodes[ib])-point(nodes[ia]));found=true;break;
                }
            }
            if(!found)return -1;
            const int n=static_cast<int>(mesh.V.size());mesh.V.push_back({{destination.x(),destination.y(),destination.z()}});
            nodes.push_back(n);xy.push_back({{p.x(),p.y()}});return static_cast<int>(nodes.size())-1;
        };
        using Loop=std::vector<int>;
        bool valid=true;
        auto clip=[&](const Loop& polygon,int a,int b,bool inside)->Loop {
            Loop result;const Vec2 start=at(a),dir=at(b)-start;const double band=eps*dir.norm();
            for(size_t k=0;k<polygon.size();++k){
                int i=polygon[k],j=polygon[(k+1)%polygon.size()];double di=cross2(dir,at(i)-start),dj=cross2(dir,at(j)-start);
                if(std::abs(di)<=band)di=0;if(std::abs(dj)<=band)dj=0;
                if(!inside){di=-di;dj=-dj;}
                if(di>=0)result.push_back(i);
                if((di<0&&dj>0)||(di>0&&dj<0)){
                    int n=intern(at(i)+(at(j)-at(i))*(di/(di-dj)));if(n<0){valid=false;return Loop();}result.push_back(n);
                }
            }
            result.erase(std::unique(result.begin(),result.end()),result.end());
            if(result.size()>1&&result.front()==result.back())result.pop_back();return result;
        };
        auto area=[&](const Loop& p){double sum=0;for(size_t k=0;k<p.size();++k)sum+=cross2(at(p[k]),at(p[(k+1)%p.size()]));return sum;};
        Loop triangles[2];
        for(int side=0;side<2;++side){triangles[side].assign(corners[side],corners[side]+3);if(area(triangles[side])<0)std::reverse(triangles[side].begin(),triangles[side].end());}
        std::vector<Loop> polygons;std::vector<int> owners;
        Loop intersection;
        for(int side=0;side<2;++side){Loop rest=triangles[side];const Loop& cutter=triangles[1-side];
            for(int k=0;k<3&&!rest.empty();++k){Loop outside=clip(rest,cutter[k],cutter[(k+1)%3],false);
                if(outside.size()>=3&&area(outside)>areaEps){polygons.push_back(outside);owners.push_back(side);}
                rest=clip(rest,cutter[k],cutter[(k+1)%3],true);
            }
            if(!valid)return false;if(side==0)intersection=rest;
        }
        if(intersection.size()<3||area(intersection)<=areaEps)return false;
        polygons.push_back(intersection);owners.push_back(0);

        // Locate all subdivisions of original edges, not only the two faces.
        std::map<Edge,std::vector<std::pair<double,int>>> splits;
        std::set<int> affected{master,slave};
        for(const auto& entry:weld)affected.insert(incident[entry.first].begin(),incident[entry.first].end());
        for(int side=0;side<2;++side)for(int k=0;k<3;++k){
            int na=original[side][k],nb=original[side][(k+1)%3];Edge edge=edgeKey(na,nb);
            const Vec2 a=project(mapped(edge.first)),e=project(mapped(edge.second))-a;
            for(int j=0;j<static_cast<int>(nodes.size());++j){
                if(nodes[j]==mapped(na)||nodes[j]==mapped(nb))continue;
                const double t=(at(j)-a).dot(e)/e.squaredNorm();
                if(t<=1e-9||t>=1-1e-9||std::abs(cross2(e,at(j)-a))>eps*e.norm())continue;
                const Vec onEdge=point(mapped(edge.first))+t*(point(mapped(edge.second))-point(mapped(edge.first)));
                const double distance=(point(nodes[j])-onEdge).norm();
                if(distance>tolerance|| (segments.count(edge)&&distance>eps))return false;
                splits[edge].push_back({t,nodes[j]});
            }
            auto split=splits.find(edge);if(split!=splits.end()){
                auto& list=split->second;std::sort(list.begin(),list.end());list.erase(std::unique(list.begin(),list.end()),list.end());
                const auto& fs=edges.find(edge)->second;affected.insert(fs.begin(),fs.end());
            }
        }
        std::vector<Face> replacement;std::vector<int> parents;
        // Ear clipping retains collinear boundary vertices, keeping neighbouring
        // facets conforming even when a polygon has more than one edge split.
        auto triangulate=[&](Loop polygon,int parent)->bool {
            polygon.erase(std::unique(polygon.begin(),polygon.end()),polygon.end());
            if(polygon.size()>1&&polygon.front()==polygon.back())polygon.pop_back();
            if(polygon.size()<3)return true;
            const Face& f=mesh.F[parent];const Vec n=normal(f).normalized(),o=point(f[0]);
            const Vec ex=(point(f[1])-o).normalized(),ey=n.cross(ex);
            auto proj=[&](int i){const Vec d=point(i)-o;return Vec2(d.dot(ex),d.dot(ey));};
            double signedArea=0;for(size_t k=0;k<polygon.size();++k)signedArea+=cross2(proj(polygon[k]),proj(polygon[(k+1)%polygon.size()]));
            if(signedArea<0)std::reverse(polygon.begin(),polygon.end());
            while(polygon.size()>3){int best=-1;double bestArea=0;
                for(int k=0;k<static_cast<int>(polygon.size());++k){
                    const int prev=(k+polygon.size()-1)%polygon.size(),next=(k+1)%polygon.size();
                    const Vec2 a=proj(polygon[prev]),b=proj(polygon[k]),c=proj(polygon[next]);const double ar=cross2(b-a,c-a);
                    if(ar<=bestArea)continue;bool occupied=false;
                    for(int j=0;j<static_cast<int>(polygon.size());++j)if(j!=prev&&j!=k&&j!=next){const Vec2 p=proj(polygon[j]);
                        if(cross2(b-a,p-a)>=-areaEps&&cross2(c-b,p-b)>=-areaEps&&cross2(a-c,p-c)>=-areaEps){occupied=true;break;}}
                    if(!occupied){best=k;bestArea=ar;}
                }
                if(best<0)return false;
                const int prev=(best+polygon.size()-1)%polygon.size(),next=(best+1)%polygon.size();
                replacement.push_back({{polygon[prev],polygon[best],polygon[next],f[3]}});parents.push_back(parent);polygon.erase(polygon.begin()+best);
            }
            replacement.push_back({{polygon[0],polygon[1],polygon[2],f[3]}});parents.push_back(parent);return true;
        };
        for(size_t i=0;i<polygons.size();++i){Loop ids;
            const Loop& polygon=polygons[i];
            for(size_t k=0;k<polygon.size();++k){int a=polygon[k],b=polygon[(k+1)%polygon.size()];ids.push_back(nodes[a]);
                const Vec2 start=at(a),edge=at(b)-start;std::vector<std::pair<double,int>> onEdge;
                for(int j=0;j<static_cast<int>(nodes.size());++j)if(j!=a&&j!=b){const double t=(at(j)-start).dot(edge)/edge.squaredNorm();
                    if(t>1e-9&&t<1-1e-9&&std::abs(cross2(edge,at(j)-start))<=eps*edge.norm())onEdge.push_back({t,nodes[j]});
                }
                std::sort(onEdge.begin(),onEdge.end());for(const auto& node:onEdge)ids.push_back(node.second);
            }
            if(!triangulate(ids,owners[i]==0?master:slave))return false;
        }
        for(int t:affected)if(t!=master&&t!=slave){const Face& f=mesh.F[t];Loop boundary;
            for(int k=0;k<3;++k){int a=f[k],b=f[(k+1)%3];boundary.push_back(mapped(a));auto split=splits.find(edgeKey(a,b));
                if(split!=splits.end()){const auto& list=split->second;
                    if(a<b)for(const auto& p:list)boundary.push_back(p.second);
                    else for(auto it=list.rbegin();it!=list.rend();++it)boundary.push_back(it->second);
                }
            }
            if(!triangulate(boundary,t))return false;
        }
        // Intersections near a corner can create a tiny connector triangle.
        // Merge that connector inside the same transaction, before judging
        // its unreliable normal. Only newly created, non-feature points may
        // be removed; all incident faces are already in this transaction.
        std::set<int> protectedNew;
        for(const auto& split:splits)if(segments.count(split.first))for(const auto& p:split.second)protectedNew.insert(p.second);
        std::map<int,int> localWeld;
        auto localRoot=[&](int n){auto it=localWeld.find(n);while(it!=localWeld.end()){n=it->second;it=localWeld.find(n);}return n;};
        bool changed=true;
        while(changed){changed=false;
            for(const Face& f:replacement){int n[3]={localRoot(f[0]),localRoot(f[1]),localRoot(f[2])};
                if(n[0]==n[1]||n[1]==n[2]||n[2]==n[0])continue;
                const Face current={{n[0],n[1],n[2],f[3]}};
                if(normal(current).norm()>=0.005*length(current)*length(current))continue;
                int a=-1,b=-1;double best=tolerance;
                for(int k=0;k<3;++k){int keep=n[k],drop=n[(k+1)%3];if(keep>drop)std::swap(keep,drop);
                    const double distance=(point(keep)-point(drop)).norm();
                    if(drop>=static_cast<int>(oldSize)&&!protectedNew.count(drop)&&distance<best){a=keep;b=drop;best=distance;}
                }
                if(a>=0){localWeld[b]=a;changed=true;break;}
            }
        }
        for(auto& f:replacement)for(int k=0;k<3;++k)f[k]=localRoot(f[k]);
        for(auto& split:splits)for(auto& p:split.second)p.second=localRoot(p.second);
        std::map<std::array<int,3>,size_t> unique;
        std::vector<unsigned char> discard(replacement.size(),0);
        for(size_t i=0;i<replacement.size();++i){const Face& f=replacement[i];const Face& before=mesh.F[parents[i]];
            if(f[0]==f[1]||f[1]==f[2]||f[2]==f[0]){discard[i]=1;continue;}
            const Vec afterNormal=normal(f),beforeNormal=normal(before);const double norm=afterNormal.norm();
            if(!(norm>1e-8*normal(mesh.F[parents[i]]).norm())||afterNormal.dot(beforeNormal)<turnCos*norm*beforeNormal.norm())return false;
            if(norm<0.001*length(f)*length(f))return false;
            auto inserted=unique.emplace(key(f),i);if(!inserted.second){if(replacement[inserted.first->second][3]!=f[3])return false;discard[i]=1;}
        }
        for(size_t i=0;i<replacement.size();++i)if(!discard[i]){
            std::vector<size_t> nearby;tree.query(box(replacement[i]),nearby);
            for(size_t id:nearby){int t=static_cast<int>(id);if(affected.count(t)||!live(t))continue;
                if(key(replacement[i])==key(mesh.F[t])){if(replacement[i][3]!=mesh.F[t][3])return false;discard[i]=1;break;}
                if(intersects(replacement[i],mesh.F[t],-1,Vec::Zero(),false,false)&&
                    !intersects(mesh.F[parents[i]],mesh.F[t],-1,Vec::Zero(),false,false))return false;
            }
            if(discard[i])continue;
            for(size_t j=0;j<i;++j)if(!discard[j]&&intersects(replacement[i],replacement[j],-1,Vec::Zero(),false,false)&&
                !intersects(mesh.F[parents[i]],mesh.F[parents[j]],-1,Vec::Zero(),false,false))return false;
        }
        std::map<int,int> remaining;
        for(int t:affected){int id=mesh.F[t][3];if(!remaining.count(id))remaining[id]=faceCounts[id];--remaining[id];}
        for(size_t i=0;i<replacement.size();++i)if(!discard[i])++remaining[replacement[i][3]];
        for(const auto& entry:remaining)if(entry.second<=0)return false;
        // A collapsed connector must actually remove coverage. Otherwise the
        // same two sheets can repeatedly create and collapse the same point.
        double beforeArea=0,afterArea=0;
        for(int t:affected)if(mesh.F[t][3]==original[0][3])beforeArea+=std::abs(normal(mesh.F[t]).dot(axis));
        for(size_t i=0;i<replacement.size();++i)if(!discard[i]&&replacement[i][3]==original[0][3])afterArea+=std::abs(normal(replacement[i]).dot(axis));
        if(beforeArea-afterArea<std::max(areaEps,0.1*area(intersection)))return false;
        // All geometric checks have finished. Update the shared topology once.
        std::vector<int> sources;for(int p:parents)sources.push_back(sourceFace[p]);
        for(int t:affected)detach(t);
        double size=DBL_MAX;if(!mesh.pointSize.empty())for(const Face& f:original)for(int k=0;k<3;++k){double s=mesh.pointSize[f[k]];if(s>0)size=std::min(size,s);}
        while(incident.size()<mesh.V.size()){
            incident.emplace_back();priority.push_back({1,1});fixed.push_back(protectedNew.count(static_cast<int>(incident.size())-1)?1:0);anchored.push_back(1);
            if(!mesh.pointSize.empty())mesh.pointSize.push_back(size==DBL_MAX?0:size);
        }
        auto slot=affected.begin();
        for(size_t i=0;i<replacement.size();++i)if(!discard[i]){int t;
            if(slot!=affected.end()){t=*slot++;mesh.F[t]=replacement[i];}else{t=static_cast<int>(mesh.F.size());mesh.F.push_back(replacement[i]);}
            if(sourceFace.size()<=static_cast<size_t>(t))sourceFace.resize(t+1);
            sourceFace[t]=sources[i];attach(t);enqueue(t);
        }
        while(slot!=affected.end()){mesh.F[*slot][0]=-1;++slot;}
        for(const auto& split:splits){const Edge edge=split.first;std::vector<int> chain{mapped(edge.first)};for(const auto& p:split.second)chain.push_back(p.second);chain.push_back(mapped(edge.second));
            if(imprintEdges.erase(edge))for(size_t k=1;k<chain.size();++k)imprintEdges.insert(edgeKey(chain[k-1],chain[k]));
            auto it=segments.find(edge);if(it==segments.end())continue;const auto indices=it->second;segments.erase(it);
            for(int index:indices){const int tag=mesh.S[index][2];bool reverse=mesh.S[index][0]!=edge.first;
                for(size_t k=1;k<chain.size();++k){int a=chain[k-1],b=chain[k];if(reverse)std::swap(a,b);int s=index;
                    if(k>1){s=static_cast<int>(mesh.S.size());mesh.S.push_back({{a,b,tag}});}else mesh.S[s]={{a,b,tag}};
                    segments[edgeKey(a,b)].push_back(s);
                }
            }
        }
        for(const auto& entry:weld)removedVertices.insert(entry.first);
        for(const auto& entry:localWeld)removedVertices.insert(entry.first);
        if(!mesh.pointSize.empty()){
            for(const auto& entry:weld){const int keep=mapped(entry.second);const double size=mesh.pointSize[entry.first];
                if(size>0&&(mesh.pointSize[keep]<=0||mesh.pointSize[keep]>size))mesh.pointSize[keep]=size;
            }
            for(const auto& entry:localWeld){const int keep=localRoot(entry.second);const double size=mesh.pointSize[entry.first];
                if(size>0&&(mesh.pointSize[keep]<=0||mesh.pointSize[keep]>size))mesh.pointSize[keep]=size;
            }
        }
        completedOverlays.insert(sourcePair);
        restore.committed=true;++patches;changedIDs.insert(original[0][3]);return true;
    }

public:
    std::set<int> removedVertices;
    int mergedSlivers=0;
    int patches=0;
    std::set<int> changedIDs;
    SurfaceImprint(dt::Mesh& m,const dt::Args& args,const std::unordered_set<int>& locked)
        :mesh(m),tree(bounds(m)),incident(m.V.size()),priority(m.V.size()),fixed(m.V.size(),0),anchored(m.V.size(),0),originalVertices(m.V.size()) {
        for(size_t i=0;i<mesh.S.size();++i){const auto& s=mesh.S[i];segments[edgeKey(s[0],s[1])].push_back(static_cast<int>(i));fixed[s[0]]=fixed[s[1]]=1;}
        for(int n:locked)if(n>=0&&n<static_cast<int>(fixed.size()))fixed[n]=1;
        sourceFace.resize(mesh.F.size());
        for(int t=0;t<static_cast<int>(mesh.F.size());++t){sourceFace[t]=t;attach(t);}
        std::vector<std::vector<int>> features(m.V.size());const double featureCos=std::cos((180-args.adpangle)*PI/180.0);
        for(const auto& entry:edges){const auto& fs=entry.second;bool feature=fs.size()!=2;
            imprintEdges.insert(entry.first);
            if(!feature){
                const int a=entry.first.first,b=entry.first.second;int opposite[2]={-1,-1};
                for(int j=0;j<2;++j)for(int k=0;k<3;++k){int n=mesh.F[fs[j]][k];if(n!=a&&n!=b)opposite[j]=n;}
                const Vec e=point(b)-point(a),na=e.cross(point(opposite[0])-point(a)),nb=(point(opposite[1])-point(a)).cross(e);
                feature=mesh.F[fs[0]][3]!=mesh.F[fs[1]][3]||!(na.norm()*nb.norm()>0)||na.dot(nb)<featureCos*na.norm()*nb.norm();
            }
            if(feature){features[entry.first.first].push_back(entry.first.second);features[entry.first.second].push_back(entry.first.first);}
        }
        for(int n=0;n<static_cast<int>(m.V.size());++n){std::set<int> ids;for(int f:incident[n])ids.insert(mesh.F[f][3]);int level=features[n].empty()?1:3;
            if(features[n].size()==2){Vec a=point(features[n][0])-point(n),b=point(features[n][1])-point(n);if(a.dot(b)<-0.996*a.norm()*b.norm())level=2;}
            priority[n]={level,static_cast<int>(ids.size())};
        }
    }
    void run() {
        for(int t=0;t<static_cast<int>(mesh.F.size());++t)if(live(t))mergeSliver(t);
        for(int t=0;t<static_cast<int>(mesh.F.size());++t)enqueue(t);
        while(!pending.empty()){
            int a=pending.front();pending.pop();queued[a]=0;if(!live(a))continue;
            if(mergeSliver(a))continue;
            std::vector<size_t> nearby;tree.query(box(mesh.F[a],nearRatio*length(mesh.F[a])),nearby);std::sort(nearby.begin(),nearby.end());
            for(size_t index:nearby){int b=static_cast<int>(index);if(a==b||!live(b)||mesh.F[a][3]!=mesh.F[b][3])continue;
                double tolerance;if(!overlapping(a,b,tolerance))continue;
                auto rank=[&](int t){std::array<std::pair<int,int>,3> r={{priority[mesh.F[t][0]],priority[mesh.F[t][1]],priority[mesh.F[t][2]]}};std::sort(r.rbegin(),r.rend());return r;};
                int master=a,slave=b;const auto ra=rank(a),rb=rank(b);
                if(ra<rb||(ra==rb&&(normal(mesh.F[a]).norm()<normal(mesh.F[b]).norm()||(normal(mesh.F[a]).norm()==normal(mesh.F[b]).norm()&&b<a))))std::swap(master,slave);
                if(attempt(master,slave,tolerance)||attempt(slave,master,tolerance)||
                    overlay(master,slave,tolerance)||overlay(slave,master,tolerance))break;
            }
        }
        mesh.F.erase(std::remove_if(mesh.F.begin(),mesh.F.end(),[](const Face& f){return f[0]<0;}),mesh.F.end());
        std::vector<unsigned char> used(mesh.V.size(),0);
        for(const Face& f:mesh.F)for(int k=0;k<3;++k)used[f[k]]=1;
        for(const auto& s:mesh.S)for(int k=0;k<2;++k)used[s[k]]=1;
        for(size_t n=originalVertices;n<mesh.V.size();++n)if(!used[n])removedVertices.insert(static_cast<int>(n));
    }
};
}

int dt::DT::imprintNearbySurfaceFaces(Mesh& mesh,Args& args) {
    if(args.constrain!=0||!Nodes.empty()||!Elems.empty()||!mesh.T.empty()||mesh.V.empty()||mesh.F.empty()||!args.periodic_P.empty()||!lockF.empty())return 0;
    if(!mesh.pointSize.empty()&&mesh.pointSize.size()!=mesh.V.size())return 0;
    for(const auto& p:mesh.V)for(double value:p)if(!std::isfinite(value))return 0;
    for(const auto& f:mesh.F)for(int j=0;j<3;++j)if(f[j]<0||f[j]>=static_cast<int>(mesh.V.size()))return 0;
    for(const auto& f:mesh.F)if(f[0]==f[1]||f[1]==f[2]||f[2]==f[0])return 0;
    for(const auto& s:mesh.S)for(int j=0;j<2;++j)if(s[j]<0||s[j]>=static_cast<int>(mesh.V.size()))return 0;
    SurfaceImprint work(mesh,args,lockV);work.run();
    if(!work.removedVertices.empty()) {
        std::vector<int> mapping(mesh.V.size(),-1);int next=0;
        for(int n=0;n<static_cast<int>(mesh.V.size());++n)if(!work.removedVertices.count(n)){
            mapping[n]=next;mesh.V[next]=mesh.V[n];if(!mesh.pointSize.empty())mesh.pointSize[next]=mesh.pointSize[n];++next;
        }
        mesh.V.resize(next);if(!mesh.pointSize.empty())mesh.pointSize.resize(next);
        for(auto& f:mesh.F)for(int k=0;k<3;++k)f[k]=mapping[f[k]];
        for(auto& e:mesh.S)for(int k=0;k<2;++k)e[k]=mapping[e[k]];
        std::unordered_set<int> updatedLocks;for(int n:lockV)if(n>=0&&n<static_cast<int>(mapping.size())&&mapping[n]>=0)updatedLocks.insert(mapping[n]);lockV.swap(updatedLocks);
    }
    if(infolevel>0&&work.mergedSlivers)meshLogger->info("Merge near-degenerate surface faces: {}",work.mergedSlivers);
    if(infolevel>0&&work.patches)meshLogger->info("Imprint nearby surface faces: {} patches",work.patches);
    if(infolevel>=2)for(int id:work.changedIDs)meshLogger->debug("Imprint surface ID: {}",id);
    return work.patches;
}
