#include "dt.h"
#include "dt_geometry_guard.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>

namespace {
using Vec = Eigen::Vector3d;
const int faceNodes[4][3]={{1,2,3},{0,3,2},{0,1,3},{0,2,1}};

double unitQuality(double value) {
    return std::isfinite(value) ? std::max(0.0,std::min(1.0,value)) : 0.0;
}

struct TetFrame {
    std::array<Vec,4> absolute, point;
    int exponent=0;
    double volume6=0, orientation=1;
    bool init(double* const vertices[4]) {
        const Vec origin(vertices[0]);double length=0;
        for(int i=1;i<4;++i)
            length=std::max(length,(Vec(vertices[i])-origin).cwiseAbs().maxCoeff());
        if(!(length>0) || !std::isfinite(length))return false;
        std::frexp(length,&exponent);
        for(int i=0;i<4;++i)for(int k=0;k<3;++k)
            absolute[i][k]=std::scalbn(vertices[i][k],-exponent);
        for(int i=0;i<4;++i)point[i]=absolute[i]-absolute[0];
        const double signedVolume=dt::GEOM_FUNC::orient3d(absolute[0].data(),absolute[1].data(),
            absolute[3].data(),absolute[2].data());
        volume6=std::abs(signedVolume);orientation=signedVolume<0?-1:1;
        return volume6>0 && std::isfinite(volume6);
    }
};

// Circumcentres of almost coplanar, almost cocircular points amplify
// cancellation in the cross-product numerator. Use two-component arithmetic
// only for these ill-conditioned cells; the common path stays in double.
struct Wide {
    double hi,lo;
    Wide(double a=0,double b=0):hi(a),lo(b){}
};
Wide exactSum(double a,double b) {
    const double sum=a+b,part=sum-a;
    return Wide(sum,(a-(sum-part))+(b-part));
}
Wide operator+(Wide a,Wide b) {
    const Wide sum=exactSum(a.hi,b.hi);
    return exactSum(sum.hi,sum.lo+a.lo+b.lo);
}
Wide operator-(Wide a,Wide b){return a+Wide(-b.hi,-b.lo);}
Wide operator*(Wide a,Wide b) {
    const double product=a.hi*b.hi;
    return exactSum(product,std::fma(a.hi,b.hi,-product)+a.hi*b.lo+a.lo*b.hi+a.lo*b.lo);
}
Vec accurateRadiusNumerator(const TetFrame& tet) {
    std::array<std::array<Wide,3>,3> edge;
    std::array<Wide,3> length2;
    for(int i=0;i<3;++i)for(int k=0;k<3;++k) {
        edge[i][k]=exactSum(tet.absolute[i+1][k],-tet.absolute[0][k]);
        length2[i]=length2[i]+edge[i][k]*edge[i][k];
    }
    Vec result;
    for(int k=0;k<3;++k) {
        const int y=(k+1)%3,z=(k+2)%3;
        const Wide value=length2[2]*(edge[0][y]*edge[1][z]-edge[0][z]*edge[1][y])+
                         length2[1]*(edge[2][y]*edge[0][z]-edge[2][z]*edge[0][y])+
                         length2[0]*(edge[1][y]*edge[2][z]-edge[1][z]*edge[2][y]);
        result[k]=value.hi+value.lo;
    }
    return result;
}

// Equilateral-volume skewness uses a regular tetrahedron with the same
// circumradius. Return V/Vregular directly: 1-(1-q) loses tiny positive q.
double volumeQuality(const TetFrame& tet) {
    const Vec& ab=tet.point[1];const Vec& ac=tet.point[2];const Vec& ad=tet.point[3];
    const Vec numerator=tet.volume6<1e-6?accurateRadiusNumerator(tet):
        (ad.squaredNorm()*ab.cross(ac)+ac.squaredNorm()*ad.cross(ab)+ab.squaredNorm()*ac.cross(ad)).eval();
    double radius=numerator.norm()/(2*tet.volume6);
    double longest=0;
    for(int i=0;i<4;++i)for(int j=i+1;j<4;++j)
        longest=std::max(longest,(tet.point[i]-tet.point[j]).norm());
    radius=std::max(radius,.5*longest);
    if(!(radius>0) || !std::isfinite(radius))return 0;
    return unitQuality((3*std::sqrt(3.0)/16)*(tet.volume6/radius)/radius/radius);
}
}

double DT::calskewness(double* a,double* b,double* c,double* d) {
    double* points[4]={a,b,c,d};TetFrame tet;
    return tet.init(points)?1-volumeQuality(tet):1;
}

double DT::orthogonal(double* a,double* b,double* c,double* d,
    double* oppositeA,double* oppositeB,double* oppositeC,double* oppositeD) {
    double* points[4]={a,b,c,d};TetFrame tet;
    if(!tet.init(points))return 0;
    const double* opposite[4]={oppositeA,oppositeB,oppositeC,oppositeD};
    const Vec center=.25*(tet.point[0]+tet.point[1]+tet.point[2]+tet.point[3]);
    double quality=volumeQuality(tet);
    for(int face=0;face<4;++face) {
        const int* f=faceNodes[face];
        const Vec area=tet.orientation*(tet.point[f[1]]-tet.point[f[0]]).cross(tet.point[f[2]]-tet.point[f[0]]);
        const double areaLength=area.norm();
        const Vec faceVector=(tet.point[f[0]]+tet.point[f[1]]+tet.point[f[2]])/3-center;
        const double denominator=areaLength*faceVector.norm();
        if(!(denominator>0))return 0;
        // n dot (faceCenter-cellCenter) = 6V/4. Robust volume avoids
        // cancellation on nearly coplanar cells at large world coordinates.
        quality=std::min(quality,unitQuality(tet.volume6/(4*denominator)));
        if(!opposite[face])continue;
        Vec neighbor;
        for(int k=0;k<3;++k)neighbor[k]=std::scalbn(opposite[face][k],-tet.exponent);
        const double neighborHeight=-tet.orientation*dt::GEOM_FUNC::orient3d(
            tet.absolute[f[0]].data(),tet.absolute[f[1]].data(),tet.absolute[f[2]].data(),neighbor.data());
        const double separation=areaLength*(neighbor-tet.absolute[face]).norm();
        if(!(neighborHeight>0) || !(separation>0))return 0;
        // Neighbor and current cell share three vertices, so Cn-Cp=(E-A)/4.
        quality=std::min(quality,unitQuality((tet.volume6+neighborHeight)/separation));
    }
    return quality;
}

double* DT::orthogonalOppositeVertex(int cell,int face) {
    const int neighbor=getNeig(cell,face);
    if(neighbor<0 || neighbor>=static_cast<int>(Elems.size()) || isDelEle(neighbor) ||
        ishulltet(neighbor) || isvirtualtet(neighbor))return nullptr;
    const int node=getoppoP(cell,face);
    return node==ghost?nullptr:Nodes[node].pt;
}

double DT::orthogonalQuality(int cell,int movedNode,const double* position) {
    if(cell<0 || cell>=static_cast<int>(Elems.size()) || isDelEle(cell) || ishulltet(cell) || isvirtualtet(cell))return 0;
    auto point=[&](int node)->double* {
        return node==movedNode && position?const_cast<double*>(position):Nodes[node].pt;
    };
    const auto& f=Elems[cell].form;
    // Keep invalid mesh orientation distinct from roundoff in a valid cosine.
    if(!(dt::GEOM_FUNC::orient3d(point(f[0]),point(f[1]),point(f[3]),point(f[2]))>0))return -1;
    double* opposite[4]={nullptr,nullptr,nullptr,nullptr};
    for(int k=0;k<4;++k) {
        const int neighbor=getNeig(cell,k);
        if(neighbor<0 || neighbor>=static_cast<int>(Elems.size()) || isDelEle(neighbor) ||
            ishulltet(neighbor) || isvirtualtet(neighbor))continue;
        const int node=getoppoP(cell,k);
        if(node!=ghost)opposite[k]=point(node);
    }
    return orthogonal(point(f[0]),point(f[1]),point(f[2]),point(f[3]),
        opposite[0],opposite[1],opposite[2],opposite[3]);
}

int DT::smoothOrthogonalPoint(int node,double minimumShapeFloor,BoundarySmoothWorkspace& boundary,bool tail) {
    if(node<0 || node>=static_cast<int>(Nodes.size()) || node==ghost || isDelNod(node) ||
        isCornerpnt(node) || lockV.count(node) || periodic_P.count(node))return 0;
    const bool onBoundary=isbndpnt(node);
    std::vector<int> star;
    if(onBoundary) {
        if(!prepareBoundarySmoothPoint(node,boundary))return 0;
        star=boundary.star;
    } else findSphere(node,star);
    star.erase(std::remove_if(star.begin(),star.end(),[&](int t){
        return isDelEle(t) || ishulltet(t) || isvirtualtet(t);
    }),star.end());
    if(star.empty())return 0;
    std::vector<int> affected=star;
    const Vec origin(Nodes[node].pt);
    double scale=DBL_MAX;
    for(int t:star) {
        for(int n:Elems[t].form)if(n!=node)scale=std::min(scale,(Vec(Nodes[n].pt)-origin).norm());
        for(int face=0;face<4;++face) {
            const int neighbor=getNeig(t,face);
            if(neighbor>=0 && neighbor<static_cast<int>(Elems.size()) && !isDelEle(neighbor) &&
                !ishulltet(neighbor) && !isvirtualtet(neighbor))affected.push_back(neighbor);
        }
    }
    if(!(scale>0) || !std::isfinite(scale))return 0;
    std::sort(affected.begin(),affected.end());
    affected.erase(std::unique(affected.begin(),affected.end()),affected.end());
    struct Score { double minimum=1,sum=0,logarithm=0,minAngle=180,minJacobian=1; int bad=0; double deficit=0; };
    auto evaluate=[&](const Vec& position,Score& score,std::vector<double>* qualities,bool probing=false)->bool {
        score=Score{};
        for(int t:star) {
            double* p[4];
            for(int k=0;k<4;++k)p[k]=Elems[t].form[k]==node?
                const_cast<double*>(position.data()):Nodes[Elems[t].form[k]].pt;
            const double q=quality_sus(p[0],p[1],p[2],p[3]);
            if(!(q>0) || (!probing && !(q>=minimumShapeFloor)))return false;
            if(tail) {
                double low=0,high=0,jacobian=0;std::vector<double> angles;
                if(!CalDihedral(p[0],p[1],p[3],p[2],low,high,angles,true,&jacobian))return false;
                score.minAngle=std::min(score.minAngle,RADIO2ANGLE(low));
                score.minJacobian=std::min(score.minJacobian,jacobian);
            }
        }
        if(qualities)qualities->resize(affected.size());
        for(size_t i=0;i<affected.size();++i) {
            const double q=orthogonalQuality(affected[i],node,position.data());
            if(!std::isfinite(q) || q<0)return false;
            score.minimum=std::min(score.minimum,q);score.sum+=tail?std::min(q,.01):q;
            score.logarithm+=std::log((tail?std::min(q,.01):q)+1e-30);
            if(tail){score.bad+=q<.01;score.deficit+=std::max(0.,.01-q);}
            if(qualities)(*qualities)[i]=q;
        }
        return true;
    };
    Score before;
    if(!evaluate(origin,before,nullptr))return 0;
    std::array<Vec,3> basis={{Vec::UnitX(),Vec::UnitY(),Vec::UnitZ()}};
    int dimension=3;
    if(onBoundary) {
        if(isFacetpnt(node)) {
            const Vec normal=Vec(boundary.normal.data())/boundary.normalLength;
            basis[0]=normal.unitOrthogonal();basis[1]=normal.cross(basis[0]);dimension=2;
        } else {
            basis[0]=(Vec(Nodes[boundary.neighbors[0]].pt)-Vec(Nodes[boundary.neighbors[1]].pt)).normalized();
            dimension=1;
        }
    }
    Vec meanGradient=Vec::Zero(),logGradient=Vec::Zero(),minimumGradient=Vec::Zero();
    // Estimate directions in the positive-volume domain. Apply the actual
    // shape floor to accepted moves, rather than clipping finite differences.
    const double h=scale*1e-5;
    for(int k=0;k<dimension;++k) {
        Score plus,minus;
        const bool validPlus=evaluate(origin+h*basis[k],plus,nullptr,tail);
        const bool validMinus=evaluate(origin-h*basis[k],minus,nullptr,tail);
        if(!validPlus && !validMinus)continue;
        if(!validPlus)plus=before;
        if(!validMinus)minus=before;
        const double divisor=(validPlus && validMinus)?2*h:h;
        if(tail)minimumGradient+=basis[k]*((plus.minimum-minus.minimum)/divisor);
        meanGradient+=basis[k]*((plus.sum-minus.sum)/divisor);
        logGradient+=basis[k]*((plus.logarithm-minus.logarithm)/divisor);
    }
    std::array<Vec,4> directions={{logGradient,meanGradient,Vec::Zero(),Vec::Zero()}};
    if(tail) {
        directions[2]=minimumGradient;
        Vec centroid=Vec::Zero();std::set<int> neighbors;
        for(int t:star)for(int n:Elems[t].form)if(n!=node)neighbors.insert(n);
        for(int n:neighbors)centroid+=Vec(Nodes[n].pt)-origin;
        Vec projected=Vec::Zero();
        for(int k=0;k<dimension;++k)projected+=basis[k]*centroid.dot(basis[k]);
        directions[3]=projected;
    }
    bool found=false;Vec bestPosition=origin;Score bestScore;
    std::vector<double> bestQualities;
    std::vector<double> qualities;
    Vec previous=Vec::Zero();
    for(int proposal=0;proposal<(tail?4:2);++proposal) {
        const Vec& gradient=directions[proposal];
        if(!gradient.allFinite() || !(gradient.norm()>0))continue;
        const Vec direction=gradient.normalized();
        if((direction-previous).squaredNorm()<1e-24)continue;
        previous=direction;
        double step=(tail?1.:.25)*scale;
        for(int trial=0;trial<10;++trial,step*=.5) {
            const Vec position=origin+step*direction;
            if(onBoundary) {
                bool valid=true;
                for(const auto& f:boundary.faces) {
                    const Vec a(Nodes[f[0]].pt),b(Nodes[f[1]].pt),c(Nodes[f[2]].pt);
                    if(!dt::pointOnBoundaryPlane(a.data(),b.data(),c.data(),position.data())){valid=false;break;}
                    const Vec x=f[0]==node?position:a,y=f[1]==node?position:b,z=f[2]==node?position:c;
                    if(!((b-a).cross(c-a).dot((y-x).cross(z-x))>0)){valid=false;break;}
                }
                if(!valid)continue;
            }
            Score after;
            if(!evaluate(position,after,&qualities) || after.minimum+1e-12*before.minimum<before.minimum)continue;
            if(!tail && (!(after.sum>before.sum+1e-12*std::max(before.sum,1e-30)) ||
                !(after.logarithm>before.logarithm+1e-12*std::max(1.0,std::abs(before.logarithm)))))continue;
            if(tail) {
                // Use the ordinary A/J safety floors; local angles need not
                // all improve while repairing the worst orthogonal quality.
                if(after.bad>before.bad || after.minAngle+1e-10*before.minAngle<std::min(before.minAngle,1.0) ||
                    after.minJacobian+1e-10*before.minJacobian<std::min(before.minJacobian,.001) ||
                    !(after.bad<before.bad || after.minimum>before.minimum*(1+1e-10) ||
                      after.deficit<before.deficit-1e-12*std::max(before.deficit,1e-30)))continue;
                if(!found || after.bad<bestScore.bad || (after.bad==bestScore.bad &&
                    (after.minimum>bestScore.minimum || (after.minimum==bestScore.minimum && after.deficit<bestScore.deficit)))) {
                    found=true;bestPosition=position;bestScore=after;bestQualities=qualities;
                }
                break;
            }
            for(int k=0;k<3;++k)Nodes[node].pt[k]=position[k];
            // Commit the complete dependency stencil, not just incident cells.
            for(size_t i=0;i<affected.size();++i)Elems[affected[i]].q=qualities[i];
            return 1;
        }
    }
    if(found) {
        for(int k=0;k<3;++k)Nodes[node].pt[k]=bestPosition[k];
        for(size_t i=0;i<affected.size();++i)Elems[affected[i]].q=bestQualities[i];
        return 1;
    }
    return 0;
}

// Preview the complete quality dependency stencil without editing host cells.
// Face neighbors are affected even when none of their vertices move.
bool DT::acceptOrthogonalPatch(const std::vector<int>& removed,
    const std::vector<std::array<int,4>>& replacement,int newNode,const double* position) {
    if(improve_Metric!=8)return true;
    if(removed.empty() || replacement.empty())return false;
    using Face=std::array<int,3>;
    auto key=[](const int* v,int face){Face x;int j=0;for(int k=0;k<4;++k)if(k!=face)x[j++]=v[k];std::sort(x.begin(),x.end());return x;};
    auto physical=[&](int t){return t>=0 && t<int(Elems.size()) && !isDelEle(t) && !ishulltet(t) && !isvirtualtet(t);};
    auto point=[&](int n)->double* {
        if(n==newNode && position)return const_cast<double*>(position);
        return n>=0 && n<int(Nodes.size()) && n!=ghost && !isDelNod(n)?Nodes[n].pt:nullptr;
    };
    const std::set<int> old(removed.begin(),removed.end());
    if(old.size()!=removed.size())return false;
    std::set<int> outside;
    std::map<Face,std::pair<int,int>> boundary;
    long double oldVolume=0,newVolume=0;
    for(int t:old) {
        if(!physical(t))return false;
        oldVolume+=calVolume(t);
        for(int f=0;f<4;++f)if(!old.count(getNeig(t,f))) {
            const int n=getNeig(t,f);
            if(!boundary.emplace(key(Elems[t].form,f),std::make_pair(n,getNeigOrd(t,f))).second)return false;
            if(physical(n))outside.insert(n);
        }
    }
    struct Link { int cell[2]={-1,-1},face[2]={-1,-1},count=0; };
    std::map<Face,Link> faces;
    std::set<std::array<int,4>> unique;
    for(int i=0;i<int(replacement.size());++i) {
        const auto& v=replacement[i];auto sorted=v;std::sort(sorted.begin(),sorted.end());
        if(std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end() || !unique.insert(sorted).second)return false;
        for(int n:v)if(!point(n))return false;
        const double volume=dt::GEOM_FUNC::orient3d(point(v[0]),point(v[1]),point(v[3]),point(v[2]));
        if(!(volume>0) || !std::isfinite(volume))return false;
        newVolume+=volume;
        for(int f=0;f<4;++f) {
            auto& link=faces[key(v.data(),f)];if(link.count==2)return false;
            link.cell[link.count]=i;link.face[link.count++]=f;
        }
    }
    if(std::abs(oldVolume-newVolume)>1e-10L*oldVolume)return false;
    size_t exterior=0;
    for(const auto& face:faces)if(face.second.count==1) {
        if(!boundary.count(face.first))return false;++exterior;
    }
    if(exterior!=boundary.size())return false;
    struct Score { double minimum=1,deficit=0,angle=180,jacobian=1;int bad=0;bool valid=true; } before,after;
    auto add=[&](Score& score,const int* v,double* const opposite[4]) {
        double* p[4]={point(v[0]),point(v[1]),point(v[2]),point(v[3])};
        const double q=orthogonal(p[0],p[1],p[2],p[3],opposite[0],opposite[1],opposite[2],opposite[3]);
        if(!std::isfinite(q) || !(q>0)){score.valid=false;return;}
        score.minimum=std::min(score.minimum,q);score.bad+=q<.01;score.deficit+=std::max(0.,.01-q);
        double low=0,high=0,j=0;std::vector<double> angles;
        if(!CalDihedral(p[0],p[1],p[3],p[2],low,high,angles,true,&j)){score.valid=false;return;}
        score.angle=std::min(score.angle,RADIO2ANGLE(low));score.jacobian=std::min(score.jacobian,j);
    };
    for(int t:removed) {
        double* opposite[4];for(int f=0;f<4;++f)opposite[f]=orthogonalOppositeVertex(t,f);
        add(before,Elems[t].form,opposite);
    }
    for(int t:outside) {
        double* opposite[4];for(int f=0;f<4;++f)opposite[f]=orthogonalOppositeVertex(t,f);
        add(before,Elems[t].form,opposite);
        for(int f=0;f<4;++f)if(old.count(getNeig(t,f))) {
            const auto it=faces.find(key(Elems[t].form,f));
            if(it==faces.end() || it->second.count!=1)return false;
            const auto& link=it->second;
            opposite[f]=point(replacement[link.cell[0]][link.face[0]]);
        }
        add(after,Elems[t].form,opposite);
    }
    for(int i=0;i<int(replacement.size());++i) {
        const auto& v=replacement[i];double* opposite[4]={};
        for(int f=0;f<4;++f) {
            const auto& link=faces.at(key(v.data(),f));
            if(link.count==2) {
                const int side=link.cell[0]==i?1:0;
                opposite[f]=point(replacement[link.cell[side]][link.face[side]]);
            } else {
                const auto& b=boundary.at(key(v.data(),f));
                if(physical(b.first))opposite[f]=point(Elems[b.first].form[b.second]);
            }
        }
        add(after,v.data(),opposite);
    }
    return before.valid && after.valid && after.minimum+1e-12*before.minimum>=before.minimum &&
        after.angle+1e-10*before.angle>=std::min(before.angle,1.0) &&
        after.jacobian+1e-10*before.jacobian>=std::min(before.jacobian,.001) && after.bad<=before.bad &&
        (after.bad<before.bad || after.minimum>before.minimum*(1+1e-10) ||
         after.deficit<before.deficit-1e-12*std::max(before.deficit,1e-30));
}

void DT::repairOrthogonalTail(double insertionAngleDegrees,int rounds) {
    MeshStageLog stage(*this,"Orthogonal tail repair",1);
    double floor=DBL_MAX;
    for(int t=0;t<int(Elems.size());++t) {
        if(isDelEle(t) || ishulltet(t) || isvirtualtet(t))continue;
        const auto& f=Elems[t].form;
        floor=std::min(floor,quality_sus(Nodes[f[0]].pt,Nodes[f[1]].pt,Nodes[f[2]].pt,Nodes[f[3]].pt));
    }
    BoundarySmoothWorkspace workspace;
    for(int round=0;round<rounds;++round) {
        std::set<int> nodes;std::vector<TopologyCandidate> candidates;
        for(int t=0;t<int(Elems.size());++t)if(!isDelEle(t) && !ishulltet(t) && !isvirtualtet(t) && Elems[t].q<.01) {
            candidates.push_back(topologyCandidate(t));for(int n:Elems[t].form)nodes.insert(n);
        }
        if(candidates.empty())break;
        int success=0;
        for(int pass=0;pass<4;++pass) {
            int moved=0;for(int n:nodes)moved+=smoothOrthogonalPoint(n,floor,workspace,true)==1;
            success+=moved;if(!moved)break;
        }
        std::stable_sort(candidates.begin(),candidates.end(),[&](const TopologyCandidate& a,const TopologyCandidate& b){return Elems[a.tet].q<Elems[b.tet].q;});
        topologyRemeshTrialBudget=std::max(topologyRemeshTrialBudget,candidates.size());
        for(const auto& c:candidates)if(topologyCandidateCurrent(c) && Elems[c.tet].q<.01)
            success+=tryTopologyRepair(c,insertionAngleDegrees)==1;
        if(!success)break;
    }
}
