#include "dt.h"
#include <cmath>

using namespace dt;
namespace {
std::array<int,4> boundaryTetKey(const Elem& e) {
    std::array<int,4> key; std::copy(e.form,e.form+4,key.begin());
    std::sort(key.begin(),key.end()); return key;
}
bool boundaryTetQuality(DT& mesh,int t,double& angle,double& sj) {
    if (mesh.isDelEle(t)||mesh.ishulltet(t)||mesh.isvirtualtet(t)) return false;
    const int* v=mesh.Elems[t].form;
    for (int n:mesh.Elems[t].form) if (!mesh.isbndpnt(n)) return false;
    double high=0,low=0; std::vector<double> unused;
    mesh.CalDihedral(mesh.Nodes[v[0]].pt,mesh.Nodes[v[1]].pt,
        mesh.Nodes[v[3]].pt,mesh.Nodes[v[2]].pt,low,high,unused,true);
    angle=RADIO2ANGLE(low);
    sj=mesh.ScaledJacobian(mesh.Nodes[v[0]].pt,mesh.Nodes[v[1]].pt,
        mesh.Nodes[v[2]].pt,mesh.Nodes[v[3]].pt);
    return true;
}
}

int DT::repairFixedBoundaryTets(double strictAngle) {
    if (omp_in_parallel()) throw std::logic_error("Boundary tet repair must run serially");
    if (modifyBnd || !AniSol.empty() || improve_Metric==8) return 0;
    MeshStageLog timer(*this,"Fixed boundary tet repair",2);
    boundaryTetRepairRejections.clear();
    struct Candidate { TopologyCandidate cell; double angle; };
    std::vector<Candidate> candidates;
    const double trigger=std::min(3.0,std::max(0.0,strictAngle));
    for (int t=0;t<static_cast<int>(Elems.size());++t) {
        double angle,sj;
        if (boundaryTetQuality(*this,t,angle,sj) && (angle<trigger || sj<.01))
            candidates.push_back({topologyCandidate(t),angle});
    }
    std::sort(candidates.begin(),candidates.end(),[](const Candidate& a,const Candidate& b) {
        return a.angle==b.angle?a.cell.tet<b.cell.tet:a.angle<b.angle;
    });
    // Separate, one-shot tail allowance. Never replenish ordinary per-run
    // insertion/remesh budgets or recursively process newly created cells.
    InteriorRepairAllowance allowance;
    allowance.plans=96; allowance.smooth=48; allowance.cells=192;
    int success=0,attempts=0;
    for (const Candidate& candidate:candidates) {
        if (!topologyCandidateCurrent(candidate.cell)) continue;
        const int seed=candidate.cell.tet;
        double angle,sj;
        if (!boundaryTetQuality(*this,seed,angle,sj) || !(angle<trigger || sj<.01)) continue;
        const auto key=boundaryTetKey(Elems[seed]);
        if (attempts==24 || success==24 || allowance.plans==0 || allowance.smooth==0) {
            boundaryTetRepairRejections[key]=1; continue;
        }
        ++attempts;
        unsigned rejected=0; bool repaired=false,interiorEdge=false;
        // The inward direction includes the seed's physical face neighbors,
        // rather than the potentially thin original BW cavity alone.
        std::set<int> centerVertices;
        for (int v:Elems[seed].form) centerVertices.insert(v);
        auto canCross=[&](int t,int f) {
            int a,b,c,d; DFC(f,a,b,c,d);
            if (isBndTri(Elems[t].form[b],Elems[t].form[c],Elems[t].form[d])) return false;
            const int n=getNeig(t,f);
            return n>=0 && !isDelEle(n) && !ishulltet(n) && !isvirtualtet(n) && Elems[n].geo==Elems[seed].geo;
        };
        for (int f=0;f<4;++f) if (canCross(seed,f))
            for (int v:Elems[getNeig(seed,f)].form) centerVertices.insert(v);
        std::array<double,3> center{};
        for (int v:centerVertices) for (int k=0;k<3;++k) center[k]+=Nodes[v].pt[k]/centerVertices.size();
        double longest2=0;
        for (int e=0;e<6;++e) longest2=std::max(longest2,distance2(
            Nodes[candidate.cell.form[Egid[e][0]]].pt,Nodes[candidate.cell.form[Egid[e][1]]].pt));
        int smoothLeft=2;
        for (int e=0;e<6 && !repaired;++e) {
            const int a=candidate.cell.form[Egid[e][0]],b=candidate.cell.form[Egid[e][1]];
            if (isBndEdg(a,b)) continue;
            interiorEdge=true;
            if (!allowance.plans) { rejected|=1; break; }
            BWRequest r; r.info=1; r.trackAccess=false; r.acceptance=BWAcceptance::VolumeOnly;
            r.seeds.push_back(seed); r.midpointEdge={{a,b}};
            r.space=.5*(Nodes[a].space+Nodes[b].space);
            for (int k=0;k<3;++k) r.point[k]=.5*(Nodes[a].pt[k]+Nodes[b].pt[k]);
            BWPlan base; --allowance.plans;
            if (planBW(*this,r,base)!=BWStatus::Ready) { rejected|=2; continue; }
            std::vector<int> expanded=base.cavity;
            std::set<int> included(expanded.begin(),expanded.end());
            for (int ring=0;ring<2 && !repaired;++ring) {
                const size_t count=expanded.size();
                for (size_t i=0;i<count;++i) for (int f=0;f<4;++f) if (canCross(expanded[i],f)) {
                    const int n=getNeig(expanded[i],f);
                    if (included.insert(n).second) expanded.push_back(n);
                }
                if (expanded.size()>48) { rejected|=2; break; }
                if (!allowance.plans) { rejected|=1; break; }
                BWRequest proposal=r; proposal.seeds=expanded;
                proposal.midpointEdge={{-1,-1}}; proposal.fixedEdgePoint=true;
                double displacement2=0;
                for (int k=0;k<3;++k) {
                    proposal.point[k]+=.05*(center[k]-r.point[k]);
                    displacement2+=(proposal.point[k]-r.point[k])*(proposal.point[k]-r.point[k]);
                }
                if (displacement2>.0025*longest2) { rejected|=2; break; }
                BWPlan plan; --allowance.plans;
                if (planBW(*this,proposal,plan)!=BWStatus::Ready) { rejected|=2; continue; }
                // An explicitly enlarged cavity must not swallow an old vertex.
                // Ordinary BW relies on its Delaunay cavity for this property.
                std::set<int> retained;
                for (const BWFace& face:plan.faces) retained.insert(face.vertices.begin(),face.vertices.end());
                bool preservesVertices=true;
                for (int t:plan.cavity) for (int v:Elems[t].form)
                    if (!retained.count(v)) preservesVertices=false;
                if (!preservesVertices) { rejected|=512; continue; }
                int offsetsLeft=0,oneSmooth=std::min(1,smoothLeft);
                const int previousSmooth=oneSmooth;
                if (!rescueInteriorPlan(plan,seed,offsetsLeft,oneSmooth,std::sqrt(longest2),&allowance,&rejected)) {
                    smoothLeft-=previousSmooth-oneSmooth; continue;
                }
                smoothLeft-=previousSmooth-oneSmooth;
                const size_t growth=plan.faces.size()>plan.cavity.size()?plan.faces.size()-plan.cavity.size():0;
                std::vector<int> slots; slots.reserve(plan.faces.size());
                const int node=addNode(plan.request.point[0],plan.request.point[1],plan.request.point[2],r.space);
                try { for (size_t i=0;i<plan.faces.size();++i) slots.push_back(addElem()); }
                catch (...) { for (int t:slots) DelEle(t); DelNod(node); throw; }
                commitBW(*this,plan,node,slots); finishBW(*this,plan);
                allowance.cells-=growth; ++success; repaired=true;
            }
        }
        if (!repaired) boundaryTetRepairRejections[key]=interiorEdge?(rejected?rejected:128):256;
    }
    if (infolevel>0) meshLogger->info("Fixed boundary tet repair: {}/{}",success,attempts);
    return success;
}

void DT::reportBoundaryTetQuality(double strictAngle) {
    MeshStageLog timer(*this,"Boundary tet census",2);
    size_t total=0,poor=0,belowStrict=0,below5=0,below3=0,below1=0,tinySJ=0;
    size_t budget=0,cavity=0,volume=0,quality=0,extreme=0,boundary=0,oldVertex=0,stalled=0,notExtreme=0,notTried=0;
    size_t printed=0; int exportedIndex=0;
    for (int t=0;t<static_cast<int>(Elems.size());++t) {
        if (isDelEle(t)||ishulltet(t)||isvirtualtet(t)) continue;
        const int index=exportedIndex++;
        double angle,sj;
        if (!boundaryTetQuality(*this,t,angle,sj)) continue;
        ++total; belowStrict+=angle<strictAngle; below5+=angle<5; below3+=angle<3; below1+=angle<1; tinySJ+=sj<.01;
        if (!(angle<strictAngle || sj<1e-5)) continue;
        ++poor;
        const auto found=boundaryTetRepairRejections.find(boundaryTetKey(Elems[t]));
        unsigned flags=found==boundaryTetRepairRejections.end()?0:found->second;
        const char* reason="tail_not_attempted";
        if (flags) {
            budget+=(flags&1)!=0; cavity+=(flags&2)!=0; volume+=(flags&4)!=0;
            quality+=(flags&24)!=0; extreme+=(flags&96)!=0; boundary+=(flags&256)!=0;
            oldVertex+=(flags&512)!=0; stalled+=(flags&128)!=0;
            reason="candidate_rejected";
        } else if (angle>=std::min(3.0,std::max(0.0,strictAngle)) && sj>=.01) {
            ++notExtreme; reason="outside_extreme_tail_trigger";
        } else { ++notTried; }
        if (infolevel>1 && printed++<128) meshLogger->debug(
            "Boundary tet: index={} vertices=[{},{},{},{}] angle={:.9g} sj={:.9g} reason={} flags={}",
            index,Elems[t].form[0],Elems[t].form[1],Elems[t].form[2],Elems[t].form[3],angle,sj,reason,flags);
    }
    meshLogger->info("Boundary-vertex tets: total={} angle<{}={} poor={} angle<5={} angle<3={} angle<1={} SJ<0.01={}",
        total,strictAngle,belowStrict,poor,below5,below3,below1,tinySJ);
    meshLogger->info("Unrepaired boundary tets: outside-tail={} not-attempted={} budget={} cavity={} volume={} quality={} extreme-star={} constrained-edges={} old-vertex={} stalled={} (rejection causes may overlap)",
        notExtreme,notTried,budget,cavity,volume,quality,extreme,boundary,oldVertex,stalled);
}
