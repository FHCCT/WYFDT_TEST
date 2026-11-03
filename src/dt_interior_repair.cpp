#include "dt.h"
#include <cmath>
#include <limits>

using namespace dt;

// A bounded fallback only for quality-repair insertion. Planning and SUS trials
// touch a small private star, never the host's nodes, cells or quality cache.
bool DT::rescueInteriorPlan(BWPlan& plan, int seed, int& plansLeft, int& smoothLeft,
    double longestEdge, InteriorRepairAllowance* allowance, unsigned* rejection) {
    size_t& searchBudget=allowance?allowance->plans:topologyInteriorSearchBudget;
    size_t& smoothBudget=allowance?allowance->smooth:topologyInteriorSmoothBudget;
    const size_t cellBudget=allowance?allowance->cells:topologyInteriorCellBudget;
    auto reject = [&](unsigned flag) { if (rejection) *rejection |= flag; return false; };
    if (improve_Metric == 8 || (searchBudget == 0 && smoothBudget == 0)) return reject(1);
    const BWRequest original = plan.request;
    std::array<double,3> center{};
    std::vector<int> vertices;
    if (plan.status == BWStatus::Ready) {
        for (int t : plan.cavity) for (int v : Elems[t].form) vertices.push_back(v);
    } else for (int v : Elems[seed].form) vertices.push_back(v);
    std::sort(vertices.begin(),vertices.end());
    vertices.erase(std::unique(vertices.begin(),vertices.end()),vertices.end());
    for (int v : vertices) {
        if (v < 0 || v == ghost) return false;
        for (int k=0;k<3;++k) center[k] += Nodes[v].pt[k] / vertices.size();
    }
    const int region = Elems[seed].geo;
    auto minima = [&](double* a,double* b,double* c,double* d,double& angle,double& sj) {
        double low=0,high=0; std::vector<double> unused;
        if (!CalDihedral(a,b,d,c,low,high,unused,true)) return false;
        const double j=ScaledJacobian(a,b,c,d);
        if (!std::isfinite(low) || !std::isfinite(j)) return false;
        angle=std::min(angle,low); sj=std::min(sj,j); return true;
    };
    auto accept = [&](const std::array<double,3>& point,double floor,double oldA,double oldJ) {
        double a=DBL_MAX,j=DBL_MAX;
        double p[3]={point[0],point[1],point[2]};
        for (const BWFace& f : plan.faces) {
            double* x=Nodes[f.vertices[0]].pt;
            double* y=Nodes[f.vertices[1]].pt;
            double* z=Nodes[f.vertices[2]].pt;
            if (!(GEOM_FUNC::orient3d(p,x,z,y)>floor)) return reject(4);
            if (!minima(p,x,y,z,a,j)) return reject(2);
        }
        const bool improves = a >= oldA*(1-1e-10) && j >= oldJ*(1-1e-10) &&
            (a>oldA+1e-8 || j>oldJ+1e-8);
        // Spend extra points only when the entire new star leaves the extreme
        // repair band (the same 3-degree / .01-SJ band used by edge collapse).
        const bool useful = a >= ANGLE2RADIO(3.0) && j >= .01;
        if (improves && infolevel > 1)
            meshLogger->debug("Interior rescue candidate: seed={} angle={:.9g}->{:.9g} sj={:.9g}->{:.9g} useful={}",
                seed,RADIO2ANGLE(oldA),RADIO2ANGLE(a),oldJ,j,useful);
        if (!improves) { if (a<oldA*(1-1e-10)) reject(8); if (j<oldJ*(1-1e-10)) reject(16); }
        if (!useful) { if (a<ANGLE2RADIO(3.0)) reject(32); if (j<.01) reject(64); }
        return improves && useful;
    };
    auto test = [&](bool useSmooth) {
        if (plan.status != BWStatus::Ready || plan.faces.empty() || plan.faces.size()>64 ||
            std::find(plan.cavity.begin(),plan.cavity.end(),seed)==plan.cavity.end()) return reject(2);
        const size_t growth=plan.faces.size()>plan.cavity.size()?plan.faces.size()-plan.cavity.size():0;
        if (growth>cellBudget) return reject(1);
        double floor=DBL_MAX,oldA=DBL_MAX,oldJ=DBL_MAX;
        for (int t:plan.cavity) {
            if (isvirtualtet(t)||ishulltet(t)||Elems[t].geo!=region) return false;
            const double v=calVolume(t);
            if (!(v>0) || !std::isfinite(v)) return false;
            floor=std::min(floor,v);
            const int* f=Elems[t].form;
            if (!minima(Nodes[f[0]].pt,Nodes[f[1]].pt,Nodes[f[2]].pt,Nodes[f[3]].pt,oldA,oldJ)) return false;
        }
        floor=std::max(minVolume_bw,floor*.1);
        for (const BWFace& f:plan.faces) {
            if (f.geo!=region || !(GEOM_FUNC::orient3d(plan.request.point.data(),
                Nodes[f.vertices[0]].pt,Nodes[f.vertices[2]].pt,Nodes[f.vertices[1]].pt)>0)) return false;
        }
        // Keep an already good candidate even when SUS does not improve it.
        const bool rawAccepted=accept(plan.request.point,floor,oldA,oldJ);
        if (useSmooth && smoothLeft>0 && smoothBudget>0) {
            --smoothLeft; --smoothBudget;
            DT local;
            local.ghost=-1; local.improve_step=true; local.improve_Metric=SUS_METRIC;
            local.infolevel=0;
            local.Nodes.emplace_back(plan.request.point[0],plan.request.point[1],plan.request.point[2]);
            std::unordered_map<int,int> mapping;
            for (const BWFace& f:plan.faces) {
                int ids[3];
                for (int k=0;k<3;++k) {
                    const int v=f.vertices[k];
                    auto found=mapping.find(v);
                    if (found==mapping.end()) {
                        const int n=static_cast<int>(local.Nodes.size());
                        mapping[v]=n; local.Nodes.push_back(Nodes[v]); ids[k]=n;
                    } else ids[k]=found->second;
                }
                local.Elems.emplace_back(0,ids[0],ids[1],ids[2]);
                Elem& e=local.Elems.back(); e.geo=region;
                for (int k=0;k<3;++k) e.neig[k+1]=(int64_t(f.adjacent[k])<<2)|f.adjacentFace[k];
            }
            local.Nodes[0].tet=0;
            if (local.smooth_sus(0,1e-12)) {
                std::array<double,3> p{{local.Nodes[0].pt[0],local.Nodes[0].pt[1],local.Nodes[0].pt[2]}};
                if (accept(p,floor,oldA,oldJ)) {
                    plan.request.point=p;
                    plan.request.midpointEdge={{-1,-1}}; plan.request.sourceFace={{-1,-1,-1}};
                    ++topologyInteriorRescued;
                    return true;
                }
            }
        }
        if (rawAccepted) ++topologyInteriorRescued;
        return rawAccepted;
    };
    if (test(true)) return true;
    // Re-plan only a few inward offsets. planBW locates coordinate-only requests;
    // no host node is allocated. Never re-plan a successfully smoothed star.
    const double alphas[]={.02,.01,.005,.04,.08,.16};
    for (double alpha:alphas) {
        if (plansLeft==0 || searchBudget==0) break;
        BWRequest r=original;
        r.midpointEdge={{-1,-1}}; r.sourceFace={{-1,-1,-1}}; r.fixedEdgePoint=true;
        double displacement2=0;
        for (int k=0;k<3;++k) { r.point[k]+=alpha*(center[k]-original.point[k]);
            displacement2+=(r.point[k]-original.point[k])*(r.point[k]-original.point[k]); }
        if (displacement2>.0025*longestEdge*longestEdge) continue;
        --plansLeft; --searchBudget;
        if (planBW(*this,r,plan)!=BWStatus::Ready) continue;
        if (test(true)) return true;
    }
    return false;
}
