#pragma once
#include "dt.h"

namespace dt { namespace detail {
// Physical A/J quantities used only for point-level smoothing acceptance.
struct JointSmoothingQuality {
    double angle=0, jacobian=0, minAngle=180, minJacobian=1, minSus=1;
    int count=0, badAngle=0, badJacobian=0;
    bool valid=true;
};
inline JointSmoothingQuality jointSmoothingQuality(DT& mesh,
    const std::vector<int>& star, int node, const double* position=nullptr) {
    JointSmoothingQuality out;
    std::vector<double> angles;
    for (int t : star) {
        double* p[4];
        for (int k=0;k<4;++k)
            p[k]=mesh.Elems[t].form[k]==node && position
                ? const_cast<double*>(position) : mesh.Nodes[mesh.Elems[t].form[k]].pt;
        double lo,hi;
        const double j=mesh.ScaledJacobian(p[0],p[1],p[2],p[3]);
        const double q=mesh.quality_sus(p[0],p[1],p[2],p[3]);
        if (!mesh.CalDihedral(p[0],p[1],p[3],p[2],lo,hi,angles,true) ||
            !(lo>0) || !std::isfinite(lo) || !(j>0) || !std::isfinite(j) || !(q>0)) {
            out.valid=false;return out;
        }
        lo*=180/std::acos(-1.0);
        out.angle+=lo;out.jacobian+=j;
        out.minAngle=std::min(out.minAngle,lo);
        out.minJacobian=std::min(out.minJacobian,j);out.minSus=std::min(out.minSus,q);
        ++out.count;out.badAngle+=lo<.1;out.badJacobian+=j<1e-5;
    }
    if (out.count) {out.angle/=out.count;out.jacobian/=out.count;}
    else out.valid=false;
    return out;
}
inline bool jointSmoothingAccept(const JointSmoothingQuality& before,
    const JointSmoothingQuality& after) {
    return before.valid && after.valid && before.count==after.count &&
        after.minAngle+1e-10>=std::min(before.minAngle,1.0) &&
        after.minJacobian+1e-12>=std::min(before.minJacobian,.001) &&
        after.badAngle<=before.badAngle && after.badJacobian<=before.badJacobian &&
        after.angle+1e-11>=before.angle && after.jacobian+1e-13>=before.jacobian;
}
} }
