#pragma once
#include "dt.h"

namespace dt { namespace detail {
using MinimaxVector = Eigen::Vector3d;
using MinimaxMatrix = Eigen::Matrix3d;
struct MinimaxCell {
    std::array<MinimaxVector,4> points;
    std::array<bool,4> moving;
};

// Each dihedral/corner is a separate smooth constraint, including ties.
// Derivatives are with respect to a common translation of the moving nodes.
// Metric 2 maximizes the minimum angle; metric 6 maximizes scaled Jacobian.
inline bool minimaxConstraints(const MinimaxCell& cell, const MinimaxVector& x,
    int metric, bool derivatives, std::array<double,6>& values,
    std::array<MinimaxVector,6>& gradients, int& count) {
    using V = MinimaxVector;
    auto p=cell.points;
    for(int k=0;k<4;++k) if(cell.moving[k]) p[k]+=x;
    const double det=GEOM_FUNC::orient3d(p[0].data(),p[1].data(),p[3].data(),p[2].data());
    if(!(det>0) || !std::isfinite(det)) return false;
    const V a=p[1]-p[0], b=p[2]-p[0], c=p[3]-p[0];
    std::array<V,4> gd;
    gd[1]=b.cross(c);gd[2]=c.cross(a);gd[3]=a.cross(b);gd[0]=-gd[1]-gd[2]-gd[3];
    V volumeGradient=V::Zero();
    if(derivatives) for(int k=0;k<4;++k) if(cell.moving[k]) volumeGradient+=gd[k];
    if(metric==6) {
        count=4;
        for(int corner=0;corner<4;++corner) {
            double product=1;V logGradient=V::Zero();
            for(int k=0;k<4;++k) if(k!=corner) {
                const V edge=p[k]-p[corner];const double length2=edge.squaredNorm();
                if(!(length2>0) || !std::isfinite(length2)) return false;
                product*=std::sqrt(length2);
                if(derivatives) logGradient+=(int(cell.moving[k])-int(cell.moving[corner]))*edge/length2;
            }
            values[corner]=std::sqrt(2.0)*det/product;
            if(derivatives) gradients[corner]=values[corner]*(volumeGradient/det-logGradient);
            if(!(values[corner]>0) || !std::isfinite(values[corner]) ||
                (derivatives && !gradients[corner].allFinite())) return false;
        }
        return true;
    }
    if(metric!=2) return false;
    count=6;
    const int faces[4][3]={{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
    std::array<V,4> normals;
    std::array<MinimaxMatrix,4> normalGradient;
    for(int f=0;f<4;++f) {
        const int i=faces[f][0],j=faces[f][1],k=faces[f][2];
        const V u=p[j]-p[i],v=p[k]-p[i];
        normals[f]=u.cross(v);
        if(derivatives) for(int axis=0;axis<3;++axis) {
            const V e=V::Unit(axis);
            normalGradient[f].col(axis)=(int(cell.moving[j])-int(cell.moving[i]))*e.cross(v)
                +(int(cell.moving[k])-int(cell.moving[i]))*u.cross(e);
        }
    }
    int edge=0;
    for(int i=0;i<3;++i) for(int j=i+1;j<4;++j,++edge) {
        const int k=i>0?0:(j>1?1:2),l=6-i-j-k;
        const V d=p[j]-p[i];const double length=d.norm();
        if(!(length>0)) return false;
        const double y=det*length,z=-normals[k].dot(normals[l]);
        values[edge]=std::atan2(y,z);
        if(derivatives) {
            const V gy=length*volumeGradient+det*(int(cell.moving[j])-int(cell.moving[i]))*d/length;
            const V gz=-normalGradient[k].transpose()*normals[l]-normalGradient[l].transpose()*normals[k];
            const double denominator=y*y+z*z;
            if(!(denominator>0) || !std::isfinite(denominator)) return false;
            gradients[edge]=(z*gy-y*gz)/denominator;
        }
        if(!(values[edge]>0) || !std::isfinite(values[edge]) ||
            (derivatives && !gradients[edge].allFinite())) return false;
    }
    return true;
}
} }
