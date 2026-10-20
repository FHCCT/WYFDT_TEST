#include "dt.h"
#include "dt_optimization_policy.h"

namespace dt {
void addOptimizationQuality(DT& mesh, OptimizationQuality& quality,
    const double* a, const double* b, const double* c, const double* d) {
    double low = 0, high = 0; std::vector<double> angles;
    auto pa = const_cast<double*>(a), pb = const_cast<double*>(b);
    auto pc = const_cast<double*>(c), pd = const_cast<double*>(d);
    const double jacobian = mesh.ScaledJacobian(pa, pb, pc, pd);
    if (!mesh.CalDihedral(pa, pb, pd, pc, low, high, angles, true) ||
        !(low > 0) || !std::isfinite(low) || !(jacobian > 0) || !std::isfinite(jacobian)) {
        quality.valid = false; return;
    }
    low = RADIO2ANGLE(low);
    ++quality.count; quality.angleSum += low; quality.jacobianSum += jacobian;
    quality.minAngle = std::min(quality.minAngle, low);
    quality.minJacobian = std::min(quality.minJacobian, jacobian);
    quality.tinyAngle += low < .1; quality.tinyJacobian += jacobian < .001;
}

OptimizationQuality measureOptimizationQuality(DT& mesh, const std::vector<int>* cells) {
    OptimizationQuality quality;
    const auto add = [&](int t) {
        if (mesh.isDelEle(t) || mesh.ishulltet(t) || mesh.isvirtualtet(t)) return;
        const auto& v = mesh.Elems[t].form;
        addOptimizationQuality(mesh, quality, mesh.Nodes[v[0]].pt, mesh.Nodes[v[1]].pt,
            mesh.Nodes[v[2]].pt, mesh.Nodes[v[3]].pt);
    };
    if (cells) for (int t : *cells) add(t);
    else for (int t = 0; t < static_cast<int>(mesh.Elems.size()); ++t) add(t);
    quality.valid = quality.valid && quality.count > 0;
    return quality;
}

bool optimizationPointAccept(const OptimizationPolicy& policy,
    const OptimizationQuality& before, const OptimizationQuality& after) {
    return before.valid && after.valid && before.count == after.count &&
        after.meanAngle() + 1e-11 >= before.meanAngle() &&
        after.meanJacobian() + 1e-13 >= before.meanJacobian() &&
        after.minAngle + 1e-10 >= policy.angleFloor &&
        after.minJacobian + 1e-12 >= policy.jacobianFloor;
}

bool optimizationBoundaryAccept(const OptimizationPolicy& policy,
    const OptimizationQuality& before, const OptimizationQuality& after) {
    if (!before.valid || !after.valid || !before.count || !after.count ||
        after.minAngle + 1e-10 < std::min(before.minAngle, .1) ||
        after.minJacobian + 1e-12 < std::min(before.minJacobian, .001) ||
        after.tinyAngle > before.tinyAngle || after.tinyJacobian > before.tinyJacobian) return false;
    // Replacing a cavity can change cell count. Protect global means using
    // the count delta, rather than comparing two unrelated local averages.
    const double delta = double(after.count) - double(before.count);
    return double(after.angleSum - before.angleSum) + 1e-10 >= policy.meanAngle * delta &&
        double(after.jacobianSum - before.jacobianSum) + 1e-12 >= policy.meanJacobian * delta;
}

double angleJacobianQuality(double angle, double jacobian) {
    return angle > 0 && jacobian > 0 && std::isfinite(angle) && std::isfinite(jacobian)
        ? std::pow(angle, .25) * std::pow(jacobian, .75) : 0;
}

}
