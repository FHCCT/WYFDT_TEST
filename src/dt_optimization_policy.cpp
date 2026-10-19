#include "dt.h"
#include "dt_optimization_policy.h"

namespace dt {
namespace {
double regionCap(const OptimizationPolicy& policy, int region) {
    const auto it = policy.regionSize.find(region);
    return it != policy.regionSize.end() && it->second > 0
        ? std::min(it->second, policy.globalMaxEdge) : policy.globalMaxEdge;
}
}

void OptimizationPolicy::captureSizes(DT& mesh) {
    pointIds.clear(); legacyLongEdges.clear();
    for (int n = 0; n < static_cast<int>(mesh.Nodes.size()); ++n)
        pointIds[mesh.Nodes[n].pt] = n;
    for (int t = 0; t < static_cast<int>(mesh.Elems.size()); ++t) {
        if (mesh.isDelEle(t) || mesh.ishulltet(t) || mesh.isvirtualtet(t)) continue;
        const auto& cell = mesh.Elems[t];
        const double cap = regionCap(*this, cell.geo);
        for (int i = 0; i < 4; ++i) for (int j = 0; j < i; ++j) {
            const int a = std::min(cell.form[i], cell.form[j]), b = std::max(cell.form[i], cell.form[j]);
            const double length = mesh.distance(mesh.Nodes[a].pt, mesh.Nodes[b].pt);
            if (length > cap * (1 + 1e-10)) legacyLongEdges[{{cell.geo, a, b}}] = length;
        }
    }
}

bool OptimizationPolicy::edgeAllowed(int region, int a, int b, const double* pa, const double* pb) const {
    if (!enabled) return true;
    double cap = regionCap(*this, region);
    if (a > b) std::swap(a, b);
    const auto old = legacyLongEdges.find({{region, a, b}});
    if (old != legacyLongEdges.end()) cap = std::max(cap, old->second);
    double squared = 0;
    for (int k = 0; k < 3; ++k) squared += (pa[k] - pb[k]) * (pa[k] - pb[k]);
    return std::isfinite(squared) && squared <= cap * cap * (1 + 2e-10);
}

bool OptimizationPolicy::moveAllowed(DT& mesh, int node, const double* point) const {
    if (!enabled) return true;
    std::vector<int> star; mesh.findSphere(node, star);
    for (int t : star) {
        if (mesh.isDelEle(t) || mesh.ishulltet(t) || mesh.isvirtualtet(t)) continue;
        for (int n : mesh.Elems[t].form) if (n != node &&
            !edgeAllowed(mesh.Elems[t].geo, node, n, point, mesh.Nodes[n].pt)) return false;
    }
    return true;
}

bool OptimizationPolicy::tetAllowed(const double* a, const double* b, const double* c, const double* d) const {
    const double* p[4] = {a, b, c, d}; int ids[4];
    for (int i = 0; i < 4; ++i) {
        const auto it = pointIds.find(p[i]);
        if (it == pointIds.end()) return false;
        ids[i] = it->second;
    }
    for (int i = 0; i < 4; ++i) for (int j = 0; j < i; ++j)
        if (!edgeAllowed(activeRegion, ids[i], ids[j], p[i], p[j])) return false;
    return true;
}

void addOptimizationQuality(DT& mesh, OptimizationQuality& quality,
    const double* a, const double* b, const double* c, const double* d) {
    double low = 0, high = 0; std::vector<double> angles;
    auto pa = const_cast<double*>(a), pb = const_cast<double*>(b);
    auto pc = const_cast<double*>(c), pd = const_cast<double*>(d);
    const double jacobian = mesh.ScaledJacobian(pa, pb, pc, pd);
    if (!mesh.CalDihedral(pa, pb, pd, pc, low, high, angles) ||
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

int DT::qualityOptimizationPass(int rounds, const std::unordered_map<int, double>& regionSize,
    double globalMaxEdge, bool smooth) {
    rounds = std::max(2, rounds);
    // Scope the guards to this optimization phase. Refinement and local trial
    // meshes retain their own acceptance policies and cannot inherit this state.
    struct RestorePolicy {
        DT& mesh; OptimizationPolicy saved; int depth, metric; bool improving;
        explicit RestorePolicy(DT& d) : mesh(d), saved(std::move(d.optimizationPolicy)),
            depth(d.optflipdeep), metric(d.improve_Metric), improving(d.improve_step) {}
        ~RestorePolicy() {
            mesh.optimizationPolicy = std::move(saved);
            mesh.optflipdeep = depth; mesh.improve_Metric = metric; mesh.improve_step = improving;
        }
    } restore(*this);
    optimizationPolicy = OptimizationPolicy{};
    auto& policy = optimizationPolicy;
    policy.enabled = true; policy.regionSize = regionSize;
    policy.globalMaxEdge = std::min(globalMaxEdge, maxEdge > 0 ? maxEdge : 1e100);
    policy.captureSizes(*this);
    improve_step = true; optflipdeep = 0;
    const int susRounds = (rounds + 1) / 2;
    for (int loop = 0; loop < rounds; ++loop) {
        if (infolevel > 0) meshLogger->info("Quality loop {}/{}: topology={}",
            loop + 1, rounds, loop < susRounds ? "SUS" : "A^0.25 J^0.75");
        updateminVolume();
        auto quality = measureOptimizationQuality(*this);
        if (!quality.valid) return 0;
        policy.meanAngle = quality.meanAngle(); policy.meanJacobian = quality.meanJacobian();
        flipEdgPass(2);

        updateminVolume();
        quality = measureOptimizationQuality(*this);
        policy.angleFloor = quality.minAngle; policy.jacobianFloor = quality.minJacobian;
        improve_Metric = loop < susRounds ? SUS_METRIC : ANGLE_JACOBIAN_METRIC;
        policy.inTopology = false; prepareQuality();
        policy.inTopology = true;
        TopologicalPass(1, 0, 1); // Pure flips: no refinement or topology insertion.
        policy.inTopology = false;

        if (smooth) {
            updateminVolume();
            improve_Metric = SUS_METRIC; prepareQuality();
            SmoothPass(1, 1);
        }
        improve_Metric = SUS_METRIC; prepareQuality();
    }
    // Restore the caller's metric with matching cached values.
    improve_Metric = restore.metric; prepareQuality();
    return 1;
}
}
