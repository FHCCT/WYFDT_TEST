#include "dt.h"
#include "dt_optimization_policy.h"
#include <limits>
#include <stdexcept>

namespace dt {
void addOptimizationQuality(DT& mesh, OptimizationQuality& quality,
    const double* a, const double* b, const double* c, const double* d) {
    double low = 0, high = 0; std::vector<double> angles;
    auto pa = const_cast<double*>(a), pb = const_cast<double*>(b);
    auto pc = const_cast<double*>(c), pd = const_cast<double*>(d);
    double jacobian = 0;
    if (!mesh.CalDihedral(pa, pb, pd, pc, low, high, angles, true, &jacobian) ||
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

namespace dt { namespace detail {
namespace {
long double tolerance(long double value) {
    return 128 * std::numeric_limits<double>::epsilon() * std::max(1e-12L, std::abs(value));
}
size_t tailCount(size_t count) { return count / 10 + (count % 10 != 0); }
long double variance(const QualityTotals& q, double center) {
    const long double shift = q.jacobian / q.count - center;
    return q.centeredSquare / q.count - shift * shift;
}
long double lowerBound(const QualityTotals& q, double cutoff) {
    const size_t k = tailCount(q.count);
    return (q.lowerAngle - (q.count - k) * static_cast<long double>(cutoff)) / k;
}
long double upperBound(const QualityTotals& q, double cutoff) {
    const size_t k = tailCount(q.count);
    return (q.upperAspect - (q.count - k) * static_cast<long double>(cutoff)) / k;
}
}

QualityBudget::QualityBudget(DT& mesh, double targetAngle) : targetAngle_(targetAngle) {
    if (!(targetAngle > 0 && targetAngle < 180) || !std::isfinite(targetAngle))
        throw std::invalid_argument("Quality tail target must be in (0,180)");
    std::vector<int> physical;
    std::vector<double> angles, aspects;
    for (int t = 0; t < static_cast<int>(mesh.Elems.size()); ++t) {
        if (mesh.isDelEle(t) || mesh.ishulltet(t) || mesh.isvirtualtet(t)) continue;
        physical.push_back(t);
        QualityTotals q;
        const auto& v = mesh.Elems[t].form;
        add(mesh, q, mesh.Nodes[v[0]].pt, mesh.Nodes[v[1]].pt, mesh.Nodes[v[2]].pt, mesh.Nodes[v[3]].pt);
        if (!q.valid) { initial_.valid = false; return; }
        angles.push_back(static_cast<double>(q.angle));
        aspects.push_back(static_cast<double>(q.aspect));
        initial_.jacobian += q.jacobian;
    }
    if (physical.empty()) return;
    jacobianCenter_ = static_cast<double>(initial_.jacobian / physical.size());
    const size_t k = tailCount(physical.size());
    std::nth_element(angles.begin(), angles.begin() + k - 1, angles.end());
    angleCutoff_ = angles[k - 1];
    std::nth_element(aspects.begin(), aspects.begin() + aspects.size() - k, aspects.end());
    aspectCutoff_ = aspects[aspects.size() - k];
    initial_ = cells(mesh, physical);
    current_ = initial_;
}

void QualityBudget::add(DT& mesh, QualityTotals& q, const double* a, const double* b,
    const double* c, const double* d) const {
    auto pa = const_cast<double*>(a), pb = const_cast<double*>(b);
    auto pc = const_cast<double*>(c), pd = const_cast<double*>(d);
    double lo = 0, hi = 0, sj = 0, aspect = 0;
    std::vector<double> angles;
    if (!mesh.CalDihedral(pa, pb, pd, pc, lo, hi, angles, true, &sj, &aspect)) { q.valid = false; return; }
    lo = RADIO2ANGLE(lo);
    if (!(lo > 0 && sj > 0) || !std::isfinite(lo) ||
        !std::isfinite(sj) || !std::isfinite(aspect)) { q.valid = false; return; }
    ++q.count;
    q.angle += lo; q.jacobian += sj; q.aspect += aspect;
    const long double centered = static_cast<long double>(sj) - jacobianCenter_;
    q.centeredSquare += centered * centered;
    q.lowerAngle += std::min(lo, angleCutoff_);
    q.upperAspect += std::max(aspect, aspectCutoff_);
    q.focus += std::min(lo, targetAngle_);
    q.minAngle = std::min(q.minAngle, lo); q.minJacobian = std::min(q.minJacobian, sj);
}

QualityTotals QualityBudget::cells(DT& mesh, const std::vector<int>& list) const {
    QualityTotals q;
    for (int t : list) {
        if (t < 0 || t >= static_cast<int>(mesh.Elems.size()) || mesh.isDelEle(t) ||
            mesh.ishulltet(t) || mesh.isvirtualtet(t)) continue;
        const auto& v = mesh.Elems[t].form;
        add(mesh, q, mesh.Nodes[v[0]].pt, mesh.Nodes[v[1]].pt, mesh.Nodes[v[2]].pt, mesh.Nodes[v[3]].pt);
    }
    return q;
}

QualityTotals QualityBudget::replacement(const QualityTotals& a, const QualityTotals& b) const {
    QualityTotals next = current_;
    if (a.count > next.count) { next.valid = false; return next; }
    next.count = next.count - a.count + b.count;
    next.angle += b.angle - a.angle; next.jacobian += b.jacobian - a.jacobian;
    next.centeredSquare += b.centeredSquare - a.centeredSquare;
    next.aspect += b.aspect - a.aspect;
    next.lowerAngle += b.lowerAngle - a.lowerAngle;
    next.upperAspect += b.upperAspect - a.upperAspect;
    next.focus += b.focus - a.focus;
    return next;
}

bool QualityBudget::accepts(const QualityTotals& a, const QualityTotals& b) const {
    if (!valid() || !a.valid || !b.valid || !a.count || !b.count) return false;
    const auto next = replacement(a, b);
    if (!next.valid || !next.count) return false;
    const long double n = next.count, refN = initial_.count;
    const long double angleMean = initial_.angle / refN, sjMean = initial_.jacobian / refN;
    const long double aspectMean = initial_.aspect / refN;
    const long double initialVariance = variance(initial_, jacobianCenter_);
    // Fixed cutoffs give rigorous bounds on the worst ceil(0.1N) means even
    // when a topology proposal changes N. This requires only local sums.
    return b.minAngle + tolerance(initial_.minAngle) >= initial_.minAngle &&
        b.minJacobian + tolerance(initial_.minJacobian) >= initial_.minJacobian &&
        next.angle / n + tolerance(angleMean) >= angleMean &&
        next.jacobian / n + tolerance(sjMean) >= sjMean &&
        variance(next, jacobianCenter_) <= initialVariance + tolerance(initialVariance) &&
        next.aspect / n <= aspectMean + tolerance(aspectMean) &&
        lowerBound(next, angleCutoff_) + tolerance(lowerBound(initial_, angleCutoff_)) >= lowerBound(initial_, angleCutoff_) &&
        upperBound(next, aspectCutoff_) <= upperBound(initial_, aspectCutoff_) + tolerance(upperBound(initial_, aspectCutoff_)) &&
        next.focus / n > current_.focus / current_.count + 1e-10L;
}

void QualityBudget::commit(const QualityTotals& a, const QualityTotals& b) {
    current_ = replacement(a, b);
}
} }
