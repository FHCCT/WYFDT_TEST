#pragma once
#include <cstddef>
#include <vector>

namespace dt {
class DT;

// State belongs to one DT run. Workers read it; only serial stage boundaries
// change it. No process-wide switches or environment variables are involved.
struct OptimizationPolicy {
    bool enabled = false, inTopology = false;
    double angleFloor = 0, jacobianFloor = 0;
    double meanAngle = 0, meanJacobian = 0;
};

struct OptimizationQuality {
    size_t count = 0, tinyAngle = 0, tinyJacobian = 0;
    long double angleSum = 0, jacobianSum = 0;
    double minAngle = 180, minJacobian = 1;
    bool valid = true;
    double meanAngle() const { return count ? double(angleSum / count) : 0; }
    double meanJacobian() const { return count ? double(jacobianSum / count) : 0; }
};

void addOptimizationQuality(DT& mesh, OptimizationQuality& quality,
    const double* a, const double* b, const double* c, const double* d);
OptimizationQuality measureOptimizationQuality(DT& mesh, const std::vector<int>* cells = nullptr);
bool optimizationPointAccept(const OptimizationPolicy& policy,
    const OptimizationQuality& before, const OptimizationQuality& after);
bool optimizationBoundaryAccept(const OptimizationPolicy& policy,
    const OptimizationQuality& before, const OptimizationQuality& after);
double angleJacobianQuality(double normalizedAngle, double jacobian);
}

namespace dt {
namespace detail {

struct QualityTotals {
    size_t count = 0;
    long double angle = 0, jacobian = 0, centeredSquare = 0, aspect = 0;
    long double lowerAngle = 0, upperAspect = 0, focus = 0;
    double minAngle = 180, minJacobian = 1;
    bool valid = true;
};

// One serial phase owns this ledger. Every proposal measures its actual old
// and new cells before writing the host, then commits only the local delta.
class QualityBudget {
public:
    QualityBudget(DT& mesh, double targetAngle);
    QualityTotals cells(DT& mesh, const std::vector<int>& cells) const;
    void add(DT& mesh, QualityTotals& total, const double* a, const double* b,
        const double* c, const double* d) const;
    bool accepts(const QualityTotals& before, const QualityTotals& after) const;
    void commit(const QualityTotals& before, const QualityTotals& after);
    double target() const { return targetAngle_; }
    bool valid() const { return initial_.valid && initial_.count > 0; }
private:
    QualityTotals initial_, current_;
    double targetAngle_, angleCutoff_ = 0, aspectCutoff_ = 0, jacobianCenter_ = 0;
    QualityTotals replacement(const QualityTotals& before, const QualityTotals& after) const;
};
}
}
