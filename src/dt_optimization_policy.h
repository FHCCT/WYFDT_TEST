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
