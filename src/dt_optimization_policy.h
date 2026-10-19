#pragma once
#include <array>
#include <cstddef>
#include <map>
#include <unordered_map>
#include <vector>

namespace dt {
class DT;

// State belongs to one DT run. Workers read it; only serial stage boundaries
// change it. No process-wide switches or environment variables are involved.
struct OptimizationPolicy {
    bool enabled = false, inTopology = false;
    double angleFloor = 0, jacobianFloor = 0;
    double meanAngle = 0, meanJacobian = 0;
    double globalMaxEdge = 1e100;
    int activeRegion = -1;
    std::unordered_map<int, double> regionSize;
    std::unordered_map<const double*, int> pointIds;
    std::map<std::array<int, 3>, double> legacyLongEdges;

    void captureSizes(DT& mesh);
    bool edgeAllowed(int region, int a, int b, const double* pa, const double* pb) const;
    bool moveAllowed(DT& mesh, int node, const double* point) const;
    bool tetAllowed(const double* a, const double* b, const double* c, const double* d) const;
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
