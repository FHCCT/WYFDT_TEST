#pragma once
#include "dt.h"
#include <limits>

namespace dt {

// Measure the exported geometric qualities, never the SUS objective or a cache.
struct SmoothingQuality {
    size_t count = 0;
    long double shapeSum = 0, angleSum = 0;
    double minShape = DBL_MAX, minAngle = DBL_MAX;
    bool valid = true;
    double meanShape() const { return count ? static_cast<double>(shapeSum / count) : 0; }
    double meanAngle() const { return count ? static_cast<double>(angleSum / count) : 0; }
};

inline SmoothingQuality measureSmoothingQuality(DT& mesh, const std::vector<int>& cells) {
    SmoothingQuality result;
    std::vector<double> angles;
    for (int t : cells) {
        if (mesh.isDelEle(t) || mesh.ishulltet(t) || mesh.isvirtualtet(t)) continue;
        const int* v = mesh.Elems[t].form;
        double *a = mesh.Nodes[v[0]].pt, *b = mesh.Nodes[v[1]].pt,
               *c = mesh.Nodes[v[2]].pt, *d = mesh.Nodes[v[3]].pt;
        const double shape = mesh.quality_sus(a, b, c, d);
        double lo, hi;
        if (!(shape > 0) || !std::isfinite(shape) ||
            !mesh.CalDihedral(a, b, d, c, lo, hi, angles) || !(lo > 0) || !std::isfinite(lo)) {
            result.valid = false;
            return result;
        }
        ++result.count;
        result.shapeSum += shape;
        result.angleSum += lo;
        result.minShape = std::min(result.minShape, shape);
        result.minAngle = std::min(result.minAngle, lo);
    }
    result.valid = result.count != 0;
    return result;
}

inline double smoothingTolerance(double value) {
    return 64 * std::numeric_limits<double>::epsilon() * std::max(1e-12, std::abs(value));
}

inline bool smoothingPreservesQuality(const SmoothingQuality& before, const SmoothingQuality& after) {
    return before.valid && after.valid && before.count == after.count &&
        after.meanShape() >= before.meanShape() - smoothingTolerance(before.meanShape()) &&
        after.meanAngle() >= before.meanAngle() - smoothingTolerance(before.meanAngle()) &&
        after.minShape >= before.minShape - smoothingTolerance(before.minShape) &&
        after.minAngle >= before.minAngle - smoothingTolerance(before.minAngle);
}

inline bool smoothingImprovesQuality(const SmoothingQuality& before, const SmoothingQuality& after) {
    return smoothingPreservesQuality(before, after) &&
        (after.meanShape() > before.meanShape() + smoothingTolerance(before.meanShape()) ||
         after.meanAngle() > before.meanAngle() + smoothingTolerance(before.meanAngle()));
}
}
