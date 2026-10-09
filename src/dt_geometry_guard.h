#pragma once
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>

namespace dt {
// Allow coordinate roundoff, not a fraction of edge length. Translate and
// scale first so thin, translated models use the same plane membership rule.
inline bool pointOnBoundaryPlane(const double* a, const double* b,
    const double* c, const double* p) {
    using Vec = Eigen::Vector3d;
    const Vec origin(a), ab = Vec(b) - origin, ac = Vec(c) - origin;
    const double scale = std::max(ab.cwiseAbs().maxCoeff(), ac.cwiseAbs().maxCoeff());
    if (!(scale > 0) || !std::isfinite(scale)) return false;
    const Vec normal = (ab / scale).cross(ac / scale);
    const double length = normal.norm();
    if (!(length > 0) || !std::isfinite(length)) return false;
    double magnitude = 0;
    for (int k = 0; k < 3; ++k)
        magnitude = std::max(magnitude, std::max(std::abs(a[k]),
            std::max(std::abs(b[k]), std::max(std::abs(c[k]), std::abs(p[k])))));
    const double tolerance = 32 * std::numeric_limits<double>::epsilon() * magnitude + 1e-12 * scale;
    const double offset = std::abs(normal.dot((Vec(p) - origin) / scale)) / length;
    return std::isfinite(offset) && offset <= tolerance / scale;
}

// Absolute volume alone cannot exclude slivers. Use dimensionless mean-ratio
// quality, and protect the poor tail of the entire old optimization cavity.
inline bool acceptableInsertionQuality(double oldMinimum, double newMinimum, bool optimizing) {
    const double numericalFloor = 1e-8;
    const double cutoff = optimizing ? std::max(numericalFloor, std::min(oldMinimum, 0.01)) : numericalFloor;
    return std::isfinite(newMinimum) && newMinimum >= cutoff;
}
}
