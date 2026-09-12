#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include "dt.h"

namespace {
    constexpr double kPi = 3.141592653589793238462643383279502884;
    constexpr double kEps = 1e-16;
    constexpr double kSqrt2 = 1.4142135623730951;

    inline double clamp(double v) {
        if (v > 1.0) {
            return 1.0;
        }
        if (v < -1.0) {
            return -1.0;
        }
        return v;
    }

    inline double vecCos(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
        const double len1 = a.norm();
        const double len2 = b.norm();
        if (len1 <= 0.0 || len2 <= 0.0 || !std::isfinite(len1) ||
            !std::isfinite(len2)) {
            return 0.0;
        }
        return clamp(a.dot(b) / (len1 * len2));
    }

    inline void triMinMaxAngle(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1,
        const Eigen::Vector3d& p2, double& minDeg,
        double& maxDeg) {
        const Eigen::Vector3d s0 = p1 - p0;
        const Eigen::Vector3d s1 = p2 - p1;
        const Eigen::Vector3d s2 = p2 - p0;
        const Eigen::Vector3d s3 = -s1;
        const double sq0 = s0.squaredNorm();
        const double sq1 = s1.squaredNorm();
        const double sq2 = s2.squaredNorm();
        if (sq0 == 0.0 || sq1 == 0.0 || sq2 == 0.0) {
            minDeg = 0.0;
            maxDeg = 0.0;
            return;
        }

        int shortSide = 0;
        if (sq1 < sq0) {
            shortSide = 1;
        }
        if (sq2 < (shortSide == 0 ? sq0 : sq1)) {
            shortSide = 2;
        }

        int longSide = 0;
        if (sq1 > sq0) {
            longSide = 1;
        }
        if (sq2 > (longSide == 0 ? sq0 : sq1)) {
            longSide = 2;
        }

        if (shortSide == 0) {
            minDeg = std::acos(vecCos(s2, s1)) * (180.0 / kPi);
        }
        else if (shortSide == 1) {
            minDeg = std::acos(vecCos(s0, s2)) * (180.0 / kPi);
        }
        else {
            minDeg = std::acos(vecCos(s0, s3)) * (180.0 / kPi);
        }

        if (longSide == 0) {
            maxDeg = std::acos(vecCos(s2, s1)) * (180.0 / kPi);
        }
        else if (longSide == 1) {
            maxDeg = std::acos(vecCos(s0, s2)) * (180.0 / kPi);
        }
        else {
            maxDeg = std::acos(vecCos(s0, s3)) * (180.0 / kPi);
        }
    }

    inline double triVolumeSkew(const Eigen::Vector3d& p0,
        const Eigen::Vector3d& p1,
        const Eigen::Vector3d& p2) {
        const double area = 0.5 * (p1 - p0).cross(p2 - p0).norm();
        if (area <= kEps || !std::isfinite(area)) {
            return 1.0;
        }

        const double aLen = (p1 - p2).norm();
        const double bLen = (p2 - p0).norm();
        const double cLen = (p0 - p1).norm();
        const double abc = aLen * bLen * cLen;
        if (!std::isfinite(abc)) {
            return 1.0;
        }

        const double radius = abc / (4.0 * area);
        const double optArea = (3.0 * std::sqrt(3.0) / 4.0) * radius * radius;
        if (optArea <= kEps || !std::isfinite(optArea)) {
            return 1.0;
        }

        const double skew = (optArea - area) / optArea;
        return std::isfinite(skew) ? skew : 1.0;
    }

    inline double tetVSkew(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
        const Eigen::Vector3d& c, const Eigen::Vector3d& d) {
        const Eigen::Vector3d ab = b - a;
        const Eigen::Vector3d ac = c - a;
        const Eigen::Vector3d ad = d - a;
        const double lengthAB2 = ab.squaredNorm();
        const double lengthAC2 = ac.squaredNorm();
        const double lengthAD2 = ad.squaredNorm();
        const Eigen::Vector3d cpBC = ab.cross(ac);
        const Eigen::Vector3d cpDB = ad.cross(ab);
        const Eigen::Vector3d cpCD = ac.cross(ad);
        const Eigen::Vector3d num =
            lengthAD2 * cpBC + lengthAC2 * cpDB + lengthAB2 * cpCD;
        const double den = 2.0 * ab.dot(cpCD);
        if (std::abs(den) < kEps || !std::isfinite(den)) {
            return 1.0;
        }

        const double circumradius = num.norm() / std::abs(den);
        const double optLength = circumradius / std::sqrt(3.0 / 8.0);
        const double optVolume =
            (1.0 / 12.0) * kSqrt2 * optLength * optLength * optLength;
        if (std::abs(optVolume) < kEps || !std::isfinite(optVolume)) {
            return 1.0;
        }

        const double volume = std::abs((d - a).dot((b - a).cross(c - a))) / 6.0;
        const double skew = (optVolume - volume) / optVolume;
        return std::isfinite(skew) ? skew : 1.0;
    }

    inline std::pair<double, double>
        tetraDihedralMinMax(const Eigen::Vector3d v[4]) {
        static constexpr int faces[4][3] = {
            {0, 2, 1},
            {0, 1, 3},
            {1, 2, 3},
            {0, 3, 2},
        };
        static constexpr int edgeFacePairs[6][2] = {
            {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3},
        };

        std::array<Eigen::Vector3d, 4> normals;
        for (int fi = 0; fi < 4; ++fi) {
            const int* face = faces[fi];
            const Eigen::Vector3d area =
                (v[face[1]] - v[face[0]]).cross(v[face[2]] - v[face[0]]);
            const double len = area.norm();
            if (len <= kEps || !std::isfinite(len)) {
                normals[fi] = Eigen::Vector3d::Zero();
            }
            else {
                normals[fi] = area / len;
            }
        }

        double minDih = std::numeric_limits<double>::infinity();
        double maxDih = -std::numeric_limits<double>::infinity();
        for (const auto& pair : edgeFacePairs) {
            const Eigen::Vector3d& n0 = normals[pair[0]];
            const Eigen::Vector3d& n1 = normals[pair[1]];
            if (n0.squaredNorm() < 1e-30 || n1.squaredNorm() < 1e-30) {
                continue;
            }
            const double angle = std::acos(clamp(-n0.dot(n1))) * (180.0 / kPi);
            minDih = std::min(minDih, angle);
            maxDih = std::max(maxDih, angle);
        }

        if (!std::isfinite(minDih) || !std::isfinite(maxDih)) {
            return { 0.0, 180.0 };
        }
        return { minDih, maxDih };
    }

    inline double faceOrthogonal(const Eigen::Vector3d& area,
        const Eigen::Vector3d& faceCenter,
        const Eigen::Vector3d& cellCenter,
        const double* oppositeCenter) {
        // fVec：当前单元质心 -> 当前面心
        const Eigen::Vector3d faceVector = faceCenter - cellCenter;
        const double areaNorm = area.norm();
        const double faceNorm = faceVector.norm();
        if (areaNorm < kEps || faceNorm < kEps || !std::isfinite(areaNorm) ||
            !std::isfinite(faceNorm) || !std::isfinite(faceVector[0]) ||
            !std::isfinite(faceVector[1]) || !std::isfinite(faceVector[2])) {
            return 0.0;
        }

        double orth = vecCos(area, faceVector);
        if (oppositeCenter != nullptr) {
            const Eigen::Vector3d neighborCenter(oppositeCenter[0], oppositeCenter[1],
                oppositeCenter[2]);
            // cVec：当前单元质心 -> 邻单元质心
            const Eigen::Vector3d cellVector = neighborCenter - cellCenter;
            const double cellNorm = cellVector.norm();
            if (cellNorm < kEps || !std::isfinite(cellNorm) ||
                !std::isfinite(cellVector[0]) || !std::isfinite(cellVector[1]) ||
                !std::isfinite(cellVector[2])) {
                return 0.0;
            }
            orth = std::min(orth, vecCos(area, cellVector));
        }
        return orth;
    }

} // namespace

double DT::calskewness(double* pa, double* pb, double* pc, double* pd) {
    const Eigen::Vector3d v[4] = {
        Eigen::Vector3d(pa[0], pa[1], pa[2]),
        Eigen::Vector3d(pb[0], pb[1], pb[2]),
        Eigen::Vector3d(pc[0], pc[1], pc[2]),
        Eigen::Vector3d(pd[0], pd[1], pd[2]),
    };

    static constexpr int faces[4][3] = {
        {0, 2, 1},
        {0, 1, 3},
        {1, 2, 3},
        {0, 3, 2},
    };

    double equiangleSkew = 0.0;
    double triVSkewMax = -std::numeric_limits<double>::infinity();

    for (const auto& face : faces) {
        double mn = 0.0;
        double mx = 0.0;
        triMinMaxAngle(v[face[0]], v[face[1]], v[face[2]], mn, mx);
        equiangleSkew = std::max(equiangleSkew,
            std::max((mx - 60.0) / 120.0, (60.0 - mn) / 60.0));
        triVSkewMax = std::max(triVSkewMax,
            triVolumeSkew(v[face[0]], v[face[1]], v[face[2]]));
    }

    const auto dihedral = tetraDihedralMinMax(v);
    if (std::isfinite(dihedral.first) && std::isfinite(dihedral.second)) {
        const double thetaDeg = std::acos(1.0 / 3.0) * (180.0 / kPi);
        const double dihedralSkewMax =
            (dihedral.second - thetaDeg) / (180.0 - thetaDeg);
        const double dihedralSkewMin = (thetaDeg - dihedral.first) / thetaDeg;
        equiangleSkew =
            std::max(equiangleSkew, std::max(dihedralSkewMin, dihedralSkewMax));
    }

    const double equivolumeSkew =
        std::max(triVSkewMax, tetVSkew(v[0], v[1], v[2], v[3]));
    return std::max(equiangleSkew, equivolumeSkew);
}

double DT::orthogonal(double* pa, double* pb, double* pc, double* pd,
    double* oppositeA, double* oppositeB, double* oppositeC,
    double* oppositeD) {
    const Eigen::Vector3d v[4] = {
        Eigen::Vector3d(pa[0], pa[1], pa[2]),
        Eigen::Vector3d(pb[0], pb[1], pb[2]),
        Eigen::Vector3d(pc[0], pc[1], pc[2]),
        Eigen::Vector3d(pd[0], pd[1], pd[2]),
    };
    const Eigen::Vector3d cellCenter = 0.25 * (v[0] + v[1] + v[2] + v[3]);
    if (!std::isfinite(cellCenter[0]) || !std::isfinite(cellCenter[1]) ||
        !std::isfinite(cellCenter[2])) {
        return 0.0;
    }

    static constexpr int faceIdx[4][3] = {
        {1, 2, 3},
        {0, 3, 2},
        {0, 1, 3},
        {0, 2, 1},
    };
    const double* oppositeCenters[4] = { oppositeA, oppositeB, oppositeC,
                                        oppositeD };

    double orth = 1.0;
    for (int fi = 0; fi < 4; ++fi) {
        const int* face = faceIdx[fi];
        // aVec：当前面的面法向
        const Eigen::Vector3d area =
            (v[face[1]] - v[face[0]]).cross(v[face[2]] - v[face[0]]);
        const Eigen::Vector3d faceCenter =
            (v[face[0]] + v[face[1]] + v[face[2]]) / 3.0;
        orth = std::min(orth, faceOrthogonal(area, faceCenter, cellCenter,
            oppositeCenters[fi]));
    }

    const double skew = calskewness(pa, pb, pc, pd);
    if (!std::isfinite(skew)) {
        return 0.0;
    }
    const double result = std::min(orth, 1.0 - skew);
    return std::isfinite(result) ? result : 0.0;
}
