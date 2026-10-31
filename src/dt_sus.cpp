#include "dt.h"
#include "dt_minimax_smoothing.h"
#include "dt_geometry_guard.h"
#include <stdexcept>

// Escobar et al., CMAME 192 (2003), 2775-2787; SUS code (2010),
// Calcula_Delta / K_Knupp_Modif_directa and TechNote.pdf.
// Independent implementation with local normalization and BFGS/Armijo.
namespace {
using Vec = Eigen::Vector3d;
using Mat = Eigen::Matrix3d;
struct SusTet {
    std::array<Vec, 4> points;
    int freeIndex;
    bool positive;
};

double sigma(std::array<Vec, 4> p) {
    return std::sqrt(2.0) * static_cast<double>(dt::GEOM_FUNC::orient3d(
        p[0].data(), p[1].data(), p[3].data(), p[2].data()));
}

double normSquared(const std::array<Vec, 4>& p) {
    double result = 0;
    for (int a = 0; a < 3; ++a)
        for (int b = a + 1; b < 4; ++b)
            result += 0.5 * (p[a] - p[b]).squaredNorm();
    return result;
}

double signedQuality(const std::array<Vec, 4>& p) {
    const double norm = normSquared(p), det = sigma(p);
    if (!(norm > 0) || !std::isfinite(det)) return 0;
    return std::copysign(std::min(1.0, 3 * std::pow(std::fabs(det), 2.0 / 3.0) / norm), det);
}

Vec qualityGradient(const SusTet& tet, const Vec& x) {
    auto p = tet.points;
    p[tet.freeIndex] = x;
    const double det = sigma(p), norm = normSquared(p);
    if (!(norm > 0)) return Vec::Zero();
    // Signed mean ratio is not differentiable at det=0. Use a symmetric
    // finite-difference direction there; actual minQ still controls acceptance.
    if (std::fabs(det) < 1e-12) {
        Vec g;
        for (int j = 0; j < 3; ++j) {
            p[tet.freeIndex][j] = x[j] + 1e-6;
            const double plus = signedQuality(p);
            p[tet.freeIndex][j] = x[j] - 1e-6;
            const double minus = signedQuality(p);
            p[tet.freeIndex][j] = x[j];
            g[j] = (plus - minus) / 2e-6;
        }
        return g;
    }
    const Vec a = p[1] - p[0], b = p[2] - p[0], c = p[3] - p[0];
    std::array<Vec, 4> gd;
    gd[1] = b.cross(c); gd[2] = c.cross(a); gd[3] = a.cross(b);
    gd[0] = -gd[1] - gd[2] - gd[3];
    Vec gn = Vec::Zero();
    for (int j = 0; j < 4; ++j)
        if (j != tet.freeIndex) gn += x - p[j];
    return 2 * std::sqrt(2.0) * gd[tet.freeIndex] / (std::cbrt(std::fabs(det)) * norm)
        - signedQuality(p) * gn / norm;
}

// For at most four gradients, inspect simplex faces directly. Degenerate
// faces are covered by their edges/vertices; larger hulls use Frank-Wolfe.
Vec smallHullNearest(const std::vector<Vec>& g) {
    Vec best = g.front();
    auto consider = [&](const Vec& p) { if (p.squaredNorm() < best.squaredNorm()) best = p; };
    const int n = static_cast<int>(g.size());
    for (int i = 0; i < n; ++i) {
        consider(g[i]);
        for (int j = i + 1; j < n; ++j) {
            const Vec u = g[j] - g[i];
            const double uu = u.squaredNorm();
            if (uu > 0) consider(g[i] + std::max(0.0, std::min(1.0, -g[i].dot(u) / uu)) * u);
            for (int k = j + 1; k < n; ++k) {
                const Vec v = g[k] - g[i];
                const double uv = u.dot(v), vv = v.squaredNorm();
                const double det = u.cross(v).squaredNorm();
                if (det <= 1e-14 * uu * vv) continue;
                const double rhsU = -g[i].dot(u), rhsV = -g[i].dot(v);
                const double b = (rhsU * vv - rhsV * uv) / det;
                const double c = (rhsV * uu - rhsU * uv) / det;
                if (b >= 0 && c >= 0 && b + c <= 1) consider(g[i] + b * u + c * v);
            }
        }
    }
    if (n == 4) {
        const Vec u=g[1]-g[0], v=g[2]-g[0], w=g[3]-g[0], rhs=-g[0];
        const double det=u.dot(v.cross(w));
        if (std::fabs(det) > 1e-12 * u.norm() * v.norm() * w.norm()) {
            const double b=rhs.dot(v.cross(w))/det;
            const double c=u.dot(rhs.cross(w))/det;
            const double d=u.dot(v.cross(rhs))/det;
            if (b >= 0 && c >= 0 && d >= 0 && b+c+d <= 1) return Vec::Zero();
        }
    }
    return best;
}

// Shared convex-hull solve for SUS and single-metric angle/SJ constraints.
Vec commonAscentDirection(std::vector<Vec>& gradients, double scale) {
    if (gradients.empty()) return Vec::Zero();
    if (!(scale > 0) || !std::isfinite(scale)) return Vec::Zero();
    for (auto& g : gradients) g /= scale;
    Vec nearest = gradients.front();
    if (gradients.size() <= 4) {
        nearest = smallHullNearest(gradients);
    }
    for (int iteration = 0; gradients.size() > 4 && iteration < 256; ++iteration) {
        size_t index = 0;
        for (size_t i = 1; i < gradients.size(); ++i)
            if (nearest.dot(gradients[i]) < nearest.dot(gradients[index])) index = i;
        const Vec toward = gradients[index] - nearest;
        // Frank-Wolfe dual gap measures remaining objective improvement.
        const double gap = -nearest.dot(toward);
        if (gap <= 1e-6 * std::max(1e-12, nearest.squaredNorm()) || toward.squaredNorm() <= 1e-30) {
            break;
        }
        const double step = std::max(0.0, std::min(1.0, -nearest.dot(toward) / toward.squaredNorm()));
        if (step < 1e-14) { break; }
        nearest += step * toward;
    }
    if (nearest.norm() <= 1e-12) return Vec::Zero();
    const Vec direction = 0.5 * nearest.normalized();
    bool ascent = true;
    for (const auto& g : gradients)
        if (g.dot(direction) <= 1e-12) { ascent = false; break; }
    if (ascent) return direction;
    return Vec::Zero();
}

// Minimum-norm point in the convex hull of the near-worst gradients.
// A common ascent direction must have positive dot product with every one.
Vec activeSetDirection(const std::vector<SusTet>& tets, const Vec& x) {
    std::vector<double> qualities;
    double worst = DBL_MAX;
    for (const auto& tet : tets) {
        auto p = tet.points; p[tet.freeIndex] = x;
        qualities.push_back(signedQuality(p));
        worst = std::min(worst, qualities.back());
    }
    // Tighten the active band if including near-worst elements blocks ascent.
    for (double relativeBand : {1e-3, 1e-5, 0.0}) {
        std::vector<Vec> gradients;
        double scale = 0;
        const double band = std::max(1e-12, relativeBand * std::max(1e-3, std::fabs(worst)));
        for (size_t i = 0; i < tets.size(); ++i) {
            if (qualities[i] > worst + band) continue;
            const Vec g = qualityGradient(tets[i], x);
            if (!g.allFinite()) return Vec::Zero();
            gradients.push_back(g); scale = std::max(scale, g.norm());
        }
        const Vec direction = commonAscentDirection(gradients, scale);
        if (direction.squaredNorm() > 0) return direction;
    }
    return Vec::Zero();
}

double energy(const std::vector<SusTet>& tets, const Vec& x, double delta,
              Vec* gradient, bool* preservesValidity) {
    const double epsilon = 1e-10;
    if (gradient) gradient->setZero();
    if (preservesValidity) *preservesValidity = true;
    double sum = 0;
    for (const auto& tet : tets) {
        auto p = tet.points;
        p[tet.freeIndex] = x;
        const double det = sigma(p), norm = normSquared(p);
        if (preservesValidity && tet.positive && det <= 0) *preservesValidity = false;
        const double shifted = det - 2 * epsilon;
        const double root = std::hypot(shifted, 2 * delta);
        // Rationalize for negative determinants to avoid cancellation.
        const double h = det >= 0 ? 0.5 * (det + root) :
            2 * (delta * delta + epsilon * epsilon - epsilon * det) / (root - det);
        if (!(h > 0) || !(norm > 0) || !std::isfinite(h))
            return std::numeric_limits<double>::infinity();
        const double k = norm / (3 * std::pow(h, 2.0 / 3.0));
        sum += k;
        if (!gradient) continue;
        const Vec a = p[1] - p[0], b = p[2] - p[0], c = p[3] - p[0];
        std::array<Vec, 4> gd;
        gd[1] = b.cross(c); gd[2] = c.cross(a); gd[3] = a.cross(b);
        gd[0] = -gd[1] - gd[2] - gd[3];
        Vec gn = Vec::Zero();
        for (int j = 0; j < 4; ++j)
            if (j != tet.freeIndex) gn += x - p[j];
        // h'/h = (1 - epsilon/h) / root, also stable for det < 0.
        *gradient += k * (gn / norm - (2.0 / 3.0) *
            ((1 - epsilon / h) / root) * std::sqrt(2.0) * gd[tet.freeIndex]);
    }
    if (gradient) *gradient /= static_cast<double>(tets.size());
    return sum / static_cast<double>(tets.size());
}
}

double dt::DT::quality_sus(double* a, double* b, double* c, double* d) {
    std::array<Vec, 4> p = { Vec(a[0], a[1], a[2]), Vec(b[0], b[1], b[2]),
        Vec(c[0], c[1], c[2]), Vec(d[0], d[1], d[2]) };
    const Vec origin = p[0];
    double scale = 0;
    for (const auto& v : p) scale = std::max(scale, (v - origin).cwiseAbs().maxCoeff());
    if (!(scale > 0) || !std::isfinite(scale)) return 0;
    // Scale the original coordinates by a power of two BEFORE any subtraction.
    // Arbitrary division (even followed by exact orient3d) can reverse the sign
    // of a nearly coplanar tetrahedron by rounding its input coordinates.
    int exponent = 0;
    std::frexp(scale, &exponent);
    if (-exponent < std::numeric_limits<double>::max_exponent) {
        const double factor = std::scalbn(1.0, -exponent);
        for (auto& v : p) v *= factor;
    } else {
        // A subnormal edge scale can require an unrepresentable reciprocal.
        for (auto& v : p)
            for (int axis = 0; axis < 3; ++axis)
                v[axis] = std::scalbn(v[axis], -exponent);
    }
    return signedQuality(p);
}

int dt::DT::smooth_sus(int iNod, double minimumQualityFloor) {
    if (iNod < 0 || iNod >= static_cast<int>(Nodes.size()) || iNod == ghost ||
        isDelNod(iNod) || isbndpnt(iNod) || isCornerpnt(iNod) ||
        lockV.count(iNod) || periodic_P.count(iNod)) return 0;

    std::vector<int> sph;
    findSphere(iNod, sph);
    if (sph.empty()) return 0;
    std::unordered_set<int> neighbors;
    for (int t : sph) {
        if (t < 0 || t >= static_cast<int>(Elems.size()) || isDelEle(t) ||
            isvirtualtet(t) || (ghost != -1 && ishulltet(t))) return 0;
        for (int n : Elems[t].form) {
            if (n < 0 || n >= static_cast<int>(Nodes.size()) || isDelNod(n)) return 0;
            if (n != iNod) neighbors.insert(n);
        }
    }
    // SmoothPass preallocates slots before workers start. Direct serial calls
    // may grow the storage; callers must not run arbitrary adjacent points concurrently.
    if (susIdleStates.size() < Nodes.size() && !omp_in_parallel())
        susIdleStates.resize(Nodes.size());
    SusIdleState* const idleState = iNod < susIdleStates.size() ? &susIdleStates[iNod] : nullptr;
    std::vector<double> neighborhood;
    std::vector<int> neighborhoodNodes;
    if (idleState) {
        // Connectivity is exact; each distinct vertex's coordinates occur once.
        neighborhoodNodes.assign(neighbors.begin(), neighbors.end());
        neighborhoodNodes.push_back(iNod);
        std::sort(neighborhoodNodes.begin(), neighborhoodNodes.end());
        neighborhood.reserve(sph.size() * 5 + neighborhoodNodes.size() * 3 + 1);
        for (int t : sph) {
            neighborhood.push_back(t);
            for (int n : Elems[t].form) neighborhood.push_back(n);
        }
        for (int n : neighborhoodNodes)
            for (int j = 0; j < 3; ++j) neighborhood.push_back(Nodes[n].pt[j]);
        // A changed pass floor must invalidate the unchanged-neighborhood cache.
        neighborhood.push_back(minimumQualityFloor);
        if (!idleState->neighborhood.empty() && idleState->neighborhood == neighborhood && idleState->skips < 3) {
            ++idleState->skips;
            return 0;
        }
        *idleState = SusIdleState{};
    }
    auto rememberIdle = [&]() {
        if (!idleState) return;
        size_t k = sph.size() * 5;
        for (int n : neighborhoodNodes)
            for (int j = 0; j < 3; ++j) neighborhood[k++] = Nodes[n].pt[j];
        idleState->neighborhood.swap(neighborhood);
        idleState->skips = 0;
    };
    const Vec origin = Eigen::Map<Vec>(Nodes[iNod].pt);
    double scale = 0;
    for (int n : neighbors) scale = std::max(scale, (Eigen::Map<Vec>(Nodes[n].pt) - origin).norm());
    if (!(scale > 0) || !std::isfinite(scale)) return 0;
    // Use the same physical-coordinate quality as Elems[].q, not the
    // regularized objective or potentially stale cached qualities.
    // Same-color nodes never share an element; neighboring coordinates stay fixed.
    std::vector<std::array<Vec, 4>> physical(sph.size());
    std::vector<int> freeIndices(sph.size());
    for (size_t k = 0; k < sph.size(); ++k)
        for (int j = 0; j < 4; ++j) {
            const int n = Elems[sph[k]].form[j];
            physical[k][j] = Eigen::Map<Vec>(Nodes[n].pt);
            if (n == iNod) freeIndices[k] = j;
        }
    size_t firstQualityCheck = 0;
    // Below the cutoff the returned value is a rejection witness, not the
    // exact minimum. Accepted candidates always scan the complete sphere.
    auto minimumQuality = [&](const Vec& position, double cutoff) {
        double result = DBL_MAX;
        size_t worstIndex = firstQualityCheck;
        const size_t first = firstQualityCheck;
        for (size_t offset = 0; offset < physical.size(); ++offset) {
            const size_t k = (first + offset) % physical.size();
            auto p = physical[k];
            p[freeIndices[k]] = position;
            const double q = quality_sus(p[0].data(), p[1].data(), p[2].data(), p[3].data());
            if (!std::isfinite(q)) return -std::numeric_limits<double>::infinity();
            if (q < result) { result = q; worstIndex = k; }
            if (q < cutoff) { firstQualityCheck = k; return q; }
        }
        firstQualityCheck = worstIndex;
        return result;
    };
    const double initialMinQuality = minimumQuality(origin, -std::numeric_limits<double>::infinity());
    if (!std::isfinite(initialMinQuality)) return 0;
    const double allowedMinimum = std::min(initialMinQuality, minimumQualityFloor);
    double currentMinQuality = initialMinQuality;
    std::vector<SusTet> tets;
    double minimum = DBL_MAX, mean = 0;
    for (int t : sph) {
        SusTet tet;
        tet.freeIndex = -1;
        for (int j = 0; j < 4; ++j) {
            const int n = Elems[t].form[j];
            tet.points[j] = (Eigen::Map<Vec>(Nodes[n].pt) - origin) / scale;
            if (n == iNod) tet.freeIndex = j;
        }
        if (tet.freeIndex < 0) return 0;
        const double det = sigma(tet.points);
        if (!std::isfinite(det)) return 0;
        tet.positive = det > 0;
        minimum = std::min(minimum, det); mean += std::fabs(det);
        tets.push_back(tet);
    }
    mean /= tets.size();
    // Match sus/src/sus.cpp (1e6), not the older TechNote value (1e11).
    const double effective = std::numeric_limits<double>::epsilon() * 1e6;
    const double delta = std::max(std::sqrt(std::max(0.0, effective * (effective - minimum))), 1e-4 * mean);
    auto evaluateEnergy = [&](const Vec& position, Vec* gradient, bool* valid) {
        return energy(tets, position, delta, gradient, valid);
    };
    Vec x = Vec::Zero(), g;
    double f = evaluateEnergy(x, &g, nullptr);
    const double initial = f;
    if (!std::isfinite(f) || !g.allFinite()) return 0;
    Mat inverse = Mat::Identity();
    const double qualityTolerance = 1e-12; // Acceptance guard; not a stopping tolerance.
    // Sphere sweep: limit depth on valid meshes; keep the untangling budget.
    const int maxIterations = initialMinQuality > 0 ? 12 : 60;
    const int maxLineSearchTrials = 24;
    const int maxSusQualityRejections = 3;
    const double smallProgressTolerance = initialMinQuality > 0 ? 1e-4 : 1e-6;
    const int smallProgressLimit = initialMinQuality > 0 ? 2 : 3;
    const int progressBlock = initialMinQuality > 0 ? 3 : 6;
    const double blockQualityTolerance = initialMinQuality > 0 ? 3e-4 : 1e-5;
    const double directionResetStep = 0.125;
    bool activeMode = false;
    int smallProgressSteps = 0;
    double activeStep = 1;
    Vec previousActiveDirection = Vec::Zero();
    double stageQuality = initialMinQuality;
    for (int iter = 0; iter < maxIterations; ++iter) {
        Vec direction = -inverse * g;
        if (!direction.allFinite() || direction.dot(g) >= 0) { inverse.setIdentity(); direction = -g; }
        if (direction.norm() > 0.5) direction *= 0.5 / direction.norm();
        const double slope = g.dot(direction);
        double step = 1, nextF = f, nextMinQuality = currentMinQuality;
        Vec nextX = x;
        bool accepted = false, blockedByQuality = false, usedActiveSet = false;
        const bool stationary = g.norm() <= 1e-8 * std::max(1.0, std::fabs(f)) || !(slope < 0);
        int qualityRejections = 0;
        if (!activeMode && !stationary) {
            for (int ls = 0; ls < maxLineSearchTrials; ++ls, step *= 0.5) {
                nextX = x + step * direction;
                bool valid = false;
                nextF = evaluateEnergy(nextX, nullptr, &valid);
                if (valid && std::isfinite(nextF) && nextF <= f + 1e-4 * step * slope) {
                    nextMinQuality = minimumQuality(origin + scale * nextX, std::min(currentMinQuality, allowedMinimum));
                    if (nextMinQuality >= std::min(currentMinQuality, allowedMinimum)) { accepted = true; break; }
                    blockedByQuality = true;
                    // Allow several smaller SUS steps before changing objectives.
                    if (++qualityRejections >= maxSusQualityRejections) break;
                }
            }
        }
        // Exhausted validity/Armijo trials also need a different direction:
        // they never reach the minQ check and do not set blockedByQuality.
        // Retain the existing fallback for an accepted but negligible step.
        if (!accepted || (blockedByQuality && (nextX - x).norm() <= 1e-10)) {
            activeMode = true;
            {
                direction = activeSetDirection(tets, x);
            }
            accepted = false;
            if (direction.squaredNorm() > 0) {
                const bool changedDirection = previousActiveDirection.squaredNorm() > 0 &&
                    direction.dot(previousActiveDirection) < 0;
                // Reuse the accepted scale even after turning; permit faster growth
                // on a new direction without retrying every large step.
                step = std::min(changedDirection ? directionResetStep : 1.0,
                    (changedDirection ? 8.0 : 2.0) * activeStep);
                for (int ls = 0; ls < maxLineSearchTrials; ++ls, step *= 0.5) {
                    nextX = x + step * direction;
                    // Reject cheap quality failures before evaluating SUS energy.
                    nextMinQuality = minimumQuality(origin + scale * nextX, currentMinQuality + qualityTolerance);
                    if (!(nextMinQuality > currentMinQuality + qualityTolerance)) {
                        continue;
                    }
                    bool valid = false;
                    nextF = evaluateEnergy(nextX, nullptr, &valid);
                    if (!valid || !std::isfinite(nextF)) continue;
                    activeStep = step;
                    previousActiveDirection = direction;
                    accepted = true; usedActiveSet = true; break;
                }
            }
        }
        if (!accepted) break;
        Vec nextG = g;
        if (!usedActiveSet) {
            evaluateEnergy(nextX, &nextG, nullptr);
            if (!nextG.allFinite()) break;
        }
        const Vec s = nextX - x, y = nextG - g;
        const double curvature = s.dot(y);
        const double qualityGain = nextMinQuality - currentMinQuality;
        const double relativeEnergyChange = std::fabs(nextF - f) / std::max(1.0, std::fabs(f));
        // Stop only after repeated negligible progress; do not relax minQ acceptance.
        const bool smallProgress = qualityGain <= smallProgressTolerance &&
            (usedActiveSet || relativeEnergyChange <= smallProgressTolerance);
        smallProgressSteps = smallProgress ? smallProgressSteps + 1 : 0;
        x = nextX; g = nextG; f = nextF;
        currentMinQuality = nextMinQuality;
        if (s.norm() < 1e-10 || smallProgressSteps >= smallProgressLimit) break;
        // Review progress in blocks rather than spending the full budget on a plateau.
        if (activeMode && (iter + 1) % progressBlock == 0) {
            if (currentMinQuality - stageQuality <= blockQualityTolerance) break;
            stageQuality = currentMinQuality;
        }
        // Do not feed a max-min step into the SUS quasi-Newton update.
        if (!usedActiveSet && curvature > 1e-12 * s.norm() * y.norm()) {
            const Mat v = Mat::Identity() - s * y.transpose() / curvature;
            inverse = v * inverse * v.transpose() + s * s.transpose() / curvature;
        } else inverse.setIdentity();
    }
    const Vec result = origin + scale * x;
    if (!result.allFinite()) return 0;
    bool preservesValidity = false;
    const double verified = evaluateEnergy((result - origin) / scale, nullptr, &preservesValidity);
    if (!preservesValidity || !std::isfinite(verified)) return 0;
    std::vector<double> qualities;
    double finalMinQuality = DBL_MAX;
    for (int t : sph) {
        double* p[4];
        double candidate[3] = {result[0], result[1], result[2]};
        for (int j = 0; j < 4; ++j) p[j] = Elems[t].form[j] == iNod ? candidate : Nodes[Elems[t].form[j]].pt;
        const double q = quality_sus(p[0], p[1], p[2], p[3]);
        if (!std::isfinite(q)) return 0;
        finalMinQuality = std::min(finalMinQuality, q);
        qualities.push_back(q);
    }
    if (finalMinQuality < allowedMinimum) return 0;
    const bool improvedQuality = finalMinQuality > initialMinQuality + qualityTolerance;
    const bool improvedEnergy = verified < initial - 1e-12 * std::max(1.0, std::fabs(initial));
    if (!improvedQuality && !improvedEnergy) { rememberIdle(); return 0; }
    for (int j = 0; j < 3; ++j) Nodes[iNod].pt[j] = result[j];
    // Store the active pass metric, while the movement guard remains SUS.
    for (size_t j = 0; j < sph.size(); ++j) {
        if (improve_Metric == SUS_METRIC) Elems[sph[j]].q = qualities[j];
        else updateQuality(sph[j]);
    }
    if (finalMinQuality - initialMinQuality <= 1e-8 && x.norm() <= 1e-8) rememberIdle();
    return 1;
}

// Quality-driven interior smoothing must enter here, including topology repair.
// Other legacy metrics retain their previous SUS fallback.
int dt::DT::smoothInteriorPoint(int iNod, double minimumQualityFloor) {
    if(improve_Metric==2 || improve_Metric==6)
        return smoothWorstQualityPoints(std::vector<int>(1,iNod),improve_Metric);
    if (optimizationPolicy.enabled) {
        if (iNod < 0 || iNod >= static_cast<int>(Nodes.size()) || iNod == ghost ||
            isDelNod(iNod) || isbndpnt(iNod)) return 0;
        std::vector<int> star; findSphere(iNod, star);
        std::array<double, 3> origin; std::copy(Nodes[iNod].pt, Nodes[iNod].pt + 3, origin.begin());
        std::vector<double> qualities; qualities.reserve(star.size());
        for (int t : star) qualities.push_back(Elems[t].q);
        const int result = smooth_sus(iNod, 1e-12);
        if (result != 1) return result;
        // An idle/rejected SUS trial needs no A/J measurement. Reconstruct the
        // old star without writing shared coordinates after a successful trial.
        OptimizationQuality before;
        for (int t : star) {
            const double* p[4];
            for (int k = 0; k < 4; ++k) {
                const int n = Elems[t].form[k];
                p[k] = n == iNod ? origin.data() : Nodes[n].pt;
            }
            addOptimizationQuality(*this, before, p[0], p[1], p[2], p[3]);
        }
        before.valid = before.valid && before.count > 0;
        const auto after = measureOptimizationQuality(*this, &star);
        if (optimizationPointAccept(optimizationPolicy, before, after)) return 1;
        std::copy(origin.begin(), origin.end(), Nodes[iNod].pt);
        for (size_t i = 0; i < star.size(); ++i) Elems[star[i]].q = qualities[i];
        susIdleStates[iNod] = SusIdleState{};
        return 0;
    }
    if (improve_Metric == 8) {
        BoundarySmoothWorkspace boundary;
        return smoothOrthogonalPoint(iNod,minimumQualityFloor,boundary);
    }
    return smooth_sus(iNod, minimumQualityFloor);
}

int dt::DT::smooth_angle(int node, double minimumQualityFloor) {
    (void)minimumQualityFloor;
    return smoothWorstQualityPoints(std::vector<int>(1,node),2);
}

int dt::DT::smoothWorstQualityPoints(const std::vector<int>& nodes, int metric,
    BoundarySmoothWorkspace* suppliedBoundary) {
    if(nodes.empty() || (metric!=2 && metric!=6)) return 0;
    std::map<int,size_t> member;
    std::vector<Vec> origins;
    std::vector<int> cells;
    std::vector<std::array<int,3>> faces;
    std::vector<Vec> normals;
    auto addNormal=[&](Vec normal) {
        for(const auto& n:normals) normal-=n*normal.dot(n);
        const double length=normal.norm();
        if(length>1e-10) normals.push_back(normal/length);
    };
    BoundarySmoothWorkspace localBoundary;
    for(int node:nodes) {
        if(node<0 || node>=static_cast<int>(Nodes.size()) || node==ghost ||
            isDelNod(node) || isCornerpnt(node) || lockV.count(node) ||
            (nodes.size()==1 && periodic_P.count(node))) return 0;
        if(!member.emplace(node,origins.size()).second) return 0;
        origins.push_back(Vec(Nodes[node].pt));
        std::vector<int> star;
        if(isbndpnt(node)) {
            auto& boundary=suppliedBoundary && nodes.size()==1 ? *suppliedBoundary : localBoundary;
            if(!prepareBoundarySmoothPoint(node,boundary)) return 0;
            star=boundary.star;
            faces.insert(faces.end(),boundary.faces.begin(),boundary.faces.end());
            if(isFacetpnt(node)) addNormal(Vec(boundary.normal.data())/boundary.normalLength);
            else {
                const Vec axis=(Vec(Nodes[boundary.neighbors[0]].pt)-Vec(Nodes[boundary.neighbors[1]].pt)).normalized();
                const Vec normal=axis.unitOrthogonal();
                addNormal(normal);addNormal(axis.cross(normal));
            }
        } else {
            findSphere(node,star);
            // Dynamic insertions can leave a stale interior classification.
            // Never move a vertex through an actual live constraint face.
            for(int t:star) for(int opposite=0;opposite<4;++opposite) {
                if(Elems[t].form[opposite]==node) continue;
                int face[3],k=0;
                for(int j=0;j<4;++j)if(j!=opposite)face[k++]=Elems[t].form[j];
                const int* entry=BndTri.find(face[0],face[1],face[2]);
                if(entry && !isDelSurTri(*entry)) return 0;
            }
        }
        // Derive feature directions from live constrained edges as well as
        // the cached vertex type. Coplanar geometry-ID seams are constraints.
        std::set<int> neighbors;
        for(int t:star)for(int other:Elems[t].form)if(other!=node)neighbors.insert(other);
        for(int other:neighbors) {
            const int* edge=BndEdg.find(node,other);
            if(!edge || isDelSurEdg(*edge))continue;
            const auto& e=SurEdgs[*edge];bool feature=e.constrain>0 || e.face.size()!=2;
            if(e.face.size()==2)feature|=SurTris[e.face[0]].parent!=SurTris[e.face[1]].parent;
            if(!feature)continue;
            const Vec offset=Vec(Nodes[other].pt)-Vec(Nodes[node].pt);
            if(!(offset.squaredNorm()>0))return 0;
            const Vec axis=offset.normalized(),normal=axis.unitOrthogonal();
            addNormal(normal);addNormal(axis.cross(normal));
        }
        cells.insert(cells.end(),star.begin(),star.end());
    }
    if(normals.size()>=3) return 0;
    std::sort(cells.begin(),cells.end());cells.erase(std::unique(cells.begin(),cells.end()),cells.end());
    if(cells.empty()) return 0;
    Mat projection=Mat::Identity();
    for(const auto& normal:normals) projection-=normal*normal.transpose();
    double scale=0;
    for(int t:cells) {
        if(t<0 || t>=static_cast<int>(Elems.size()) || isDelEle(t) || ishulltet(t) || isvirtualtet(t)) return 0;
        for(int i=0;i<3;++i) for(int j=i+1;j<4;++j)
            scale=std::max(scale,(Vec(Nodes[Elems[t].form[i]].pt)-Vec(Nodes[Elems[t].form[j]].pt)).norm());
    }
    if(!(scale>0) || !std::isfinite(scale)) return 0;
    std::vector<detail::MinimaxCell> tets(cells.size());
    for(size_t i=0;i<cells.size();++i) {
        const auto& form=Elems[cells[i]].form;
        Vec anchor(Nodes[form[0]].pt);
        for(int k=0;k<4;++k) if(member.count(form[k])) {anchor=Vec(Nodes[form[k]].pt);break;}
        for(int k=0;k<4;++k) {
            tets[i].points[k]=(Vec(Nodes[form[k]].pt)-anchor)/scale;
            tets[i].moving[k]=member.count(form[k])!=0;
        }
    }
    auto point=[&](int node,const Vec& x)->Vec {
        const auto found=member.find(node);
        return found==member.end() ? Vec(Nodes[node].pt) : Vec(origins[found->second]+scale*x);
    };
    auto boundaryValid=[&](const Vec& x) {
        // A nearly stationary projected hull can amplify roundoff when its
        // direction is normalized. Verify actual displacements independently.
        for(size_t i=0;i<origins.size();++i) {
            const Vec displacement=Vec(origins[i]+scale*x)-origins[i];
            const double tolerance=32*std::numeric_limits<double>::epsilon()*origins[i].cwiseAbs().maxCoeff()+1e-12*scale;
            for(const auto& normal:normals)
                if(std::abs(normal.dot(displacement))>tolerance)return false;
        }
        for(const auto& face:faces) {
            const Vec a(Nodes[face[0]].pt),b(Nodes[face[1]].pt),c(Nodes[face[2]].pt);
            const Vec na=point(face[0],x),nb=point(face[1],x),nc=point(face[2],x);
            if(!pointOnBoundaryPlane(a.data(),b.data(),c.data(),na.data()) ||
                !pointOnBoundaryPlane(a.data(),b.data(),c.data(),nb.data()) ||
                !pointOnBoundaryPlane(a.data(),b.data(),c.data(),nc.data()) ||
                !((b-a).cross(c-a).dot((nb-na).cross(nc-na))>0)) return false;
        }
        return true;
    };
    auto evaluate=[&](const Vec& x,bool derivatives,std::vector<double>& values,std::vector<Vec>& gradients) {
        values.clear();gradients.clear();double worst=DBL_MAX;
        std::array<double,6> q;std::array<Vec,6> g;int count=0;
        for(const auto& tet:tets) {
            if(!detail::minimaxConstraints(tet,x,metric,derivatives,q,g,count)) return -DBL_MAX;
            for(int k=0;k<count;++k) {
                values.push_back(q[k]);worst=std::min(worst,q[k]);
                if(derivatives) gradients.push_back(projection*g[k]);
            }
        }
        return worst;
    };
    // Real-coordinate verification makes acceptance match the exported mesh,
    // including after normalized displacements are rounded back to doubles.
    auto actualMinimum=[&](const Vec& x) {
        double worst=DBL_MAX;std::vector<double> angles;
        for(int t:cells) {
            const auto& form=Elems[t].form;
            std::array<Vec,4> p={{point(form[0],x),point(form[1],x),point(form[2],x),point(form[3],x)}};
            if(!(GEOM_FUNC::orient3d(p[0].data(),p[1].data(),p[3].data(),p[2].data())>0)) return -DBL_MAX;
            double q=0;
            if(metric==2) {
                double high;
                if(!CalDihedral(p[0].data(),p[1].data(),p[3].data(),p[2].data(),q,high,angles,true)) return -DBL_MAX;
            } else q=ScaledJacobian(p[0].data(),p[1].data(),p[2].data(),p[3].data());
            if(!(q>0) || !std::isfinite(q)) return -DBL_MAX;
            worst=std::min(worst,q);
        }
        return worst;
    };
    Vec x=Vec::Zero();
    std::vector<double> values,candidateValues;
    std::vector<Vec> gradients,unused,active;
    values.reserve(tets.size()*6);gradients.reserve(tets.size()*6);candidateValues.reserve(tets.size()*6);
    const double initial=actualMinimum(x);
    if(!(initial>0)) return 0;
    double worst=initial;int smallProgress=0;
    // Klingner/Stellar-style nonsmooth max-min smoothing: active gradients,
    // closest convex-hull point, inactive-function intersection, then line search.
    // The 12 local iterations and 24 backtracks match the prior angle solver.
    for(int iteration=0;iteration<12;++iteration) {
        const double current=evaluate(x,true,values,gradients);
        if(!(current>0)) break;
        Vec direction=Vec::Zero();double rate=0,activeBand=0;
        for(double relativeBand:{.03,1e-3,1e-5,0.0}) {
            const double band=std::max(1e-12,relativeBand*current);
            active.clear();double gradientScale=0;
            for(size_t k=0;k<values.size();++k) if(values[k]<=current+band) {
                active.push_back(gradients[k]);gradientScale=std::max(gradientScale,gradients[k].norm());
            }
            direction=projection*commonAscentDirection(active,gradientScale);
            if(direction.norm()<=1e-12) {direction.setZero();continue;}
            direction.normalize();rate=DBL_MAX;
            for(size_t k=0;k<values.size();++k) if(values[k]<=current+band)
                rate=std::min(rate,gradients[k].dot(direction));
            if(rate>0) {activeBand=band;break;}
            direction.setZero();
        }
        if(direction.squaredNorm()==0) break;
        double step=.25;
        for(size_t k=0;k<values.size();++k) if(values[k]>current+activeBand) {
            const double slope=gradients[k].dot(direction);
            if(slope<rate) {
                const double intersection=(values[k]-current)/(rate-slope);
                if(intersection>0 && std::isfinite(intersection)) step=std::min(step,intersection);
            }
        }
        const double tolerance=std::max(1e-12,1e-10*std::abs(current));
        bool accepted=false;Vec next=x;double nextWorst=current;
        for(int trial=0;trial<24 && step>1e-14;++trial,step*=.5) {
            next=projection*(x+step*direction);
            if(!boundaryValid(next)) continue;
            nextWorst=evaluate(next,false,candidateValues,unused);
            if(nextWorst>current+tolerance && nextWorst-current>=.1*step*rate) {accepted=true;break;}
        }
        if(!accepted) break;
        smallProgress=nextWorst-current<1e-5 ? smallProgress+1 : 0;
        x=next;worst=nextWorst;
        if(smallProgress>=2 || step<1e-10) break;
    }
    const double tolerance=std::max(1e-12,1e-10*std::abs(initial));
    if(!(worst>initial+tolerance) || !boundaryValid(x) || !(actualMinimum(x)>initial+tolerance)) return 0;
    for(size_t i=0;i<nodes.size();++i) {
        const Vec result=origins[i]+scale*x;
        for(int axis=0;axis<3;++axis) Nodes[nodes[i]].pt[axis]=result[axis];
    }
    for(int t:cells) updateQuality(t);
    return 1;
}
