#include "dt.h"

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
        if (nearest.norm() <= 1e-12) continue;
        const Vec direction = 0.5 * nearest.normalized();
        bool ascent = true;
        for (const auto& g : gradients)
            if (g.dot(direction) <= 1e-12) { ascent = false; break; }
        if (ascent) return direction;
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
    for (auto& v : p) { v -= origin; scale = std::max(scale, v.cwiseAbs().maxCoeff()); }
    if (!(scale > 0) || !std::isfinite(scale)) return 0;
    for (auto& v : p) v /= scale;
    return signedQuality(p);
}

int dt::DT::smooth_sus(int iNod) {
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
    if (susReuseIdle && susIdleStates.size() < Nodes.size() && !omp_in_parallel())
        susIdleStates.resize(Nodes.size());
    SusIdleState* idleState = susReuseIdle && iNod < susIdleStates.size() ? &susIdleStates[iNod] : nullptr;
    std::vector<double> neighborhood;
    if (idleState) {
        neighborhood.reserve(sph.size() * 17);
        for (int t : sph) {
            neighborhood.push_back(t);
            for (int n : Elems[t].form) {
                neighborhood.push_back(n);
                for (int j = 0; j < 3; ++j) neighborhood.push_back(Nodes[n].pt[j]);
            }
        }
        if (!idleState->neighborhood.empty() && idleState->neighborhood == neighborhood && idleState->skips < 3) {
            ++idleState->skips;
            return 0;
        }
        *idleState = SusIdleState{};
    }
    auto rememberIdle = [&]() {
        if (!idleState) return;
        size_t k = 0;
        for (int t : sph) {
            ++k;
            for (int n : Elems[t].form) {
                ++k;
                for (int j = 0; j < 3; ++j) neighborhood[k++] = Nodes[n].pt[j];
            }
        }
        *idleState = SusIdleState{neighborhood, 0};
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
                    nextMinQuality = minimumQuality(origin + scale * nextX, currentMinQuality);
                    if (nextMinQuality >= currentMinQuality) { accepted = true; break; }
                    blockedByQuality = true;
                    // Allow several smaller SUS steps before changing objectives.
                    if (++qualityRejections >= maxSusQualityRejections) break;
                }
            }
        }
        // Try smaller SUS steps first. If minQ blocks progress (including
        // a numerically tiny step), switch objectives at the current position.
        if (activeMode || stationary || (blockedByQuality && (!accepted || (nextX - x).norm() <= 1e-10))) {
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
    if (finalMinQuality < initialMinQuality) return 0;
    const bool improvedQuality = finalMinQuality > initialMinQuality + qualityTolerance;
    const bool improvedEnergy = verified < initial - 1e-12 * std::max(1.0, std::fabs(initial));
    if (!improvedQuality && !improvedEnergy) { rememberIdle(); return 0; }
    for (int j = 0; j < 3; ++j) Nodes[iNod].pt[j] = result[j];
    for (size_t j = 0; j < sph.size(); ++j) Elems[sph[j]].q = qualities[j];
    if (finalMinQuality - initialMinQuality <= 1e-8 && x.norm() <= 1e-8) rememberIdle();
    return 1;
}
