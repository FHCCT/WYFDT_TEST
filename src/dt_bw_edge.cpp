#include "dt.h"
#include "dt_bw.h"
#include <cmath>

namespace dt {
namespace {
struct EdgeQuality {
    int count = 0, tinyAngle = 0, tinyJacobian = 0;
    double angle = 180, jacobian = 1;
};

EdgeQuality tetQuality(DT& mesh, const double* a, const double* b,
    const double* c, const double* d) {
    EdgeQuality q;
    double low = 0, high = 0;
    std::vector<double> angles;
    mesh.CalDihedral(const_cast<double*>(a), const_cast<double*>(b),
        const_cast<double*>(d), const_cast<double*>(c), low, high, angles);
    q.angle = RADIO2ANGLE(low);
    q.jacobian = mesh.ScaledJacobian(const_cast<double*>(a), const_cast<double*>(b),
        const_cast<double*>(c), const_cast<double*>(d));
    q.count = 1;
    q.tinyAngle = q.angle < 0.1;
    q.tinyJacobian = q.jacobian < 0.001;
    return q;
}

void addQuality(EdgeQuality& total, const EdgeQuality& q) {
    total.count += q.count;
    total.tinyAngle += q.tinyAngle;
    total.tinyJacobian += q.tinyJacobian;
    total.angle = std::min(total.angle, q.angle);
    total.jacobian = std::min(total.jacobian, q.jacobian);
}

EdgeQuality fillQuality(DT& mesh, const BWPlan& plan) {
    EdgeQuality q;
    for (const BWFace& face : plan.faces) {
        const auto& v = face.vertices;
        if (face.geo == DT::virtualID ||
            std::find(v.begin(), v.end(), mesh.ghost) != v.end()) continue;
        addQuality(q, tetQuality(mesh, plan.request.point.data(),
            mesh.Nodes[v[0]].pt, mesh.Nodes[v[1]].pt, mesh.Nodes[v[2]].pt));
    }
    return q;
}

void addRetainedCell(DT& mesh, EdgeQuality& q, int tet) {
    if (mesh.ishulltet(tet) || mesh.isvirtualtet(tet)) return;
    const int* v = mesh.Elems[tet].form;
    addQuality(q, tetQuality(mesh, mesh.Nodes[v[0]].pt, mesh.Nodes[v[1]].pt,
        mesh.Nodes[v[2]].pt, mesh.Nodes[v[3]].pt));
}

bool betterPosition(DT& mesh, const BWPlan& best, EdgeQuality a,
    const BWPlan& trial, EdgeQuality b) {
    // Different positions can produce different cavities. Compare both results
    // on their union, including old cells retained by only one of the plans.
    for (int t : trial.cavity)
        if (std::find(best.cavity.begin(), best.cavity.end(), t) == best.cavity.end())
            addRetainedCell(mesh, a, t);
    for (int t : best.cavity)
        if (std::find(trial.cavity.begin(), trial.cavity.end(), t) == trial.cavity.end())
            addRetainedCell(mesh, b, t);
    if (!b.count || !(b.angle > 0 && b.jacobian > 0) ||
        !std::isfinite(b.angle) || !std::isfinite(b.jacobian) ||
        b.tinyAngle > a.tinyAngle || b.tinyJacobian > a.tinyJacobian ||
        b.angle + 1e-10 < std::min(a.angle, 0.1) ||
        b.jacobian + 1e-12 < std::min(a.jacobian, 0.001)) return false;
    // A geometric mean treats both dimensionless tail measures symmetrically.
    return std::log(b.angle / 70.5287793655093) + std::log(b.jacobian) >
        std::log(a.angle / 70.5287793655093) + std::log(a.jacobian) + 1e-7;
}
}

BWStatus planBWOnEdge(DT& mesh, const BWRequest& request, BWPlan& plan) {
    const int a = request.midpointEdge[0], b = request.midpointEdge[1];
    if (request.fixedEdgePoint || request.info == 0 || request.anisotropic || a < 0 || b < 0 || a == b ||
        a >= static_cast<int>(mesh.Nodes.size()) || b >= static_cast<int>(mesh.Nodes.size()) ||
        a == mesh.ghost || b == mesh.ghost || mesh.isDelNod(a) || mesh.isDelNod(b))
        return planBW(mesh, request, plan);
    // Own the request even if the caller borrows it from this plan's workspace.
    const BWRequest original = request;
    BWStatus status = planBW(mesh, original, plan);
    EdgeQuality bestQuality;
    if (status == BWStatus::Ready) {
        bestQuality = fillQuality(mesh, plan);
        if (bestQuality.angle >= 0.1 && bestQuality.jacobian >= 0.001) return status;
    }
    // Only poor/failed edge candidates pay for at most four extra BW plans.
    // Keep the original valid point when no better valid position is found.
    const double positions[] = {0.49, 0.51, 0.48, 0.52};
    BWPlan trial;
    for (double fraction : positions) {
        BWRequest candidate = original;
        for (int k = 0; k < 3; ++k)
            candidate.point[k] = (1 - fraction) * mesh.Nodes[a].pt[k] + fraction * mesh.Nodes[b].pt[k];
        candidate.space = (1 - fraction) * mesh.Nodes[a].space + fraction * mesh.Nodes[b].space;
        if (planBW(mesh, candidate, trial) != BWStatus::Ready) continue;
        const EdgeQuality q = fillQuality(mesh, trial);
        if (status != BWStatus::Ready || betterPosition(mesh, plan, bestQuality, trial, q)) {
            std::swap(plan, trial);
            bestQuality = q;
            status = BWStatus::Ready;
        }
    }
    return status;
}
}
