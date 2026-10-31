#include "dt.h"
#include <stdexcept>

namespace dt {

int DT::insertPlanarBoundaryPoint(int face, const std::array<double, 3>& weights) {
    if (face < 0 || face >= static_cast<int>(SurTris.size()) || isDelSurTri(face) || lockF.count(face)) return 0;
    const SurTri original = SurTris[face];
    const int* constraint = BndTri.find(original.form[0], original.form[1], original.form[2]);
    if (!constraint || *constraint != face) return 0;
    std::array<int, 3> edges;
    for (int k = 0; k < 3; ++k) {
        const int a = original.form[k], b = original.form[(k + 1) % 3];
        const int* edge = BndEdg.find(a, b);
        if (!edge || isDelSurEdg(*edge) || lockE.count(*edge) || periodic_P.count(a)) return 0;
        edges[k] = *edge;
    }
    int owner = -1;
    if (!isMeshFace(original.form[0], original.form[1], original.form[2], &owner)) return 0;
    int opposite = -1;
    for (int k = 0; k < 4; ++k)
        if (std::find(original.form, original.form + 3, Elems[owner].form[k]) == original.form + 3) opposite = k;
    if (opposite < 0) return 0;
    const int neighbor = getNeig(owner, opposite);
    if (neighbor < 0 || isDelEle(neighbor)) return 0;
    std::array<double, 3> point = {{0, 0, 0}};
    double space = 0;
    for (int k = 0; k < 3; ++k) {
        space += weights[k] * Nodes[original.form[k]].space;
        for (int j = 0; j < 3; ++j) point[j] += weights[k] * Nodes[original.form[k]].pt[j];
    }
    const int node = addNode(point[0], point[1], point[2], space);
    setbndpnt(node);
    // Only this triangle may be replaced. Other surface triangles, feature
    // edges and all region interfaces remain barriers to BW expansion.
    struct RestoreFace {
        DT& mesh; int face, node; SurTri original; bool committed = false;
        RestoreFace(DT& m, int f, int n, const SurTri& t) : mesh(m), face(f), node(n), original(t) {}
        ~RestoreFace() {
            if (!committed) {
                mesh.BndTri.add(original.form[0], original.form[1], original.form[2], face);
                mesh.DelNod(node);
            }
        }
    } restore{*this, face, node, original};
    BndTri.erase(original.form[0], original.form[1], original.form[2]);
    BWRequest request = makeBWRequest(*this, node, {owner, neighbor}, 3);
    request.trackAccess = false;
    request.initShell = 2;
    request.acceptance = BWAcceptance::VolumeOnly;
    std::copy(original.form, original.form + 3, request.sourceFace.begin());
    BWPlan& plan = serialBW.plan;
    if (planBW(*this, request, plan) != BWStatus::Ready) return 0;
    if (std::find(plan.cavity.begin(), plan.cavity.end(), owner) == plan.cavity.end() ||
        std::find(plan.cavity.begin(), plan.cavity.end(), neighbor) == plan.cavity.end()) return 0;
    const int regions[2] = {Elems[owner].geo, Elems[neighbor].geo};
    double oldMinimumVolume = DBL_MAX;
    for (int cell : plan.cavity) {
        if (Elems[cell].geo != regions[0] && Elems[cell].geo != regions[1]) return 0;
        if (!ishulltet(cell) && !isvirtualtet(cell)) {
            const double volume = calVolume(cell);
            if (!std::isfinite(volume) || !(volume > 0)) return 0;
            oldMinimumVolume = std::min(oldMinimumVolume, volume);
        }
    }
    if (oldMinimumVolume == DBL_MAX) return 0;
    const double floor = std::max(minVolume_bw, 0.1 * oldMinimumVolume);
    std::array<int, 3> childUses = {{0, 0, 0}};
    for (const BWFace& f : plan.faces) {
        if (f.geo != regions[0] && f.geo != regions[1]) return 0;
        for (int k = 0; k < 3; ++k)
            if (std::find(f.vertices.begin(), f.vertices.end(), original.form[k]) != f.vertices.end() &&
                std::find(f.vertices.begin(), f.vertices.end(), original.form[(k + 1) % 3]) != f.vertices.end()) ++childUses[k];
        if (f.geo == virtualID || std::find(f.vertices.begin(), f.vertices.end(), ghost) != f.vertices.end()) continue;
        const double volume = GEOM_FUNC::orient3d(point.data(), Nodes[f.vertices[0]].pt,
            Nodes[f.vertices[2]].pt, Nodes[f.vertices[1]].pt);
        if (!std::isfinite(volume) || !(volume > floor)) return 0;
    }
    // Every proposed child surface must have a cell on each side (one side
    // may be the exterior hull), not just a positive-volume tetrahedron.
    if (std::find_if(childUses.begin(), childUses.end(), [](int n) { return n != 2; }) != childUses.end()) return 0;
    SurTris.reserve(SurTris.size() + 3);
    SurEdgs.reserve(SurEdgs.size() + 3);
    auto& slots = serialBW.slots;
    slots.clear(); slots.reserve(plan.faces.size());
    try {
        for (size_t k = 0; k < plan.faces.size(); ++k) slots.push_back(addElem());
    } catch (...) {
        for (int cell : slots) DelEle(cell);
        throw;
    }
    commitBW(*this, plan, node, slots);
    finishBW(*this, plan);
    restore.committed = true;
    setDelSurTri(face);
    std::array<int, 3> children;
    for (int k = 0; k < 3; ++k) {
        children[k] = static_cast<int>(SurTris.size());
        SurTri triangle = {};
        triangle.form[0] = original.form[k]; triangle.form[1] = original.form[(k + 1) % 3]; triangle.form[2] = node;
        triangle.parent = original.parent; triangle.info = 0;
        SurTris.push_back(triangle);
        BndTri.add(triangle.form[0], triangle.form[1], triangle.form[2], children[k]);
        for (int& incident : SurEdgs[edges[k]].face) if (incident == face) { incident = children[k]; break; }
    }
    for (int k = 0; k < 3; ++k) {
        const int edge = static_cast<int>(SurEdgs.size());
        SurEdg spoke;
        spoke.iStart = original.form[k]; spoke.iEnd = node;
        spoke.face = {children[k], children[(k + 2) % 3]};
        SurEdgs.push_back(std::move(spoke));
        BndEdg.add(original.form[k], node, edge);
    }
    FacetNum[original.parent] += 2;
    return 1;
}

int DT::refineHighDegreePlanarPoints(int maxDegree) {
    if (omp_in_parallel()) throw std::logic_error("Planar degree refinement must run outside a parallel region");
    if (maxDegree < 0) throw std::invalid_argument("Planar degree threshold must be nonnegative");
    if (!maxDegree || !modifyBnd || !AniSol.empty()) return 0;
    MeshStageLog stage(*this, "Planar degree refinement", 2);
    // Freeze the original high-degree points. New vertices and points pushed
    // over the threshold by this pass do not start a chain of subdivisions.
    std::vector<int> degree(Nodes.size(), 0), candidates;
    for (int f = 0; f < static_cast<int>(SurTris.size()); ++f) if (!isDelSurTri(f))
        for (int n : SurTris[f].form) if (n != ghost && isFacetpnt(n)) ++degree[n];
    for (int n = 0; n < static_cast<int>(degree.size()); ++n)
        if (degree[n] > maxDegree && !isDelNod(n) && !lockV.count(n) && !periodic_P.count(n)) candidates.push_back(n);
    std::stable_sort(candidates.begin(), candidates.end(), [&](int a, int b) { return degree[a] > degree[b]; });
    int inserted = 0;
    for (int n : candidates) {
        std::unordered_set<int> star;
        findSphere_tri(n, star);
        std::vector<int> faces(star.begin(), star.end());
        std::sort(faces.begin(), faces.end());
        int largest = -1, geometry = -1;
        double maximumArea = 0, scale = 0;
        bool eligible = true;
        for (int f : faces) {
            if (isDelSurTri(f)) continue;
            if (lockF.count(f)) { eligible = false; break; }
            if (geometry == -1) geometry = SurTris[f].parent;
            if (SurTris[f].parent != geometry) { eligible = false; break; }
            const auto& v = SurTris[f].form;
            const double area = calArea(Nodes[v[0]].pt, Nodes[v[1]].pt, Nodes[v[2]].pt);
            if (area > maximumArea) { maximumArea = area; largest = f; }
            for (int p : v) scale = std::max(scale, distance(Nodes[n].pt, Nodes[p].pt));
        }
        if (!eligible || largest < 0 || !(scale > 0)) continue;
        const auto& triangle = SurTris[largest].form;
        const Eigen::Vector3d a(Nodes[triangle[0]].pt), b(Nodes[triangle[1]].pt), c(Nodes[triangle[2]].pt);
        const Eigen::Vector3d normal = (b - a).cross(c - a);
        if (!(normal.norm() > 1e-12 * scale * scale)) continue;
        for (int f : faces) if (!isDelSurTri(f)) for (int p : SurTris[f].form)
            if (std::abs(normal.dot(Eigen::Vector3d(Nodes[p].pt) - a)) > 1e-10 * normal.norm() * scale) eligible = false;
        if (!eligible) continue;
        const std::array<double, 3> weights = {{1.0 / 3, 1.0 / 3, 1.0 / 3}};
        inserted += insertPlanarBoundaryPoint(largest, weights) > 0;
    }
    if (infolevel > 1) meshLogger->debug("Planar degree insertion: {}/{} (degree > {})", inserted, candidates.size(), maxDegree);
    return inserted;
}

} // namespace dt
