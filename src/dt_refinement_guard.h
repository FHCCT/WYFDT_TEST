#pragma once
#include "dt.h"

namespace dt {

namespace detail {
struct BoundaryRefinementSector {
    int faces = 0;
    double angle = 0;
    bool bounded = false;
};

inline double boundaryRefinementCornerAngle(DT& mesh, int face, int vertex) {
    const auto& points = mesh.SurTris[face].form;
    int position = 0;
    while (position < 3 && points[position] != vertex) ++position;
    if (position == 3) return -1;
    const Eigen::Vector3d origin(mesh.Nodes[vertex].pt);
    const Eigen::Vector3d a = Eigen::Vector3d(mesh.Nodes[points[(position + 1) % 3]].pt) - origin;
    const Eigen::Vector3d b = Eigen::Vector3d(mesh.Nodes[points[(position + 2) % 3]].pt) - origin;
    if (!(a.squaredNorm() > 0 && b.squaredNorm() > 0)) return -1;
    return std::atan2(a.cross(b).norm(), a.dot(b)) * 180.0 / std::acos(-1.0);
}

inline BoundaryRefinementSector boundaryRefinementSector(DT& mesh, int face, int vertex) {
    BoundaryRefinementSector result;
    const int geometry = mesh.SurTris[face].parent;
    std::vector<int> faces(1, face);
    std::unordered_set<int> seen{face};
    std::set<int> rays;
    for (size_t i = 0; i < faces.size(); ++i) {
        // An ambiguous/large fan supplies no restriction; do not guess its angle.
        if (faces.size() > 128) return BoundaryRefinementSector{};
        const int current = faces[i];
        const double angle = boundaryRefinementCornerAngle(mesh, current, vertex);
        if (!(angle >= 0) || !std::isfinite(angle)) return BoundaryRefinementSector{};
        result.angle += angle;
        for (int other : mesh.SurTris[current].form) if (other != vertex) {
            const int* edge = mesh.BndEdg.find(vertex, other);
            if (!edge || mesh.isDelSurEdg(*edge)) return BoundaryRefinementSector{};
            const auto& boundary = mesh.SurEdgs[*edge];
            int count = 0, next = -1;
            for (int neighbor : boundary.face)
                if (!mesh.isDelSurTri(neighbor) && mesh.SurTris[neighbor].parent == geometry) {
                    ++count;
                    if (neighbor != current) next = neighbor;
                }
            // Feature edges and nonmanifold branches bound distinct sectors.
            if (boundary.constrain > 0 || count != 2) {
                rays.insert(other);
                continue;
            }
            if (next >= 0 && seen.insert(next).second) faces.push_back(next);
        }
    }
    result.faces = static_cast<int>(faces.size());
    result.bounded = rays.size() == 2;
    return result;
}

// Bisection preserves the angle sum at an existing vertex and cannot reduce
// its triangle count. An already overfull feature sector cannot meet the same
// boundary-angle target by adding more points. Flips/smoothing can still act;
// re-evaluate their current result rather than permanently marking geometry.
inline bool boundaryQualityRefinementStagnates(DT& mesh, int face, double target) {
    for (int vertex : mesh.SurTris[face].form) {
        if (!(boundaryRefinementCornerAngle(mesh, face, vertex) < target)) continue;
        const auto sector = boundaryRefinementSector(mesh, face, vertex);
        const double required = target * sector.faces;
        if (sector.bounded && sector.angle + 1e-8 * std::max(sector.angle, required) < required)
            return true;
    }
    return false;
}
} // namespace detail

// Quality repair stays within the requested edge's incident triangles.
// Following an increasing-edge chain can split remote, otherwise good areas.
inline int longestAdjacentQualityBoundaryEdge(DT& mesh, int edge) {
    if (edge < 0 || mesh.isDelSurEdg(edge)) return -1;
    int best = edge;
    double maximum = mesh.distance2(mesh.Nodes[mesh.SurEdgs[edge].iStart].pt,
        mesh.Nodes[mesh.SurEdgs[edge].iEnd].pt);
    for (int face : mesh.SurEdgs[edge].face) {
        if (mesh.isDelSurTri(face) || mesh.lockF.count(face)) continue;
        const auto& points = mesh.SurTris[face].form;
        for (int k = 0; k < 3; ++k) {
            const int* candidate = mesh.BndEdg.find(points[k], points[(k + 1) % 3]);
            if (!candidate || mesh.isDelSurEdg(*candidate) || mesh.lockE.count(*candidate)) continue;
            bool locked = false;
            for (int neighbor : mesh.SurEdgs[*candidate].face)
                locked |= mesh.lockF.count(neighbor) != 0;
            if (locked) continue;
            const double length = mesh.distance2(mesh.Nodes[points[k]].pt,
                mesh.Nodes[points[(k + 1) % 3]].pt);
            if (length > maximum * (1 + 1e-12)) {
                maximum = length;
                best = *candidate;
            }
        }
    }
    return best;
}

// Follow increasing edge lengths before bisection. Splitting a short edge
// repeatedly while leaving a long side intact creates a thin triangle fan.
inline int longestRefinementBoundaryEdge(DT& mesh, int edge) {
    const auto lengthSquared = [&](int e) {
        return mesh.distance2(mesh.Nodes[mesh.SurEdgs[e].iStart].pt,
            mesh.Nodes[mesh.SurEdgs[e].iEnd].pt);
    };
    for (;;) {
        if (edge < 0 || edge >= static_cast<int>(mesh.SurEdgs.size()) ||
            mesh.isDelSurEdg(edge) || mesh.lockE.count(edge)) return -1;
        int longest = edge;
        double maximum = lengthSquared(edge);
        for (int face : mesh.SurEdgs[edge].face) {
            if (mesh.isDelSurTri(face)) return -1;
            if (mesh.lockF.count(face)) return -1;
            const auto& p = mesh.SurTris[face].form;
            for (int k = 0; k < 3; ++k) {
                const int* found = mesh.BndEdg.find(p[k], p[(k+1)%3]);
                if (!found || mesh.isDelSurEdg(*found)) return -1;
                const double length = lengthSquared(*found);
                if (length > maximum * (1 + 1e-12)) {
                    maximum = length;
                    longest = *found;
                }
            }
        }
        if (longest == edge) return edge;
        edge = longest;
    }
}

// Used only by the serial size-refinement pass. Rebuild after smoothing; during
// refinement, append new edges and ignore deleted ones without changing the BW API.
class RefinementBoundaryEdges {
public:
    // Redirect only within 0.1% of the shorter insertion/boundary edge.
    static double proximityRatio() { return 0.001; }

    explicit RefinementBoundaryEdges(DT& mesh) : mesh_(mesh) {
        typename Tree::BoundingBox box;
        for (int k = 0; k < 3; ++k) { box[k] = DBL_MAX; box[k + 3] = -DBL_MAX; }
        bool found = false;
        for (int e = 0; e < static_cast<int>(mesh_.SurEdgs.size()); ++e) {
            if (mesh_.isDelSurEdg(e)) continue;
            const auto& edge = mesh_.SurEdgs[e];
            for (int n : {edge.iStart, edge.iEnd}) for (int k = 0; k < 3; ++k) {
                box[k] = std::min(box[k], mesh_.Nodes[n].pt[k]);
                box[k + 3] = std::max(box[k + 3], mesh_.Nodes[n].pt[k]);
            }
            found = true;
        }
        if (!found) return;
        tree_.reset(new Tree(box));
        appendNewEdges();
    }

    void appendNewEdges() {
        if (!tree_) return;
        for (; nextEdge_ < mesh_.SurEdgs.size(); ++nextEdge_) {
            const int e = static_cast<int>(nextEdge_);
            if (mesh_.isDelSurEdg(e)) continue;
            const auto& edge = mesh_.SurEdgs[e];
            Segment segment;
            segment.row(0) = Eigen::Map<const Eigen::Vector3d>(mesh_.Nodes[edge.iStart].pt);
            segment.row(1) = Eigen::Map<const Eigen::Vector3d>(mesh_.Nodes[edge.iEnd].pt);
            if (!(segment.row(1) - segment.row(0)).squaredNorm()) continue;
            // BinaryTree requires dense insertion IDs, independently of mesh edge IDs.
            const size_t id = edgeIds_.size();
            edgeIds_.push_back(e);
            tree_->insert(segment, id);
        }
    }

    int nearby(const double* point, double insertionEdgeLength, int region) {
        if (!tree_ || !(insertionEdgeLength > 0) || !std::isfinite(insertionEdgeLength)) return -1;
        appendNewEdges();
        const Eigen::Vector3d p = Eigen::Map<const Eigen::Vector3d>(point);
        if (!p.allFinite()) return -1;
        typename Tree::BoxVector queryPoint = p.transpose();
        tree_->query(queryPoint, proximityRatio() * insertionEdgeLength, candidates_);
        int best = -1;
        double bestDistance = DBL_MAX;
        for (size_t id : candidates_) {
            const int e = edgeIds_[id];
            if (mesh_.isDelSurEdg(e)) continue;
            const auto& edge = mesh_.SurEdgs[e];
            const Eigen::Vector3d a = Eigen::Map<const Eigen::Vector3d>(mesh_.Nodes[edge.iStart].pt);
            const Eigen::Vector3d ab = Eigen::Map<const Eigen::Vector3d>(mesh_.Nodes[edge.iEnd].pt) - a;
            const double lengthSquared = ab.squaredNorm();
            if (!(lengthSquared > 0)) continue;
            const double projection = (p - a).dot(ab) / lengthSquared;
            // Proximity to an endpoint alone is a vertex issue, not edge encroachment.
            if (!(projection > 0 && projection < 1)) continue;
            const double distanceSquared = (p - a - projection * ab).squaredNorm();
            const double limit = proximityRatio() * std::min(insertionEdgeLength, std::sqrt(lengthSquared));
            if (!(distanceSquared <= limit * limit) || distanceSquared > bestDistance) continue;
            // Do not split an edge belonging only to a different solid/region.
            mesh_.findSphere(edge.iStart, star_);
            bool touchesRegion = false;
            for (int t : star_) {
                if (mesh_.isDelEle(t) || mesh_.ishulltet(t) || mesh_.isvirtualtet(t)) continue;
                if (mesh_.Elems[t].geo == region && mesh_.isNod_in_Tet(edge.iEnd, t) >= 0) {
                    touchesRegion = true;
                    break;
                }
            }
            if (!touchesRegion) continue;
            if (distanceSquared == bestDistance && best >= 0) {
                const bool feature = edge.constrain > 0;
                const bool bestFeature = mesh_.SurEdgs[best].constrain > 0;
                if (feature != bestFeature ? !feature : e > best) continue;
            }
            best = e;
            bestDistance = distanceSquared;
        }
        return best;
    }

    bool protectedEdge(int edge) const {
        if (mesh_.lockE.count(edge)) return true;
        for (int face : mesh_.SurEdgs[edge].face) if (mesh_.lockF.count(face)) return true;
        return false;
    }

private:
    using Segment = Eigen::Matrix<double, 2, 3>;
    using Tree = BinaryTree<Segment>;
    DT& mesh_;
    std::unique_ptr<Tree> tree_;
    size_t nextEdge_ = 0;
    std::vector<int> edgeIds_, star_;
    std::vector<size_t> candidates_;
};
}
