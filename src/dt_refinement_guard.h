#pragma once
#include "dt.h"

namespace dt {

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
