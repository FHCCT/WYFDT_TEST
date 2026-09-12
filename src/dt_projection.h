#pragma once
#include <Eigen/Dense>
#include <vector>
#include "../extern/tree/dt_binary_tree.hpp"

namespace dt {
// Reference geometry belongs to one mesh. Build/reset between algorithm calls;
// parallel projection queries only read these arrays and trees.
struct FineMeshProjection {
    using Edge = Eigen::Matrix<double, 2, 3>;
    using EdgeTree = BinaryTree<Edge>;
    bool enabled = false;
    std::vector<Eigen::Matrix3d> triangles;
    BinaryTree<Eigen::Matrix3d> tree;
    bool edgesEnabled = false;
    std::vector<Edge> edges;
    EdgeTree edgeTree;
};
}
