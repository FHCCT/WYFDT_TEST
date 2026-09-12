#ifndef SRC_WYFDT_DT_BINARY_TREE_HPP_
#define SRC_WYFDT_DT_BINARY_TREE_HPP_
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <array>
#include <bitset>
#include <iostream>
#include <limits>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
#define DIM DerivedV::ColsAtCompileTime
// #pragma optimize("", off)
namespace dt {
		/**
		 * @brief The axis aligned spatial tree for all
		 * DerivedV::ColsAtCompileTimeensional simplex.
		 *
		 * For instance:
		 * simplex: triangle , tree: octree
		 * simplex: edge,      tree: quadtree
		 *
		 *
		 * @tparam DerivedV simplex, decribed by a m*n matrix, while m=the number of
		 * points, and n = DerivedV::ColsAtCompileTime
		 */
		template <typename DerivedV>
		class BinaryTree {
		public:
			// using DerivedV::ColsAtCompileTime =
			// DerivedV::ColsAtCompileTime;
			using BoxVector = Eigen::Matrix<double, 1, DerivedV::ColsAtCompileTime>;
			using Vertices = Eigen::Matrix<double, -1, DerivedV::ColsAtCompileTime>;
			using Topos = Eigen::Matrix<int, -1, DerivedV::RowsAtCompileTime>;
			using BoundingBox = std::array<double, 2 * DerivedV::ColsAtCompileTime>;
			using HashArray = std::array<size_t, 1 << DerivedV::ColsAtCompileTime>;

			BinaryTree() {}
			explicit BinaryTree(
				const std::array<double, DerivedV::ColsAtCompileTime * 2>& initial_box) {
				tree_node_array_.emplace_back(0, std::numeric_limits<size_t>::max(),
					initial_box);
				generateChildHash();
				diag_len_of_box_ = 0;
				for (int i = 0; i < DIM; i++) {
					diag_len_of_box_ += (initial_box[i] - initial_box[DIM + i]) *
						(initial_box[i] - initial_box[DIM + i]);
				}
				diag_len_of_box_ = std::sqrt(diag_len_of_box_);
			}

			explicit BinaryTree(const Vertices& V, const Topos& F) {
				std::array<double, 2 * DerivedV::ColsAtCompileTime> box_array;
				for (size_t i = 0; i < DerivedV::ColsAtCompileTime; i++) {
					box_array[i] = V(0, i);
					box_array[i + DerivedV::ColsAtCompileTime] = V(0, i);
				}
				for (size_t i = 1; i < V.rows(); i++) {
					for (size_t k = 0; k < DerivedV::ColsAtCompileTime; k++) {
						box_array[k] = std::min(box_array[k], V(i, k));
						box_array[k + DerivedV::ColsAtCompileTime] =
							std::max(box_array[k + DerivedV::ColsAtCompileTime], V(i, k));
					}
				}

				init(box_array);
				size_t index = std::numeric_limits<size_t>::max();
				for (int i = 0; i < F.rows(); i++) {
					DerivedV M;
					for (size_t k = 0; k < DerivedV::RowsAtCompileTime; k++) {
						M.row(k) = V.row(F(i, k));
					}

					index = insert(M, i, index);
				}
			}

			explicit BinaryTree(const DerivedV& box) {
				std::array<double, DerivedV::ColsAtCompileTime * 2> initial_box;
				for (int i = 0; i < DerivedV::ColsAtCompileTime; i++) {
					initial_box[i] = box(0, i);
					initial_box[i + DerivedV::ColsAtCompileTime] = box(1, i);
				}
				tree_node_array_.emplace_back(0, std::numeric_limits<size_t>::max(),
					initial_box);
				generateChildHash();
				for (int i = 0; i < DIM; i++) {
					diag_len_of_box_ += (box[i] - box[DIM + i]) * (box[i] - box[DIM + i]);
				}
				diag_len_of_box_ = std::sqrt(diag_len_of_box_);
			}

			explicit BinaryTree(BoxVector& v1, BoxVector& v2) {
				std::array<double, DerivedV::ColsAtCompileTime * 2> initial_box;
				for (int i = 0; i < DerivedV::ColsAtCompileTime; i++) {
					initial_box[i] = v1(i);
					initial_box[i + DerivedV::ColsAtCompileTime] = v2(i);
				}
				tree_node_array_.emplace_back(0, std::numeric_limits<size_t>::max(),
					initial_box);
				generateChildHash();
				for (int i = 0; i < DIM; i++) {
					diag_len_of_box_ += (initial_box[i] - initial_box[DIM + i]) *
						(initial_box[i] - initial_box[DIM + i]);
				}
				diag_len_of_box_ = std::sqrt(diag_len_of_box_);
			}
			struct Data {
				Data() : simplex(DerivedV(DerivedV())), local_hint(0) {}

				explicit Data(const DerivedV& s)
					: simplex(s),
					local_hint(0) {  ///> modify std::double::max() -> 0 due to the most
					/// greatest common ancester is zero.
					for (int i = 0; i < DIM; i++) {
						box[i] = s(0, i);
						box[i + DIM] = s(0, i);
					}
					for (int j = 1; j < s.rows(); j++) {
						for (int i = 0; i < DIM; i++) {
							box[i] = std::min(box[i], s(j, i));
							box[i + DIM] = std::max(box[i + DIM], s(j, i));
						}
					}
				}
				DerivedV simplex;
				BoundingBox box;
				size_t local_hint;
				size_t id;
			};
			struct TreeNode {
				TreeNode(const size_t& current_id, const size_t& father_id,
					const BoundingBox& bbox)
					: id(current_id), fatherid(father_id), depth(0) {
					setbox(bbox);
				}
				size_t id;             ///< the index of current node
				std::size_t fatherid;  ///< the index of father node of current node
				std::vector<size_t>
					children;  ///> the index array of the children node within the tree
				std::vector<size_t>
					datas;  ///> the arrravector_to_matrixy that store the indexes of data
				BoundingBox box;
				BoundingBox coordcentor;
				size_t depth;
				inline void setbox(const BoundingBox& bbox) {
					box = bbox;
					for (int i = 0; i < DIM; i++) {
						coordcentor[i] = (box[i] + box[DIM + i]) * 0.5;
					}
				}
				inline bool isLeafNode() const { return children.empty(); }
			};

			void init(const Vertices& V, const Topos& F) {
				std::array<double, 2 * DerivedV::ColsAtCompileTime> box_array;
				for (size_t i = 0; i < DerivedV::ColsAtCompileTime; i++) {
					box_array[i] = V(0, i);
					box_array[i + DerivedV::ColsAtCompileTime] = V(0, i);
				}
				for (size_t i = 1; i < V.rows(); i++) {
					for (size_t k = 0; k < DerivedV::ColsAtCompileTime; k++) {
						box_array[k] = std::min(box_array[k], V(i, k));
						box_array[k + DerivedV::ColsAtCompileTime] =
							std::max(box_array[k + DerivedV::ColsAtCompileTime], V(i, k));
					}
				}

				init(box_array);
				size_t index = std::numeric_limits<size_t>::max();
				for (int i = 0; i < F.rows(); i++) {
					DerivedV M;
					for (size_t k = 0; k < DerivedV::RowsAtCompileTime; k++) {
						M.row(k) = V.row(F(i, k));
					}

					index = insert(M, i, index);
				}
			}

			void init(const BoundingBox& box) {
				if (tree_node_array_.size()) {
					throw std::runtime_error("The binary tree have be initialized.");
					return;
				}
				tree_node_array_.emplace_back(0, std::numeric_limits<size_t>::max(), box);
				generateChildHash();
				diag_len_of_box_ = 0;
				for (int i = 0; i < DIM; i++) {
					diag_len_of_box_ += (box[i] - box[DIM + i]) * (box[i] - box[DIM + i]);
				}
				diag_len_of_box_ = std::sqrt(diag_len_of_box_);
			}

			/**
			 * @brief
			 *
			 * @param M means the simplex n * dim, n is the number of point.
			 *
			 * @return the index of the inserted element
			 * @note Assuming that M must intersect with the box of node[id]
			 */
			size_t insert(const DerivedV& M, const size_t& simplex_id,
				const size_t& hint_id = std::numeric_limits<size_t>::max()) {
				Data data(M);
				data.id = simplex_id;
				size_t index = tree_data_array_.size();
				tree_data_array_.push_back(data);
				if (tree_data_array_.size() <= hint_id) {
					insert(simplex_id, 0);
				}
				else {
					size_t id = tree_data_array_[hint_id].local_hint;
					while (id != 0 && !beFullyIn(tree_node_array_[id].box, data.box)) {
						id = tree_node_array_[id].fatherid;
					}
					insert(simplex_id, id);
				}
				return index;
			}

			/**
			 * @brief
			 *
			 * @param M means the
			 *
			 * @return the index of the inserted element
			 * @note Assuming that M must intersect with the box of node[id]
			 */
			size_t reinsert(const DerivedV& M, const size_t& simplex_id,
				const size_t& hint_id = std::numeric_limits<size_t>::max()) {
				Data data(M);
				data.id = simplex_id;
				size_t index = simplex_id;
				tree_data_array_[simplex_id] = (data);
				if (tree_data_array_.size() <= hint_id) {
					insert(simplex_id, 0);
				}
				else {
					size_t id = tree_data_array_[hint_id].local_hint;
					while (id != 0 && !beFullyIn(tree_node_array_[id].box, data.box)) {
						id = tree_node_array_[id].fatherid;
					}
					insert(simplex_id, id);
				}
				return simplex_id;
			}

			void remove(const size_t& simplex_id,
				const size_t& hint_id = std::numeric_limits<size_t>::max()) {
				if (tree_data_array_.size() <= hint_id) {
					removeRecursion(simplex_id, 0);
				}
				else {
					size_t id = tree_data_array_[hint_id].local_hint;
					while (!beFullyIn(tree_node_array_[id].box,
						tree_data_array_[simplex_id].box)) {
						id = tree_node_array_[id].fatherid;
					}
					removeRecursion(simplex_id, id);
				}
				//  tree_data_array_.erase(simplex_id); //todo set deleted
			}

			void query(const BoundingBox& box, std::vector<DerivedV>& result,
				const size_t& hint_id = std::numeric_limits<size_t>::max()) const {
				std::vector<size_t> simplex_id_result;
				query(box, simplex_id_result, hint_id);
				if (DerivedV::RowsAtCompileTime > 1) {
					std::sort(simplex_id_result.begin(), simplex_id_result.end());
					auto index_it =
						std::unique(simplex_id_result.begin(), simplex_id_result.end());
					simplex_id_result.erase(index_it, simplex_id_result.end());
				}
				result.reserve(simplex_id_result.size());
				for (auto& id : simplex_id_result) {
					result.push_back(tree_data_array_[id].simplex);
				}
			}

			void query(const BoxVector& query_xyz, const double& radius,
				std::vector<size_t>& result,
				const size_t& hint_id = std::numeric_limits<size_t>::max()) const {
				BoundingBox box;
				for (int i = 0; i < DerivedV::ColsAtCompileTime; ++i) {
					box[i] = query_xyz[i] - radius;
					box[i + DerivedV::ColsAtCompileTime] = query_xyz[i] + radius;
				}
				result.clear();
				query(box, result, hint_id);
				if (DerivedV::RowsAtCompileTime > 1) {
					std::sort(result.begin(), result.end());
					auto index_it = std::unique(result.begin(), result.end());
					result.erase(index_it, result.end());
				}
			}

			void query(const BoundingBox& box, std::vector<size_t>& result,
				const size_t& hint_id = std::numeric_limits<size_t>::max()) const {
				if (tree_data_array_.empty()) return;
				if (tree_data_array_.size() <= hint_id) {
					queryRecursion(box, 0, result);
				}
				else {
					size_t id = tree_data_array_[hint_id].local_hint;
					while (id != 0 && !beFullyIn(tree_node_array_[id].box, box)) {
						id = tree_node_array_[id].fatherid;
					}
					queryRecursion(box, id, result);
				}
				if (DerivedV::RowsAtCompileTime > 1) {
					/*if (result.size() > 10000) {
					  std::cout << result.size() << std::endl;
					  std::sort(result.begin(), result.end());
					  auto index_it = std::unique(result.begin(), result.end());
					  result.erase(index_it, result.end());
					  while (1) {
						result.clear();
						if (tree_data_array_.size() <= hint_id) {
						  queryRecursion(box, 0, result);
						} else {
						  size_t id = tree_data_array_[hint_id].local_hint;
						  while (id != 0 && !beFullyIn(tree_node_array_[id].box, box)) {
							id = tree_node_array_[id].fatherid;
						  }
						  queryRecursion(box, id, result);
						}
					  }
					}*/
					std::sort(result.begin(), result.end());
					auto index_it = std::unique(result.begin(), result.end());
					result.erase(index_it, result.end());
				}
			}

			size_t projection(
				const BoxVector& xyz, BoxVector& ans_xyz,
				const double reference_box_len = -1,
				const size_t& hint_id = std::numeric_limits<size_t>::max()) {
				size_t simplex_id;
				double distance;
				Eigen::Vector3d
					barycentric_coordinates;  /// TODO : Not Vector3d maybe Vector2d.
				return queryNearestTriangle(xyz, simplex_id, distance, ans_xyz,
					barycentric_coordinates, reference_box_len,
					hint_id);
			}

			/// return hint id
			size_t queryNearestTriangle(
				const BoxVector& xyz, size_t& simplex_id, double& distance,
				BoxVector& closet_point, Eigen::Vector3d& barycentric_coordinates,
				const double reference_box_len = -1,
				const size_t& hint_id = std::numeric_limits<size_t>::max()) {
				if (tree_data_array_.empty() || tree_node_array_.empty()) {
					throw std::runtime_error("Tree is empty!");
					return 0;
				}

				distance = std::numeric_limits<double>::max();
				double query_length = reference_box_len;
				if (query_length <= 0) {
					query_length = diag_len_of_box_ * 1e-4;
				};

				assert(query_length > 0);

				std::vector<size_t> result;
				while (distance > 2 * query_length) {
					query_length *= 2;
					std::array<double, 6> box;
					for (int i = 0; i < 3; i++) {
						box[i] = xyz(0, i) - query_length;
						box[3 + i] = xyz(0, i) + query_length;
					}
					result.clear();
					this->query(box, result, hint_id);
					for (auto tsimplex_id : result) {
						double dis;
						BoxVector tcloset_point;
						Eigen::Vector3d tbarycentric_coordinates;
						Eigen::Matrix<double, 1, 3> P;
						P(0, 0) = xyz.x();
						P(0, 1) = xyz.y();
						P(0, 2) = xyz.z();

						this->point_simplex_distance(P, tree_data_array_[tsimplex_id].simplex,
							dis, tcloset_point,
							tbarycentric_coordinates);
						if (distance > dis) {
							distance = dis;
							barycentric_coordinates = tbarycentric_coordinates;
							simplex_id = tsimplex_id;
							closet_point = tcloset_point;
						}
						if (distance == 0) break;
					}
					if (result.size() > 0) break;
				}
				return simplex_id;
			}

			size_t queryNearestEdge(
				const BoxVector& xyz, size_t& simplex_id, double& distance,
				BoxVector& closet_point, Eigen::Vector2d& barycentric_coordinates,
				const double reference_box_len = -1,
				const size_t& hint_id = std::numeric_limits<size_t>::max()) {
				static_assert(DerivedV::RowsAtCompileTime == 2,
					"queryNearestEdge requires an edge simplex with two vertices.");
				if (tree_data_array_.empty() || tree_node_array_.empty()) {
					throw std::runtime_error("Tree is empty!");
					return 0;
				}

				distance = std::numeric_limits<double>::max();
				double query_length = reference_box_len;
				if (query_length <= 0) {
					query_length = diag_len_of_box_ * 1e-4;
				};

				assert(query_length > 0);

				std::vector<size_t> result;
				while (distance > 2 * query_length) {
					query_length *= 2;
					std::array<double, 6> box;
					for (int i = 0; i < 3; i++) {
						box[i] = xyz(0, i) - query_length;
						box[3 + i] = xyz(0, i) + query_length;
					}
					result.clear();
					this->query(box, result, hint_id);
					for (auto tsimplex_id : result) {
						double dis;
						BoxVector tcloset_point;
						Eigen::Vector2d tbarycentric_coordinates;
						Eigen::Matrix<double, 1, 3> P;
						P(0, 0) = xyz.x();
						P(0, 1) = xyz.y();
						P(0, 2) = xyz.z();

						this->point_simplex_distance(P, tree_data_array_[tsimplex_id].simplex,
							dis, tcloset_point,
							tbarycentric_coordinates);
						if (distance > dis) {
							distance = dis;
							barycentric_coordinates = tbarycentric_coordinates;
							simplex_id = tsimplex_id;
							closet_point = tcloset_point;
						}
						if (distance == 0) break;
					}
					if (result.size() > 0) break;
				}
				return simplex_id;
			}

			void test() {
				for (auto i : tree_node_array_) {
					if (i.isLeafNode()) {
						std::cout << "id=" << i.id << " father id=" << i.fatherid
							<< " size=" << i.datas.size() << " " << std::endl;
					}
				}
				std::cout << "hello world!" << std::endl;
			}

			std::vector<std::vector<DerivedV> > getTriInstance() {
				std::vector<std::vector<DerivedV> > ans;
				for (int i = 0; i < tree_node_array_.size(); i++) {
					if (tree_node_array_[i].isLeafNode() &&
						tree_node_array_[i].datas.size()) {
						ans.push_back(std::vector<DerivedV>());
						for (int j = 0; j < tree_node_array_[i].datas.size(); j++) {
							ans.back().push_back(
								tree_data_array_[tree_node_array_[i].datas[j]].simplex);
						}
					}
				}
				return ans;
			}

			std::vector<std::vector<std::pair<size_t, BoundingBox> > > getInstance() {
				std::vector<std::vector<std::pair<size_t, BoundingBox> > > ans;
				for (int i = 0; i < tree_node_array_.size(); i++) {
					if (tree_node_array_[i].isLeafNode() &&
						tree_node_array_[i].datas.size()) {
						ans.push_back(std::vector<std::pair<size_t, BoundingBox> >());
						for (int j = 0; j < tree_node_array_[i].datas.size(); j++) {
							ans.back().push_back(std::make_pair(
								tree_node_array_[i].datas[j],
								tree_data_array_[tree_node_array_[i].datas[j]].box));
						}
					}
				}
				return ans;
			}

			template <typename Derivedp, typename DerivedPV, typename DerivedEle,
				typename Derivedsqr_d, typename Derivedc, typename Derivedb>
			void point_simplex_distance(const Eigen::MatrixBase<Derivedp>& p,
				const std::vector<DerivedPV>& V,
				const std::vector<DerivedEle>& Ele,
				const typename DerivedEle::Index primitive,
				Derivedsqr_d& sqr_d,
				Eigen::MatrixBase<Derivedc>& c,
				Eigen::PlainObjectBase<Derivedb>& bary) {
				typedef typename Derivedp::Scalar Scalar;
				typedef typename Eigen::Matrix<Scalar, 1, Derivedp::ColsAtCompileTime>
					Vector;
				typedef Vector Point;
				// typedef Derivedb BaryPoint;
				typedef Eigen::Matrix<typename Derivedb::Scalar, 1, 3> BaryPoint;

				const auto& Dot = [](const Point& a, const Point& b) -> Scalar {
					return a.dot(b);
					};
				// Real-time collision detection, Ericson, Chapter 5
				const auto& ClosestBaryPtPointTriangle =
					[&Dot](Point p, Point a, Point b, Point c,
						BaryPoint& bary_out) -> Point {
							// Check if P in vertex region outside A
							Vector ab = b - a;
							Vector ac = c - a;
							Vector ap = p - a;
							Scalar d1 = Dot(ab, ap);
							Scalar d2 = Dot(ac, ap);
							if (d1 <= 0.0 && d2 <= 0.0) {
								// barycentric coordinates (1,0,0)
								bary_out << 1, 0, 0;
								return a;
							}
							// Check if P in vertex region outside B
							Vector bp = p - b;
							Scalar d3 = Dot(ab, bp);
							Scalar d4 = Dot(ac, bp);
							if (d3 >= 0.0 && d4 <= d3) {
								// barycentric coordinates (0,1,0)
								bary_out << 0, 1, 0;
								return b;
							}
							// Check if P in edge region of AB, if so return projection of P onto AB
							Scalar vc = d1 * d4 - d3 * d2;
							if (a != b) {
								if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
									Scalar v = d1 / (d1 - d3);
									// barycentric coordinates (1-v,v,0)
									bary_out << 1 - v, v, 0;
									return a + v * ab;
								}
							}
							// Check if P in vertex region outside C
							Vector cp = p - c;
							Scalar d5 = Dot(ab, cp);
							Scalar d6 = Dot(ac, cp);
							if (d6 >= 0.0 && d5 <= d6) {
								// barycentric coordinates (0,0,1)
								bary_out << 0, 0, 1;
								return c;
							}
							// Check if P in edge region of AC, if so return projection of P onto AC
							Scalar vb = d5 * d2 - d1 * d6;
							if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
								Scalar w = d2 / (d2 - d6);
								// barycentric coordinates (1-w,0,w)
								bary_out << 1 - w, 0, w;
								return a + w * ac;
							}
							// Check if P in edge region of BC, if so return projection of P onto BC
							Scalar va = d3 * d6 - d5 * d4;
							if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
								Scalar w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
								// barycentric coordinates (0,1-w,w)
								bary_out << 0, 1 - w, w;
								return b + w * (c - b);
							}
							// P inside face region. Compute Q through its barycentric coordinates
							// (u,v,w)
							Scalar denom = 1.0 / (va + vb + vc);
							Scalar v = vb * denom;
							Scalar w = vc * denom;
							bary_out << 1.0 - v - w, v, w;
							return a + ab * v +
								ac * w;  // = u*a + v*b + w*c, u = va * denom = 1.0-v-w
					};

				assert(Ele[0].cols() <= 3 &&
					"Only simplices up to triangles are considered");

				assert((Derivedb::RowsAtCompileTime == 1 ||
					Derivedb::ColsAtCompileTime == 1) &&
					"bary must be Eigen Vector or Eigen RowVector");
				assert(((Derivedb::RowsAtCompileTime == -1 ||
					Derivedb::ColsAtCompileTime == -1) ||
					(Derivedb::RowsAtCompileTime == Ele[0].rows() ||
						Derivedb::ColsAtCompileTime == Ele[0].rows())) &&
					"bary must be Dynamic or size of Ele.cols()");

				BaryPoint tmp_bary;
				c = (const Derivedc&)ClosestBaryPtPointTriangle(
					p, V[Ele[primitive](0)],
					// modulo is a HACK to handle points, segments and triangles. Because of
					// this, we need 3d buffer for bary
					V[Ele[primitive](1 % Ele[0].cols())],
					V[Ele[primitive](2 % Ele[0].cols())], tmp_bary);
				bary.resize(Derivedb::RowsAtCompileTime == 1 ? 1 : Ele[0].cols(),
					Derivedb::ColsAtCompileTime == 1 ? 1 : Ele[0].cols());
				bary.head(Ele[0].cols()) = tmp_bary.head(Ele[0].cols());
				sqr_d = (p - c).norm();
			}

			template <typename Derivedp, typename DerivedPV, typename Derivedsqr_d,
				typename Derivedc, typename Derivedb>
			void point_simplex_distance(const Eigen::MatrixBase<Derivedp>& p,
				const DerivedPV& V, Derivedsqr_d& sqr_d,
				Eigen::MatrixBase<Derivedc>& c,
				Eigen::PlainObjectBase<Derivedb>& bary) {
				typedef typename Derivedp::Scalar Scalar;
				typedef typename Eigen::Matrix<Scalar, 1, Derivedp::ColsAtCompileTime>
					Vector;
				typedef Vector Point;
				// typedef Derivedb BaryPoint;
				typedef Eigen::Matrix<typename Derivedb::Scalar, 1, 3> BaryPoint;

				const auto& Dot = [](const Point& a, const Point& b) -> Scalar {
					return a.dot(b);
					};
				// Real-time collision detection, Ericson, Chapter 5
				const auto& ClosestBaryPtPointTriangle =
					[&Dot](const Point p, const Point a, const Point b, const Point c,
						BaryPoint& bary_out) -> Point {
							// Check if P in vertex region outside A
							Vector ab = b - a;
							Vector ac = c - a;
							Vector ap = p - a;
							Scalar d1 = Dot(ab, ap);
							Scalar d2 = Dot(ac, ap);
							if (d1 <= 0.0 && d2 <= 0.0) {
								// barycentric coordinates (1,0,0)
								bary_out << 1, 0, 0;
								return a;
							}
							// Check if P in vertex region outside B
							Vector bp = p - b;
							Scalar d3 = Dot(ab, bp);
							Scalar d4 = Dot(ac, bp);
							if (d3 >= 0.0 && d4 <= d3) {
								// barycentric coordinates (0,1,0)
								bary_out << 0, 1, 0;
								return b;
							}
							// Check if P in edge region of AB, if so return projection of P onto AB
							Scalar vc = d1 * d4 - d3 * d2;
							if (a != b) {
								if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
									Scalar v = d1 / (d1 - d3);
									// barycentric coordinates (1-v,v,0)
									bary_out << 1 - v, v, 0;
									return a + v * ab;
								}
							}
							// Check if P in vertex region outside C
							Vector cp = p - c;
							Scalar d5 = Dot(ab, cp);
							Scalar d6 = Dot(ac, cp);
							if (d6 >= 0.0 && d5 <= d6) {
								// barycentric coordinates (0,0,1)
								bary_out << 0, 0, 1;
								return c;
							}
							// Check if P in edge region of AC, if so return projection of P onto AC
							Scalar vb = d5 * d2 - d1 * d6;
							if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
								Scalar w = d2 / (d2 - d6);
								// barycentric coordinates (1-w,0,w)
								bary_out << 1 - w, 0, w;
								return a + w * ac;
							}
							// Check if P in edge region of BC, if so return projection of P onto BC
							Scalar va = d3 * d6 - d5 * d4;
							if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
								Scalar w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
								// barycentric coordinates (0,1-w,w)
								bary_out << 0, 1 - w, w;
								return b + w * (c - b);
							}
							// P inside face region. Compute Q through its barycentric coordinates
							// (u,v,w)
							Scalar denom = 1.0 / (va + vb + vc);
							Scalar v = vb * denom;
							Scalar w = vc * denom;
							bary_out << 1.0 - v - w, v, w;
							return a + ab * v +
								ac * w;  // = u*a + v*b + w*c, u = va * denom = 1.0-v-w
					};

				assert(DerivedV::RowsAtCompileTime <= 3 &&
					"Only simplices up to triangles are considered");

				assert((Derivedb::RowsAtCompileTime == 1 ||
					Derivedb::ColsAtCompileTime == 1) &&
					"bary must be Eigen Vector or Eigen RowVector");
				assert(((Derivedb::RowsAtCompileTime == -1 ||
					Derivedb::ColsAtCompileTime == -1) ||
					(Derivedb::RowsAtCompileTime == DerivedV::RowsAtCompileTime ||
						Derivedb::ColsAtCompileTime == DerivedV::RowsAtCompileTime)) &&
					"bary must be Dynamic or size of Ele.cols()");

				BaryPoint tmp_bary;
				c = (const Derivedc&)ClosestBaryPtPointTriangle(
					p, V.row(0),
					// modulo is a HACK to handle points, segments and triangles. Because of
					// this, we need 3d buffer for bary
					V.row(1 % DerivedV::RowsAtCompileTime),
					V.row(2 % DerivedV::RowsAtCompileTime), tmp_bary);
				bary.resize(
					Derivedb::RowsAtCompileTime == 1 ? 1 : DerivedV::RowsAtCompileTime,
					Derivedb::ColsAtCompileTime == 1 ? 1 : DerivedV::RowsAtCompileTime);
				bary.head(DerivedV::RowsAtCompileTime) =
					tmp_bary.head(DerivedV::RowsAtCompileTime);
				sqr_d = (p - c).norm();
			}

		protected:
			inline size_t getHashNum(const Data& data,
				const BoundingBox& box_centor) const {
				return getHashNum(data.box, box_centor);
			}
			inline size_t getHashNum(const BoundingBox& box,
				const BoundingBox& box_centor) const {
				size_t hash_num = 0;
				size_t hash_id = 1;
				for (int i = 0; i < DIM; i++) {
					if (box[i] < box_centor[i]) {
						hash_num |= hash_id;
					}
					hash_id <<= 1;
					if (box[i + DIM] >= box_centor[i]) {
						hash_num |= hash_id;
					}
					hash_id <<= 1;
				}
				return hash_num;
			}
			void queryRecursion(const BoundingBox& box, const size_t& node_id,
				std::vector<size_t>& result) const {
				auto& coord_center = tree_node_array_[node_id].coordcentor;
				if (tree_node_array_[node_id].isLeafNode()) {
					for (auto i : tree_node_array_[node_id].datas) {
						auto& currunt_box = tree_data_array_[i].box;
						size_t j = 0;

						for (; j < DerivedV::ColsAtCompileTime; j++) {
							if (currunt_box[j] > box[j + DIM] || currunt_box[j + DIM] < box[j])
								break;
						}
						if (j == DerivedV::ColsAtCompileTime) result.push_back(i);
					}
				}
				else {
					for (int i = 0; i < hash_array.size(); i++) {
						auto data_hash = getHashNum(box, coord_center);

						if (!((~data_hash) & hash_array[i])) {  // if intersected
							queryRecursion(box, tree_node_array_[node_id].children[i], result);
						}
					}
				}
			}
			void removeRecursion(const size_t& simplex_id, const size_t& node_id) {
				Data& data = (tree_data_array_[simplex_id]);
				auto& coord_center = tree_node_array_[node_id].coordcentor;
				if (tree_node_array_[node_id].isLeafNode()) {
					for (auto it = tree_node_array_[node_id].datas.begin();
						it != tree_node_array_[node_id].datas.end(); it++) {
						if (*it == simplex_id) {
							tree_node_array_[node_id].datas.erase(it);
							break;
						}
					}
				}
				else {
					for (int i = 0; i < hash_array.size(); i++) {
						auto data_hash = getHashNum(data, coord_center);

						if (!((~data_hash) & hash_array[i])) {  // if intersected
							removeRecursion(simplex_id, tree_node_array_[node_id].children[i]);
						}
					}
				}
			}

		protected:
			void insert(const size_t& simplex_id, const size_t& node_id) {
				if (tree_node_array_[node_id].isLeafNode()) {
					tree_node_array_[node_id].datas.push_back(simplex_id);

					if (tree_node_array_[node_id].datas.size() > kMaxAllowedNodeDataSize &&
						tree_node_array_[node_id].depth < kMaxAllowedTreeHeight &&
						tree_node_array_.size() <
						tree_data_array_.size() / 2) {  // if decomposition is needed
						size_t child_id = tree_node_array_.size();
						std::vector<size_t> data_hashs;
						for (auto j = 0; j < tree_node_array_[node_id].datas.size(); j++) {
							data_hashs.push_back(
								getHashNum(tree_data_array_[tree_node_array_[node_id].datas[j]],
									tree_node_array_[node_id].coordcentor));
						}
						for (int i = 0; i < hash_array.size(); i++) {
							tree_node_array_[node_id].children.push_back(child_id + i);
							tree_node_array_.push_back(
								TreeNode(child_id + i, node_id, BoundingBox()));  ///> add box
							tree_node_array_.back().depth = tree_node_array_[node_id].depth + 1;
						}
						for (int i = 0; i < hash_array.size(); i++) {
							BoundingBox box;
							size_t hash_num = 1;
							for (int j = 0; j < DIM; j++) {
								if (hash_array[i] & hash_num) {  ///> construct the box for children
									box[j] = tree_node_array_[node_id].box[j];
									box[j + DIM] = tree_node_array_[node_id].coordcentor[j];
								}
								else {
									box[j] = tree_node_array_[node_id].coordcentor[j];
									box[j + DIM] = tree_node_array_[node_id].box[j + DIM];
								}
								hash_num <<= 2;
							}
							///> insert the child node into array
							tree_node_array_[child_id + i].setbox(box);  ///> add box
							///> insert element

							for (auto j = 0; j < tree_node_array_[node_id].datas.size(); j++) {
								if (!((~data_hashs[j]) & hash_array[i])) {
									insert(tree_node_array_[node_id].datas[j], child_id + i);
								}
							}
						}
						tree_node_array_[node_id].datas.clear();
					}
				}
				else {  ///> if the node has children
					auto data_hash = getHashNum(tree_data_array_[simplex_id].box,
						tree_node_array_[node_id].coordcentor);
					for (int i = 0; i < hash_array.size(); i++) {
						if (!((~data_hash) & hash_array[i])) {  // if intersected
							if (beFullyIn(
								tree_node_array_[tree_node_array_[node_id].children[i]].box,
								tree_data_array_[simplex_id].box)) {
								tree_data_array_[simplex_id].local_hint =
									tree_node_array_[node_id].children[i];
							}
							insert(simplex_id, tree_node_array_[node_id].children[i]);
						}
					}
				}
			}
			/**
			 * @brief caculate the hash_array for all octants. For instance, in
			 * in Quadtree, 0101 means the southwestern octants, while 1010 means the
			 * northeastern octant. in Octree, 010101 means the lower southeastern
			 * octants. the example of one octant decomposition scheme:
			 * ---------------------
			 * |  1001   |   1010  |
			 * |-------------------|
			 * |  0101   |   0110  |
			 * ---------------------
			 * @note this function only called in constructive function for once
			 */
		public:
			void generateChildHash() {
				hash_array.fill(0);
				for (int i = 0; i < DerivedV::ColsAtCompileTime; i++) {
					size_t current_size = 1 << i;
					for (int j = 0; j < current_size; j++) {
						hash_array[j + current_size] = hash_array[j] | (1 << (2 * i + 1));
						hash_array[j] = hash_array[j] | (1 << (2 * i));
					}
				}
			}
			/**
			 * @ brief whether box1 contain the box2
			 */
			inline bool beFullyIn(const BoundingBox& box1,
				const BoundingBox& box2) const {
				for (size_t i = 0; i < DerivedV::ColsAtCompileTime; i++) {
					if (box1[i] > box2[i] || box1[i + DIM] < box2[i + DIM]) return false;
				}
				return true;
			}

		protected:
			static const size_t kMaxAllowedTreeHeight = 15;
			static const size_t kMaxAllowedNodeDataSize = 80;

			HashArray hash_array;
			std::vector<TreeNode> tree_node_array_;
			std::vector<Data> tree_data_array_;
			double diag_len_of_box_ = 0.0;
		};
}  // namespace dt
#endif  // SRC_WYFDT_DT_BINARY_TREE_HPP_
		// #pragma optimize("", on)
