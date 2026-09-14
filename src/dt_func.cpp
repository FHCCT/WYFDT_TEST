#include "dt.h"
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include "spdlog/logger.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace {
	using FineMeshEdge = Eigen::Matrix<double, 2, 3>;
	using FineMeshEdgeTree = dt::BinaryTree<FineMeshEdge>;


	static Eigen::Vector3d closest_point_on_triangle(
		const Eigen::Vector3d& p,
		const Eigen::Vector3d& a,
		const Eigen::Vector3d& b,
		const Eigen::Vector3d& c
	) {
		const Eigen::Vector3d ab = b - a;
		const Eigen::Vector3d ac = c - a;
		const Eigen::Vector3d ap = p - a;
		const double d1 = ab.dot(ap);
		const double d2 = ac.dot(ap);
		if (d1 <= 0.0 && d2 <= 0.0) {
			return a;
		}

		const Eigen::Vector3d bp = p - b;
		const double d3 = ab.dot(bp);
		const double d4 = ac.dot(bp);
		if (d3 >= 0.0 && d4 <= d3) {
			return b;
		}

		const double vc = d1 * d4 - d3 * d2;
		if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
			const double v = d1 / (d1 - d3);
			return a + v * ab;
		}

		const Eigen::Vector3d cp = p - c;
		const double d5 = ab.dot(cp);
		const double d6 = ac.dot(cp);
		if (d6 >= 0.0 && d5 <= d6) {
			return c;
		}

		const double vb = d5 * d2 - d1 * d6;
		if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
			const double w = d2 / (d2 - d6);
			return a + w * ac;
		}

		const double va = d3 * d6 - d5 * d4;
		if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
			const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
			return b + w * (c - b);
		}

		const double denom = 1.0 / (va + vb + vc);
		const double v = vb * denom;
		const double w = vc * denom;
		return a + v * ab + w * ac;
	}


	static bool valid_vertex_index(int index, size_t vertex_count) {
		return index >= 0 && static_cast<size_t>(index) < vertex_count;
	}

	static bool valid_triangle(const dt::Mesh& mesh, const std::array<int, 3>& face) {
		const size_t vertex_count = mesh.V.size();
		if (!valid_vertex_index(face[0], vertex_count) ||
			!valid_vertex_index(face[1], vertex_count) ||
			!valid_vertex_index(face[2], vertex_count)) {
			return false;
		}

		const auto& va = mesh.V[face[0]];
		const auto& vb = mesh.V[face[1]];
		const auto& vc = mesh.V[face[2]];
		const Eigen::Vector3d a(va[0], va[1], va[2]);
		const Eigen::Vector3d b(vb[0], vb[1], vb[2]);
		const Eigen::Vector3d c(vc[0], vc[1], vc[2]);
		return std::isfinite(a.squaredNorm()) &&
			std::isfinite(b.squaredNorm()) &&
			std::isfinite(c.squaredNorm()) &&
			((b - a).cross(c - a)).squaredNorm() > 1e-30;
	}


	static std::array<int, 2> sorted_edge_key(int a, int b) {
		return a < b ? std::array<int, 2>{{ a, b }} : std::array<int, 2>{{ b, a }};
	}

	static void append_feature_edge(
		dt::FineMeshProjection& projection,
		const dt::Mesh& mesh,
		int a,
		int b,
		std::set<std::array<int, 2>>& seen_edges
	) {
		if (!valid_vertex_index(a, mesh.V.size()) ||
			!valid_vertex_index(b, mesh.V.size()) ||
			a == b) {
			return;
		}

		std::array<int, 2> key = sorted_edge_key(a, b);
		if (seen_edges.find(key) != seen_edges.end()) {
			return;
		}

		const auto& va = mesh.V[a];
		const auto& vb = mesh.V[b];
		Eigen::Vector3d p0(va[0], va[1], va[2]);
		Eigen::Vector3d p1(vb[0], vb[1], vb[2]);
		if ((p1 - p0).squaredNorm() <= 1e-30) {
			return;
		}

		FineMeshEdge edge;
		edge.row(0) = p0;
		edge.row(1) = p1;
		seen_edges.insert(key);
		projection.edges.push_back(edge);
	}

	static bool is_fine_mesh_crease_edge(
		const dt::Mesh& mesh,
		const std::array<int, 2>& edge,
		const std::array<int, 3>& face0,
		const std::array<int, 3>& face1,
		double angle
	) {
		int p1 = edge[0];
		int p2 = edge[1];
		int p3 = -1;
		int p4 = -1;

		for (int i = 0; i < 3; ++i) {
			if (face0[i] != p1 && face0[i] != p2) {
				p3 = face0[i];
				break;
			}
		}
		for (int i = 0; i < 3; ++i) {
			if (face1[i] != p1 && face1[i] != p2) {
				p4 = face1[i];
				break;
			}
		}
		if (p3 < 0 || p4 < 0) {
			return false;
		}

		const auto& v1 = mesh.V[p1];
		const auto& v2 = mesh.V[p2];
		const auto& v3 = mesh.V[p3];
		const auto& v4 = mesh.V[p4];
		Eigen::Vector3d a(v1[0], v1[1], v1[2]);
		Eigen::Vector3d b(v2[0], v2[1], v2[2]);
		Eigen::Vector3d c(v3[0], v3[1], v3[2]);
		Eigen::Vector3d d(v4[0], v4[1], v4[2]);
		Eigen::Vector3d n1 = (b - a).cross(c - a);
		Eigen::Vector3d n2 = (d - a).cross(b - a);
		const double n1_len = n1.norm();
		const double n2_len = n2.norm();
		if (n1_len <= 1e-30 || n2_len <= 1e-30) {
			return false;
		}

		double aa = n1.dot(n2) / (n1_len * n2_len);
		if (std::abs(aa) > 1.0) {
			return false;
		}
		return 180.0 - RADIO2ANGLE(std::acos(aa)) < angle;
	}

	static void extract_fine_mesh_feature_edges(dt::FineMeshProjection& projection, const dt::Mesh& mesh) {
		projection.edgesEnabled = false;
		projection.edges.clear();
		std::set<std::array<int, 2>> seen_edges;

		for (const auto& edge : mesh.S) {
			append_feature_edge(projection, mesh, edge[0], edge[1], seen_edges);
		}

		std::vector<std::array<int, 3>> faces;
		std::vector<int> face_geo;
		for (const auto& face_with_id : mesh.F) {
			std::array<int, 3> face = {{face_with_id[0], face_with_id[1], face_with_id[2]}};
			if (valid_triangle(mesh, face)) {
				faces.push_back(face);
				face_geo.push_back(face_with_id[3]);
			}
		}
		if (faces.empty()) return;

		struct EdgeFaceRecord {
			std::array<int, 3> face;
			int geo;
		};
		std::map<std::array<int, 2>, std::vector<EdgeFaceRecord>> edge_faces;
		for (size_t i = 0; i < faces.size(); ++i) {
			const auto& face = faces[i];
			EdgeFaceRecord record;
			record.face = face;
			record.geo = face_geo[i];
			edge_faces[sorted_edge_key(face[0], face[1])].push_back(record);
			edge_faces[sorted_edge_key(face[1], face[2])].push_back(record);
			edge_faces[sorted_edge_key(face[2], face[0])].push_back(record);
		}

		const double angle_check = 120.0;
		for (const auto& item : edge_faces) {
			const auto& adjacent_faces = item.second;
			bool is_feature = adjacent_faces.size() != 2;
			if (!is_feature && adjacent_faces[0].geo != adjacent_faces[1].geo) {
				is_feature = true;
			}
			if (!is_feature) {
				is_feature = is_fine_mesh_crease_edge(
					mesh,
					item.first,
					adjacent_faces[0].face,
					adjacent_faces[1].face,
					angle_check
				);
			}

			if (is_feature) {
				append_feature_edge(projection, mesh, item.first[0], item.first[1], seen_edges);
			}
		}
	}

	static bool build_fine_mesh_feature_edge_tree(dt::FineMeshProjection& projection) {
		projection.edgesEnabled = false;
		projection.edgeTree = FineMeshEdgeTree();
		if (projection.edges.empty()) {
			return false;
		}

		FineMeshEdgeTree::Vertices V(projection.edges.size() * 2, 3);
		FineMeshEdgeTree::Topos E(projection.edges.size(), 2);
		for (int i = 0; i < static_cast<int>(projection.edges.size()); ++i) {
			const int base = i * 2;
			V.row(base) = projection.edges[i].row(0);
			V.row(base + 1) = projection.edges[i].row(1);
			E(i, 0) = base;
			E(i, 1) = base + 1;
		}

		try {
			projection.edgeTree.init(V, E);
		}
		catch (...) {
			projection.edgeTree = FineMeshEdgeTree();
			return false;
		}

		projection.edgesEnabled = true;
		return true;
	}

	static Eigen::Vector3d project_boundary_point_to_fine_mesh_bruteforce(
		const dt::FineMeshProjection& projection,
		const Eigen::Vector3d& p,
		double max_projection_distance
	) {
		double best_distance = std::numeric_limits<double>::max();
		Eigen::Vector3d best_point = p;

		for (const auto& tri : projection.triangles) {
			const Eigen::Vector3d q = closest_point_on_triangle(
				p,
				tri.row(0).transpose(),
				tri.row(1).transpose(),
				tri.row(2).transpose()
			);
			const double distance = (q - p).norm();
			if (distance < best_distance) {
				best_distance = distance;
				best_point = q;
			}
		}

		return best_distance <= max_projection_distance ? best_point : p;
	}
}

bool DT::build_fine_mesh_projection_tree(const std::string& fine_mesh_file) {
	auto& projection = fineMeshProjection;
	projection.enabled = false;
	projection.triangles.clear();
	projection.tree = dt::BinaryTree<Eigen::Matrix3d>();
	projection.edgesEnabled = false;
	projection.edges.clear();
	projection.edgeTree = FineMeshEdgeTree();

	if (fine_mesh_file.empty()) {
		return false;
	}

	dt::Mesh fine_mesh;
	std::string filename = fine_mesh_file;
	if (!dt::readMesh(filename, fine_mesh)) {
		return false;
	}

	for (const auto& face_with_id : fine_mesh.F) {
		const std::array<int, 3> face = {{face_with_id[0], face_with_id[1], face_with_id[2]}};
		if (!valid_triangle(fine_mesh, face)) continue;
		Eigen::Matrix3d triangle;
		for (int i = 0; i < 3; ++i) {
			const auto& vertex = fine_mesh.V[face[i]];
			triangle.row(i) = Eigen::Vector3d(vertex[0], vertex[1], vertex[2]);
		}
		projection.triangles.push_back(triangle);
	}
	extract_fine_mesh_feature_edges(projection, fine_mesh);
	build_fine_mesh_feature_edge_tree(projection);
	if (projection.triangles.empty()) {
		return false;
	}

	Eigen::MatrixXd V(projection.triangles.size() * 3, 3);
	Eigen::MatrixXi F(projection.triangles.size(), 3);
	for (int i = 0; i < static_cast<int>(projection.triangles.size()); ++i) {
		const int base = i * 3;
		for (int j = 0; j < 3; ++j) {
			V.row(base + j) = projection.triangles[i].row(j);
			F(i, j) = base + j;
		}
	}

	try {
		projection.tree.init(V, F);
	}
	catch (...) {
		projection.tree = dt::BinaryTree<Eigen::Matrix3d>();
		return false;
	}

	projection.enabled = true;
	return true;
}

bool DT::project_boundary_point_to_fine_mesh(
	double* in,
	double* out,
	double max_projection_distance
) {
	auto& projection = fineMeshProjection;
	if (in == nullptr || out == nullptr) {
		return false;
	}

	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];

	if (!projection.enabled ||
		projection.triangles.empty() ||
		max_projection_distance < 0.0) {
		return false;
	}

	const Eigen::Vector3d p(in[0], in[1], in[2]);
	Eigen::Vector3d projected = p;

	try {
		Eigen::Matrix<double, 1, 3> query(p.x(), p.y(), p.z());
		Eigen::Matrix<double, 1, 3> closest;
		Eigen::Vector3d barycentric;
		size_t faceidx = 0;
		double distance = std::numeric_limits<double>::max();
		projection.tree.queryNearestTriangle(query, faceidx, distance, closest, barycentric);
		if (distance <= max_projection_distance) {
			projected = Eigen::Vector3d(closest(0, 0), closest(0, 1), closest(0, 2));
		}
	}
	catch (...) {
		projected = project_boundary_point_to_fine_mesh_bruteforce(projection, p, max_projection_distance);
	}

	out[0] = projected.x();
	out[1] = projected.y();
	out[2] = projected.z();
	
	return true;
}

bool DT::project_segment_point_to_fine_mesh(
	double* in,
	double* out,
	double max_projection_distance
) {
	auto& projection = fineMeshProjection;
	if (in == nullptr || out == nullptr) {
		return false;
	}

	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];

	if (!projection.edgesEnabled ||
		projection.edges.empty() ||
		max_projection_distance < 0.0) {
		return false;
	}

	try {
		Eigen::Matrix<double, 1, 3> query(in[0], in[1], in[2]);
		Eigen::Matrix<double, 1, 3> closest;
		Eigen::Vector2d barycentric;
		size_t edgeidx = 0;
		double distance = std::numeric_limits<double>::max();

		projection.edgeTree.queryNearestEdge(
			query,
			edgeidx,
			distance,
			closest,
			barycentric
		);

		if (distance > max_projection_distance) {
			return false;
		}

		out[0] = closest(0, 0);
		out[1] = closest(0, 1);
		out[2] = closest(0, 2);
		return true;
	}
	catch (...) {
		out[0] = in[0];
		out[1] = in[1];
		out[2] = in[2];
		return false;
	}
}
////////////////////////////////////////////////////////////////
/**************************  function *************************/
////////////////////////////////////////////////////////////////

int DT::addElem(bool UseVacancy) {
	int newE = -1;
	if (UseVacancy && !Evacancy.empty()) {
		newE = Evacancy.front();
		Evacancy.pop();
		Elems[newE].form[0] = 0;
		Elems[newE].form[1] = 0;
		Elems[newE].form[2] = 0;
		Elems[newE].form[3] = 0;
		Elems[newE].neig[0] = -1;
		Elems[newE].neig[1] = -1;
		Elems[newE].neig[2] = -1;
		Elems[newE].neig[3] = -1;
		Elems[newE].info = 0;
		Elems[newE].geo = -1;
		Elems[newE].q = -1;
	}
	else {
		newE = Elems.size();
		Elems.push_back(Elem());
	}
	return newE;
}

int DT::addNode(bool UseVacancy) {
	int newN = -1;
	if (UseVacancy && !Nvacancy.empty()) {
		newN = Nvacancy.front();
		Nvacancy.pop();
		Nodes[newN].pt[0] = 0;
		Nodes[newN].pt[1] = 0;
		Nodes[newN].pt[2] = 0;
		Nodes[newN].space = 0;
		Nodes[newN].info = 0;
		Nodes[newN].tet = 0;
		Nodes[newN].type = 0;
		Nodes[newN].occupying.store(-1);
	}
	else {
		newN = Nodes.size();
		Nodes.push_back(Node());
	}

	return newN;
}

int DT::addElem(int pa, int pb, int pc, int pd, bool UseVacancy) {
	int newE = -1;
	if (UseVacancy && !Evacancy.empty()) {
		newE = Evacancy.front();
		Evacancy.pop();
		Elems[newE].form[0] = pa;
		Elems[newE].form[1] = pb;
		Elems[newE].form[2] = pc;
		Elems[newE].form[3] = pd;
		Elems[newE].neig[0] = -1;
		Elems[newE].neig[1] = -1;
		Elems[newE].neig[2] = -1;
		Elems[newE].neig[3] = -1;
		Elems[newE].info = 0;
		Elems[newE].geo = -1;
		Elems[newE].q = -1;
	}
	else {
		newE = Elems.size();
		Elems.push_back(Elem(pa, pb, pc, pd));
	}
	return newE;
}

int DT::addNode(double x, double y, double z, double space, bool UseVacancy) {
	int newN = -1;
	if (UseVacancy && !Nvacancy.empty()) {
		newN = Nvacancy.front();
		Nvacancy.pop();
		Nodes[newN].pt[0] = x;
		Nodes[newN].pt[1] = y;
		Nodes[newN].pt[2] = z;
		Nodes[newN].space = space;
		Nodes[newN].info = 0;
		Nodes[newN].tet = 0;
		Nodes[newN].type = 0;
		Nodes[newN].occupying.store(-1);
	}
	else {
		newN = Nodes.size();
		Nodes.push_back(Node(x, y, z, space));
	}

	return newN;
}

//adjacent connect and decode
bool DT::bond(int t1, int t2) {
	int nig1 = -1, nig2 = -1;

	for (int i = 0; i < 4; ++i) {
		bool found = false;
		for (int j = 0; j < 4; ++j) {
			if (Elems[t1].form[i] == Elems[t2].form[j]) {
				found = true;
				break;
			}
		}
		if (!found) {
			if (nig1 == -1) nig1 = i;
			else return false;
		}
	}

	for (int j = 0; j < 4; ++j) {
		bool found = false;
		for (int i = 0; i < 4; ++i) {
			if (Elems[t2].form[j] == Elems[t1].form[i]) {
				found = true;
				break;
			}
		}
		if (!found) {
			if (nig2 == -1) nig2 = j;
			else return false;
		}
	}

	bond(t1, nig1, t2, nig2);
	return true;
}

void DT::bond(int t1, int nig1, int t2, int nig2) {
	if (t1 != -1 && t2 != -1) {
		Elems[t1].neig[nig1] = ((int64_t)t2 << 2) | nig2;
		Elems[t2].neig[nig2] = ((int64_t)t1 << 2) | nig1;
	}
	else if (t1 == -1) {
		Elems[t2].neig[nig2] = -1;
	}
	else if (t2 == -1) {
		Elems[t1].neig[nig1] = -1;
	}
	return;
}

/*
* During point insert ,it don't ues.and may destroy
* So should build before
* Mesh improvement and boundary recover
*/
void DT::setAllP2T() {
	for (int i = 0; i < Elems.size(); i++) {
		if (isDelEle(i))
			continue;
		if (isvirtualtet(i))
			continue;
		if (ishulltet(i))
			continue;
		//Elems[i].info = 0;
		for (int j = 0; j < 4; j++) {
			setP2T(Elems[i].form[j], i);
		}
	}
	return;
}

int DT::isNod_in_Tet(int iNod, int tet)
{
	if (tet<0 || tet>Elems.size()) 
		return-1;
	for (int m = 0; m < 4; m++)
	{
		if (Elems[tet].form[m] == iNod)
			return m;//id
	}
	return -1;//Node don't int tet
}

bool DT::isDelEle(int E) {
	return Elems[E].info < 0;
}

void DT::DelEle(int E) {
	if (isDelEle(E))return;
	Elems[E].form[0] = Elems[E].form[1] = Elems[E].form[2] = Elems[E].form[3] = -1;
	Elems[E].neig[0] = Elems[E].neig[1] = Elems[E].neig[2] = Elems[E].neig[3] = -1;
	Elems[E].info = -1; Elems[E].geo = -1, Elems[E].q = -1;
	Evacancy.push(E);//Collect delete id
	return;
}

bool DT::isDelNod(int N) {
	bool result;
#pragma omp critical
	{
		result = Nodes[N].info < 0;
	}
	return result;
}

void DT::DelNod(int N) {
	if (isDelNod(N))
		return;
#pragma omp critical
	{
		Nodes[N].pt[0] = Nodes[N].pt[1] = Nodes[N].pt[2] = 0;
		Nodes[N].space = 0;
		Nodes[N].tet = -1;
		Nodes[N].info = -1;
		Nodes[N].type = 0;
		Nvacancy.push(N);//Collect delete id
	}
	return;
}

bool DT::isDelSurTri(int i) {
	return SurTris[i].info == -100;
}
bool DT::isDelSurEdg(int i) {
	return SurEdgs[i].info == -100;
}
bool DT::setDelSurTri(int i) {
	return SurTris[i].info = -100;
}
bool DT::setDelSurEdg(int i) {
	return SurEdgs[i].info = -100;
}
bool DT::isRecBndEdg(int i) {
	return SurEdgs[i].info > 0;
}
bool DT::isRecBndTri(int i) {
	return SurTris[i].info > 0;
}

int DT::findTriParent(int i) {
	while (i >= nSurTris) {
		int parent = SurTris[i].parent;
		if (parent > SurTris.size())
			return parent;
		int parentinfo = SurTris[parent].info;
		if (parentinfo <= i && i <= parentinfo + 3) {
			i = parent;
		}
		else {
			return parent;
		}
	}
	return SurTris[i].parent;
}

bool DT::isBndEdg(const int p1, const int p2) {
	return BndEdg.find(p1, p2) != nullptr;
}
bool DT::isBndTri(const int p1, const int p2, const int p3) {
	return BndTri.find(p1, p2, p3) != nullptr;
}
int DT::getoppoP(int tet, int i) {
	return Elems[getNeig(tet, i)].form[getNeigOrd(tet, i)];
}
int DT::setP2T(int p, int t) {
	return Nodes[p].tet = t;
}
int DT::getP2T(int p) {
	return Nodes[p].tet;
}
bool DT::isbndpnt(int i) {
	return Nodes[i].type > 0;
}
void DT::setbndpnt(int i) {
	Nodes[i].type = 1;
	return;
}
void DT::clearbndpnt(int i) {
	Nodes[i].type = 0;
	return;
}
bool DT::isFacetpnt(int i) {
	return Nodes[i].type == 1;
}

bool DT::isSegmentpnt(int i) {
	return Nodes[i].type == 3;
}

bool DT::isCornerpnt(int i) {
	if (Nodes[i].type == 2 || Nodes[i].type > 3)
		return true;
	return false;
}

//Clean up empty bodies and points
//in this step ,we clean all info,and ignore everything
//so this step should used in last
void DT::clearNodesElems() {
	int i, j = 0;
	//Nodes info store node's new address,and info will clean
	for (i = 0; i < Nodes.size(); i++) {
		if (isDelNod(i))
			continue;
		Nodes[i].info = j++;
	}

	//reset Elems's form,don't reset Elems's neig,it's info,so after this step,neig is useless
	for (i = 0; i < Elems.size(); i++) {
		if (isDelEle(i)) //this step have no hull tet must
			continue;
		for (j = 0; j < 4; j++) {
			Elems[i].form[j] = Nodes[Elems[i].form[j]].info;
		}
	}
	//clean SurTri
	for (i = 0; i < SurTris.size(); i++) {
		if (isDelSurTri(i)) //this step have no hull tet must
			continue;
		for (j = 0; j < 3; j++) {
			if (Nodes[SurTris[i].form[j]].info == -1) {
				setDelSurTri(i);
				break;
			}
			SurTris[i].form[j] = Nodes[SurTris[i].form[j]].info;
		}
	}
	return;
}

void DT::getMeshEdgebyGeo(Mesh& mesh) {
	EdgeHasher<int> Edgindex;
	int te = 1;

	for (int i = 0; i < SurEdgs.size(); i++) {
		if (isDelSurEdg(i))
			continue;
		if (SurEdgs[i].face.size() != 2) {
			mesh.S.push_back({ Nodes[SurEdgs[i].iStart].info,Nodes[SurEdgs[i].iEnd].info,SurEdgs[i].geo });
		}
		else if (SurTris[SurEdgs[i].face[0]].parent != SurTris[SurEdgs[i].face[1]].parent) {
			mesh.S.push_back({ Nodes[SurEdgs[i].iStart].info,Nodes[SurEdgs[i].iEnd].info,SurEdgs[i].geo });
		}
	}
	return;
}
//--------------------- Mesh Quality ---------------------
double DT::tetquality(double v1[3], double v2[3], double v3[3], double v4[3], double* AniMetric, int qualmeasure) {
	if (qualmeasure == SUS_METRIC) return quality_sus(v1, v2, v3, v4);
	double quality = 0;
	if (qualmeasure == 1) {
		//efective in removing large dihedral angles.
		quality = vlrms3ratio(v1, v2, v4, v3);
	}
	else if (qualmeasure == 2) {
		//use minimum dihedral angle is available.
		double minangle, maxangle;
		std::vector<double> Dihedral;
		if (CalDihedral(v1, v2, v4, v3, minangle, maxangle, Dihedral)) {
			quality = std::sin(minangle);
		}
		else {
			quality = 0;//invert
		}
	}
	else if (qualmeasure == 3 || improve_Metric == 4) {
		//use Volume.
		quality = dt::GEOM_FUNC::orient3d(v1, v2, v4, v3) / 6.0;
	}
	else if (qualmeasure == 5) {
		quality = caltet33_ani(v1, v2, v3, v4, AniMetric);
	}
	else if (qualmeasure == 6) {
		quality = ScaledJacobian(v1, v2, v3, v4);
	}
	return quality;
}

//calculate the dihedral angle in radians by volume
/* compute the tangent of the angle using the tangent formula:

	tan(theta_ij) = - 6 * V * l_ij
					--------------
					dot(n_k, n_l)

	because this formula is accurate in the entire range.
*/
int  DT::CalDihedral(double v1[3], double v2[3], double v3[3], double v4[3],
	double& minangle, double& maxangle, std::vector<double>& Dihedral) {
	double point[4][3];      /* tet vertices */
	double edgelength[3][4]; /* the lengths of each of the edges of the tet */
	double facenormal[4][3]; /* the normals of each face of the tet */
	double dx, dy, dz;       /* intermediate values of edge lengths */
	double pyrvolume;        /* volume of tetrahedron */
	int i, j, k, l;          /* loop indices */
	double E[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
	double angle, tantheta;
	double dotproduct;

	minangle = DBL_MAX;
	maxangle = 0.0;
	Dihedral.resize(6);

	for (i = 0; i < 3; i++) {
		point[0][i] = v1[i];
		point[1][i] = v2[i];
		point[2][i] = v3[i];
		point[3][i] = v4[i];
	}

	/* calculate the volume*6 of the tetrahedron */
	pyrvolume = (double)dt::GEOM_FUNC::orient3d(point[0], point[1], point[2], point[3]);

	/* if the volume is zero, the quality is zero, no reason to continue */
	if (pyrvolume <= 0.0)
		return 0;

	/* for each vertex/face of the tetrahedron */
	for (i = 0; i < 4; i++) {
		j = (i + 1) & 3;
		if ((i & 1) == 0) {
			k = (i + 3) & 3;
			l = (i + 2) & 3;
		}
		else {
			k = (i + 2) & 3;
			l = (i + 3) & 3;
		}

		/* compute the normal for each face */
		facenormal[i][0] =
			(point[k][1] - point[j][1]) * (point[l][2] - point[j][2]) -
			(point[k][2] - point[j][2]) * (point[l][1] - point[j][1]);
		facenormal[i][1] =
			(point[k][2] - point[j][2]) * (point[l][0] - point[j][0]) -
			(point[k][0] - point[j][0]) * (point[l][2] - point[j][2]);
		facenormal[i][2] =
			(point[k][0] - point[j][0]) * (point[l][1] - point[j][1]) -
			(point[k][1] - point[j][1]) * (point[l][0] - point[j][0]);

		/* compute edge lengths (squared) */
		for (j = i + 1; j < 4; j++) {
			dx = point[i][0] - point[j][0];
			dy = point[i][1] - point[j][1];
			dz = point[i][2] - point[j][2];
			edgelength[i][j] = dx * dx + dy * dy + dz * dz;
		}
	}

	/* for each edge in the tetrahedron */
	for (i = 0; i < 3; i++) {
		for (j = i + 1; j < 4; j++) {
			k = (i > 0) ? 0 : (j > 1) ? 1 : 2;
			l = 6 - i - j - k;

			dotproduct = facenormal[k][0] * facenormal[l][0] + facenormal[k][1] * facenormal[l][1] + facenormal[k][2] * facenormal[l][2];

			if (dotproduct != 0.0) {
				tantheta = (-pyrvolume * sqrt(edgelength[i][j])) / dotproduct;
				/* now compute the actual angle */
				angle = atan(tantheta);
			}
			else angle = PI / 2.0;

			/* adjust angle for sign of dot product */
			if (dotproduct > 0) angle += PI;

			/* make negative angles positive */
			if (angle < 0) angle += 2.0 * PI;

			if (dotproduct == 0.0) angle = PI / 2.0;

			//Because will change 2 3,so idx will change
			int idx = -1;
			if (i == 0) {
				if (j == 1)idx = 0;
				else if (j == 2)idx = 2;
				else idx = 1;
			}
			else if (i == 1) {
				if (j == 2)idx = 4;
				else idx = 3;
			}
			else {
				idx = 5;
			}
			Dihedral[idx] = angle;
			if (angle < minangle) minangle = angle;
			if (angle > maxangle) maxangle = angle;
		}
	}
	return 1;
}

//volme edge ratio
double DT::vlrms3ratio(double v1[3], double v2[3], double v3[3], double v4[3]) {
	double point[4][3];			 /* tet vertices */
	// double edgelength[3][4];  /* the lengths of each of the edges of the tet */
	double volume;			 /* volume of tetrahedron */
	int i, j;				 /* loop indices */
	double edgelengthsum = 0.0;

	for (i = 0; i < 3; i++) {
		point[0][i] = v1[i];
		point[1][i] = v2[i];
		point[2][i] = v3[i];
		point[3][i] = v4[i];
	}
	/* calculate the volume*6 of the tetrahedron */
	volume = (double)dt::GEOM_FUNC::orient3d(v1, v2, v3, v4);

	/* for each edge in the tetrahedron */
	for (i = 0; i < 3; i++) {
		for (j = i + 1; j < 4; j++) {
			edgelengthsum += distance2(point[i], point[j]);
		}
	}
	/* compute the root mean square */
	const double SQ_2 = 12 * sqrt(3.0);
	return (SQ_2 * volume) / (edgelengthsum * sqrt(edgelengthsum));
}

// quality based on shortest height / longest edge
// 1: best, equilateral tetrahedron
// 0: worst, degenerate tetrahedron
double DT::AspectRatio(double v1[3], double v2[3], double v3[3], double v4[3]) {
	double point[4][3];

	for (int i = 0; i < 3; i++) {
		point[0][i] = v1[i];
		point[1][i] = v2[i];
		point[2][i] = v3[i];
		point[3][i] = v4[i];
	}

	// signed 6 * volume
	double volume6 = (double)dt::GEOM_FUNC::orient3d(v1, v2, v3, v4);

	if (volume6 <= 0.0) {
		return volume6;
	}

	// longest edge length squared
	double maxedgelen2 = 0.0;

	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 4; j++) {
			double len2 = distance2(point[i], point[j]);

			if (len2 > maxedgelen2) {
				maxedgelen2 = len2;
			}
		}
	}

	if (maxedgelen2 <= 0.0) {
		return 0.0;
	}

	double maxedgelen = sqrt(maxedgelen2);

	// shortest vertex-to-opposite-face height
	double minheight = DBL_MAX;

	for (int i = 0; i < 4; i++) {
		int ids[3];
		int cnt = 0;

		for (int j = 0; j < 4; j++) {
			if (j != i) {
				ids[cnt++] = j;
			}
		}

		double a[3], b[3], cross[3];

		for (int k = 0; k < 3; k++) {
			a[k] = point[ids[1]][k] - point[ids[0]][k];
			b[k] = point[ids[2]][k] - point[ids[0]][k];
		}

		cross[0] = a[1] * b[2] - a[2] * b[1];
		cross[1] = a[2] * b[0] - a[0] * b[2];
		cross[2] = a[0] * b[1] - a[1] * b[0];

		// |cross| = 2 * area of the opposite face
		double face2area = sqrt(
			cross[0] * cross[0] +
			cross[1] * cross[1] +
			cross[2] * cross[2]
		);

		if (face2area <= 0.0) {
			return 0.0;
		}

		// h = 6V / (2A)
		double height = volume6 / face2area;

		if (height < minheight) {
			minheight = height;
		}
	}

	if (minheight <= 0.0) {
		return 0.0;
	}

	// normalized quality:
	// equilateral tetrahedron gives 1
	// degenerate tetrahedron approaches 0
	double quality = sqrt(3.0 / 2.0) * minheight / maxedgelen;

	// clamp to [0, 1] to avoid tiny floating-point overshoot
	if (quality < 0.0) {
		quality = 0.0;
	}

	if (quality > 1.0) {
		quality = 1.0;
	}

	return quality;
}

int DT::getmm(int p1, int p2, int p3, int p4, double* mm) {
	int form[4] = { p1,p2,p3,p4 };
	int n = 0;

	for (int k = 0; k < 6; k++) mm[k] = 0;
	for (int j = 0; j < 4; j++) {
		//if (isSegmentpnt(form[j])) continue;
		//n++;
		for (int k = 0; k < 6; k++) {
			mm[k] += AniSol[form[j]][k];
		}
	}

	//if (n == 0)
	//	return 0;
	//double dd = 1.0 / n;
	//for (int k = 0; k < 6; k++) mm[k] *= dd;
	for (int k = 0; k < 6; k++) mm[k] /= 4.0;
	return 1;
}

double DT::ScaledJacobian(double v1[3],double v2[3],double v3[3],double v4[3], double eps) {
	// 当前程序中的正向四面体顺序：
	// orient3d(v1, v2, v4, v3) = 6 * signed volume
	const double jacobian =
		static_cast<double>(
			dt::GEOM_FUNC::orient3d(v1, v2, v4, v3)
			);

	if (jacobian == 0.0) {
		return 0.0;
	}

	// 六条边长
	const double l12 = distance(v1, v2);
	const double l13 = distance(v1, v3);
	const double l14 = distance(v1, v4);
	const double l23 = distance(v2, v3);
	const double l24 = distance(v2, v4);
	const double l34 = distance(v3, v4);

	// 四个顶点处三条关联边长度的乘积
	const double lambda1 = l12 * l13 * l14;
	const double lambda2 = l12 * l23 * l24;
	const double lambda3 = l13 * l23 * l34;
	const double lambda4 = l14 * l24 * l34;

	// Scaled Jacobian 取四个顶点归一化 Jacobian 的最小值。
	// 因为四个顶点的 Jacobian 绝对值相同，
	// 等价于除以最大的边长乘积。
	const double denominator = std::max(
		std::max(lambda1, lambda2),
		std::max(lambda3, lambda4)
	);

	if (denominator == 0.0) {
		return 0.0;
	}

	double quality = std::sqrt(2.0) * jacobian / denominator;

	// 防止浮点误差导致结果略微超出理论范围
	if (quality > 1.0) {
		quality = 1.0;
	}
	else if (quality < -1.0) {
		quality = -1.0;
	}

	return quality;
}

//double DT::ScaledJacobian(double v1[3], double v2[3],double v3[3], double v4[3],double eps){
//	const double o = dt::GEOM_FUNC::orient3d(v1, v2, v4, v3);
//	if (o == 0.0) 
//		return 0.0; // 退化
//	const double sgn = (o > 0.0) ? 1.0 : -1.0;
//
//	// 2) Jacobian（用标准三重积的“大小”，符号交给 orient）
//	double e12[3], e13[3], e14[3];
//	vsub3(v2, v1, e12);
//	vsub3(v3, v1, e13);
//	vsub3(v4, v1, e14);
//
//	double c[3];
//	cross(e12, e13, c);
//	const double Jabs = std::fabs(dot(c, e14)); // 只取大小
//	const double J = sgn * Jabs;             // 符号对齐 orient3d
//
//	if (J == 0)
//		return o;
//	// 3) VERDICT 分母：四个顶点处“三条 incident 边长乘积”的最大值
//	// vertex v1: |v2-v1| |v3-v1| |v4-v1|
//	// vertex v2: |v1-v2| |v3-v2| |v4-v2|
//	// vertex v3: |v1-v3| |v2-v3| |v4-v3|
//	// vertex v4: |v1-v4| |v2-v4| |v3-v4|
//	double a[3], b[3], c3[3];
//
//	vsub3(v2, v1, a); vsub3(v3, v1, b); vsub3(v4, v1, c3);
//	const double lam1 = lenvec(a) * lenvec(b) * lenvec(c3);
//
//	vsub3(v1, v2, a); vsub3(v3, v2, b); vsub3(v4, v2, c3);
//	const double lam2 = lenvec(a) * lenvec(b) * lenvec(c3);
//
//	vsub3(v1, v3, a); vsub3(v2, v3, b); vsub3(v4, v3, c3);
//	const double lam3 = lenvec(a) * lenvec(b) * lenvec(c3);
//
//	vsub3(v1, v4, a); vsub3(v2, v4, b); vsub3(v3, v4, c3);
//	const double lam4 = lenvec(a) * lenvec(b) * lenvec(c3);
//
//	const double denom = std::max(std::max(lam1, lam2), std::max(lam3, lam4));
//	if (denom <= eps) return 0.0;
//
//	// 4) scaled jacobian
//	double q = std::sqrt(2.0) * (J / denom);
//
//	// 5) 强制范围 [-1, 1]
//	if (q > 1.0) q = 1.0;
//	if (q < -1.0) q = -1.0;
//	return q;
//}

/**
 * Compute the quality of the tet pt with respect to the anisotropic metric \a
 * met. \f$ Q = V_met(K) / (sum(len(edge_K)^2)^(3/2) \f$ and for a calssic
 * storage of metrics at ridges.
 */
//#pragma optimize("",off)
double DT::caltet33_ani(double v1[3], double v2[3], double v3[3], double v4[3], double* mm) {
	double ab[3], ac[3], ad[3], bc[3], bd[3], cd[3];
	double  v1c, v2c, v3c, vol, det, h1, h2, h3, h4, h5, h6, rap, num, cal;
	double minV = 1e-7;

	// Edge vectors
	for (int i = 0; i < 3; i++) {
		ab[i] = v2[i] - v1[i];
		ac[i] = v3[i] - v1[i];
		ad[i] = v4[i] - v1[i];
		bc[i] = v3[i] - v2[i];
		bd[i] = v4[i] - v2[i];
		cd[i] = v4[i] - v3[i];
	}

	// Geometric volume
	vol = dt::GEOM_FUNC::orient3d(v1, v2, v4, v3);
	
	if (vol <= minV)
		return vol;

	// Determinant of metric matrix
	det = mm[0] * (mm[3] * mm[5] - mm[4] * mm[4])
		- mm[1] * (mm[1] * mm[5] - mm[2] * mm[4])
		+ mm[2] * (mm[1] * mm[4] - mm[2] * mm[3]);

	if (det <= 1e-14)
		return vol;

	det = sqrt(det) * vol;

	// Edge metric lengths squared
#define LEN_SQ(u) (mm[0]*u[0]*u[0] + mm[3]*u[1]*u[1] + mm[5]*u[2]*u[2] + 2.0*(mm[1]*u[0]*u[1] + mm[2]*u[0]*u[2] + mm[4]*u[1]*u[2]))

	h1 = LEN_SQ(ab);
	h2 = LEN_SQ(ac);
	h3 = LEN_SQ(ad);
	h4 = LEN_SQ(bc);
	h5 = LEN_SQ(bd);
	h6 = LEN_SQ(cd);

#undef LEN_SQ

	// Quality metric
	rap = h1 + h2 + h3 + h4 + h5 + h6;
	if (rap < 0)
		return vol;

	num = sqrt(rap) * rap;

	cal = det / num;

	if (!(cal >= 0 && cal <= 1))
		return vol;

	//double angle = tetquality(v1, v2, v3, v4, mm, 2);

	//if (angle < 2e-4) {
	//	return std::min({ vol, angle, cal * angle });
	//}

	return cal;// std::min(cal, vol);
}
//#pragma optimize("",on)

double DT::caltri33_ani(double v1[3], double v2[3], double v3[3], double* AniMetric) {
	double ab[3], ac[3], bc[3], mm[6];
	double area, det, h1, h2, h3, rap, num, cal;

	// 平均各向异性度量张量 (3个点)
	for (int i = 0; i < 6; i++) {
		mm[i] = (AniMetric[i] + AniMetric[6 + i] + AniMetric[12 + i]) / 3.0;
	}

	// 边向量
	for (int i = 0; i < 3; i++) {
		ab[i] = v2[i] - v1[i];
		ac[i] = v3[i] - v1[i];
		bc[i] = v3[i] - v2[i];
	}

	// 几何面积（向量叉积的一半）
	double nx = ab[1] * ac[2] - ab[2] * ac[1];
	double ny = ab[2] * ac[0] - ab[0] * ac[2];
	double nz = ab[0] * ac[1] - ab[1] * ac[0];
	area = 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
	if (area <= 0.0) return 0.0;

	// 各向异性度量张量的行列式
	det = mm[0] * (mm[3] * mm[5] - mm[4] * mm[4])
		- mm[1] * (mm[1] * mm[5] - mm[2] * mm[4])
		+ mm[2] * (mm[1] * mm[4] - mm[2] * mm[3]);

	double minor2 = mm[0] * mm[3] - mm[1] * mm[1];
	if (det <= 1e-14 || minor2<1e-14) return 0.0;

	det = sqrt(det) * area;

	// 在度量 mm 下的边长平方
#define LEN_SQ(u) (mm[0]*u[0]*u[0] + mm[3]*u[1]*u[1] + mm[5]*u[2]*u[2] + 2.0*(mm[1]*u[0]*u[1] + mm[2]*u[0]*u[2] + mm[4]*u[1]*u[2]))

	h1 = LEN_SQ(ab);
	h2 = LEN_SQ(ac);
	h3 = LEN_SQ(bc);

#undef LEN_SQ

	// 质量度量
	rap = h1 + h2 + h3;
	if (rap <= 0) return 0.0;

	num = rap * sqrt(rap);

	cal = det / num;

	return cal;
}

void DT::printfDihedral(double& minD, double& minAvgD, double& maxD, double& maxAvgD) {
    if (infolevel <= 0) return;
    const int minID = calculateDihedral(minD, minAvgD, maxD, maxAvgD);
    meshLogger->info("Dihedral: min: {:.8f} minAvg: {:.3f} max: {:.3f} maxAvg: {:.3f} minID: {}",
        minD, minAvgD, maxD, maxAvgD, minID);
}

int DT::calculateDihedral(double& minD, double& minAvgD, double& maxD, double& maxAvgD) {
    DTParallelScope parallelScope;
    struct alignas(64) Statistics {
        double minAngle = DBL_MAX, maxAngle = DBL_MIN, sumMin = 0, sumMax = 0;
        int count = 0, minID = -1;
    };
    const int workers = activeDTThreads(*this, 0, Elems.size());
    std::vector<Statistics> partial(workers);
#pragma omp parallel num_threads(workers) if(workers > 1)
    {
        Statistics local;
        std::vector<double> angles(6);
#pragma omp for schedule(static)
        for (int t = 0; t < static_cast<int>(Elems.size()); ++t) {
            if (isDelEle(t) || isvirtualtet(t) || ishulltet(t)) continue;
            const auto& p = Elems[t].form;
            double lo = 0, hi = 0;
            CalDihedral(Nodes[p[0]].pt, Nodes[p[1]].pt, Nodes[p[3]].pt, Nodes[p[2]].pt, lo, hi, angles);
            lo = RADIO2ANGLE(lo); hi = RADIO2ANGLE(hi);
            ++local.count; local.sumMin += lo; local.sumMax += hi;
            if (lo < local.minAngle) { local.minAngle = lo; local.minID = t; }
            local.maxAngle = std::max(local.maxAngle, hi);
        }
        partial[omp_get_thread_num()] = local;
    }
    Statistics total;
    for (const auto& local : partial) {
        total.count += local.count; total.sumMin += local.sumMin; total.sumMax += local.sumMax;
        if (local.minID >= 0 && (local.minAngle < total.minAngle ||
            (local.minAngle == total.minAngle && (total.minID < 0 || local.minID < total.minID)))) {
            total.minAngle = local.minAngle; total.minID = local.minID;
        }
        total.maxAngle = std::max(total.maxAngle, local.maxAngle);
    }
    minD = total.count ? total.minAngle : 0; maxD = total.count ? total.maxAngle : 0;
    minAvgD = total.count ? total.sumMin / total.count : 0;
    maxAvgD = total.count ? total.sumMax / total.count : 0;
    return total.minID;
}
//--------------------- Hilbert sort ---------------------
void DT::Hilbert(const std::vector<std::array<double, 3>>& V, std::vector<int>& order) {
    MeshStageLog stageLog(*this, "Hilbert ordering", 2);
	int  ngroup = 0, Vsize = V.size();
	std::vector<std::array<double, 4>> Varray;
	Varray.resize(Vsize);

	for (int i = 0; i < Vsize; i++) {
		for (int j = 0; j < 3; j++)
			Varray[i][j] = V[i][j];
		Varray[i][3] = order[i];
	}
	multiscale_sort(Varray, Vsize, ngroup);

	for (int i = 0; i < Vsize; i++) {
		order[i] = Varray[i][3];
	}
	return;
}

void DT::multiscale_sort(std::vector<std::array<double, 4>>& Varray, int Asize, int& depth)
{
	int middle = 0;
	if (Asize >= 64) {
		depth++;
		middle = Asize * 0.125;
		multiscale_sort(Varray, middle, depth);
	}
	hilbert_sort(Varray, middle, Asize - middle, 0, 0,
		minW[0] * 1.01, maxW[0] * 1.01, minW[1] * 1.01,
		maxW[1] * 1.01, minW[2] * 1.01, maxW[2] * 1.01, 0);
	return;
}

void DT::hilbert_sort(std::vector<std::array<double, 4>>& Varray,
	int bg, int Asize, int e, int d, double bxmin, double bxmax,
	double bymin, double bymax, double bzmin, double bzmax, int depth)
{
	double x1, x2, y1, y2, z1, z2;
	int p[9], w, e_w, d_w, k, ei, di;
	int n = 3, mask = 7;

	p[0] = 0;
	p[8] = Asize;
	p[4] = hilbert_split(Varray, bg, p[8], Trans[e][d][3], Trans[e][d][4],
		bxmin, bxmax, bymin, bymax, bzmin, bzmax);
	p[2] = hilbert_split(Varray, bg, p[4], Trans[e][d][1], Trans[e][d][2],
		bxmin, bxmax, bymin, bymax, bzmin, bzmax);
	p[1] = hilbert_split(Varray, bg, p[2], Trans[e][d][0], Trans[e][d][1],
		bxmin, bxmax, bymin, bymax, bzmin, bzmax);
	p[3] = hilbert_split(Varray, bg + p[2], p[4] - p[2],
		Trans[e][d][2], Trans[e][d][3], bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[2];
	p[6] = hilbert_split(Varray, bg + p[4], p[8] - p[4],
		Trans[e][d][5], Trans[e][d][6], bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[4];
	p[5] = hilbert_split(Varray, bg + p[4], p[6] - p[4],
		Trans[e][d][4], Trans[e][d][5], bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[4];
	p[7] = hilbert_split(Varray, bg + p[6], p[8] - p[6],
		Trans[e][d][6], Trans[e][d][7], bxmin, bxmax, bymin, bymax, bzmin, bzmax) + p[6];

	if ((depth + 1) == 52)
		return;

	for (w = 0; w < 8; w++) {
		if ((p[w + 1] - p[w]) > 8) {
			if (w == 0) {
				e_w = 0;
			}
			else {
				k = 2 * ((w - 1) / 2);
				e_w = k ^ (k >> 1); // = gc(k).
			}
			k = e_w;
			e_w = ((k << (d + 1)) & mask) | ((k >> (n - d - 1)) & mask);
			ei = e ^ e_w;
			if (w == 0) {
				d_w = 0;
			}
			else {
				d_w = ((w % 2) == 0) ? HbTab[w - 1] : HbTab[w];
			}
			di = (d + d_w + 1) % n;
			if (Trans[e][d][w] & 1) { // x-axis
				x1 = 0.5 * (bxmin + bxmax);
				x2 = bxmax;
			}
			else {
				x1 = bxmin;
				x2 = 0.5 * (bxmin + bxmax);
			}
			if (Trans[e][d][w] & 2) { // y-axis
				y1 = 0.5 * (bymin + bymax);
				y2 = bymax;
			}
			else {
				y1 = bymin;
				y2 = 0.5 * (bymin + bymax);
			}
			if (Trans[e][d][w] & 4) { // z-axis
				z1 = 0.5 * (bzmin + bzmax);
				z2 = bzmax;
			}
			else {
				z1 = bzmin;
				z2 = 0.5 * (bzmin + bzmax);
			}
			hilbert_sort(Varray, bg + p[w], p[w + 1] - p[w], ei, di,
				x1, x2, y1, y2, z1, z2, depth + 1);
		} // if (p[w+1] - p[w] > 1)
	}     // w
	return;
}

int DT::hilbert_split(std::vector<std::array<double, 4>>& Varray, int bg, int Asize, int gc0,
	int gc1, double bxmin, double bxmax, double bymin, double bymax, double bzmin, double bzmax)
{
	int axis, d, i, j;
	double split;

	axis = (gc0 ^ gc1) >> 1;

	if (axis == 0) {
		split = 0.5 * (bxmin + bxmax);
	}
	else if (axis == 1) {
		split = 0.5 * (bymin + bymax);
	}
	else { // == 2
		split = 0.5 * (bzmin + bzmax);
	}

	d = ((gc0 & (1 << axis)) == 0) ? 1 : -1;
	i = 0;
	j = Asize - 1;

	if (d > 0) {
		do {
			for (; i < Asize; i++) {
				if (Varray[i + bg][axis] >= split)break;
			}
			for (; j >= 0; j--) {
				if (Varray[j + bg][axis] < split)break;
			}
			if (i == (j + 1))break;
			std::swap(Varray[i + bg], Varray[j + bg]);
		} while (true);
	}
	else {
		do {
			for (; i < Asize; i++) {
				if (Varray[i + bg][axis] <= split)break;
			}
			for (; j >= 0; j--) {
				if (Varray[j + bg][axis] > split)break;
			}
			if (i == (j + 1))break;
			std::swap(Varray[i + bg], Varray[j + bg]);
		} while (true);
	}
	return i;
}

// Linear algebra operators.
double DT::dot(double* v1, double* v2) {
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// cross() computes the cross product: n = v1 cross v2.
void DT::cross(double* v1, double* v2, double* n)
{
	n[0] = v1[1] * v2[2] - v1[2] * v2[1];
	n[1] = v1[2] * v2[0] - v1[0] * v2[2];
	n[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

// distance() computes the Euclidean distance between two points.
double DT::distance(double* p1, double* p2) {
	return sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) +
		(p2[1] - p1[1]) * (p2[1] - p1[1]) +
		(p2[2] - p1[2]) * (p2[2] - p1[2]));
}

double DT::distance2(double* p1, double* p2) {
	return norm2(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
}
double DT::calArea(double* p1, double* p2, double* p3) {
	double Cross[3] = { 0 };
	double p1p2[3] = { p2[0] - p1[0] , p2[1] - p1[1] ,p2[2] - p1[2] };
	double p1p3[3] = { p3[0] - p1[0] , p3[1] - p1[1], p3[2] - p1[2] };

	cross(p1p2, p1p3, Cross);

	return lenvec(Cross) / 2.0;
}
double DT::lenvec(double* x) {
	return sqrt(norm2(x[0], x[1], x[2]));
}
double DT::norm2(double x, double y, double z) {
	return (x) * (x)+(y) * (y)+(z) * (z);
}
double DT::cal_ani_length(double* v1, double* v2, double* m1, double* m2) {
	// Compute edge vector
	double ux = v2[0] - v1[0];
	double uy = v2[1] - v1[1];
	double uz = v2[2] - v1[2];

	// Compute squared length in metric m1
	double dd1 = m1[0] * ux * ux + m1[3] * uy * uy + m1[5] * uz * uz
		+ 2.0 * (m1[1] * ux * uy + m1[2] * ux * uz + m1[4] * uy * uz);
	if (dd1 <= 0.0) dd1 = 0.0;

	// Compute squared length in metric m2
	double dd2 = m2[0] * ux * ux + m2[3] * uy * uy + m2[5] * uz * uz
		+ 2.0 * (m2[1] * ux * uy + m2[2] * ux * uz + m2[4] * uy * uz);
	if (dd2 <= 0.0) dd2 = 0.0;

	// Evaluate the length based on similarity or weighted average
	double length;
	if (fabs(dd1 - dd2) < 0.05) {
		length = sqrt(0.5 * (dd1 + dd2));
	}
	else {
		length = (sqrt(dd1) + sqrt(dd2) + 4.0 * sqrt(0.5 * (dd1 + dd2))) / 6.0;
	}

	return length;
}

int DT::Interpolate_met33_ani(double* m1, double* m2, double* m, double s) {
	constexpr double EPS = 1e-12;

	Eigen::Matrix3d M1, M2;
	M1 << m1[0], m1[1], m1[2],
		m1[1], m1[3], m1[4],
		m1[2], m1[4], m1[5];

	M2 << m2[0], m2[1], m2[2],
		m2[1], m2[3], m2[4],
		m2[2], m2[4], m2[5];

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solverM1(M1);
	if (solverM1.info() != Eigen::Success) return 0;

	Eigen::Matrix3d P = solverM1.eigenvectors();
	Eigen::Vector3d lambda = solverM1.eigenvalues();

	for (int i = 0; i < 3; ++i) {
		if (lambda[i] < EPS) return 0;
	}

	Eigen::Matrix3d LambdaSqrt = lambda.cwiseSqrt().asDiagonal();
	Eigen::Matrix3d LambdaInvSqrt = lambda.cwiseSqrt().cwiseInverse().asDiagonal();

	Eigen::Matrix3d sqrtM1 = P * LambdaSqrt * P.transpose();
	Eigen::Matrix3d is = P * LambdaInvSqrt * P.transpose();

	Eigen::Matrix3d isnis = is * M2 * is;
	isnis = 0.5 * (isnis + isnis.transpose());

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solverIsnis(isnis);
	if (solverIsnis.info() != Eigen::Success) return 0;

	Eigen::Matrix3d Q = solverIsnis.eigenvectors();
	Eigen::Vector3d mu = solverIsnis.eigenvalues();

	for (int i = 0; i < 3; ++i) {
		if (mu[i] < -EPS) return 0;
		mu[i] = std::max(mu[i], EPS);

		double dd = s * std::sqrt(mu[i]) + (1.0 - s);
		dd = dd * dd;
		if (dd < EPS) return 0;

		mu[i] = mu[i] / dd;
	}

	Eigen::Matrix3d M = sqrtM1 * Q * mu.asDiagonal() * Q.transpose() * sqrtM1;
	M = 0.5 * (M + M.transpose());

	m[0] = M(0, 0);
	m[1] = M(0, 1);
	m[2] = M(0, 2);
	m[3] = M(1, 1);
	m[4] = M(1, 2);
	m[5] = M(2, 2);

	return 1;
}

//int DT::Interpolate_met33_ani(double* m1, double* m2, double* m, double s) {
//	constexpr double EPS = 1e-12;
//
//	// 构造两个 3x3 对称张量矩阵
//	Eigen::Matrix3d M1, M2;
//	M1 << m1[0], m1[1], m1[2],
//		m1[1], m1[3], m1[4],
//		m1[2], m1[4], m1[5];
//
//	M2 << m2[0], m2[1], m2[2],
//		m2[1], m2[3], m2[4],
//		m2[2], m2[4], m2[5];
//
//	// M1 的特征分解
//	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solverM1(M1);
//	if (solverM1.info() != Eigen::Success) return 0;
//
//	Eigen::Matrix3d P = solverM1.eigenvectors();
//	Eigen::Vector3d lambda = solverM1.eigenvalues();
//
//	for (int i = 0; i < 3; ++i)
//		if (lambda[i] < EPS) return 0;
//
//	// 计算 is = P * diag(1/sqrt(lambda)) * P^T
//	Eigen::Matrix3d LambdaInvSqrt = lambda.cwiseSqrt().cwiseInverse().asDiagonal();
//	Eigen::Matrix3d is = P * LambdaInvSqrt * P.transpose();
//
//	// isnis = is * M2 * is
//	Eigen::Matrix3d isnis = is * M2 * is;
//
//	// 对 isnis 进行特征分解
//	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solverIsnis(isnis);
//	if (solverIsnis.info() != Eigen::Success) return 0;
//	Eigen::Matrix3d Q = solverIsnis.eigenvectors();
//	Eigen::Vector3d mu = solverIsnis.eigenvalues();
//
//	for (int i = 0; i < 3; ++i) {
//		if (mu[i] < 0.0) return 0;
//		double dd = s * std::sqrt(mu[i]) + (1.0 - s);
//		dd = dd * dd;
//		if (dd < EPS) return 0;
//		mu[i] = mu[i] / dd;
//	}
//
//	// 构造 P = is * Q
//	Eigen::Matrix3d Pfull = is * Q;
//
//	// 求逆
//	Eigen::Matrix3d PfullInv;
//	bool invertible;
//	double det;
//	Pfull.computeInverseWithCheck(PfullInv, invertible, det);
//	if (!invertible) return 0;
//
//	// 插值后的张量：P^{-T} * diag(mu) * P^{-1}
//	Eigen::Matrix3d M = PfullInv.transpose() * mu.asDiagonal() * PfullInv;
//
//	// 转换为对称张量的压缩格式
//	m[0] = M(0, 0);
//	m[1] = M(0, 1);
//	m[2] = M(0, 2);
//	m[3] = M(1, 1);
//	m[4] = M(1, 2);
//	m[5] = M(2, 2);
//
//	return 1;
//}

int DT::Interpolate_met(int p1, int p2, int newp, double s) {
	if (AniSol.size() != 0) {
		double m1[6] = { 0 }, m2[6] = { 0 }, m_inter[6] = { 0 };
		for (int mm = 0; mm < 6; mm++) {
			m1[mm] = AniSol[p1][mm];
			m2[mm] = AniSol[p2][mm];
		}

		if (AniSol.size() <= newp) {
			AniSol.resize(std::max(AniSol.size() * 2, static_cast<size_t>(newp + 1)));
		}

		if (Interpolate_met33_ani(m1, m2, m_inter, s)) {
			for (int mm = 0; mm < 6; mm++) {
				AniSol[newp][mm] = m_inter[mm];
			}
		}
		else {
			for (int mm = 0; mm < 6; mm++) {
				AniSol[newp][mm] = (m1[mm] + m2[mm]) / 2.0;
			}
		}
	}
	return 0;
}

/*
* Input four point
* Calculate the centre and radius of a tetrahedron
* Input can be four points (tetrahedron) / three points (triangle)
* Calculate the outer sphere / outer circle
* Return false if degenerate, coplanar/collinear
*/
bool DT::calCircum(double* pa, double* pb, double* pc, double* pd, double* cent, double* radius) {
	double A[4][4], rhs[4], D;
	int indx[4];

	// Compute the coefficient matrix A (3x3).
	A[0][0] = pb[0] - pa[0];
	A[0][1] = pb[1] - pa[1];
	A[0][2] = pb[2] - pa[2];
	A[1][0] = pc[0] - pa[0];
	A[1][1] = pc[1] - pa[1];
	A[1][2] = pc[2] - pa[2];
	if (pd != NULL) {
		A[2][0] = pd[0] - pa[0];
		A[2][1] = pd[1] - pa[1];
		A[2][2] = pd[2] - pa[2];
	}
	else {
		cross(A[0], A[1], A[2]);
	}

	// Compute the right hand side vector b (3x1).
	rhs[0] = 0.5 * dot(A[0], A[0]);
	rhs[1] = 0.5 * dot(A[1], A[1]);
	if (pd != NULL) {
		rhs[2] = 0.5 * dot(A[2], A[2]);
	}
	else {
		rhs[2] = 0.0;
	}

	// Solve the 3 by 3 equations use LU decomposition with partial pivoting
	//   and backward and forward substitute..
	if (!lu_decmp(A, 3, indx, &D, 0)) {
		if (radius != (double*)NULL) *radius = 0.0;
		return false;
	}
	lu_solve(A, 3, indx, rhs, 0);
	if (cent != (double*)NULL) {
		cent[0] = pa[0] + rhs[0];
		cent[1] = pa[1] + rhs[1];
		cent[2] = pa[2] + rhs[2];
	}
	if (radius != (double*)NULL) {
		*radius = sqrt(rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2]);
	}
	return true;
}
/*
* Check if three vertices (from left to right):
* left, mid, and right are collinear.
*/
bool DT::P_in_Line(double* mid, double* left, double* right) {
	double eps = distance(left, right) / 1000.0;
	if (mid[0] < std::min(left[0], right[0]) - eps || mid[0] > std::max(left[0], right[0]) + eps ||
		mid[1] < std::min(left[1], right[1]) - eps || mid[1] > std::max(left[1], right[1]) + eps ||
		mid[2] < std::min(left[2], right[2]) - eps || mid[2] > std::max(left[2], right[2]) + eps)
	{
		return false;
	}

	double seg[3] = { right[0] - left[0], right[1] - left[1], right[2] - left[2] };
	double v[3] = { mid[0] - left[0],   mid[1] - left[1],   mid[2] - left[2] };

	double seg_len2 = dot(seg, seg);
	double c[3];
	cross(v, seg, c);
	if (dot(c, c) > eps * eps * seg_len2) {
		return false;
	}

	double proj = dot(v, seg);
	if (proj <= 0 || proj >= seg_len2) {
		return false;
	}

	return true;
}

bool DT::isParallel(double* p1, double* p2, double* p3, double* p4)
{
	const double eps = 1e-12;

	// 两条线方向向量
	double v1[3] = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };
	double v2[3] = { p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2] };

	// 叉积
	double c[3] = {
		v1[1] * v2[2] - v1[2] * v2[1],
		v1[2] * v2[0] - v1[0] * v2[2],
		v1[0] * v2[1] - v1[1] * v2[0]
	};

	// 若叉积长度接近0，则平行
	double norm2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];

	return norm2 < eps;
}
/*
* lu_decmp()    Compute the LU decomposition of a matrix.Reference tetgen
*/
bool DT::lu_decmp(double lu[4][4], int n, int* ps, double* d, int N) {
	double scales[4] = {0};
	double pivot, biggest, mult, tempf;
	int pivotindex = 0;
	int i, j, k;

	*d = 1.0;                                      // No row interchanges yet.

	for (i = N; i < n + N; i++) {                  // For each row.
		// Find the largest element in each row for row equilibration
		biggest = 0.0;
		for (j = N; j < n + N; j++)
			if (biggest < (tempf = fabs(lu[i][j])))
				biggest = tempf;
		if (biggest != 0.0)
			scales[i] = 1.0 / biggest;
		else {
			scales[i] = 0.0;
			return false;                          // Zero row: singular matrix.
		}
		ps[i] = i;                                 // Initialize pivot sequence.
	}

	for (k = N; k < n + N - 1; k++) {                      // For each column.
		// Find the largest element in each column to pivot around.
		biggest = 0.0;
		for (i = k; i < n + N; i++) {
			if (biggest < (tempf = fabs(lu[ps[i]][k]) * scales[ps[i]])) {
				biggest = tempf;
				pivotindex = i;
			}
		}
		if (biggest == 0.0) {
			return false;                         // Zero column: singular matrix.
		}
		if (pivotindex != k) {                    // Update pivot sequence.
			j = ps[k];
			ps[k] = ps[pivotindex];
			ps[pivotindex] = j;
			*d = -(*d);                           // ...and change the parity of d.
		}

		// Pivot, eliminating an extra variable  each time
		pivot = lu[ps[k]][k];
		for (i = k + 1; i < n + N; i++) {
			lu[ps[i]][k] = mult = lu[ps[i]][k] / pivot;
			if (mult != 0.0) {
				for (j = k + 1; j < n + N; j++)
					lu[ps[i]][j] -= mult * lu[ps[k]][j];
			}
		}
	}
	// (lu[ps[n + N - 1]][n + N - 1] == 0.0) ==> A is singular.
	return lu[ps[n + N - 1]][n + N - 1] != 0.0;
}
/*
* lu_solve()    Solves the linear equation:  Ax = b,Reference tetgen
*/
void DT::lu_solve(double lu[4][4], int n, int* ps, double* b, int N) {
	int i, j;
	double X[4]={0}, dot;

	for (i = N; i < n + N; i++) X[i] = 0.0;

	// Vector reduction using U triangular matrix.
	for (i = N; i < n + N; i++) {
		dot = 0.0;
		for (j = N; j < i + N; j++)
			dot += lu[ps[i]][j] * X[j];
		X[i] = b[ps[i]] - dot;
	}

	// Back substitution, in L triangular matrix.
	for (i = n + N - 1; i >= N; i--) {
		dot = 0.0;
		for (j = i + 1; j < n + N; j++)
			dot += lu[ps[i]][j] * X[j];
		X[i] = (X[i] - dot) / lu[ps[i]][i];
	}

	for (i = N; i < n + N; i++) b[i] = X[i];
	return;
}

double DT::insphere_s(int a, int b, int c, int d, int e) {
	double sign, * pa, * pb, * pc, * pd, * pe;

	pa = Nodes[a].pt;
	pb = Nodes[b].pt;
	pc = Nodes[c].pt;
	pd = Nodes[d].pt;
	pe = Nodes[e].pt;
	sign = dt::GEOM_FUNC::insphere(pa, pb, pc, pd, pe);
	if (sign != 0.0) {
		return sign;
	}

	// Symbolic perturbation.
	double pt[5][3];
	double oriA, oriB;
	int swaps, count, n, i, mp[5];
	mp[0] = a; mp[1] = b; mp[2] = c; mp[3] = d; mp[4] = e;

	// Sort the five points such that their indices are in the increasing order.
	swaps = 0; // Record the total number of swaps.
	n = 5;
	do {
		count = 0;
		n = n - 1;
		for (i = 0; i < n; i++) {
			if (mp[i] > mp[i + 1]) {
				std::swap(mp[i], mp[i + 1]);
				count++;
			}
		}
		swaps += count;
	} while (count > 0); // Continue if some points are swapped.

	pa = Nodes[mp[0]].pt;
	pb = Nodes[mp[1]].pt;
	pc = Nodes[mp[2]].pt;
	pd = Nodes[mp[3]].pt;
	pe = Nodes[mp[4]].pt;

	oriA = dt::GEOM_FUNC::orient3d(pb, pc, pd, pe);
	if (oriA != 0.0) {
		// Flip the sign if there are odd number of swaps.
		if ((swaps % 2) != 0) oriA = -oriA;
		return oriA;
	}

	oriB = -dt::GEOM_FUNC::orient3d(pa, pc, pd, pe);
	// Flip the sign if there are odd number of swaps.
	if ((swaps % 2) != 0) oriB = -oriB;
	return oriB;
}

void DT::calnormal(int b, int c, int d, double noraml[]) {
	//normal [b,c,d]->a
	double* pb = Nodes[b].pt;
	double* pc = Nodes[c].pt;
	double* pd = Nodes[d].pt;
	double bd[3] = { pd[0] - pb[0] ,pd[1] - pb[1] ,pd[2] - pb[2] };
	double bc[3] = { pc[0] - pb[0] ,pc[1] - pb[1] ,pc[2] - pb[2] };

	cross(bc, bd, noraml);
	return;
}

void DT::tensorproduct33(double* v1, double* v2, double* ans) {
	ans[0] = v1[0] * v2[0]; ans[1] = v1[1] * v2[0]; ans[2] = v1[2] * v2[0];
	ans[3] = v1[0] * v2[1]; ans[4] = v1[1] * v2[1]; ans[5] = v1[2] * v2[1];
	ans[6] = v1[0] * v2[2]; ans[7] = v1[1] * v2[2]; ans[8] = v1[2] * v2[2];
	return;
}
bool DT::inverseM(double* in, double* out) {
	double det = in[0] * (in[4] * in[8] - in[5] * in[7]) -
		in[1] * (in[3] * in[8] - in[5] * in[6]) +
		in[2] * (in[3] * in[7] - in[4] * in[6]);

	if (det == 0) {
		//std::cerr << "Error: The in is singular, its inverse does not exist." << std::endl;
		return false;
	}

	double adj[9] = { 0 };
	adj[0] = in[4] * in[8] - in[5] * in[7];
	adj[1] = in[2] * in[7] - in[1] * in[8];
	adj[2] = in[1] * in[5] - in[2] * in[4];
	adj[3] = in[5] * in[6] - in[3] * in[8];
	adj[4] = in[0] * in[8] - in[2] * in[6];
	adj[5] = in[2] * in[3] - in[0] * in[5];
	adj[6] = in[3] * in[7] - in[4] * in[6];
	adj[7] = in[1] * in[6] - in[0] * in[7];
	adj[8] = in[0] * in[4] - in[1] * in[3];

	for (int i = 0; i < 9; ++i) {
		out[i] = adj[i] / det;
	}
	return true;
}
void DT::vecTimesMatrix13_33(double* vec, double* Matrix, double* vecans) {
	for (int i = 0; i < 3; ++i) {
		vecans[i] = 0.0;
		for (int j = 0; j < 3; ++j) {
			//vecans[i] += vec[j] * Matrix[i * 3 + j];
			vecans[i] += vec[j] * Matrix[j * 3 + i];
		}
	}
	return;
}
void  DT::projectPointToPlane(double* v1, double* v2, double* v3, double* p, double* projection) {
	// Prepare data
	double v12[3] = { v2[0] - v1[0] ,v2[1] - v1[1] ,v2[2] - v1[2] }; // Vector from v1 to v2
	double v13[3] = { v3[0] - v1[0] ,v3[1] - v1[1] ,v3[2] - v1[2] }; // Vector from v1 to v3
	double p1[3] = { p[0] - v1[0] ,p[1] - v1[1] ,p[2] - v1[2] }; // Vector from v1 to point p

	double nor[3];
	// Calculate normal of the plane
	cross(v12, v13, nor);
	double lennor = lenvec(nor);
	// Compute the perpendicular distance from point p to the plane
	double  d = dot(nor, p1) / lennor;

	// Calculate the projection point on the plane
	for (int i = 0; i < 3; i++)
		projection[i] = p[i] - d * nor[i] / lennor;

	return;
}
void DT::calFactor(int n, int& a, int& b, int& c) {
	int minSum = INT_MAX;
	for (int i = 1; i * i * i <= n; ++i) {
		if (n % i == 0) {
			for (int j = i; j * j <= n / i; ++j) {
				if ((n / i) % j == 0) {
					int k = n / (i * j);
					int sum = i + j + k;
					if (sum < minSum) {
						minSum = sum;
						a = i;
						b = j;
						c = k;
					}
				}
			}
		}
	}
	return;
}

void DT::calBarycenter(int i, double* pnt) {
	//get barycenter of this tet
	for (int j = 0; j < 3; j++) {
		pnt[j] = 0;
		for (int k = 0; k < 4; k++)
			pnt[j] += Nodes[Elems[i].form[k]].pt[j];
		pnt[j] /= 4.0;
	}
	return;
}
void DT::distanceToPlane(double p1[3], double p2[3], double p3[3], double p4[3], double& dis) {
	// 计算平面法向量（p1p2和p1p3的叉积）
	double v1[3] = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] }; // 向量p1p2
	double v2[3] = { p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2] }; // 向量p1p3
	double normal[3]; // 平面法向量

	// 叉积计算
	cross(v1, v2, normal);

	// 计算平面方程 Ax + By + Cz + D = 0 中的 A, B, C, D
	double A = normal[0];
	double B = normal[1];
	double C = normal[2];
	double D = -(A * p1[0] + B * p1[1] + C * p1[2]);

	// 计算点p4到平面的距离
	dis = std::fabs(A * p4[0] + B * p4[1] + C * p4[2] + D) / std::sqrt(A * A + B * B + C * C);
	return;
}

double DT::segmentSegmentDistance(double* p1, double* q1, double* p2, double* q2){
	const double EPS = 1e-12;

	double u[3], v[3], w[3];
	vsub3(q1, p1, u);   // u = q1 - p1
	vsub3(q2, p2, v);   // v = q2 - p2
	vsub3(p1, p2, w);   // w = p1 - p2

	double a = dot(u, u); // always >= 0
	double b = dot(u, v);
	double c = dot(v, v); // always >= 0
	double d = dot(u, w);
	double e = dot(v, w);
	double D = a * c - b * b; // always >= 0

	double sc, sN, sD = D;
	double tc, tN, tD = D;

	// 两条线段都退化成点
	if (a < EPS && c < EPS) {
		return distance(p1, p2);
	}

	// 第一条退化成点
	if (a < EPS) {
		sc = 0.0;
		tc = (c < EPS ? 0.0 : e / c);
		tc = std::max(0.0, std::min(1.0, tc));
	}
	// 第二条退化成点
	else if (c < EPS) {
		tc = 0.0;
		sc = -d / a;
		sc = std::max(0.0, std::min(1.0, sc));
	}
	else {
		// 一般情况
		if (D < EPS) {
			// 近乎平行
			sN = 0.0;
			sD = 1.0;
			tN = e;
			tD = c;
		}
		else {
			sN = (b * e - c * d);
			tN = (a * e - b * d);

			if (sN < 0.0) {
				sN = 0.0;
				tN = e;
				tD = c;
			}
			else if (sN > sD) {
				sN = sD;
				tN = e + b;
				tD = c;
			}
		}

		if (tN < 0.0) {
			tN = 0.0;

			if (-d < 0.0) {
				sN = 0.0;
			}
			else if (-d > a) {
				sN = sD;
			}
			else {
				sN = -d;
				sD = a;
			}
		}
		else if (tN > tD) {
			tN = tD;

			if ((-d + b) < 0.0) {
				sN = 0.0;
			}
			else if ((-d + b) > a) {
				sN = sD;
			}
			else {
				sN = (-d + b);
				sD = a;
			}
		}

		sc = (std::abs(sN) < EPS ? 0.0 : sN / sD);
		tc = (std::abs(tN) < EPS ? 0.0 : tN / tD);
	}

	double dP[3];
	dP[0] = w[0] + sc * u[0] - tc * v[0];
	dP[1] = w[1] + sc * u[1] - tc * v[1];
	dP[2] = w[2] + sc * u[2] - tc * v[2];

	return lenvec(dP);
}

void DT::vsub3( double a[3],  double b[3], double r[3]) {
	r[0] = a[0] - b[0]; r[1] = a[1] - b[1]; r[2] = a[2] - b[2];
	return;
}
double DT::calculateTriangleAngles(const double* A, const double* B, const double* C, std::vector<double>& angles) {
	angles.resize(3);
	angles[0] = angles[1] = angles[2] = 181;

	//Cal A
	double AB[3] = { 0 }, AC[3] = { 0 };
	for (int i = 0; i < 3; i++) {
		AB[i] = B[i] - A[i];
		AC[i] = C[i] - A[i];
	}
	double aa = dot(AB, AC) / (lenvec(AB) * lenvec(AC));
	if (abs(aa) > 1)
		return 0;
	double angleA = RADIO2ANGLE(acos(aa));
	angles[0] = angleA;

	//Cal B
	double BA[3] = { 0 }, BC[3] = { 0 };
	for (int i = 0; i < 3; i++) {
		BA[i] = A[i] - B[i];
		BC[i] = C[i] - B[i];
	}
	aa = dot(BA, BC) / (lenvec(BA) * lenvec(BC));
	if (abs(aa) > 1)
		return 0;
	double angleB = RADIO2ANGLE(acos(aa));
	angles[1] = angleB;
	//Cal C
	double  angleC = 180 - angleA - angleB;
	angles[2] = angleC;

	return std::min(std::min(angleA, angleB), angleC);
}
double DT::calVolume(int i) {
	if (isDelEle(i))
		return 0;
	double* p0 = Nodes[Elems[i].form[0]].pt;
	double* p1 = Nodes[Elems[i].form[1]].pt;
	double* p2 = Nodes[Elems[i].form[2]].pt;
	double* p3 = Nodes[Elems[i].form[3]].pt;
	return dt::GEOM_FUNC::orient3d(p0, p1, p3, p2);
}

void DT::addRandomP(double x, double y, double z, int n) {
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<double> distr(0.0, 1.0);
	double dis = 1e-3;
	for (int i = 0; i < n; i++) {

		int newp = addNode(x+ dis * distr(eng),y+ dis * distr(eng),z+ dis * distr(eng),0);

		int iSrch = 0;
		while (isDelEle(iSrch))
			iSrch++;

		std::vector<int> srchtet = { iSrch };
		int ret = BW_insert_vertex(newp, srchtet, 1);
		if (ret <= 0) {
			DelNod(newp);
			i--;
		}
	}
	return;
}

//void DT::export_ring_csv_min(int iNod, const double oldpos[3], int mode /*整型标志，便于扩展*/) {
//	return;
//	// 0) 新旧位移（用当前节点坐标做 newpos）
//	const double newpos[3] = { Nodes[iNod].pt[0], Nodes[iNod].pt[1], Nodes[iNod].pt[2] };
//	const double disp[3] = { newpos[0] - oldpos[0], newpos[1] - oldpos[1], newpos[2] - oldpos[2] };
//
//	// 1) 取以 iNod 为顶点的四面体列表
//	std::vector<int> sph;
//	findSphere(iNod, sph);
//	const int T = (int)sph.size();
//
//	// 2) 提取对面三点（全局 id），并按“首次出现”收集唯一邻点（不排序）
//	std::vector<std::array<int, 3>> tris_global;
//	tris_global.reserve(T);
//
//	std::unordered_set<int> seen;  seen.reserve(T * 3 * 2);
//	std::vector<int>        nbrs;  nbrs.reserve(T * 3);
//
//	for (int ei : sph) {
//		int j = 0; while (j < 4 && Elems[ei].form[j] != iNod) ++j;
//		int a, b, c, d; DFC(j, d, a, b, c); // d==j
//
//		int ga = Elems[ei].form[a];
//		int gb = Elems[ei].form[b];
//		int gc = Elems[ei].form[c];
//
//		tris_global.push_back({ ga, gb, gc });
//
//		if (seen.insert(ga).second) nbrs.push_back(ga);
//		if (seen.insert(gb).second) nbrs.push_back(gb);
//		if (seen.insert(gc).second) nbrs.push_back(gc);
//	}
//
//	// （可选）如果以后需要几何稳定性，可按距离排序：
//	// auto dist2 = [&](int nid){ double dx=Nodes[nid].pt[0]-oldpos[0], dy=..., dz=...; return dx*dx+dy*dy+dz*dz; };
//	// std::sort(nbrs.begin(), nbrs.end(), [&](int a,int b){ return dist2(a) < dist2(b); });
//
//	const int P = (int)nbrs.size();
//
//	// 3) 局部尺度 s = 邻距中位数（欧氏），下限保护
//	auto dist = [&](int nid)->double {
//		double dx = Nodes[nid].pt[0] - oldpos[0];
//		double dy = Nodes[nid].pt[1] - oldpos[1];
//		double dz = Nodes[nid].pt[2] - oldpos[2];
//		double r2 = dx * dx + dy * dy + dz * dz;
//		return std::sqrt(r2 > 1e-30 ? r2 : 1e-30);
//		};
//	double s = 1.0;
//	if (P > 0) {
//		std::vector<double> d; d.reserve(P);
//		for (int nid : nbrs) d.push_back(dist(nid));
//		std::nth_element(d.begin(), d.begin() + d.size() / 2, d.end());
//		s = std::max(d[d.size() / 2], 1e-12);
//	}
//
//	if (s < 1e-4) return;  // 避免归一化爆炸
//	if (s > 1.0)  return;  // 通常 s 不应大于 1
//	// 4) 归一化位移幅度 ndisp，作为唯一导出判据
//	const double disp_norm = std::sqrt(disp[0] * disp[0] + disp[1] * disp[1] + disp[2] * disp[2]);
//	const double ndisp = disp_norm / s;
//	if (!(ndisp >= 1e-3 && ndisp <= 1)) {
//		return;  // 不导出
//	}
//
//	// 4) 建立本地索引映射（按“首次出现顺序”）
//	std::unordered_map<int, int> lid; lid.reserve(P * 2);
//	for (int idx = 0; idx < P; ++idx) lid[nbrs[idx]] = idx + 1;
//
//	// 5) 写 CSV：meta + 点表 + 三角面表
//	std::string fname = "E:\\SmoothTest\\1\\vol_" + std::to_string(smid++) + ".csv";
//	std::ofstream fout(fname);
//	if (!fout.is_open()) return;
//
//	fout << std::fixed << std::setprecision(10);
//
//	// meta：位移用归一化标签；原点视为 (0,0,0)
//	const double ndx = disp[0] / s, ndy = disp[1] / s, ndz = disp[2] / s;
//	fout << "type,mode,dx,dy,dz,s,P,T\n";
//	fout << "meta," << mode << "," << ndx << "," << ndy << "," << ndz << ","
//		<< s << "," << P << "," << T << "\n";
//
//	// 点表：相对 oldpos 并除以 s；索引 1..P
//	fout << "type,idx,vx,vy,vz\n";
//	for (int idx = 0; idx < P; ++idx) {
//		int nid = nbrs[idx];
//		double vx = (Nodes[nid].pt[0] - oldpos[0]) / s;
//		double vy = (Nodes[nid].pt[1] - oldpos[1]) / s;
//		double vz = (Nodes[nid].pt[2] - oldpos[2]) / s;
//		fout << "pt," << (idx + 1) << "," << vx << "," << vy << "," << vz << "\n";
//	}
//
//	// 三角面表：用本地点索引表示，拓扑保持与遍历顺序一致
//	fout << "type,ia,ib,ic\n";
//	for (const auto& tri : tris_global) {
//		int ia = lid.at(tri[0]);
//		int ib = lid.at(tri[1]);
//		int ic = lid.at(tri[2]);
//		fout << "tri," << ia << "," << ib << "," << ic << "\n";
//	}
//	fout.close();
//
//	return;
//}


std::chrono::high_resolution_clock::time_point DT::getTime_now() {
	return std::chrono::high_resolution_clock::now();
}

double DT::getTime(std::chrono::high_resolution_clock::time_point t1,
	std::chrono::high_resolution_clock::time_point t2) {
	return std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() * 1e-9;
}

void DT::prinfPnt(int iNod) {
	return;
}

void DT::outHullTri(std::string filename) {
	int i, j;
	Mesh mesh;
	mesh.V.reserve(Nodes.size());
	//mesh.S.resize(SurEdgs.size());
	mesh.F.reserve(Elems.size());
	for (i = 0; i < Nodes.size(); i++) {//have ghost
		mesh.V.push_back({ Nodes[i].pt[0],Nodes[i].pt[1],Nodes[i].pt[2] });
	}

	for (i = 0; i < Elems.size(); i++) {
		if (ishulltet(i))
			mesh.F.push_back({ Elems[i].form[0],Elems[i].form[1],Elems[i].form[2],0 });
	}

	dt::writeVTK(filename, mesh, true);
	return;
}

void DT::outTempMesh(std::string filename) {
	int i, j;
	Mesh mesh;
	mesh.V.reserve(Nodes.size());
	//mesh.S.resize(SurEdgs.size());
	mesh.F.reserve(SurTris.size());
	mesh.T.reserve(Elems.size());
	for (i = 0; i < Nodes.size(); i++) {//have ghost
		mesh.V.push_back({ Nodes[i].pt[0],Nodes[i].pt[1],Nodes[i].pt[2] });
	}
	//for (i = 0; i < SurEdgs.size(); i++) {//have ghost
	//	if (isDelSurEdg(i))
	//		continue;
	//	mesh.S.push_back({ SurEdgs[i].iStart,SurEdgs[i].iEnd });
	//}
	for (i = 0; i < SurTris.size(); i++) {//have ghost
		if (isDelSurTri(i) || SurTris[i].info >= nSurTris)
			continue;
		mesh.F.push_back({ SurTris[i].form[0],SurTris[i].form[1],SurTris[i].form[2] });
	}
	for (i = 0; i < Elems.size(); i++) {
		if (ishulltet(i) || isDelEle(i) /*|| isvirtualtet(i)*/)
			continue;
		mesh.T.push_back({ Elems[i].form[0],Elems[i].form[1],Elems[i].form[2],Elems[i].form[3],Elems[i].geo });
	}
	
	//out mesh
	auto& V = mesh.V;
	auto& S = mesh.S;
	auto& F = mesh.F;
	auto& T = mesh.T;
	auto& PS = mesh.pointSize;

	FILE* outFile = fopen(filename.c_str(), "w");
	if (outFile == nullptr) {
		std::cout << "Write VTK file failed. - " << filename << "\n";
		return;
	}

	std::cout << "Writing mesh to - " << filename << "\n";
	fprintf(outFile, "# vtk DataFile Version 2.0\n");
	fprintf(outFile, "TetWild Mesh\n");
	fprintf(outFile, "ASCII\n");
	fprintf(outFile, "DATASET UNSTRUCTURED_GRID\n");

	// POINTS
	fprintf(outFile, "POINTS %d double\n", (int)V.size());
	for (int i = 0; i < (int)V.size(); i++)
		fprintf(outFile, "%.16lf %.16lf %.16lf\n", V[i][0], V[i][1], V[i][2]);

	// CELLS
	fprintf(outFile, "CELLS %ld %ld\n",
		S.size() + F.size() + T.size(),
		S.size() * 3 + F.size() * 4 + T.size() * 5);

	for (int i = 0; i < (int)S.size(); i++)
		fprintf(outFile, "2 %d %d \n", S[i][0], S[i][1]);

	for (int i = 0; i < (int)F.size(); i++)
		fprintf(outFile, "3 %d %d %d\n", F[i][0], F[i][1], F[i][2]);


	for (int i = 0; i < (int)T.size(); i++)
		fprintf(outFile, "4 %d %d %d %d\n", T[i][0], T[i][1], T[i][2], T[i][3]);

	// CELL_TYPES
	fprintf(outFile, "CELL_TYPES  %ld\n", S.size() + F.size() + T.size());
	for (int i = 0; i < (int)S.size(); i++)
		fprintf(outFile, "3\n");   // VTK_LINE
	for (int i = 0; i < (int)F.size(); i++)
		fprintf(outFile, "5\n");   // VTK_TRIANGLE

	for (int i = 0; i < (int)T.size(); i++)
		fprintf(outFile, "10\n");      // VTK_TETRA

	// -----------------------------
	// POINT_DATA
	// -----------------------------
	//bool hasPointData = (!PS.empty()) || (!AniSol.empty());
	//if (hasPointData) {
	//	fprintf(outFile, "POINT_DATA  %ld\n", V.size());

	//	// PointSize
	//	if (!PS.empty()) {
	//		fprintf(outFile, "SCALARS PointSize float 1\n");
	//		fprintf(outFile, "LOOKUP_TABLE default\n");
	//		for (int i = 0; i < (int)PS.size(); i++)
	//			fprintf(outFile, "%lf\n", PS[i]);
	//	}

	//	// AniSol: 输出原始 6 分量
	//	if (!AniSol.empty()) {
	//			// 方式1：输出为 6 分量 field
	//		fprintf(outFile, "FIELD FieldData 1\n");
	//		fprintf(outFile, "AniSol 6 %ld double\n", V.size());
	//		for (int i = 0; i < (int)mesh.V.size(); i++) {
	//			fprintf(outFile, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf\n",
	//				AniSol[i][0], AniSol[i][1], AniSol[i][2],
	//				AniSol[i][3], AniSol[i][4], AniSol[i][5]);
	//		}

	//		// 方式2：输出为完整 3x3 tensor，便于 ParaView 使用
	//		fprintf(outFile, "TENSORS AniMetric double\n");
	//		for (int i = 0; i < (int)mesh.V.size(); i++) {
	//			const double m11 = AniSol[i][0];
	//			const double m12 = AniSol[i][1];
	//			const double m13 = AniSol[i][2];
	//			const double m22 = AniSol[i][3];
	//			const double m23 = AniSol[i][4];
	//			const double m33 = AniSol[i][5];

	//			fprintf(outFile, "%.16lf %.16lf %.16lf\n", m11, m12, m13);
	//			fprintf(outFile, "%.16lf %.16lf %.16lf\n", m12, m22, m23);
	//			fprintf(outFile, "%.16lf %.16lf %.16lf\n", m13, m23, m33);
	//		}
	//		
	//	}
	//}

	// -----------------------------
	// CELL_DATA
	// -----------------------------
	fprintf(outFile, "CELL_DATA  %ld\n", S.size() + F.size() + T.size());
	

	fprintf(outFile, "SCALARS GeoId float 1\n");
	fprintf(outFile, "LOOKUP_TABLE default\n");

	for (int i = 0; i < (int)S.size(); i++)
		fprintf(outFile, "%d\n", S[i][2]);
	for (int i = 0; i < (int)F.size(); i++)
		fprintf(outFile, "%d\n", F[i][3]);
	for (int i = 0; i < (int)T.size(); i++)
		fprintf(outFile, "%d\n", T[i][4]);
		

	fclose(outFile);

	return;
}

void DT::outUnRecvEdge(std::string filename) {
	FILE* outFile = fopen(filename.c_str(), "w");
	if (outFile == nullptr) {
		std::cout << "Write VTK file failed. - " << filename << "\n";
		return;
	}
	std::cout << "Writing mesh to - " << filename << "\n";
	fprintf(outFile, "# vtk DataFile Version 2.0\n");
	fprintf(outFile, "TetWild Mesh\n");
	fprintf(outFile, "ASCII\n");
	fprintf(outFile, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(outFile, "POINTS %d double\n", (int)Nodes.size());
	for (int i = 0; i < Nodes.size(); i++)
		fprintf(outFile, "%.16lf %.16lf %.16lf\n", Nodes[i].pt[0], Nodes[i].pt[1], Nodes[i].pt[2]);

	int noRecEdg = 0;
	for (int i = 0; i < SurEdgs.size(); i++)
		if (!isRecBndEdg(i))
			noRecEdg++;

	fprintf(outFile, "CELLS %ld %ld\n", noRecEdg, noRecEdg * 3);

	for (int i = 0; i < SurEdgs.size(); i++)
		if (!isRecBndEdg(i))
			fprintf(outFile, "2 %d %d\n", SurEdgs[i].iStart, SurEdgs[i].iEnd);

	fprintf(outFile, "CELL_TYPES  %ld\n", noRecEdg);

	for (int i = 0; i < noRecEdg; i++)
		fprintf(outFile, "3\n");

	fclose(outFile);
	return;
}
void DT::printSph_VTK(const std::vector<int> Sph, std::string filename) {
	int i, j;
	Mesh mesh;
	mesh.V.reserve(Nodes.size());
	mesh.T.reserve(Sph.size());
	for (i = 0; i < Nodes.size(); i++) {//have ghost
		mesh.V.push_back({ Nodes[i].pt[0],Nodes[i].pt[1],Nodes[i].pt[2] });
	}
	for (int i = 0; i < Sph.size(); i++) {
		if (Sph[i] == -1 || ishulltet(Sph[i]) || isDelEle(Sph[i]))
			continue;
		mesh.T.push_back({ Elems[Sph[i]].form[0],Elems[Sph[i]].form[1],Elems[Sph[i]].form[2],Elems[Sph[i]].form[3] });
	}
	dt::writeVTK(filename, mesh);
	return;
}

void DT::checkMeshError() {
#ifndef _DEBUG
	//Don't start in test
	return;
#endif
	int ia, ib, ic, id, a, b, c, d, i, j, k, minidx = -1;
	double minVolum = DBL_MAX;
	//outTempMesh("./temp.vtk");
	for (i = 0; i < Elems.size(); i++) {
		if (isDelEle(i))
			continue;
		for (j = 0; j < 4; j++) {
			if (Elems[i].form[j] < 0 || Elems[i].form[j] >= Nodes.size()) {
				meshLogger->error("Elems:{} has wrong Nod", i);
				meshLogger->error("Elems:{} ,info {}, {} {} {} {}", i, Elems[i].info, Elems[i].form[0], Elems[i].form[1], Elems[i].form[2], Elems[i].form[3]);
				throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
			}
			if (isDelNod(Elems[i].form[j])) {
				meshLogger->error("Elems:{} has Del Nod", i);
				throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
			}
			if (ghost != -1 && j != 3 && Elems[i].form[j] == ghost) {
				meshLogger->error("Elems:{} has wrong ghost", i);
				throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
			}
			int neig = getNeig(i, j);
			int neigOrd = getNeigOrd(i, j);
			if (neig == -1) {
				meshLogger->error("Elems:{} 's neig {} is -1 !", i, j);
				throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
			}
			if (Elems[i].form[j] == Elems[neig].form[neigOrd]) {
				meshLogger->error("Non-Convex! Elems:{} {}!", i, neig);
				meshLogger->error("{}: {} {} {} {}   {}: {} {} {} {}",
					i, Elems[i].form[0], Elems[i].form[1], Elems[i].form[2], Elems[i].form[3],
					neig, Elems[neig].form[0], Elems[neig].form[1], Elems[neig].form[2], Elems[neig].form[3]);
				throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
			}
			if (neig != -1 && isDelEle(neig)) {
				meshLogger->error("Elems:{} 's {} neig {} is delElm!", i, j, neig);
			}
			if (i != getNeig(neig, neigOrd) || j != getNeigOrd(neig, neigOrd)) {
				meshLogger->error("Adjacent between wrong:{} {} != {} {}", i, j, neig, neigOrd);
			}
			DNC(j, ia, ib, ic, id);
			b = Elems[i].form[ib];
			c = Elems[i].form[ic];
			d = Elems[i].form[id];
			int js = 0;
			for (k = 0; k < 4; k++) {
				if (Elems[neig].form[k] == b || Elems[neig].form[k] == c || Elems[neig].form[k] == d)
					js++;
			}
			if (js != 3) {
				meshLogger->error("Neighbor wrong:{} {}", i, neig);
			}
		}
		if (!ishulltet(i)&&!isvirtualtet(i)) {//check volume
			double* pa, * pb, * pc, * pd, ori;
			pa = Nodes[Elems[i].form[0]].pt;
			pb = Nodes[Elems[i].form[1]].pt;
			pc = Nodes[Elems[i].form[2]].pt;
			pd = Nodes[Elems[i].form[3]].pt;
			ori = dt::GEOM_FUNC::orient3d(pa, pb, pd, pc);
			if (ori < minVolum) {
				minVolum = ori;
				minidx = i;
			}
			if (ori <= 0) {
				meshLogger->error("Negtive Volume {}:{}", i, ori);
			}
		}
	}
	for (i = 0; i < Nodes.size(); i++) {
		if (isDelNod(i))
			continue;
		if (ghost != -1 && i == ghost)
			continue;
		if (Nodes[i].tet != 0 && isNod_in_Tet(i, Nodes[i].tet) == -1) {
			meshLogger->error("Error point {} to tet {}", i, Nodes[i].tet);
		}
	}
	for (i = 0; i < SurTris.size(); i++) {
		if (isDelSurTri(i) || SurTris[i].info > 1)
			continue;
		for (int j = 0; j < 3; j++) {
			int p1 = SurTris[i].form[j];
			int p2 = SurTris[i].form[(j + 1) % 3];
			int ee0 = BndEdg.get(p1, p2);
			if (ee0 == -1) {
				assert(0);
			}
		}
	}

	if (infolevel >= 2) {
		meshLogger->debug("Min Volume:{} {} : {},{},{},{}", minidx, minVolum,
			Elems[minidx].form[0], Elems[minidx].form[1], Elems[minidx].form[2], Elems[minidx].form[3]);
		//std::vector<int> worst = { minidx };
		//printSph(worst);
		//if (minVolum < 1e-30)
		//	outTempMesh("./temp.vtk");
	}
	return;
}

void DT::spdlogoutfile(bool outlogfile) {
	if (!outlogfile) {
        if (fileLogging) meshLogger.reset();
        fileLogging = false;
        return;
    }

	std::time_t now = std::time(nullptr);
	std::tm local_time{};

#ifdef _WIN32
	localtime_s(&local_time, &now);
#else
	localtime_r(&now, &local_time);
#endif

	std::ostringstream filename;
	filename << "dt_"
		<< std::put_time(&local_time, "%m_%d_%H_%M_%S")
        << "_" << std::hex << reinterpret_cast<uintptr_t>(this)
		<< ".log";

	// 终端输出
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

	// 文件输出，true 表示覆盖同名文件
	auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
		filename.str(),
		true
	);

	std::vector<spdlog::sink_ptr> sinks{ console_sink, file_sink };

	auto logger = std::make_shared<spdlog::logger>(
		"DT_logger",
		sinks.begin(),
		sinks.end()
	);

	meshLogger = logger;
    fileLogging = true;

	meshLogger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
	meshLogger->set_level(infolevel == 0 ? spdlog::level::err :
        infolevel == 1 ? spdlog::level::info : spdlog::level::debug);
	meshLogger->flush_on(spdlog::level::info);

	meshLogger->debug("Log file created: {}", filename.str());
}

void DT::checkEdgeDegree() {
	EdgeHasher<int> EdgeDegree;
	std::vector<int> Shell, shell_point, count(50);
	int ia, ib, pa, pb, sum = 0;

	for (int i = 0; i < Elems.size(); i++) {
		if (isDelEle(i))
			continue;
		for (int j = 0; j < 6; j++) {
			ia = Egid[j][0];
			ib = Egid[j][1];
			pa = Elems[i].form[ia];
			pb = Elems[i].form[ib];
			if (isBndEdg(pa, pb))
				continue;
			auto degreeEntry = EdgeDegree.try_emplace(pa, pb, 0);
			if (degreeEntry.second) {
				findShell(i, ia, ib, Shell, shell_point);
				*degreeEntry.first = static_cast<int>(Shell.size());
				count[Shell.size()]++;
			}
		}
	}
	for (int i = 3; i < 10; i++) {
		sum += count[i];
		meshLogger->debug("{}:{}", i, count[i]);
	}
	meshLogger->debug("sum:{}", sum);
	return;
}

int DT::checkEdgeLen() {
	int  pa, pb, maxEid = -1, numE = 0;
	double sumE = 0, maxE = 0, minE = DBL_MAX, dt;
	EdgeHasher<int> EdgeDegree;
	for (int i = 0; i < Elems.size(); i++) {
		if (isDelEle(i) || ishulltet(i))
			continue;
		for (int j = 0; j < 6; j++) {
			pa = Elems[i].form[Egid[j][0]];
			pb = Elems[i].form[Egid[j][1]];
			if (EdgeDegree.try_emplace(pa, pb, i * 6 + j).second) {
				double edgelen = distance(Nodes[pa].pt, Nodes[pb].pt);
				sumE += edgelen;
				maxE = std::max(maxE, edgelen);
				minE = std::min(minE, edgelen);
				numE++;
			}
		}
	}
	meshLogger->debug("TE_min  : {:.3f}  TE_max: {:.3f}  AvgE: {:.3f}   MaxE: {:.3f}   MinE: {:.3f}", minEdge, maxEdge, sumE / numE, maxE, minE);
	maxE = 0;
	for (int i = 0; i < SurEdgs.size(); i++) {
		if (isDelSurEdg(i))
			continue;

		pa = SurEdgs[i].iStart;
		pb = SurEdgs[i].iEnd;
		double edgelen = distance(Nodes[pa].pt, Nodes[pb].pt);
		if (edgelen > maxE) {
			maxE = edgelen;
			maxEid = i;
		}
	}
	return maxEid;
}
void DT::printSph(const std::vector<int> Sph) {
	int i, j;
	meshLogger->debug("Sph size: {}", Sph.size());
	for (i = 0; i < Sph.size(); i++) {
		meshLogger->debug("{}: {} {} {} {}", Sph[i], Elems[Sph[i]].form[0],
			Elems[Sph[i]].form[1], Elems[Sph[i]].form[2], Elems[Sph[i]].form[3]);
	}
	return;
}

uint64_t DT::getFreeMemory() {
#ifdef _WIN32
	MEMORYSTATUSEX memStatus;
	memStatus.dwLength = sizeof(memStatus);
	if (GlobalMemoryStatusEx(&memStatus)) {
		return static_cast<uint64_t>(memStatus.ullAvailPhys);
	}
	else {
		// 获取内存信息失败，可以记录日志或处理错误
		return 0;
	}
#else
	std::ifstream meminfo("/proc/meminfo");
	if (!meminfo.is_open()) {
		// 无法打开 /proc/meminfo，可能不是 Linux 系统
		return 0;
	}

	std::string line;
	uint64_t availableMemoryKB = 0;

	while (std::getline(meminfo, line)) {
		if (line.find("MemAvailable:") == 0) {
			// 提取数字部分
			size_t pos = line.find_first_of("0123456789");
			if (pos != std::string::npos) {
				availableMemoryKB = std::stoull(line.substr(pos));
			}
			break;
		}
	}

	meminfo.close();

	// 将 KB 转换为 Bytes
	return availableMemoryKB * 1024;
#endif
}

void DT::printMemoryUsage() {
#ifdef _WIN32
	// Windows 系统获取内存使用情况
	PROCESS_MEMORY_COUNTERS_EX pmc;
	if (GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc))) {
		SIZE_T usedMemory = pmc.WorkingSetSize;
		MEMORYSTATUSEX memStatus;
		memStatus.dwLength = sizeof(memStatus);
		if (GlobalMemoryStatusEx(&memStatus)) {
			SIZE_T totalMemory = memStatus.ullTotalPhys;
			SIZE_T freeMemory = memStatus.ullAvailPhys;

			double usedMemoryGB = usedMemory / (1024.0 * 1024 * 1024);
			meshLogger->debug("Used Memory : {:.3f} GB   Free Memory: {:.3f} GB", usedMemoryGB, freeMemory / (1024.0 * 1024 * 1024));

			if (usedMemoryGB > 0) {
				double elemsPerGB = Elems.size() / usedMemoryGB / 10000.0;
				meshLogger->debug("Elems Num   : {}   Elems per 1 GB: {:.3f} W", Elems.size(), elemsPerGB);
			}
		}
	}
#else
	// 获取当前程序使用的内存
	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	double usedMemoryGB = usage.ru_maxrss / 1024.0 / 1024.0; // ru_maxrss 以 KB 为单位

	// 获取系统的总内存和可用内存
	std::ifstream meminfo("/proc/meminfo");
	std::string line;
	long long totalMemory = 0, availableMemory = 0;

	while (std::getline(meminfo, line)) {
		if (line.find("MemTotal:") == 0) {
			totalMemory = std::stoll(line.substr(line.find_first_of("0123456789")));
		}
		else if (line.find("MemAvailable:") == 0) {
			availableMemory = std::stoll(line.substr(line.find_first_of("0123456789")));
			break;
		}
	}
	meminfo.close();

	double totalMemoryGB = totalMemory / 1024.0 / 1024.0;
	double freeMemoryGB = availableMemory / 1024.0 / 1024.0;

	meshLogger->debug("Program Used Memory: {:.3f} GB", usedMemoryGB);
	meshLogger->debug("Total System Memory: {:.3f} GB", totalMemoryGB);
	meshLogger->debug("Free System Memory: {:.3f} GB", freeMemoryGB);

	if (usedMemoryGB > 0) {
		double elemsPerGB = Elems.size() / usedMemoryGB / 10000.0;
		meshLogger->debug("Number of Elems: {}  Elems per 1 GB: {:.3f} W", Elems.size(), elemsPerGB);
	}
#endif
	return;
}

double DT::getPeakMegabytesUsed()
{
#ifdef _WIN32
	HANDLE hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, GetCurrentProcessId());
	if (NULL == hProcess) return 0;

	PROCESS_MEMORY_COUNTERS pmc;
	double mem = 0;
	if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc)))
	{
		mem = pmc.PeakWorkingSetSize / 1048576.0;
	}

	CloseHandle(hProcess);
	return mem;
#else
	return 0;
#endif
}
namespace dt {
MeshStageLog::MeshStageLog(DT& owner, const char* stage, int verbosity, MeshStageSummary report)
    : mesh(owner), name(stage), level(verbosity > 1 ? spdlog::level::debug : spdlog::level::info),
      active(owner.infolevel >= verbosity && owner.meshLogger->should_log(level)), summary(report) {
    if (!active) return;
    if (summary == MeshStageSummary::MeshChange) initialTets = countTets();
    mesh.meshLogger->log(level, "{} begin", name);
    started = std::chrono::steady_clock::now();
}

size_t MeshStageLog::countTets() const {
    size_t count = 0;
    for (int t = 0; t < static_cast<int>(mesh.Elems.size()); ++t)
        if (!mesh.isDelEle(t) && !mesh.isvirtualtet(t) && !mesh.ishulltet(t)) ++count;
    return count;
}

MeshStageLog::~MeshStageLog() noexcept { finish(); }

void MeshStageLog::finish(size_t finalTets) noexcept {
    if (!active) return;
    active = false;
    try {
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        if (summary == MeshStageSummary::Time) {
            mesh.meshLogger->log(level, "{}: time={:.6f}s", name, seconds);
            return;
        }
        if (finalTets == size_t(-1)) finalTets = countTets();
        const double speed = seconds > 0 ? static_cast<double>(finalTets) / seconds / 10000.0 : 0;
        if (summary == MeshStageSummary::MeshChange)
            mesh.meshLogger->log(level, "{}: time={:.6f}s tets={}->{} speed={:.3f} W/s",
                name, seconds, initialTets, finalTets, speed);
        else
            mesh.meshLogger->log(level, "{}: time={:.6f}s tets={} speed={:.3f} W/s",
                name, seconds, finalTets, speed);
    } catch (...) {
        // Diagnostics must not replace an algorithm exception during unwinding.
    }
}
}
