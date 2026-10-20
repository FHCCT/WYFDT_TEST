#include "dt.h"

DECL_VOLTET int API_Tetrahedralize(
	dt::Mesh& mesh,
	dt::Args& args)
{
	dt::DT d;
	return d.tetrahedralize(mesh, args);
}

DECL_VOLTET int API_Optimization(
	dt::Mesh& mesh, dt::Args& args
) {
	dt::DT d;
	return d.optimization(mesh, args);
}

DECL_VOLTET int API_adaptation_by_pSize(
	dt::Mesh& mesh, dt::Args& args
) {
	dt::DT d;
	return d.adaptation_by_pSize(mesh, args);
}

DECL_VOLTET int API_adaptation_by_pError(
	dt::Mesh& mesh,
	dt::Args& args,
	std::vector<std::array<double, 4>>& addVertex,
	double GrowRatio
)
{
	dt::DT d;
	return d.adaptation_by_pError(mesh, args, addVertex, GrowRatio);
}

DECL_VOLTET int API_adaptation_by_Tetid(
	dt::Mesh& mesh,
	dt::Args& args,
	std::vector<int> refine_tri_id,
	std::vector<int> refine_tet_id,
	double GrowRatio
)
{
	dt::DT d;
	return d.adaptation_by_Tetid(mesh, args, refine_tri_id ,refine_tet_id, GrowRatio);
}

DECL_VOLTET int API_adaptation_Coarse(
	dt::Mesh& mesh,
	dt::Args& args
)
{
	dt::DT d;
	return d.adaptation_Coarse(mesh, args);
}

DECL_VOLTET int API_adaptation_Anisotropic(
	dt::Mesh& mesh,
	dt::Args& args,
	std::vector<std::array<double, 6>>& anisol,
	std::vector<int> lockFactes,
	std::vector<int> lockVertex
) {
	dt::DT d;
	return d.adaptation_by_ani(mesh, args, anisol, lockFactes, lockVertex);
}

DECL_VOLTET int API_adaptation_YunBoSzieControl(
	dt::Mesh& mesh, dt::Args& args, double lamasize,
	std::unordered_map<int, double>& facetSize,
	std::unordered_map<int, int>& facetNum,
	std::unordered_map<int, double>& elementSize,
	std::unordered_map<int, int>& elementNum
) {
	dt::DT d;
	return d.adaptation_by_YunBoSzieControl(mesh, args, lamasize, facetSize,facetNum, elementSize, elementNum);
}

DECL_VOLTET int API_ReadMesh_File(
	std::string in_filename,
	dt::Mesh& mesh)
{
	if (!dt::readVTK(in_filename, mesh))
	{
		spdlog::error("readVtk failed!");
		return 0;
	}
	return 1;
}

DECL_VOLTET int API_WriteMesh(
	std::string out_filename,
	dt::Mesh& mesh, bool addSurTri)
{
	return dt::writeVTK(out_filename, mesh, addSurTri);
}
