#include "dt.h"
#include "mesh_io.h"
#include <CLI11.hpp>
using namespace std;

#ifdef DT_EXEC
int main(int argc, char* argv[]) {
	std::string in_filename;
	std::string in_adp_filename;
	std::string in_per_filename;
	std::string out_filename;
	bool outwithsur = true;
	int adptype =0;
	dt::Args args;
	CLI::App app{ "RobustTetMeshing" };
	app.add_option("--input", in_filename, "Input surface mesh INPUT in .vtk format. (string, required)")->required();
	app.add_option("--fine_mesh_file", args.fine_mesh_file, "Input dense/reference mesh for boundary point projection. (string, optional)");
	app.add_option("--adpin", in_adp_filename, "Input Vertex. (string, required)");
	app.add_option("--perin", in_per_filename, "Input Periodic Vertex. (string, required)");
	app.add_option("--adptype", adptype, "Adp type,1 P_error,2 Tid,3 Coares. (int, required)");
	app.add_option("--output", out_filename, "Output tetmesh OUTPUT in .vtk format.  (string, required)");
	app.add_option("--constrain", args.constrain, "Keep boundary constrain. (int, optional, int:1)");
	app.add_option("--size", args.size, "Mesh size,The smaller the value, the larger the number. (double, optional, default:1)");
	app.add_option("--minEdge", args.minEdge, "Min Edge Length.(-1 is off) (double, optional, default:-1)");
	app.add_option("--maxEdge", args.maxEdge, "Max Edge Length.(-1 is off) (double, optional, default:-1)");
	app.add_option("--growsize", args.growsize, "Edge grow size.(-1 is off) (double, optional, default:>1)");
	app.add_option("--infolevel", args.infolevel, "Information out level. (int, optional, int:1)");
	app.add_option("--watertightcheck", args.watertightcheck, "Check Mesh Hole. (int, optional, int:0)");
	app.add_option("--dighole", args.hole, "Dig hole's coordinates. (double,double,double)");
	app.add_option("--diglayer", args.layer, "Dig layer idx,1  (int)");
	app.add_option("--digbnd", args.hole_bnd_vec, "Dig hole boundary idx,1  (int)");
	app.add_option("--refine", args.refine, "If refine mesh. (int, optional, int:1)");
	app.add_option("--optlevel", args.optlevel, "Mesh improvement level. (int, optional, int:2)");
	app.add_option("--optloop", args.optloop, "Mesh improvement loop. (int, optional, int:3)");
	app.add_option("--optanglestrict", args.optanglestrict, "Insertion threshold for minimum dihedral angle in degrees (0 disables insertion, default:20)");
	app.add_option("--optjacob", args.optTh, "Mesh improvement target opt jacobian. (double, optional, default:0.01)");
	app.add_option("--adpangle", args.adpangle, "Mesh Adptation Angle. (double, optional, default:160)");
	app.add_option("--nthread", args.nthread, "Maximum threads for every stage; active teams adapt to mesh size (-1: automatic).");
	app.add_option("--extrashell", args.extrashell, "Extraction Shell. (int, optional, int:-1)");
	app.add_option("--autoflip", args.autoflip, "use new auto flip. (int, optional, int:1)");
	app.add_option("--getmeshedg", args.getmeshedg, "get mesh bndedge by input Tris' geo id. (int, optional, int:1)");
	app.add_option("--ignoreIntersect", args.ignoreIntersect, "Ingore intersect input for autoGrid. (int, optional, int:0)");
	app.add_option("--outworsttet", args.outworsttet, "Out worst quality tet to VTK. (int, optional, int:0)");
	app.add_option("--outwithsur", outwithsur, "Out vtk don't with surface.(int, optional, int:1)");
	try {
		app.parse(argc, argv);
	}
	catch (const CLI::ParseError& e) {
		return app.exit(e);
	}
	args.filename = in_filename.substr(0, in_filename.length() - 4);
	if (out_filename.empty()) {
		out_filename = in_filename.substr(0, in_filename.length() - 4) + "_out.vtk";
	}

	dt::Mesh mesh;
	dt::DT d;
	for (int i = 2 ; i < 500; i++)
		args.layer.push_back(i);

	if (dt::readMesh(in_filename, mesh)) {
		if (in_per_filename.size() > 0) {
			readPeriodicP(in_per_filename, args.periodic_P);
		}

		if (adptype == 0) {
			if (d.tetrahedralize(mesh, args)){
				//for (auto it : args.periodic_P)
				//	printf("%d,", it);
				dt::writeVTK(out_filename, mesh, outwithsur);
				//mesh.T.resize(0);
				//out_filename = "./newin.vtk";
				//dt::writeVTK(out_filename, mesh, outwithsur);
			}
			else
				spdlog::info("tetrahedralize failed!");
		}
		else {
			int ret;
			if (adptype == 1) {
				std::vector<std::array<double, 4>> addVertex;
				readVertex(in_adp_filename, addVertex);
				ret = d.adaptation_by_pError(mesh, args, addVertex, 1.25);
			}
			else if (adptype == 2) {
				std::vector<int> refine_tet_id;// = { 1,2,3,4,5,6,7,8,9,10 };
				std::vector<int> refine_tri_id;// = { 137,35 };
				dt::readRefineT(in_adp_filename, refine_tet_id);
				//args.periodic_P = { 73,74,9,7,81,4,1,67,3,6,67,1,68,0,4,81,74,73,71,69,7,9,5,80,69,71,0,68,2,82,82,2,6,3,83,84,84,83,80,5,79,78,78,79,77,75,75,77,87,88,88,87,8,91,91,8,91,90,90,91,93,94,94,93 };
				args.adpangle = 160;
				args.optloop = 3;
				args.optTh =0.01;
				// Keep the user-provided minimum dihedral insertion angle.
				ret = d.adaptation_by_Tetid(mesh, args, refine_tri_id, refine_tet_id, 1.254);
			}
			else if (adptype == 3) {
				ret = d.adaptation_Coarse(mesh, args);
			}
			else if (adptype == 4) {
				dt::readNodesSize(in_adp_filename, mesh.pointSize);
				ret = d.adaptation_by_pSize(mesh, args);
			}
			else if (adptype == 5) {
				std::vector<std::array<double, 6>> anisol;
				dt::readAniSol(in_adp_filename, anisol);
				std::vector<int> lockF, lockV;
				ret = d.adaptation_by_ani(mesh, args, anisol, lockF, lockV);
			}
			else if (adptype == 6) {
				std::unordered_map<int, double> facetSize;// = { {1, 6.66}, {2, 3.66}, {3, 1.66}, {4, 10.66}, {5, 16.66}, };
				std::unordered_map<int, int> facetNum;// = { {1, 10000} };
				std::unordered_map<int, double> elementSize = { {85932,19.9842},{80800,19.9842},{81532,19.9842},{81662,19.9842},{98274,19.9842},{82332,19.9842},{110684,19.9842},{83184,19.9842},{97181,19.9842},{89818,19.9842},{85184,19.9842},{82901,19.9842},{90228,19.9842},{145720,19.9842},{82462,19.9842},{82771,19.9842},{83314,19.9842},{83444,19.9842},{97965,19.9842},{94263,19.9842},{111351,19.9842},{83828,19.9842},{84238,19.9842},{84368,19.9842},{102519,19.9842},{144687,19.9842},{89034,19.9842},{84498,19.9842},{110375,19.9842},{143517,9.52708},{95665,19.9842},{84628,19.9842},{84924,19.9842},{85054,19.9842},{93479,19.9842},{85506,19.9842},{85802,19.9842},{90898,19.9842},{86215,19.9842},{88452,19.9842},{104409,19.9842},{99006,19.9842},{88582,19.9842},{108267,19.9842},{88904,19.9842},{94933,19.9842},{96397,19.9842},{99738,19.9842},{99868,19.9842},{100652,19.9842},{100782,19.9842},{110967,19.9842},{101452,19.9842},{109295,19.9842},{102236,19.9842},{103303,19.9842},{109965,19.9842},{104087,19.9842},{105079,19.9842},{105749,19.9842},{106419,19.9842},{107151,19.9842},{107535,19.9842},{108999,19.9842},{111761,19.9842},{112145,19.9842},{112467,19.9842},{112763,19.9842},{113072,19.9842},{113118,9.52708},{113177,9.52708},{137177,9.52708},{138130,19.9842},{139020,19.9842},{140479,19.9842},{141193,19.9842},{142652,19.9842},{142782,19.9842},{142912,19.9842},{143042,19.9842},{143172,19.9842},{143302,19.9842},{143647,19.9842},{143777,19.9842},{143907,19.9842},{144037,19.9842},{144167,19.9842},{144297,19.9842},{144427,19.9842},{144557,19.9842},{144817,19.9842},{144947,19.9842},{160413,19.9842},{164693,9.52708} };

				std::unordered_map<int, int> elementNum;// = { {1641, 60} };
				//args.periodic_P = { 1406,22,1534,1536,1534,1535,1536,1534,1536,1537,0,23,0,1408,1420,1418,1420,8,1412,1416,23,0,23,1417,1416,1412,1418,1420,1418,30,9,28,1546,1542,1546,1549,28,9,1535,1537,1535,1534,1537,1535,1537,1536,1408,1417,1408,0,30,8,30,1418,1417,1408,1417,23,8,30,8,1420,1419,1414,1419,13,1414,1419,1414,34,1545,1539,1545,1544,19,36,36,19,1265,6,13,34,13,1419,1413,1405,1413,37,34,13,34,1414,20,37,20,1405,1538,1544,1538,1539,1544,1538,1544,1545,1539,1545,1539,1538,1540,1548,1540,1541,1548,1540,1548,1547,1541,1547,1541,1540,1,1410,1547,1541,1547,1548,1542,1546,1542,1543,1410,1,37,20,37,1413,1543,1549,1543,1542,1549,1543,1549,1546,875,737,737,875,1405,1413,1405,20,6,1265,1266,15,15,1266,22,1406 };
				
				//double portSize = 0.25;
				//facetSize[1] = 0.45;
				//facetSize[2] = 0.45;
				//facetSize[114731] = portSize;
				//facetSize[114709] = portSize;

				args.optloop = 10;
				args.optTh = 0.3;
				args.optanglestrict = 1;
				// Keep the user-provided minimum dihedral insertion angle.
				args.adpangle = 160;
				double lamasize = 19.9842;// 45.418557386999999;
				ret = d.adaptation_by_YunBoSzieControl(mesh, args, lamasize, facetSize, facetNum, elementSize, elementNum);
			}
			else if (adptype == 7) {
				ret = d.optimization(mesh, args);
			}

			if (ret == 1)
				dt::writeVTK(out_filename, mesh, outwithsur);
			else
				spdlog::info("Mesh Adaptation failed!");
		}
	}
	else
		spdlog::error("readVtk failed!");

	return 0;
}
#endif

#ifdef DT_LIBRARY
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
#endif