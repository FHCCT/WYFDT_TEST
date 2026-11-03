#include "dt.h"
#include "mesh_io.h"
#include <CLI11.hpp>
#include "dt_sizecontrol_profile.h"
using namespace std;
//--input ../test\20260915_thread\adp/thread2_adpin.vtk --adpin ../test\20260915_thread\adp/thread2_refine_id.txt --adptype 2
//--input ../test\20260916_dtAdp\pcb2\in.vtk --adpin ../test\20260916_dtAdp\pcb2\in.txt --adptype 2
//--input ../test\Periodic\0928\14\in.vtk --adpin ..\test\Periodic\0928\14\in.txt --adptype 2
int main(int argc, char* argv[]) {
	std::string in_filename;
	std::string in_adp_filename;
	std::string in_per_filename;
	std::string out_filename;
    std::string sizeControlProfile;
	bool outwithsur = true;
	int adptype =0;
	dt::Args args;
	CLI::App app{ "RobustTetMeshing" };
	app.add_option("--input", in_filename, "Input surface mesh INPUT in .vtk format. (string, required)")->required();
	app.add_option("--fine_mesh_file", args.fine_mesh_file, "Input dense/reference mesh for boundary point projection. (string, optional)");
	app.add_option("--adpin", in_adp_filename, "Input Vertex. (string, required)");
	app.add_option("--perin", in_per_filename, "Input Periodic Vertex. (string, required)");
    app.add_option("--sizecontrol-profile", sizeControlProfile, "Replay size-control map and branch parameters (adptype 6)");
	app.add_option("--adptype", adptype, "Adp type,1 P_error,2 Tid,3 Coares. (int, required)");
	app.add_option("--output", out_filename, "Output tetmesh OUTPUT in .vtk format.  (string, required)");
	app.add_option("--constrain", args.constrain, "Keep boundary constrain. (int, optional, int:1)");
	app.add_option("--size", args.size, "Mesh size,The smaller the value, the larger the number. (double, optional, default:1)");
	app.add_option("--minEdge", args.minEdge, "Min Edge Length.(-1 is off) (double, optional, default:-1)");
	app.add_option("--maxEdge", args.maxEdge, "Max Edge Length.(-1 is off) (double, optional, default:-1)");
	app.add_option("--growsize", args.growsize, "Edge grow size.(-1 is off) (double, optional, default:>1)");
	app.add_option("--infolevel", args.infolevel, "Log detail: 0=errors only, 1=main stages/timing/quality, 2=all stage details (default:1)");
	app.add_option("--watertightcheck", args.watertightcheck, "Check Mesh Hole. (int, optional, int:0)");
	app.add_option("--dighole", args.hole, "Dig hole's coordinates. (double,double,double)");
	app.add_option("--diglayer", args.layer, "Dig layer idx,1  (int)");
	app.add_option("--digbnd", args.hole_bnd_vec, "Dig hole boundary idx,1  (int)");
	app.add_option("--refine", args.refine, "If refine mesh. (int, optional, int:1)");
	app.add_option("--optlevel", args.optlevel, "Mesh improvement level. (int, optional, int:2)");
	app.add_option("--optloop", args.optloop, "Mesh improvement rounds (default:5, adptype 2:3)");
	app.add_option("--optanglestrict", args.optanglestrict, "Minimum dihedral insertion threshold in degrees (default:1). optlevel 5: >0 enables Q<0.01 tail repair; 0 keeps pure flips/smoothing. Other levels: 0 disables only the angle trigger (SJ<1e-5 still triggers).");
	app.add_option("--optTh", args.optTh, "Worst-quality fraction per optimization round, in [0,1] (default:0.2, adptype 2:0.5)")->check(CLI::Range(0.0, 1.0));
	app.add_option("--adpangle", args.adpangle, "Adaptation boundary angle (default:170, adptype 2:160)");
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
    d.infolevel = std::max(0, std::min(2, args.infolevel));
    d.meshLogger->set_level(d.infolevel == 0 ? spdlog::level::err :
        d.infolevel == 1 ? spdlog::level::info : spdlog::level::debug);
	//for (int i = 2 ; i < 500; i++)
	//	args.layer.push_back(i);

    int readResult;
    { dt::MeshStageLog readLog(d, "Read input"); readResult = dt::readMesh(in_filename, mesh); }
	if (readResult) {
		if (in_per_filename.size() > 0) {
			readPeriodicP(in_per_filename, args.periodic_P);
		}

		if (adptype == 0) {
			if (d.tetrahedralize(mesh, args)) {
				if (d.infolevel > 0) d.meshLogger->info("Writing mesh to: {}", out_filename);
				dt::writeVTK(out_filename, mesh, outwithsur);
			}
			else
				{ d.meshLogger->error("tetrahedralize failed!"); return 1; }
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
                if (!dt::readRefineT(in_adp_filename, refine_tri_id, refine_tet_id)) return 1;
                if (std::any_of(refine_tri_id.begin(), refine_tri_id.end(), [&](int id) { return id >= mesh.F.size(); }) ||
                    std::any_of(refine_tet_id.begin(), refine_tet_id.end(), [&](int id) { return id >= mesh.T.size(); })) {
                    d.meshLogger->error("Refinement ID is outside the input mesh range");
                    return 1;
                }
                if (d.infolevel > 0) d.meshLogger->info("Refinement selection: triangles={} tets={}",
                    refine_tri_id.size(), refine_tet_id.size());
				//args.periodic_P = { 115,0,117,201,122,7,123,195,131,9,134,113,140,197,141,196,198,184,200,187,204,223,205,222,206,221,207,224,208,225,209,227,210,226,211,228,212,229,213,231,214,230,215,232,216,235,217,234,218,233,219,237,220,236,558,559,46,175,47,27,49,30,105,146,107,26,136,53,150,22,151,25,155,156,157,50,159,158,161,164,163,147,176,177,182,181,193,194,202,162,203,100,9,7,11,8,46,175,47,27,49,30,105,146,107,26,113,0,114,2,129,128,131,122,133,127,134,115,136,53,150,22,151,25,155,156,157,50,159,158,161,164,163,147,176,177,182,181,189,190,192,191,193,194,202,162,203,100,208,207,209,214,210,213,211,215,212,204,219,218,220,217,225,224,226,231,227,230,228,232,229,223,236,234,237,233 };
				// Shared isotropic quality profile; explicit CLI values take precedence.
				if (app.get_option("--adpangle")->count() == 0) args.adpangle = 160;
				if (app.get_option("--optloop")->count() == 0) args.optloop = 10;
				if (app.get_option("--nthread")->count() == 0) args.nthread = 4;
				if (app.get_option("--optTh")->count() == 0) args.optTh = 0.5;
				if (app.get_option("--optanglestrict")->count() == 0) args.optanglestrict = 1;
				ret = d.adaptation_by_Tetid(mesh, args, refine_tri_id, refine_tet_id, 1.268);
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
				std::unordered_map<int, double> elementSize = { {173124,4.99604},{295413,4.99604},{309565,4.99604},{170801,2.83756},{183943,4.99604},{191332,4.99604},{253008,4.99604},{173052,4.99604},{181614,2.83756},{249383,4.99604},{191222,4.99604},{215787,4.99604},{199468,4.99604},{188891,4.99604},{188169,4.99604},{196990,4.99604},{192747,4.99604},{250519,4.99604},{299078,4.99604},{194525,4.99604},{195264,4.99604},{198716,4.99604},{199578,4.99604},{265313,4.99604},{249540,4.99604},{213265,2.83756},{215043,4.99604},{243544,4.99604},{246618,4.99604},{272401,4.99604},{304643,4.99604},{246715,4.99604},{247480,4.99604},{249455,4.99604},{250604,4.99604},{251583,4.99604},{252234,4.99604},{310563,4.99604},{252344,4.99604},{256581,4.99604},{297904,4.99604},{257289,4.99604},{274551,4.99604},{257997,4.99604},{258692,4.99604},{259387,4.99604},{260750,2.83756},{310871,4.99604},{272304,4.99604},{261485,4.99604},{261557,4.99604},{262261,4.99604},{266214,4.99604},{264512,4.99604},{266129,4.99604},{265216,4.99604},{265385,4.99604},{295303,4.99604},{267193,4.99604},{269100,4.99604},{274648,4.99604},{269878,4.99604},{272207,4.99604},{272511,4.99604},{274344,4.99604},{274454,4.99604},{274720,2.83756},{293470,2.83756},{295498,4.99604},{296477,4.99604},{296562,4.99604},{297541,4.99604},{297638,4.99604},{297735,4.99604},{309434,4.99604},{297807,4.99604},{309768,4.99604},{298014,4.99604},{298099,4.99604},{299163,4.99604},{300142,4.99604},{300227,4.99604},{301206,4.99604},{308465,4.99604},{308850,4.99604},{308922,4.99604},{308994,4.99604},{309066,4.99604},{309138,4.99604},{309210,4.99604},{309282,4.99604},{309358,4.99604},{309506,4.99604},{309637,4.99604},{309696,4.99604},{309840,4.99604},{309916,4.99604},{309992,4.99604},{310068,4.99604},{310144,4.99604},{310216,4.99604},{310275,4.99604},{310347,4.99604},{310419,4.99604},{310491,4.99604},{310635,4.99604},{310694,4.99604},{310753,4.99604},{310812,4.99604},{310943,4.99604},{311002,4.99604},{311424,4.99604} };

				std::unordered_map<int, int> elementNum;// = { {1641, 60} };
				//args.periodic_P = { 1406,22,1534,1536,1534,1535,1536,1534,1536,1537,0,23,0,1408,1420,1418,1420,8,1412,1416,23,0,23,1417,1416,1412,1418,1420,1418,30,9,28,1546,1542,1546,1549,28,9,1535,1537,1535,1534,1537,1535,1537,1536,1408,1417,1408,0,30,8,30,1418,1417,1408,1417,23,8,30,8,1420,1419,1414,1419,13,1414,1419,1414,34,1545,1539,1545,1544,19,36,36,19,1265,6,13,34,13,1419,1413,1405,1413,37,34,13,34,1414,20,37,20,1405,1538,1544,1538,1539,1544,1538,1544,1545,1539,1545,1539,1538,1540,1548,1540,1541,1548,1540,1548,1547,1541,1547,1541,1540,1,1410,1547,1541,1547,1548,1542,1546,1542,1543,1410,1,37,20,37,1413,1543,1549,1543,1542,1549,1543,1549,1546,875,737,737,875,1405,1413,1405,20,6,1265,1266,15,15,1266,22,1406 };
				
				//double portSize = 0.25;
				//facetSize[1] = 0.45;
				//facetSize[2] = 0.45;
				//facetSize[114731] = portSize;
				//facetSize[114709] = portSize;

				args.optloop = 10;
				args.optTh = 0.4;
				args.optanglestrict = 0.01;
				// Keep the user-provided minimum dihedral insertion angle.
				args.adpangle = 170;
				double lamasize = 14.9896;
                if (!sizeControlProfile.empty()) {
                    try { dt::readSizeControlProfile(sizeControlProfile,args,lamasize,facetSize,facetNum,elementSize,elementNum); }
                    catch (const std::exception& error) { d.meshLogger->error("{}",error.what()); return 1; }
                }
				ret = d.adaptation_by_YunBoSzieControl(mesh, args, lamasize, facetSize, facetNum, elementSize, elementNum);
			}
			else if (adptype == 7) {
				ret = d.optimization(mesh, args);
			}

			if (ret == 1) {
				if (d.infolevel > 0) d.meshLogger->info("Writing mesh to: {}", out_filename);
				dt::writeVTK(out_filename, mesh, outwithsur);
			}
			else
				spdlog::error("Mesh Adaptation failed!");
		}
	}
	else
		spdlog::error("readVtk failed!");

	return 0;
}
