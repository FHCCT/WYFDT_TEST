#include "dt.h"
#include "mesh_io.h"
#include <CLI11.hpp>
using namespace std;
//--input ../test\20260915_thread\adp/thread2_adpin.vtk --adpin ../test\20260915_thread\adp/thread2_refine_id.txt --adptype 2
//--input ../test\20260916_dtAdp\pcb2\in.vtk --adpin ../test\20260916_dtAdp\pcb2\in.txt --adptype 2
//--input ../test\Periodic\0928\14\in.vtk --adpin ..\test\Periodic\0928\14\in.txt --adptype 2
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
	app.add_option("--infolevel", args.infolevel, "Log detail: 0=errors only, 1=main stages/timing/quality, 2=all stage details (default:1)");
	app.add_option("--watertightcheck", args.watertightcheck, "Check Mesh Hole. (int, optional, int:0)");
	app.add_option("--dighole", args.hole, "Dig hole's coordinates. (double,double,double)");
	app.add_option("--diglayer", args.layer, "Dig layer idx,1  (int)");
	app.add_option("--digbnd", args.hole_bnd_vec, "Dig hole boundary idx,1  (int)");
	app.add_option("--refine", args.refine, "If refine mesh. (int, optional, int:1)");
	app.add_option("--optlevel", args.optlevel, "Mesh improvement level. (int, optional, int:2)");
	app.add_option("--optloop", args.optloop, "Mesh improvement rounds (default:5, adptype 2:18)");
	app.add_option("--optanglestrict", args.optanglestrict, "Insertion threshold for minimum dihedral angle in degrees (0 disables insertion; default:1, adptype 2:0)");
	app.add_option("--optTh", args.optTh, "Quality candidate threshold in the active metric (default:0.2, adptype 2:0.8)");
	app.add_option("--adpangle", args.adpangle, "Adaptation boundary angle (default:170, adptype 2:180)");
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
				ret = d.adaptation_by_pError(mesh, args, addVertex, 1.243);
			}
			else if (adptype == 2) {
				std::vector<int> refine_tet_id;// = { 1,2,3,4,5,6,7,8,9,10 };
				std::vector<int> refine_tri_id;// = { 137,35 };
				dt::readRefineT(in_adp_filename, refine_tet_id);
				//args.periodic_P = { 115,0,117,201,122,7,123,195,131,9,134,113,140,197,141,196,198,184,200,187,204,223,205,222,206,221,207,224,208,225,209,227,210,226,211,228,212,229,213,231,214,230,215,232,216,235,217,234,218,233,219,237,220,236,558,559,46,175,47,27,49,30,105,146,107,26,136,53,150,22,151,25,155,156,157,50,159,158,161,164,163,147,176,177,182,181,193,194,202,162,203,100,9,7,11,8,46,175,47,27,49,30,105,146,107,26,113,0,114,2,129,128,131,122,133,127,134,115,136,53,150,22,151,25,155,156,157,50,159,158,161,164,163,147,176,177,182,181,189,190,192,191,193,194,202,162,203,100,208,207,209,214,210,213,211,215,212,204,219,218,220,217,225,224,226,231,227,230,228,232,229,223,236,234,237,233 };
				// Shared isotropic quality profile; explicit CLI values take precedence.
				args.adpangle = 160;
				args.optloop =3;
				args.nthread = 4;
				args.optTh = 0.3;
				args.optanglestrict =0.01;
				ret = d.adaptation_by_Tetid(mesh, args, refine_tri_id, refine_tet_id, 1.243);
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
				std::unordered_map<int, double> elementSize = { {331613,2.94395},{226131,6.24505},{332699,2.94395},{239142,2.94395},{227512,6.24505},{505147,6.24505},{227675,6.24505},{507264,6.24505},{230627,6.24505},{240153,2.94395},{290912,2.94395},{229246,6.24505},{396422,6.24505},{331254,2.94395},{227865,6.24505},{443703,6.24505},{257547,6.24505},{368526,6.24505},{247198,6.24505},{244640,2.94395},{231188,6.24505},{240861,2.94395},{703186,6.24505},{232264,2.94395},{438249,6.24505},{230790,6.24505},{519522,6.24505},{345862,6.24505},{337702,6.24505},{231537,2.94395},{240946,2.94395},{231726,2.94395},{410066,6.24505},{245231,2.94395},{519902,6.24505},{445766,6.24505},{231915,2.94395},{444273,6.24505},{232453,2.94395},{232802,2.94395},{245303,2.94395},{679538,3.08421},{507101,6.24505},{244291,2.94395},{239416,2.94395},{238793,6.24505},{697847,6.24505},{276395,2.94395},{233151,2.94395},{244725,2.94395},{233500,2.94395},{234705,2.94395},{242418,2.94395},{235054,2.94395},{243753,2.94395},{239227,2.94395},{235970,2.94395},{235243,2.94395},{239879,2.94395},{235592,2.94395},{442507,6.24505},{235781,2.94395},{336398,2.94395},{236319,2.94395},{502103,6.24505},{237700,6.24505},{273882,6.24505},{238630,6.24505},{559990,6.24505},{239501,2.94395},{275404,2.94395},{239690,2.94395},{645719,6.24505},{257384,6.24505},{240068,2.94395},{446526,6.24505},{240502,2.94395},{671301,3.08421},{240587,2.94395},{240776,2.94395},{283089,2.94395},{241295,2.94395},{490650,6.24505},{241635,2.94395},{241380,2.94395},{241465,2.94395},{241550,2.94395},{241984,2.94395},{242333,2.94395},{244102,2.94395},{244810,2.94395},{370477,6.24505},{244882,2.94395},{364057,6.24505},{333690,2.94395},{246684,6.24505},{342546,6.24505},{287319,2.94395},{248579,6.24505},{282470,2.94395},{249960,6.24505},{256986,2.94395},{257737,6.24505},{259646,6.24505},{261027,6.24505},{699376,6.24505},{521232,6.24505},{262408,6.24505},{370097,6.24505},{263789,6.24505},{705148,6.24505},{702996,6.24505},{270408,2.94395},{271789,6.24505},{271979,6.24505},{272389,6.24505},{273770,6.24505},{295732,2.94395},{274384,6.24505},{274625,2.94395},{274814,2.94395},{449442,6.24505},{347551,6.24505},{275003,2.94395},{275593,2.94395},{275994,2.94395},{276796,2.94395},{399886,6.24505},{278001,2.94395},{461207,6.24505},{421206,6.24505},{278402,2.94395},{471376,6.24505},{278803,2.94395},{279044,2.94395},{279233,2.94395},{279634,2.94395},{519142,6.24505},{421396,6.24505},{280035,2.94395},{401457,6.24505},{361040,6.24505},{280498,6.24505},{347361,6.24505},{281879,6.24505},{282069,6.24505},{282711,2.94395},{514399,6.24505},{282900,2.94395},{283330,2.94395},{620212,2.94395},{283731,2.94395},{704958,6.24505},{468614,6.24505},{439630,6.24505},{283920,2.94395},{284109,2.94395},{368336,6.24505},{284510,2.94395},{284911,2.94395},{703349,6.24505},{285312,2.94395},{286517,2.94395},{295142,2.94395},{286918,2.94395},{405980,6.24505},{287720,2.94395},{421586,6.24505},{289101,6.24505},{484612,6.24505},{400076,6.24505},{290482,6.24505},{290723,2.94395},{291101,2.94395},{724625,6.24318},{443068,6.24505},{436868,6.24505},{291502,2.94395},{291691,2.94395},{292092,2.94395},{643792,2.94395},{504984,6.24505},{292493,2.94395},{559800,6.24505},{292894,2.94395},{294099,2.94395},{424897,6.24505},{294500,2.94395},{456874,6.24505},{294901,2.94395},{295331,2.94395},{296133,2.94395},{374973,6.24505},{303933,3.08421},{695931,6.24505},{304018,2.94395},{560180,6.24505},{323740,2.94395},{324141,2.94395},{330612,2.94395},{331013,2.94395},{518952,6.24505},{422343,6.24505},{331313,2.94395},{331554,2.94395},{331802,2.94395},{419831,6.24505},{331861,2.94395},{332050,2.94395},{628362,2.94395},{332239,2.94395},{374620,6.24505},{332640,2.94395},{333100,2.94395},{705637,6.24505},{407304,6.24505},{333289,2.94395},{660181,3.20364},{442005,6.24505},{334091,2.94395},{334492,2.94395},{335697,2.94395},{702806,6.24505},{336098,2.94395},{336339,2.94395},{336799,2.94395},{425060,6.24505},{337200,2.94395},{338100,6.24505},{338290,6.24505},{698363,6.24505},{370287,6.24505},{339671,6.24505},{341052,6.24505},{341242,6.24505},{520282,6.24505},{357337,6.24505},{343100,6.24505},{488755,6.24505},{344481,6.24505},{346052,6.24505},{346242,6.24505},{353663,6.24505},{346981,6.24505},{347171,6.24505},{468112,6.24505},{418937,2.94395},{349330,6.24505},{705474,6.24505},{518593,6.24505},{350711,6.24505},{352092,6.24505},{352282,6.24505},{435297,6.24505},{355044,6.24505},{355442,6.24505},{355956,6.24505},{357735,6.24505},{359659,6.24505},{357925,6.24505},{358088,6.24505},{463969,6.24505},{358278,6.24505},{361295,6.24505},{362676,6.24505},{365438,6.24505},{365940,6.24505},{398149,6.24505},{367321,6.24505},{367511,6.24505},{585876,3.20364},{503674,6.24505},{367701,6.24505},{367891,6.24505},{531029,6.24505},{368081,6.24505},{412318,6.24505},{369907,6.24505},{371858,6.24505},{490136,6.24505},{373239,6.24505},{374810,6.24505},{483231,6.24505},{445386,6.24505},{399696,6.24505},{376354,6.24505},{377735,6.24505},{378237,6.24505},{379618,6.24505},{380999,6.24505},{382380,6.24505},{393660,3.08421},{395041,6.24505},{396768,6.24505},{432535,6.24505},{398651,6.24505},{399153,6.24505},{399343,6.24505},{452568,6.24505},{399533,6.24505},{440551,6.24505},{401647,6.24505},{453949,6.24505},{433916,6.24505},{403028,6.24505},{403218,6.24505},{404599,6.24505},{519332,6.24505},{408685,6.24505},{410996,6.24505},{411510,6.24505},{540556,6.24505},{411920,6.24505},{698815,6.24505},{421841,6.24505},{419127,6.24505},{419641,6.24505},{420021,6.24505},{420951,6.24505},{422506,6.24505},{422916,6.24505},{423326,6.24505},{504586,6.24505},{424707,6.24505},{425250,6.24505},{425440,6.24505},{426821,6.24505},{427011,6.24505},{428392,6.24505},{520852,6.24505},{429773,6.24505},{431154,6.24505},{444653,6.24505},{436678,6.24505},{442905,6.24505},{487374,6.24505},{443323,6.24505},{446336,6.24505},{443513,6.24505},{443893,6.24505},{444083,6.24505},{491829,6.24505},{444463,6.24505},{444843,6.24505},{445033,6.24505},{445223,6.24505},{445576,6.24505},{445956,6.24505},{446146,6.24505},{701702,6.24505},{447681,6.24505},{447871,6.24505},{448061,6.24505},{450823,6.24505},{451325,6.24505},{451827,6.24505},{451990,6.24505},{452180,6.24505},{455330,6.24505},{455493,6.24505},{457064,6.24505},{458445,6.24505},{559610,6.24505},{459826,6.24505},{462588,6.24505},{465350,6.24505},{466731,6.24505},{469995,6.24505},{472757,6.24505},{472920,6.24505},{573592,2.94395},{485993,6.24505},{491164,6.24505},{491666,6.24505},{492906,6.24505},{501913,6.24505},{502293,6.24505},{504176,6.24505},{505557,6.24505},{506938,6.24505},{507427,6.24505},{514046,3.08421},{518443,6.24505},{514236,6.24505},{514562,6.24505},{514725,6.24505},{514888,6.24505},{516888,6.24505},{519712,6.24505},{520092,6.24505},{520472,6.24505},{520662,6.24505},{521042,6.24505},{550083,6.24505},{560370,6.24505},{594687,3.08421},{603422,3.08421},{612252,3.08421},{620402,6.24505},{637173,3.08421},{695605,2.94395},{695768,6.24505},{696094,6.24505},{698010,6.24505},{698200,6.24505},{699213,6.24505},{699539,6.24505},{701865,6.24505},{702055,6.24505},{702218,6.24505},{702408,6.24505},{703539,6.24505},{703702,6.24505},{703865,6.24505},{704028,6.24505},{705311,6.24505} };

				std::unordered_map<int, int> elementNum;// = { {1641, 60} };
				//args.periodic_P = { 1406,22,1534,1536,1534,1535,1536,1534,1536,1537,0,23,0,1408,1420,1418,1420,8,1412,1416,23,0,23,1417,1416,1412,1418,1420,1418,30,9,28,1546,1542,1546,1549,28,9,1535,1537,1535,1534,1537,1535,1537,1536,1408,1417,1408,0,30,8,30,1418,1417,1408,1417,23,8,30,8,1420,1419,1414,1419,13,1414,1419,1414,34,1545,1539,1545,1544,19,36,36,19,1265,6,13,34,13,1419,1413,1405,1413,37,34,13,34,1414,20,37,20,1405,1538,1544,1538,1539,1544,1538,1544,1545,1539,1545,1539,1538,1540,1548,1540,1541,1548,1540,1548,1547,1541,1547,1541,1540,1,1410,1547,1541,1547,1548,1542,1546,1542,1543,1410,1,37,20,37,1413,1543,1549,1543,1542,1549,1543,1549,1546,875,737,737,875,1405,1413,1405,20,6,1265,1266,15,15,1266,22,1406 };
				
				//double portSize = 0.25;
				//facetSize[1] = 0.45;
				//facetSize[2] = 0.45;
				//facetSize[114731] = portSize;
				//facetSize[114709] = portSize;

				args.optloop = 10;
				args.optTh = 0.3;
				args.optanglestrict = 0.01;
				// Keep the user-provided minimum dihedral insertion angle.
				args.adpangle = 170;
				constexpr double lamasize = 33.3334;
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
#endif

#if defined(DT_LIBRARY) && !defined(DT_EXEC)
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