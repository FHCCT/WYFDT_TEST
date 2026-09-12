#ifndef dt_API_h
#define dt_API_h

#ifdef __linux__
#define DECL_VOLTET
#elif _WIN32
#define DECL_VOLTET __declspec(dllexport)
#else
#error "Unsupported platform"
#endif

#include <functional>
#include <vector>
#include <array>
#include <string>
#include <set>
#include <map>
#include <unordered_map>

namespace dt {
	struct Mesh {//for io easily
		std::vector<std::array<double, 3>> V;
		std::vector<std::array<int, 3>> S;
		std::vector<std::array<int, 4>> F;
		std::vector<std::array<int, 5>> T;
		std::vector<double> pointSize;
	};

	struct Args {
		int infolevel = 1;			//输出信息等级，设置为0输出最少，作为库，建议设置为0
		int constrain = 1;			//0，conforming;1,constrain
		int refine = 1;				//是否开启网格细化，默认打开
		int optlevel = 3;			//优化等级，0~4，0不优化，默认最高4，质量最慢，实测3更好插点
		int optloop = 5;			//优化轮数，默认8
		int outworsttet = 0;		//打印最差的四面体单元，用于debug
		int nthread = 4; // Thread upper bound for all stages; <=0 uses available CPUs.
		int ignoreIntersect = 0;	//Ignore intersect input, may result in non-watertight
		int watertightcheck = 0;	//Find non-mainfold edge，检测水密
		int extrashell = 0;
		int autoflip = 1;			//use new auto flip
		int getmeshedg = 0;	
		int outlogfile = 0;	
		double size = 1;			//网格过渡尺寸，默认1，越小越密
		double minEdge = -1;		//最小网格边长，默认不启动（-1），可以主动设置
		double maxEdge = -1;		//最大网格边长，默认不启动（-1），可以主动设置
		double growsize = 1.01;     //边界增长率，默认1.01，大于1
		double optanglestrict = 1; // Minimum dihedral angle in degrees for insertion; 0 disables insertion.
		double optTh = 0.2;		    //通用优化目标
		double adpangle = 170;
		std::string filename;
		std::string fine_mesh_file;		// optional dense/reference mesh for boundary point projection
		std::vector<double> hole;	//挖洞的坐标，可以是连续的x1,y1,z1;x2,y2,z2;... ...
		std::vector<int> layer;		//挖洞的层，适用于流场.
		std::vector<int> hole_bnd_vec; //挖洞的边界面ID
		std::map< int, std::set<int>> bnd_bodyid;
		std::vector<int> periodic_P;
		std::function<double(const double&, const double&, const double&)> sizingFunc = nullptr; //尺寸场接口
	};
}

#ifdef DT_LIBRARY
/*
 *  @brief	Main Tetrahedralize
 *
 *  @param[in]		mesh		input mesh
 *  @param[in]		args		input arguments
 *  @param[out]		mesh		output mesh
 *  @return
 *	1				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_Tetrahedralize(
	dt::Mesh& mesh, dt::Args& args
);

/*
 *  @brief	Main optimization
 *
 *  @param[in]		mesh		input mesh
 *  @param[in]		args		input arguments
 *  @param[out]		mesh		output mesh
 *  @return
 *	1				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_Optimization(
	dt::Mesh& mesh, dt::Args& args
);

/*
 *  @brief	Main Adaptation
 *
 *  @param[in]		mesh		input mesh
 *  @param[in]		args		input arguments
 *  @param[out]		mesh		output mesh
 *  @return
 *	1				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_adaptation_by_pSize(
	dt::Mesh& mesh, dt::Args& args
);

DECL_VOLTET int API_adaptation_by_pError(
	dt::Mesh& mesh, dt::Args& args, std::vector<std::array<double, 4>>& addVertex, double GrowRatio
);

DECL_VOLTET int API_adaptation_by_Tetid(
	dt::Mesh& mesh, dt::Args& args, std::vector<int> refine_tri_id, std::vector<int> refine_tet_id, double GrowRatio
);

DECL_VOLTET int API_adaptation_Coarse(
	dt::Mesh& mesh, dt::Args& args
);

DECL_VOLTET int API_adaptation_Anisotropic(
	dt::Mesh& mesh,
	dt::Args& args,
	std::vector<std::array<double, 6>>& anisol,
	std::vector<int> lockFactes,
	std::vector<int> lockVertex
);

DECL_VOLTET int API_adaptation_YunBoSzieControl(
	dt::Mesh& mesh, dt::Args& args,
	double lamasize,
	std::unordered_map<int, double>& facetSize,
	std::unordered_map<int, int>& facetNum,
	std::unordered_map<int, double>& elementSize,
	std::unordered_map<int, int>& elementNum
);

/*
 *  @brief	Create Mesh from vtk file
 *
 *  @param[in]		filename		input filename
 *  @param[out]		mesh			output mesh
 *  @return
 *	1				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_ReadMesh_File(
	std::string in_filename, dt::Mesh& mesh
);

/*
 *  @brief	Create Mesh from vtk file
 *
 *  @param[in]		filename		output filename
 *  @param[in]		mesh			output mesh
 */
DECL_VOLTET int API_WriteMesh(
	std::string out_filename, dt::Mesh& mesh, bool addSurTri = false
);
#endif
#endif
