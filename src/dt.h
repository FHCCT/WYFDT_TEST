#pragma once
#ifndef dt_h
#define dt_h

#include <set>
#include <map>
#include <array>
#include <cmath>
#include <string>
#include <mutex>
#include <cfloat>
#include <queue>
#include <bitset>
#include <chrono>
#include <random>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "dt_API.h"
#include "spdlog/spdlog.h"
#include "dt_define.h"
#include "dt_hash.h"
#include "geom_func.h"
#include "mesh_io.h"
#include "../extern/tree/dt_binary_tree.hpp"
#include "dt_mem.h"
#include "dt_bw.h"
#include "dt_bw_parallel.h"
#include "dt_parallel.h"
#include "dt_projection.h"
#include "dt_logger.h"

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <unistd.h> // 用于 sysconf
#include <sys/resource.h> // 用于 getrusage
#endif

using namespace dt;
//============================================================================//
//
// dt
//
// Tetrahedron-based mesh data structure.
// Elementary flip operations
// Boundary recovery
// Control Steiner points
// Mesh improvement
//============================================================================//
namespace dt {
    // Distinct instances may run concurrently. Mutating operations on one DT
    // require exclusive ownership; its parallel entry points manage worker access.
	class DT {
	public:
		DT();
        InstanceLogger meshLogger;
		~DT();

		/****************************  Data ****************************/
#ifdef MEMORY_POOL
		SmallVector<Node> Nodes;
		SmallVector<Elem> Elems;
		SmallVector<double> qual;
#else
		std::vector<Node> Nodes;
		std::vector<Elem> Elems;
#endif // DEBUG
		std::vector<SurTri> SurTris;
		std::vector<SurEdg> SurEdgs;
		std::queue<int> Evacancy;
		std::queue<int> Evacancy_thread[128];
		std::queue<int> Nvacancy;
		std::vector<int> SteinerOrd;
		std::vector<std::pair<int, int>> EdgSteiner;//point idx;parent edge idx
		std::vector<std::pair<int, int>> TriSteiner;//point idx;parent face idx
		std::vector<std::array<double, 6>> AniSol;
		std::unordered_set<int> lockF;
		std::unordered_set<int> lockV;
		std::unordered_set<int> lockE;
		
		bool QuantityControl = false;
		std::unordered_map<int, int> GeoNum;
		std::unordered_map<int, int> FacetNum;
		std::vector<std::array<int, 3>> flipnmRecll;
        std::array<std::vector<std::array<int, 3>>, 128> parallelFlipHistory;
        std::array<int, 128> parallelFlipCount{};
        bool parallelTopologyBatch = false; // No container growth inside a topology batch.

		std::unordered_map<int, std::vector<int>> periodic_P;
		enum { SUS_METRIC = 9 }; // Signed Knupp mean-ratio quality; SUS regularized smoothing.
		
		int ghost;
		int nSurNodes;
		int nSurTris;
		int seg[2];//some global information for boundary recover
		int fac[3];
		int fliplevel;
		int fliplevel_face;
        struct SusIdleState { std::vector<double> neighborhood; int skips = 0; };
        std::vector<SusIdleState> susIdleStates;
        bool susReuseIdle = true;
		int improve_Metric;
		int num_threads = 1; // Single thread-count upper bound for every stage.
        size_t parallelMinPointsPerThread = 8192;
        bool threadsInitialized = false;
		int virtualID = -2;
		int tempfliptime = 0;
		int maxfliptime;
		int optflipdeep = 3;
		bool improve_step;
		bool ignoreIntersect = false;
		bool modifyBnd = false;
		bool addBoxFlag = false;
		double AvgEdgLen = 0;
		double meshSize = 1;
		double cos_collinear_ang_tol;
		double maxW[3], minW[3], dist_max;
		double minEdge = -1;
		double maxEdge = -1;
		double growsize = -1;
		int initBWshell = 0;
		int addst = 0;
		int addstbnd = 0;
		double minVolume_bw = 0;
		double ani_upper = 2.0;
		double ani_lower = 0.3;

		//hash for quakily check bnd
		TriHasher<int> BndTri;
		EdgeHasher<int> BndEdg;
		
		//int mp4time = 0;

		//for boundary size transition
		std::vector<int> surTri_mapping;
		/***************** user interaction parameters *****************/
		int infolevel;
		/**************************  process **************************/
		int dt_init(Mesh& mesh, Args& args);
		int tetrahedralize(Mesh& mesh, Args& args);
		int BndPntInst(Mesh& mesh, Args& args);
		int BoudaryRecover(Mesh& mesh, Args& args);
		int ColorVirtualTet(Args& args);
		int MeshRefine(Args& args);
		int MeshImprove(Args& args);
		int RemoveTet(Args& args);
		int outMesh(Mesh& mesh, Args& args);
		int optimization(Mesh& mesh, Args& args);
		int adaptation_by_pError(Mesh& mesh, Args& args, std::vector<std::array<double, 4>>& addVertex, double GrowRatio);
		int adaptation_by_Tetid(Mesh& mesh, Args& args, std::vector<int> refine_tri_id, std::vector<int> refine_tet_id, double GrowRatio);
		int adaptation_by_pSize(Mesh& mesh, Args& args);
		int adaptation_by_YunBoSzieControl(Mesh& mesh, Args& args, double lamasize, std::unordered_map<int, double>& facetSize,std::unordered_map<int, int>& facetNum,std::unordered_map<int, double>& elementSize, std::unordered_map<int, int>& elementNum);
		int adaptation_Coarse(Mesh& mesh, Args& args);
		int adaptation_by_ani(Mesh& mesh, Args& args, std::vector<std::array<double, 6>>& anisol, std::vector<int> lockFactes, std::vector<int> lockVertex);
		/************************** algorithm ************************/
		void buildPntInfo(Mesh& mesh);
		void buildTetInfo(Mesh& mesh, Args& args);
		int BW_insert_vertex(int iNod, std::vector<int>& srchtet, int info, int thread_n = -1);
		int locate_pnt(int iNod, int& searchtet);
		void AddBox(double scaled);
		int findShell(const int t, const  int p2, const int p3, std::vector<int>& Shell, std::vector<int>& shell_point, int thread_n = -1);
		int findSphere(const  int p, std::vector<int>& Sphere);
		int findSphere_pnt(const  int p, std::unordered_set<int>& Sphere_pnt);
		int findSphere_global(const  int p, std::vector<int>& Sphere);
		int findSphere_tri(const  int p, std::unordered_set<int>& Sphere_tri);
		//int findSphere_tri_p(const  int p, std::unordered_set<int>& neig_p);
		//boundary recover
		void buildBndInfo(Mesh& mesh, Args& args, bool buildSize = true);
		void SurMeshClean(Mesh& mesh, Args& args);
		void AttachSeg2Pnt(int ie);
		void DelSingleEdge(int ie);
		int AttachPnt2Seg(int iNod, int targetE);
		bool isMeshEdge(const int p1, const int p2, int* tet = NULL);
		bool isMeshFace(const int p1, const int p2, const int p3, int* tet = NULL);
		int AutorecoverEdges(Args& args);
		void updateFliptype(std::unordered_map<int, int>& N_lostE_pre, std::queue<int>& lost);
		int recoverEdgesPass(Args& args);
		int recoverFacesPass(Args& args);
		int removeStPass(Args& args);
		int recoverEdges(std::queue<int>& lost, int fullsearch, int info);
		int recoverFaces(std::queue<int>& lost, int info);
		int recoverEdge(const int targetE, int fullsearch, int info);
		int findIntersectwithEdgs(const int targetE, std::vector<std::array<int, 3>>& Vid);
		int recoverEdgebyFlip(const int targetE, int dirflag, int info);
		int recoverFace(const int targetE, int info);
		int recoverFacebyFlip_Split(const int targetF, int info);
		int recoverFacebyaddinSt(const int targetF, std::vector<int>& newN, int info);
		int addinnerSteiner_Edge(const int targetE, std::vector<int>& newN, int info);
		int addinnerSteiner_Edge2(const int targetE, std::vector<int>& newN, int info);
		int addinnerSteiner_Face(const int targetF, std::vector<int>& newN, double dis[2]);
		int removeEdgStiner(const int idx, int level);
		int removeTriStiner(const int idx);
		int flipBndEdgPass(int nloop);
		int flipBndEdge(int index);
		int splitBndEdge(const int targetE, int info);
		int splitBndTri(const int targetF, std::vector<int>& shell, double intPnt[], int info);
		int removePnt(const int iNod,int tryTime=10);
		int ignoreE(int targetE);
		int ignoreF(int targetF);
		int finddirection(const int p1, const int p2, int& srctet);
		int finddirection_global(const int p1, const int p2, int& srctet);
		int DealIntersect(int pa, int pb, int p1, int p2);
		//main flip
		int removeface(std::vector<int>& oldtet, int ia, int info = 0, int thread_n = -1);
		int removeEdge(std::vector<int>& oldtet, int ia, int ib, int info = 0, int thread_n = -1);
		int flipnm(std::vector<int>& oldtet, const int ia, const int ib, int level, int maxlevel, double minq, int thread_n = -1);
		int flip32(std::vector<int>& oldtet, int ia, int ib, int thread_n = -1);
		int flip23(std::vector<int>& oldtet, int a, int thread_n = -1);
		int flip41(std::vector<int>& oldtet, int iNod);
		bool flipintersectcheck(int fliptype, int a, int b, int c, int d, int e);
		//MeshRefine
		void updateSize(int iNod, Args& args);
		void buildspace(Mesh& mesh, Args& args);
		int creatNewV_Grav(int i, Args& args);
		int ColorTets();
		int shellextra();
		int ColorTetNeig(int iElm, int color);
        int colorTetComponent(int iElm, int color, std::vector<unsigned>& marks, unsigned stamp, std::vector<int>& queue);
		/****************** Mesh Improvement **************/
		int improve_init(Args& args);
		int TraditionalOptPass(Args& args);
		int OrthogonalityOptPass(Args& args);
		int QuicklyOptPass(Args& args);
		int QuicklyOptCoarsePass(Args& args);
		int OptSizeControl_yunbo(double lamasize, std::unordered_map<int, double>& facetSize, std::unordered_map<int, int>& facetNum, std::unordered_map<int, double>& elementSize, std::unordered_map<int, int>& elementNum, Args& args);
		int OptbyAnisotropicPass(Args& args);
		int VolumeImprovePass(Args& args);
		int prepareQuality();
		void printQuality(double  improve_goal, int& nbad, double& minq, bool outWorst = false);
		void updateQuality(int i);
		void updateminVolume(void);
		///****************** Smoothing **************/
		int SmoothPass(int nloop, double improve_goal);
		// smooth

		int smooth_sus(int iNod);
		double quality_sus(double* a, double* b, double* c, double* d);

		// Volume Control
		int SmoothPassForVolume(int nloop, int smooth_type);
		int smooth_diff(int iNod, int type);
		int smooth_volume(int iNod, bool equalAngle = false);
		int getVolGrad(int iNod, std::vector<int> sph, double* VolGrad, std::vector<double> area);
		int getHessian(int iNod, std::vector<int> sph, double* Hessian, std::vector<double> area);
		void calGlobalEnergy(double& minV, double& maxV, int& nTet, double& Avg, double& Variance, double& Energy, int Energy_tpye);
		double getVolEnergy(std::vector<double> volume, std::vector<double> area);
		// Anisotropic
		int smooth_ani(int iNod);
		int smoothBndPnt_ani(int iNod);
		// Bnd
		int smoothBndPntPass(int loop);
		int smoothBndPnt(int iNod);
		// Disturbance
		int disturbPnt(const int iNod);
		int disturbPnt_search(const int iNod, double norm[]);
		
		///****************** Topological **************/
        struct TopologyCandidate { int tet; std::array<int, 4> form; };
        TopologyCandidate topologyCandidate(int t);
        bool topologyCandidateCurrent(const TopologyCandidate& candidate);
        bool needsTopologyInsertion(int t, double angleDegrees);
        bool hasBadDihedral(double angleDegrees);
        // Quality threshold for flips; minimum dihedral in degrees for insertion (0 disables it).
        int TopologicalPass(double improve_goal, double insert_angle_degrees, int nloop);
        int TopologicalPass_serial(const std::vector<TopologyCandidate>& candidates, double improve_goal, double insert_angle_degrees);
        int TopologicalPass_parallel(const std::vector<TopologyCandidate>& candidates, double improve_goal, double insert_angle_degrees, int workers);
        int tryTopologyInsertion(const TopologyCandidate& candidate, double angleDegrees);
        int removebadtet(int& iElm, int thread_n = -1);
        int removebadtet_addPnt(int iElm);
		int findtet(int p[], std::vector<int> sph);
		int matchtet(int p[], int t);

		///****************** Anisotropic adaptation **************/
		int lockingPass(std::vector<int> lockFactes, std::vector<int> lockVertex);
		void Smooth_size_ani(Mesh& mesh, std::vector<std::array<double, 6>>& anisol, std::vector<int> lockFactes, std::vector<int> lockVertex);
		
		// Anisotropy measure
		int getmm(int p1, int p2, int p3, int p4, double* mm);
		double caltet33_ani(double v1[3], double v2[3], double v3[3], double v4[3], double* AniMetric);
		double caltri33_ani(double v1[3], double v2[3], double v3[3], double* AniMetric);
		double cal_ani_length(double* v1, double* v2, double* m1, double* m2);
		int Interpolate_met33_ani(double* m1, double* m2, double* m, double s);
		int Interpolate_met(int m1, int m2, int m, double s);

		///****************** Size Control **************/
		int sizeControlPass(double lower, double upper);
		// Split 
		bool ifSplitBetter(int index);
		int splitBndEdgPass(double upper);
		int splitLongEdgPass(double upper);
		int splitLongEdgPass_noParallel(double upper);
		int splitEdg(int index, int deep = 0);
		// Contract
		int collapsePointLevel(int iNod);
		double contractionEdgeLength(int p1, int p2);
		int contractEdgPass(double lower);
		int contshortEdgPass(double lower);
		int canDestroyShortEdge(int iElm, int ia, int ib, double Threshold = 0);
		int destroyShortEdge(int iElm, int ia, int ib, double Threshold = 0);
		int tryDestroyShortEdge(int iElm, int ia, int ib, double Threshold = 0);
		int collapseEdg(int index);
		int collapseEdg(int index, int expectedDelete, int expectedKeep);
		int canCollapseBoundaryEdgeDirected(int edgeIndex, int deletePoint, int keepPoint);
		int canTryDestroyShortEdge(int iElm, int ia, int ib, double Threshold, int* destroyIa = nullptr, int* keepIb = nullptr);
		// filp
		int flipEdgPass(int nloop);
		void Type_Vertex_Edg(double Angle, const Mesh& mesh);;
		int ifflipEdg(int index);
		int checkFlipNormal(int p1, int p2, int p3, int p4);
		int flipEdg(int index, int deep = 0);

		///****************** OpenMP **************/
		//omp algorithm
		void evalNodesToSmooth(std::vector<int>& nodesToSmooth, double improve_goal);
        void colorBadQualityNodes(std::vector<std::vector<int>>& colors, double improve_goal);
	    // atomic operation
		inline bool tryOccupying(std::atomic<int>& occupying, int thread_n) {
			int current = occupying.load(std::memory_order_relaxed);

			if (current == thread_n) {
				return true;
			}

			if (current == -1) {
				int expected = -1;
				return occupying.compare_exchange_strong(expected, thread_n);
			}

			return false;
		}
		inline void clearOccupying(std::atomic<int>& occupying, int thread_n) {
			int expected = thread_n;
			occupying.compare_exchange_strong(expected, -1);
		}

		/**************************  function *************************/
		int addElem(bool UseVacancy = true);
		int addElem(int thread_n, bool UseVacancy = true);
		int addElem(int pa, int pb, int pc, int pd, bool UseVacancy = true);
		int addElem(int pa, int pb, int pc, int pd, int thread_n, bool UseVacancy = true);
		int addNode(bool UseVacancy = true);
		int addNode(double x, double y, double z, double space, bool UseVacancy = true);
		int addNode(double x, double y, double z, double space, int thread, bool UseVacancy = true);
		bool isDelEle(int E);
		bool isDelNod(int N);
		void DelEle(int E);
		void DelEle(int E, int thread_n);
		void DelNod(int N);
		//-------------------- Mesh info set and read --------------------
		int isNod_in_Tet(int iNod, int tet);
		int setP2T(int p, int t);
		int getP2T(int p);
		void setAllP2T();//set all point to tet
		bool isDelSurTri(int i);
		bool isDelSurEdg(int i);
		bool setDelSurTri(int i);
		bool setDelSurEdg(int i);
		bool isRecBndEdg(int i);
		bool isRecBndTri(int i);
		int findTriParent(int i);
		bool isBndEdg(const int p1, const int p2);
		bool isBndTri(const int p1, const int p2, const int p3);
		void clearbndpnt(int i);
		bool isbndpnt(int i);
		void setbndpnt(int i);
		bool isFacetpnt(int i);
		bool isSegmentpnt(int i);
		bool isCornerpnt(int i);
		inline int set_bit(int& info, int x) { return info = info | (1 << x); }
		inline int clear_bit(int& info, int x) { return info = info & ~(1 << x); }
		inline int get_bit(int info, int x) { return (info >> x) & 1; }
		inline int set_bit(int64_t& info, int x) { return info = info | (1 << x); }
		inline int clear_bit(int64_t& info, int x) { return info = info & ~(1 << x); }
		inline int get_bit(int64_t info, int x) { return (info >> x) & 1; }
		inline bool ishulltet(int tet) { return Elems[tet].form[3] == ghost; }
		inline bool isvirtualtet(int tet) { return Elems[tet].geo == virtualID; }
		inline bool setvirtualtet(int tet) { return Elems[tet].geo = virtualID; }
		void clearNodesElems();
		void getMeshEdgebyGeo(Mesh& mesh);
		//---------------- Adjacent connect and decode ----------------
		bool bond(int t1, int t2);
		void bond(int t1, int nig1, int t2, int nig2);
		inline int  getNeig(int t1, int i) {
			if (Elems[t1].neig[i] == -1) return -1;
			return Elems[t1].neig[i] >> 2;
		}
		inline int  getNeigOrd(int t1, int i) {
			if (Elems[t1].neig[i] == -1) return -1;
			return Elems[t1].neig[i] & 3;
		}
		int getoppoP(int tet, int i);
		//--------------------- Mesh Quality ---------------------
		void printfDihedral(double& minD, double& minAvgD, double& maxD, double& maxAvgD);
		int CalDihedral(double v1[3], double v2[3], double v3[3], double v4[3], double& minangle, double& maxangle, std::vector<double>& Dihedral);
		double tetquality(double v1[3], double v2[3], double v3[3], double v4[3], double* AniMetric, int qualmeasure);
		double vlrms3ratio(double v1[3], double v2[3], double v3[3], double v4[3]);
		double ScaledJacobian( double p0[3],  double p1[3],  double p2[3],  double p3[3], double eps= 1e-30);
		double AspectRatio(double v1[3], double v2[3], double v3[3], double v4[3]);
		//--------------------- Orthogonality ---------------------
		double calskewness(double* pa, double* pb, double* pc, double* pd);
		double orthogonal(double* pa, double* pb, double* pc, double* pd, double* oppositeA, double* oppositeB, double* oppositeC, double* oppositeD);

		//--------------------- Hilbert sort ---------------------
		void Hilbert(const std::vector<std::array<double, 3>>& V, std::vector<int>& order);
		void multiscale_sort(std::vector<std::array<double, 4>>& Varray, int, int& depth);
		void hilbert_sort(std::vector<std::array<double, 4>>& Varray, int bg, int Asize, int, int d, double, double, double, double, double, double, int depth);
		int  hilbert_split(std::vector<std::array<double, 4>>& Varray, int bg, int Asize, int gc0, int gc1, double, double, double, double, double, double);
		//----------------- Linear algebra operators ----------------
		bool calCircum(double* pa, double* pb, double* pc, double* pd, double* cent, double* radius);
		bool P_in_Line(double* a, double* b, double* p);
		bool isParallel(double* p1, double* p2, double* p3, double* p4);
		bool lu_decmp(double lu[4][4], int n, int* ps, double* d, int N);
		void lu_solve(double lu[4][4], int n, int* ps, double* b, int N);
		void cross(double* v1, double* v2, double* n);
		void calnormal(int b, int c, int d, double noraml[]);
		void tensorproduct33(double* v1, double* v2, double* ans);
		bool inverseM(double* in, double* out);
		void vecTimesMatrix13_33(double* vec, double* Matric, double* vecans);
		void  projectPointToPlane(double* v1, double* v2, double* v3, double* p, double* projection);
		void calFactor(int n, int& a, int& b, int& c);
		void calBarycenter(int i, double* pnt);
		void distanceToPlane(double p1[3], double p2[3], double p3[3], double p4[3], double& dis);
		void vsub3( double a[3],  double b[3], double r[3]);
		double segmentSegmentDistance(double* p1, double* q1, double* p2, double* q2);
		double calculateTriangleAngles(const double* A, const double* B, const double* C, std::vector<double>& angles);
		double calVolume(int i);
		double lenvec(double* x);
		double dot(double* v1, double* v2);
		double distance(double* p1, double* p2);
		double distance2(double* p1, double* p2);
		double calArea(double* p1, double* p2, double* p3);
		double norm2(double x, double y, double z);
		double insphere_s(int a, int b, int c, int d, int e);
		//------------------------ Special aim tool ------------------------
		//void randQueue(std::queue<int>& q);
		void addRandomP(double x, double y, double z, int n);
		///void export_ring_csv_min(int iNod, const double oldpos[3], int mode);
		//------------------------- Debug info ------------------------
		void spdlogoutfile(bool outlogfile);
		int checkEdgeLen();
		double getPeakMegabytesUsed();
		uint64_t getFreeMemory();
		void printMemoryUsage();
		void prinfPnt(int iNod);
		void checkMeshError();
		void checkEdgeDegree();
		void outHullTri(std::string filename);
		void outTempMesh(std::string filename);
		void outUnRecvEdge(std::string filename);
		void printSph(const std::vector<int>);
		void printSph_VTK(const std::vector<int>, std::string filename);
		double getTime(std::chrono::high_resolution_clock::time_point t1, std::chrono::high_resolution_clock::time_point t2);
		std::chrono::high_resolution_clock::time_point getTime_now();
		bool build_fine_mesh_projection_tree(const std::string& fine_mesh_file);
		bool project_boundary_point_to_fine_mesh(double* in, double* out, double max_projection_distance);
		bool project_segment_point_to_fine_mesh(double* in, double* out, double max_projection_distance);
    private:
        void resetMeshState();
        int createRefineCandidate(int tet, Args& args, bool& rejected);
        bool fileLogging = false;
        FineMeshProjection fineMeshProjection;
        // Only the serial BW compatibility entry uses this workspace. Parallel
        // algorithms own their BWPlans explicitly and never borrow this buffer.
        BWSerialWorkspace serialBW;
	};
}
#endif
