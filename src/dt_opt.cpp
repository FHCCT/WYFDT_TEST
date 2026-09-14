#include "dt.h"
#include <exception>
#include "spdlog/sinks/null_sink.h"

namespace {
constexpr double topologyShortEdgeRatio = 1.0 / 200.0;

// A tiny edge may join two feature vertices created by a narrow surface strip.
// Bound any geometric relaxation by both its own length and the local shell
// scale. Locked vertices, normal changes and positive volumes remain checked.
double shortBoundaryCollapseTolerance(DT& mesh, int pa, int pb, const std::vector<int>& shell) {
    if (!mesh.improve_step || !mesh.modifyBnd || !mesh.isbndpnt(pa) || !mesh.isbndpnt(pb)) return 0;
    double longestSquared = 0;
    for (int t : shell) {
        if (mesh.ishulltet(t) || mesh.isvirtualtet(t)) continue;
        for (int e = 0; e < 6; ++e)
            longestSquared = std::max(longestSquared, mesh.distance2(
                mesh.Nodes[mesh.Elems[t].form[Egid[e][0]]].pt,
                mesh.Nodes[mesh.Elems[t].form[Egid[e][1]]].pt));
    }
    const double edgeSquared = mesh.distance2(mesh.Nodes[pa].pt, mesh.Nodes[pb].pt);
    if (!(longestSquared > 0) || !std::isfinite(longestSquared) ||
        !(edgeSquared < longestSquared * topologyShortEdgeRatio * topologyShortEdgeRatio)) return 0;
    return std::min(std::sqrt(edgeSquared), 1e-3 * std::sqrt(longestSquared));
}
}

//-------------------------mesh improvement-------------------------
/* optlevel info
*  0:Don't improvement
*  Reduction of maximum dihedral
*  1: volume-Side ratio
*  2: volume-Side ratio + Angle
*  3: volume-Side ratio + Angle + Jacobian
*  4:
*  5:
*  7:quickly opt
*  10:Vbase
*  11:Abase
*  12:Lap
*  13:CPT
*  14:ODT
*  15:UVAT
*/
int DT::MeshImprove(Args& args) {
    infolevel = std::max(0, std::min(2, args.infolevel));
    meshLogger->set_level(infolevel == 0 ? spdlog::level::err :
        infolevel == 1 ? spdlog::level::info : spdlog::level::debug);
	if (args.optlevel == 0)
		return 0;
	if (infolevel > 0)
		meshLogger->debug("Mesh improve start");

    MeshStageLog stageLog(*this, "Optimization", 1, MeshStageSummary::MeshChange);
	if (!improve_init(args))
		return 0;

	if (args.optlevel == 1 || args.optlevel == 2 || args.optlevel == 3) {
		//Traditional boundary constrain optimization
		TraditionalOptPass(args);
	}
	else if (args.optlevel == 5) {
		//quickly opt,modify bnd
		OrthogonalityOptPass(args);
	}
	else if (args.optlevel == 7) {
		//quickly opt,modify bnd
		QuicklyOptPass(args);
	}
	else if (args.optlevel == 8) {
		//quickly flip, Coarse
		QuicklyOptCoarsePass(args);
	}
	else if (args.optlevel == 9) {
		//sizecontrol by anisotropic size
		OptbyAnisotropicPass(args);
	}
	else if (args.optlevel >= 10 && args.optlevel <= 15) {
		//for volume&angle uniform paper
		return VolumeImprovePass(args);
	}
	else {
		meshLogger->warn("This optimization method is awaiting development... ...");
	}

	return 0;
}

int DT::improve_init(Args& args) {
	//"Initialize Quality Field,3 times should be sufficient.
	improve_step = 1;//tell dt it is improving now;

	for (int i = 0; i < EdgSteiner.size(); i++)
		setbndpnt(EdgSteiner[i].first);
	for (int i = 0; i < TriSteiner.size(); i++)
		setbndpnt(TriSteiner[i].first);

	if (!threadsInitialized) initializeDTThreads(*this, args);
	flipnmRecll.reserve(100000);

	//pre Thread pools
	int EvacancySize = 100;
	const int poolThreads = activeDTThreads(*this);
	for (int i = 0; i < poolThreads; i++) {
		while (Evacancy_thread[i].size() < EvacancySize) {
			int loc = addElem();
			DelEle(loc, i);
		}
	}

	//printf opt info
	if (infolevel > 0) {
		if (args.optlevel == 8) meshLogger->debug("Opt_level         : Coarsening");
		else meshLogger->debug("Opt_level          : {}", args.optlevel);
		meshLogger->debug("Opt_Loop_num       : {}", args.optloop);
		meshLogger->debug("Opt_Threshold : {}", args.optTh);
		meshLogger->debug("Opt_Angle_strict: {}", args.optanglestrict);
	}

	return 1;
}

int DT::TraditionalOptPass(Args& args) {
	int nLoop = 3, nbad = 0;
	double improve_goal = 1, minD = 0, minAvgD = 0, maxD = 0, maxAvgD = 0, minq;
	nLoop = std::max(args.optloop, nLoop);
	improve_goal = args.optTh;

	improve_Metric = SUS_METRIC;
	prepareQuality();
	printQuality(improve_goal, nbad, minq, false);
	minVolume_bw = 1e-16;

	for (int loop = 0; loop < nLoop; loop++) {
		if (infolevel > 0)
			meshLogger->info("Loop {}:", loop);

		if (args.constrain == 0) {
			updateminVolume();
			improve_step = false;
			flipEdgPass(1);
			improve_step = true;
		}
        else if (infolevel > 0) meshLogger->info("Flip bndEdg: {}/{}", 0, 0);

		prepareQuality();
		TopologicalPass(improve_goal, args.optlevel > 2 ? args.optanglestrict : 0.0, 1);
		if (args.optlevel > 1) SmoothPass(1, improve_goal);
		printQuality(improve_goal, nbad, minq, false);

		if (nbad == 0) break;
	}

	printfDihedral(minD, minAvgD, maxD, maxAvgD);
	if (args.outworsttet)
		printQuality(improve_goal, nbad, minq, true);

	return 1;
}

int DT::OptbyAnisotropicPass(Args& args) {
	int nLoop = 10, nbad = 0, splitcontract_old = INT_MAX;
	double improve_goal = 1, minD = 0, minAvgD = 0, maxD = 0, maxAvgD = 0, minq;
	optflipdeep = 0;
	improve_Metric = 5;//Anisotropic quality check
	improve_goal = 0.0288675;
	nLoop = args.optloop;// std::max(args.optloop, 3);
	minVolume_bw = 1e-16;

	for (int loop = 0; loop < nLoop; loop++) {
		if (infolevel > 0) meshLogger->info("Loop {}:", loop);

		int splitcontract = sizeControlPass(0.65, 1.3);

		improve_Metric = 5;//Anisotropic quality check
		improve_goal = 0.0288675;
		prepareQuality();

		TopologicalPass(improve_goal, 0.0, 1);

		SmoothPass(1, 1);

		improve_Metric = 3;//Volume
		improve_goal = 1e-12;
		prepareQuality();
		for (int loop0 = 0; loop0 < 2; loop0++) {
			TopologicalPass(improve_goal, args.optanglestrict, 1);
			SmoothPass(1, improve_goal);
			printQuality(improve_goal, nbad, minq, false);
			if (nbad == 0) break;
		}
	}

	improve_Metric = 3;//Volume
	improve_goal = 1e-10;
	prepareQuality();
	for (int loop = 0; loop < 2; loop++) {
		TopologicalPass(improve_goal, args.optanglestrict, 1);
		SmoothPass(1, improve_goal);
		printQuality(improve_goal, nbad, minq, false);
		if (nbad == 0) break;
	}

	improve_Metric = 2;//
	improve_goal = std::sin(ANGLE2RADIO(1e-2));
	prepareQuality();
	for (int loop = 0; loop < 2; loop++) { 
		TopologicalPass(improve_goal, args.optanglestrict, 1);
		SmoothPass(1, improve_goal);
		printQuality(improve_goal, nbad, minq, false);
		if (nbad == 0) break;
	}

	return 1;
}

int DT::OrthogonalityOptPass(Args& args) {
	int nLoop = 3, nbad = 0;
	double  minD = 0, minAvgD = 0, maxD = 0, maxAvgD = 0, minq;

	double improve_goal = 0.1;

	nLoop = std::max(args.optloop, nLoop);
	improve_Metric = 8; //Orthogonality
	prepareQuality();
	printQuality(0.01, nbad, minq, false);

	for (int loop = 0; loop < nLoop; loop++) {
		TopologicalPass(improve_goal, 0.0, 1);
		SmoothPass(1, improve_goal);
		printQuality(0.01, nbad, minq, false);
	}
	return 1;
}

int DT::QuicklyOptPass(Args& args) {
	int nLoop = 1, nbad = 0;
	double  minD = 0, minAvgD = 0, maxD = 0, maxAvgD = 0, minq;
	double improve_goal = 1;// // (2.0 * std::sin(ANGLE2RADIO(args.optangle)) * args.optratio) / (std::sin(ANGLE2RADIO(args.optangle)) + args.optratio);
	optflipdeep = 0;

	nLoop = std::max(args.optloop, nLoop);
	for (int loop = 0; loop < nLoop ; loop++) {
		updateminVolume();

		//Get min volume,yunbo size control don't use minvolume
		improve_step = false;
		flipEdgPass(1);
		improve_step = true;

		improve_Metric = SUS_METRIC;
		improve_goal = args.optTh;
		prepareQuality();
		//printQuality(improve_goal, nbad, minq, false);
		TopologicalPass(improve_goal, args.optanglestrict, 1);
		SmoothPass(1, improve_goal);
		printQuality(improve_goal, nbad, minq, false);

		if (!hasBadDihedral(args.optanglestrict))
			break;
	}
	return 1;
}

int DT::OptSizeControl_yunbo(double lamasize, 
	std::unordered_map<int, double>& facetSize,
	std::unordered_map<int, int>& facetNum,
	std::unordered_map<int, double>& elementSize,
	std::unordered_map<int, int>& elementNum, Args& args) {
	int nLoop = 3, nbad = 0, nstop = 0;
	double minD = 0, minAvgD = 0, maxD = 0, maxAvgD = 0, minq;
	double improve_goal = 1;// = args.optTh;// (2.0 * std::sin(ANGLE2RADIO(args.optangle)) * args.optratio) / (std::sin(ANGLE2RADIO(args.optangle)) + args.optratio);
	optflipdeep = 0;

	QuantityControl = true;//for BW/flip geoNum Quantity Control

	nLoop = std::max(args.optloop, nLoop);
	for (int loop = 0; loop < nLoop; loop++) {
		//Get min volume,yunbo size control don't use minvolume

		updateminVolume();

		flipEdgPass(1);

		improve_Metric = SUS_METRIC;
		improve_goal = args.optTh;
		prepareQuality();
		//printQuality(improve_goal, nbad, minq, false);
		TopologicalPass(improve_goal, args.optanglestrict, 1);
		SmoothPass(1, improve_goal);
		printQuality(improve_goal, nbad, minq, false);

		// flip boundary Edges,it is not correct, need to do
		improve_step = false;
		GeoNum.clear();
		for (int i = 0; i < Elems.size(); i++) {
			if (isDelEle(i) || ishulltet(i) || isvirtualtet(i))
				continue;
			GeoNum[Elems[i].geo]++;
		}

		FacetNum.clear();
		for (int i = 0; i < SurTris.size(); i++) {
			if (isDelSurTri(i))
				continue;
			int fid = SurTris[i].parent;
			FacetNum[fid]++;
		}

		struct splitSort {
			int EdgeID;
			int geoID;
			double len;
		};

		std::vector<splitSort> waitsplit;

		for (int i = 0; i < Elems.size(); i++) {
			if (isDelEle(i) || ishulltet(i) || isvirtualtet(i))
				continue;

			bool ifsplit = false;
			int geoID = Elems[i].geo;

			double maxEdgLen = 0;
			int maxEdgIdx = -1;
			for (int j = 0; j < 6; j++) {
				int p1 = Elems[i].form[Egid[j][0]];
				int p2 = Elems[i].form[Egid[j][1]];

				if (!isbndpnt(p1) || !isbndpnt(p2) || !isBndEdg(p1, p2))
					continue;

				////Old Edg idndex may be destroy
				double dis = distance(Nodes[p1].pt, Nodes[p2].pt);
				if (dis > maxEdgLen) {
					maxEdgLen = dis;
					maxEdgIdx = j;
				}
			}

			if (lamasize > 0 && maxEdgLen > lamasize * 4.0 / 3.0)
				ifsplit = true;

			if (!ifsplit) {
				if (elementSize[geoID] > 0 && maxEdgLen > elementSize[geoID]) {
					ifsplit = true;
				}
				else if (elementNum[geoID] > 0 && GeoNum[geoID] < elementNum[geoID]) {
					ifsplit = true;
				}
			}

			if (ifsplit && maxEdgIdx != -1) {
				int p1 = Elems[i].form[Egid[maxEdgIdx][0]];
				int p2 = Elems[i].form[Egid[maxEdgIdx][1]];
				if (auto* boundaryEntry = BndEdg.find(p1, p2)) {
					const int boundaryIndex = *boundaryEntry;
					int Edgid = boundaryIndex;
					if (isDelSurEdg(Edgid))
						continue;
					waitsplit.push_back({ Edgid,geoID ,maxEdgLen });
				}
			}
		}
		
		for (int i = 0; i < SurEdgs.size() ; i++) {
			if (isDelSurEdg(i) || lockE.count(i))
				continue;
			// not find
			double dis = distance(Nodes[SurEdgs[i].iStart].pt, Nodes[SurEdgs[i].iEnd].pt);
			for (int j = 0; j < SurEdgs[i].face.size(); j++) {
				auto it = std::find_if(waitsplit.begin(), waitsplit.end(),
					[i](const splitSort& s) { return s.EdgeID == i; });
				if (it != waitsplit.end()) {
					break;
				}
				int fid = SurTris[SurEdgs[i].face[j]].parent;
				if (facetSize.count(fid) && facetSize[fid] > 0) {
					if (dis > facetSize[fid]) {
						if (facetNum.count(fid) && facetNum[fid] > 0 && FacetNum[fid] >= facetNum[fid]) {
							continue;
						}
						waitsplit.push_back({ i,-fid,dis });
					}
				}
				else if (facetNum.count(fid) && facetNum[fid] > 0 && FacetNum[fid] < facetNum[fid]) {
					waitsplit.push_back({ i,-fid,dis });
				}
			}
		}

		// 按 len 递增排序
		std::sort(waitsplit.begin(), waitsplit.end(),
			[](const splitSort& a, const splitSort& b) {
				return a.len > b.len;
			});

		int stop = 0;

		for (auto it : waitsplit) {
			if (it.geoID>0 && elementNum[it.geoID]>0 && GeoNum[it.geoID] > elementNum[it.geoID])
				continue;
			if (it.geoID < 0) {
				int fid = -it.geoID;
				if (facetNum[fid]>0 && FacetNum[fid] > facetNum[fid])
					continue;
			}
			splitEdg(it.EdgeID);
			stop++;
		}

		GeoNum.clear();
		for (int i = 0; i < Elems.size(); i++) {
			if (isDelEle(i) || ishulltet(i) || isvirtualtet(i))
				continue;
			GeoNum[Elems[i].geo]++;
		}

		for (int i = 0; i < Elems.size(); i++) {
			if (isDelEle(i) || ishulltet(i) || isvirtualtet(i))
				continue;

			bool ifsplit = false;
			int geoID = Elems[i].geo;

			double maxEdgLen = 0;
			int maxEdgIdx = -1;
			for (int j = 0; j < 6; j++) {
				int p1 = Elems[i].form[Egid[j][0]];
				int p2 = Elems[i].form[Egid[j][1]];

				if (isBndEdg(p1, p2))
					continue;
				////Old Edg idndex may be destroy
				double dis = distance(Nodes[p1].pt, Nodes[p2].pt);
				if (dis > maxEdgLen) {
					maxEdgLen = dis;
					maxEdgIdx = j;
				}
			}

			if (lamasize > 0 && maxEdgLen > lamasize * 4.0 / 3.0)
				ifsplit = true;

			if (!ifsplit) {
				if (elementSize[geoID] > 0 && maxEdgLen > elementSize[geoID]) {
					ifsplit = true;
				}
				else if (elementNum[geoID] > 0 && GeoNum[geoID] < elementNum[geoID]) {
					ifsplit = true;
				}
			}

			if (ifsplit) {
				int p1 = Elems[i].form[Egid[maxEdgIdx][0]];
				int p2 = Elems[i].form[Egid[maxEdgIdx][1]];

				double addV[3] = { 0.5 * (Nodes[p1].pt[0] + Nodes[p2].pt[0]), 0.5 * (Nodes[p1].pt[1] + Nodes[p2].pt[1]), 0.5 * (Nodes[p1].pt[2] + Nodes[p2].pt[2]) };

				int iNod = addNode(addV[0], addV[1], addV[2], Nodes[p1].space * 0.5 + Nodes[p2].space * 0.5);

				int searchtet = i;
				locate_pnt(iNod, searchtet);
				if (isvirtualtet(searchtet)) {
					DelNod(iNod);
					continue;
				}

				//add inner point
				std::vector<int> nearTets = { searchtet };

				int ret = BW_insert_vertex(iNod, nearTets, 1);
				if (ret == 1) {//success
					stop++;
					continue;
				}
				else {//failed
					DelNod(iNod);
				}
			}
		}

		if (stop < 10 && minq>improve_goal)
			break;
		improve_step = true;
	}

	return 1;
}

int DT::QuicklyOptCoarsePass(Args& args) {
	int nLoop = 3, nbad = 0;
	double improve_goal = 1, minD = 0, minAvgD = 0, maxD = 0, maxAvgD = 0, minq;
	improve_Metric = 1;
	improve_goal = SUS_METRIC;
	nLoop = std::max(args.optloop, 1);
	prepareQuality();
	for (int loop = 0; loop < nLoop; loop++) {
		printQuality(improve_goal, nbad, minq, false);
		//sizeControlPass(1e10, 1e10);
		improve_step = false;

		int contract1 = contractEdgPass(1e10);
		if (infolevel > 0) meshLogger->debug("Contract  BndEdg: {}", contract1);
		int contract2 = contshortEdgPass(1e10);
		if (infolevel > 0) meshLogger->debug("Contract MeshEdg: {}", contract2);

		flipEdgPass(1);

		improve_step = true;

		prepareQuality();
		TopologicalPass(improve_goal, 0.0, 1);
		SmoothPass(1, improve_goal);

		if (contract1 + contract2 < 100)
			break;
	}

	printQuality(improve_goal, nbad, minq, false);
	return 1;
}

int DT::VolumeImprovePass(Args& args) {
	//Volume uniformity
	meshLogger->debug("Volume uniformity");
	int nLoop = std::max(args.optloop, 100), ntet, fail = 0, nbad;
	double minV, maxV, Avg, Energy, Variance, oldminAvgD, max_minAvgD = 0;
	double minq = 0, minD = 0, minAvgD = 0, maxD = 0, maxAvgD = 0, improve_goal;
    // The stopping criterion needs angles even when diagnostic scans are disabled.
    const auto updateDihedral = [&]() {
        if (infolevel > 0) printfDihedral(minD, minAvgD, maxD, maxAvgD);
        else calculateDihedral(minD, minAvgD, maxD, maxAvgD);
    };

	int smooth_type = args.optlevel < 15 ? args.optlevel - 10 : 0;
	int Energy_tpye = 2;
	if (smooth_type == 0) Energy_tpye = 0;
	if (smooth_type == 1) Energy_tpye = 1;
	improve_Metric = 3;
	prepareQuality();
	meshLogger->debug("Initial:");
	if (infolevel >= 2) calGlobalEnergy(minV, maxV, ntet, Avg, Variance, Energy, Energy_tpye);
	updateDihedral();
	oldminAvgD = minAvgD;
	max_minAvgD = std::max(max_minAvgD, std::floor(minAvgD * (smooth_type == 0 ? 10.0 : 100.0)) / (smooth_type == 0 ? 10.0 : 100.0));
	for (int loop = 0; loop < nLoop; loop++) {
		meshLogger->info("Loop {}:", loop);
		//sizeControlPass(minEdge, maxEdge);
		improve_Metric = 1;
		improve_goal = 0.4;
		//improve_goal = std::sin(ANGLE2RADIO(40));
		prepareQuality();
		TopologicalPass(improve_goal, 0.0, 10);
		//printQuality(improve_goal, nbad, minq, false);
		//flipEdgPass(3);
		improve_Metric = 3;
		prepareQuality();
		SmoothPassForVolume(5, smooth_type);
		updateDihedral();
		if (infolevel >= 2) calGlobalEnergy(minV, maxV, ntet, Avg, Variance, Energy, Energy_tpye);
		//outTempMesh("./temp_" + std::to_string(loop) + ".vtk");
		if (std::floor(minAvgD * (smooth_type == 0 ? 10.0 : 100.0)) / (smooth_type == 0 ? 10.0 : 100.0) <= max_minAvgD) {
			fail++;
		}
		else {
			max_minAvgD = std::floor(minAvgD * (smooth_type == 0 ? 10.0 : 100.0)) / (smooth_type == 0 ? 10.0 : 100.0);
			fail = 0;
		}
		if (fail == 10)
			break;
	}
	if (args.optlevel == 15) {
		fail = 0;
		for (int loop = 0; loop < nLoop; loop++) {
			meshLogger->info("Loop {}", loop);
			//sizeControlPass(minEdge, maxEdge);
			improve_Metric = 1;
			//improve_goal = std::sin(ANGLE2RADIO(40));
			improve_goal = 0.4;
			prepareQuality();
			TopologicalPass(improve_goal, 0.0, 10);
			//flipEdgPass(3);
			improve_Metric = 3;
			prepareQuality();
			SmoothPassForVolume(5, 1);
			updateDihedral();
			if (infolevel >= 2) calGlobalEnergy(minV, maxV, ntet, Avg, Variance, Energy, 1);
			if (std::floor(minAvgD * 100.0) / 100.0 <= max_minAvgD) {
				fail++;
			}
			else {
				max_minAvgD = std::floor(minAvgD * 100.0) / 100.0;
				fail = 0;
			}
			if (fail == 10)
				break;
		}
	}
	return 0;
}

int DT::prepareQuality() {
    MeshStageLog stageLog(*this, "Evaluate quality", 2);
    DTParallelScope parallelScope;
    const int workers = activeDTThreads(*this, 0, Elems.size());
#pragma omp parallel for num_threads(workers) schedule(static) if(workers > 1)
	for (int i = 0; i < Elems.size(); i++)
	{
		if (isDelEle(i) || isvirtualtet(i)  || (ghost != -1 && ishulltet(i)))
			continue;
		int a = Elems[i].form[0];
		int b = Elems[i].form[1];
		int c = Elems[i].form[2];
		int d = Elems[i].form[3];

		double q = 0;

		if (improve_Metric != 8) {
			double AniMetric[6] = { 0 };
			if (improve_Metric == 5) {
				getmm(a, b, c, d, AniMetric);
			}
			q = tetquality(Nodes[a].pt, Nodes[b].pt, Nodes[c].pt, Nodes[d].pt, AniMetric, improve_Metric);
			if (q < 0) {
				double ori = calVolume(i);
				meshLogger->warn("Warning:Inverted tet: {} {}", i, ori);
				q = tetquality(Nodes[a].pt, Nodes[b].pt, Nodes[c].pt, Nodes[d].pt, AniMetric, improve_Metric);
			}
		}
		else {
			q = orthogonal(Nodes[a].pt, Nodes[b].pt, Nodes[c].pt, Nodes[d].pt,
				getoppoP(i, 0) == ghost ? NULL : Nodes[getoppoP(i, 0)].pt,
				getoppoP(i, 1) == ghost ? NULL : Nodes[getoppoP(i, 1)].pt,
				getoppoP(i, 2) == ghost ? NULL : Nodes[getoppoP(i, 2)].pt,
				getoppoP(i, 3) == ghost ? NULL : Nodes[getoppoP(i, 3)].pt);
			if (q < 0) {
				double ori = calVolume(i);
				q = orthogonal(Nodes[a].pt, Nodes[b].pt, Nodes[c].pt, Nodes[d].pt,
					getoppoP(i, 0) == ghost ? NULL : Nodes[getoppoP(i, 0)].pt,
					getoppoP(i, 1) == ghost ? NULL : Nodes[getoppoP(i, 1)].pt,
					getoppoP(i, 2) == ghost ? NULL : Nodes[getoppoP(i, 2)].pt,
					getoppoP(i, 3) == ghost ? NULL : Nodes[getoppoP(i, 3)].pt);
			}
		}
		Elems[i].q = q;
	}
	return 1;
}

void DT::printQuality(double  improve_goal, int& nbad, double& minq, bool outWorst) {
	// Quality drives stopping conditions at every verbosity level.
	int nElem = Elems.size(), nSum = 0;
	double SumQ = 0;
	minq = DBL_MAX;
	nbad = 0;

	for (int i = 0; i < nElem; i++) {
		if (isDelEle(i) || isvirtualtet(i) || ishulltet(i))
			continue;

		double q = Elems[i].q;

		SumQ += Elems[i].q;
		nSum++;
		if (Elems[i].q < improve_goal)
			nbad++;

		if (Elems[i].q < minq) {
			minq = Elems[i].q;
		}
	}

	if (infolevel > 0) {
		meshLogger->info("Quality : bad:{} minQ:{:3e} avg:{:.3e}", nbad, minq, nSum ? SumQ / nSum : 0.0);
	}

	return;
}

void DT::updateQuality(int i) {
	if (ishulltet(i) || isvirtualtet(i) || (ghost != -1 && ishulltet(i)))
		return;
	int a = Elems[i].form[0];
	int	b = Elems[i].form[1];
	int	c = Elems[i].form[2];
	int	d = Elems[i].form[3];
	double q = 0;
	if (improve_Metric != 8) {
		double AniMetric[6] = { 0 };
		if (improve_Metric == 5) {
			getmm(a, b, c, d, AniMetric);
		}
		q = tetquality(Nodes[a].pt, Nodes[b].pt, Nodes[c].pt, Nodes[d].pt, AniMetric, improve_Metric);
		//if (q < 0) {
		//	double ori = calVolume(i);
		//	meshLogger->warn("Warning:Inverted tet: {} {}", i, ori);
		//	q = tetquality(Nodes[a].pt, Nodes[b].pt, Nodes[c].pt, Nodes[d].pt, AniMetric, improve_Metric);
		//}
	}
	else {
		q = orthogonal(Nodes[a].pt, Nodes[b].pt, Nodes[c].pt, Nodes[d].pt,
			getoppoP(i, 0) == ghost ? NULL : Nodes[getoppoP(i, 0)].pt,
			getoppoP(i, 1) == ghost ? NULL : Nodes[getoppoP(i, 1)].pt,
			getoppoP(i, 2) == ghost ? NULL : Nodes[getoppoP(i, 2)].pt,
			getoppoP(i, 3) == ghost ? NULL : Nodes[getoppoP(i, 3)].pt);
	}
	Elems[i].q = q;
	return;
}

void DT::updateminVolume(void) {
    DTParallelScope parallelScope;
    const int workers = activeDTThreads(*this, 0, Elems.size());
    std::vector<double> minima(workers, 1e-10);
#pragma omp parallel num_threads(workers) if(workers > 1)
    {
        double localMin = 1e-10;
#pragma omp for schedule(static)
        for (int t = 0; t < static_cast<int>(Elems.size()); ++t) {
            if (isDelEle(t) || ishulltet(t) || isvirtualtet(t)) continue;
            localMin = std::min(calVolume(t), localMin);
        }
        minima[omp_get_thread_num()] = localMin;
    }
    minVolume_bw = *std::min_element(minima.begin(), minima.end());
}

//Smooth loop_parallel
//if checkQ open ,keep every smooth is improvement
int DT::SmoothPass(int nloop, double improve_goal)
{
    MeshStageLog stageLog(*this, "Smooth", 2);
    DTParallelScope parallelScope;
    if (susReuseIdle) susIdleStates.resize(Nodes.size());
    std::vector<std::vector<int>> colors;
    colorBadQualityNodes(colors, improve_goal);
    size_t maxGroup = 0;
    for (const auto& group : colors) maxGroup = std::max(maxGroup, group.size());
    const int workers = activeDTThreads(*this, 0, maxGroup);
    long long success = 0, attempts = 0;
    // One team per pass; each color ends with a barrier before the next color.
#pragma omp parallel num_threads(workers) if(workers > 1) reduction(+:success, attempts)
    {
        for (int loop = 0; loop < nloop; ++loop) {
            for (size_t color = 0; color < colors.size(); ++color) {
                const auto& group = colors[color];
#pragma omp for schedule(dynamic, 16)
                for (int i = 0; i < static_cast<int>(group.size()); ++i)
                {
                    ++attempts;
                    success += smooth_sus(group[i]) == 1;
                }
            }
        }
    }
    if (infolevel > 0) meshLogger->info("Smooth: {}/{}", success, attempts);
    return 1;
}

//conbine laplace smooth and numerical smooth
int DT::smooth_ani(int iNod) {
	if (AniSol.size() == 0 || isbndpnt(iNod))
		return 0;
	double minq = DBL_MAX, sum_vol = 0, mm[6], vol, det, wcen[3] = { 0 }, d[3], newpos[3], alpha = 1;
	double* verts[4];
	std::vector<int> sph;
	std::vector<double> q;

	findSphere(iNod, sph);
	int n = sph.size();
	//get worse quality
	q.resize(n);

	//smooth by Anisotropic weighted centroid
	for (int i = 0; i < n; i++)
	{
		int ie = sph[i];
		int* form = Elems[ie].form;

		minq = std::min(minq, Elems[sph[i]].q);

		// Compute average metric (6 components)
		for (int j = 0; j < 6; j++) {
			mm[j] = 0.25 * (AniSol[form[0]][j] + AniSol[form[1]][j] + AniSol[form[2]][j] + AniSol[form[3]][j]);
		}

		vol = calVolume(ie);

		// Determinant of metric matrix
		det = mm[0] * (mm[3] * mm[5] - mm[4] * mm[4])
			- mm[1] * (mm[1] * mm[5] - mm[2] * mm[4])
			+ mm[2] * (mm[1] * mm[4] - mm[2] * mm[3]);

		if (det <= 1e-30)
			return 0.0;

		vol *= sqrt(det);

		sum_vol += vol;

		//get weighted centroid
		for (int m = 0; m < 3; m++)
			wcen[m] += 0.25 * vol * (Nodes[form[0]].pt[m] + Nodes[form[1]].pt[m] + Nodes[form[2]].pt[m] + Nodes[form[3]].pt[m]);
	}

	if (sum_vol <= 1e-30)
		return 0.0;

	for (int m = 0; m < 3; m++)
	{
		wcen[m] /= sum_vol;					//weighted
		d[m] = wcen[m] - Nodes[iNod].pt[m];				//displacement
		newpos[m] = Nodes[iNod].pt[m] + alpha * d[m];	//new position
	}

	bool moveflag = true;
	int iter = 0;
	while (iter < 6) {
		moveflag = true;
		for (int i = 0; i < n; i++)
		{
			for (int m = 0; m <= 3; m++)
			{
				int iElemNd = Elems[sph[i]].form[m];
				verts[m] = (iElemNd != iNod) ? Nodes[iElemNd].pt : newpos;
			}

			double AniMetric[6] = { 0 };
			getmm(Elems[sph[i]].form[0], Elems[sph[i]].form[1], Elems[sph[i]].form[2], Elems[sph[i]].form[3], AniMetric);

			q[i] = tetquality(verts[0], verts[1], verts[2], verts[3], AniMetric, improve_Metric);
			if (q[i] < minq) {
				moveflag = false;
				break; // This tet becomes invalid.
			}
		}
		if (moveflag) {
			break;
		}
		else {
			alpha /= 2.0;
			for (int j = 0; j < 3; j++) {
				newpos[j] = Nodes[iNod].pt[j] + alpha * d[j];
			}
			iter++;
		}
	} // while (iter < 3)
	if (moveflag) {
		//update qual
		for (int i = 0; i < n; i++)
			Elems[sph[i]].q = q[i];

		//update point position
		for (int j = 0; j < 3; j++)
			Nodes[iNod].pt[j] = newpos[j];

		return 1;
	}
	return 0;
}


DT::TopologyCandidate DT::topologyCandidate(int t) {
    TopologyCandidate candidate; candidate.tet = t;
    std::copy(Elems[t].form, Elems[t].form + 4, candidate.form.begin());
    return candidate;
}

bool DT::topologyCandidateCurrent(const TopologyCandidate& candidate) {
    const int t = candidate.tet;
    return t >= 0 && t < static_cast<int>(Elems.size()) && !isDelEle(t) &&
        !isvirtualtet(t) && !ishulltet(t) &&
        std::equal(candidate.form.begin(), candidate.form.end(), Elems[t].form);
}

bool DT::needsTopologyInsertion(int t, double angleDegrees) {
    if (!(angleDegrees > 0) || t < 0 || t >= static_cast<int>(Elems.size()) ||
        isDelEle(t) || isvirtualtet(t) || ishulltet(t)) return false;
    double minAngle = 0, maxAngle = 0;
    std::vector<double> angles;
    const auto& p = Elems[t].form;
    // CalDihedral expects the opposite orientation to DT's stored ordering.
    if (!CalDihedral(Nodes[p[0]].pt, Nodes[p[1]].pt, Nodes[p[3]].pt, Nodes[p[2]].pt,
        minAngle, maxAngle, angles)) return true;
    return !std::isfinite(minAngle) || RADIO2ANGLE(minAngle) < angleDegrees;
}

bool DT::hasBadDihedral(double angleDegrees) {
    if (!(angleDegrees > 0)) return false;
    for (int t = 0; t < static_cast<int>(Elems.size()); ++t)
        if (needsTopologyInsertion(t, angleDegrees)) return true;
    return false;
}

int DT::tryTopologyRepair(const TopologyCandidate& candidate, double angleDegrees) {
    if (parallelTopologyBatch || omp_in_parallel())
        throw std::logic_error("Topology repair must run outside a parallel region");
    if (!topologyCandidateCurrent(candidate)) return 0;

    // Collapse is independent of the angle threshold, including when insertion
    // is disabled. Lengths here are physical lengths, not sizing-field lengths.
    if (improve_step && modifyBnd) {
        std::array<double, 6> squaredLengths;
        double longestSquared = 0;
        for (int e = 0; e < 6; ++e) {
            squaredLengths[e] = distance2(Nodes[candidate.form[Egid[e][0]]].pt,
                Nodes[candidate.form[Egid[e][1]]].pt);
            longestSquared = std::max(longestSquared, squaredLengths[e]);
        }
        if (longestSquared > 0 && std::isfinite(longestSquared)) {
            const double limit = longestSquared * topologyShortEdgeRatio * topologyShortEdgeRatio;
            for (int e = 0; e < 6; ++e) {
                if (!(squaredLengths[e] < limit)) continue;
                const int p1 = candidate.form[Egid[e][0]], p2 = candidate.form[Egid[e][1]];
                const int* entry = BndEdg.find(p1, p2);
                if (!entry || isDelSurEdg(*entry) || SurEdgs[*entry].info > 1 || lockE.count(*entry)) continue;
                // Preserve distinct feature curves, as in contractEdgPass.
                if (SurEdgs[*entry].constrain == 0 && collapsePointLevel(p1) >= 2 && collapsePointLevel(p2) >= 2)
                    continue;
                // The master entry selects face -> segment -> corner, tries
                // both directions at equal levels, and handles periodic pairs.
                if (collapseEdg(*entry) == 1) return 1;
                if (!topologyCandidateCurrent(candidate)) return 0;
            }
        }
    }
    if (!needsTopologyInsertion(candidate.tet, angleDegrees)) return 0;
    return removebadtet_addPnt(candidate.tet);
}

int DT::TopologicalPass(double improve_goal, double insert_angle_degrees, int nloop) {
    MeshStageLog stageLog(*this, "Topology", 2);
    if (!std::isfinite(insert_angle_degrees) || insert_angle_degrees < 0 || insert_angle_degrees > 180)
        throw std::invalid_argument("Topology insertion angle must be in [0, 180] degrees");
    size_t totalSuccess = 0, totalCandidates = 0;
    for (int loop = 0; loop < nloop; ++loop) {
        std::vector<TopologyCandidate> candidates;
        for (int t = 0; t < static_cast<int>(Elems.size()); ++t) {
            if (isDelEle(t) || isvirtualtet(t) || ishulltet(t) || Elems[t].q < 0) continue;
            // Only the quality threshold selects candidates. Check the angle
            // of an individual candidate only after its flip attempt fails.
            if (Elems[t].q <= improve_goal)
                candidates.push_back(topologyCandidate(t));
        }
        const int workers = activeDTThreads(*this, 0, candidates.size());
        int success;
        if (workers > 1 && improve_step)
            success = TopologicalPass_parallel(candidates, improve_goal, insert_angle_degrees, workers);
        else
            success = TopologicalPass_serial(candidates, improve_goal, insert_angle_degrees);
        totalSuccess += success;
        totalCandidates += candidates.size();
        if (success < 10) break;
    }
    if (infolevel > 0) meshLogger->info("Topology: {}/{}", totalSuccess, totalCandidates);
    return 1;
}

int DT::TopologicalPass_serial(const std::vector<TopologyCandidate>& candidates,
    double improve_goal, double insert_angle_degrees) {
    MeshStageLog stageLog(*this, "Serial topology", 2);
    int success = 0;
    for (const auto& candidate : candidates) {
        if (!topologyCandidateCurrent(candidate) || Elems[candidate.tet].q < 0 || Elems[candidate.tet].q > improve_goal) continue;
        int t = candidate.tet;
        const int result = removebadtet(t, -1);
        if (result == 1) ++success;
        else if (result == 0 && t >= 0 && !isDelEle(t))
            success += tryTopologyRepair(topologyCandidate(t), insert_angle_degrees) == 1;
    }
    return success;
}

namespace {
// Scheduler-only scratch storage. Dense membership avoids per-element hash
// allocation on large passes; small passes avoid allocating mesh-sized arrays.
class TopologyIndexSet {
    bool dense_;
    std::vector<unsigned char> marked_;
    std::unordered_set<int> sparse_;
    std::vector<int> values_;
public:
    explicit TopologyIndexSet(size_t capacity) : dense_(capacity != 0), marked_(capacity, 0) {
        values_.reserve(128);
    }
    void clear() {
        if (dense_) for (int value : values_) marked_[value] = 0;
        else sparse_.clear();
        values_.clear();
    }
    bool count(int value) const {
        return dense_ ? value >= 0 && static_cast<size_t>(value) < marked_.size() && marked_[value] != 0
                      : sparse_.count(value) != 0;
    }
    bool insert(int value) {
        if (dense_) {
            if (static_cast<size_t>(value) >= marked_.size()) marked_.resize(static_cast<size_t>(value) + 1, 0);
            if (marked_[value]) return false;
            marked_[value] = 1;
        } else if (!sparse_.insert(value).second) return false;
        values_.push_back(value);
        return true;
    }
    template<class Iterator> void insert(Iterator first, Iterator last) {
        for (; first != last; ++first) insert(*first);
    }
    size_t size() const { return values_.size(); }
    std::vector<int>::const_iterator begin() const { return values_.begin(); }
    std::vector<int>::const_iterator end() const { return values_.end(); }
};
}

int DT::TopologicalPass_parallel(const std::vector<TopologyCandidate>& candidates,
    double improve_goal, double insert_angle_degrees, int workers) {
    MeshStageLog stageLog(*this, "Parallel topology", 2);
    DTParallelScope parallelScope;
    workers = std::min(workers, activeDTThreads(*this, 0, candidates.size()));
    if (workers <= 1) return TopologicalPass_serial(candidates, improve_goal, insert_angle_degrees);
    std::queue<TopologyCandidate> pending;
    for (const auto& candidate : candidates) pending.push(candidate);
    int success = 0;
    const bool dense = candidates.size() >= 256;
    const size_t tetCapacity = dense ? Elems.size() : 0, nodeCapacity = dense ? Nodes.size() : 0;
    TopologyIndexSet batchReadTets(tetCapacity), batchWriteTets(tetCapacity), batchWriteNodes(nodeCapacity);
    TopologyIndexSet seedTets(tetCapacity), readTets(tetCapacity), writeTets(tetCapacity), shellVertices(nodeCapacity);
    TopologyIndexSet visited(tetCapacity);
    std::vector<int> star, points;
    std::vector<unsigned char> starOrdinals;
    std::vector<int> batch, result;
    std::vector<std::exception_ptr> errors;
    std::vector<TopologyCandidate> failed;
    batch.reserve(workers); result.reserve(workers); errors.reserve(workers); failed.reserve(workers);
    // All callers have already checked the tetrahedron index. Keep its four
    // vertex reads together in this hot, read-only traversal.
    auto pointOrdinal = [&](int node, int tet) {
        const auto& cell = Elems[tet];
        if (cell.info < 0) return -1;
        for (int j = 0; j < 4; ++j) if (cell.form[j] == node) return j;
        return -1;
    };
    // Same traversal as findSphere, using reusable LOCAL visit marks. This
    // lambda runs only during serial planning, never alongside mesh mutation.
    auto findPlanningStar = [&](int node) {
        star.clear(); starOrdinals.clear(); visited.clear();
        if (node < 0 || node >= static_cast<int>(Nodes.size())) return;
        int seed = getP2T(node);
        int seedOrdinal = seed >= 0 && seed < static_cast<int>(Elems.size())
            ? pointOrdinal(node, seed) : -1;
        if (seedOrdinal < 0) {
            seed = -1;
            for (int t = 0; t < static_cast<int>(Elems.size()); ++t)
                if (!isDelEle(t) && isNod_in_Tet(node, t) >= 0) { seed = t; break; }
        }
        if (seed < 0) return;
        if (seedOrdinal < 0) seedOrdinal = pointOrdinal(node, seed);
        visited.insert(seed); star.push_back(seed); starOrdinals.push_back(static_cast<unsigned char>(seedOrdinal));
        for (size_t k = 0; k < star.size(); ++k) {
            const int t = star[k], ord = starOrdinals[k];
            for (int f = 0; f < 4; ++f) {
                if (f == ord) continue;
                const int next = getNeig(t, f);
                // Planning never mutates the mesh: an already visited cell
                // remains live and contains this node until traversal finishes.
                if (next < 0 || next >= static_cast<int>(Elems.size()) || visited.count(next)) continue;
                const int nextOrdinal = pointOrdinal(node, next);
                if (nextOrdinal < 0) continue;
                visited.insert(next); star.push_back(next);
                starOrdinals.push_back(static_cast<unsigned char>(nextOrdinal));
            }
        }
    };
    while (!pending.empty()) {
        batch.clear();
        batchReadTets.clear(); batchWriteTets.clear(); batchWriteNodes.clear();
        size_t requiredSlots = 0;
        const size_t lookahead = std::min(pending.size(), static_cast<size_t>(workers) * 8);
        for (size_t attempt = 0; attempt < lookahead && batch.size() < workers; ++attempt) {
            TopologyCandidate c = pending.front(); pending.pop();
            if (!topologyCandidateCurrent(c) || Elems[c.tet].q < 0 || Elems[c.tet].q > improve_goal) continue;
            // Depth zero changes only seed-edge shells. Rollback searches the
            // full stars of their vertices; findSphere also reads face neighbors.
            shellVertices.clear(); seedTets.clear(); readTets.clear(); writeTets.clear();
            seedTets.insert(c.tet);
            for (int edge = 0; edge < 6; ++edge) {
                const int a = Egid[edge][0], b = Egid[edge][1];
                if (isBndEdg(c.form[a], c.form[b])) continue;
                findShell(c.tet, a, b, star, points);
                seedTets.insert(star.begin(), star.end());
            }
            for (int f = 0; f < 4; ++f) {
                int face[3], k = 0;
                for (int j = 0; j < 4; ++j) if (j != f) face[k++] = c.form[j];
                const int n = getNeig(c.tet, f);
                if (!isBndTri(face[0], face[1], face[2]) && n >= 0 && !isDelEle(n)) seedTets.insert(n);
            }
            writeTets.insert(seedTets.begin(), seedTets.end());
            for (int t : seedTets) {
                for (int v : Elems[t].form) shellVertices.insert(v);
                for (int f = 0; f < 4; ++f) {
                    const int n = getNeig(t, f);
                    if (n >= 0) writeTets.insert(n); // bond() writes the outside neighbor.
                }
            }
            // Reject an overlapping write neighborhood before constructing the
            // more expensive read neighborhood. The acceptance test is unchanged.
            bool conflict = false;
            for (int v : shellVertices) if (batchWriteNodes.count(v)) { conflict = true; break; }
            if (!conflict) for (int t : writeTets) if (batchReadTets.count(t)) { conflict = true; break; }
            if (conflict) { pending.push(c); continue; }
            for (int p : shellVertices) {
                findPlanningStar(p);
                for (int t : star) {
                    readTets.insert(t);
                    // findSphere tests the neighbor before testing membership.
                    for (int f = 0; f < 4; ++f) {
                        const int n = getNeig(t, f);
                        if (n >= 0) readTets.insert(n);
                    }
                }
            }
            readTets.insert(writeTets.begin(), writeTets.end());
            for (int t : readTets) if (batchWriteTets.count(t)) { conflict = true; break; }
            if (conflict) { pending.push(c); continue; }
            batchWriteNodes.insert(shellVertices.begin(), shellVertices.end());
            batchReadTets.insert(readTets.begin(), readTets.end());
            batchWriteTets.insert(writeTets.begin(), writeTets.end());
            batch.push_back(c.tet);
            // At depth zero each forward 2->3 reduces the remaining shell by
            // one. Rollback returns slots; this exceeds the peak live growth.
            requiredSlots = std::max(requiredSlots, seedTets.size() * 12 + 32);
        }
        if (batch.empty()) continue;
        // Allocate outside OpenMP: even SmallVector's block table and size
        // are shared mutable state. Each worker only reuses its own free slots.
        for (int w = 0; w < workers; ++w)
            while (Evacancy_thread[w].size() < requiredSlots) {
                int t = addElem(); DelEle(t, w);
            }
        result.assign(batch.size(), 0);
        errors.assign(batch.size(), std::exception_ptr{});
        parallelTopologyBatch = true;
#pragma omp parallel for num_threads(workers) schedule(static, 1)
        for (int j = 0; j < static_cast<int>(batch.size()); ++j) {
            try { result[j] = removebadtet(batch[j], omp_get_thread_num()); }
            catch (...) { errors[j] = std::current_exception(); }
        }
        parallelTopologyBatch = false;
        for (const auto& error : errors) if (error) std::rethrow_exception(error);
        // Finish this batch's failed attempts before starting another batch.
        // Snapshot identities first: a collapse or insertion can invalidate another.
        failed.clear();
        for (size_t j = 0; j < batch.size(); ++j) {
            if (result[j] == 1) ++success;
            else if (result[j] == 0 && batch[j] >= 0 && !isDelEle(batch[j]))
                failed.push_back(topologyCandidate(batch[j]));
        }
        for (const auto& candidate : failed)
            success += tryTopologyRepair(candidate, insert_angle_degrees) == 1;
    }
    return success;
}

// Pure flips: parallel calls use depth zero; serial calls retain optflipdeep.
int DT::removebadtet(int& iElm, int thread_n) {
	if (iElm < 0 || iElm >= static_cast<int>(Elems.size())) return 0;
	if (isDelEle(iElm))
		return 1;//remove success

	const int deepth = thread_n >= 0 ? 0 : (improve_step ? optflipdeep : 10);

	int p[4] = { Elems[iElm].form[0] ,Elems[iElm].form[1] ,Elems[iElm].form[2] ,Elems[iElm].form[3] };
	if (thread_n != -1) {
		for (int j = 0; j < 4; j++) {
			if (!tryOccupying(Nodes[p[j]].occupying, thread_n)) {
				for (int k = 0; k < j; k++)
					clearOccupying(Nodes[p[k]].occupying, thread_n);
				return -1;
			}
		}
	}

    struct SeedRelease {
        DT* mesh; const int* points; int owner;
        ~SeedRelease() {
            if (owner >= 0) for (int j = 0; j < 4; ++j)
                mesh->clearOccupying(mesh->Nodes[points[j]].occupying, owner);
        }
    } seedRelease{this, p, thread_n};

	for (int i = 0; i < 6; i++) {
		//find a edge not try
		int pa = Elems[iElm].form[Egid[i][0]];
		int pb = Elems[iElm].form[Egid[i][1]];
		if (isBndEdg(pa, pb)) {
			continue;
		}
		std::vector<int> wait_remove = { iElm };

		int ret = removeEdge(wait_remove, Egid[i][0], Egid[i][1], deepth, thread_n);

		if (ret == 1) {
			return 1;
		}

		//remove fail but tet still idx change
		iElm = findtet(p, wait_remove);
		if (iElm == -1) {
			return 0;
		}
	}

	for (int i = 0; i < 4; i++) {
		std::vector<int> wait_remove = { iElm };

		int ret = removeface(wait_remove, i, deepth, thread_n);
		if (ret == 1) {
			return 1;//remove success
		}
		//remove fail but tet still idx change
		iElm = findtet(p, wait_remove);
		if (iElm == -1) {
			return 0;
		}
	}


	return 0;
}

int DT::removebadtet_addPnt(int iElm) {
    if (parallelTopologyBatch || omp_in_parallel())
        throw std::logic_error("Topology insertion must run outside a parallel region");
    if (iElm < 0 || iElm >= static_cast<int>(Elems.size()) || isDelEle(iElm) ||
        isvirtualtet(iElm) || ishulltet(iElm)) return 0;
	int p0 = Elems[iElm].form[0];
	int p1 = Elems[iElm].form[1];
	int p2 = Elems[iElm].form[2];
	int p3 = Elems[iElm].form[3];

	double pnt[3] = { 0.0 }, space = 0;
	int a, b, c, d;
	std::vector<int> shell, shellp;

	int nBndpnt = isbndpnt(p0) + isbndpnt(p1) + isbndpnt(p2) + isbndpnt(p3);
	int nBndedg = 0;
	for (int j = 0; j < 6; j++) {
		int pa = Elems[iElm].form[Egid[j][0]];
		int pb = Elems[iElm].form[Egid[j][1]];

		if (isBndEdg(pa, pb)) {
			nBndedg++;
		}
	}

	//printf("%d %d\n", nBndpnt, nBndedg);
	//std::vector<int> worst = { iElm };
	//printSph_VTK(worst, "./" + std::to_string(iElm) + "_" + std::to_string(nBndpnt) + "_" + std::to_string(nBndedg) + "_" + ".vtk");

	if (nBndpnt >= 4) {
		{
			int maxIdx = -1;
			double maxDis = 0;
			/******** Split Long Edge ********/
			std::set<int> alltet;
			for (int i = 0; i < 4; i++) {
				std::vector<int> sph;
				findSphere(Elems[iElm].form[i], sph);
				for (auto it : sph) {
					if (ishulltet(it) || isvirtualtet(it))
						continue;
					alltet.insert(it);
				}
			}

			int tp1 = -1, tp2 = -1;
			for (auto it : alltet) {
				for (int j = 0; j < 6; j++) {
					int pa = Elems[it].form[Egid[j][0]];
					int pb = Elems[it].form[Egid[j][1]];

					double dis = distance(Nodes[pa].pt, Nodes[pb].pt);

					if (auto* boundaryEntry = BndEdg.find(pa, pb)) {
						const int boundaryIndex = *boundaryEntry;
						if (!modifyBnd)
							continue;

						int Edgid = boundaryIndex;
						if (SurEdgs[Edgid].constrain > 0)
							dis *= 1.12;//1.118 1.25
					}

					if (dis > maxDis) {
						maxDis = dis;
						tp1 = pa;
						tp2 = pb;
					}
				}
			}

			if (tp1 != -1 && tp2 != -1) {
				if (auto* boundaryEntry = BndEdg.find(tp1, tp2)) {
					const int boundaryIndex = *boundaryEntry;
					int Edgid = boundaryIndex;
					int ret = splitEdg(Edgid);
				}
				else {
					int tempt = iElm;
					isMeshEdge(tp1, tp2, &tempt);
					//findShell(tempt, isNod_in_Tet(tp1, tempt), isNod_in_Tet(tp2, tempt), shell, shellp);

					space = (Nodes[tp1].space + Nodes[tp2].space) * 0.5;
					for (int k = 0; k < 3; k++)
						pnt[k] = (Nodes[tp1].pt[k] + Nodes[tp2].pt[k]) * 0.5;

					int newp = addNode(pnt[0], pnt[1], pnt[2], space);

					tempt = iElm;
					int loc = locate_pnt(newp, tempt);
					if (loc < 1 || isvirtualtet(tempt) || ishulltet(tempt)) {//iNod is same point,and it's idx = -sameP
						DelNod(newp);
					}
					else {
						if (AniSol.size() > 0)
							Interpolate_met(tp1, tp2, newp, 0.5);

						shell.clear();
						shell.push_back(tempt);
						//B_W
						int ret = BW_insert_vertex(newp, shell, 1);
						if (ret == 1) {
							smooth_sus(newp);
						}
						else {
							DelNod(newp);
						}
					}
				}
			}
		}

		if (!isDelEle(iElm) && Elems[iElm].form[0] == p0 && Elems[iElm].form[1] == p1 && Elems[iElm].form[2] == p2 && Elems[iElm].form[3] == p3) {
			double len[6];
			int order[6] = { 0,1,2,3,4,5 };

			for (int j = 0; j < 6; j++) {
				int pa = Elems[iElm].form[Egid[j][0]];
				int pb = Elems[iElm].form[Egid[j][1]];
				len[j] = distance(Nodes[pa].pt, Nodes[pb].pt);
				if (auto* boundaryEntry = BndEdg.find(pa, pb)) {
					const int boundaryIndex = *boundaryEntry;
					int Edgid = boundaryIndex;
					if (SurEdgs[Edgid].constrain > 0)
						len[j] *= 1.12;
				}
			}

			// sort index
			std::sort(order, order + 6, [&](int a, int b) {
				return len[a] > len[b];
				});

			for (int idx = 0; idx < 6; idx++) {
				if (isDelEle(iElm) || Elems[iElm].form[0] != p0 || Elems[iElm].form[1] != p1 || Elems[iElm].form[2] != p2 || Elems[iElm].form[3] != p3)
					break;

				int j = order[idx];
				int pa = Elems[iElm].form[Egid[j][0]];
				int pb = Elems[iElm].form[Egid[j][1]];

				if (auto* boundaryEntry = BndEdg.find(pa,pb)) {
					const int boundaryIndex = *boundaryEntry;
					if (modifyBnd) {
						int Edgid = boundaryIndex;
						int ret = splitEdg(Edgid);
					}
					continue;
				}

				for (int k = 0; k < 3; k++)
					pnt[k] = Nodes[pa].pt[k] + Nodes[pb].pt[k];

				space = Nodes[pa].space + Nodes[pb].space;

				if (nBndpnt == 4 && nBndedg == 5) {

					findShell(iElm, Egid[j][0], Egid[j][1], shell, shellp);

					for (int t = 0; t < shellp.size(); t++) {
						for (int k = 0; k < 3; k++)
							pnt[k] += Nodes[shellp[t]].pt[k];
						space += Nodes[shellp[t]].space;
					}

					space /= 1.0 * (shellp.size() + 2);
					for (int k = 0; k < 3; k++)
						pnt[k] /= 1.0 * (shellp.size() + 2);
				}
				else {
					space /= 2.0;
					for (int k = 0; k < 3; k++)
						pnt[k] /= 2.0;
				}

				int newp = addNode(pnt[0], pnt[1], pnt[2], space);

				int tempt = iElm;
				int loc = locate_pnt(newp, tempt);
				if (loc < 1 || isvirtualtet(tempt) || ishulltet(tempt)) {//iNod is same point,and it's idx = -sameP
					DelNod(newp);
				}
				else {
					shell.clear();
					shell.push_back(tempt);
					//B_W
					int ret = BW_insert_vertex(newp, shell, 1);
					if (ret == 1) {
						smooth_sus(newp);
					}
					else {
						DelNod(newp);
					}
				}
			}
		}
	
		if (!isDelEle(iElm) && Elems[iElm].form[0] == p0 && Elems[iElm].form[1] == p1 && Elems[iElm].form[2] == p2 && Elems[iElm].form[3] == p3) {
			for (int k = 0; k < 3; k++)
				pnt[k] = (Nodes[p0].pt[k] + Nodes[p1].pt[k] + Nodes[p2].pt[k] + Nodes[p3].pt[k]) / 4.0;

			space = (Nodes[p0].space + Nodes[p1].space + Nodes[p2].space + Nodes[p3].space) / 4.0;

			int newp = addNode(pnt[0], pnt[1], pnt[2], space);

			int tempt = iElm;
			int loc = locate_pnt(newp, tempt);
			if (loc < 1 || isvirtualtet(tempt) || ishulltet(tempt) || tempt != iElm) {
				DelNod(newp);
			}
			else {
				shell.clear();
				shell.push_back(tempt);

				//B_W
				int ret = BW_insert_vertex(newp, shell, 1);
				if (ret == 1) {
					smooth_sus(newp);
				}
				else {
					DelNod(newp);
				}
			}
		}
	}
    else {
        // The current insertion strategy handles cells with four boundary vertices.
        return 0;
    }

	if (isDelEle(iElm) || Elems[iElm].form[0] != p0 || Elems[iElm].form[1] != p1 || Elems[iElm].form[2] != p2 || Elems[iElm].form[3] != p3)
		return 1;

	return 0;
}

//find tet have p0,p1,p2,p3
int DT::findtet(int p[], std::vector<int> sph) {
	int i;
	for (i = 0; i < sph.size(); i++)
		if (isNod_in_Tet(p[0], sph[i]) != -1)
			if (isNod_in_Tet(p[1], sph[i]) != -1)
				if (isNod_in_Tet(p[2], sph[i]) != -1)
					if (isNod_in_Tet(p[3], sph[i]) != -1)
						return matchtet(p, sph[i]);
	findSphere(p[0], sph);
	for (i = 0; i < sph.size(); i++)
		if (isNod_in_Tet(p[0], sph[i]) != -1)
			if (isNod_in_Tet(p[1], sph[i]) != -1)
				if (isNod_in_Tet(p[2], sph[i]) != -1)
					if (isNod_in_Tet(p[3], sph[i]) != -1)
						return matchtet(p, sph[i]);
	return -1;
}

//Adjust the order of the 4 points of the tet
int DT::matchtet(int p[], int t) {
	int oldp[4] = { 0 }, neig[4] = { 0 }, neigOrd[4] = { 0 };
	for (int i = 0; i < 4; i++) {
		if (Elems[t].form[i] == p[0]) { oldp[0] = i; }
		else if (Elems[t].form[i] == p[1]) { oldp[1] = i; }
		else if (Elems[t].form[i] == p[2]) { oldp[2] = i; }
		else if (Elems[t].form[i] == p[3]) { oldp[3] = i; }
		neig[i] = getNeig(t, i);
		neigOrd[i] = getNeigOrd(t, i);
	}
	for (int i = 0; i < 4; i++) {
		Elems[t].form[i] = p[i];
		bond(t, i, neig[oldp[i]], neigOrd[oldp[i]]);
	}
	return t;
}

int DT::SmoothPassForVolume(int nloop, int smooth_type) {
    MeshStageLog stageLog(*this, "Volume smooth", 2);
	std::vector<int> nodes;
	evalNodesToSmooth(nodes, 1);
    long long success = 0, attempts = 0;

	for (int loop = 0; loop < nloop; loop++) {
		//if constrain,there will change
		int iSuccess = 0;
		int iTried = 0;
		//Volume uniformity

		//Don't use openmp
		for (int i = 0; i < nodes.size(); ++i)
		{
			const int iNode = nodes[i];
			if (isbndpnt(iNode) || isDelNod(iNode) || iNode == ghost)
				continue;
			iTried++;
			if (smooth_type == 0) {//volume equal
				if (smooth_volume(iNode, false) == 1)
					iSuccess++;
			}
			else if (smooth_type == 1) {//angle equal
				if (smooth_volume(iNode, true) == 1)
					iSuccess++;
			}
			else if (smooth_type == 2) {//laplace
				if (smooth_diff(iNode, 0) == 1)
					iSuccess++;
			}
			else if (smooth_type == 3) {//weighted centre of gravity
				if (smooth_diff(iNode, 1) == 1)
					iSuccess++;
			}
			else if (smooth_type == 4) {//NODT,weighted external center
				if (smooth_diff(iNode, 2) == 1)
					iSuccess++;
			}
		}

        success += iSuccess;
        attempts += iTried;
	}
    if (infolevel > 0) meshLogger->info("Smooth: {}/{}", success, attempts);
	return 0;
}
/*
* Different smoothing methods
* type==0: laplace
* type==1: weighted centre of gravity
* type==2: NODT,weighted external center
*/
//Use with caution when having negative volumes
int DT::smooth_diff(int iNod, int type) {
	int i, m, k, iElem, iSecond;
	double newpos[3] = { 0 }, aimpos[3] = { 0 }, d[3] = { 0 }, alpha = 1;
	double* verts[4];
	std::vector<int> sph;
	std::unordered_map<int, bool> seen;
	std::vector<double> newVolume;

	findSphere(iNod, sph);
	int n = sph.size();
	newVolume.resize(n);
	for (i = 0; i < n; i++) {
		if (ishulltet(sph[i]))
			return 0;
	}

	if (type == 0) {
		//smooth by laplace
		for (i = 0; i < n; i++) {
			iElem = sph[i];
			for (m = 0; m < 4; m++) {
				iSecond = Elems[iElem].form[m];
				if (iSecond != iNod) {
					if (seen.find(iSecond) == seen.end()) {
						seen[iSecond] = true;
						aimpos[0] += Nodes[iSecond].pt[0];
						aimpos[1] += Nodes[iSecond].pt[1];
						aimpos[2] += Nodes[iSecond].pt[2];
					}
				}
			}
		}
		for (int m = 0; m < 3; m++) {
			aimpos[m] /= seen.size();
		}
	}
	else if (type == 1) {
		//CPT
		double volSum = 0;
		for (i = 0; i < n; i++) {
			iElem = sph[i];
			double vol = calVolume(iElem);
			volSum += vol;
			for (m = 0; m < 4; m++) {
				iSecond = Elems[iElem].form[m];
				aimpos[0] += Nodes[iSecond].pt[0] * vol;
				aimpos[1] += Nodes[iSecond].pt[1] * vol;
				aimpos[2] += Nodes[iSecond].pt[2] * vol;
			}
		}
		for (int m = 0; m < 3; m++) {
			aimpos[m] /= 4.0 * volSum;
		}
	}
	else if (type == 2) {
		//ODT
		double volSum = 0;
		for (i = 0; i < n; i++) {
			iElem = sph[i];
			double vol = calVolume(iElem);
			double cen[3] = { 0 }, rad = 0;
			volSum += vol;
			if (!calCircum(
				Nodes[Elems[iElem].form[0]].pt,
				Nodes[Elems[iElem].form[1]].pt,
				Nodes[Elems[iElem].form[2]].pt,
				Nodes[Elems[iElem].form[3]].pt, cen, &rad)) {
				//coplanar,can't calculate outer sphere
				return 0;
			}
			aimpos[0] += cen[0] * vol;
			aimpos[1] += cen[1] * vol;
			aimpos[2] += cen[2] * vol;
		}
		for (int m = 0; m < 3; m++) {
			aimpos[m] /= volSum;
		}
	}
	else if (type == 3) {
		//A_base
		for (i = 0; i < n; i++) {
			iElem = sph[i];
			double p[3] = { 0 };
			int ia = -1, ib = -1, ic = -1, id = -1;
			DFC(isNod_in_Tet(iNod, iElem), ia, ib, ic, id);
			projectPointToPlane(
				Nodes[Elems[iElem].form[ib]].pt,
				Nodes[Elems[iElem].form[ic]].pt,
				Nodes[Elems[iElem].form[id]].pt,
				Nodes[Elems[iElem].form[ia]].pt, p);
			aimpos[0] += p[0];
			aimpos[1] += p[1];
			aimpos[2] += p[2];
		}
		for (int m = 0; m < 3; m++) {
			aimpos[m] /= n;
		}
	}

	for (int m = 0; m < 3; m++) {
		d[m] = aimpos[m] - Nodes[iNod].pt[m];				//displacement
		newpos[m] = Nodes[iNod].pt[m] + alpha * d[m];		//new position
	}

	bool moveflag = true;
	int iter = 0;
	while (iter < 6) {
		moveflag = true;
		for (i = 0; i < n; i++) {
			for (int m = 0; m <= 3; m++)
			{
				int iElemNd = Elems[sph[i]].form[m];
				verts[m] = (iElemNd != iNod) ? Nodes[iElemNd].pt : newpos;
			}
			newVolume[i] = tetquality(verts[0], verts[1], verts[2], verts[3], NULL, improve_Metric);
			if (newVolume[i] <= 0) {
				moveflag = false;
				break; // This tet becomes invalid.
			}
		}
		if (moveflag) {
			break;
		}
		else {
			alpha /= 2.0;
			for (int j = 0; j < 3; j++) {
				newpos[j] = Nodes[iNod].pt[j] + alpha * d[j];
			}
			iter++;
		}
	} // while (iter < 6)

	if (moveflag) {
		for (i = 0; i < n; i++) {
			Elems[sph[i]].q = newVolume[i];//use volume
		}

		//update point position
		for (int j = 0; j < 3; j++)
			Nodes[iNod].pt[j] = newpos[j];

		if (type == 3 && lenvec(d) > 1e-6) {
			smooth_diff(iNod, 3);
		}
		return 1;
	}
	return 0;
}

//Smooth for Volume uniformity
int DT::smooth_volume(int iNod, bool equalAngle) {
	int i, m, k, iElem, iSecond;
	double newpos[3], oldEnergy, newEnergy = DBL_MAX, minq = 1e10;
	double* verts[4];
	std::vector<int> sph;
	std::vector<double> area, oldVolume, newVolume;
	double oldpos[3] = { Nodes[iNod].pt[0],Nodes[iNod].pt[1],Nodes[iNod].pt[2] };

	findSphere(iNod, sph);
	int n = sph.size();

	if (!improve_step) {
		for (i = 0; i < n; i++) {
			if (ishulltet(sph[i])) {
				return 0;
			}

			Elems[sph[i]].q = calVolume(sph[i]);
		}
	}

	//get worse quality
	if (!equalAngle) {
		area.resize(0);
	}
	else {
		area.resize(n);
	}
	newVolume.resize(n);
	oldVolume.resize(n);

	for (i = 0; i < n; i++) {
		iElem = sph[i];
		oldVolume[i] = Elems[iElem].q;

		if (equalAngle) {
			for (m = 0; m < 4; m++) {
				iSecond = Elems[iElem].form[m];
				if (iSecond == iNod) {
					area[i] = calArea(Nodes[Elems[iElem].form[(m + 1) % 4]].pt,
						Nodes[Elems[iElem].form[(m + 2) % 4]].pt, Nodes[Elems[iElem].form[(m + 3) % 4]].pt);
				}
			}
		}
	}
	oldEnergy = getVolEnergy(oldVolume, area);

	int descentNum = 0;
	while (1) {
		descentNum++;//Number of descents
		double VolGrad_initial[3] = { 0 };
		double VolGrad[3] = { 0 };
		double Hessian[9] = { 0 };
		double HessianT[9] = { 0 };
		getVolGrad(iNod, sph, VolGrad_initial, area);
		getHessian(iNod, sph, Hessian, area);
		if (!inverseM(Hessian, HessianT)) {
			return 0;
		}
		vecTimesMatrix13_33(VolGrad_initial, HessianT, VolGrad);

		double alpha = 1;
		double beta = 0.8;
		double gamma = 1e-4;//0.01;
		while (1) {
			//new position
			for (i = 0; i < 3; i++)
				newpos[i] = Nodes[iNod].pt[i] - VolGrad[i] * alpha;
			for (i = 0; i < n; i++) {
				for (int m = 0; m <= 3; m++)
				{
					int iElemNd = Elems[sph[i]].form[m];
					verts[m] = (iElemNd != iNod) ? Nodes[iElemNd].pt : newpos;
				}
				newVolume[i] = tetquality(verts[0], verts[1], verts[2], verts[3], NULL, improve_Metric);
				if (newVolume[i] < 0)
					break;
			}

			if (i == n) {
				//all tet positive now
				newEnergy = getVolEnergy(newVolume, area);
				if (newEnergy <= oldEnergy + gamma * alpha * dot(VolGrad, VolGrad_initial)) {
					break;
				}
			}

			alpha *= beta;

			if (alpha < /*1e-8*/1e-10)
				break;
		}
		if (newEnergy < oldEnergy) {
			//can't improvement more
			for (i = 0; i < n; i++) {
				Elems[sph[i]].q = newVolume[i];
			}

			for (int i = 0; i < 3; i++)
				Nodes[iNod].pt[i] = newpos[i];
			//if (newEnergy / oldEnergy > 0.9999) {
			if (fabs((newEnergy - oldEnergy) / oldEnergy) < 1e-5){
				return 1;
			}
			else {
				oldEnergy = newEnergy;
			}
		}
		else {

			if (descentNum > 1) {
				//export_ring_csv_min(iNod, oldpos, 1);
				return 1;
			}
			else
				return 0;
		}
	}
	return 0;
}

int DT::getVolGrad(int iNod, std::vector<int> sph, double* VolGrad, std::vector<double> area) {
	VolGrad[0] = VolGrad[1] = VolGrad[2] = 0;
	for (int i = 0; i < sph.size(); i++) {
		int j = 0, a, b, c, d;
		for (; j < 4; j++) {
			if (Elems[sph[i]].form[j] == iNod)
				break;
		}
		DFC(j, d, a, b, c);
		double* pa = Nodes[Elems[sph[i]].form[a]].pt;
		double* pb = Nodes[Elems[sph[i]].form[b]].pt;
		double* pc = Nodes[Elems[sph[i]].form[c]].pt;
		double papb[3] = { pb[0] - pa[0],  pb[1] - pa[1],  pb[2] - pa[2] };
		double papc[3] = { pc[0] - pa[0],  pc[1] - pa[1],  pc[2] - pa[2] };
		double ans[3] = { 0 };
		cross(papb, papc, ans);
		if (area.size() == 0) {
			VolGrad[0] += Elems[sph[i]].q * ans[0] / 3.0;
			VolGrad[1] += Elems[sph[i]].q * ans[1] / 3.0;
			VolGrad[2] += Elems[sph[i]].q * ans[2] / 3.0;
		}
		else {
			VolGrad[0] += Elems[sph[i]].q * ans[0] / area[i] / 3.0;
			VolGrad[1] += Elems[sph[i]].q * ans[1] / area[i] / 3.0;
			VolGrad[2] += Elems[sph[i]].q * ans[2] / area[i] / 3.0;
		}
	}
	return 0;
}

int DT::getHessian(int iNod, std::vector<int> sph, double* Hessian, std::vector<double> area) {
	Hessian[0] = Hessian[1] = Hessian[2] = Hessian[3] = Hessian[4] = Hessian[5] = Hessian[6] = Hessian[7] = Hessian[8] = 0;
	for (int i = 0; i < sph.size(); i++) {
		int j = 0, a, b, c, d;
		for (; j < 4; j++) {
			if (Elems[sph[i]].form[j] == iNod)
				break;
		}
		DFC(j, d, a, b, c);
		double* pa = Nodes[Elems[sph[i]].form[a]].pt;
		double* pb = Nodes[Elems[sph[i]].form[b]].pt;
		double* pc = Nodes[Elems[sph[i]].form[c]].pt;
		double papb[3] = { pb[0] - pa[0],  pb[1] - pa[1],  pb[2] - pa[2] };
		double papc[3] = { pc[0] - pa[0],  pc[1] - pa[1],  pc[2] - pa[2] };
		double ans[3] = { 0 };
		double hes[9] = { 0 };
		cross(papb, papc, ans);
		tensorproduct33(ans, ans, hes);
		if (area.size() == 0) {
			for (j = 0; j < 9; j++) {
				Hessian[j] += hes[j] / 18.0;
			}
		}
		else {
			for (j = 0; j < 9; j++) {
				Hessian[j] += hes[j] / area[i] / 18.0;
			}
		}
	}

	for (int r = 0; r < 3; ++r) {
		for (int c = r + 1; c < 3; ++c) {
			double sym = 0.5 * (Hessian[r * 3 + c] + Hessian[c * 3 + r]);
			Hessian[r * 3 + c] = Hessian[c * 3 + r] = sym;
		}
	}
	return 0;
}

int  DT::disturbPnt(const int iNod) {
	if (isbndpnt(iNod))
		return 0;
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<double> distr(0.0, 1.0);
	std::vector<int> sph;
	findSphere(iNod, sph);
	double pold[3] = { Nodes[iNod].pt[0],Nodes[iNod].pt[1],Nodes[iNod].pt[2] };
	for (int i = 0; i < 10; i++) {
		Nodes[iNod].pt[0] = pold[0] + 1e-6 * distr(eng);
		Nodes[iNod].pt[1] = pold[1] + 1e-6 * distr(eng);
		Nodes[iNod].pt[2] = pold[2] + 1e-6 * distr(eng);
		int j = 0;
		for (j = 0; j < sph.size(); j++) {
			if (ishulltet(sph[j]))
				continue;
			double vol = calVolume(sph[j]);
			if (vol <= 0) {
				break;
			}
		}
		if (j == sph.size()) {
			return 1;
		}
	}

	Nodes[iNod].pt[0] = pold[0];
	Nodes[iNod].pt[1] = pold[1];
	Nodes[iNod].pt[2] = pold[2];
	return 0;
}
int DT::disturbPnt_search(const int iNod, double norm[]) {
	std::vector<int> sph;
	findSphere(iNod, sph);
	double pold[3] = { Nodes[iNod].pt[0],Nodes[iNod].pt[1],Nodes[iNod].pt[2] };
	for (int i = 0; i < 2; i++) {
		double step = 1e-3;
		while (step > 1e-10) {
			Nodes[iNod].pt[0] = pold[0] + step * norm[0];
			Nodes[iNod].pt[1] = pold[1] + step * norm[1];
			Nodes[iNod].pt[2] = pold[2] + step * norm[2];
			int j = 0;
			for (j = 0; j < sph.size(); j++) {
				if (ishulltet(sph[j]))
					continue;
				double vol = calVolume(sph[j]);
				if (vol <= 0) {
					break;
				}
			}
			step /= 2.0;
			if (j == sph.size()) {
				return 1;
			}
		}
		for (int j = 0; j < 3; j++) {
			norm[j] *= -1.0;
		}
	}

	Nodes[iNod].pt[0] = pold[0];
	Nodes[iNod].pt[1] = pold[1];
	Nodes[iNod].pt[2] = pold[2];
	return 0;
}

void DT::calGlobalEnergy(double& minV, double& maxV, int& nTet,
	double& Avg, double& Variance, double& Energy, int Energy_tpye) {
	Avg = Energy = Variance = nTet = 0;
	minV = 1e10, maxV = -1e10;
	double Sum = 0;
	int minidx = 0, maxidx = 0, n;
	for (int i = 0; i < Elems.size(); i++) {
		if (isDelEle(i) || isvirtualtet(i))
			continue;

		if (minV > Elems[i].q) {
			minV = Elems[i].q;
			minidx = i;
		}
		if (maxV < Elems[i].q) {
			maxV = Elems[i].q;
			maxidx = i;
		}
		Sum += Elems[i].q;
		nTet++;
	}
	Avg = Sum / nTet;
	for (int i = 0; i < Elems.size(); i++) {
		if (isDelEle(i) || isvirtualtet(i))
			continue;
		Variance += std::pow((Elems[i].q - Avg), 2);
	}
	Variance /= nTet;

	if (Energy_tpye < 2) {
		std::vector<int> sph;
		std::vector<double> area, Volume;
		for (int iNod = 0; iNod < Nodes.size(); iNod++) {
			if (isbndpnt(iNod) || isDelNod(iNod) || iNod == ghost)
				continue;
			findSphere(iNod, sph);
			n = sph.size();
			Volume.resize(n);
			area.resize(n);
			if (Energy_tpye == 0) {
				area.resize(0);
			}
			for (int i = 0; i < n; i++) {
				int iElem = sph[i];
				if (Elems[iElem].q < 0)
					continue;
				Volume[i] = Elems[iElem].q;

				if (Energy_tpye == 1) {
					for (int m = 0; m < 4; m++) {
						int iSecond = Elems[iElem].form[m];
						if (iSecond == iNod) {
							area[i] = calArea(Nodes[Elems[iElem].form[(m + 1) % 4]].pt,
								Nodes[Elems[iElem].form[(m + 2) % 4]].pt, Nodes[Elems[iElem].form[(m + 3) % 4]].pt);
						}
					}
				}
			}

			Energy += getVolEnergy(Volume, area);
			//printf("%d %lf\n", iNod, Energy);
		}
	}

	meshLogger->debug("minV:{:.6e} maxV:{:.6e} Num:{} Avg:{:.6e} Variance:{:.6e} Energy:{:.6e}",
		minV, maxV, nTet, Avg, Variance, Energy);
	return;
}

double DT::getVolEnergy(std::vector<double> volume, std::vector<double> area) {
	double energy = 0;
	if (area.size() == 0) {
		for (int i = 0; i < volume.size(); i++)
			energy += volume[i] * volume[i];
	}
	else {
		for (int i = 0; i < volume.size(); i++)
			energy += volume[i] * volume[i] / area[i];
	}
	return energy;
}


//#pragma optimize("",off)
int DT::sizeControlPass(double lower, double upper) {
    MeshStageLog stageLog(*this, "Size control", 2);
	int  splitNum = 0, contract = 0;
	improve_step = false;

	//split bndEdge
	if (modifyBnd) {
		for (int i = 0; i < 1; i++) {
			splitNum = splitBndEdgPass(upper);
			if (infolevel > 0) meshLogger->debug("Split    BndEdge : {}", splitNum);

			contract = contractEdgPass(lower);
			if (infolevel > 0) meshLogger->debug("Contract BndEdge : {}", contract);
			
			flipEdgPass(3);

		}
	}

	//split long Edge
	splitNum = splitLongEdgPass(upper);
	if (infolevel > 0) meshLogger->debug("Split    : {}", splitNum);

	contract = contshortEdgPass(lower);
	if (infolevel > 0) meshLogger->debug("Contract : {}", contract);

	improve_step = true;

	return splitNum + contract;
}
//#pragma optimize("",on)

int DT::splitBndEdgPass(double upper) {
    MeshStageLog stageLog(*this, "Split boundary edges", 2);
	int n = SurEdgs.size(), splitNum = 0;
	for (int i = 0; i < n; i++) {
		if (isDelSurEdg(i))
			continue;
		if (lockE.count(i))
			continue;
		int p1 = SurEdgs[i].iStart;
		int p2 = SurEdgs[i].iEnd;
		bool splitflag = false;

		if (AniSol.size() == 0) {
			double s1 = Nodes[p1].space;
			double s2 = Nodes[p2].space;
			double s_aim = (s1 + s2) / 2.0;
			double dt = distance(Nodes[p1].pt, Nodes[p2].pt);
			splitflag = dt > s_aim * upper;
		}
		else {
			double m1[6] = { 0 }, m2[6] = { 0 };
			for (int mm = 0; mm < 6; mm++) {
				m1[mm] = AniSol[p1][mm];
				m2[mm] = AniSol[p2][mm];
			}
			double  anilen = cal_ani_length(Nodes[p1].pt, Nodes[p2].pt, m1, m2);

			splitflag = anilen > upper;
			if (!splitflag && SurEdgs[i].face.size() == 2) {
				if (ifSplitBetter(i))
					splitflag = true;
			}
		}
		if (splitflag) {
			int ret = splitEdg(i);
			if (ret > 0) {
				splitNum++;;
			}
		}
	}
	return splitNum;
}

bool DT::ifSplitBetter(int index) {
	int p1 = SurEdgs[index].iStart;
	int p2 = SurEdgs[index].iEnd;

	int f1 = SurEdgs[index].face[0];
	int f2 = SurEdgs[index].face[1];
	int p3 = -1, p4 = -1;

	double pnt[3];
	for (int i = 0; i < 3; i++) {
		pnt[i] = (Nodes[p1].pt[i] + Nodes[p2].pt[i]) * 0.5;
	}

	int newp = addNode(pnt[0], pnt[1], pnt[2], 0);
	Interpolate_met(p1, p2, newp, 0.5);

	for (int j = 0; j < 3; j++) {
		if (SurTris[f1].form[j] != p1 && SurTris[f1].form[j] != p2) {
			p3 = SurTris[f1].form[j];
			break;
		}
	}
	for (int j = 0; j < 3; j++) {
		if (SurTris[f2].form[j] != p1 && SurTris[f2].form[j] != p2) {
			p4 = SurTris[f2].form[j];
			break;
		}
	}

	double AniMetric[18] = { 0 };
	for (int k = 0; k < 6; k++) {
		AniMetric[0 * 6 + k] = AniSol[p1][k];
		AniMetric[1 * 6 + k] = AniSol[p2][k];
		AniMetric[2 * 6 + k] = AniSol[p3][k];
	}

	double q1 = caltri33_ani(Nodes[p1].pt, Nodes[p2].pt, Nodes[p3].pt, AniMetric);

	for (int k = 0; k < 6; k++) {
		AniMetric[0 * 6 + k] = AniSol[p1][k];
		AniMetric[1 * 6 + k] = AniSol[p2][k];
		AniMetric[2 * 6 + k] = AniSol[p4][k];
	}

	double q2 = caltri33_ani(Nodes[p1].pt, Nodes[p2].pt, Nodes[p4].pt, AniMetric);

	for (int k = 0; k < 6; k++) {
		AniMetric[0 * 6 + k] = AniSol[p1][k];
		AniMetric[1 * 6 + k] = AniSol[newp][k];
		AniMetric[2 * 6 + k] = AniSol[p3][k];
	}

	double q11 = caltri33_ani(Nodes[p1].pt, pnt, Nodes[p3].pt, AniMetric);

	for (int k = 0; k < 6; k++) {
		AniMetric[0 * 6 + k] = AniSol[newp][k];
		AniMetric[1 * 6 + k] = AniSol[p2][k];
		AniMetric[2 * 6 + k] = AniSol[p3][k];
	}

	double q12 = caltri33_ani(pnt, Nodes[p2].pt, Nodes[p3].pt, AniMetric);

	for (int k = 0; k < 6; k++) {
		AniMetric[0 * 6 + k] = AniSol[p1][k];
		AniMetric[1 * 6 + k] = AniSol[newp][k];
		AniMetric[2 * 6 + k] = AniSol[p4][k];
	}

	double q21 = caltri33_ani(Nodes[p1].pt, pnt, Nodes[p4].pt, AniMetric);

	for (int k = 0; k < 6; k++) {
		AniMetric[0 * 6 + k] = AniSol[newp][k];
		AniMetric[1 * 6 + k] = AniSol[p2][k];
		AniMetric[2 * 6 + k] = AniSol[p4][k];
	}

	double q22 = caltri33_ani(pnt, Nodes[p2].pt, Nodes[p4].pt, AniMetric);

	DelNod(newp);

	if (std::min(q1, q2) > std::min({ q11, q12, q21, q22 }))
		return false;

	return true;
}

int DT::splitLongEdgPass_noParallel(double upper) {
	int splitNum = 0, n = Elems.size();
	for (int i = 0; i < n; i++) {
		if (isvirtualtet(i) || isDelEle(i) || ishulltet(i))
			continue;

		int longEdgidx = -1;
		double longEedgLen = 0, s_aim;
		//find longest edge
		for (int j = 0; j < 6; j++) {
			int p1 = Elems[i].form[Egid[j][0]];
			int p2 = Elems[i].form[Egid[j][1]];
			if (isBndEdg(p1, p2)) {
				continue;
			}
			double s1 = Nodes[p1].space;
			double s2 = Nodes[p2].space;
			s_aim = (s1 + s2) / 2.0;
			if (AniSol.size() == 0) {
				double dt = distance(Nodes[p1].pt, Nodes[p2].pt);
				if (dt > s_aim * upper) {
					if (dt > longEedgLen) {
						longEedgLen = dt;
						longEdgidx = j;
					}
				}
			}
			else {
				double m1[6] = { 0 }, m2[6] = { 0 };
				for (int mm = 0; mm < 6; mm++) {
					m1[mm] = AniSol[p1][mm];
					m2[mm] = AniSol[p2][mm];
				}
				double anilen = cal_ani_length(Nodes[p1].pt, Nodes[p2].pt, m1, m2);
				//if (isCornerpnt(p1) + isCornerpnt(p2) + isSegmentpnt(p1) + isSegmentpnt(p2) == 1)
				//	anilen *=1.5;
				if (anilen > upper) {
					if (anilen > longEedgLen) {
						longEedgLen = anilen;
						longEdgidx = j;
					}
				}
			}
		}
		//split it
		if (longEdgidx != -1) {
			int p1 = Elems[i].form[Egid[longEdgidx][0]];
			int p2 = Elems[i].form[Egid[longEdgidx][1]];
			double pnt[3] = { 0 };
			for (int k = 0; k < 3; k++) {
				pnt[k] = (Nodes[p1].pt[k] + Nodes[p2].pt[k]) / 2;
			}
			int newp = addNode(pnt[0], pnt[1], pnt[2], s_aim);

			Interpolate_met(p1, p2, newp, 0.5);

			std::vector<int> Shell, Shell_point;
			findShell(i, Egid[longEdgidx][0], Egid[longEdgidx][1], Shell, Shell_point);

			int ret = BW_insert_vertex(newp, Shell, 2, -1);
			if (ret <= 0) {
				DelNod(newp);
			}
			else {
				splitNum++;
				continue;
			}
		}
	}

	return splitNum;
}

//#pragma optimize("",off)
int DT::splitLongEdgPass(double upper) {
    MeshStageLog stageLog(*this, "Split long edges", 2);
    // Insertion grows shared containers and can expand its cavity beyond a
    // partition. Use the serial insertion path, as in TopologicalPass.
    return splitLongEdgPass_noParallel(upper);
}
//#pragma optimize("",on)

int DT::collapsePointLevel(int iNod) {
	if (iNod == ghost)
		return 4;
	if (isCornerpnt(iNod))
		return 3;
	if (isSegmentpnt(iNod))
		return 2;
	if (isFacetpnt(iNod))
		return 1;
	return 0;
}

double DT::contractionEdgeLength(int p1, int p2) {
	if (AniSol.empty()) {
		const double target = 0.5 * (Nodes[p1].space + Nodes[p2].space);
		if (target <= 0.0)
			return DBL_MAX;
		return distance(Nodes[p1].pt, Nodes[p2].pt) / target;
	}

	double m1[6] = { 0.0 }, m2[6] = { 0.0 };
	for (int i = 0; i < 6; ++i) {
		m1[i] = AniSol[p1][i];
		m2[i] = AniSol[p2][i];
	}
	return cal_ani_length(Nodes[p1].pt, Nodes[p2].pt, m1, m2);
}

int DT::contractEdgPass(double lower) {
    MeshStageLog stageLog(*this, "Collapse boundary edges", 2);
	const int edgeCount = static_cast<int>(SurEdgs.size());
	int contract = 0;

	for (int i = 0; i < edgeCount; ++i) {
		if (isDelSurEdg(i) || lockE.count(i))
			continue;

		const int p1 = SurEdgs[i].iStart;
		const int p2 = SurEdgs[i].iEnd;
		if (isCornerpnt(p1) && isCornerpnt(p2))
			continue;
		if (SurEdgs[i].constrain == 0 &&
			(isSegmentpnt(p1) || isCornerpnt(p1)) &&
			(isSegmentpnt(p2) || isCornerpnt(p2)))
			continue;

		const double physicalLength = distance(Nodes[p1].pt, Nodes[p2].pt);
		const bool forceByMinEdge = minEdge != -1 && physicalLength < minEdge * lower;
		if (!forceByMinEdge && contractionEdgeLength(p1, p2) >= lower)
			continue;

		if (collapseEdg(i) == 1)
			++contract;
	}

	return contract;
}

int DT::contshortEdgPass(double lower) {
    MeshStageLog stageLog(*this, "Collapse short edges", 2);
	const int elemCount = static_cast<int>(Elems.size());
	int contract = 0;

	for (int i = 0; i < elemCount; ++i) {
		if (isvirtualtet(i) || isDelEle(i) || ishulltet(i))
			continue;

		int shortestEdge = -1;
		double shortestLength = DBL_MAX;
		for (int j = 0; j < 6; ++j) {
			const int p1 = Elems[i].form[Egid[j][0]];
			const int p2 = Elems[i].form[Egid[j][1]];
			if (isBndEdg(p1, p2))
				continue;

			const double edgeLength = contractionEdgeLength(p1, p2);
			if (edgeLength < lower && edgeLength < shortestLength) {
				shortestLength = edgeLength;
				shortestEdge = j;
			}
		}

		if (shortestEdge == -1)
			continue;

		const int ia = Egid[shortestEdge][0];
		const int ib = Egid[shortestEdge][1];
		if (tryDestroyShortEdge(i, ia, ib, 1e-14) == 1)
			++contract;
	}

	return contract;
}

int DT::flipEdgPass(int nloop) {
    MeshStageLog stageLog(*this, "Boundary flips", 2);
	if (!modifyBnd) {
        if (infolevel > 0) meshLogger->info("Flip bndEdg: {}/{}", 0, 0);
		return 0;
    }

	int nsuccess = 0;
    size_t attempts = 0;
	for (int loop = 0; loop < nloop; loop++) {
		int tpsuc = 0, fail = 0;
		for (int i = 0; i < SurEdgs.size(); i++) {
			if (isDelSurEdg(i) || SurEdgs[i].info > 1 || SurEdgs[i].constrain > 0)
				continue;
			if (lockE.count(i))
				continue;
			if (!ifflipEdg(i))
				continue;

			if (periodic_P.size() != 0) {
				int p1 = SurEdgs[i].iStart;
				int p2 = SurEdgs[i].iEnd;
				if (periodic_P.find(p1) != periodic_P.end() && periodic_P.find(p2) != periodic_P.end()) {
					int eid = -1;
					for (auto it1 : periodic_P[p1]) {
						for (auto it2 : periodic_P[p2]) {
							const auto* boundaryEntry = BndEdg.find(it1, it2);
							if (boundaryEntry && isParallel(Nodes[p1].pt, Nodes[p2].pt, Nodes[it1].pt, Nodes[it2].pt)) {
								const int boundaryIndex = *boundaryEntry;
								eid = boundaryIndex;
								break;
							}
						}
						if (eid != -1) break;
					}
					if (eid != -1) {
						if (!ifflipEdg(eid))
							continue;
					}
				}
			}

			if (flipEdg(i)) {
				tpsuc++;
			}
			else {
				fail++;
			}
		}
		nsuccess += tpsuc;
        attempts += tpsuc + fail;
		if (infolevel > 0)
			meshLogger->debug("Boundary flip round: {} of {}, fail: {}", tpsuc, SurEdgs.size(), fail);
		if (tpsuc == 0)
			break;
	}

	// need to update virtual tet
	int hullidx = 0;
	for (int i = 0; i < Elems.size(); i++) {
		if (ishulltet(i)) {
			hullidx = i;
			break;
		}
	}

	int nvirtual = ColorTetNeig(hullidx, virtualID);

    if (infolevel > 0) meshLogger->info("Flip bndEdg: {}/{}", nsuccess, attempts);
	return nsuccess;
}

int DT::canDestroyShortEdge(int iElm, int ia, int ib, double Threshold) {
	if (iElm < 0 || iElm >= static_cast<int>(Elems.size()) || isDelEle(iElm) ||
		ia < 0 || ia >= 4 || ib < 0 || ib >= 4 || ia == ib)
		return 0;

	const int pa = Elems[iElm].form[ia];
	const int pb = Elems[iElm].form[ib];
	if (pa == ghost || pa == pb || isDelNod(pa) || lockV.count(pa) || (pb != ghost && isDelNod(pb)))
		return 0;

	const int levelA = collapsePointLevel(pa);
	const int levelB = collapsePointLevel(pb);
	if (levelA > levelB)
		return 0;
	if (isbndpnt(pa) && isbndpnt(pb)) {
		if (!isBndEdg(pa, pb))
			return 0;
	}
	else if (isbndpnt(pa) && !isbndpnt(pb)) {
		return 0;
	}

	double* target = Nodes[pb].pt;
	std::vector<int> shell, shellPoints, sphereA, sphereB;
	findShell(iElm, ia, ib, shell, shellPoints);
	if (shell.empty())
		return 0;
    const double shortEdgeTolerance = shortBoundaryCollapseTolerance(*this, pa, pb, shell);
    if (isCornerpnt(pa) && isCornerpnt(pb) && !(shortEdgeTolerance > 0)) return 0;


	for (int tempE : shell) {
		int ta = -1, tb = -1;
		for (int j = 0; j < 4; ++j) {
			if (Elems[tempE].form[j] == pa)
				ta = j;
			else if (Elems[tempE].form[j] == pb)
				tb = j;
		}
		if (ta == -1 || tb == -1)
			return 0;

		const int neiga = getNeig(tempE, ta);
		const int neigb = getNeig(tempE, tb);
		const int neigaOrd = getNeigOrd(tempE, ta);
		const int neigbOrd = getNeigOrd(tempE, tb);
		if (Elems[neiga].form[neigaOrd] == Elems[neigb].form[neigbOrd])
			return 0;
	}

	findSphere(pa, sphereA);
	findSphere(pb, sphereB);

	if (pb != ghost) {
		std::unordered_set<int> neighbours;
		findSphere_pnt(pa, neighbours);
		for (int neighbour : neighbours) {
			if (neighbour == pb || neighbour == ghost)
				continue;
			if (isMeshEdge(pb, neighbour))
				continue;
			if (contractionEdgeLength(pb, neighbour) > ani_upper)
				return 0;
		}
	}

	if (isbndpnt(pa)) {
		std::unordered_set<int> surfaceTriangles;
		findSphere_tri(pa, surfaceTriangles);
		const double edgeLength = distance(Nodes[pa].pt, target);
		const double deviationTolerance = std::max(0.01 * edgeLength, shortEdgeTolerance);
		const double minimumNormalCos = std::cos(5.0 * PI / 180.0);
		int previousSegmentPoint = -1;

		for (int tri : surfaceTriangles) {
			int p3 = -1, p4 = -1;
			bool containsPb = false;
			for (int j = 0; j < 3; ++j) {
				const int point = SurTris[tri].form[j];
				if (point == pa)
					continue;
				if (point == pb)
					containsPb = true;
				if (p3 == -1)
					p3 = point;
				else
					p4 = point;
			}
			if (containsPb)
				continue;
			if (p3 == -1 || p4 == -1)
				return 0;

			double oldA[3] = { Nodes[p3].pt[0] - Nodes[pa].pt[0], Nodes[p3].pt[1] - Nodes[pa].pt[1], Nodes[p3].pt[2] - Nodes[pa].pt[2] };
			double oldB[3] = { Nodes[p4].pt[0] - Nodes[pa].pt[0], Nodes[p4].pt[1] - Nodes[pa].pt[1], Nodes[p4].pt[2] - Nodes[pa].pt[2] };
			double newA[3] = { Nodes[p3].pt[0] - target[0], Nodes[p3].pt[1] - target[1], Nodes[p3].pt[2] - target[2] };
			double newB[3] = { Nodes[p4].pt[0] - target[0], Nodes[p4].pt[1] - target[1], Nodes[p4].pt[2] - target[2] };
			double oldNormal[3] = { 0.0 }, newNormal[3] = { 0.0 };
			cross(oldB, oldA, oldNormal);
			cross(newB, newA, newNormal);
			const double oldLength = lenvec(oldNormal);
			const double newLength = lenvec(newNormal);
			if (oldLength <= 1e-14 * dist_max * dist_max || newLength <= 1e-14 * dist_max * dist_max)
				return 0;
			if (dot(oldNormal, newNormal) / (oldLength * newLength) < minimumNormalCos)
				return 0;

			double planeA[3] = { Nodes[p3].pt[0] - Nodes[pa].pt[0], Nodes[p3].pt[1] - Nodes[pa].pt[1], Nodes[p3].pt[2] - Nodes[pa].pt[2] };
			double planeB[3] = { Nodes[p4].pt[0] - Nodes[pa].pt[0], Nodes[p4].pt[1] - Nodes[pa].pt[1], Nodes[p4].pt[2] - Nodes[pa].pt[2] };
			double planeNormal[3] = { 0.0 };
			cross(planeA, planeB, planeNormal);
			const double planeLength = lenvec(planeNormal);
			double displacement[3] = { target[0] - Nodes[pa].pt[0], target[1] - Nodes[pa].pt[1], target[2] - Nodes[pa].pt[2] };
			if (planeLength <= 1e-14 * dist_max * dist_max ||
				std::fabs(dot(planeNormal, displacement)) / planeLength > deviationTolerance)
				return 0;

			if (isSegmentpnt(pa)) {
				for (int point : { p3, p4 }) {
					const int edge = BndEdg.get(pa, point);
					if (edge != -1 && SurEdgs[edge].constrain > 0 && point != pb) {
						if (previousSegmentPoint == -1)
							previousSegmentPoint = point;
						else if (previousSegmentPoint != point)
							return 0;
					}
				}
			}
		}

		if (isSegmentpnt(pa) && previousSegmentPoint != -1 &&
			calArea(Nodes[previousSegmentPoint].pt, Nodes[pa].pt, target) > 1e-2)
			return 0;
	}

	if (pb != ghost) {
		for (int tempE : sphereA) {
			if (ishulltet(tempE) || std::find(shell.begin(), shell.end(), tempE) != shell.end())
				continue;
			const double orientation = dt::GEOM_FUNC::orient3d(
				Elems[tempE].form[0] == pa ? target : Nodes[Elems[tempE].form[0]].pt,
				Elems[tempE].form[1] == pa ? target : Nodes[Elems[tempE].form[1]].pt,
				Elems[tempE].form[3] == pa ? target : Nodes[Elems[tempE].form[3]].pt,
				Elems[tempE].form[2] == pa ? target : Nodes[Elems[tempE].form[2]].pt);
			if (orientation <= Threshold)
				return 0;
		}
	}

	return 1;
}

int DT::canTryDestroyShortEdge(
	int iElm,
	int ia,
	int ib,
	double Threshold,
	int* destroyIa,
	int* keepIb
) {
	if (iElm < 0 || iElm >= static_cast<int>(Elems.size()) ||
		isDelEle(iElm) ||
		ia < 0 || ia >= 4 ||
		ib < 0 || ib >= 4 ||
		ia == ib) {
		return 0;
	}

	const int pa = Elems[iElm].form[ia];
	const int pb = Elems[iElm].form[ib];

	if (pa == pb || isDelNod(pa) || isDelNod(pb))
		return 0;

	const int levelA = collapsePointLevel(pa);
	const int levelB = collapsePointLevel(pb);

	auto accept = [&](int da, int kb) -> int {
		if (!canDestroyShortEdge(iElm, da, kb, Threshold))
			return 0;

		if (destroyIa != nullptr)
			*destroyIa = da;

		if (keepIb != nullptr)
			*keepIb = kb;

		return 1;
		};

	if (levelA < levelB)
		return accept(ia, ib);

	if (levelB < levelA)
		return accept(ib, ia);

	if (accept(ia, ib))
		return 1;

	return accept(ib, ia);
}

int DT::tryDestroyShortEdge(int iElm, int ia, int ib, double Threshold) {
	if (iElm < 0 || iElm >= static_cast<int>(Elems.size()) || ia < 0 || ia >= 4 || ib < 0 || ib >= 4)
		return 0;

	const int pa = Elems[iElm].form[ia];
	const int pb = Elems[iElm].form[ib];
	const int levelA = collapsePointLevel(pa);
	const int levelB = collapsePointLevel(pb);

	if (levelA < levelB)
		return destroyShortEdge(iElm, ia, ib, Threshold);
	if (levelB < levelA)
		return destroyShortEdge(iElm, ib, ia, Threshold);

	int ret = destroyShortEdge(iElm, ia, ib, Threshold);
	if (ret == 0)
		ret = destroyShortEdge(iElm, ib, ia, Threshold);
	return ret;
}

// Directional collapse: delete the vertex at ia and keep the vertex at ib.
int DT::destroyShortEdge(int iElm, int ia, int ib, double Threshold) {
	if (iElm < 0 || iElm >= static_cast<int>(Elems.size()) || isDelEle(iElm) ||
		ia < 0 || ia >= 4 || ib < 0 || ib >= 4 || ia == ib)
		return 0;

	const int pa = Elems[iElm].form[ia];
	const int pb = Elems[iElm].form[ib];
	if (pa == ghost || pa == pb || isDelNod(pa) || lockV.count(pa) || (pb != ghost && isDelNod(pb)))
		return 0;

	const int levelA = collapsePointLevel(pa);
	const int levelB = collapsePointLevel(pb);
	if (levelA > levelB)
		return 0;
	if (isbndpnt(pa) && isbndpnt(pb)) {
		if (!isBndEdg(pa, pb))
			return 0;
	}
	else if (isbndpnt(pa) && !isbndpnt(pb)) {
		return 0;
	}

	double* target = Nodes[pb].pt;
	std::vector<int> shell, shellPoints, sphereA, sphereB;
	findShell(iElm, ia, ib, shell, shellPoints);
	if (shell.empty())
		return 0;
    const double shortEdgeTolerance = shortBoundaryCollapseTolerance(*this, pa, pb, shell);
    if (isCornerpnt(pa) && isCornerpnt(pb) && !(shortEdgeTolerance > 0)) return 0;


	// Reject non-manifold reconnections before changing any topology.
	for (int tempE : shell) {
		int ta = -1, tb = -1;
		for (int j = 0; j < 4; ++j) {
			if (Elems[tempE].form[j] == pa)
				ta = j;
			else if (Elems[tempE].form[j] == pb)
				tb = j;
		}
		if (ta == -1 || tb == -1)
			return 0;

		const int neiga = getNeig(tempE, ta);
		const int neigb = getNeig(tempE, tb);
		const int neigaOrd = getNeigOrd(tempE, ta);
		const int neigbOrd = getNeigOrd(tempE, tb);
		if (Elems[neiga].form[neigaOrd] == Elems[neigb].form[neigbOrd])
			return 0;
	}

	findSphere(pa, sphereA);
	findSphere(pb, sphereB);

	// Only pa moves. Check edges that are newly created between pb and pa's neighbours.
	if (pb != ghost) {
		std::unordered_set<int> neighbours;
		findSphere_pnt(pa, neighbours);
		for (int neighbour : neighbours) {
			if (neighbour == pb || neighbour == ghost)
				continue;
			if (isMeshEdge(pb, neighbour))
				continue;
			if (contractionEdgeLength(pb, neighbour) > ani_upper)
				return 0;
		}
	}

	// Boundary validation is directional, so a failed pa -> pb attempt can safely retry pb -> pa.
	if (isbndpnt(pa)) {
		std::unordered_set<int> surfaceTriangles;
		findSphere_tri(pa, surfaceTriangles);
		const double edgeLength = distance(Nodes[pa].pt, target);
		const double deviationTolerance = std::max(0.01 * edgeLength, shortEdgeTolerance);
		const double minimumNormalCos = std::cos(5.0 * PI / 180.0);
		int previousSegmentPoint = -1;

		for (int tri : surfaceTriangles) {
			int p3 = -1, p4 = -1;
			bool containsPb = false;
			for (int j = 0; j < 3; ++j) {
				const int point = SurTris[tri].form[j];
				if (point == pa)
					continue;
				if (point == pb)
					containsPb = true;
				if (p3 == -1)
					p3 = point;
				else
					p4 = point;
			}
			if (containsPb)
				continue;
			if (p3 == -1 || p4 == -1)
				return 0;

			double oldA[3] = { Nodes[p3].pt[0] - Nodes[pa].pt[0], Nodes[p3].pt[1] - Nodes[pa].pt[1], Nodes[p3].pt[2] - Nodes[pa].pt[2] };
			double oldB[3] = { Nodes[p4].pt[0] - Nodes[pa].pt[0], Nodes[p4].pt[1] - Nodes[pa].pt[1], Nodes[p4].pt[2] - Nodes[pa].pt[2] };
			double newA[3] = { Nodes[p3].pt[0] - target[0], Nodes[p3].pt[1] - target[1], Nodes[p3].pt[2] - target[2] };
			double newB[3] = { Nodes[p4].pt[0] - target[0], Nodes[p4].pt[1] - target[1], Nodes[p4].pt[2] - target[2] };
			double oldNormal[3] = { 0.0 }, newNormal[3] = { 0.0 };
			cross(oldB, oldA, oldNormal);
			cross(newB, newA, newNormal);
			const double oldLength = lenvec(oldNormal);
			const double newLength = lenvec(newNormal);
			if (oldLength <= 1e-14 * dist_max * dist_max || newLength <= 1e-14 * dist_max * dist_max)
				return 0;
			if (dot(oldNormal, newNormal) / (oldLength * newLength) < minimumNormalCos)
				return 0;

			double planeA[3] = { Nodes[p3].pt[0] - Nodes[pa].pt[0], Nodes[p3].pt[1] - Nodes[pa].pt[1], Nodes[p3].pt[2] - Nodes[pa].pt[2] };
			double planeB[3] = { Nodes[p4].pt[0] - Nodes[pa].pt[0], Nodes[p4].pt[1] - Nodes[pa].pt[1], Nodes[p4].pt[2] - Nodes[pa].pt[2] };
			double planeNormal[3] = { 0.0 };
			cross(planeA, planeB, planeNormal);
			const double planeLength = lenvec(planeNormal);
			double displacement[3] = { target[0] - Nodes[pa].pt[0], target[1] - Nodes[pa].pt[1], target[2] - Nodes[pa].pt[2] };
			if (planeLength <= 1e-14 * dist_max * dist_max ||
				std::fabs(dot(planeNormal, displacement)) / planeLength > deviationTolerance)
				return 0;

			if (isSegmentpnt(pa)) {
				for (int point : { p3, p4 }) {
					const int edge = BndEdg.get(pa, point);
					if (edge != -1 && SurEdgs[edge].constrain > 0 && point != pb) {
						if (previousSegmentPoint == -1)
							previousSegmentPoint = point;
						else if (previousSegmentPoint != point)
							return 0;
					}
				}
			}
		}

		if (isSegmentpnt(pa) && previousSegmentPoint != -1 &&
			calArea(Nodes[previousSegmentPoint].pt, Nodes[pa].pt, target) > 1e-2)
			return 0;
	}

	// Moving pa exactly onto pb must preserve positive volume in every tetrahedron not deleted with the shell.
	if (pb != ghost) {
		for (int tempE : sphereA) {
			if (ishulltet(tempE) || std::find(shell.begin(), shell.end(), tempE) != shell.end())
				continue;
			const double orientation = dt::GEOM_FUNC::orient3d(
				Elems[tempE].form[0] == pa ? target : Nodes[Elems[tempE].form[0]].pt,
				Elems[tempE].form[1] == pa ? target : Nodes[Elems[tempE].form[1]].pt,
				Elems[tempE].form[3] == pa ? target : Nodes[Elems[tempE].form[3]].pt,
				Elems[tempE].form[2] == pa ? target : Nodes[Elems[tempE].form[2]].pt);
			if (orientation <= Threshold)
				return 0;
		}
	}

	for (int tempE : shell) {
		int ta = -1, tb = -1;
		for (int j = 0; j < 4; ++j) {
			if (Elems[tempE].form[j] == pa)
				ta = j;
			else if (Elems[tempE].form[j] == pb)
				tb = j;
		}
		const int neiga = getNeig(tempE, ta);
		const int neigb = getNeig(tempE, tb);
		bond(neiga, getNeigOrd(tempE, ta), neigb, getNeigOrd(tempE, tb));
		DelEle(tempE);
	}

	for (int tempE : sphereA) {
		if (isDelEle(tempE))
			continue;
		for (int j = 0; j < 4; ++j) {
			if (Elems[tempE].form[j] == pa) {
				Elems[tempE].form[j] = pb;
				if (pb == ghost) {
					int points[4] = { Elems[tempE].form[(j + 1) % 4], Elems[tempE].form[(j + 2) % 4], Elems[tempE].form[(j + 3) % 4], Elems[tempE].form[j] };
					matchtet(points, tempE);
				}
				break;
			}
		}
	}

	for (int tempE : sphereA) {
		if (isDelEle(tempE))
			continue;
		for (int j = 0; j < 4; ++j)
			setP2T(Elems[tempE].form[j], tempE);
	}
	for (int tempE : sphereB) {
		if (isDelEle(tempE))
			continue;
		for (int j = 0; j < 4; ++j)
			setP2T(Elems[tempE].form[j], tempE);
	}

	for (int point : shellPoints) {
		if (isNod_in_Tet(point, Nodes[point].tet) != -1)
			continue;
		bool found = false;
		for (int tempE = 0; tempE < static_cast<int>(Elems.size()); ++tempE) {
			if (!isDelEle(tempE) && isNod_in_Tet(point, tempE) != -1) {
				Nodes[point].tet = tempE;
				found = true;
				break;
			}
		}
		if (!found)
			DelNod(point);
	}

	DelNod(pa);
	if (improve_step) {
		for (int tempE : sphereA) {
			if (!isDelEle(tempE)) updateQuality(tempE);
		}
		for (int tempE : sphereB) {
			if (!isDelEle(tempE))
				updateQuality(tempE);
		}
	}

	return 1;
}

void DT::Type_Vertex_Edg(double Angle, const Mesh& mesh) {
	// Suspended surface inspection
	for (int i = 0; i < Nodes.size(); i++) {
		if (isDelNod(i) || i == ghost)
			continue;
		if (isNod_in_Tet(i, getP2T(i)) == -1) {
			meshLogger->error("Points {} in a suspended surface, unconnected tetrahedrons!", i);
			//throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
		}
	}

	if (infolevel > 0) 
		meshLogger->debug("Conforming adaptive angle: {}", Angle);
	// Currently, we consider Angle as the ridge line.
	double angleCheck = Angle;
	for (int index = 0; index < SurEdgs.size(); index++) {
		if (isDelSurEdg(index)) {
			continue;
		}
		if (SurEdgs[index].face.size() != 2) {
			SurEdgs[index].constrain = 1;
			Nodes[SurEdgs[index].iStart].type++;
			Nodes[SurEdgs[index].iEnd].type++;
		}
		else {
			int f1 = SurEdgs[index].face[0];
			int f2 = SurEdgs[index].face[1];

			if (SurTris[f1].parent != SurTris[f2].parent)
			{
				SurEdgs[index].constrain = 1;
				Nodes[SurEdgs[index].iStart].type++;
				Nodes[SurEdgs[index].iEnd].type++;
			}
			else {
				//base on angle
				int p1 = SurEdgs[index].iStart;
				int p2 = SurEdgs[index].iEnd;
				int p3 = -1, p4 = -1;

				for (int i = 0; i < 3; i++) {
					if (SurTris[f1].form[i] != p1 && SurTris[f1].form[i] != p2) {
						p3 = SurTris[f1].form[i];
						break;
					}
				}
				for (int i = 0; i < 3; i++) {
					if (SurTris[f2].form[i] != p1 && SurTris[f2].form[i] != p2) {
						p4 = SurTris[f2].form[i];
						break;
					}
				}

				double n1[3] = { 0 }, n2[3] = { 0 }, Dihedral = 0;
				calnormal(p1, p2, p3, n1);
				calnormal(p1, p4, p2, n2);

				double aa = dot(n1, n2) / (lenvec(n1) * lenvec(n2));
				if (abs(aa) <= 1 && (180 - RADIO2ANGLE(acos(aa)) < Angle)) {
					SurEdgs[index].constrain = 1;
					Nodes[SurEdgs[index].iStart].type++;
					Nodes[SurEdgs[index].iEnd].type++;
				}
			}
		}
	}

	//Determine by input BndEdg
	if (mesh.S.size() != 0) {
		for (int i = 0; i < mesh.S.size(); i++) {
			int p1 = mesh.S[i][0];
			int p2 = mesh.S[i][1];
			int index = BndEdg.get(p1, p2);
			if (index == -1)
				continue;
			if (isDelSurEdg(index)) {
				continue;
			}

			SurEdgs[index].geo = mesh.S[i][2];
			if (SurEdgs[index].constrain == 1) {
				SurEdgs[index].constrain = 2;
				//is constrained bnd edge
				continue;
			}

			SurEdgs[index].constrain = 2;
			Nodes[p1].type++;
			Nodes[p2].type++;
		}
	}

	std::vector<std::unordered_set<int>> pp_neig;
	pp_neig.resize(Nodes.size());
	for (int index = 0; index < SurEdgs.size(); index++) {
		int p1 = SurEdgs[index].iStart;
		int p2 = SurEdgs[index].iEnd;
		pp_neig[p1].insert(p2);
		pp_neig[p2].insert(p1);
	}

	//determine by angle
	for (int iNod = 0; iNod < Nodes.size(); iNod++) {
		if (isSegmentpnt(iNod)) {
			std::unordered_map<int, bool> seen;
			const auto& neig_p = pp_neig[iNod];
			/*std::unordered_set<int> neig_p;
			findSphere_tri_p(iNod, neig_p);*/

			for (auto it : neig_p) {
				if (isCornerpnt(it) || isSegmentpnt(it))
				{
if (auto* boundaryEntry = BndEdg.find(iNod, it))
					{
						const int boundaryIndex = *boundaryEntry;
						int Edgidx = boundaryIndex;
						if (SurEdgs[Edgidx].constrain>0) {
							if (seen.find(it) == seen.end())
							{
								seen[it] = true;
							}
						}
					}
				}
			}
			int p1 = -1, p2 = -1;
			for (auto it : seen) {
				if (p1 == -1) {
					p1 = it.first;
				}
				else {
					p2 = it.first;
				}
			}

			double iNodp1[3] = { Nodes[p1].pt[0] - Nodes[iNod].pt[0],  Nodes[p1].pt[1] - Nodes[iNod].pt[1],  Nodes[p1].pt[2] - Nodes[iNod].pt[2] };
			double iNodp2[3] = { Nodes[p2].pt[0] - Nodes[iNod].pt[0],  Nodes[p2].pt[1] - Nodes[iNod].pt[1],  Nodes[p2].pt[2] - Nodes[iNod].pt[2] };
			double angle = dot(iNodp1, iNodp2) / (lenvec(iNodp1) * lenvec(iNodp2));
			if (angle > -0.996/*std::cos(ANGLE2RADIO(Angle))*/) {
				Nodes[iNod].type++;
			}
		}
	}

	return;
}

#pragma optimize("",off)
int DT::splitEdg(int index,int deep) {
	//split Bnd Edge
	int i, j, newp, srchtet = -1, p1, p2;
	double pnt[3] = { 0 }, space = 0;
	p1 = SurEdgs[index].iStart;
	p2 = SurEdgs[index].iEnd;

	std::vector<int> sph;
	findSphere(p1, sph);
	for (i = 0; i < sph.size(); i++) {
		if (isNod_in_Tet(p2, sph[i]) != -1) {
			srchtet = sph[i];
			break;
		}
	}

	if (srchtet == -1)
		return 0;

	int ia, ib;
	for (int m = 0; m < 4; m++) {
		if (Elems[srchtet].form[m] == p1) {
			ia = m;
		}
		if (Elems[srchtet].form[m] == p2) {
			ib = m;
		}
	}

	std::vector<int> shell, shellp;
	findShell(srchtet, ia, ib, shell, shellp);

	int  manifold;
	std::vector<int> f;

	manifold = SurEdgs[index].face.size();//might non mainfold happen
	f.resize(manifold * 3);
	for (i = 0; i < manifold; i++) {
		f[i] = SurEdgs[index].face[i];
		FacetNum[SurTris[f[i]].parent]--;
		int p3 = -1;
		for (j = 0; j < 3; j++) {
			if (SurTris[f[i]].form[j] != p1 && SurTris[f[i]].form[j] != p2) {
				p3 = SurTris[f[i]].form[j];
				break;
			}
		}
		if (std::find(shellp.begin(), shellp.end(),p3) == shellp.end()) {
			//printf("%d %d %d\n", BndEdg.get(p1, p2), BndEdg.get(p1, p3), BndEdg.get(p3, p2));
			//printf("%d %d %d %d\n", isMeshEdge(p1, p2) ? 1 : 0, isMeshEdge(p1, p3) ? 1 : 0, isMeshEdge(p3, p2) ? 1 : 0, isMeshFace(p1,p3, p2) ? 1 : 0);
			//checkMeshError();
			meshLogger->warn("Warning: Split Edge, {} {} {}", p1, p2,p3);
			return 0;
		}
	}

	initBWshell = shell.size();
	double oldminVolume_bw = minVolume_bw;
	double shellVol = DBL_MAX;
	for (int i = 0; i < shell.size(); i++) {
		if (ishulltet(shell[i]))
			continue;
		shellVol = std::min(calVolume(shell[i]), shellVol);
	}
	minVolume_bw = std::max(shellVol / 10.0, minVolume_bw);

	for (i = 0; i < 3; i++) {
		pnt[i] = (Nodes[p1].pt[i] + Nodes[p2].pt[i]) * 0.5;
	}
	space = Nodes[p1].space * 0.5 + Nodes[p2].space * 0.5;
	newp = addNode(pnt[0], pnt[1], pnt[2], space);

	setbndpnt(newp);
	Interpolate_met(p1, p2, newp, 0.5);

	//before add steiner point,we need to deltet BndTri connect to Edge
	for (i = 0; i < SurEdgs[index].face.size(); i++) {
		int Fidx = SurEdgs[index].face[i];
		BndTri.erase(SurTris[Fidx].form[0], SurTris[Fidx].form[1], SurTris[Fidx].form[2]);
	}
	BndEdg.erase(p1, p2);

	//Boundary edge splitting allows for lower constraints.
	if (BW_insert_vertex(newp, shell, 3) != 1) {
		minVolume_bw = oldminVolume_bw;
		//recover edge min constraints
		DelNod(newp);
		for (i = 0; i < SurEdgs[index].face.size(); i++) {
			int Fidx = SurEdgs[index].face[i];
			BndTri.add(SurTris[Fidx].form[0], SurTris[Fidx].form[1], SurTris[Fidx].form[2], Fidx);
		}
		BndEdg.add(p1, p2, index);
		//meshLogger->debug("Split Bnd Edge {} failed!", index);
		return 0;
	}
	minVolume_bw = oldminVolume_bw;

	//set new point type
	// F=1 S=3 C=else
	// F F = F
	// F S = F
	// F C = F
	// S S = S
	// S C = S
	// C C = S
	if (SurEdgs[index].constrain>0) {
		Nodes[newp].type = 3;
	}

	//recover edge min constraints
	for (i = 0; i < SurEdgs[index].face.size(); i++) {
		int Fidx = SurEdgs[index].face[i];
		setDelSurTri(Fidx);
	}
	setDelSurEdg(index);
	/********************** update bnd info ***********************/
	int  e1, e2;
	std::vector<int> e;

	e.resize(manifold);
	for (i = manifold; i < manifold * 3; i++) {
		f[i] = SurTris.size();
		SurTris.emplace_back(SurTri());//add new surTri
	}

	//add [p1,newp]
	e1 = SurEdgs.size();
	BndEdg.add(p1, newp, e1);
	SurEdgs.emplace_back(SurEdg());
	SurEdgs[e1].iStart = p1;
	SurEdgs[e1].iEnd = newp;
	SurEdgs[e1].constrain = SurEdgs[index].constrain;
	SurEdgs[e1].geo = SurEdgs[index].geo;
	SurEdgs[e1].info = 0;//don't recover
	SurEdgs[e1].face.resize(manifold);
	for (i = 0; i < manifold; i++)
		SurEdgs[e1].face[i] = f[manifold + i * 2];//first face

	//add [p2,newp]
	e2 = SurEdgs.size();
	BndEdg.add(p2, newp, e2);
	SurEdgs.emplace_back(SurEdg());
	SurEdgs[e2].iStart = p2;
	SurEdgs[e2].iEnd = newp;
	SurEdgs[e2].constrain = SurEdgs[index].constrain;
	SurEdgs[e2].geo = SurEdgs[index].geo;
	SurEdgs[e2].info = 0;//don't recover
	SurEdgs[e2].face.resize(manifold);
	for (i = 0; i < manifold; i++)
		SurEdgs[e2].face[i] = f[manifold + i * 2 + 1];//first face

	//add new edge [newp,f[i].form[x]!=p1,p2]
	for (i = 0; i < manifold; i++) {
		int p3 = -1;
		for (j = 0; j < 3; j++) {
			if (SurTris[f[i]].form[j] != p1 && SurTris[f[i]].form[j] != p2) {
				p3 = SurTris[f[i]].form[j];
				break;
			}
		}
		e[i] = SurEdgs.size();
		BndEdg.add(p3, newp, e[i]);
		SurEdgs.emplace_back(SurEdg());
		SurEdgs[e[i]].iStart = p3;
		SurEdgs[e[i]].iEnd = newp;
		//SurEdgs[e[i]].parent = index;
		SurEdgs[e[i]].info = 0;//don't recover
		SurEdgs[e[i]].face.resize(2);
		SurEdgs[e[i]].face[0] = f[manifold + i * 2];//first face
		SurEdgs[e[i]].face[1] = f[manifold + i * 2 + 1];//first face

		if (!isMeshEdge(p3, newp) && isNod_in_Tet(p3, getP2T(p3)) != -1) {
			int ret = recoverEdge(e[i], 1, 0);
			if (ret != 1)
				throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
		}
	}
	//updata surTri connect p1
	for (i = manifold; i < manifold * 3; i = i + 2) {
		for (j = 0; j < 3; j++) {
			int parentform = SurTris[f[(i - manifold) / 2]].form[j];
			SurTris[f[i]].form[j] = parentform == p2 ? newp : parentform;
		}
		SurTris[f[i]].parent = f[(i - manifold) / 2];
		for (j = 0; j < 3; j++) {
			int tempE = BndEdg.get(SurTris[f[i]].form[(j + 1) % 3], SurTris[f[i]].form[(j + 2) % 3]);
			//SurTris[f[i]].edgs[j] = tempE;
			if (SurTris[f[i]].form[(j + 1) % 3] != newp && SurTris[f[i]].form[(j + 2) % 3] != newp) {
				for (int k = 0; k < SurEdgs[tempE].face.size(); k++) {
					if (SurEdgs[tempE].face[k] == SurTris[f[i]].parent) {
						SurEdgs[tempE].face[k] = f[i];
						break;
					}
				}
			}
		}
		SurTris[f[i]].info = 0;
		SurTris[f[i]].parent = SurTris[SurTris[f[i]].parent].parent;
		FacetNum[SurTris[f[i]].parent]++;
		BndTri.add(SurTris[f[i]].form[0], SurTris[f[i]].form[1], SurTris[f[i]].form[2], f[i]);//add this sub bnd Tris
	}
	//updata surTri connect p2
	for (i = manifold + 1; i < manifold * 3; i = i + 2) {
		for (j = 0; j < 3; j++) {
			int parentform = SurTris[f[(i - manifold - 1) / 2]].form[j];
			SurTris[f[i]].form[j] = parentform == p1 ? newp : parentform;
		}
		SurTris[f[i]].parent = f[(i - manifold - 1) / 2];
		for (j = 0; j < 3; j++) {
			int tempE = BndEdg.get(SurTris[f[i]].form[(j + 1) % 3], SurTris[f[i]].form[(j + 2) % 3]);
			//SurTris[f[i]].edgs[j] = tempE;
			if (SurTris[f[i]].form[(j + 1) % 3] != newp && SurTris[f[i]].form[(j + 2) % 3] != newp) {
				for (int k = 0; k < SurEdgs[tempE].face.size(); k++) {
					if (SurEdgs[tempE].face[k] == SurTris[f[i]].parent) {
						SurEdgs[tempE].face[k] = f[i];
						break;
					}
				}
			}
		}
		SurTris[f[i]].info = 0;
		SurTris[f[i]].parent = SurTris[SurTris[f[i]].parent].parent;
		FacetNum[SurTris[f[i]].parent]++;
		BndTri.add(SurTris[f[i]].form[0], SurTris[f[i]].form[1], SurTris[f[i]].form[2], f[i]);//add this sub bnd Tris
	}

	//return new point idx
	//	smoothBndPnt_ani(newp);
	//}

	if (/*deep == 0 &&*/periodic_P.size() != 0) {
		if (periodic_P.find(p1) != periodic_P.end() && periodic_P.find(p2) != periodic_P.end()) {
			const auto partners1 = periodic_P[p1];
			const auto partners2 = periodic_P[p2];
			std::unordered_set<int> processedEdges;
			for (auto it1 : partners1) {
				if (it1 == p2)
					continue;
				for (auto it2 : partners2) {
					if (it2 == p1)
						continue;
					if (isParallel(Nodes[p1].pt, Nodes[p2].pt, Nodes[it1].pt, Nodes[it2].pt)) {
						if (auto* boundaryEntry = BndEdg.find(it1, it2)) {
							const int boundaryIndex = *boundaryEntry;
							int eid = boundaryIndex;
							if (!processedEdges.insert(eid).second)
								continue;
							int ret = splitEdg(eid, 1);
							if (ret != 0) {
								if (std::find(periodic_P[ret].begin(), periodic_P[ret].end(), newp) == periodic_P[ret].end())
									periodic_P[ret].push_back(newp);

								if (std::find(periodic_P[newp].begin(), periodic_P[newp].end(), ret) == periodic_P[newp].end())
									periodic_P[newp].push_back(ret);

								continue; // Synchronize every matching periodic edge.
							}
							else {

								int i, j, srchtet = -1, p1, p2;
								double pnt[3] = { 0 }, space = 0;
								p1 = SurEdgs[eid].iStart;
								p2 = SurEdgs[eid].iEnd;

								std::vector<int> sph;
								findSphere(p1, sph);
								for (i = 0; i < sph.size(); i++) {
									if (isNod_in_Tet(p2, sph[i]) != -1) {
										srchtet = sph[i];
										break;
									}
								}

								if (srchtet == -1)
									return 0;

								int ia, ib;
								for (int m = 0; m < 4; m++) {
									if (Elems[srchtet].form[m] == p1) {
										ia = m;
									}
									if (Elems[srchtet].form[m] == p2) {
										ib = m;
									}
								}

								std::vector<int> shell, shellp;
								findShell(srchtet, ia, ib, shell, shellp);

								initBWshell = shell.size();
								double oldminVolume_bw = minVolume_bw;
								double shellVol = DBL_MAX;
								for (int i = 0; i < shell.size(); i++) {
									if (ishulltet(shell[i]))
										continue;
								}
								// Explicit mesh export is left to the caller.
								meshLogger->warn("The periodicity may be disrupted during splitEdg, {} {}", p1, p2);
							}
						}
						else {
							std::unordered_set<int> Sphere_pnt;
							findSphere_pnt(it1, Sphere_pnt);
							for (auto midp : Sphere_pnt) {
								double mid12[3] = {(Nodes[it1].pt[0] + Nodes[it2].pt[0]) / 2.0,(Nodes[it1].pt[1] + Nodes[it2].pt[1]) / 2.0,(Nodes[it1].pt[2] + Nodes[it2].pt[2]) / 2.0};
								if (distance2(Nodes[midp].pt, mid12)<1e-10) {
									if (std::find(periodic_P[midp].begin(), periodic_P[midp].end(), newp) == periodic_P[midp].end())
										periodic_P[midp].push_back(newp);

									if (std::find(periodic_P[newp].begin(), periodic_P[newp].end(), midp) == periodic_P[newp].end()) 
										periodic_P[newp].push_back(midp);
								
									break;
								}
							}
						}
					}
				}
			}
		}
	}

	return newp;
}
#pragma optimize("",on)

int DT::collapseEdg(int index) {
	return collapseEdg(index, -1, -1);
}

int DT::canCollapseBoundaryEdgeDirected(int edgeIndex, int deletePoint, int keepPoint) {
	if (edgeIndex < 0 || edgeIndex >= static_cast<int>(SurEdgs.size()) || isDelSurEdg(edgeIndex) || lockE.count(edgeIndex))
		return 0;

	if (deletePoint < 0 || keepPoint < 0 || deletePoint == keepPoint)
		return 0;

	if (isDelNod(deletePoint) || isDelNod(keepPoint))
		return 0;

	if (!isBndEdg(deletePoint, keepPoint))
		return 0;

	int edgeTet = -1;
	std::vector<int> edgeSphere;
	findSphere(deletePoint, edgeSphere);

	for (int tet : edgeSphere) {
		if (isDelEle(tet) || isvirtualtet(tet) || ishulltet(tet))
			continue;

		if (isNod_in_Tet(keepPoint, tet) != -1) {
			edgeTet = tet;
			break;
		}
	}

	if (edgeTet == -1)
		return 0;

	int eia = -1;
	int eib = -1;

	for (int j = 0; j < 4; ++j) {
		if (Elems[edgeTet].form[j] == deletePoint)
			eia = j;
		else if (Elems[edgeTet].form[j] == keepPoint)
			eib = j;
	}

	if (eia == -1 || eib == -1)
		return 0;

	// Directional check: deletePoint -> keepPoint only.
	return canDestroyShortEdge(edgeTet, eia, eib, 1e-14);
}

int DT::collapseEdg(int index, int expectedDelete, int expectedKeep) {
	if (index < 0 || index >= static_cast<int>(SurEdgs.size()) || isDelSurEdg(index) || lockE.count(index))
		return 0;

	int p1 = SurEdgs[index].iStart;
	int p2 = SurEdgs[index].iEnd;

	const int oldP1 = p1;
	const int oldP2 = p2;

	double oldP1Pt[3] = {
		Nodes[p1].pt[0],
		Nodes[p1].pt[1],
		Nodes[p1].pt[2]
	};

	double oldP2Pt[3] = {
		Nodes[p2].pt[0],
		Nodes[p2].pt[1],
		Nodes[p2].pt[2]
	};

	const bool directionalCollapse = (expectedDelete != -1 || expectedKeep != -1);

	if (directionalCollapse) {
		if (expectedDelete == -1 || expectedKeep == -1)
			return 0;

		if (expectedDelete == expectedKeep)
			return 0;

		if (isDelNod(expectedDelete) || isDelNod(expectedKeep))
			return 0;

		if (!isBndEdg(expectedDelete, expectedKeep))
			return 0;
	}

	// ------------------------------------------------------------
	// 1. Locate one tetrahedron containing the boundary edge
	// ------------------------------------------------------------
	int sourceTet = -1;
	std::vector<int> sphere;
	findSphere(p1, sphere);

	for (int tet : sphere) {
		if (isDelEle(tet) || isvirtualtet(tet) || ishulltet(tet))
			continue;

		if (isNod_in_Tet(p2, tet) != -1) {
			sourceTet = tet;
			break;
		}
	}

	if (sourceTet == -1)
		return 0;

	int ia = -1;
	int ib = -1;

	for (int j = 0; j < 4; ++j) {
		if (Elems[sourceTet].form[j] == p1)
			ia = j;
		else if (Elems[sourceTet].form[j] == p2)
			ib = j;
	}

	if (ia == -1 || ib == -1)
		return 0;

	// ------------------------------------------------------------
	// 2. Decide collapse direction
	// ------------------------------------------------------------
	int destroyIa = -1;
	int keepIb = -1;

	if (directionalCollapse) {
		for (int j = 0; j < 4; ++j) {
			if (Elems[sourceTet].form[j] == expectedDelete)
				destroyIa = j;
			else if (Elems[sourceTet].form[j] == expectedKeep)
				keepIb = j;
		}

		if (destroyIa == -1 || keepIb == -1)
			return 0;

		// Directional check: expectedDelete -> expectedKeep only.
		if (!canDestroyShortEdge(sourceTet, destroyIa, keepIb, 1e-14))
			return 0;
	}
	else {
		if (!canTryDestroyShortEdge(sourceTet, ia, ib, 1e-14, &destroyIa, &keepIb))
			return 0;
	}

	const int collapseDelete = Elems[sourceTet].form[destroyIa];
	const int collapseKeep = Elems[sourceTet].form[keepIb];

	// ------------------------------------------------------------
	// 3. Find required periodic paired collapse
	//    Only the master/original collapse searches its periodic pair.
	//    Directional slave collapse must not recurse again.
	// ------------------------------------------------------------
	// Store every paired edge and its required delete -> keep direction.
	std::vector<std::array<int, 3>> requiredPeriodicEdges;
	std::unordered_set<int> periodicEdgeIds;

	if (!directionalCollapse && periodic_P.size() != 0) {
		auto deletePeriodic = periodic_P.find(collapseDelete);
		auto keepPeriodic = periodic_P.find(collapseKeep);

		const bool deleteHasPeriodic = (deletePeriodic != periodic_P.end());
		const bool keepHasPeriodic = (keepPeriodic != periodic_P.end());

		// One endpoint has periodic partner but the other does not:
		// reject to avoid breaking periodic consistency.
		if (deleteHasPeriodic != keepHasPeriodic)
			return 0;

		if (deleteHasPeriodic && keepHasPeriodic) {
			for (int it1 : deletePeriodic->second) {
				if (it1 == collapseKeep || isDelNod(it1))
					continue;

				for (int it2 : keepPeriodic->second) {
					if (it2 == collapseDelete || it1 == it2 || isDelNod(it2))
						continue;

					if (!isParallel(
						Nodes[collapseDelete].pt,
						Nodes[collapseKeep].pt,
						Nodes[it1].pt,
						Nodes[it2].pt
					)) {
						continue;
					}

					const auto* boundaryEntry = BndEdg.find(it1, it2);
					if (!boundaryEntry)
						continue;

					int eid = *boundaryEntry;

					if (!periodicEdgeIds.insert(eid).second)
						continue;

					// Every matching replica must pass before the original is changed.
					if (!canCollapseBoundaryEdgeDirected(eid, it1, it2))
						return 0;
					requiredPeriodicEdges.push_back({ eid, it1, it2 });
				}
			}

			if (requiredPeriodicEdges.empty())
				return 0;
		}
	}

	// ------------------------------------------------------------
	// 4. Collect boundary faces incident to a point
	// ------------------------------------------------------------
	auto collectBoundaryFaces = [&](int point) {
		std::vector<int> faces;
		std::vector<int> pointSphere;

		findSphere(point, pointSphere);

		for (int tet : pointSphere) {
			if (isDelEle(tet))
				continue;

			for (int j = 0; j < 4; ++j) {
				if (Elems[tet].form[j] == point)
					continue;

				int a, b, c, d;
				DFC(j, a, b, c, d);

				b = Elems[tet].form[b];
				c = Elems[tet].form[c];
				d = Elems[tet].form[d];

				if (auto* boundaryEntry = BndTri.find(b, c, d)) {
					const int boundaryIndex = *boundaryEntry;
					int fid = boundaryIndex;
					if (fid >= 0 && !isDelSurTri(fid))
						faces.push_back(fid);
				}
			}
		}

		return faces;
		};

	std::vector<int> surfaceFacesP1 = collectBoundaryFaces(p1);
	std::vector<int> surfaceFacesP2 = collectBoundaryFaces(p2);

	// ------------------------------------------------------------
	// 5. Execute collapse
	//    Critical change:
	//    use destroyShortEdge(sourceTet, destroyIa, keepIb),
	//    not tryDestroyShortEdge(sourceTet, ia, ib).
	// ------------------------------------------------------------
	if (destroyShortEdge(sourceTet, destroyIa, keepIb, 1e-14) == 0)
		return 0;

	std::vector<int> surfaceFaces;

	// After destroyShortEdge, exactly one endpoint should be deleted.
	// Keep the original convention below:
	// p1 = deleted point, p2 = survivor.
	if (isDelNod(p2)) {
		surfaceFaces.swap(surfaceFacesP2);
		std::swap(p1, p2);
	}
	else {
		surfaceFaces.swap(surfaceFacesP1);
	}

	if (!isDelNod(p1) || isDelNod(p2))
		return 0;

	double collapseP1Pt[3] = { 0.0, 0.0, 0.0 };
	double collapseP2Pt[3] = { 0.0, 0.0, 0.0 };

	for (int k = 0; k < 3; ++k) {
		collapseP1Pt[k] = (p1 == oldP1) ? oldP1Pt[k] : oldP2Pt[k];
		collapseP2Pt[k] = (p2 == oldP2) ? oldP2Pt[k] : oldP1Pt[k];
	}

	// ------------------------------------------------------------
	// 6. Remove the two boundary faces incident to the collapsed edge
	// ------------------------------------------------------------
	for (int fid : SurEdgs[index].face) {
		if (isDelSurTri(fid))
			continue;

		BndTri.erase(
			SurTris[fid].form[0],
			SurTris[fid].form[1],
			SurTris[fid].form[2]
		);

		setDelSurTri(fid);

		int leftEdge = -1;
		int rightEdge = -1;
		std::vector<int> neighbouringFaces;

		for (int j = 0; j < 3; ++j) {
			if (SurTris[fid].form[j] == p2) {
				rightEdge = BndEdg.get(
					SurTris[fid].form[(j + 1) % 3],
					SurTris[fid].form[(j + 2) % 3]
				);

				setDelSurEdg(rightEdge);
				BndEdg.erase(SurEdgs[rightEdge].iStart, SurEdgs[rightEdge].iEnd);

				for (int face : SurEdgs[rightEdge].face) {
					if (face != fid)
						neighbouringFaces.push_back(face);
				}
			}
			else if (SurTris[fid].form[j] == p1) {
				leftEdge = BndEdg.get(
					SurTris[fid].form[(j + 1) % 3],
					SurTris[fid].form[(j + 2) % 3]
				);

				for (int face : SurEdgs[leftEdge].face) {
					if (face != fid)
						neighbouringFaces.push_back(face);
				}
			}
		}

		if (leftEdge == -1 || rightEdge == -1) {
			throw EXCEPTIONSTRING(
				std::string("error exit in") +
				std::string(__FILE__) +
				std::to_string(__LINE__)
			);
		}

		if (SurEdgs[leftEdge].constrain == 0 && SurEdgs[rightEdge].constrain > 0)
			SurEdgs[leftEdge].constrain = SurEdgs[rightEdge].constrain;

		SurEdgs[leftEdge].face = neighbouringFaces;
	}

	setDelSurEdg(index);
	BndEdg.erase(SurEdgs[index].iStart, SurEdgs[index].iEnd);

	// ------------------------------------------------------------
	// 7. Redirect remaining boundary faces and edges
	//    from deleted point p1 to survivor p2
	// ------------------------------------------------------------
	for (int fid : surfaceFaces) {
		if (isDelSurTri(fid))
			continue;

		BndTri.erase(
			SurTris[fid].form[0],
			SurTris[fid].form[1],
			SurTris[fid].form[2]
		);

		int pa = -1;
		int pb = -1;

		for (int j = 0; j < 3; ++j) {
			if (SurTris[fid].form[j] == p1) {
				SurTris[fid].form[j] = p2;
				pa = SurTris[fid].form[(j + 1) % 3];
				pb = SurTris[fid].form[(j + 2) % 3];
				break;
			}
		}

		BndTri.add(
			SurTris[fid].form[0],
			SurTris[fid].form[1],
			SurTris[fid].form[2],
			fid
		);

		for (int neighbour : { pa, pb }) {
			if (neighbour == -1)
				continue;

			const auto* boundaryEntry = BndEdg.find(p1, neighbour);
			if (!boundaryEntry)
				continue;

			const int edge = *boundaryEntry;

			BndEdg.erase(p1, neighbour);

			if (SurEdgs[edge].iStart == p1)
				SurEdgs[edge].iStart = p2;
			else if (SurEdgs[edge].iEnd == p1)
				SurEdgs[edge].iEnd = p2;

			BndEdg.add(SurEdgs[edge].iStart, SurEdgs[edge].iEnd, edge);
		}
	}

	// ------------------------------------------------------------
	// 8. Update periodic map
	//    Only the master/original collapse updates periodic_P.
	//    Directional paired collapse should not update it independently.
	// ------------------------------------------------------------
	if (!directionalCollapse && periodic_P.size() != 0) {
		auto erasePeriodicValue = [](std::vector<int>& values, int value) {
			values.erase(
				std::remove(values.begin(), values.end(), value),
				values.end()
			);
			};

		auto removePeriodicPoint = [&](int point) {
			auto it = periodic_P.find(point);

			if (it != periodic_P.end()) {
				std::vector<int> partners = it->second;

				for (int partner : partners) {
					auto partnerIt = periodic_P.find(partner);

					if (partnerIt == periodic_P.end())
						continue;

					erasePeriodicValue(partnerIt->second, point);

					if (partnerIt->second.empty())
						periodic_P.erase(partnerIt);
				}

				periodic_P.erase(point);
			}

			for (auto mapIt = periodic_P.begin(); mapIt != periodic_P.end();) {
				erasePeriodicValue(mapIt->second, point);

				if (mapIt->second.empty())
					mapIt = periodic_P.erase(mapIt);
				else
					++mapIt;
			}
			};

		auto addPeriodicPair = [&](int pa, int pb) {
			if (pa == pb || isDelNod(pa) || isDelNod(pb))
				return;

			auto& va = periodic_P[pa];
			if (std::find(va.begin(), va.end(), pb) == va.end())
				va.push_back(pb);

			auto& vb = periodic_P[pb];
			if (std::find(vb.begin(), vb.end(), pa) == vb.end())
				vb.push_back(pa);
			};

		std::vector<int> periodicDeletedPoints;
		std::vector<int> periodicSurvivors;
		for (const auto& pairedEdge : requiredPeriodicEdges) {
			const int requiredPeriodicEdge = pairedEdge[0];
			const int requiredPeriodicP1 = pairedEdge[1];
			const int requiredPeriodicP2 = pairedEdge[2];
			int ret = collapseEdg(requiredPeriodicEdge, requiredPeriodicP1, requiredPeriodicP2);

			if (ret != 0 && isDelNod(requiredPeriodicP1) && !isDelNod(requiredPeriodicP2)) {
				periodicDeletedPoints.push_back(requiredPeriodicP1);
				periodicSurvivors.push_back(requiredPeriodicP2);
			}
			else {
				meshLogger->warn(
					"The periodicity may be disrupted during collapseEdg, {} {}",
					requiredPeriodicP1,
					requiredPeriodicP2
				);
			}
		}

		// Keep the map intact until all directional collapses have been attempted.
		removePeriodicPoint(p1);
		for (int point : periodicDeletedPoints)
			removePeriodicPoint(point);
		for (int survivor : periodicSurvivors)
			addPeriodicPair(p2, survivor);
	}

	return 1;
}

int DT::ifflipEdg(int index) {
	int mainfold = SurEdgs[index].face.size();
	if (mainfold != 2)
		return 0;

	int p1 = SurEdgs[index].iStart;
	int p2 = SurEdgs[index].iEnd;
	int f1 = SurEdgs[index].face[0];
	int f2 = SurEdgs[index].face[1];
	int p3 = -1, p4 = -1;

	for (int i = 0; i < 3; i++) {
		if (SurTris[f1].form[i] != p1 && SurTris[f1].form[i] != p2) {
			p3 = SurTris[f1].form[i];
			break;
		}
	}
	for (int i = 0; i < 3; i++) {
		if (SurTris[f2].form[i] != p1 && SurTris[f2].form[i] != p2) {
			p4 = SurTris[f2].form[i];
			break;
		}
	}

	if (!checkFlipNormal(p1, p2, p3, p4))
		return 0;

	double disedge = segmentSegmentDistance(Nodes[p1].pt, Nodes[p2].pt, Nodes[p3].pt, Nodes[p4].pt);
	if (disedge > 1e-1 * distance(Nodes[p1].pt,Nodes[p2].pt))
		return 0;
	
	// angle check
	std::vector<double> init1, init2, after1, after2;
	double min1 = calculateTriangleAngles(Nodes[p1].pt, Nodes[p2].pt, Nodes[p3].pt, init1);
	double min2 = calculateTriangleAngles(Nodes[p1].pt, Nodes[p2].pt, Nodes[p4].pt, init2);

	double min3 = calculateTriangleAngles(Nodes[p1].pt, Nodes[p3].pt, Nodes[p4].pt, after1);
	double min4 = calculateTriangleAngles(Nodes[p2].pt, Nodes[p3].pt, Nodes[p4].pt, after2);

	if (init1[0] + init2[0] >= 170 || init1[1] + init2[1] >= 170) //179.5
		return 0;

	if (AniSol.size() == 0) {
		if (std::min(min3, min4) <= std::min(min1, min2) + 0.1)
			return 0;
	}
	else {
		double AniMetric[18] = { 0 };
		for (int k = 0; k < 6; k++) {
			AniMetric[0 * 6 + k] = AniSol[p1][k];
			AniMetric[1 * 6 + k] = AniSol[p2][k];
			AniMetric[2 * 6 + k] = AniSol[p3][k];
		}

		double q1 = caltri33_ani(Nodes[p1].pt, Nodes[p2].pt, Nodes[p3].pt, AniMetric);

		for (int k = 0; k < 6; k++) {
			AniMetric[0 * 6 + k] = AniSol[p1][k];
			AniMetric[1 * 6 + k] = AniSol[p2][k];
			AniMetric[2 * 6 + k] = AniSol[p4][k];
		}

		double q2 = caltri33_ani(Nodes[p1].pt, Nodes[p2].pt, Nodes[p4].pt, AniMetric);

		for (int k = 0; k < 6; k++) {
			AniMetric[0 * 6 + k] = AniSol[p3][k];
			AniMetric[1 * 6 + k] = AniSol[p4][k];
			AniMetric[2 * 6 + k] = AniSol[p1][k];
		}

		double q3 = caltri33_ani(Nodes[p3].pt, Nodes[p4].pt, Nodes[p1].pt, AniMetric);

		for (int k = 0; k < 6; k++) {
			AniMetric[0 * 6 + k] = AniSol[p3][k];
			AniMetric[1 * 6 + k] = AniSol[p4][k];
			AniMetric[2 * 6 + k] = AniSol[p2][k];
		}

		double q4 = caltri33_ani(Nodes[p3].pt, Nodes[p4].pt, Nodes[p2].pt, AniMetric);

		if (std::min(q3, q4) <= std::min(q1, q2) + 1e-3)
			return 0;
	}

	return 1;
}

int DT::checkFlipNormal(int p1, int p2, int p3, int p4) {
    using Vec = Eigen::Vector3d;
    const int ids[4] = {p1, p2, p3, p4};
    for (int id : ids)
        if (id < 0 || id >= static_cast<int>(Nodes.size())) return 0;
    const Vec origin(Nodes[p1].pt[0], Nodes[p1].pt[1], Nodes[p1].pt[2]);
    std::array<Vec, 4> points;
    double scale = 0;
    for (int i = 0; i < 4; ++i) {
        points[i] = Vec(Nodes[ids[i]].pt[0], Nodes[ids[i]].pt[1], Nodes[ids[i]].pt[2]) - origin;
        if (!points[i].allFinite()) return 0;
        scale = std::max(scale, points[i].cwiseAbs().maxCoeff());
    }
    if (!(scale > 0)) return 0;
    for (auto& point : points) point /= scale;

    // Relative area detects unreliable normals independently of mesh size.
    auto normal = [&](int a, int b, int c, Vec& n, double& length) {
        const Vec ab = points[b] - points[a], ac = points[c] - points[a];
        const double edgeSquared = std::max(ab.squaredNorm(),
            std::max(ac.squaredNorm(), (points[c] - points[b]).squaredNorm()));
        n = ab.cross(ac);
        length = n.norm();
        return length > 1e-12 * edgeSquared;
    };

    Vec old1, old2, next1, next2;
    double oldLength1, oldLength2, nextLength1, nextLength2;
    const bool old1Valid = normal(0, 1, 2, old1, oldLength1);
    const bool old2Valid = normal(1, 0, 3, old2, oldLength2);
    // Always validate the replacement, even if the old patch is degenerate.
    if (!normal(2, 3, 0, next1, nextLength1) || !normal(3, 2, 1, next2, nextLength2)) return 0;
    next1 /= nextLength1;
    next2 /= nextLength2;
    const double newCos = next1.dot(next2);
    if (!old1Valid || !old2Valid) {
        // Consistent normals exclude overlapping/oppositely oriented faces.
        return newCos >= std::cos(ANGLE2RADIO(5.0));
    }

    const double minFaceCos = std::cos(ANGLE2RADIO(30.0));
    const Vec oldPatch = old1 + old2;
    old1 /= oldLength1;
    old2 /= oldLength2;
    if (old1.dot(old2) < minFaceCos || newCos < minFaceCos) return 0;
    // The reversed candidate has exactly the opposite normals; no need to
    // recompute both triangles to choose the common orientation.
    if ((next1 + next2).dot(oldPatch) < 0) {
        next1 = -next1;
        next2 = -next2;
    }
    return std::max(next1.dot(old1), next1.dot(old2)) >= minFaceCos
        && std::max(next2.dot(old1), next2.dot(old2)) >= minFaceCos;
}
//waiting TODO:virtual tet seting
int DT::flipEdg(int index , int deep) {
	if (deep > 1)
		return 0;
	int p1 = SurEdgs[index].iStart;
	int p2 = SurEdgs[index].iEnd;
	int f1 = SurEdgs[index].face[0];
	int f2 = SurEdgs[index].face[1];
	int p3, p4, p2_f1_idx, p1_f2_idx;

	for (int i = 0; i < 3; i++) {
		if (SurTris[f1].form[i] != p1 && SurTris[f1].form[i] != p2) {
			p3 = SurTris[f1].form[i];
		}
		else if (SurTris[f1].form[i] == p2) {
			p2_f1_idx = i;
		}
	}
	for (int i = 0; i < 3; i++) {
		if (SurTris[f2].form[i] != p1 && SurTris[f2].form[i] != p2) {
			p4 = SurTris[f2].form[i];
		}
		else if (SurTris[f2].form[i] == p1) {
			p1_f2_idx = i;
		}
	}

	//if ((p1 == 93872 || p1 == 14013) && (p2 == 93872 || p2 == 14013))
	//{
	//	printf("%d %d %d %d %d\n", index, p1, p2, p3, p4);
	//}

	int ia = -1, ib = -1, ic = -1, id = -1;
	int tet = -1, i1 = -1, i2 = -1;
	if (isMeshEdge(p1, p2, &tet)) {
		for (int i = 0; i < 4; i++) {
			if (Elems[tet].form[i] == p1) {
				i1 = i;
			}
			else if (Elems[tet].form[i] == p2) {
				i2 = i;
			}
		}
	}
	else {
		return 0;
	}

	std::vector<int> shell, shell_p, shell1, shell_p1;

	findShell(tet, i1, i2, shell, shell_p);

	std::map<int, int> oldgeoid;
	double minVolume_flip = DBL_MAX;
	for (int i = 0; i < shell.size(); i++) {
		oldgeoid[Elems[shell[i]].geo] = 1;
		if (!isvirtualtet(shell[i]) && !ishulltet(shell[i])) {
			double tpvol = calVolume(shell[i]);
			minVolume_flip = std::min(minVolume_flip, tpvol);
		}
	}
	minVolume_flip /= 10.0;
	dt::Mesh mesh;

	int gp = -1;
	double gpnt[3] = { 0,0,0 };
	std::unordered_map<int, int> nodemp;
	std::unordered_map<int, int> mpn;
	for (int i = 0; i < shell.size(); ++i) {
		for (int j = 0; j < 4; j++) {
			int ip = Elems[shell[i]].form[j];
			if (nodemp.count(ip))
				continue;
			if (ip != ghost) {
				mesh.V.push_back({ Nodes[ip].pt[0], Nodes[ip].pt[1], Nodes[ip].pt[2] });
				gpnt[0] += Nodes[ip].pt[0];
				gpnt[1] += Nodes[ip].pt[1];
				gpnt[2] += Nodes[ip].pt[2];
			}
			else {
				gp = mesh.V.size();
				mesh.V.push_back({ 0,0,0 });
			}
			nodemp[ip] = mesh.V.size() - 1;
			mpn[mesh.V.size() - 1] = ip;
		}

	}

	//update ghost point place
	//if (gp != -1) {
	//	double mid[3] = { 0,0,0 };
	//	for (int i = 0; i < 3; ++i) {
	//		mid[i] = (Nodes[p1].pt[i] + Nodes[p2].pt[i] + Nodes[p3].pt[i] + Nodes[p4].pt[i]) / 4.0;
	//	}

	//	mesh.V[gp][0] = 2 * mid[0] - gpnt[0] / (nodemp.size() - 1);
	//	mesh.V[gp][1] = 2 * mid[1] - gpnt[1] / (nodemp.size() - 1);
	//	mesh.V[gp][2] = 2 * mid[2] - gpnt[2] / (nodemp.size() - 1);
	//}

	if (gp != -1) {
		double mid[3] = { 0.0, 0.0, 0.0 };
		for (int i = 0; i < 3; ++i) {
			mid[i] = 0.25 * (
				Nodes[p1].pt[i] +
				Nodes[p2].pt[i] +
				Nodes[p3].pt[i] +
				Nodes[p4].pt[i]
				);
		}

		auto sameFace = [](int a, int b, int c, int x, int y, int z) -> bool {
			return (a == x || a == y || a == z) &&
				(b == x || b == y || b == z) &&
				(c == x || c == y || c == z);
			};

		auto normalizeVec = [&](double n[3]) -> bool {
			double l = lenvec(n);
			if (l <= 1e-14 * dist_max * dist_max) {
				return false;
			}
			n[0] /= l;
			n[1] /= l;
			n[2] /= l;
			return true;
			};

		auto getSolidFaceNormal = [&](int a, int b, int c, double n[3]) -> bool {
			for (int ie : shell) {
				if (ie < 0) continue;
				if (isDelEle(ie)) continue;
				if (ishulltet(ie)) continue;
				if (isvirtualtet(ie)) continue;

				for (int lf = 0; lf < 4; ++lf) {
					int ia0 = -1, ib0 = -1, ic0 = -1, id0 = -1;
					DFC(lf, ia0, ib0, ic0, id0);

					int v0 = Elems[ie].form[ib0];
					int v1 = Elems[ie].form[ic0];
					int v2 = Elems[ie].form[id0];

					if (!sameFace(v0, v1, v2, a, b, c)) {
						continue;
					}

					calnormal(v0, v1, v2, n);

					return normalizeVec(n);
				}
			}

			return false;
			};

		double n1[3] = { 0.0, 0.0, 0.0 };
		double n2[3] = { 0.0, 0.0, 0.0 };

		bool ok1 = getSolidFaceNormal(p1, p2, p3, n1);
		bool ok2 = getSolidFaceNormal(p1, p2, p4, n2);

		double dir[3] = { 0.0, 0.0, 0.0 };

		if (ok1 && ok2) {
			dir[0] = -(n1[0] + n2[0]);
			dir[1] = -(n1[1] + n2[1]);
			dir[2] = -(n1[2] + n2[2]);

			double ldir = lenvec(dir);

			// 如果两个法向几乎抵消，退化为只用其中一个面法向的反向。
			if (ldir <= 1e-14) {
				dir[0] = -n1[0];
				dir[1] = -n1[1];
				dir[2] = -n1[2];
				ldir = lenvec(dir);
			}

			if (ldir <= 1e-14) {
				// 极端失败情况：不要制造不可控 ghost 点。
				mesh.V[gp][0] = mid[0];
				mesh.V[gp][1] = mid[1];
				mesh.V[gp][2] = mid[2];
			}
			else {
				dir[0] /= ldir;
				dir[1] /= ldir;
				dir[2] /= ldir;

				double d12 = distance(Nodes[p1].pt, Nodes[p2].pt);
				double d34 = distance(Nodes[p3].pt, Nodes[p4].pt);
				double offset = 0.1 * std::min(d12, d34);

				mesh.V[gp][0] = mid[0] + offset * dir[0];
				mesh.V[gp][1] = mid[1] + offset * dir[1];
				mesh.V[gp][2] = mid[2] + offset * dir[2];
			}
		}
		else {
			return 0;
		}
	}


	for (int i = 0; i < shell.size(); i++) {
		for (int j = 0; j < 4; j++) {
			int neig = getNeig(shell[i], j);
			if (std::find(shell.begin(), shell.end(), neig) == shell.end()) {
				DNC(j, ia, ib, ic, id);
				int pa = Elems[shell[i]].form[ib];
				int pb = Elems[shell[i]].form[ic];
				int pc = Elems[shell[i]].form[id];
				mesh.F.push_back({ nodemp[pa],nodemp[pb],nodemp[pc] });
			}

		}
	}
	mesh.F.push_back({ nodemp[p3],nodemp[p4],nodemp[p1] });
	mesh.F.push_back({ nodemp[p3],nodemp[p4],nodemp[p2] });

	//printSph_VTK(shell, "./shell_" + std::to_string(index) + ".vtk");
	//if (index == 242504) {
	//	printSph_VTK(shell,"./shell.vtk");
	//	//// Explicit mesh export is left to the caller.
	//	checkMeshError();
	//	std::string load = "./"+std::to_string(index) + "_out.vtk";
	//	writeVTK(load, mesh, true);
	//}

	dt::Args tempargs;
	tempargs.infolevel = 0;
	tempargs.refine = 0;
	tempargs.optlevel = 1;
	tempargs.optloop = 1;
	tempargs.nthread = num_threads;
	tempargs.constrain = 1;
	tempargs.outlogfile = 0;

	dt::DT d;
    d.parallelMinPointsPerThread = parallelMinPointsPerThread;
    // This is a speculative local mesh. A rejected trial is a recoverable
    // optimization detail, not an error of the caller's mesh generation.
    d.meshLogger = std::make_shared<spdlog::logger>("topology_trial",
        std::make_shared<spdlog::sinks::null_sink_mt>());
    try {
        if (!d.tetrahedralize(mesh, tempargs)) {
            meshLogger->debug("Local topology trial rejected");
            return 0;
        }
    } catch (...) {
        meshLogger->debug("Local topology trial rejected");
        return 0;
    }

	if (mesh.V.size() > nodemp.size()/* || shell.size()< mesh.T.size()*/)
		return 0;

	std::set<int> color;
	for (int i = 0; i < mesh.T.size(); i++) {
		color.insert(mesh.T[i][4]);
	}
	if (color.size() != 2 && mesh.T.size() != 2)
		return 0;

	for (int i = 0; i < mesh.T.size(); i++) {
		if (mpn[mesh.T[i][0]] != ghost && mpn[mesh.T[i][1]] != ghost && mpn[mesh.T[i][2]] != ghost && mpn[mesh.T[i][3]] != ghost) {
			double pa[3] = { mesh.V[mesh.T[i][0]][0],mesh.V[mesh.T[i][0]][1],mesh.V[mesh.T[i][0]][2] };
			double pb[3] = { mesh.V[mesh.T[i][1]][0],mesh.V[mesh.T[i][1]][1],mesh.V[mesh.T[i][1]][2] };
			double pc[3] = { mesh.V[mesh.T[i][2]][0],mesh.V[mesh.T[i][2]][1],mesh.V[mesh.T[i][2]][2] };
			double pd[3] = { mesh.V[mesh.T[i][3]][0],mesh.V[mesh.T[i][3]][1],mesh.V[mesh.T[i][3]][2] };
			double ori = dt::GEOM_FUNC::orient3d(pa, pb, pd, pc);
			if (ori <= minVolume_flip)
				return 0;
		}
	}

	struct Int3 {
		int x, y, z;
		Int3(int a, int b, int c) {
			int vals[3] = { a, b, c };
			std::sort(vals, vals + 3);
			x = vals[0]; y = vals[1]; z = vals[2];
		}
		bool operator==(const Int3& other) const {
			return x == other.x && y == other.y && z == other.z;
		}
	};

	struct Int3Hasher {
		std::size_t operator()(const Int3& k) const {
			std::size_t h1 = std::hash<int>()(k.x);
			std::size_t h2 = std::hash<int>()(k.y);
			std::size_t h3 = std::hash<int>()(k.z);
			return ((h1 ^ (h2 << 1)) >> 1) ^ (h3 << 1);
		}
	};

	// --- 创建临时映射表 ---
	std::unordered_map<Int3, int64_t, Int3Hasher> Tri;
	std::unordered_map<Int3, int, Int3Hasher> checkTri;

	for (int i = 0; i < shell.size(); i++) {
		for (int j = 0; j < 4; j++) {
			int neig = getNeig(shell[i], j);
			int neigord = getNeigOrd(shell[i], j);
			if (std::find(shell.begin(), shell.end(), neig) == shell.end()) {
				DNC(j, ia, ib, ic, id);
				int pa = Elems[shell[i]].form[ib];
				int pb = Elems[shell[i]].form[ic];
				int pc = Elems[shell[i]].form[id];
				Tri[Int3(pa, pb, pc)] = ((int64_t)neig << 2) | neigord;
				checkTri[Int3(pa, pb, pc)] = 1;
			}
		}
	}

	//check if dt is right
	for (int i = 0; i < mesh.T.size(); i++) {
		for (int j = 0; j < 4; j++) {
			DNC(j, ia, ib, ic, id);
			int pa = mpn[mesh.T[i][ib]];
			int pb = mpn[mesh.T[i][ic]];
			int pc = mpn[mesh.T[i][id]];
			Int3 query_key(pa, pb, pc);
			auto it = checkTri.find(query_key);
			if (it != checkTri.end()) {
				it->second++;
			}
			else {
				checkTri[Int3(pa, pb, pc)] = 1;
			}
		}
	}
	for (auto it : checkTri) {
		if (it.second != 2)
			return 0;
	}

	if (deep == 0 && periodic_P.size() != 0) {
		if (periodic_P.find(p1) != periodic_P.end() && periodic_P.find(p2) != periodic_P.end()) {
			for (auto it1 : periodic_P[p1]) {
				for (auto it2 : periodic_P[p2]) {
					const auto* boundaryEntry = BndEdg.find(it1, it2);
					if (boundaryEntry && isParallel(Nodes[p1].pt, Nodes[p2].pt, Nodes[it1].pt, Nodes[it2].pt)) {
						const int boundaryIndex = *boundaryEntry;
						int eid = boundaryIndex;
						if (eid == index)
							continue;
						if (flipEdg(eid,1)) {
							break;
						}
						else {
							return 0;
						}
					}
				}
			}
		}
	}

	//if it will be correct
	for (int i = 0; i < shell.size(); i++) {
		DelEle(shell[i]);
	}

	std::queue<int> newEvec;
	std::vector<int> newEvec2;

	for (int i = 0; i < mesh.T.size(); i++) {
		int forghost = 0;
		for (; forghost < 3; forghost++) {
			if (mpn[mesh.T[i][forghost]] == ghost) {
				break;
			}
		}
		DFC(forghost, ia, ib, ic, id);
		int	newE = addElem(mpn[mesh.T[i][ib]], mpn[mesh.T[i][ic]], mpn[mesh.T[i][id]], mpn[mesh.T[i][ia]]);

		newEvec.push(newE);
		newEvec2.push_back(newE);
		for (int j = 0; j < 4; j++) {
			setP2T(Elems[newE].form[j], newE);
			DNC(j, ia, ib, ic, id);
			int pa = Elems[newE].form[ib];
			int pb = Elems[newE].form[ic];
			int pc = Elems[newE].form[id];
			Int3 query_key(pa, pb, pc);
			auto it = Tri.find(query_key);
			if (it != Tri.end()) {
				int64_t neiginfo = it->second;
				int neig = neiginfo >> 2;
				int neigo = neiginfo & 3;
				bond(newE, j, neig, neigo);
			}
			else {
				Tri[Int3(pa, pb, pc)] = ((int64_t)newE << 2) | j;
			}
		}
	}

	//update f1,f2
	BndTri.erase(SurTris[f1].form[0], SurTris[f1].form[1], SurTris[f1].form[2]);
	BndTri.erase(SurTris[f2].form[0], SurTris[f2].form[1], SurTris[f2].form[2]);


	int newf1 = f1;
	int newf2 = f2;
	auto faceEntry1 = BndTri.try_emplace(p3, p4, p1, newf1);
	if (!faceEntry1.second) {
		const int boundaryIndex = *faceEntry1.first;
		setDelSurTri(f1);
		newf1 = boundaryIndex;
	}
	else {
		SurTris[newf1].form[p2_f1_idx] = p4;
	}


	auto faceEntry2 = BndTri.try_emplace(p3, p4, p2, newf2);
	if (!faceEntry2.second) {
		const int boundaryIndex = *faceEntry2.first;
		setDelSurTri(f2);
		newf2 = boundaryIndex;
	}
	else
	{
		SurTris[newf2].form[p1_f2_idx] = p3;
	}

	//update e0,e1,e2
	BndEdg.erase(p1, p2);
	auto edgeEntry = BndEdg.try_emplace(p3, p4, index);
	if (!edgeEntry.second) {
		const int boundaryIndex = *edgeEntry.first;
		setDelSurEdg(index);
		index = boundaryIndex;
		SurEdgs[index].face.push_back(newf1);
		SurEdgs[index].face.push_back(newf2);
	}
	else {
		SurEdgs[index].iStart = p3;
		SurEdgs[index].iEnd = p4;
	}
	SurEdgs[index].info = 0;


	int neigEdg[4] = { -1 };
	neigEdg[0] = BndEdg.get(p1, p4);
	for (int i = 0; i < SurEdgs[neigEdg[0]].face.size(); i++) {
		if (SurEdgs[neigEdg[0]].face[i] == f2) {
			SurEdgs[neigEdg[0]].face[i] = newf1;
			break;
		}
	}

	neigEdg[1] = BndEdg.get(p2, p3);
	for (int i = 0; i < SurEdgs[neigEdg[1]].face.size(); i++) {
		if (SurEdgs[neigEdg[1]].face[i] == f1) {
			SurEdgs[neigEdg[1]].face[i] = newf2;
			break;
		}
	}
	SurTris[f1].info = 0;
	SurTris[f2].info = 0;

	neigEdg[2] = BndEdg.get(p1, p3);
	for (int i = 0; i < SurEdgs[neigEdg[2]].face.size(); i++) {
		if (SurEdgs[neigEdg[2]].face[i] == f1) {
			SurEdgs[neigEdg[2]].face[i] = newf1;
			break;
		}
	}
	neigEdg[3] = BndEdg.get(p2, p4);
	for (int i = 0; i < SurEdgs[neigEdg[3]].face.size(); i++) {
		if (SurEdgs[neigEdg[3]].face[i] == f2) {
			SurEdgs[neigEdg[3]].face[i] = newf2;
			break;
		}
	}

	//for geoid
	int pushfailtime = 0;
	while (!newEvec.empty()) {
		int ie = newEvec.front();
		newEvec.pop();

		if (ishulltet(ie)) {
			setvirtualtet(ie);
			continue;
		}

		int bndnum = 0;
		for (int k = 0; k < 4; k++) {
			DFC(k, ia, ib, ic, id);
			int pb = Elems[ie].form[ib];
			int pc = Elems[ie].form[ic];
			int pd = Elems[ie].form[id];
			if (isBndTri(pb, pc, pd)) {
				bndnum++;
				continue;
			}
			int neig = getNeig(ie, k);

			if (Elems[neig].geo == -1) {
				continue;
			}
			else {
				Elems[ie].geo = Elems[neig].geo;
				oldgeoid[Elems[ie].geo] = 0;
			}
		}
		if (Elems[ie].geo == -1) {
			if (pushfailtime++ > 1000) {
				auto it1 = oldgeoid.begin();
				auto it2 = std::next(it1);
				if (it1->second == 0) Elems[ie].geo = it2->first;
				if (it2->second == 0) Elems[ie].geo = it1->first;
			}
			else if (bndnum == 4) {
				for (int k = 0; k < 4; k++) {
					int neig = getNeig(ie, k);
					if (Elems[neig].geo > 0) {
						Elems[ie].geo = Elems[neig].geo;
						break;
					}
				}
			}
			else newEvec.push(ie);

		}
	}

	//for (int i = 0; i < 4; i++) {
	//	int ie = neigEdg[i];
	//	if (isDelSurEdg(ie) || SurEdgs[ie].info > 1 || SurEdgs[ie].constrain > 0)
	//		continue;
	//	if (ifflipEdg(ie)) {
	//		flipEdg(ie, deep+1);
	//	}
	//}

	return 1;
}

int DT::smoothBndPntPass(int loop) {
	for (int i = 0; i < loop; i++)
		for (int j = 0; j < Nodes.size(); j++)
			smoothBndPnt(j);
	return 0;
}

// constrain boundary vertex
int DT::smoothBndPnt(int iNod) {
	if (improve_Metric == SUS_METRIC) return 0; // SUS keeps boundary nodes fixed.
	if (periodic_P.size() != 0) {
		// Temporarily not smooth
		if (periodic_P.find(iNod) != periodic_P.end())
			return 0;
	}
	if (improve_Metric == 5) {
		// Anisotropic
		return smoothBndPnt_ani(iNod);
	}
	if (!isbndpnt(iNod)) {
		// if it's not boundary vertex
		return 0;
	}
	else if (isCornerpnt(iNod)) {
		// Don't smooth corner vertex
		return 0;
	}
	else if (isSegmentpnt(iNod)) {
		// Smooth Segment vertex
		// find SurEdgs connect to iNod
		int i, m, k, iElem, iSecond;
		double newpos[3], pernewpos[3], avg[3] = { 0 }, d[3], alpha = 0.1, minq = DBL_MAX;
		double* verts[4];

		std::vector<int> sph, persph;
		std::vector<double> q, perq;
		std::unordered_map<int, bool> seen;

		// ------------------------------------------------------------
		// Periodic paired node
		// ------------------------------------------------------------
		int periNod = -1;
		if (periodic_P.size() != 0) {
			if (periodic_P.find(iNod) != periodic_P.end()) {
				periNod = periodic_P[iNod][0];
			}
		}
		if (periNod == iNod) {
			periNod = -1;
		}

		const int thread_n = omp_get_thread_num();
		if (!tryOccupying(Nodes[iNod].occupying, thread_n)) {
			return -1;
		}
		std::vector<int> vecn;
		vecn.push_back(iNod);

		findSphere(iNod, sph);
		int n = sph.size();

		if (periNod != -1) {
			findSphere(periNod, persph);
		}

		int pern = persph.size();

		if (n <= 1) {
			for (auto it : vecn)
				clearOccupying(Nodes[it].occupying, thread_n);

			return 0;
		}

		if (periNod != -1 && pern <= 1) {
			for (auto it : vecn)
				clearOccupying(Nodes[it].occupying, thread_n);

			return 0;
		}

		// ------------------------------------------------------------
		// Get current worst quality from both cavities
		// ------------------------------------------------------------
		q.resize(n);

		if (periNod != -1) {
			perq.resize(pern);
		}

		for (i = 0; i < n; i++) {
			if (Elems[sph[i]].q >= 0) {
				minq = std::min(minq, Elems[sph[i]].q);
			}
		}

		if (periNod != -1) {
			for (i = 0; i < pern; i++) {
				if (Elems[persph[i]].q >= 0) {
					minq = std::min(minq, Elems[persph[i]].q);
				}
			}
		}

		if (minq == DBL_MAX) {
			for (auto it : vecn)
				clearOccupying(Nodes[it].occupying, thread_n);

			return 0;
		}

		// ------------------------------------------------------------
		// Smooth by Laplace on the original segment node
		// Occupying logic is kept exactly in the original style.
		// ------------------------------------------------------------
		for (i = 0; i < n; i++)
		{
			iElem = sph[i];
			for (m = 0; m < 4; m++)
			{
				iSecond = Elems[iElem].form[m];
				vecn.push_back(iSecond);

				if (!tryOccupying(Nodes[iSecond].occupying, thread_n)) {
					for (auto it : vecn)
						clearOccupying(Nodes[it].occupying, thread_n);
					return -1;
				}

				if (iSecond != iNod)
				{
					if (isCornerpnt(iSecond) || isSegmentpnt(iSecond))
					{
if (auto* boundaryEntry = BndEdg.find(iNod, iSecond))
						{
							const int boundaryIndex = *boundaryEntry;
							int Edgidx = boundaryIndex;
							if (SurEdgs[Edgidx].constrain > 0) {
								if (seen.find(iSecond) == seen.end())
								{
									seen[iSecond] = true;
									avg[0] += Nodes[iSecond].pt[0];
									avg[1] += Nodes[iSecond].pt[1];
									avg[2] += Nodes[iSecond].pt[2];
								}
							}
						}
					}
				}
			}
		}

		// If Error SurEdg FOUND
		if (seen.size() != 2) {
			for (auto it : vecn)
				clearOccupying(Nodes[it].occupying, thread_n);
			return 0;
		}

		int p1 = seen.begin()->first;
		int p2 = std::next(seen.begin())->first;

		if (!P_in_Line(Nodes[iNod].pt, Nodes[p1].pt, Nodes[p2].pt))
		{
			for (auto it : vecn)
				clearOccupying(Nodes[it].occupying, thread_n);
			return 0;
		}

		// ------------------------------------------------------------
		// Optional periodic-side segment consistency check
		// Does not touch Occupying.
		// ------------------------------------------------------------
		if (periNod != -1) {
			std::unordered_map<int, bool> perSeen;

			for (i = 0; i < pern; i++) {
				iElem = persph[i];

				for (m = 0; m < 4; m++) {
					iSecond = Elems[iElem].form[m];

					if (iSecond != periNod)
					{
						if (isCornerpnt(iSecond) || isSegmentpnt(iSecond))
						{
if (auto* boundaryEntry = BndEdg.find(periNod, iSecond))
							{
								const int boundaryIndex = *boundaryEntry;
								int Edgidx = boundaryIndex;
								if (SurEdgs[Edgidx].constrain > 0) {
									if (perSeen.find(iSecond) == perSeen.end())
									{
										perSeen[iSecond] = true;
									}
								}
							}
						}
					}
				}
			}

			if (perSeen.size() != 2) {
				for (auto it : vecn)
					clearOccupying(Nodes[it].occupying, thread_n);
				return 0;
			}

			int pp1 = perSeen.begin()->first;
			int pp2 = std::next(perSeen.begin())->first;

			if (!P_in_Line(Nodes[periNod].pt, Nodes[pp1].pt, Nodes[pp2].pt))
			{
				for (auto it : vecn)
					clearOccupying(Nodes[it].occupying, thread_n);
				return 0;
			}
		}

		// ------------------------------------------------------------
		// Candidate movement
		// ------------------------------------------------------------
		for (int m = 0; m < 3; m++)
		{
			avg[m] /= seen.size();							// center of gravity
			d[m] = avg[m] - Nodes[iNod].pt[m];				// displacement
			newpos[m] = Nodes[iNod].pt[m] + alpha * d[m];	// new position

			if (periNod != -1) {
				// This assumes translational periodicity.
				// If periodicity is rotational or reflective,
				// d should be transformed before being applied.
				pernewpos[m] = Nodes[periNod].pt[m] + alpha * d[m];
			}
		}

		// ------------------------------------------------------------
		// Quality-protected line search
		// ------------------------------------------------------------
		bool moveflag = true;
		int iter = 0;
		while (iter < 6) {
			moveflag = true;

			// Check original cavity
			for (i = 0; i < n; i++)
			{
				if (Elems[sph[i]].q < 0)
				{
					continue;
				}

				for (int m = 0; m <= 3; m++)
				{
					int iElemNd = Elems[sph[i]].form[m];
					verts[m] = (iElemNd != iNod) ? Nodes[iElemNd].pt : newpos;
				}

				q[i] = tetquality(
					verts[0], verts[1], verts[2], verts[3],
					NULL,
					improve_Metric
				);

				if (q[i] <= minq) {
					moveflag = false;
					break;
				}
			}

			// Check periodic paired cavity
			if (moveflag && periNod != -1) {
				for (i = 0; i < pern; i++)
				{
					if (Elems[persph[i]].q < 0)
					{
						continue;
					}

					for (int m = 0; m <= 3; m++)
					{
						int iElemNd = Elems[persph[i]].form[m];
						verts[m] = (iElemNd != periNod) ? Nodes[iElemNd].pt : pernewpos;
					}

					perq[i] = tetquality(
						verts[0], verts[1], verts[2], verts[3],
						NULL,
						improve_Metric
					);

					if (perq[i] <= minq) {
						moveflag = false;
						break;
					}
				}
			}

			if (moveflag) {
				break;
			}
			else {
				alpha /= 2.0;

				for (int j = 0; j < 3; j++) {
					newpos[j] = Nodes[iNod].pt[j] + alpha * d[j];

					if (periNod != -1) {
						pernewpos[j] = Nodes[periNod].pt[j] + alpha * d[j];
					}
				}

				iter++;
			}
		}

		if (moveflag) {
			// ------------------------------------------------------------
			// Update qualities of original cavity
			// ------------------------------------------------------------
			for (i = 0; i < n; i++) {
				if (Elems[sph[i]].q != -1)
					Elems[sph[i]].q = q[i];
			}

			// ------------------------------------------------------------
			// Update qualities of periodic paired cavity
			// Do not assume persph.size() == sph.size().
			// ------------------------------------------------------------
			if (periNod != -1) {
				for (i = 0; i < pern; i++) {
					if (Elems[persph[i]].q != -1)
						Elems[persph[i]].q = perq[i];
				}
			}

			// ------------------------------------------------------------
			// Update point positions
			// ------------------------------------------------------------
			for (int j = 0; j < 3; j++)
				Nodes[iNod].pt[j] = newpos[j];

			if (periNod != -1) {
				for (int j = 0; j < 3; j++)
					Nodes[periNod].pt[j] = pernewpos[j];
			}

			for (auto it : vecn)
				clearOccupying(Nodes[it].occupying, thread_n);

			return 1;
		}
		else {
			for (auto it : vecn)
				clearOccupying(Nodes[it].occupying, thread_n);

			return 0;
		}
	}
else if (isFacetpnt(iNod)) {
	// Smooth Facet vertex
	int i, m, k, iElem, iSecond;
	double newpos[3], pernewpos[3], avg[3] = { 0 }, d[3], alpha = 0.1, minq = DBL_MAX;
	double* verts[4];

	std::vector<int> sph, persph;
	std::vector<double> q, perq;
	std::unordered_map<int, bool> seen;

	// ------------------------------------------------------------
	// Periodic paired node
	// ------------------------------------------------------------
	int periNod = -1;
	if (periodic_P.size() != 0) {
		if (periodic_P.find(iNod) != periodic_P.end()) {
			periNod = periodic_P[iNod][0];
		}
	}
	if (periNod == iNod) {
		periNod = -1;
	}

	const int thread_n = omp_get_thread_num();
	if (!tryOccupying(Nodes[iNod].occupying, thread_n)) {
		return -1;
	}

	std::vector<int> vecn;
	vecn.push_back(iNod);

	findSphere(iNod, sph);
	int n = sph.size();

	if (periNod != -1) {
		findSphere(periNod, persph);
	}
	int pern = persph.size();

	if (n <= 1) {
		for (auto it : vecn)
			clearOccupying(Nodes[it].occupying, thread_n);
		return 0;
	}

	if (periNod != -1 && pern <= 1) {
		for (auto it : vecn)
			clearOccupying(Nodes[it].occupying, thread_n);
		return 0;
	}

	// ------------------------------------------------------------
	// Get current worst quality from both cavities
	// ------------------------------------------------------------
	q.resize(n);

	if (periNod != -1) {
		perq.resize(pern);
	}

	for (i = 0; i < n; i++) {
		if (Elems[sph[i]].q >= 0) {
			minq = std::min(minq, Elems[sph[i]].q);
		}
	}

	if (periNod != -1) {
		for (i = 0; i < pern; i++) {
			if (Elems[persph[i]].q >= 0) {
				minq = std::min(minq, Elems[persph[i]].q);
			}
		}
	}

	if (minq == DBL_MAX) {
		for (auto it : vecn)
			clearOccupying(Nodes[it].occupying, thread_n);
		return 0;
	}

	// ------------------------------------------------------------
	// Smooth by Laplace on the original facet node
	// Occupying logic is kept in the original style.
	// ------------------------------------------------------------
	for (i = 0; i < n; i++)
	{
		iElem = sph[i];
		for (m = 0; m < 4; m++)
		{
			iSecond = Elems[iElem].form[m];
			vecn.push_back(iSecond);

			if (!tryOccupying(Nodes[iSecond].occupying, thread_n)) {
				for (auto it : vecn)
					clearOccupying(Nodes[it].occupying, thread_n);
				return -1;
			}

			if (iSecond != iNod)
			{
				if (isbndpnt(iSecond)) {
					if (isBndEdg(iNod, iSecond)) {
						if (seen.find(iSecond) == seen.end())
						{
							seen[iSecond] = true;
							avg[0] += Nodes[iSecond].pt[0];
							avg[1] += Nodes[iSecond].pt[1];
							avg[2] += Nodes[iSecond].pt[2];
						}
					}
				}
			}
		}
	}

	// Avoid division by zero.
	if (seen.size() == 0) {
		for (auto it : vecn)
			clearOccupying(Nodes[it].occupying, thread_n);
		return 0;
	}

	// ------------------------------------------------------------
	// Candidate movement
	// ------------------------------------------------------------
	for (int m = 0; m < 3; m++)
	{
		avg[m] /= seen.size();							// center of gravity
		d[m] = avg[m] - Nodes[iNod].pt[m];				// displacement
		newpos[m] = Nodes[iNod].pt[m] + alpha * d[m];	// new position

		if (periNod != -1) {
			// This assumes translational periodicity.
			// If periodicity is rotational or reflective,
			// d should be transformed before being applied.
			pernewpos[m] = Nodes[periNod].pt[m] + alpha * d[m];
		}
	}

	// ------------------------------------------------------------
	// Quality-protected line search
	// ------------------------------------------------------------
	bool moveflag = true;
	int iter = 0;

	while (iter < 6) {
		moveflag = true;

		// Check original cavity
		for (i = 0; i < n; i++)
		{
			if (Elems[sph[i]].q < 0)
			{
				continue;
			}

			for (int m = 0; m <= 3; m++)
			{
				int iElemNd = Elems[sph[i]].form[m];
				verts[m] = (iElemNd != iNod) ? Nodes[iElemNd].pt : newpos;
			}

			q[i] = tetquality(
				verts[0], verts[1], verts[2], verts[3],
				NULL,
				improve_Metric
			);

			if (q[i] <= minq) {
				moveflag = false;
				break; // This tet becomes invalid.
			}
		}

		// Check periodic paired cavity
		if (moveflag && periNod != -1) {
			for (i = 0; i < pern; i++)
			{
				if (Elems[persph[i]].q < 0)
				{
					continue;
				}

				for (int m = 0; m <= 3; m++)
				{
					int iElemNd = Elems[persph[i]].form[m];
					verts[m] = (iElemNd != periNod) ? Nodes[iElemNd].pt : pernewpos;
				}

				perq[i] = tetquality(
					verts[0], verts[1], verts[2], verts[3],
					NULL,
					improve_Metric
				);

				if (perq[i] <= minq) {
					moveflag = false;
					break;
				}
			}
		}

		if (moveflag) {
			break;
		}
		else {
			alpha /= 2.0;

			for (int j = 0; j < 3; j++) {
				newpos[j] = Nodes[iNod].pt[j] + alpha * d[j];

				if (periNod != -1) {
					pernewpos[j] = Nodes[periNod].pt[j] + alpha * d[j];
				}
			}

			iter++;
		}
	}

	if (moveflag) {
		// ------------------------------------------------------------
		// Update qualities of original cavity
		// ------------------------------------------------------------
		for (i = 0; i < n; i++) {
			if (Elems[sph[i]].q != -1)
				Elems[sph[i]].q = q[i];
		}

		// ------------------------------------------------------------
		// Update qualities of periodic paired cavity
		// Do not assume persph.size() == sph.size().
		// ------------------------------------------------------------
		if (periNod != -1) {
			for (i = 0; i < pern; i++) {
				if (Elems[persph[i]].q != -1)
					Elems[persph[i]].q = perq[i];
			}
		}

		// ------------------------------------------------------------
		// Update point positions
		// ------------------------------------------------------------
		for (int j = 0; j < 3; j++)
			Nodes[iNod].pt[j] = newpos[j];

		if (periNod != -1) {
			for (int j = 0; j < 3; j++)
				Nodes[periNod].pt[j] = pernewpos[j];
		}

		for (auto it : vecn)
			clearOccupying(Nodes[it].occupying, thread_n);

		return 1;
	}
	else {
		for (auto it : vecn)
			clearOccupying(Nodes[it].occupying, thread_n);

		return 0;
	}
	}
	return 0;
}

// Legacy anisotropic boundary-point smoothing. Kept for reference.
int DT::smoothBndPnt_ani(int iNod) {
	if (!isbndpnt(iNod)) {
		// if it's not boundary vertex
		return 0;
	}
	else if (isCornerpnt(iNod)) {
		// Don't smooth corner vertex
		return 0;
	}
	else if (isSegmentpnt(iNod)) {
		// Smooth Segment vertex
		// find SurEdgs connect to iNod
		int i, m, k, iElem, iSecond;
		double pernewpos[3],newpos[3], d[3], alpha = 1, minq = DBL_MAX;//alpha is smooth Gradient
		double* verts[4];
		std::vector<int> sph,persph;
		std::vector<double> q,perq;
		std::unordered_map<int, bool> seen;

		int periNod = -1;
		if (periodic_P.size() != 0) {
			if (periodic_P.find(iNod) != periodic_P.end()) {
				periNod = periodic_P[iNod][0];
			}
		}

		findSphere(iNod, sph);
		int n = sph.size(),pern=-1;
		//get worse quality
		q.resize(n);
		for (i = 0; i < n; i++) {
			if (Elems[sph[i]].q >= 0) {
				minq = std::min(minq, Elems[sph[i]].q);
			}
		}
		if (periNod != -1) {
			findSphere(periNod, persph);
			pern = persph.size();
			perq.resize(pern);
			for (i = 0; i < pern; i++) {
				if (Elems[persph[i]].q >= 0) {
					minq = std::min(minq, Elems[persph[i]].q);
				}
			}
		}

		//smooth by laplace
		for (i = 0; i < n; i++)
		{
			iElem = sph[i];
			for (m = 0; m < 4; m++)
			{
				iSecond = Elems[iElem].form[m];
				if (iSecond != iNod)
				{
					if (isCornerpnt(iSecond) || isSegmentpnt(iSecond))
					{
if (auto* boundaryEntry = BndEdg.find(iNod, iSecond))
						{
							const int boundaryIndex = *boundaryEntry;
							int Edgidx = boundaryIndex;
							if (SurEdgs[Edgidx].constrain>0) {
								if (seen.find(iSecond) == seen.end())
								{
									seen[iSecond] = true;
								}
							}
						}
					}
				}
			}
		}

		int p1 = -1, p2 = -1;
		for (auto it : seen) {
			if (p1 == -1)
				p1 = it.first;
			else
				p2 = it.first;
		}

		double m1[6] = { 0 }, m2[6] = { 0 }, m0[6] = { 0 };
		for (int mm = 0; mm < 6; mm++) {
			m0[mm] = AniSol[iNod][mm];
			m1[mm] = AniSol[p1][mm];
			m2[mm] = AniSol[p2][mm];
		}

		// 初始二分区间
		double left = 0.0, right = 1.0;
		double tol = 1e-6;
		int max_iter = 30;

		double  pos[3];
		double* A = Nodes[p1].pt;
		double* B = Nodes[p2].pt;
		for (int iter = 0; iter < max_iter; iter++) {
			double mid = 0.5 * (left + right);

			// 插值得到候选点
			for (int i = 0; i < 3; i++)
				pos[i] = (1 - mid) * A[i] + mid * B[i];

			double len1 = cal_ani_length(A, pos, m1, m0);
			double len2 = cal_ani_length(B, pos, m2, m0);

			double diff = len1 - len2;
			if (fabs(diff) < tol)
				break;

			if (diff > 0)
				right = mid;
			else
				left = mid;
		}

		double projected[3] = { pos[0], pos[1], pos[2] };
		double proj_tol = 0.3 * distance(A, B);
		bool segment_proj_ok = project_segment_point_to_fine_mesh(pos, projected, proj_tol);
		if (segment_proj_ok) {
			minq = std::min(minq, 1e-10);
		}

		for (int m = 0; m < 3; m++)
		{				//center of gravity
			d[m] = projected[m] - Nodes[iNod].pt[m];				//displacement
			newpos[m] = Nodes[iNod].pt[m] + alpha * d[m];	//new position
			if (periNod != -1) {
				pernewpos[m] = Nodes[periNod].pt[m] + alpha * d[m];
			}
		}

		int e1 = BndEdg.get(p1, iNod);
		int e2 = BndEdg.get(p2, iNod);
		std::vector<int> chkF;

		for (auto it : SurEdgs[e1].face)
			chkF.push_back(it);
		for (auto it : SurEdgs[e2].face)
			chkF.push_back(it);

		double oldNormal[4][3] = { 0.0 };
		double oldNormalLen[4] = { 0.0 };

		for (int fi = 0; fi < (int)chkF.size(); fi++) {
			int iTri = chkF[fi];

			int p0 = SurTris[iTri].form[0];
			int p1 = SurTris[iTri].form[1];
			int p2 = SurTris[iTri].form[2];

			double* v0 = Nodes[p0].pt;
			double* v1 = Nodes[p1].pt;
			double* v2 = Nodes[p2].pt;

			double e01[3] = {
				v1[0] - v0[0],
				v1[1] - v0[1],
				v1[2] - v0[2]
			};

			double e02[3] = {
				v2[0] - v0[0],
				v2[1] - v0[1],
				v2[2] - v0[2]
			};

			cross(e01, e02, oldNormal[fi]);
			oldNormalLen[fi] = lenvec(oldNormal[fi]);

			if (oldNormalLen[fi] <= 1e-30) return 0;
		}

		bool moveflag = true;
		int iter = 0;
		while (iter < 6) {
			moveflag = true;
			for (i = 0; i < n; i++)
			{
				if (Elems[sph[i]].q < 0)
				{
					continue;
				}
				for (int m = 0; m <= 3; m++)
				{
					int iElemNd = Elems[sph[i]].form[m];
					verts[m] = (iElemNd != iNod) ? Nodes[iElemNd].pt : newpos;
				}

				double AniMetric[6] = { 0 };
				getmm(Elems[sph[i]].form[0], Elems[sph[i]].form[1], Elems[sph[i]].form[2], Elems[sph[i]].form[3], AniMetric);

				q[i] = tetquality(verts[0], verts[1], verts[2], verts[3], AniMetric, improve_Metric);

				if (q[i] <= minq) {
					moveflag = false;
					break; // This tet becomes invalid.
				}
			}
			if(moveflag && periNod != -1) {
				for (i = 0; i < pern; i++)
				{
					if (Elems[persph[i]].q < 0)
					{
						continue;
					}
					for (int m = 0; m <= 3; m++)
					{
						int iElemNd = Elems[persph[i]].form[m];
						verts[m] = (iElemNd != periNod) ? Nodes[iElemNd].pt : pernewpos;
					}

					double AniMetric[6] = { 0 };
					getmm(Elems[persph[i]].form[0], Elems[persph[i]].form[1], Elems[persph[i]].form[2], Elems[persph[i]].form[3], AniMetric);

					perq[i] = tetquality(verts[0], verts[1], verts[2], verts[3], AniMetric, improve_Metric);

					if (perq[i] <= minq) {
						moveflag = false;
						break; // This tet becomes invalid.
					}
				}
			}
			if (moveflag) {
				for (int fi = 0; fi < (int)chkF.size(); fi++) {
					int iTri = chkF[fi];

					int p0 = SurTris[iTri].form[0];
					int p1 = SurTris[iTri].form[1];
					int p2 = SurTris[iTri].form[2];

					double* old0 = Nodes[p0].pt;
					double* old1 = Nodes[p1].pt;
					double* old2 = Nodes[p2].pt;

					double* v0 = (p0 == iNod) ? newpos : old0;
					double* v1 = (p1 == iNod) ? newpos : old1;
					double* v2 = (p2 == iNod) ? newpos : old2;

					double e01[3] = {
						v1[0] - v0[0],
						v1[1] - v0[1],
						v1[2] - v0[2]
					};

					double e02[3] = {
						v2[0] - v0[0],
						v2[1] - v0[1],
						v2[2] - v0[2]
					};

					double newNormal[3] = { 0.0, 0.0, 0.0 };
					cross(e01, e02, newNormal);

					double newNormalLen = lenvec(newNormal);
					if (newNormalLen <= 1e-30) {
						moveflag = false;
						break;
					}

					double cosAngle =
						dot(oldNormal[fi], newNormal) /
						(oldNormalLen[fi] * newNormalLen);

					if (cosAngle > 1.0) cosAngle = 1.0;
					if (cosAngle < -1.0) cosAngle = -1.0;

					if (cosAngle < 0.99) {
						moveflag = false;
						break;
					}
				}
			}
			if (moveflag) {
				break;
			}
			else {
				alpha /= 2.0;
				for (int j = 0; j < 3; j++) {
					newpos[j] = Nodes[iNod].pt[j] + alpha * d[j];
					if (periNod != -1) {
						pernewpos[j] = Nodes[periNod].pt[j] + alpha * d[j];
					}
				}
				iter++;
			}
		} // while (iter < 3)
		if (moveflag) {
			//update qual
			for (i = 0; i < n; i++) {
				if (Elems[sph[i]].q != -1) {
					Elems[sph[i]].q = q[i];
				}
			}

			if (periNod != -1) {
				for (i = 0; i < pern; i++) {
					if (Elems[persph[i]].q != -1) {
						Elems[persph[i]].q = perq[i];
					}
				}
			}

			//update point position
			for (int j = 0; j < 3; j++)
				Nodes[iNod].pt[j] = newpos[j];

			if (periNod != -1) {
				for (int j = 0; j < 3; j++)
					Nodes[periNod].pt[j] = pernewpos[j];
			}
			return 1;
		}
		else {
			return 0;
		}
	}
	else if (isFacetpnt(iNod)) {
		int i, m, a, b, c, d;
		double* verts[4];
		double newpos[3],pernewpos[3], target[3], dis[3];
		double alpha = 1.0;
		double minq = DBL_MAX;

		std::vector<int> sph,persph;
		std::set<int> surF;
		std::set<int> surNod;
		std::vector<double> q,perq;

		findSphere(iNod, sph);
		int n = (int)sph.size(),pern=-1;
		if (n == 0) return 0;

		q.resize(n);

		int periNod = -1;
		if (periodic_P.size() != 0) {
			if (periodic_P.find(iNod) != periodic_P.end()) {
				periNod = periodic_P[iNod][0];
			}
		}
		if (periNod != -1) {
			findSphere(periNod, persph);
			pern = (int)persph.size();
	
			perq.resize(pern);
		}


		// ------------------------------------------------------------
		// 1. Current worst quality in the local cavity
		// ------------------------------------------------------------
		for (i = 0; i < n; i++) {
			if (Elems[sph[i]].q >= 0) {
				minq = std::min(minq, Elems[sph[i]].q);
			}
		}

		if (periNod != -1) {
			for (i = 0; i < pern; i++) {
				if (Elems[persph[i]].q >= 0) {
					minq = std::min(minq, Elems[persph[i]].q);
				}
			}
		}


		// ------------------------------------------------------------
		// 2. Collect incident boundary triangles around iNod
		// ------------------------------------------------------------
		for (i = 0; i < n; i++) {
			int iElem = sph[i];
			if (isDelEle(iElem)) continue;

			for (m = 0; m < 4; m++) {
				if (Elems[iElem].form[m] == iNod) continue;

				DFC(m, a, b, c, d);

				b = Elems[iElem].form[b];
				c = Elems[iElem].form[c];
				d = Elems[iElem].form[d];

				if (auto* boundaryEntry = BndTri.find(b, c, d)) {
					const int boundaryIndex = *boundaryEntry;
					int iTri = boundaryIndex;
					if (iTri >= 0 && !isDelSurTri(iTri)) {
						surF.insert(iTri);
					}
				}
			}
		}

		if (surF.empty()) return 0;

		// ------------------------------------------------------------
		// 3. Collect surface one-ring neighboring nodes
		// ------------------------------------------------------------
		for (auto iTri : surF) {
			int p0 = SurTris[iTri].form[0];
			int p1 = SurTris[iTri].form[1];
			int p2 = SurTris[iTri].form[2];

			if (p0 != iNod) surNod.insert(p0);
			if (p1 != iNod) surNod.insert(p1);
			if (p2 != iNod) surNod.insert(p2);
		}

		if (surNod.size() < 2) return 0;

		// ------------------------------------------------------------
		// 4. Compute anisotropic weighted surface target
		// ------------------------------------------------------------
		target[0] = target[1] = target[2] = 0.0;
		double* p0 = Nodes[iNod].pt;
		double wsum = 0.0;
		double avgLen = 0.0;
		int avgCnt = 0;

		for (auto iTri : surF) {
			int fp0 = SurTris[iTri].form[0];
			int fp1 = SurTris[iTri].form[1];
			int fp2 = SurTris[iTri].form[2];

			double* v0 = Nodes[fp0].pt;
			double* v1 = Nodes[fp1].pt;
			double* v2 = Nodes[fp2].pt;

			double area = calArea(v0, v1, v2);
			if (area <= 1e-30) continue;

			double mm[6] = { 0.0 };
			for (int j = 0; j < 6; j++) {
				mm[j] = (AniSol[fp0][j] + AniSol[fp1][j] + AniSol[fp2][j]) / 3.0;
			}

			double det = mm[0] * (mm[3] * mm[5] - mm[4] * mm[4])
				- mm[1] * (mm[1] * mm[5] - mm[2] * mm[4])
				+ mm[2] * (mm[1] * mm[4] - mm[2] * mm[3]);
			if (det <= 1e-30) continue;

			double w = std::sqrt(det) * area;
			target[0] += w * (v0[0] + v1[0] + v2[0]) / 3.0;
			target[1] += w * (v0[1] + v1[1] + v2[1]) / 3.0;
			target[2] += w * (v0[2] + v1[2] + v2[2]) / 3.0;
			wsum += w;
		}

		for (auto pj : surNod) {
			double elen = distance(p0, Nodes[pj].pt);
			if (elen < 1e-30) continue;
			avgLen += elen;
			avgCnt++;
		}

		if (wsum <= 1e-30 || avgCnt == 0) return 0;

		target[0] /= wsum;
		target[1] /= wsum;
		target[2] /= wsum;
		avgLen /= (double)avgCnt;

		// ------------------------------------------------------------
		// 5. Project the anisotropic target onto the fine boundary mesh
		// ------------------------------------------------------------
		double projected[3] = { target[0], target[1], target[2] };
		project_boundary_point_to_fine_mesh(target, projected, avgLen);

		dis[0] = projected[0] - p0[0];
		dis[1] = projected[1] - p0[1];
		dis[2] = projected[2] - p0[2];

		double dlen = lenvec(dis);
		if (dlen < 1e-30) return 0;

		newpos[0] = p0[0] + alpha * dis[0];
		newpos[1] = p0[1] + alpha * dis[1];
		newpos[2] = p0[2] + alpha * dis[2];

		if (periNod != -1) {
			pernewpos[0] = Nodes[periNod].pt[0] + alpha * dis[0];
			pernewpos[1] = Nodes[periNod].pt[1] + alpha * dis[1];
			pernewpos[2] = Nodes[periNod].pt[2] + alpha * dis[2];
		}

		// ------------------------------------------------------------
		// 6. Quality-protected line search
		// ------------------------------------------------------------
		bool moveflag = true;
		int iter = 0;

		while (iter < 8) {
			moveflag = true;

			// --------------------------------------------------------
			// 7.1 Check boundary surface triangles do not flip
			// --------------------------------------------------------
			for (auto iTri : surF) {
				int p[3] = {
					SurTris[iTri].form[0],
					SurTris[iTri].form[1],
					SurTris[iTri].form[2]
				};

				double* old0 = Nodes[p[0]].pt;
				double* old1 = Nodes[p[1]].pt;
				double* old2 = Nodes[p[2]].pt;

				double* v0 = (p[0] == iNod) ? newpos : old0;
				double* v1 = (p[1] == iNod) ? newpos : old1;
				double* v2 = (p[2] == iNod) ? newpos : old2;

				// 1. Area preservation check
				double oldArea = calArea(old0, old1, old2);
				double newArea = calArea(v0, v1, v2);

				if (oldArea <= 1e-30 || newArea < 0.3 * oldArea) {
					moveflag = false;
					break;
				}

				// 2. Normal direction / flipping check
				double old_e01[3] = {
					old1[0] - old0[0],
					old1[1] - old0[1],
					old1[2] - old0[2]
				};

				double old_e02[3] = {
					old2[0] - old0[0],
					old2[1] - old0[1],
					old2[2] - old0[2]
				};

				double new_e01[3] = {
					v1[0] - v0[0],
					v1[1] - v0[1],
					v1[2] - v0[2]
				};

				double new_e02[3] = {
					v2[0] - v0[0],
					v2[1] - v0[1],
					v2[2] - v0[2]
				};

				double old_n[3] = { 0.0, 0.0, 0.0 };
				double new_n[3] = { 0.0, 0.0, 0.0 };

				cross(old_e01, old_e02, old_n);
				cross(new_e01, new_e02, new_n);

				double old_n_len = lenvec(old_n);
				double new_n_len = lenvec(new_n);

				if (old_n_len <= 1e-30 || new_n_len <= 1e-30) {
					moveflag = false;
					break;
				}

				double cosOldNew = dot(old_n, new_n) / (old_n_len * new_n_len);

				if (cosOldNew > 1.0) cosOldNew = 1.0;
				if (cosOldNew < -1.0) cosOldNew = -1.0;

				// cos(60 degrees) = 0.5.
				// cos (10) =0.98
				// This also prevents flipping because flipped triangles have cosOldNew < 0.
				if (cosOldNew < 0.98) {
					moveflag = false;
					break;
				}
			}

			// --------------------------------------------------------
			// 7.2 Check volume and tet quality
			// --------------------------------------------------------
			for (i = 0; i < n; i++) {
				if (Elems[sph[i]].q < 0) continue;

				for (m = 0; m < 4; m++) {
					int iElemNd = Elems[sph[i]].form[m];
					verts[m] = (iElemNd != iNod) ? Nodes[iElemNd].pt : newpos;
				}

				double AniMetric[6] = { 0.0 };

				getmm(Elems[sph[i]].form[0], Elems[sph[i]].form[1], Elems[sph[i]].form[2], Elems[sph[i]].form[3], AniMetric);

				q[i] = tetquality(verts[0], verts[1], verts[2], verts[3], AniMetric, improve_Metric);

				// Do not degrade the current local worst quality.
				if (q[i] < minq) {
					moveflag = false;
					break;
				}
			}
			if (moveflag && periNod != -1) {
				for (i = 0; i < pern; i++) {
					if (Elems[persph[i]].q < 0) continue;

					for (m = 0; m < 4; m++) {
						int iElemNd = Elems[persph[i]].form[m];
						verts[m] = (iElemNd != periNod) ? Nodes[iElemNd].pt : pernewpos;
					}

					double AniMetric[6] = { 0.0 };

					getmm(Elems[persph[i]].form[0], Elems[persph[i]].form[1], Elems[persph[i]].form[2], Elems[persph[i]].form[3], AniMetric);

					perq[i] = tetquality(verts[0], verts[1], verts[2], verts[3], AniMetric, improve_Metric);

					// Do not degrade the current local worst quality.
					if (perq[i] < minq) {
						moveflag = false;
						break;
					}
				}
			}

			if (moveflag) break;

			alpha *= 0.5;

			newpos[0] = p0[0] + alpha * dis[0];
			newpos[1] = p0[1] + alpha * dis[1];
			newpos[2] = p0[2] + alpha * dis[2];
			if (periNod != -1) {
				pernewpos[0] = Nodes[periNod].pt[0] + alpha * dis[0];
				pernewpos[1] = Nodes[periNod].pt[1] + alpha * dis[1];
				pernewpos[2] = Nodes[periNod].pt[2] + alpha * dis[2];
			}

			iter++;
		}

		if (!moveflag) return 0;

		// ------------------------------------------------------------
		// 8. Require at least one local element to improve
		// ------------------------------------------------------------
		bool improve = false;

		for (i = 0; i < n; i++) {
			if (Elems[sph[i]].q >= 0 && q[i] > Elems[sph[i]].q + 1e-14) {
				improve = true;
				break;
			}
		}

		if (periNod != -1) {
			for (i = 0; i < pern; i++) {
				if (Elems[persph[i]].q >= 0 && perq[i] > Elems[persph[i]].q + 1e-14) {
					improve = true;
					break;
				}
			}

		}

		if (!improve) return 0;

		// ------------------------------------------------------------
		// 9. Update qualities and point position
		// ------------------------------------------------------------
		for (i = 0; i < n; i++) {
			if (Elems[sph[i]].q != -1) {
				Elems[sph[i]].q = q[i];
			}
		}

		Nodes[iNod].pt[0] = newpos[0];
		Nodes[iNod].pt[1] = newpos[1];
		Nodes[iNod].pt[2] = newpos[2];

		if (periNod != -1) {
			for (i = 0; i < pern; i++) {
				if (Elems[persph[i]].q != -1) {
					Elems[persph[i]].q = perq[i];
				}
			}

			Nodes[periNod].pt[0] = pernewpos[0];
			Nodes[periNod].pt[1] = pernewpos[1];
			Nodes[periNod].pt[2] = pernewpos[2];
		}
		return 1;
	}
	return 0;
}

void DT::Smooth_size_ani(Mesh& mesh, std::vector<std::array<double, 6>>& anisol, std::vector<int> lockFactes, std::vector<int> lockVertex) {
	AniSol = anisol;
	return;
}
