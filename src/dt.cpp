#include "dt.h"

using namespace std;

constexpr int DT::virtualID; // C++14 definition for reference-taking callers.

DT::DT() {};
DT::~DT() {};

int DT::tetrahedralize(Mesh& mesh, Args& args)
{
	if (!dt_init(mesh, args))
		return 0;
    MeshStageLog totalLog(*this, "Tetrahedralize");

	if (!BndPntInst(mesh, args)) {
        meshLogger->error("Boundary point insertion failed");
        return 0;
    }

	BoudaryRecover(mesh, args);

	ColorVirtualTet(args);

	MeshRefine(args);


	MeshImprove(args);
	RemoveTet(args);
	outMesh(mesh, args);
	return 1;
}

// Start a new mesh run, retaining allocated storage and independent tuning
// parameters. No worker may be using this same instance while it is reset.
void DT::resetMeshState() {
    fineMeshProjection.enabled = false; fineMeshProjection.edgesEnabled = false;
    Nodes.resize(0); Elems.resize(0);
#ifdef MEMORY_POOL
    qual.resize(0);
#endif
    SurTris.clear(); SurEdgs.clear();
    while (!Evacancy.empty()) Evacancy.pop();
    while (!Nvacancy.empty()) Nvacancy.pop();
    for (auto& queue : Evacancy_thread) while (!queue.empty()) queue.pop();
    SteinerOrd.clear(); EdgSteiner.clear(); TriSteiner.clear(); AniSol.clear();
    BndTri.clear(); BndEdg.clear(); periodic_P.clear(); surTri_mapping.clear();
    lockF.clear(); lockE.clear(); lockV.clear();
    GeoNum.clear(); FacetNum.clear(); QuantityControl = false;
    flipnmRecll.clear();
    for (auto& history : parallelFlipHistory) history.clear();
    parallelFlipCount.fill(0); parallelTopologyBatch = false;
    susIdleStates.clear(); serialBW.slots.clear();
    ignoreIntersect = false; modifyBnd = false; addBoxFlag = false;
    tempfliptime = 0; initBWshell = 0; addst = 0; addstbnd = 0; AvgEdgLen = 0;
    minVolume_bw = 0; // Derived from the mesh/algorithm run.
}

int DT::dt_init(Mesh& mesh, Args& args)
{
	dt::GEOM_FUNC::exactinit();
    resetMeshState();
	ghost = -1;
	infolevel = std::max(0, std::min(2, args.infolevel));
	spdlogoutfile(args.outlogfile);
    meshLogger->set_level(infolevel == 0 ? spdlog::level::err :
        infolevel == 1 ? spdlog::level::info : spdlog::level::debug);
	if (infolevel > 0) meshLogger->info("Version 2026.09.14");
    MeshStageLog initLog(*this, "Initialize");
	improve_step = false;
	cos_collinear_ang_tol = cos(179.9999 / 180. * PI);
	seg[0] = seg[1] = -1;
	fac[0] = fac[1] = fac[2] = -1;
	if (args.ignoreIntersect) {
		ignoreIntersect = true;
		args.refine = 0;
		args.optlevel = 0;
	}
	meshSize = args.size;
	minEdge = args.minEdge <= 0 ? -1 : args.minEdge;
	maxEdge = args.maxEdge <= 0 ? -1 : args.maxEdge;
	if (minEdge != -1 && maxEdge != -1 && maxEdge < minEdge)
		maxEdge = -1;
	growsize = args.growsize < 1 ? 1 : args.growsize;
	improve_Metric = 3;
	maxfliptime = 10000;
	if (args.size <= 0) {
		meshLogger->error("args.size = {}", args.size);
		return 0;
	}
	if (mesh.V.size() == 0) {
		meshLogger->error("Mesh'Vertex is empty");
		return 0;
	}
	if (mesh.F.size() == 0) {
		meshLogger->error("Mesh'Facet is empty");
		return 0;
	}

	if (!build_fine_mesh_projection_tree(args.fine_mesh_file) &&
		!args.fine_mesh_file.empty()) {
		meshLogger->warn("Failed to build fine mesh projection tree: {}", args.fine_mesh_file);
	}

	threadsInitialized = false; // dt_init starts a new run, even on a reused DT object.
	initializeDTThreads(*this, args, mesh.V.size());

    if (infolevel > 0) {
        meshLogger->info("constrain: {}", args.constrain);
        if (args.refine == 1)
            meshLogger->info("refine: {} size={} minEdge={} maxEdge={} growsize={}",
                args.refine, args.size, args.minEdge, args.maxEdge, args.growsize);
        else
            meshLogger->info("refine: {}", args.refine);
        meshLogger->info("nthread: {} (requested={}, initial={})", num_threads,
            args.nthread, activeDTThreads(*this, mesh.V.size()));
        meshLogger->info("optlevel: {} optloop={} optanglestrict={} optTh={} adpangle={}",
            args.optlevel, args.optloop, args.optanglestrict, args.optTh, args.adpangle);
    }

	for (int i = 0; i < args.periodic_P.size() / 2; i++) {
		int a = args.periodic_P[i * 2];
		int b = args.periodic_P[i * 2 + 1];

		auto& va = periodic_P[a];
		if (std::find(va.begin(), va.end(), b) == va.end()) {
			va.push_back(b);
		}

		auto& vb = periodic_P[b];
		if (std::find(vb.begin(), vb.end(), a) == vb.end()) {
			vb.push_back(a);
		}
	}
	return 1;
}

int DT::BndPntInst(Mesh& mesh, Args& args)
{
    MeshStageLog stageLog(*this, "Boundary points", 1);
	int i, j;
	double v1[3], v2[3], n[3];
	// Read input Pnts
	buildPntInfo(mesh);
	if (infolevel > 0)
		meshLogger->debug("Delaunizing boundary points.");
	/**set ghost at nSurNodes**/
	//init ghost nodes;store at nSurNodes
	ghost = addNode();
	/********************* Hilbert sort input nodes ******************/
	if (infolevel > 0) meshLogger->debug("Hilbert sort");
	std::vector<int> order;
	order.resize(nSurNodes);
	for (i = 0; i < nSurNodes; i++)
		order[i] = i;
	Hilbert(mesh.V, order);
	/************************ init first tet *************************/
	{
		if (infolevel > 0) meshLogger->debug("Create first tet");
		double epsilon = 1e-30;
		// Calculate the diagonal size of its bounding box.
		double boxsize = sqrt(norm2(maxW[0] - minW[0], maxW[0] - minW[0], maxW[2] - minW[0]));

		// Make sure the second vertex is not identical with the first one.
		i = 1;
		while ((distance(Nodes[order[0]].pt, Nodes[order[i]].pt) / boxsize) < epsilon && i < nSurNodes - 1) {
			i++;
		}
		if (i > 1) { // Swap i to index 1.
			std::swap(order[i], order[1]);
		}
		// Make sure the third vertex is not collinear with the first two.
		i = 2;
		for (j = 0; j < 3; j++) {
			v1[j] = Nodes[order[1]].pt[j] - Nodes[order[0]].pt[j];
			v2[j] = Nodes[order[i]].pt[j] - Nodes[order[0]].pt[j];
		}
		cross(v1, v2, n);
		while ((lenvec(n) / boxsize * boxsize) < epsilon && i < nSurNodes - 1) {
			i++;
			for (j = 0; j < 3; j++) {
				v2[j] = Nodes[order[i]].pt[j] - Nodes[order[0]].pt[j];
			}
			cross(v1, v2, n);
		}
		if (i > 2) {
			std::swap(order[i], order[2]);// Swap i to index 2.
		}
		// Make sure the fourth vertex is not coplanar with the first three.
		i = 3;
		double ori = dt::GEOM_FUNC::orient3d(Nodes[order[0]].pt, Nodes[order[1]].pt, Nodes[order[2]].pt, Nodes[order[i]].pt);
		while (fabs(ori) < epsilon && i++ < nSurNodes - 1) {
			ori = dt::GEOM_FUNC::orient3d(Nodes[order[0]].pt, Nodes[order[1]].pt, Nodes[order[2]].pt, Nodes[order[i]].pt);
		}

		if (i == nSurNodes){//plane
			if (infolevel > 0) meshLogger->debug("It is a 2D model.");
			return 0;
		}

		if (i > 3) {// Swap i to index 3.
			std::swap(order[i], order[3]);
		}

		//   right-hand rule.
		if (ori > 0.0) {// Swap the first two vertices.
			std::swap(order[0], order[1]);
		}
		// Create first tet.
		int t_first = addElem(order[0], order[1], order[2], order[3]);

		// Create four hull tet.
		int t_0 = addElem(order[1], order[2], order[3], ghost);
		int t_1 = addElem(order[2], order[0], order[3], ghost);
		int t_2 = addElem(order[0], order[1], order[3], ghost);
		int t_3 = addElem(order[0], order[2], order[1], ghost);

		//connect first tet with four hull
		bond(t_first, 0, t_0, 3);
		bond(t_first, 1, t_1, 3);
		bond(t_first, 2, t_2, 3);
		bond(t_first, 3, t_3, 3);

		//connect four hull
		bond(t_3, 0, t_0, 2);
		bond(t_3, 2, t_1, 2);
		bond(t_3, 1, t_2, 2);
		bond(t_0, 0, t_1, 1);
		bond(t_1, 0, t_2, 1);
		bond(t_2, 0, t_0, 1);
	}

	insertDelaunayPoints(*this, order);
	
	AddBox(2.0);
	return 1;
}

void DT::AddBox(double scaled) {
	addBoxFlag = true;

	double C[3], H[3], Hx[3];
	for (int i = 0; i < 3; i++) {
		C[i] = 0.5 * (minW[i] + maxW[i]);      // center
		H[i] = 0.5 * (maxW[i] - minW[i]);      // half size
		Hx[i] = H[i] * scaled;                 // scaled half size
	}

	// 放大后的 min / max
	double minX[3], maxX[3];
	for (int i = 0; i < 3; i++) {
		minX[i] = C[i] - Hx[i];
		maxX[i] = C[i] + Hx[i];
	}

	// 8 个角点（AABB 标准顺序）
	double box[8][3] = {
		{minX[0], minX[1], minX[2]},
		{maxX[0], minX[1], minX[2]},
		{maxX[0], maxX[1], minX[2]},
		{minX[0], maxX[1], minX[2]},
		{minX[0], minX[1], maxX[2]},
		{maxX[0], minX[1], maxX[2]},
		{maxX[0], maxX[1], maxX[2]},
		{minX[0], maxX[1], maxX[2]}
	};

	int searchtet = 0;
	while (isDelEle(searchtet)||ishulltet(searchtet))
		searchtet++;

	for (int i = 0; i < 8; i++) {
		int iNod = addNode(box[i][0], box[i][1], box[i][2], 0);
		std::vector<int> firsttet = { searchtet };
		int BW_recall = BW_insert_vertex(iNod, firsttet, 0);
		searchtet = firsttet[0];
		if (BW_recall == 1) {//success
			setbndpnt(iNod);
			continue;
		}
		else {
			BW_insert_vertex(iNod, firsttet, 0);
			meshLogger->error("Inserting box point failed.");
		}
	}
	return;
}

int DT::BoudaryRecover(Mesh& mesh, Args& args) {
    MeshStageLog stageLog(*this, "Boundary recovery", 1);
	if (infolevel > 0) meshLogger->debug("Recovering boundaries.");

	// Read input Faces
	buildBndInfo(mesh, args);

	if (args.constrain == 0) {
		//allow modify bnd
		modifyBnd = true;	//set Bnd Modification
		Type_Vertex_Edg(args.adpangle, mesh); //determine the constrained edge by geoinfo
		//flipBndEdgPass(3);
	}

	/************************** Recover Edges **************************/
	if (args.ignoreIntersect != 1 && args.autoflip == 1) {
		AutorecoverEdges(args);
	}
	else {
		//for AutoGrid, use traditional boundary edge recover
		recoverEdgesPass(args);
	}

	/************************** Recover Faces **************************/
	recoverFacesPass(args);

	if (!ignoreIntersect && args.constrain == 1) {
		removeStPass(args);
	}
    else if (!ignoreIntersect && addst > 0) {
        removeInteriorSteiner();
    }


	if (infolevel > 0) {
		meshLogger->debug("Add steiner point: {}", addst);
		meshLogger->debug("Add Boundary point: {}", addstbnd);
	}

	return 0;
}

int DT::AutorecoverEdges(Args& args) {
    MeshStageLog stageLog(*this, "Recover edges by flips", 2);
	std::queue<int> lost;
	std::unordered_map<int, int> N_lostE_pre;

	//find all lost edges
	for (int i = 0; i < Elems.size(); i++) {
		if (isDelEle(i))
			continue;
		for (int j = 0; j < 6; j++) {
			int* edgeid = BndEdg.find(Elems[i].form[Egid[j][0]], Elems[i].form[Egid[j][1]]);
			if (edgeid) {
				SurEdgs[*edgeid].info = 1;
			}
		}
	}

	//if (args.constrain == 0) {
	//	flipBndEdgPass(2);
	//}

	//lost.push(2);
	for (int i = 0; i < SurEdgs.size(); i++) {
		if (!isRecBndEdg(i)) {
			lost.push(i);
			N_lostE_pre[SurEdgs[i].iStart]++;
			N_lostE_pre[SurEdgs[i].iEnd]++;
		}
	}

	if (infolevel > 0)
		meshLogger->debug("Lost Edges: {} / {}", (int)lost.size(), (int)SurEdgs.size());

	/**************************** Start Recover Edge *******************************/
	int tempfliplevel = 1;
	while (lost.size() != 0) {
		int  nlost = lost.size();
		int success = 0;
	
		for (int i = 0; i < nlost; i++) {
			int te = lost.front();
			lost.pop();
			if (isDelSurEdg(te))
				continue;

			//if (periodic_P.size() == 0 && args.constrain == 0) {
			//	if (!SurEdgs[te].constrain && !isMeshEdge(SurEdgs[te].iStart, SurEdgs[te].iEnd) && ifflipEdg(te)) {
			//		flipBndEdge(te);
			//	}
			//}

			int fullsearch = 0;
			int stflag = 0;
			fliplevel = tempfliplevel - SurEdgs[te].info * 10;

			if (SurEdgs[te].info <= -3) {
				fliplevel = std::max(1000, fliplevel);
				fullsearch = 1;
			}

			if (SurEdgs[te].info <= -4) stflag = 1;
			if (SurEdgs[te].info <= -5) stflag = 2;

			int ret = recoverEdge(te, fullsearch, stflag);

			if (ret == 0) {//recover edge fail
				lost.push(te);
			}
			else if (ret > 1) {//recover edge fail,and split it
				for (int j = ret; j < SurEdgs.size(); j++) {
					//if (periodic_P.size() == 0 && args.constrain == 0) {
					//	if (!isMeshEdge(SurEdgs[j].iStart, SurEdgs[j].iEnd) && ifflipEdg(j)) {
					//		flipBndEdge(j);
					//	}
					//}
					if (recoverEdge(j, 1, 0) == 0){
						SurEdgs[j].info = -3;
						lost.push(j);//add new edge wait for recover
					}
				}
				success++;
			}
			else if (ret == 1) {
				success++;
			}
		}

		updateFliptype(N_lostE_pre, lost);

		if (infolevel >= 2)
			meshLogger->debug("level:{} lost edges:{}", tempfliplevel, (int)lost.size());
		tempfliplevel++;
		if (tempfliplevel > 1000) {
			meshLogger->error("Boundary recovery failed!");
			throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
		}
	}

	seg[0] = seg[1] = -1;

	if (lost.empty() && infolevel > 0)
		meshLogger->debug("Finish Recover Edges.");

	return 0;
}

void DT::updateFliptype(std::unordered_map<int, int>& N_lostE_pre, std::queue<int>& lost) {
	std::queue<int> tempQ = lost;
	std::unordered_map<int, int> N_lostE_now;
	while (!tempQ.empty()) {
		int x = tempQ.front();
		tempQ.pop();
		N_lostE_now[SurEdgs[x].iStart]++;
		N_lostE_now[SurEdgs[x].iEnd]++;
	}

	tempQ = lost;
	while (!tempQ.empty()) {
		int x = tempQ.front();
		tempQ.pop();
		int pa = SurEdgs[x].iStart;
		int pb = SurEdgs[x].iEnd;
		int olda = N_lostE_pre[pa];
		int oldb = N_lostE_pre[pb];
		int newa = N_lostE_now[pa];
		int newb = N_lostE_now[pb];

		std::swap(SurEdgs[x].iStart, SurEdgs[x].iEnd);

		if (olda <= newa && oldb <= newb) {
			SurEdgs[x].info = std::max(SurEdgs[x].info - 1, -90);
			if (newa + newb > 100) {
				SurEdgs[x].info = std::min(SurEdgs[x].info, -4);
			}
		}
		else {
			SurEdgs[x].info = std::min(SurEdgs[x].info + 1, 0); //0
		}
	}

	N_lostE_pre.clear();
	for (auto it : N_lostE_now) {
		N_lostE_pre[it.first] = it.second;
	}

	return;
}

int DT::recoverEdgesPass(Args& args) {
    MeshStageLog stageLog(*this, "Recover edges", 2);
	std::queue<int> lost;
	int success = 0;
	//find all lost edges
	for (int i = 0; i < Elems.size(); i++) {
		if (isDelEle(i))
			continue;
		for (int j = 0; j < 6; j++) {
			int* edgeid = BndEdg.find(Elems[i].form[Egid[j][0]], Elems[i].form[Egid[j][1]]);
			if (edgeid) {
				SurEdgs[*edgeid].info = 1;
			}
		}
	}

	for (int i = 0; i < SurEdgs.size(); i++) {
		if (!isRecBndEdg(i)) {
			lost.push(i);
		}
	}

	if (infolevel > 0)
		meshLogger->debug("Lost Edges: {} / {}", (int)lost.size(), (int)SurEdgs.size());

	//-------------------------------- only easy flip -----------------------------
	fliplevel = 1;
	while (lost.size() != 0) {
		success = recoverEdges(lost, 0, 0);

		if (infolevel >= 2)
			meshLogger->debug("level:{} lost edges:{}", fliplevel, (int)lost.size());
		if (success == 0) {
			break;
		}
		fliplevel++;
	}

	fliplevel = std::max(1000, fliplevel);
	while (lost.size() != 0) {
		success = recoverEdges(lost, 1, 0);
		if (infolevel >= 2)
			meshLogger->debug("level:{} lost edges:{}", fliplevel, (int)lost.size());
		if (success == 0) {
			break;
		}
	}

	if (args.ignoreIntersect == 1) {
		while (lost.size() != 0) {
			int lostE = lost.front();
			ignoreE(lostE);
			lost.pop();
			continue;
		}
	}

	//---------- easy flip + fullsearch flip + add inner/middle steiner point -------------
	int trybdrc = 0;
	while (lost.size() != 0) {
		int oldlost = lost.size();
		if (success == 0) {
			if (trybdrc++ > 100) {
				meshLogger->error("Boundary recovery failed!");
				throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
			}
			success = recoverEdges(lost, 1, trybdrc > 1 ? 2 : 1);
		}

		success = recoverEdges(lost, 1, 0);
		if (infolevel >= 2)
			meshLogger->debug("Recover by add steiner, lost edges:{}", (int)lost.size());
	}


	if (lost.empty() && infolevel > 0)
		meshLogger->debug("Finish Recover Edges.");

	return 0;
}

int DT::recoverFacesPass(Args& args) {
    MeshStageLog stageLog(*this, "Recover faces", 2);
	int i, j, success, a, b, c, d;
	std::queue<int> lost;
				
	//find all lost edges
	for (i = 0; i < Elems.size(); i++) {
		if (isDelEle(i))
			continue;
		for (j = 0; j < 4; j++) {
			DNC(j, a, b, c, d);
			int* faceid = BndTri.find(Elems[i].form[b], Elems[i].form[c], Elems[i].form[d]);
			if (faceid) {
				SurTris[*faceid].info = 1;
			}
		}
	}
	for (i = 0; i < SurTris.size(); i++) {
		if (!isRecBndTri(i) && !isDelSurTri(i)) {
			lost.push(i);//find a lost edge
		}
	}
	if (infolevel > 0) meshLogger->debug("Lost Faces: {} / {}", (int)lost.size(), (int)SurTris.size());

	//-------------------------------  flip ---------------------------
	fliplevel_face = 0;
	fliplevel = 1000;

	while (lost.size() != 0) {
		success = recoverFaces(lost, 0);

		if (infolevel >= 2)
			meshLogger->debug("level:{} lost faces:{}", fliplevel_face, (int)lost.size());
		if (fliplevel_face++ > 100)
			break;
		if (success == 0) {
			fliplevel_face = std::max(fliplevel_face, 1000);//set a big flip level
		}
	}

	if (args.ignoreIntersect == 1) {
		while (lost.size() != 0) {
			int lostF = lost.front();
			ignoreF(lostF);
			lost.pop();
		}
	}

    // Keep boundary connectivity while the interior-only attempt is useful.
    if (!args.ignoreIntersect && !lost.empty()) {
        recoverFaces(lost, 1);
        recoverFaces(lost, 0);
    }

	//-----------------------  flip +  + split ---------------------
	int tryloop = 0;
	while (lost.size() != 0) {
		success = recoverFaces(lost, 2);
		success = recoverFaces(lost, 0);

		if (infolevel >= 2)
			meshLogger->debug("level:{} lost faces:{}", fliplevel_face, (int)lost.size());
		if (tryloop++ > 1000) {
			throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
		}
	}


	if (lost.empty() && infolevel > 0)
		meshLogger->debug("Finish Recover Faces.");

	return 0;
}

// Boundary points are retained here; only disposable interior points are tried.
void DT::removeInteriorSteiner() {
    MeshStageLog stageLog(*this, "Remove interior Steiner points", 2);
	for (int i = nSurNodes; i < Nodes.size(); i++) {
		if (isDelNod(i) || isbndpnt(i) || i == ghost) {
			continue;
		}
		if (!removePnt(i, 100)){
			smooth_volume(i, true);
		}
	}
}

int DT::removeStPass(Args& args) {
    MeshStageLog stageLog(*this, "Remove boundary Steiner points", 2);
	if (!args.constrain)
		return 1;

	int nRmvS = 0;
    std::vector<int> remainingOrder;
    remainingOrder.reserve(SteinerOrd.size());
	int it = TriSteiner.size() - 1;
	int ie = EdgSteiner.size() - 1;
	for (int j = SteinerOrd.size() - 1; j >= 0; j--) {
		if (SteinerOrd[j] == 2) {
			if (!removeTriStiner(it)) {//can't remove
                remainingOrder.push_back(2);
				if (infolevel >= 2) meshLogger->warn("Can't remove Tri steiner point:{}", TriSteiner[it].first);
			}
			else {//success remove steiner point
				nRmvS++;
				if (infolevel >= 2) meshLogger->debug("Remove Tri steiner point:{}", TriSteiner[it].first);
				vector<std::pair<int, int>>::iterator iter = TriSteiner.begin() + it;
				TriSteiner.erase(iter);//erase this steiner info
			}
			it--;
		}
		else if (SteinerOrd[j] == 1) {
			int  ret = 1;

			if (EdgSteiner[ie].first != -1)
				ret = removeEdgStiner(ie, 0);
			//checkMeshError();
			if (ret != 1) {//can't remove
                remainingOrder.push_back(1);
				if (infolevel >= 2) meshLogger->warn("Can't remove Edge steiner point:{} in {} {}", EdgSteiner[ie].first,
					SurEdgs[EdgSteiner[ie].second].iStart, SurEdgs[EdgSteiner[ie].second].iEnd);
			}
			else if (ret == 1) {//success remove steiner point
				nRmvS++;
				if (infolevel >= 2) meshLogger->debug("Remove steiner point:{} in {} {}", EdgSteiner[ie].first,
					SurEdgs[EdgSteiner[ie].second].iStart, SurEdgs[EdgSteiner[ie].second].iEnd);
				vector<std::pair<int, int>>::iterator iter = EdgSteiner.begin() + ie;
				EdgSteiner.erase(iter);//erase this steiner info
			}
			ie--;
		}
	}
	std::reverse(remainingOrder.begin(), remainingOrder.end());
    SteinerOrd.swap(remainingOrder);

    removeInteriorSteiner();

	if (infolevel > 0) {
		if (EdgSteiner.size() + TriSteiner.size() == 0) meshLogger->debug("Keep constrain!");
		meshLogger->debug("Remove Bnd steiner: {}", nRmvS);
		meshLogger->debug("Leave Bnd steiner: {}", (int)EdgSteiner.size() + (int)TriSteiner.size());
	}
	return 0;
}
/*
* remove hulltet
* dig hole
*/
int DT::ColorVirtualTet(Args& args) {
    MeshStageLog stageLog(*this, "Classify regions", 1);
	if (infolevel > 0)
		meshLogger->debug("Color outer tet & hole!");


	//color all tet
	int subdomain = ColorTets();
	if (infolevel > 0)
		meshLogger->debug("Have {} SubDomain.", subdomain);

	if (args.extrashell) {
		shellextra();
	}

	//box bnd type clean
	for (int i = ghost; i < ghost + 9; i++) {
		Nodes[i].type = 0;
	}

    // These counts are diagnostic only; do not scan the mesh when silent.
    if (infolevel > 0) {
        size_t innerCount = 0;
        int nEdgSteiner = 0;
        std::vector<unsigned char> counted(Nodes.size(), 0);
        for (int i = ghost; i < static_cast<int>(Nodes.size()); ++i)
            if (!isDelNod(i) && isbndpnt(i)) ++nEdgSteiner;
        for (int i = 0; i < static_cast<int>(Elems.size()); ++i) {
            if (isvirtualtet(i) || isDelEle(i)) continue;
            for (int n : Elems[i].form)
                if (!isbndpnt(n) && !counted[n]) { counted[n] = 1; ++innerCount; }
        }
        meshLogger->debug("{} steiner in volume", innerCount);
        meshLogger->debug("{} steiner in Boundary", nEdgSteiner);
    }

    // Lazy index includes all slots, matching the original geo-based scans.
    // Apply renames sequentially so body-ID collisions retain their old order.
    std::map<int, std::vector<int>> regionElements;
    bool indexed = false;
    auto changeRegion = [&](int from, int to) {
        if (!indexed) {
            for (int t = 0; t < static_cast<int>(Elems.size()); ++t)
                regionElements[Elems[t].geo].push_back(t);
            indexed = true;
        }
        auto found = regionElements.find(from);
        if (found == regionElements.end() || found->second.empty()) return false;
        if (from == to) return true;
        auto& destination = regionElements[to];
        destination.reserve(destination.size() + found->second.size());
        for (int t : found->second) { Elems[t].geo = to; destination.push_back(t); }
        regionElements.erase(found);
        return true;
    };

	//dig hole
	if (args.hole.size() != 0) {
        int firstSearch = 0;
        while (firstSearch < static_cast<int>(Elems.size()) && isDelEle(firstSearch)) ++firstSearch;
        if (firstSearch == static_cast<int>(Elems.size())) return -1;
		for (int i = 0; i < args.hole.size() / 3; i++) {
            int iSrch = firstSearch;
			if (infolevel > 0)
				meshLogger->debug("Hole: {} {} {}", args.hole[i * 3], args.hole[i * 3 + 1], args.hole[i * 3 + 2]);
			int findpoint = addNode(args.hole[i * 3], args.hole[i * 3 + 1], args.hole[i * 3 + 2], 0);
			locate_pnt(findpoint, iSrch);
			DelNod(findpoint);
			int iColor = Elems[iSrch].geo;
            changeRegion(iColor, virtualID);
		}
	}

	//dig layer
    if (!args.layer.empty()) {
        for (int layer : args.layer) {
            if (layer > subdomain) continue;
            changeRegion(layer, virtualID);
        }
    }

	setAllP2T();

    const bool needRegionBoundary = infolevel >= 2 || !args.hole_bnd_vec.empty() || !args.bnd_bodyid.empty();
	std::set<int> setTriGeo;

	//clean outer steiner point,only do in constrain
	for (int i = nSurNodes + 1; i < Nodes.size(); i++) {
		if (isDelNod(i))
			continue;
		int iElm = getP2T(i);
		if (isDelEle(iElm))
			DelNod(i);
	}

	for (int i = 0; i < SurTris.size(); i++) {
		if (isDelSurTri(i))
			continue;
		for (int j = 0; j < 3; j++) {
			if (isDelNod(SurTris[i].form[j])) {
				setDelSurTri(i);
				break;
			}
		}
		if (needRegionBoundary && !isDelSurTri(i))
			setTriGeo.insert(SurTris[i].parent);
	}

	if (needRegionBoundary) {
		int ia, ib, ic, id;
		std::map<int, std::set<int>> t2s;
		std::map<int, int> usedtetid;
        const std::unordered_set<int> holeBoundary(args.hole_bnd_vec.begin(), args.hole_bnd_vec.end());

		for (int i = 0; i < Elems.size(); i++) {
			if (isDelEle(i) || isvirtualtet(i) || ishulltet(i)) {
				continue;
			}
			for (int j = 0; j < 4; j++) {
				if (isbndpnt(Elems[i].form[j])) {
					int tetgeo = Elems[i].geo;
					usedtetid[tetgeo] = 1;
					for (int k = 0; k < 4; k++) {
						DFC(k, ia, ib, ic, id);
						int pa = Elems[i].form[ib];
						int pb = Elems[i].form[ic];
						int	pc = Elems[i].form[id];
						if (auto* boundaryEntry = BndTri.find(pa, pb, pc)) {
							const int boundaryIndex = *boundaryEntry;
							t2s[tetgeo].insert(findTriParent(boundaryIndex));
						}
					}
					break;
				}
			}
		}
		//meshLogger->debug("TetGeo connect to TriGeo:");
		std::map<int, std::set<int>> geo2tet;
		for (const auto& tid : t2s) {
			//std::printf("%d : ", tid.first);
			//for (auto fid : tid.second) {
			//	std::printf("%d ", fid);
			//}
			//std::printf("\n");
			if (args.bnd_bodyid.size() != 0) {
				for (const auto& temptid : args.bnd_bodyid) {
					int belong = 0;
					for (auto tempfid : temptid.second) {
						if (tid.second.find(tempfid) != tid.second.end()) {
							belong++;
						}
					}
					if ((1.0) * belong / temptid.second.size() > 0.6) {
						geo2tet[temptid.first].insert(tid.first);
					}
				}
			}

			//dig hole by tri ID
			if (args.hole_bnd_vec.size() != 0) {
				bool ifdeltet = true;
				for (auto bd : tid.second) {
					if (holeBoundary.count(bd) == 0) {
						ifdeltet = false;
						break;
					}
				}
				if (ifdeltet) {
                    changeRegion(tid.first, virtualID);
				}
			}
		}

		//update body ID
		if (args.bnd_bodyid.size() != 0) {
			//for (auto tid : args.bnd_bodyid) {
				//std::printf("%d : ", tid.first);
				//for (auto fid : tid.second) {
				//	std::printf("%d ", fid);
				//}
				//std::printf("\n");
			//}
			//for (auto it : geo2tet) {
			//	printf("body id: %d connect to ", it.first);
			//	for (auto itt : it.second) {
			//		printf("%d ", itt);
			//	}
			//	printf("\n");
			//}

			for (auto it = geo2tet.begin(); it != geo2tet.end();) {
				if (it->second.size() == 1) {
					int oldbodyid = *it->second.begin();
					int newbodyid = it->first;
					//std::printf("old body id: %d new body id:%d\n", oldbodyid, newbodyid);
					if (usedtetid[newbodyid] == 1) {
						int tempbodyid = newbodyid + 1;
						while (usedtetid[tempbodyid] == 1) {
							tempbodyid++;
						}
						usedtetid[tempbodyid] = 1;

                        const bool changeid = changeRegion(newbodyid, tempbodyid);
						if (changeid) {
							for (auto& pair : geo2tet) {
								// get set
								auto& mySet = pair.second;
								// clean x
								auto it2 = mySet.find(newbodyid); // find x
								if (it2 != mySet.end()) {
									mySet.erase(it2);
									mySet.insert(tempbodyid);
								}
							}
						}
					}

                    changeRegion(oldbodyid, newbodyid);
					// clean geo2tet
					for (auto& pair : geo2tet) {
						// get set
						auto& mySet = pair.second;
						// clean x
						auto it2 = mySet.find(oldbodyid); // find x
						if (it2 != mySet.end()) {
							mySet.erase(it2);
						}
					}
					it = geo2tet.begin();
				}
				else {
					it++;
				}
			}
		}

	}


	for (int i = 0; i < SurTris.size(); i++) {
		//SurTris[i].info > 1 ,this SurTri is divided
		if (isDelSurTri(i) || SurTris[i].info > 1)
			continue;
		SurTris[i].parent = findTriParent(i);
	}

	return 0;
}

int DT::RemoveTet(Args& args) {
    MeshStageLog stageLog(*this, "Remove exterior cells", 1);
	if (infolevel > 0)
		meshLogger->debug("Remove outer tet & Dig hole!");

	if (ignoreIntersect && (args.optlevel!=1)) {
		//for Autogrid, Maintain the suspension point
		for (int i = 0; i < Elems.size(); i++) {
			if (isDelEle(i)) {
				continue;
			}
			if (ishulltet(i)) {
				DelEle(i);
			}
		}
	}
	else {
		//remove hulltet
		for (int i = 0; i < Elems.size(); i++) {
			if (isvirtualtet(i)) {
				DelEle(i);
			}
		}
	}

	setAllP2T();

	for (int i = nSurNodes; i < Nodes.size(); i++) {
		if (isDelNod(i) || isbndpnt(i))
			continue;
		int iElm = getP2T(i);
		if (isDelEle(iElm))
			DelNod(i);
	}

	//check Geo
	//int ia, ib, ic, id, pa, pb, pc;
	//for (int i = 0; i < Elems.size(); i++) {
	//	if (isDelEle(i)) continue;
	//	for (int j = 0; j < 4; j++) {
	//		DNC(j, id, ia, ib, ic);
	//		pa = Elems[i].form[ia];
	//		pb = Elems[i].form[ib];
	//		pc = Elems[i].form[ic];
	//		if (isBndTri(pa, pb, pc)) {//chcek if it is bnd Tri
	//			continue;
	//		}
	//		int neig = getNeig(i, j);
	//		if (Elems[i].geo != Elems[neig].geo)
	//			printf("%d %d\n", i, neig);
	//	}
	//}
	//del ghost
	DelNod(ghost);
	ghost = -1;
	return 0;
}

//add inner point,using default size now
int DT::MeshRefine(Args& args) {
    MeshStageLog stageLog(*this, "Refinement", 1);
	if (!args.refine) {
		return 0;
	}

	if (infolevel > 0) meshLogger->debug("MeshRefine Start");
	if (args.sizingFunc) {
		if (infolevel > 0) meshLogger->debug("Size function control!");
	}
	else if (args.growsize != -1) {//use Boundary transition control
		if (infolevel > 0) meshLogger->debug("Edge grow ratio control: {}", args.growsize);
	}
	else {
		if (infolevel > 0) meshLogger->debug("Uniform Size control: {}", args.size);
	}
	if (minEdge != -1) {
		if (infolevel > 0) meshLogger->debug("Target Min Edge Length: {}", minEdge);
	}
	if (maxEdge != -1) {
		if (infolevel > 0) meshLogger->debug("Target Max Edge Length: {}", maxEdge);
	}
	/*********************** incremental insert inner point *******************/
	if (infolevel > 0) meshLogger->debug("Before mesh refine have tet: {}", (int)Elems.size());

	if (/*num_threads > 1*/ 0) {
		refineBWParallel(*this, args);
	} else {
		int loop = 0;
        std::vector<uint8_t> rejected(Elems.size(), 0);
		while (1) {
			int success = 0;
			int n = Elems.size();
			for (int i = 0; i < n; i++) {
				if (isDelEle(i) || isvirtualtet(i) || ishulltet(i))
					continue;
				bool wasRejected = rejected[i] != 0;
                int iNod = createRefineCandidate(i, args, wasRejected);
                rejected[i] = wasRejected;
				if (iNod == -1)
					continue;//don't need add point

				std::vector<int> tempS = { i };

				if (BW_insert_vertex(iNod, tempS, 2) == 1) {

					success++;
                    rejected.resize(Elems.size(), 0);
                    for (int t : serialBW.plan.newElements) rejected[t] = 0;
					updateSize(iNod, args);
					if (infolevel >= 2) {
						if (Nodes.size() % 1000000 == 0) {
							printMemoryUsage();
						}
					}
				}
				else {
					rejected[i] = 1;
					DelNod(iNod);
				}
			}
			//
			if (infolevel > 0)
				meshLogger->debug("Loop:{: <2} add:{: <5} Elem:{: <8} Node:{: <8}", loop++, success, (int)Elems.size(), (int)Nodes.size());
			if (success == 0)
				break;//can't create new tet economically
		}
	}
	//need clear Elems'info and set point to tet
	for (int i = 0; i < Elems.size(); i++) {
		if (isDelEle(i) || isvirtualtet(i) || ishulltet(i))
			continue;
		Elems[i].info = 0;
		for (int j = 0; j < 4; j++) {
			setP2T(Elems[i].form[j], i);
		}
	}
	if (infolevel > 0) meshLogger->debug("After mesh refine have tet: {}", (int)Elems.size());
	return 1;
}

////////////////////////////////////////////////////////////////
/**************************  algorithm ************************/
////////////////////////////////////////////////////////////////
void DT::buildPntInfo(Mesh& mesh) {
    MeshStageLog stageLog(*this, "Build point data", 2);
	nSurNodes = mesh.V.size();
	nSurTris = mesh.F.size();
	int ntet = mesh.T.size();
	if (infolevel > 0) meshLogger->debug("Input nods: {}", nSurNodes);
	if (infolevel > 0) meshLogger->debug("Input Tris: {}", nSurTris);
	if (ntet != 0 && infolevel > 0)
		meshLogger->debug("Input Tets: {}", mesh.T.size());
	if (mesh.S.size() != 0 && infolevel > 0)
		meshLogger->debug("Input Segs: {}", mesh.S.size());
	//alloc memory
	if (ntet == 0) {
		uint64_t freeM = getFreeMemory();
		if (infolevel > 0)
			meshLogger->debug("Free Memory: {:.3f} GB", 1.0 * freeM / 1024 / 1024 / 1024);

		uint64_t MaxN = freeM / sizeof(Elem) / 10;
		uint64_t MaxE = freeM / sizeof(Node) * 0.8;

		uint64_t ReserveN = std::min<uint64_t>(nSurNodes * 5, MaxN);
		uint64_t ReserveE = std::min<uint64_t>(nSurTris * 10, MaxE);

		Nodes.reserve(static_cast<size_t>(ReserveN));
		Elems.reserve(static_cast<size_t>(ReserveE));
	}
	else {
		Nodes.reserve(nSurNodes * 1.5);
		Elems.reserve(ntet * 1.5);
	}
	maxW[0] = maxW[1] = maxW[2] = -DBL_MAX;
	minW[0] = minW[1] = minW[2] = DBL_MAX;

	//input SurNodes,from 0
	for (int i = 0; i < nSurNodes; i++)
	{
		int iNode = addNode();
		for (int j = 0; j < 3; j++) {
			Nodes[iNode].pt[j] = mesh.V[i][j];
			//Nodes[iNode].pt[j] = std::round(mesh.V[i][j] * 1e10) / 1e10;
			maxW[j] = max(maxW[j], Nodes[iNode].pt[j]);
			minW[j] = min(minW[j], Nodes[iNode].pt[j]);
		}
	}
	dist_max = std::sqrt(std::pow((maxW[0] - minW[0]), 2) +
		std::pow((maxW[1] - minW[1]), 2) +
		std::pow((maxW[2] - minW[2]), 2));
	if (infolevel > 0) meshLogger->debug("Diagonal: {}", dist_max);
	return;
}

// Compatibility entry point for serial callers. Parallel insertion uses frozen
// planning and conflict-free batch commits in dt_bw_parallel.cpp.
int DT::BW_insert_vertex(int iNod, std::vector<int>& srchtet, int info, int thread_n) {
    if (srchtet.empty()) return 0;
    // Node occupation alone cannot protect concurrent location or container growth.
    if (omp_in_parallel()) return -1;
    // This entry is serial and non-reentrant for a given mesh. Reuse its
    // mesh-owned workspace without retaining another mesh's state on the thread.
    BWPlan& plan = serialBW.plan;
    BWRequest request = makeBWRequest(*this, iNod, srchtet, info);
    request.trackAccess = false;
    const BWStatus status = planBW(*this, request, plan);
    if (status == BWStatus::Duplicate) return -plan.duplicateNode;
    if (status != BWStatus::Ready) return -1;

    // Preserve the legacy contention check for threaded callers outside a team.
    // Release only locks acquired here; an enclosing operation may own others.
    struct Occupation {
        DT& mesh;
        int owner;
        std::vector<int> acquired;
        ~Occupation() {
            for (int node : acquired) mesh.clearOccupying(mesh.Nodes[node].occupying, owner);
        }
    } occupation{*this, thread_n, {}};
    if (thread_n >= 0) {
        occupation.acquired.reserve(plan.lockNodes.size());
        for (int node : plan.lockNodes) {
            int expected = -1;
            if (Nodes[node].occupying.compare_exchange_strong(expected, thread_n))
                occupation.acquired.push_back(node);
            else if (expected != thread_n) return -1;
        }
    }
    std::vector<int>& slots = serialBW.slots;
    slots.clear();
    slots.reserve(plan.faces.size());
    try {
        for (size_t i = 0; i < plan.faces.size(); ++i)
            slots.push_back(thread_n < 0 ? addElem() : addElem(thread_n));
    } catch (...) {
        for (int slot : slots) {
            if (thread_n < 0) DelEle(slot); else DelEle(slot, thread_n);
        }
        throw;
    }
    commitBW(*this, plan, iNod, slots);
    finishBW(*this, plan, thread_n);
    if (info == 0) srchtet[0] = slots.back();
    if (thread_n >= 0) clearOccupying(Nodes[iNod].occupying, thread_n);
    return 1;
}

//locate the point
//negtive: have same point,return -idx [-,0]
//1:int tet or hull
//[4,14]:on edge (01,00)~(11,10),two point idx
//[16,32,48,64]:on face,(idx>>4)-1
//#pragma optimize("",off)
int DT::locate_pnt(int iNod, int& searchTet)
{
	int  i, ia, ib, ic, id, p[3], mOrtMin, prevElem = -1, findtry = 0;
	double ori, ortd, ortMin, init = searchTet;
	Node* SrchPt = &Nodes[iNod];
	int zero[2];

	if (ishulltet(searchTet)) {//init search tet don't in hull
		searchTet = getNeig(searchTet, 3);
	}

	// Walk through tetrahedra to locate the point.
	while (1)
	{
		if (iNod > ghost && findtry++ > 10000)
			return 100;
		ortMin = DBL_MAX;
		mOrtMin = -1;
		zero[0] = zero[1] = -1;
		if (ishulltet(searchTet) || (prevElem != -1 && searchTet == init))
			return 1;
		// Let searchtet be the face such that 'searchpt' lies above to it.
		for (i = 0; i < 4; i++) {
			if (getNeig(searchTet, i) == prevElem)
				continue;
			DNC(i, ia, ib, ic, id);
			p[0] = Elems[searchTet].form[ib];
			p[1] = Elems[searchTet].form[ic];
			p[2] = Elems[searchTet].form[id];
			ori = dt::GEOM_FUNC::orient3d(Nodes[p[0]].pt, Nodes[p[1]].pt, Nodes[p[2]].pt, SrchPt->pt);
			if (ori < /*-1e-14*/0.0) {
				// Get vertical distance
				ortd = ori / std::max(calArea(Nodes[p[0]].pt, Nodes[p[1]].pt, Nodes[p[2]].pt), 1e-30);
				if (ortd < ortMin) {
					mOrtMin = i;
					ortMin = ortd;
				}
			}
			else if (ori == 0) {
				if (zero[0] == -1) zero[0] = i;
				else if (zero[1] == -1) zero[1] = i;
				else { //three face ori==3,have same point
					if (ib != zero[0] && ib != zero[1])return p[0] * -1;
					if (ic != zero[0] && ic != zero[1])return p[1] * -1;
					if (id != zero[0] && id != zero[1])return p[2] * -1;
				}
			}
		}

		if (mOrtMin != -1) { /* find least neigtive,it is closer to right loc */
			prevElem = searchTet;
			searchTet = getNeig(searchTet, mOrtMin); /* find the next */
		}
		else
			break;//no negtive value,point in tet
	}
	if (zero[0] != -1 && zero[1] == -1) {//point in face
		return (zero[0] + 1) << 4;//point in face,tet's (x>>4)-1 ,()
	}
	if (zero[0] != -1 && zero[1] != -1) {//point in edge
		ib = ic = -1;
		for (i = 3; i >= 0; i--) {
			if (i != zero[0] && i != zero[1]) {
				if (ib == -1)ib = i;
				else if (ic == -1)ic = i;
			}
		}
		return ib << 2 | ic; //[(1,0),(3,2)]=[0100,1110]=[4,14]
	}
	return 1;//return 1,in tet
}
//#pragma optimize("",on)

//#pragma optimize("",off)
int DT::findShell(const int t, const int p2, const int p3, 
	std::vector<int>& Shell, std::vector<int>& shell_point, int thread_n) {
	int pStart, pEnd, pa, pb, p0 = 0, p1 = 0, i, neigOrd, nowt = t, n = 0;
	shell_point.clear();
	Shell.clear();

	DDNC(p0, p1, p2, p3);
	pStart = Elems[t].form[p0];
	pEnd = Elems[t].form[p1];
	pa = Elems[t].form[p2];
	pb = Elems[t].form[p3];
	shell_point.emplace_back(pStart);
	shell_point.emplace_back(pEnd);
	Shell.emplace_back(t);
	n++;
	while (pEnd != pStart) {
		neigOrd = getNeigOrd(nowt, p0);//neig tet's direction
		int oldnowt = nowt;
		nowt = getNeig(nowt, p0);
		if (nowt == -1) {
			checkMeshError();
			throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
		}
		pEnd = Elems[nowt].form[neigOrd];//update end pEnd
		if (thread_n != -1 && pEnd != pStart) {
			if (!tryOccupying(Nodes[pEnd].occupying, thread_n)) {
				for (i = 2; i < shell_point.size(); i++)
					clearOccupying(Nodes[shell_point[i]].occupying, thread_n);
				return -1;
			}
		}
		shell_point.emplace_back(pEnd);
		for (i = 0; i < 4; i++) {
			if (i != neigOrd && Elems[nowt].form[i] != pa && Elems[nowt].form[i] != pb) {
				p0 = i;//update p0
				break;
			}
		}
		Shell.emplace_back(nowt);
		n++;
		if (n > 10000) {
			checkMeshError();
			throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
		}
	}
	Shell.resize(n);
	shell_point.resize(n);
	return Shell.size();
}
//#pragma optimize("",on)

// Read-only traversal: no shared element visit bits and no P2T repairs.
int DT::findSphere(int node, std::vector<int>& sphere) {
    sphere.clear();
    if (node < 0 || node >= static_cast<int>(Nodes.size())) return 0;
    int seed = getP2T(node);
    if (seed < 0 || seed >= static_cast<int>(Elems.size()) || isDelEle(seed) || isNod_in_Tet(node, seed) < 0) {
        seed = -1;
        for (int t=0; t<static_cast<int>(Elems.size()); ++t)
            if (!isDelEle(t) && isNod_in_Tet(node,t)>=0) { seed=t; break; }
    }
    if (seed < 0) return 0;
    std::unordered_set<int> visited;
    visited.insert(seed); sphere.push_back(seed);
    for (size_t k=0; k<sphere.size(); ++k) {
        const int t=sphere[k];
        const int ord=isNod_in_Tet(node,t);
        for (int face=0; face<4; ++face) {
            if (face==ord) continue;
            const int next=getNeig(t,face);
            if (next<0 || next>=static_cast<int>(Elems.size()) || isDelEle(next) || isNod_in_Tet(node,next)<0) continue;
            if (visited.insert(next).second) sphere.push_back(next);
        }
    }
    return static_cast<int>(sphere.size());
}


int DT::findSphere_pnt(const int p, std::unordered_set<int>& Sphere_pnt) {
	int iElm, ord, src[4], i, nig;
	std::queue<int> que;
	Sphere_pnt.clear();
	std::vector<int> Sphere;
	iElm = getP2T(p);
	set_bit(Elems[iElm].info, 31);//o_bit mean this tet has been visited
	que.push(iElm);
	while (!que.empty()) {
		iElm = que.front();
		que.pop();
		Sphere.emplace_back(iElm);

		if ((ord = isNod_in_Tet(p, iElm)) == -1) {//p don't in iElm
			for (i = 0; i < Sphere.size(); i++) {
				clear_bit(Elems[Sphere[i]].info, 31);//Clear position 0 bit
			}
			return -1;
		}
		DNC(ord, src[0], src[1], src[2], src[3]);
		for (i = 1; i < 4; i++) {
			Sphere_pnt.insert(Elems[iElm].form[src[i]]);
			nig = getNeig(iElm, src[i]);
			if (nig == -1)
				continue;
			if (get_bit(Elems[nig].info, 31))//neig has been visited
				continue;
			set_bit(Elems[nig].info, 31);//o_bit mean this tet has been visited
			que.push(nig);
		}
	}
	//clear info
	for (i = 0; i < Sphere.size(); i++) {
		clear_bit(Elems[Sphere[i]].info, 31);//Clear position 0 bit
	}
	return Sphere_pnt.size();
}

int DT::findSphere_global(const  int p, std::vector<int>& Sphere) {
	Sphere.clear();
	for (int i = 0; i < Elems.size(); i++) {
		if (isDelEle(i))
			continue;
		if (isNod_in_Tet(p, i) != -1)
			Sphere.push_back(i);
	}
	return Sphere.size();
}

int DT::findSphere_tri(const  int p, std::unordered_set<int>& Sphere_tri) {
	Sphere_tri.clear();

	//find first tri contain p
	int iElm, src[4], it = -1;//it is first triangle contain  p
	std::vector<int> visitE;
	std::queue<int> que;
	iElm = getP2T(p);
	if (isNod_in_Tet(p, iElm) == -1) {
		meshLogger->error("Wrong point to tet! P:{} T:{}", p, iElm);
		throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
	}

	set_bit(Elems[iElm].info, 31);//o_bit mean this tet has been visited
	que.push(iElm);
	visitE.emplace_back(iElm);

	while (!que.empty()) {
		iElm = que.front();
		que.pop();
		int ord = isNod_in_Tet(p, iElm);

		if (ord == -1) {//p don't in iElm
			meshLogger->error("Wrong point to tet! P:{} T:{}", p, iElm);
			throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
		}

		DNC(ord, src[0], src[1], src[2], src[3]);
		for (int i = 1; i < 4; i++) {
			if (auto* boundaryEntry = BndTri.find(Elems[iElm].form[(src[i] + 1) % 4], Elems[iElm].form[(src[i] + 2) % 4], Elems[iElm].form[(src[i] + 3) % 4])) {
				const int boundaryIndex = *boundaryEntry;
				it = boundaryIndex;
			}
			int nig = getNeig(iElm, src[i]);
			if (nig == -1)
				continue;
			if (get_bit(Elems[nig].info, 31))//neig has been visited
				continue;
			set_bit(Elems[nig].info, 31);//o_bit mean this tet has been visited
			que.push(nig);
			visitE.emplace_back(nig);
		}

		if (it != -1)
			break;
	}
	//clear info
	for (int i = 0; i < visitE.size(); i++) {
		clear_bit(Elems[visitE[i]].info, 31);//Clear position 0 bit
	}

	if (it == -1)
		return 0;

	while (!que.empty()) {
		que.pop();
	}

	Sphere_tri.insert(it);

	que.push(it);
	while (!que.empty()) {
		it = que.front();
		que.pop();

		for (int i = 0; i < 3; i++) {
			int p1 = SurTris[it].form[i];
			if (p1 == p)
				continue;

			for (auto neigf : SurEdgs[BndEdg.get(p, p1)].face) {
				if (Sphere_tri.find(neigf) != Sphere_tri.end())
					continue;
				que.push(neigf);
				Sphere_tri.insert(neigf);
			}
		}
	}

	return Sphere_tri.size();
}

//int DT::findSphere_tri_p(const  int p, std::unordered_set<int>& neig_p) {
//	neig_p.clear();
//
//	//find first tri contain p
//	int iElm, src[4], it = -1;//it is first triangle contain  p
//	std::vector<int> visitE;
//	std::queue<int> que;
//	iElm = getP2T(p);
//	if (isNod_in_Tet(p, iElm) == -1) {
//		if (meshLogger->level() != spdlog::level::off) printf("Wrong point to tet! P:%d T:%d\n", p, iElm);
//		meshLogger->error("Wrong point to tet! P:{} T:{}", p, iElm);
//		throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
//	}
//
//	set_bit(Elems[iElm].info, 31);//o_bit mean this tet has been visited
//	que.push(iElm);
//	visitE.emplace_back(iElm);
//
//	while (!que.empty()) {
//		iElm = que.front();
//		que.pop();
//		int ord = isNod_in_Tet(p, iElm);
//
//		if (ord == -1) {//p don't in iElm
//			if (meshLogger->level() != spdlog::level::off) printf("Wrong point to tet! P:%d T:%d\n", p, iElm);
//			meshLogger->error("Wrong point to tet! P:{} T:{}", p, iElm);
//			throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
//		}
//
//		DNC(ord, src[0], src[1], src[2], src[3]);
//		for (int i = 1; i < 4; i++) {
//			if (isBndTri(Elems[iElm].form[(src[i] + 1) % 4], Elems[iElm].form[(src[i] + 2) % 4], Elems[iElm].form[(src[i] + 3) % 4])) {
//				it = BndTri.get(Elems[iElm].form[(src[i] + 1) % 4], Elems[iElm].form[(src[i] + 2) % 4], Elems[iElm].form[(src[i] + 3) % 4]);
//			}
//			int nig = getNeig(iElm, src[i]);
//			if (nig == -1)
//				continue;
//			if (get_bit(Elems[nig].info, 31))//neig has been visited
//				continue;
//			set_bit(Elems[nig].info, 31);//o_bit mean this tet has been visited
//			que.push(nig);
//			visitE.emplace_back(nig);
//		}
//
//		if (it != -1)
//			break;
//	}
//	//clear info
//	for (int i = 0; i < visitE.size(); i++) {
//		clear_bit(Elems[visitE[i]].info, 31);//Clear position 0 bit
//	}
//
//	if (it == -1)
//		return 0;
//
//	while (!que.empty()) {
//		que.pop();
//	}
//
//	std::unordered_map<int, bool> visitf;
//	visitf[it] = true;
//
//	que.push(it);
//	while (!que.empty()) {
//		it = que.front();
//		que.pop();
//
//		for (int i = 0; i < 3; i++) {
//			int p1 = SurTris[it].form[i];
//			if (p1 == p)
//				continue;
//			neig_p.insert(p1);
//			for (auto neigf : SurEdgs[BndEdg.get(p, p1)].face) {
//				if (visitf[neigf])
//					continue;
//				que.push(neigf);
//				visitf[neigf] = true;
//			}
//		}
//	}
//
//	return neig_p.size();
//}

/*
* if p1,p2 is MeshEdge,the first tet have p1,p2,will destroy BW
*/
bool DT::isMeshEdge(const int p1, const int p2, int* tet) {
	if (isDelNod(p1) || isDelNod(p2))
		return false;
	int iElm = getP2T(p1);
	if (isNod_in_Tet(p1, iElm) == -1) {
		return 0;
		meshLogger->error("Wrong point to tet! P:{} T:{}", p1, iElm);
		throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
	}

	std::vector<int> sph;
	findSphere(p1, sph);
	for (int i = 0; i < sph.size(); i++) {
		if (isNod_in_Tet(p2, sph[i]) != -1) {
			if (tet) *tet = sph[i];//need to store info
			return true;
		}
	}
	return false;
}

bool DT::isMeshFace(const int p1, const int p2, const int p3, int* tet) {
	if (isDelNod(p1) || isDelNod(p2) || isDelNod(p3))
		return false;
	int iElm = getP2T(p1);
	if (isNod_in_Tet(p1, iElm) == -1) {
		meshLogger->error("Wrong point to tet! P:{} T:{}", p1, iElm);
		throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
	}

	std::vector<int> sph;
	findSphere(p1, sph);
	for (int i = 0; i < sph.size(); i++) {
		if (isNod_in_Tet(p2, sph[i]) != -1) {
			if (isNod_in_Tet(p3, sph[i]) != -1) {
				if (tet)
					*tet = sph[i];//need to store info
				return true;
			}
		}
	}
	return false;
}

int DT::recoverEdges(std::queue<int>& lost, int fullsearch, int info) {
	//randQueue(lost);
	int  nlost = lost.size();
	int success = 0;
	for (int i = 0; i < nlost; i++) {
		int targetEdge = lost.front();
		lost.pop();
		int ret = recoverEdge(targetEdge, fullsearch, info);
		if (ret == 0) {//recover edge fail
			lost.push(targetEdge);
		}
		else if (ret > 1) {//recover edge fail,and split it
			for (int j = ret; j < SurEdgs.size(); j++) {
				if (recoverEdge(j, fullsearch, 0) == 0)
					lost.push(j);//add new edge wait for recover
			}
			success++;
		}
		else if (ret == 1) {
			success++;
		}
	}
	seg[0] = seg[1] = -1;
	return success;
}

int DT::recoverFaces(std::queue<int>& lost, int info) {
	//randQueue(lost);
	int  nlost = lost.size();
	int success = 0;
	for (int i = 0; i < nlost; i++) {
		int targetTri = lost.front();
		lost.pop();
		int ret = recoverFace(targetTri, info);
		fac[0] = fac[1] = fac[2] = -1;
		if (ret == 0) {//recover Tri fail
			lost.push(targetTri);
		}
		else if (ret > 1) {//recover edge fail,and split it
			for (int j = ret; j < SurTris.size(); j++) {
				lost.push(j);//add new edge wait for recover
			}
			success++;
		}
		else if (ret == 1) {
			success++;
		}
	}
	return success;
}

int DT::recoverEdge(const int targetE, int fullsearch, int info) {
	if (isDelSurEdg(targetE) || isRecBndEdg(targetE))
		return 1;

	SurEdg* lostE = &SurEdgs[targetE];
	seg[0] = lostE->iStart;
	seg[1] = lostE->iEnd;
	int  filpdepth = fliplevel, ret = 0;
	std::vector<int> newN;

	ret = recoverEdgebyFlip(targetE, 0, filpdepth << 1 | 0);

	if (ret == 1) {
		if(!isDelSurEdg(targetE))
			SurEdgs[targetE].info = 1;//set this edge is recovered
		return 1;
	}
	else if (ret < 0) 
		return ignoreE(targetE);
	else if (ret > 1)
		return ret;

	ret = recoverEdgebyFlip(targetE, 1, filpdepth << 1 | 0);

	if (ret == 1) {
		if (!isDelSurEdg(targetE))
			SurEdgs[targetE].info = 1;//set this edge is recovered
		return 1;
	}
	else if (ret < 0) 
		return ignoreE(targetE);
	else if (ret > 1)
		return ret;

	//fullsearch
	if (fullsearch == 1) {
		ret = recoverEdgebyFlip(targetE, 0, filpdepth << 1 | 1);

		if (ret == 1) {
			SurEdgs[targetE].info = 1;//set this edge is recovered
			return 1;
		}
		else if (ret < 0) 
			return ignoreE(targetE);
		else if (ret > 1)
			return ret;
	}

    // info: 0 flip only; 1 interior only; 2 interior then boundary; 3 boundary only.
    if (info == 3) return splitBndEdge(targetE, 1);
    if (info > 0) {
        ret = addinnerSteiner_Edge(targetE, newN, 0);
        if (ret == 0 && info > 1) ret = splitBndEdge(targetE, 1);
        // Retain useful interior points until both edges and faces are recovered.
    }

	return ret;
}

//find all Intersect with target edges
int DT::findIntersectwithEdgs(const int targetE, std::vector<std::array<int, 3>>& Vid) {
	int srctet, intTyp, intCod, ia, ib, ic, id, pa, pb, pc, pd;
	double linep[2][3] = { 0 }, facept[3][3] = { 0 };
	int p1 = SurEdgs[targetE].iStart;
	int p2 = SurEdgs[targetE].iEnd;

	int dir = finddirection(p1, p2, srctet);
	if (dir == -20) {
		//some error happen
		return 0;
	}

	if (isNod_in_Tet(p2, srctet) != -1)
		return 1;

	for (int k = 0; k < 3; k++) {
		linep[0][k] = Nodes[p2].pt[k];//some reback close to point start
		linep[1][k] = Nodes[p1].pt[k];
	}

	while (1) {
		/***************find next tetand next intersect type, update dir******************/
		if (4 <= dir && dir <= 7) {//Across Face
			dir -= 4;
			DNC(dir, id, ia, ib, ic);
			Vid.push_back({ Elems[srctet].form[ia] ,Elems[srctet].form[ib] ,Elems[srctet].form[ic] });
			int neigdir = getNeigOrd(srctet, dir);//To reduce a calculate
			srctet = getNeig(srctet, dir);//found next tet,oppo srctet
			int i;
			for (i = 0; i < 4; i++) {
				if (i == neigdir)
					continue;
				DNC(i, id, ia, ib, ic);
				for (int k = 0; k < 3; k++) {
					facept[0][k] = Nodes[Elems[srctet].form[ia]].pt[k];
					facept[1][k] = Nodes[Elems[srctet].form[ib]].pt[k];
					facept[2][k] = Nodes[Elems[srctet].form[ic]].pt[k];
				}
				dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, nullptr);
				if (intTyp == 0) {//don't intersect
					continue;
				}
				else if (intTyp == 1) {//intersect a point
					int pintersect = intCod == 0 ? ia : (intCod == 1 ? ib : ic);
					int IntersectPnt = Elems[srctet].form[pintersect];
					if (IntersectPnt == p2)
						return 0;//find p2
					else if (IntersectPnt != p1) {//round in differt direction,don choose p1 as ic
						if (!ignoreIntersect) {
							if (!isbndpnt(IntersectPnt)) {
								if (removePnt(IntersectPnt)) {
									Vid.clear();
									return findIntersectwithEdgs(targetE, Vid);
								}
								if (disturbPnt(IntersectPnt)) {
									Vid.clear();
									return findIntersectwithEdgs(targetE, Vid);
								}
							}
							return 0;
							meshLogger->error("Colline happen,{} located in {} {}", IntersectPnt, p1, p2);
							throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
						}
						else {
							if (infolevel > 0) meshLogger->warn("Colline happen,{} located in {} {}", IntersectPnt, p1, p2);
							return -1;
						}
					}
				}
				else if (intTyp == 2) {//intersect with edge
					int nextia = intCod == 0 ? ia : (intCod == 1 ? ib : ic);
					int nextib = (intCod + 1) % 3 == 0 ? ia : ((intCod + 1) % 3 == 1 ? ib : ic);
					//conect nextia,nextib to next dir
					dir = -(nextia << 2 | nextib);
					break;
				}
				else if (intTyp == 3) {//intersect with face
					dir = i + 4;
					break;
				}
			}//for (i = 0; i < 4; i++)
			if (i == 4)
				return 0;
		}
		else if (-14 <= dir && dir <= -1) {//Across Edge
			ia = ((-dir) >> 2) & 3;
			ib = (-dir) & 3;
			DDNC(ic, id, ia, ib);
			pa = Elems[srctet].form[ia];
			pb = Elems[srctet].form[ib];
			pc = Elems[srctet].form[ic];
			for (auto it : Vid) {
				if (it[0] == Elems[srctet].form[ia] && it[1] == Elems[srctet].form[ib])
					return 0;
			}

			Vid.push_back({ Elems[srctet].form[ia] ,Elems[srctet].form[ib] ,-1 });
			double papb_p1p2[3];//Intersection of pa,pb with p1,p2
			for (int k = 0; k < 3; k++) {
				facept[0][k] = Nodes[pa].pt[k];
				facept[1][k] = Nodes[pb].pt[k];
				facept[2][k] = Nodes[pc].pt[k];
			}
			dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, papb_p1p2);
			double papb_p1p2_p2 = distance2(papb_p1p2, Nodes[p2].pt);//Distance from intersection to p2
			if (papb_p1p2_p2 < 1e-10) {
				//due to accuray
				return 0;
			}
			int istart = Elems[srctet].form[id];
			while (1) {
				id = getNeigOrd(srctet, ic);//neig tet's direction
				srctet = getNeig(srctet, ic);
				for (int i = 0; i < 4; i++) {
					if (Elems[srctet].form[i] == pa) { ia = i; }
					else if (Elems[srctet].form[i] == pb) { ib = i; }
					else if (i != id && Elems[srctet].form[i] != pa && Elems[srctet].form[i] != pb) { ic = i; }
				}
				pc = Elems[srctet].form[ic];
				pd = Elems[srctet].form[id];

				if (pc == ghost || pd == ghost)
					continue;
				if (pd == p2) {//found p2
					return 0;
				}
				double ori1 = dt::GEOM_FUNC::orient3d(Nodes[pa].pt, Nodes[pb].pt, Nodes[p1].pt, Nodes[pc].pt);
				double ori2 = dt::GEOM_FUNC::orient3d(Nodes[pa].pt, Nodes[pb].pt, Nodes[p1].pt, Nodes[pd].pt);
				if ((ori1 > 0 && ori2 > 0) || (ori1 < 0 && ori2 < 0))
					continue;
				//check if [p1,p2] intersect with [pa,pb,pc] or [pa,pb,pd]
				bool findnext = false;
				double next_intersect[3];//Intersection of pa,pb with p1,p2
				for (int i = 0; i < 2; i++) {
					int ix = i == 0 ? ia : ib;
					for (int k = 0; k < 3; k++) {
						facept[0][k] = Nodes[Elems[srctet].form[ix]].pt[k];
						facept[1][k] = Nodes[pd].pt[k];
						facept[2][k] = Nodes[pc].pt[k];
					}
					dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, next_intersect);
					if (intTyp == 0) {//don't intersect
						continue;
					}
					else if (intTyp == 1) {//intersect a point
						int pintersect = intCod == 0 ? ix : (intCod == 1 ? id : ic);
						int IntersectPnt = Elems[srctet].form[pintersect];
						if (IntersectPnt == p2)
							return 0;//find p2
						if (Elems[srctet].form[pintersect] != p1) {//round in differt direction,don choose p1 as ic
							if (!ignoreIntersect) {
								if (!isbndpnt(IntersectPnt)) {
									if (removePnt(IntersectPnt)) {
										Vid.clear();
										return findIntersectwithEdgs(targetE, Vid);
									}
									if (disturbPnt(IntersectPnt)) {
										Vid.clear();
										return findIntersectwithEdgs(targetE, Vid);
									}
								}
								return 0;
								meshLogger->error("Colline happen,{} located in {} {}", IntersectPnt, p1, p2);
								throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
							}
							else {
								if (infolevel > 0) meshLogger->warn("Colline happen,{} located in {} {}", IntersectPnt, p1, p2);
								return -1;
							}
						}
					}
					else if (intTyp == 2) {//intersect with edge
						double tempdis2 = distance2(next_intersect, Nodes[p2].pt);
						if (tempdis2 > papb_p1p2_p2) {
							//need to update,because isMeshEdge don't find right dir
							//papb_p1p2_p2 = tempdis2;
							break;
						}
						findnext = true;
						int nextia = intCod == 0 ? ix : (intCod == 1 ? id : ic);
						int nextib = (intCod + 1) % 3 == 0 ? ix : ((intCod + 1) % 3 == 1 ? id : ic);
						//conect nextia,nextib to next dir
						dir = -(nextia << 2 | nextib);
						break;
					}
					else if (intTyp == 3) {//intersect with face
						double tempdis2 = distance2(next_intersect, Nodes[p2].pt);
						if (tempdis2 > papb_p1p2_p2) {
							//need to update,because isMeshEdge don't find right dir
							//papb_p1p2_p2 = tempdis2;
							break;
						}
						findnext = true;
						dir = (i == 0 ? ib : ia) + 4;//oppo ia/ib
						break;
					}
					else {
						return 0;//Too complex,ignore,
					}
				}
				if (findnext)
					break;
				if (pd == istart)
					return 0;
			}
		}
		else if (0 <= dir && dir <= 3) {//Across Vertex
			int IntersectPnt = Elems[srctet].form[dir];
			if (IntersectPnt == p2)
				return 0;//Success recover
			else {
				if (!isbndpnt(IntersectPnt)) {
					if (removePnt(IntersectPnt)) {
						Vid.clear();
						return findIntersectwithEdgs(targetE, Vid);
					}
					if (disturbPnt(IntersectPnt)) {
						Vid.clear();
						return findIntersectwithEdgs(targetE, Vid);
					}
				}
				if (!ignoreIntersect) {
					return 0;
					meshLogger->error("Colline happen,{} located in {} {}", IntersectPnt, p1, p2);
					throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
				}
				else {
					if (infolevel > 0) meshLogger->warn("Colline happen,{} located in {} {}", IntersectPnt, p1, p2);
					return -1;
				}
			}
		}
	}
	return 0;
}

//info's 0's bit mean if goto full search
int DT::recoverEdgebyFlip(const int targetE, int dirflag, int info)
{
	int dir, srctet, ia, ib, ic, id, pa, pb, pc, pd, i, tried = 0;
	SurEdg* lostE = &SurEdgs[targetE];
	int p1 = lostE->iStart;
	int p2 = lostE->iEnd;
	if (dirflag == 1)
		std::swap(p1, p2);

	//std::vector<std::array<int, 3>> Vid;
	//findIntersectwithEdgs(targetE, Vid);

	//Mesh mesh;
	//mesh.V.reserve(Nodes.size());
	//for (i = 0; i < Nodes.size(); i++) {//have ghost
	//	mesh.V.push_back({ Nodes[i].pt[0],Nodes[i].pt[1],Nodes[i].pt[2] });
	//}

	//for (auto it : Vid) {
	//	if (it[2] != -1) {	// Face
	//		mesh.F.push_back({ it[0], it[1], it[2] });
	//	}
	//	else {// Edge
	//		mesh.F.push_back({ it[0], it[1], it[1] });
	//	}
	//}

	//std::ostringstream filename;
	//filename << R"(C:\Users\yfwan\Desktop\BndRcv\MP4\bndrec\mesh_)"
	//	<< std::setw(5) << std::setfill('0') << mp4time
	//	<< ".vtk";

	//dt::writeVTK(filename.str(), mesh, true);
	//filename.str("");
	//filename.clear();
	//filename << R"(C:\Users\yfwan\Desktop\BndRcv\MP4\bndrec\Back_)"
	//	<< std::setw(5) << std::setfill('0') << mp4time++
	//	<< ".vtk";

	//outTempMesh(filename.str());

	while (1) {
		if (isDelSurEdg(targetE))
			return 2;
		//find search direction
		dir = finddirection(p1, p2, srctet);

		if (dir == -20 || tried++ > 1000) {
			dir = finddirection(p1, p2, srctet);
			//some error happen
			return 0;
		}
		else if (0 <= dir && dir <= 3) {//Across Vertex
			int IntersectPnt = Elems[srctet].form[dir];
			if (IntersectPnt == p2)
				return 1;//Success recover
			else {
				if (!isbndpnt(IntersectPnt)) {
					if (removePnt(IntersectPnt)) {
						//success
						continue;
					}
					if (disturbPnt(IntersectPnt)) {
						continue;
					}
					else {
						return splitBndEdge(targetE, -IntersectPnt);
					}
				}
				if (!ignoreIntersect) {
					//if (meshLogger->level() != spdlog::level::off)
					//printf("Colline happen, %d located in %d,%d\n", IntersectPnt, p1, p2);
					meshLogger->warn("Colline happen,{} located in {} {}", IntersectPnt, p1, p2);
					return AttachPnt2Seg(IntersectPnt, targetE);
					//throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
				}
				else {
					//this is bnd point,remove all tris connect to this tets
					if (infolevel > 0) meshLogger->warn("Colline happen,{} located in {} {}", IntersectPnt, p1, p2);
					return -1;
				}
			}
		}
		else if (4 <= dir && dir <= 7) {//Across Face
			dir -= 4;
			DNC(dir, id, ia, ib, ic);
			pa = Elems[srctet].form[ia];
			pb = Elems[srctet].form[ib];
			pc = Elems[srctet].form[ic];
			if (auto* boundaryEntry = BndTri.find(pa, pb, pc)) {
				const int boundaryIndex = *boundaryEntry;//chcek if it is bnd Tri
				if (!ignoreIntersect) {
					meshLogger->error("Intersection Face: {}:{} {} {} | Edge:{} {}", boundaryIndex, pa, pb, pc, p1, p2);
					//outTempMesh("./temp.vtk");
					ignoreE(targetE);
					return 1;
					throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
				}
				else {
					if (infolevel >= 2) meshLogger->warn("Intersection Face: {}:{} {} {} | Edge:{} {}", boundaryIndex, pa, pb, pc, p1, p2);
					return -1;
				}
			}
			std::vector<int> rF = { srctet };//store init and new tet

			if (removeface(rF, dir, info >> 1) == 1) {
				continue;//success remove a face,go to next loop
			}
		}
		else if (-14 <= dir && dir <= -1) {//Across Edge
			ia = ((-dir) >> 2) & 3;
			ib = (-dir) & 3;
			pa = Elems[srctet].form[ia];
			pb = Elems[srctet].form[ib];
			if (isBndEdg(pa, pb)) { //chcek if it is bnd Tri
				if (!ignoreIntersect) {
					meshLogger->error("Intersection Edge: {} {} | {} {}", pa, pb, p1, p2);
					//outTempMesh("./temp.vtk");
					int ret;
					
					try {
						ret = DealIntersect(pa, pb, p1, p2);
					}
					catch (const std::exception& e) {
						throw EXCEPTIONSTRING(
							std::string("error exit in ") +
							__FILE__ + ":" +
							std::to_string(__LINE__) +
							", reason: " + e.what()
						);
					}

					if (ret != -1) {
						return ret;
					}
					else {
						throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
					}
				}
				else {
					if (infolevel >= 2) meshLogger->warn("Intersection Edge: {} {} | {} {}", pa, pb, p1, p2);
					return -1;
				}
			}

			std::vector<int> rE = { srctet };//store init and new tet
			if (removeEdge(rE, ia, ib, info >> 1) == 1) {
				continue;
			}
		}
		else {
			meshLogger->error("Error happen when recover Edge:{} {}", p1, p2);
		}

		/*if info&1==1 ,found all face and edge intersect and remove them*/
		if ((info & 1) == 1) {
			std::vector<std::array<int, 3>> Vid;
			int ret = findIntersectwithEdgs(targetE, Vid);
			int sus = 0;
            // A failed neighboring flip does not rule out other intersections.
            const int searchDepth = std::min(info >> 1, 32);

			for (const auto& it : Vid) {

				if (it[2] != -1) {
					//try to remove face
					if (isMeshFace(it[0], it[1], it[2], &srctet)) {
						for (int t = 0; t < 4; t++) {
							if (Elems[srctet].form[t] != it[0] && Elems[srctet].form[t] != it[1] && Elems[srctet].form[t] != it[2]) {
								dir = t;
								break;
							}
						}
						std::vector<int> rf = { srctet };
						if (removeface(rf, dir, searchDepth) == 1)
							sus++;

					}
				}
				else {
					//try to remove edge
					if (isMeshEdge(it[0], it[1], &srctet)) {
						for (int t = 0; t < 4; t++) {
							if (Elems[srctet].form[t] == it[0]) { ia = t; }
							else if (Elems[srctet].form[t] == it[1]) { ib = t; }
						}
						std::vector<int> re = { srctet };
						ret = removeEdge(re, ia, ib, searchDepth);

						if (ret == 1)
							sus++;

					}
				}
			}
			if (sus > 0)
				return recoverEdgebyFlip(targetE, dirflag, info);
		}
		break;
	}
	return 0;
}


namespace {
// One recovery attempt owns all traversal marks and snapshots. No mesh marks
// are borrowed, and a flip never leaves cached tetrahedron indices in use.
struct FaceRecoveryPatch {
    std::vector<std::array<int, 2>> edges;
    std::vector<std::array<int, 4>> topology;
    std::vector<int> pending, sphere;
    std::unordered_set<int> visited;
    std::unordered_map<uint64_t, int> intersections;

    void collect(DT& mesh, int face) {
        edges.clear(); topology.clear(); pending.clear(); visited.clear(); intersections.clear();
        const std::array<int, 3> target = {mesh.SurTris[face].form[0],
            mesh.SurTris[face].form[1], mesh.SurTris[face].form[2]};
        double tri[3][3];
        for (int i=0;i<3;++i) for (int k=0;k<3;++k) tri[i][k]=mesh.Nodes[target[i]].pt[k];
        auto enqueueTet = [&](int tet) {
            if (tet>=0 && !mesh.isDelEle(tet) && !mesh.ishulltet(tet) && visited.insert(tet).second)
                pending.push_back(tet);
        };
        // Include vertex contacts and coplanar contacts, not only transverse
        // intersections. These connect the cut cells in degenerate positions.
        for (int p:target) {
            sphere.clear(); mesh.findSphere(p,sphere);
            for (int tet:sphere) enqueueTet(tet);
        }
        for (size_t next=0;next<pending.size();++next) {
            const int tet=pending[next];
            std::array<int,4> form;
            for (int k=0;k<4;++k) form[k]=mesh.Elems[tet].form[k];
            auto canonical=form; std::sort(canonical.begin(),canonical.end());
            topology.push_back(canonical);
            for (int a=0;a<4;++a) for (int b=a+1;b<4;++b) {
                const int u=std::min(form[a],form[b]),v=std::max(form[a],form[b]);
                const uint64_t key=(uint64_t(uint32_t(u))<<32)|uint32_t(v);
                auto hit=intersections.emplace(key,0);
                if (hit.second) {
                    double line[2][3]; int code=-1;
                    for (int k=0;k<3;++k) {line[0][k]=mesh.Nodes[u].pt[k];line[1][k]=mesh.Nodes[v].pt[k];}
                    dt::GEOM_FUNC::lin_tri_intersect3d(line,tri,&hit.first->second,&code,nullptr);
                    const int type=hit.first->second;
                    if ((type==dt::GEOM_FUNC::LTI_INTERSECT_FAC || type==dt::GEOM_FUNC::LTI_INTERSECT_EDG)
                        && std::find(target.begin(),target.end(),u)==target.end()
                        && std::find(target.begin(),target.end(),v)==target.end()) edges.push_back({u,v});
                }
                if (hit.first->second!=dt::GEOM_FUNC::LTI_INTERSECT_NUL) {
                    // Traverse the two faces containing this intersecting edge.
                    for (int k=0;k<4;++k) if (k!=a && k!=b) enqueueTet(mesh.getNeig(tet,k));
                }
            }
        }
        std::sort(topology.begin(),topology.end());
    }
};
}

int DT::recoverFacebyLocalFlips(const int targetF) {
    const std::array<int,3> target={SurTris[targetF].form[0],SurTris[targetF].form[1],SurTris[targetF].form[2]};
    FaceRecoveryPatch patch;
    std::set<std::vector<std::array<int,4>>> tried;
    std::vector<int> shell;
    for (;;) {
        if (isMeshFace(target[0],target[1],target[2])) return 1;
        patch.collect(*this,targetF);
        if (patch.edges.empty() || !tried.insert(patch.topology).second) return 0;
        for (const auto& edge:patch.edges) {
            int tet=-1;
            if (isBndEdg(edge[0],edge[1]) || !isMeshEdge(edge[0],edge[1],&tet)) continue;
            shell.assign(1,tet);
            removeEdge(shell,isNod_in_Tet(edge[0],tet),isNod_in_Tet(edge[1],tet),fliplevel_face);
            if (isMeshFace(target[0],target[1],target[2])) return 1;
        }
        // Even a failed edge removal can retain intermediate flips during
        // recovery. Compare actual topology, rather than just success counts.
    }
}

//info == Insertion point method
int DT::recoverFace(const int targetF, int info) {
    if (isDelSurTri(targetF) || isRecBndTri(targetF)) return 1;
    const std::array<int,3> target = {SurTris[targetF].form[0], SurTris[targetF].form[1], SurTris[targetF].form[2]};
    for (int j=0;j<3;++j) fac[j]=target[j];
    for (int j=0;j<3;++j) {
        const int* edge = BndEdg.find(target[j],target[(j+1)%3]);
        if (!edge) return 0;
        const int e=*edge;
        if (isMeshEdge(SurEdgs[e].iStart,SurEdgs[e].iEnd)) continue;
        SurEdgs[e].info=-1;
        const int firstNewFace=static_cast<int>(SurTris.size());
        const int ret=recoverEdge(e,1,info);
        for(int k=0;k<3;++k) fac[k]=target[k];
        seg[0]=seg[1]=-1;
        if (static_cast<int>(SurTris.size())>firstNewFace) return firstNewFace;
        if (ret!=1 || !isMeshEdge(SurEdgs[e].iStart,SurEdgs[e].iEnd)) return 0;
    }
    // Try every anchor edge with flips before considering a surface split.
    int ret=recoverFacebyFlip_Split(targetF,0);
    if (ret==0 && (info>0 || fliplevel_face>=1000)) ret=recoverFacebyLocalFlips(targetF);
    if (ret==0 && (info==1 || info==2)) {
        std::vector<int> newNodes;
        ret=recoverFacebyaddinSt(targetF,newNodes,info);
        if(ret==0 && !newNodes.empty()) {
            ret=recoverFacebyFlip_Split(targetF,0);
            if(ret==0) ret=recoverFacebyLocalFlips(targetF);
        }
    }
    if(ret==0 && info>=2) ret=recoverFacebyFlip_Split(targetF,2);
    if(ret==1) { if(!isDelSurTri(targetF)) SurTris[targetF].info=1; return 1; }
    if(ret>1) return ret;
    if(ret==-1) return ignoreF(targetF);
    return 0;
}

int DT::recoverFacebyaddinSt(const int targetF, std::vector<int>& newN, int info) {
    const int* target=SurTris[targetF].form;
    if (isMeshFace(target[0],target[1],target[2])) return 1;
    double tri[3][3],center[3]={};
    for (int i=0;i<3;++i) for (int k=0;k<3;++k) {
        tri[i][k]=Nodes[target[i]].pt[k]; center[k]+=tri[i][k]/3.0;
    }
    const double area=calArea(tri[0],tri[1],tri[2]);
    if (!(area>0)) return 0;
    FaceRecoveryPatch patch;
    patch.collect(*this,targetF);
    double best=0,dis[2]={},base[3]={};
    for (const auto& edge:patch.edges) {
        if (isBndEdg(edge[0],edge[1])) continue;
        double line[2][3],intersection[3]; int type=0,code=-1;
        for (int i=0;i<2;++i) for (int k=0;k<3;++k) line[i][k]=Nodes[edge[i]].pt[k];
        dt::GEOM_FUNC::lin_tri_intersect3d(line,tri,&type,&code,intersection);
        if (type!=dt::GEOM_FUNC::LTI_INTERSECT_FAC) continue;
        double heights[2];
        for (int i=0;i<2;++i) heights[i]=-dt::GEOM_FUNC::orient3d(tri[0],tri[1],tri[2],line[i])/(2.0*area);
        if (heights[0]<0) std::swap(heights[0],heights[1]);
        if (!(heights[0]>0 && heights[1]<0)) continue;
        // Favor a residual crossing with room on both sides. Blend its foot
        // towards the centroid to keep the new points away from facet edges.
        const double clearance=std::min(heights[0],-heights[1]);
        if (clearance>best) {
            best=clearance; dis[0]=heights[0]; dis[1]=heights[1];
            for (int k=0;k<3;++k) base[k]=center[k]+0.5*(intersection[k]-center[k]);
        }
    }
    if (!(best>0)) return 0;
    const size_t before=newN.size();
    int ret=addInteriorFacePoints(targetF,newN,dis,base);
    // Keep the established centroid fallback when the residual-based position
    // cannot be inserted at all; do not add another pair after successful BW.
    if (ret==0 && newN.size()==before && distance2(base,center)>0)
        ret=addInteriorFacePoints(targetF,newN,dis,center);
    return ret;
}

//info mean when to split this face,default:7
int DT::recoverFacebyFlip_Split(const int targetF, int info)
{
	int dir = -1, srctet = -1, pa = -1, pb = -1, pc = -1, pd = -1, pe = -1;
	int i = -1, j = -1, k = -1, a = -1, b = -1, c = -1, d = -1, intTyp = -1, intCod = -1;
	double linep[2][3] = { 0 }, facept[3][3] = { 0 }, intPnt[3] = { 0 };//, middle[3], mindis;
	int p1 = SurTris[targetF].form[0];
	int p2 = SurTris[targetF].form[1];
	int p3 = SurTris[targetF].form[2];

	for (k = 0; k < 3; k++) {
		facept[0][k] = Nodes[p1].pt[k];
		facept[1][k] = Nodes[p2].pt[k];
		facept[2][k] = Nodes[p3].pt[k];
	}
	int flipdeapth = fliplevel_face;
	for (i = 0; i < 3; i++) {
		int startPidx = -1;
		while (1) {
			if (i == 0) { pa = p1; pb = p2; pe = p3; }
			else if (i == 1) { pa = p2; pb = p3; pe = p1; }
			else if (i == 2) { pa = p3; pb = p1; pe = p2; }
			//find a tet contain p1,p2;
			if (!isMeshEdge(pa, pb, &srctet)) return 0;
			//dir = finddirection(pa, pb, srctet);
			a = isNod_in_Tet(pa, srctet);
			b = isNod_in_Tet(pb, srctet);

			DDNC(c, d, a, b);
			while (1) {
				pc = Elems[srctet].form[c];
				pd = Elems[srctet].form[d];
				if (pc == pe || pd == pe)
					return 1;//face has been recoverd
				if (startPidx == -1) {
					startPidx = pc;
				}
				else if (startPidx != -1 && startPidx == pc) {
					break;
				}

				if (pc == ghost || pd == ghost) {
					//Non-convex
					std::vector<int> shell, shellp;
					findShell(srctet, c, d, shell, shellp);
					if (shellp.size() == 3 
						&& std::find(shellp.begin(), shellp.end(), p1) != shellp.end()
						&& std::find(shellp.begin(), shellp.end(), p2) != shellp.end()
						&& std::find(shellp.begin(), shellp.end(), p3) != shellp.end()){
						flip32(shell, c, d, -1);
						return isMeshFace(p1,p2,p3) ? 1 : 0;
					}
				}
				else {
					//check if pc,pd intersect with pa,pb,pc
					for (k = 0; k < 3; k++) {
						linep[0][k] = Nodes[pc].pt[k];
						linep[1][k] = Nodes[pd].pt[k];
					}
					dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, intPnt);
					if (intTyp == 3 || intTyp == 2) {//intersect with face
						if (isBndEdg(pc, pd)) {//intersect input mesh
							if (!ignoreIntersect) {
								//outTempMesh("./temp.vtk");
								meshLogger->error("Intersection Face: {} : {} {} {} | Edge:{} {}", targetF, p1, p2, p3, pc, pd);
								ignoreF(targetF);
								return 1;
								throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
							}
							else {
								if (infolevel >= 2) meshLogger->warn("Intersection Face: {} {} {} | Edge:{} {}", p1, p2, p3, pc, pd);
								return -1;
							}
						}

						std::vector<int> shell, shell_pnt;
						findShell(srctet, c, d, shell, shell_pnt);

						//Line's point in face,report intersect face too...
						for (int j = 0; j < 2; j++) {
							double dis = distance(intPnt, linep[j]);
							if (dis == 0) {
								int intsectpnt = j == 0 ? pc : pd;
								if (intsectpnt < ghost) {
									;
								}
								else {
									if (!removePnt(intsectpnt)) {
										double normal[3] = { 0 };
										calnormal(p1, p2, p3, normal);
										if (disturbPnt_search(intsectpnt, normal) == 1) {
											return 0;
										}
										else {
											//only split, no insert
											if (info > 1) {
												std::vector<int> shell;
												return splitBndTri(targetF, shell, intPnt, intsectpnt);
											}
											return 0;
										}
									}
									return 0;
								}
							}
						}

						if (info > 1) {
							for (int j = 0; j < 3; j++) {
								double dis = distance(intPnt, facept[j]);
								if (dis == 0) {
									//Mesh edge intersect with tri point, due to add inner steiner point
									if (disturbPnt(pc) || disturbPnt(pd))
										return recoverFacebyFlip_Split(targetF, 0);
									else {
										//split gravity, but should not be here
										double vv[3] = { 0 };
										return splitBndTri(targetF, shell, vv, -1);
									}
								}
							}
							//spilt this Tris,will return new sub Tri's idx
							return splitBndTri(targetF, shell, intPnt, -2);
						}

						//remove this edge
						std::vector<int> oldtet = { srctet };
						int ret = removeEdge(oldtet, c, d, flipdeapth);
						if (ret == 1)
							i = -1;
						break;//again found pa,pb
					}//if (intTyp == 3)
					else if (intTyp != 0) {
						//outTempMesh("./temp.vtk");
						meshLogger->error("A strange intersect happen:{}", intTyp);
						throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
					}
				}//if (pc != ghost && pd != ghost)
				d = getNeigOrd(srctet, c);
				srctet = getNeig(srctet, c);
				for (j = 0; j < 4; j++) {
					if (j != d && Elems[srctet].form[j] != pa && Elems[srctet].form[j] != pb) {
						c = j;//update p0
						break;
					}
				}
			}//while (1)
			break;
		}
	}//for (i = 0; i < 3; i++)
	return 0;
}
//#pragma optimize("",on)

int DT::ignoreE(int targetE) {
	setDelSurEdg(targetE);
	BndEdg.erase(SurEdgs[targetE].iStart, SurEdgs[targetE].iEnd);
	for (int i = 0; i < SurEdgs[targetE].face.size(); i++) {
		ignoreF(SurEdgs[targetE].face[i]);
	}

	return -1;
}
int DT::ignoreF(int targetF) {
	BndTri.erase(SurTris[targetF].form[0], SurTris[targetF].form[1], SurTris[targetF].form[2]);
	setDelSurTri(targetF);
	for (int k = 0; k < 3; k++) {
		int edgidx = BndEdg.get(SurTris[targetF].form[k], SurTris[targetF].form[(k + 1) % 3]);
		if (edgidx == -1) {
			continue;
		}
		if (SurEdgs[edgidx].face.size() == 1) {
			SurEdgs[edgidx].face.erase(SurEdgs[edgidx].face.begin());
			setDelSurEdg(edgidx);
			BndEdg.erase(SurTris[targetF].form[k], SurTris[targetF].form[(k + 1) % 3]);
		}
		else {
			for (std::vector<int>::iterator it = SurEdgs[edgidx].face.begin(); it != SurEdgs[edgidx].face.end(); it++) {
				if ((*it) == targetF) {
					SurEdgs[edgidx].face.erase(it);
					break;
				}
			}
		}
	}
	return -1;
}
/*
* return 0~3,p1p2 collinear
* return 4+(0~3),intersect with sph[0].nig[0~3]
* return -(a<<2|b)ntersect with srctet' edge_a_b
*/
int DT::finddirection(const int p1, const int p2, int& srctet) {
	int i, ord, ia = 0, ib = 0, ic = 0, id = 0, pb, pc, pd, oldtet, tried = 0;
	double ori_bc, ori_cd, ori_db, startpt[3], endpt[3];
	std::vector<int> visited;

	for (i = 0; i < 3; i++) {
		startpt[i] = Nodes[p1].pt[i];
		endpt[i] = Nodes[p2].pt[i];
	}

	srctet = getP2T(p1);
	if (isNod_in_Tet(p1, srctet) == -1) {
		throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
	}
	if (ishulltet(srctet))
		srctet = getNeig(srctet, 3);
	visited.emplace_back(srctet);
	set_bit(Elems[srctet].info, 0);
	while (1) {
		if (tried++ > 10000)
			return -20;
		if (ishulltet(srctet)) {
			//it must be convex!,find correct tet
			int tord = isNod_in_Tet(p1, srctet);
			oldtet = srctet;
			for (i = 0; i < 4; i++) {
				if (i == tord)
					continue;
				if (i < 3)
					srctet = getNeig(getNeig(oldtet, i), 3);
				else
					srctet = getNeig(oldtet, 3);
				if (!get_bit(Elems[srctet].info, 0))
					break;
			}
			if (i == 4) {
				srctet = oldtet;
				//for (i = 0; i < visited.size(); i++)
				//	clear_bit(Elems[visited[i]].info, 0);

				//double ori = dt::GEOM_FUNC::orient3d(Nodes[Elems[srctet].form[0]].pt,
				//	Nodes[Elems[srctet].form[1]].pt, Nodes[Elems[srctet].form[2]].pt, Nodes[p2].pt);
				//if (ori < 0) {
				//	//non-convex
				//	return 4 + tord;
				//}
				//else if (ori == 0) {
				//	if (tord == 0) return -(1 << 2 | 2);
				//	if (tord == 1) return -(0 << 2 | 2);
				//	if (tord == 2) return -(0 << 2 | 1);
				//}
				return finddirection_global(p1, p2, srctet);

				//srctet = oldtet;
				//return 4 + tord;
				//checkMeshError();
				//return finddirection_global(p1, p2, srctet);
			}
			visited.emplace_back(srctet);
			set_bit(Elems[srctet].info, 0);
		}
		ord = isNod_in_Tet(p1, srctet);//ord p1 of search tet
		DNC(ord, ia, ib, ic, id);
		pb = Elems[srctet].form[ib];
		pc = Elems[srctet].form[ic];
		pd = Elems[srctet].form[id];
		if (pb == p2 || pc == p2 || pd == p2) {
			for (i = 0; i < visited.size(); i++)
				clear_bit(Elems[visited[i]].info, 0);
			if (pb == p2) return ib;
			if (pc == p2) return ic;
			if (pd == p2) return id;
		}

		ori_bc = dt::GEOM_FUNC::orient3d(startpt, Nodes[pb].pt, Nodes[pc].pt, endpt);
		ori_cd = dt::GEOM_FUNC::orient3d(startpt, Nodes[pc].pt, Nodes[pd].pt, endpt);
		ori_db = dt::GEOM_FUNC::orient3d(startpt, Nodes[pd].pt, Nodes[pb].pt, endpt);

		oldtet = srctet;
		// Decide the move direction.
		if (ori_bc > 0) {
			if (ori_cd > 0) {
				if (ori_db > 0) {
					// move to oppo b,c,d
					srctet = getNeig(oldtet, ib);
					if (get_bit(Elems[srctet].info, 0)) {
						srctet = getNeig(oldtet, ic);
						if (get_bit(Elems[srctet].info, 0))
							srctet = getNeig(oldtet, id);
					}
				}
				else {//ori_db < 0
					// move to oppo b,d
					srctet = getNeig(oldtet, ib);
					if (get_bit(Elems[srctet].info, 0))
						srctet = getNeig(oldtet, id);
				}
			}
			else {//ori_cd < 0
				if (ori_db > 0) {
					// move to oppo c,d
					srctet = getNeig(oldtet, ic);
					if (get_bit(Elems[srctet].info, 0))
						srctet = getNeig(oldtet, id);
				}
				else //ori_cd < 0 && ori_db < 0
					// move to oppo d
					srctet = getNeig(oldtet, id);
			}
			visited.emplace_back(srctet);
			set_bit(Elems[srctet].info, 0);
		}
		else {//ori_bc < 0
			if (ori_cd > 0) {
				if (ori_db > 0) {
					// move to oppo b,c
					srctet = getNeig(oldtet, ib);
					if (get_bit(Elems[srctet].info, 0))
						srctet = getNeig(oldtet, ic);
				}
				else {//ori_bc < 0 && ori_db < 0
					// move to oppo b
					srctet = getNeig(oldtet, ib);
				}
				visited.emplace_back(srctet);
				set_bit(Elems[srctet].info, 0);
			}
			else {//ori_bc < 0 && ori_cd < 0
				if (ori_db > 0) {
					// move to oppo c
					srctet = getNeig(oldtet, ic);
					visited.emplace_back(srctet);
					set_bit(Elems[srctet].info, 0);
				}
				else {//ori_bc <= 0 && ori_cd <= 0 && ori_db <= 0
					for (i = 0; i < visited.size(); i++)
						clear_bit(Elems[visited[i]].info, 0);
					if (ori_bc == 0) {
						if (ori_cd == 0) {//ori_bc == 0
							return ic;//coline ac
						}
						if (ori_db == 0) {//ori_bc == 0
							return ib;//coline ab
						}
						return -(ib << 2 | ic);//across edge bc
					}
					if (ori_cd == 0) {
						if (ori_db == 0) {//ori_cd == 0
							return id;//coline ad
						}
						return -(ic << 2 | id);//across edge cd
					}
					if (ori_db == 0) {
						return -(ib << 2 | id);//across edge bd
					}
					// pa->'endpt' crosses the face bcd.
					return 4 + ia;
				}
			}
		}
	}
	for (i = 0; i < visited.size(); i++)
		clear_bit(Elems[visited[i]].info, 0);
	return 0;
}

/*
* return 0~3,p1p2 collinear
* return 4+(0~3),intersect with sph[0].nig[0~3]
* return -(a<<2|b)ntersect with srctet' edge_a_b
* now, it will be non convex
*/
int DT::finddirection_global(const int p1, const int p2, int& srctet) {
	int ia, ib, ic, id, pb, pc, pd, intTyp, intCod;
	double linep[2][3], facept[3][3], intPnt[3];
	std::vector<int> sph;

	findSphere(p1, sph);
	for (int i = 0; i < 3; i++) {
		linep[0][i] = Nodes[p1].pt[i];
		linep[1][i] = Nodes[p2].pt[i];
	}

	for (int i = 0; i < sph.size(); i++) {
		srctet = sph[i];
		if (ishulltet(srctet))
			continue;
		int ord = isNod_in_Tet(p1, srctet);
		DFC(ord, ia, ib, ic, id);
		pb = Elems[srctet].form[ib];
		pc = Elems[srctet].form[ic];
		pd = Elems[srctet].form[id];

		for (int k = 0; k < 3; k++) {
			facept[0][k] = Nodes[pb].pt[k];
			facept[1][k] = Nodes[pc].pt[k];
			facept[2][k] = Nodes[pd].pt[k];
		}

		dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, intPnt);

		if (intTyp == 1) {
			//intersect with Node
			if (intCod == 0) {
				return ib;
			}
			else if (intCod == 1) {
				return ic;
			}
			else if (intCod == 2) {
				return id;
			}
		}
		else if (intTyp == 2) {
			//intersect with Edge,-(a<<2|b)
			if (intCod == 0) {
				return -(ib << 2 | ic);
			}
			else if (intCod == 1) {
				return  -(ic << 2 | id);
			}
			else if (intCod == 2) {
				return  -(id << 2 | ib);
			}
		}
		else if (intTyp == 3) {
			//intersect with Face
			return ia;
		}
	}

	double dismin = DBL_MAX;
	int tempsrc = -1, finaldir = -1;
	for (int i = 0; i < sph.size(); i++) {
		tempsrc = sph[i];
		if (!ishulltet(tempsrc))
			continue;
		int dir = isNod_in_Tet(p1, tempsrc);
		double dis = 0;
		int p3 = -1, p4 = -1;
		for (int j = 0; j < 3; j++) {
			if (dir == j)
				continue;
			if (p3 == -1) p3 = Elems[tempsrc].form[j];
			else p4 = Elems[tempsrc].form[j];
		}

		double v1[3], v2[3], v3[3], c12[3], c23[3];
		for (int j = 0; j < 3; j++) {
			v1[j] = Nodes[p3].pt[j] - Nodes[p1].pt[j];
			v2[j] = Nodes[p2].pt[j] - Nodes[p1].pt[j];
			v3[j] = Nodes[p4].pt[j] - Nodes[p1].pt[j];
		}
		cross(v1, v2, c12);
		cross(v2, v3, c23);

		if (dot(c12,c23)<0)
			continue;

		double sin12 = std::sqrt(c12[0] * c12[0] + c12[1] * c12[1] + c12[2] * c12[2]);
		double sin23 = std::sqrt(c23[0] * c23[0] + c23[1] * c23[1] + c23[2] * c23[2]);

		double cos12 = DT::dot(v1, v2);
		double cos23 = DT::dot(v2, v3);

		double angle1 = RADIO2ANGLE(std::atan2(sin12, cos12));
		double angle2 = RADIO2ANGLE(std::atan2(sin23, cos23));
		
		if ((angle1 < 0 && angle2>0) || (angle1 > 0 && angle2 < 0) || angle1 + angle2>180.0)
			continue;

		srctet = tempsrc;
		double ori = dt::GEOM_FUNC::orient3d(Nodes[p1].pt, Nodes[p2].pt, Nodes[p3].pt, Nodes[p4].pt);
		if (std::fabs(ori) < 1e-20) {
			return -(isNod_in_Tet(p3, tempsrc) << 2 | isNod_in_Tet(p4, tempsrc));
		}
		else {
			return 4 + dir;
		}
	}

	return -20;
}

int DT::DealIntersect(int pa, int pb, int p1, int p2) {
	double linep[2][2][3];

	for (int d = 0; d < 3; ++d) {
		linep[0][0][d] = Nodes[pa].pt[d];
		linep[0][1][d] = Nodes[pb].pt[d];
		linep[1][0][d] = Nodes[p1].pt[d];
		linep[1][1][d] = Nodes[p2].pt[d];
	}

	int intTyp[2] = { 0, 0 };
	int intIdx[2] = { -1, -1 };
	double pnt[3];

	int isInt = dt::GEOM_FUNC::lin_lin_intersect3d_exact(
		linep, intTyp, intIdx, pnt
	);

	if (isInt != 1) {
		// 0: 不相交
		// 2: 共线重叠，不是唯一交点，不适合直接 addNode
		return -1;
	}

	// 如果只想处理两条边内部相交，而不处理端点接触，可以加这个判断
	if (intTyp[0] != dt::GEOM_FUNC::LTI_INTERSECT_EDG ||
		intTyp[1] != dt::GEOM_FUNC::LTI_INTERSECT_EDG) {
		return -1;
	}

	double space = 0.25 * (
		Nodes[pa].space +
		Nodes[pb].space +
		Nodes[p1].space +
		Nodes[p2].space
		);

	int e1 = BndEdg.get(pa, pb);
	if (e1 == -1)
		return -1;
	BndEdg.erase(pa, pb);
	for (int i = 0; i < SurEdgs[e1].face.size(); i++) {
		int Fidx = SurEdgs[e1].face[i];
		BndTri.erase(SurTris[Fidx].form[0], SurTris[Fidx].form[1], SurTris[Fidx].form[2]);
	}

	int e2 = BndEdg.get(p1, p2);
	if (e2 == -1)
		return -1;
	BndEdg.erase(p1, p2);
	for (int i = 0; i < SurEdgs[e2].face.size(); i++) {
		int Fidx = SurEdgs[e2].face[i];
		BndTri.erase(SurTris[Fidx].form[0], SurTris[Fidx].form[1], SurTris[Fidx].form[2]);
	}

	int newp = addNode(pnt[0], pnt[1], pnt[2], space);
	std::vector<int> srchtet = { Nodes[pa].tet };
	int addbw = BW_insert_vertex(newp, srchtet, 3);
	if (addbw <= 0) {//BW_1,add bnd point
		DelNod(newp);
		return -1;
	}

	setbndpnt(newp);

	int ret = AttachPnt2Seg(newp, e1);
	if (ret < 0) {
		return ret;
	}

	int ret2 = AttachPnt2Seg(newp, e2);
	if (ret2 < 0) {
		return ret2;
	}

	return ret;
}

int DT::flipBndEdgPass(int nloop) {
	if (!modifyBnd)
		return 0;

	for (int loop = 0; loop < nloop; loop++) {
		int nsuccess = 0, fail = 0;
		for (int i = 0; i < SurEdgs.size(); i++) {
			if (isDelSurEdg(i) || SurEdgs[i].info >= 1 || SurEdgs[i].constrain > 0 || SurEdgs[i].face.size()!=2)
				continue;
			if (ifflipEdg(i)) {
				if (flipBndEdge(i)) {
					nsuccess++;
				}
				else {
					fail++;
				}
			}
		}
		//outTempMesh("./temp_" + std::to_string(loop) + ".vtk");
		if (infolevel > 0) 
			meshLogger->debug("Flip bndEdg: {}/{}", nsuccess, nsuccess + fail);
		if (nsuccess < 5) 
			break;
	}

	return 0;
}

int DT::flipBndEdge(int index) {
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

	SurTris[f1].info = 0;
	SurTris[f2].info = 0;

	int e1 = BndEdg.get(p1, p4);
	for (int i = 0; i < SurEdgs[e1].face.size(); i++) {
		if (SurEdgs[e1].face[i] == f2) {
			SurEdgs[e1].face[i] = newf1;
			break;
		}
	}

	e1 = BndEdg.get(p1, p3);
	for (int i = 0; i < SurEdgs[e1].face.size(); i++) {
		if (SurEdgs[e1].face[i] == f1) {
			SurEdgs[e1].face[i] = newf1;
			break;
		}
	}

	e1 = BndEdg.get(p2, p3);
	for (int i = 0; i < SurEdgs[e1].face.size(); i++) {
		if ( SurEdgs[e1].face[i] == f1) {
			SurEdgs[e1].face[i] = newf2;
			break;
		}
	}


	 e1 = BndEdg.get(p2, p4);
	for (int i = 0; i < SurEdgs[e1].face.size(); i++) {
		if (SurEdgs[e1].face[i] == f2) {
			SurEdgs[e1].face[i] = newf2;
			break;
		}
	}

	if (periodic_P.size() != 0) {
		if (periodic_P.find(p1) != periodic_P.end() && periodic_P.find(p2) != periodic_P.end()) {
			for (auto it1 : periodic_P[p1]) {
				for (auto it2 : periodic_P[p2]) {
					const auto* boundaryEntry = BndEdg.find(it1, it2);
					if (boundaryEntry && isParallel(Nodes[p1].pt, Nodes[p2].pt, Nodes[it1].pt, Nodes[it2].pt)) {
						const int boundaryIndex = *boundaryEntry;
						int eid = boundaryIndex;
						int ret=flipBndEdge(eid);
						if (ret == 1) {
							return 1;
						}
						else {
							meshLogger->warn("The periodicity may be disrupted during flipBndEdge, {} {}", p1, p2);
						}
					}
				}
			}
		}
	}
	return 1;
}

//return e1:e1 to SurEdgs.size() is new bnd edg
int DT::splitBndEdge(const int lostE, int info) {
	int i, j, newp, p1, p2;
	double pnt[3] = { 0 }, space = 0;
	long double dis1 = 0;
	p1 = SurEdgs[lostE].iStart;
	p2 = SurEdgs[lostE].iEnd;
	std::vector<int> srchtet = { getP2T(p1) };

	double edgelen = distance(Nodes[p1].pt, Nodes[p2].pt);

	if (info == 0) {//midpoint insertion point
		for (i = 0; i < 3; i++) {
			pnt[i] = (Nodes[p1].pt[i] + Nodes[p2].pt[i]) * 0.5;
		}
		space = Nodes[p1].space * 0.5 + Nodes[p2].space * 0.5;
	}
	else if (info == 1) {//insert point at intersection
		std::vector<std::array<int, 3>> Vid;
		if (findIntersectwithEdgs(lostE, Vid) == 1) {
			//have been recover
			return 1;
		}
		if (Vid.size() == 0) {
			return splitBndEdge(lostE, 0);
		}
		std::vector<int> tempvec, shellp;

		double linep[2][3], facept[3][3];
		int srctet, intTyp, intCod;
		for (int k = 0; k < 3; k++) {
			linep[0][k] = Nodes[p1].pt[k];
			linep[1][k] = Nodes[p2].pt[k];
		}

		double best = 2;
		double tempnt[3];

		for (auto it : Vid) {
			if (it[2] == -1) {
				//try to remove edge
				if (isMeshEdge(it[0], it[1], &srctet)) {
					findShell(srctet, isNod_in_Tet(it[0], srctet), isNod_in_Tet(it[1], srctet), tempvec, shellp);
					int j;
					for (j = 0; j < shellp.size(); j++) {
						if (shellp[j] == ghost || shellp[j] == p1 || shellp[j] == p2)
							continue;
						for (int k = 0; k < 3; k++) {
							facept[0][k] = Nodes[it[0]].pt[k];
							facept[1][k] = Nodes[it[1]].pt[k];
							facept[2][k] = Nodes[shellp[j]].pt[k];
						}
						dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, tempnt);

						if (intCod == 0) {
							break;
						}
						else {
							continue;
						}
					}
					if (j == shellp.size()) {
						srchtet.clear();
						srchtet.push_back(getP2T(p1));
						continue;
					}
				}
				else continue;

				double d1 = distance(tempnt, Nodes[p1].pt);
				double d2 = distance(tempnt, Nodes[p2].pt);

				const long double candidateRatio = d1 / (d2 + d1);

				if (std::fabs(0.5 - candidateRatio) < best) {
					srchtet = tempvec;
					best = std::fabs(0.5 - candidateRatio);
					dis1 = candidateRatio;
					for (int k = 0; k < 3; k++)
						pnt[k] = tempnt[k];
				}
				else
					break;
			}
		}

		if(best>1)
		for (auto it : Vid) {
			if (it[2] != -1) {
				if (isMeshFace(it[0], it[1], it[2], &srctet)) {
					//try to remove face
					for (int k = 0; k < 3; k++) {
						facept[0][k] = Nodes[it[0]].pt[k];
						facept[1][k] = Nodes[it[1]].pt[k];
						facept[2][k] = Nodes[it[2]].pt[k];
					}
					dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, tempnt);
					if (intTyp != 3) {
						return splitBndEdge(lostE, 0);
					}

					double d1 = distance(tempnt, Nodes[p1].pt);
					double d2 = distance(tempnt, Nodes[p2].pt);
					const long double candidateRatio = d1 / (d2 + d1);

					if (std::fabs(0.5 - candidateRatio) < best) {
						srchtet.clear();
						srchtet.push_back(srctet);
						int forthp = Elems[srctet].form[0] + Elems[srctet].form[1] + Elems[srctet].form[2] + Elems[srctet].form[3] - it[0] - it[1] - it[2];
						srchtet.push_back(getNeig(srctet, isNod_in_Tet(forthp, srctet)));

						best = std::fabs(0.5 - candidateRatio);
					dis1 = candidateRatio;
						for (int k = 0; k < 3; k++)
							pnt[k] = tempnt[k];
					}
					else
						break;
				}
			}
		}
	
		if (dis1 < 1e-19 || dis1 > 1- 1e-19) {
			return splitBndEdge(lostE, 0);
		}

		space = Nodes[p1].space * (1 - dis1) + Nodes[p2].space * dis1;
	}

	//before add steiner point,we need to deltet BndTri connect to Edge
	for (i = 0; i < SurEdgs[lostE].face.size(); i++) {
		int Fidx = SurEdgs[lostE].face[i];
		BndTri.erase(SurTris[Fidx].form[0], SurTris[Fidx].form[1], SurTris[Fidx].form[2]);
	}

	BndEdg.erase(p1, p2);

	if (info >= 0) {
		newp = addNode(pnt[0], pnt[1], pnt[2], space);

		int ret = BW_insert_vertex(newp, srchtet, 3);
		if (ret <= 0) {//BW_1,add bnd point
			DelNod(newp);
			if (infolevel > 0)
				meshLogger->warn("Add Edge steiner point failed");
			for (i = 0; i < SurEdgs[lostE].face.size(); i++) {
				int Fidx = SurEdgs[lostE].face[i];
				BndTri.add(SurTris[Fidx].form[0], SurTris[Fidx].form[1], SurTris[Fidx].form[2], Fidx);
			}
			BndEdg.add(p1, p2, lostE);
			if (info == 0) {
				return ignoreE(lostE);
			}
			else {
				return splitBndEdge(lostE, 0);
			}
			return 0;
		}
	}
	else {
		newp = -info;
	}

	setbndpnt(newp);
	addstbnd++;
	EdgSteiner.push_back(std::make_pair(newp, lostE));
	SteinerOrd.push_back(1);//edge

	/********************** update bnd info ***********************/
	int  e1, e2, manifold;
	std::vector<int> f;
	std::vector<int> e;

	manifold = SurEdgs[lostE].face.size();//might non mainfold happen
	f.resize(manifold * 3);
	e.resize(manifold);
	for (i = 0; i < manifold; i++) {
		f[i] = SurEdgs[lostE].face[i];
	}
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
	//SurEdgs[e1].parent = lostE;
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
	//SurEdgs[e2].parent = lostE;
	SurEdgs[e2].info = 0;//don't recover
	SurEdgs[e2].face.resize(manifold);
	for (i = 0; i < manifold; i++)
		SurEdgs[e2].face[i] = f[manifold + i * 2 + 1];//first face
	//update lostE's info,it's son is e1
	SurEdgs[lostE].info = e1;
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
		//SurEdgs[e[i]].parent = lostE;
		SurEdgs[e[i]].info = 0;//don't recover
		SurEdgs[e[i]].face.resize(2);
		SurEdgs[e[i]].face[0] = f[manifold + i * 2];//first face
		SurEdgs[e[i]].face[1] = f[manifold + i * 2 + 1];//first face
	}
	//updata old surTri
	for (i = 0; i < manifold; i++) {
		SurTris[f[i]].info = f[manifold + i * 2];
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
		BndTri.add(SurTris[f[i]].form[0], SurTris[f[i]].form[1], SurTris[f[i]].form[2], f[i]);//add this sub bnd Tris
	}

	if (infolevel >= 2)
		meshLogger->debug("Split Edge:{} [{},{}],add steiner:{}", lostE, p1, p2, newp);
	return e1;
}

/*
* info=-1,center of gravity insertion point
* info=-2,intersection insertion
* else, point index
*/
int DT::splitBndTri(const int targetF, std::vector<int>& shell, double intPnt[], int info) {
	int newp = 0;
	const std::array<int, 3> originalForm = {SurTris[targetF].form[0], SurTris[targetF].form[1], SurTris[targetF].form[2]};

	int p1 = originalForm[0];//one of old Tri point
	int p2 = originalForm[1];//one of old Tri point
	int p3 = originalForm[2];//one of old Tri point
	if (info == -1) {//use gravity point
		double newpt[3];
		for (int i = 0; i < 3; i++) {
			newpt[i] = (Nodes[p1].pt[i] + Nodes[p2].pt[i] + Nodes[p3].pt[i]) / 3;
		}
		double space = (Nodes[p1].space + Nodes[p2].space + Nodes[p3].space) / 3.0;
		newp = addNode(newpt[0], newpt[1], newpt[2], space);//new point idx
		shell.clear();
		shell.push_back(getP2T(p1));
	}
	else if (info == -2) {//use intersection point
		double space = (Nodes[p1].space + Nodes[p2].space + Nodes[p3].space) / 3.0;//not right
		newp = addNode(intPnt[0], intPnt[1], intPnt[2], space);//new point idx
	}
	else {
		newp = info;
	}

	//before add steiner point,we need to deltet BndTri connect to Edge
	BndTri.erase(SurTris[targetF].form[0], SurTris[targetF].form[1], SurTris[targetF].form[2]);

	if (info ==-1 || info ==-2) {
		int ret = BW_insert_vertex(newp, shell, 3);
		if (ret <= 0) {//BW_1,add bnd point
			DelNod(newp);
			if (infolevel > 0)
				meshLogger->warn("Add Face steiner point failed");
			BndTri.add(SurTris[targetF].form[0], SurTris[targetF].form[1], SurTris[targetF].form[2], targetF);
			return 0;
			return splitBndTri(targetF, shell, intPnt, -1);
		}
	}

	addstbnd++;
	setbndpnt(newp);
	TriSteiner.push_back(std::make_pair(newp, targetF));//build steiner point info
	SteinerOrd.push_back(2);//Tri
	/*********************** Update BndTris info ***********************/
	int f[3], e[3];
	for (int i = 0; i < 3; i++) {
		f[i] = SurTris.size();
		SurTris.emplace_back(SurTri());//add new surTri
		SurTris[f[i]].form[0] = originalForm[i];
		SurTris[f[i]].form[1] = originalForm[(i + 1) % 3];
		SurTris[f[i]].form[2] = newp;
		SurTris[f[i]].parent = targetF;
		SurTris[f[i]].info = 0;
		BndTri.add(SurTris[f[i]].form[0], SurTris[f[i]].form[1], SurTris[f[i]].form[2], f[i]);
	}
	for (int i = 0; i < 3; i++) {
		e[i] = SurEdgs.size();
		SurEdgs.emplace_back(SurEdg());//add new surEdg
		SurEdgs[e[i]].iStart = newp;
		SurEdgs[e[i]].iEnd = originalForm[i];
		SurEdgs[e[i]].face.resize(2);
		SurEdgs[e[i]].face[0] = f[i];
		SurEdgs[e[i]].face[1] = f[(i + 2) % 3];
		//SurEdgs[e[i]].parent = 0;
		SurEdgs[e[i]].info = 0;
		BndEdg.add(SurEdgs[e[i]].iStart, SurEdgs[e[i]].iEnd, e[i]);
	}
	for (int i = 0; i < 3; i++) {
		p1 = originalForm[i];
		p2 = originalForm[(i + 1) % 3];
		int edgid = BndEdg.get(p1, p2);
		for (int j = 0; j < SurEdgs[edgid].face.size(); j++) {
			if (SurEdgs[edgid].face[j] == targetF) {
				SurEdgs[edgid].face[j] = f[i];
				break;
			}
		}
	}
	//link old Tri to sub Tri
	SurTris[targetF].info = f[0];
	if (infolevel >= 2)
		meshLogger->debug("Split Tri:{}:{} {} {},add steiner:{}", targetF, originalForm[0], originalForm[1], originalForm[2], newp);

	return f[0];
}

int DT::addinnerSteiner_Edge(const int lostE, std::vector<int>& newN, int info) {
	std::vector<std::array<int, 3>> Vid;

	if (findIntersectwithEdgs(lostE, Vid)) { //have been recovered
		return 1;
	}

	if ((Vid.size() >= info && info != 0) || Vid.size() == 0)
		return 0;

	int pa = SurEdgs[lostE].iStart;
	int pb = SurEdgs[lostE].iEnd;
	double p1[3] = {Nodes[pa].pt[0],Nodes[pa].pt[1],Nodes[pa].pt[2]};
    double p2[3] = {Nodes[pb].pt[0],Nodes[pb].pt[1],Nodes[pb].pt[2]};
	double LenEdg = distance(p1, p2);
	double space = (Nodes[pa].space + Nodes[pb].space) / 2.0;
	double linep[2][3], facept[3][3], intPnt[3], pnt[3];
	int intTyp, intCod, srctet, nVid = Vid.size();

	for (int k = 0; k < 3; k++) {
		linep[0][k] = p1[k];
		linep[1][k] = p2[k];
	}

	//int nume = 0, numf = 0;
	//for (auto it : Vid) {
	//	if (it[2] == -1) nume++;
	//	else numf++;
	//}
	//printf("lostE:%d  Vid:%d  info:%d  e:%d  f:%d\n",lostE, Vid.size(),info, nume,numf);
	
	//printf("%d %d %d\n",Vid.size(), pa, pb);
	//std::set<int> alltet;
	//for (auto it : Vid) {
	//	if (it[2] != -1) {	// Face
	//		if (isMeshFace(it[0], it[1], it[2], &srctet)) {
	//			alltet.insert(srctet);
	//			int p4 = Elems[srctet].form[0] + Elems[srctet].form[1] + Elems[srctet].form[2] + Elems[srctet].form[3] - it[0] - it[1] - it[2];
	//			int ie4 = isNod_in_Tet(p4, srctet);
	//			int neig = getNeig(srctet, ie4);
	//			alltet.insert(neig);
	//		}
	//	}
	//	else {// Edge
	//		if (isMeshEdge(it[0], it[1], &srctet)) {
	//			std::vector<int>  shell, shellp;
	//			findShell(srctet, isNod_in_Tet(it[0], srctet), isNod_in_Tet(it[1], srctet), shell, shellp);
	//			for (auto it : shell)
	//				alltet.insert(it);
	//		}//if (isMeshEdge(it[0], it[1], &srctet))
	//	}
	//}//for(auto it:Vid
	//std::vector<int> sph;
	//for (auto it : alltet)
	//	sph.push_back(it);
	//printSph_VTK(sph, "./temp.vtk");
	//throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));

	for(auto it:Vid){
		if (it[2] != -1) {	// Face
			if (isMeshFace(it[0], it[1], it[2], &srctet)) {
				std::vector<int> oldtet = { srctet };
				int dir = isNod_in_Tet(Elems[srctet].form[0] + Elems[srctet].form[1] + Elems[srctet].form[2] + Elems[srctet].form[3] - it[0] - it[1] - it[2], srctet);
				int a, b, c, d;
				DFC(dir, a, b, c, d);
				b = Elems[srctet].form[b];
				c = Elems[srctet].form[c];
				d = Elems[srctet].form[d];
				int adddir = removeface(oldtet, dir, 10, -1);

				if (adddir < 0) {

					int rmvEdg = -adddir - 1, pc = -1, pd = -1;
					if (rmvEdg == 0) { pc = b, pd = d; }
					if (rmvEdg == 1) { pc = d, pd = c; }
					if (rmvEdg == 2) { pc = c, pd = b; }

					//try to remove face
					for (int k = 0; k < 3; k++) {
						facept[0][k] = Nodes[it[0]].pt[k];
						facept[1][k] = Nodes[it[1]].pt[k];
						facept[2][k] = Nodes[it[2]].pt[k];
					}

					dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, intPnt);
					int newp = addNode();
					Nodes[newp].space = space;

					for (int j = 0; j < 3; j++) 
						Nodes[newp].pt[j] = (Nodes[pc].pt[j] + Nodes[pd].pt[j] + intPnt[j]) / 3.0;
				
					isMeshFace(it[0], it[1], it[2], &srctet);
					dir = isNod_in_Tet(Elems[srctet].form[0] + Elems[srctet].form[1] + Elems[srctet].form[2] + Elems[srctet].form[3] - it[0] - it[1] - it[2], srctet);

					int BW_tet = getP2T(pa);
					int loc = locate_pnt(newp, BW_tet);
					if (loc < 1) { DelNod(newp); continue; }

					std::vector<int> BW_vec = { srctet,getNeig(srctet,dir)};

					int ret = BW_insert_vertex(newp, BW_vec, 3);

					if (ret <= 0)
						DelNod(newp);
					else { ++addst; newN.push_back(newp); }
				}			
				if (recoverEdge(lostE, 0, 0) == 1) {
					return 1;
				}
			}
		}
		else {// Edge
			if (isMeshEdge(it[0], it[1], &srctet)) {
				std::vector<int>  shell, shellp;
				findShell(srctet, isNod_in_Tet(it[0], srctet), isNod_in_Tet(it[1], srctet), shell, shellp);

				bool hullflag = false;
				for (int i = 0; i < shellp.size(); i++) {
					if (shellp[i] == ghost) {
						hullflag = true;
						break;
					}
				}

				for (int i = 0; i < shellp.size(); i++) {
					if (shellp[i] == ghost || shellp[i] == pa || shellp[i] == pb)
						continue;
					for (int j = 0; j < 3; j++) {
						facept[0][j] = Nodes[it[0]].pt[j];
						facept[1][j] = Nodes[it[1]].pt[j];
						facept[2][j] = Nodes[shellp[i]].pt[j];
					}
					dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, intPnt);

					if (intCod != 0)
						continue;

					int px = it[1];
					if (distance2(intPnt, Nodes[it[0]].pt) > distance2(intPnt, Nodes[it[1]].pt))
						px = it[0];

					int newp = addNode();
					Nodes[newp].space = space;

					for (int k = 0; k < 3; k++)
						Nodes[newp].pt[k] = (intPnt[k] + Nodes[px].pt[k]) / 2.0;

					int BW_tet = getP2T(pa);
					int loc = locate_pnt(newp, BW_tet);

					if (loc < 1 || loc == 100 || ishulltet(BW_tet)) {
						DelNod(newp);
						continue;
					}

					if (BW_insert_vertex(newp, shell, 3) <= 0) {
						DelNod(newp);
						continue;
					}
					else { ++addst; newN.push_back(newp); }

					// direction
					double norm[3] = { 0 };
					calnormal(pa, pb, newp, norm);
					double normLen = lenvec(norm);
					for (int j = 0; j < 3; j++) norm[j] /= normLen;

					double SA[3] = { Nodes[shellp[i]].pt[0] - p1[0], Nodes[shellp[i]].pt[1] - p1[1], Nodes[shellp[i]].pt[2] - p1[2] };
					double nrom_sa = dot(norm, SA);
					if (nrom_sa > 0.0) {
						norm[0] = -norm[0]; norm[1] = -norm[1]; norm[2] = -norm[2];
					}

					double area = calArea(p1, p2, Nodes[newp].pt);
					double d = area / LenEdg * 2.0;
					double newpos[3] = { 0 }, alpha = 0.5;
					double* verts[4];

					std::vector<int> sph;
					findSphere(newp, sph);

					int iter = 0;
					while (iter++ < 16) {
						bool moveflag = true;
						for (int j = 0; j < 3; j++) newpos[j] = Nodes[newp].pt[j] + d * alpha * norm[j];
						alpha /= 2.0;
						for (int j = 0; j < sph.size(); j++) {
							if (ishulltet(sph[j])) continue;
							for (int m = 0; m < 4; m++) {
								int iElemNd = Elems[sph[j]].form[m];
								verts[m] = (iElemNd != newp) ? Nodes[iElemNd].pt : newpos;
							}
							double ori = dt::GEOM_FUNC::orient3d(verts[0], verts[1], verts[2], verts[3]);
							if (ori >= 0) {
								moveflag = false;
								break; // This tet becomes invalid.
							}
						}
						if (moveflag) {
							for (int j = 0; j < 3; j++)
								Nodes[newp].pt[j] = newpos[j];
							break;
						}
					} // while (iter < 16)
					break;
				} //for (int i = 0; i < shellp.size(); i++)
				if (recoverEdge(lostE, 0, 0) == 1) {
					return 1;
				}
			}//if (isMeshEdge(it[0], it[1], &srctet))
		}
	}//for(auto it:Vid

	return  addinnerSteiner_Edge(lostE, newN, Vid.size());
	return  0;
}

int DT::addinnerSteiner_Edge2(const int lostE, std::vector<int>& newN, int info) {
	std::vector<std::array<int, 3>> Vid;
	int ret = findIntersectwithEdgs(lostE, Vid);

	if (ret == 1) { //have been recovered
		return 1;
	}
	if (Vid.size() == 0) { //have been recovered
		return 0;
	}

	int pa = SurEdgs[lostE].iStart;
	int pb = SurEdgs[lostE].iEnd;
	double* p1 = Nodes[pa].pt;
	double* p2 = Nodes[pb].pt;
	int nf = SurEdgs[lostE].face.size();
	double lenEdg = distance(p1,p2);
	double space = (Nodes[pa].space + Nodes[pb].space) / 2.0;
	double pnt[3] = { 0 };
	double linep[2][3], facept[3][3], intPnt[3];
	int intTyp, intCod, srctet;
	std::vector < std::array<double, 3>> Normal;

	std::set<int> P_all;
	for (auto it : Vid) {
		for (int j = 0; j < 3; j++) {
			if (it[j] == ghost || it[j] == -1)
				continue;
			P_all.insert(it[j]);
		}
	}

	double minlen = lenEdg;
	for (auto pit : P_all) {
		double* p3 = Nodes[pit].pt;
		double area = calArea(p1, p2, p3);
		double height = area / lenEdg * 2.0;
		minlen = std::min(minlen, height);
	}

	double len = minlen < 1e-8 ? 1e-8 : minlen;

	for (int i = 0; i < nf; i++) {
		int sf = SurEdgs[lostE].face[i];

		int pc = SurTris[sf].form[0] + SurTris[sf].form[1] + SurTris[sf].form[2] - pa - pb;
		double  norm[3] = { 0 };
		calnormal(pa, pb, pc, norm);
		double normLen = lenvec(norm);
		if (normLen == 0) {
			if (i == 1) {
				for (int j = 0; j < 3; j++)
					norm[j] = -Normal[0][j];
			}
		}
		else {
			for (int j = 0; j < 3; j++)
				norm[j] /= normLen;
		}

		Normal.push_back({ norm[0] ,norm[1] ,norm[2] });
	}

	//for (auto it : Vid) {
	//	printf("%d %d %d\n", it[0], it[1],it[2]);
	//	//if (isMeshEdge(it[0], it[1], &srctet)) {
	//	//	findShell(srctet, isNod_in_Tet(it[0], srctet), isNod_in_Tet(it[1], srctet), tempvec, shellp);
	//	//	printSph_VTK(tempvec, "./shell.vtk");
	//	//}
	//}
	//printf("%d -> ", Vid.size());


	//printf("%d: ", Vid.size());
	//for (auto it : Vid) {
	//	if(it[2]==-1) printf("e ");
	//	else printf("f ");
	//}
	//printf("\n");

	//根据交点
	//std::unordered_map<int, int> mp;
	//for (auto it : Vid) {
	//	std::vector<int>  tempvec, shellp;

	//	if (mp[it[0]] == 1 || mp[it[1]] == 1 || mp[it[2]] == 1);
	//		//continue;

	//	for (int k = 0; k < 3; k++) {
	//		linep[0][k] = p1[k];
	//		linep[1][k] = p2[k];
	//	}

	//	if (it[2] != -1) {
	//		//try to remove face
	//		for (int k = 0; k < 3; k++) {
	//			facept[0][k] = Nodes[it[0]].pt[k];
	//			facept[1][k] = Nodes[it[1]].pt[k];
	//			facept[2][k] = Nodes[it[2]].pt[k];
	//		}
	//		dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, pnt);
	//		mp[it[0]] = 1;
	//		mp[it[1]] = 1;
	//		mp[it[2]] = 1;
	//	}
	//	else {
	//		if (isMeshEdge(it[0], it[1], &srctet)) {
	//			findShell(srctet, isNod_in_Tet(it[0], srctet), isNod_in_Tet(it[1], srctet), tempvec, shellp);
	//			int j;
	//			for (j = 0; j < shellp.size(); j++) {
	//				if (shellp[j] == ghost || shellp[j] == pa || shellp[j] == pb)
	//					continue;
	//				for (int k = 0; k < 3; k++) {
	//					facept[0][k] = Nodes[it[0]].pt[k];
	//					facept[1][k] = Nodes[it[1]].pt[k];
	//					facept[2][k] = Nodes[shellp[j]].pt[k];
	//				}
	//				dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, pnt);

	//				if (intCod != 0)
	//					continue;

	//				mp[it[0]] = 1;
	//				mp[it[1]] = 1;

	//				break;
	//			}
	//		}
	//	}

	for (int j = 0; j < 3; j++) {
		pnt[j] = (p1[j] + p2[j]) * 0.5;
	}

	for (int i = 0; i < nf; i++) {
		double tempNorm[3] = { 0 };
		for (int j = 0; j < 3; j++) {
			tempNorm[j] = (Normal[i][j] - Normal[(i + 1) % nf][j]) / 2.0;
		}

		double ratio = 0.5;
		int newp = addNode();
		Nodes[newp].space = space;
		while (1) {
			if (ratio < 0.01) {
				DelNod(newp);
				break;
			}
			for (int j = 0; j < 3; j++) {
				Nodes[newp].pt[j] = pnt[j] + len * tempNorm[j] * ratio;
			}
			ratio *= 0.5;

			int loc, srchtet;
			for (int j = 0; j < 3; j++) {
				srchtet = getP2T(pa);
				loc = locate_pnt(newp, srchtet);
				if (loc != 100)
					break;
			}

			if (loc < 1 || loc == 100 || ishulltet(srchtet)) {
				break;
			}
			else if (loc > 15) {//point in face
				int ord = (loc >> 4) - 1;
				int a = 0, b = 0, c = 0, d = 0;
				DNC(ord, a, b, c, d);
				b = Elems[srchtet].form[b];
				c = Elems[srchtet].form[c];
				d = Elems[srchtet].form[d];
				if (isBndTri(b, c, d))
					continue;
			}
			else if (4 <= loc && loc <= 14) {//point in edge
				int a = (loc & 12) >> 2, b = loc & 3;
				a = Elems[srchtet].form[a];
				b = Elems[srchtet].form[b];
				if (isBndEdg(a, b))
					continue;
			}
			std::vector<int> tempS = { srchtet };

			int ret = BW_insert_vertex(newp, tempS, 3);

			if (ret <= 0) {
				DelNod(newp);
			}
			else {
				addst++;
				newN.push_back(newp);
				if (recoverEdge(lostE, 0, 0) == 1) {
					//printf("success\n");
					return 1;
				}
			}
			break;
		}//while (1)
	}
	//}
	//printf("fail\n");
	return 0;
}

//dis[0],for <p1,p2,p3> is positive
//dis[1],for <p1,p2,p3> is negative
int DT::addinnerSteiner_Face(const int lostF, std::vector<int>& newN, double dis[2]) {
    double center[3]={};
    for (int i=0;i<3;++i) for (int k=0;k<3;++k) center[k]+=Nodes[SurTris[lostF].form[i]].pt[k]/3.0;
    return addInteriorFacePoints(lostF,newN,dis,center);
}

int DT::addInteriorFacePoints(const int lostF, std::vector<int>& newN, const double dis[2], const double base[3]) {
	int i, j, newp, srchtet, p1, p2, p3;
	double norm[3] = { 0 };
	p1 = SurTris[lostF].form[0];
	p2 = SurTris[lostF].form[1];
	p3 = SurTris[lostF].form[2];

	//cal normal
	calnormal(p1, p2, p3, norm);
	double normLen = lenvec(norm);
    if (!(normLen > 0)) return 0;
	for (i = 0; i < 3; i++)
		norm[i] /= normLen;

	double space = (Nodes[p1].space + Nodes[p2].space + Nodes[p3].space) / 3.0;//space


	for (i = 0; i < 2; i++) {
		double len = std::abs(dis[i]); // Both signed distances are lengths on opposite sides.

		//update normal
		if (i == 1) {
			for (j = 0; j < 3; j++)
				norm[j] *= -1;
		}

		double ratio = 0.5;
		newp = addNode();
		Nodes[newp].space = space;
		while (1) {
			if (ratio < 0.01 || len < 1e-16) {
				DelNod(newp);
				break;
			}
			for (j = 0; j < 3; j++) {
				Nodes[newp].pt[j] = base[j] + len * norm[j] * ratio;
			}
			ratio *= 0.5;

			int loc;
			for (int j = 0; j < 3; j++) {
				srchtet = getP2T(SurTris[lostF].form[j]);
				loc = locate_pnt(newp, srchtet);
				if (loc != 100)
					break;
			}

			if (loc < 1 || loc == 100 || ishulltet(srchtet))
				continue;
			else if (loc > 15) {//point in face
				int ord = (loc >> 4) - 1;
				int a = 0, b = 0, c = 0, d = 0;
				DNC(ord, a, b, c, d);
				b = Elems[srchtet].form[b];
				c = Elems[srchtet].form[c];
				d = Elems[srchtet].form[d];
				if (isBndTri(b, c, d))
					continue;
			}
			else if (4 <= loc && loc <= 14) {//point in edge
				int a = (loc & 12) >> 2, b = loc & 3;
				a = Elems[srchtet].form[a];
				b = Elems[srchtet].form[b];
				if (isBndEdg(a, b))
					continue;
			}
			std::vector<int> tempS = { srchtet };

			int ret = BW_insert_vertex(newp, tempS, 3);
			if (ret <= 0) {
                // BW failure is transactional; reuse this uninserted node for
                // a shorter step, and release it when this side is exhausted.
				continue;
			}
			else {
				addst++;
				newN.push_back(newp);
				if (isMeshEdge(SurTris[lostF].form[0], SurTris[lostF].form[1]) &&
					isMeshEdge(SurTris[lostF].form[1], SurTris[lostF].form[2]) &&
					isMeshEdge(SurTris[lostF].form[2], SurTris[lostF].form[0])) {
					if (recoverFace(lostF, 0) == 1)
						return 1;
				}
			}
			break;
		}
	}

	return 0;
}
//try to remove bnd edges steiner point
//#pragma optimize("",off)
int DT::removeEdgStiner(const int idx, int level) {
	if (level > 3) return 0;
	int iNod = EdgSteiner[idx].first;			//steiner point
	if (iNod == -1) return 1;
	int lostedg = EdgSteiner[idx].second;		//lost edges
	int p1 = SurEdgs[lostedg].iStart;
	int p2 = SurEdgs[lostedg].iEnd;
	int subedg0 = SurEdgs[lostedg].info;		//sub edge 0
	int subedg1 = SurEdgs[lostedg].info + 1;	//sub edge 1
	int mainfold = SurEdgs[lostedg].face.size();
	int a, b, c, d, i, j, k, m;
	double edgelen = distance(Nodes[p1].pt, Nodes[p2].pt);
	std::vector<int> newN;
	std::vector<double> volvec;

	newN.push_back(iNod);
	//If this edge has been recovered due to accuracy issues
	int tempE, trytmvst = 0;
	if (isMeshEdge(p1, p2, &tempE)) {
		std::vector<int> oldtet = { tempE };
		if (removeEdge(oldtet, isNod_in_Tet(p1, tempE), isNod_in_Tet(p2, tempE), 1000, -1) == 0) {
			trytmvst = 1;
		}
	}

	if (trytmvst == 1) {
		clearbndpnt(iNod);
		smooth_volume(iNod, true); //smooth_diff(iNod, 0);
	}
	else {
		if (mainfold == 1) {
			if (infolevel >= 2)
				meshLogger->debug("Can't remove single edge steiner now!");
			return 2;
		}
		/************************** Build subTri Hash ************************/
		TriHasher<int> subTri;
		for (i = 0; i < SurEdgs[subedg0].face.size(); i++) {
			SurTri* s = &SurTris[SurEdgs[subedg0].face[i]];
			subTri.add(s->form[0], s->form[1], s->form[2], SurEdgs[subedg0].face[i]);
			if (!isMeshFace(s->form[0], s->form[1], s->form[2])) {
				int ret = recoverFace(SurEdgs[subedg0].face[i], 0);
			}
		}
		for (i = 0; i < SurEdgs[subedg1].face.size(); i++) {
			SurTri* s = &SurTris[SurEdgs[subedg1].face[i]];
			subTri.add(s->form[0], s->form[1], s->form[2], SurEdgs[subedg1].face[i]);
			if (!isMeshFace(s->form[0], s->form[1], s->form[2])) {
				int ret = recoverFace(SurEdgs[subedg1].face[i], 0);
			}
		}
		/**************************** color sphere **************************/
		std::vector<int> sph;
		findSphere(iNod, sph);

		std::vector<int> sphColor(sph.size(), -1);
		std::map<int, int> mp2sph;
		for (i = 0; i < sph.size(); i++)
			mp2sph[sph[i]] = i;
		int color = 0;
		for (i = 0; i < sph.size(); i++) {
			if (sphColor[i] != -1)
				continue;
			std::queue<int> q;//for BFS
			q.push(i);
			sphColor[i] = color;
			while (!q.empty()) {
				int src = sph[q.front()];
				q.pop();
				for (m = 0; m < 4; m++) {
					int neig = getNeig(src, m);
					if (mp2sph.find(neig) == mp2sph.end())
						continue;//neig don't in sph
					if (sphColor[mp2sph[neig]] != -1)
						continue;
					DNC(m, a, b, c, d);
					b = Elems[src].form[b];
					c = Elems[src].form[c];
					d = Elems[src].form[d];
					if (subTri.find(b, c, d))
						continue;
					q.push(mp2sph[neig]);
					sphColor[mp2sph[neig]] = color;
				}
			}
			color++;
		}
		if (color != mainfold) {
			if (infolevel > 0)
				meshLogger->debug("Remove Edge steiner point {} error,classify sph fail!", idx);
			return 0;
		}
		/*********************** determine norm and new pt ***********************/
		std::vector<std::vector<int>> sphcla(color);
		std::vector<std::vector<double>> normal(color, std::vector<double>(3));
		std::vector<std::vector<double>> newPcoord(color, std::vector<double>(3));
		std::vector<std::vector<int>> bndElem(mainfold * 2, std::vector<int>(3, -1));
		for (i = 0; i < mainfold; i++) {
			//continue hull tet,ignore volume
			bool hullflag = false;
			//------------- Classifiy sphere ------------
			for (j = 0; j < sph.size(); j++) {
				if (sphColor[j] == i) {
					sphcla[i].push_back(sph[j]);
					if (ishulltet(sph[j]))
						hullflag = true;
				}
			}

			//--------------- get sub Tri and bnd ---------------
			int nsubsph = sphcla[i].size();
			std::vector<std::vector<int>> bndpt(nsubsph + 2, std::vector<int>(3, -1));
			std::map<int, int> PTV;//parent bnd tri have been visited?
			int SubTriord = 0;
			for (k = 0; k < 3; k++)
				normal[i][k] = 0;
			for (j = 0; j < nsubsph; j++) {
				int t0 = sphcla[i][j];
				for (m = 0; m < 4; m++) {
					int neig = getNeig(t0, m);
					DFC(m, a, b, c, d);
					b = Elems[t0].form[b];
					c = Elems[t0].form[c];
					d = Elems[t0].form[d];

					if (mp2sph.find(neig) == mp2sph.end()) {
						bndpt[j][0] = b;//store bnd point
						bndpt[j][1] = c;
						bndpt[j][2] = d;
						continue;//neig don't in sph
					}
					//In it's sphere
					if (sphColor[mp2sph[neig]] == i)
						continue;

					//determine ord of last two tri
					int tt = subTri.get(b, c, d);
					int Tid = SurTris[tt].parent;//parent tri's idx
					if (std::find(SurEdgs[lostedg].face.begin(), SurEdgs[lostedg].face.end(), Tid) == SurEdgs[lostedg].face.end()) {
						continue;
					}
					if (PTV.find(Tid) == PTV.end()) {
						PTV[Tid] = SubTriord++;
						//first bnd
						bndpt[nsubsph + PTV[Tid]][0] = b;
						bndpt[nsubsph + PTV[Tid]][1] = c;
						bndpt[nsubsph + PTV[Tid]][2] = d;
						for (k = 0; k < 3; k++) {
							if (bndpt[nsubsph + PTV[Tid]][k] == iNod) {
								for (int k1 = 0; k1 < 3; k1++) {
									if (SurTris[Tid].form[k1] != b && SurTris[Tid].form[k1] != c && SurTris[Tid].form[k1] != d) {
										bndpt[nsubsph + PTV[Tid]][k] = SurTris[Tid].form[k1];
										break;
									}
								}
								break;
							}
						}
						//normal [b,c,d]->a
						double f[3];
						calnormal(bndpt[nsubsph + PTV[Tid]][0], bndpt[nsubsph + PTV[Tid]][1], bndpt[nsubsph + PTV[Tid]][2], f);
						double flen = lenvec(f);
						for (k = 0; k < 3; k++)
							f[k] /= flen;
						for (k = 0; k < 3; k++)
							normal[i][k] += f[k];
					}
				}
			}
			if (PTV.size() != 2) {
				return 0;
			}

			//store bnd infomation
			for (j = 0; j < 2; j++)
				for (k = 0; k < 3; k++)
					bndElem[i * 2 + j][k] = bndpt[nsubsph + j][k];

			//Normalization
			for (k = 0; k < 3; k++)
				normal[i][k] /= 2.0;

			//----------------- determine the new point ----------------
			double oldp[3] = { Nodes[iNod].pt[0],Nodes[iNod].pt[1],Nodes[iNod].pt[2] };
			double newp[3] = { oldp[0] + edgelen * normal[i][0] ,oldp[1] + edgelen * normal[i][1] ,oldp[2] + edgelen * normal[i][2] };
			double ori, len = edgelen/10.0;
			double ShortestDistance = 1e-16;

			for (j = bndpt.size() - 1; j >= bndpt.size() - 2; j--) {
				ori = dt::GEOM_FUNC::orient3d(Nodes[bndpt[j][0]].pt, Nodes[bndpt[j][1]].pt, newp, Nodes[bndpt[j][2]].pt);
				if (ori <= 0) {
					std::swap(bndpt[j][1], bndpt[j][2]);
				}
			}

			//try to find a position,let all tet's volume is positive
			while (len > ShortestDistance) {
				bool allpositive = true;
				double newp[3] = { oldp[0] + len * normal[i][0] ,oldp[1] + len * normal[i][1] ,oldp[2] + len * normal[i][2] };
				for (j = bndpt.size() - 1; j >= 0; j--) {
					if (bndpt[j][0] == ghost || bndpt[j][1] == ghost || bndpt[j][2] == ghost)
						continue;
					ori = dt::GEOM_FUNC::orient3d(Nodes[bndpt[j][0]].pt, Nodes[bndpt[j][1]].pt, newp, Nodes[bndpt[j][2]].pt);
					if (ori <= 0) {
						if (j == bndpt.size() - 1 || j == bndpt.size() - 2) {
							len = ShortestDistance;
						}
						len *= 0.5;
						allpositive = false;
						break;
					}
				}
				if (allpositive) {
					for (k = 0; k < 3; k++)
						newPcoord[i][k] = newp[k];//store new position in normal for easily
					break;
				}
			}

			//can't find a position let all tet's volume is positive
			if (len < ShortestDistance) {
				//return 0;
				if (level < 1 && !hullflag) {
					int ret = -1;
					for (auto it : sph) {
						if (ishulltet(it)||isDelEle(it))
							continue;
						double vol = calVolume(it);
						if (vol < 1e-14 && vol >=0) {
							improve_step = true;
							improve_Metric = 3;
  							ret = removebadtet(it, -1);
                            if (ret == 0) ret = removebadtet_addPnt(it);
							improve_step = false;
						}
					}
					
					return removeEdgStiner(idx, level + 1);
				}

				for (k = 0; k < 3; k++)
					newPcoord[i][k] = oldp[k];//store new position in normal for easily
			}
		}

		/************************ creat new tet and neig info **************************/
		std::vector<int> newElemIdx(2 * mainfold);
		for (i = 0; i < mainfold; i++) {
			if (i == 0) {
				for (k = 0; k < 3; k++)
					Nodes[newN[i]].pt[k] = newPcoord[i][k];//obtain new position
			}//iNod used for first class
			else {
				newN.push_back(addNode(newPcoord[i][0], newPcoord[i][1], newPcoord[i][2], Nodes[iNod].space));
				//---------- update elems connect to iNod -----------
				for (j = 0; j < sphcla[i].size(); j++) {
					int t0 = sphcla[i][j];
					for (m = 0; m < 4; m++) {
						if (Elems[t0].form[m] == iNod) {
							Elems[t0].form[m] = newN[i];
							break;
						}
					}
				}
			}
			clearbndpnt(newN[i]);
			//reset point to tet
			setP2T(newN[i], sphcla[i][0]);
			//create new tet
			for (j = 0; j < 2; j++) {
				newElemIdx[i * 2 + j] = addElem(bndElem[i * 2 + j][0], bndElem[i * 2 + j][1], bndElem[i * 2 + j][2], newN[i]);
			}

			//connect 2 neig information
			for (j = 0; j < 2; j++) {
				int t0 = newElemIdx[i * 2 + j];
				for (k = 0; k < 3; k++) {//neig[0] connect to another class
					if (getNeig(t0, k) != -1)
						continue;
					DFC(k, a, b, c, d);
					b = Elems[t0].form[b];
					c = Elems[t0].form[c];
					d = Elems[t0].form[d];

					for (m = 0; m < sphcla[i].size(); m++) {
						int t1 = sphcla[i][m];
						int* form = Elems[t1].form;
						int ia = -1, ib = -1, ic = -1, id = -1;
						for (int k1 = 0; k1 < 4; k1++) {
							if (form[k1] == b)ib = k1;
							else if (form[k1] == c)ic = k1;
							else if (form[k1] == d)id = k1;
							else ia = k1;
						}
						if (ib == -1 || ic == -1 || id == -1)
							continue;
						//find neig
						bond(t0, k, t1, ia);
						break;
					}
				}
			}
		}
		//conect last one neig
		for (i = 0; i < mainfold * 2; i++) {
			int t0 = newElemIdx[i];
			for (k = 0; k < 4; k++) {//neig[0] connect to another class
				if (getNeig(t0, k) != -1)
					continue;
				DNC(k, a, b, c, d);
				b = Elems[t0].form[b];
				c = Elems[t0].form[c];
				d = Elems[t0].form[d];
				for (j = i + 1; j < mainfold * 2; j++) {
					int t1 = newElemIdx[j];
					int* form = Elems[t1].form;
					int ia = -1, ib = -1, ic = -1, id = -1;
					for (int k1 = 0; k1 < 4; k1++) {
						if (form[k1] == b)ib = k1;
						else if (form[k1] == c)ic = k1;
						else if (form[k1] == d)id = k1;
						else ia = k1;
					}
					if (ib == -1 || ic == -1 || id == -1)
						continue;
					//find neig
					bond(t0, k, t1, ia);
					break;
				}
			}
		}
	}
	/************************ update SurTris and SurEdgs **************************/
    // Two half-edges plus one spoke per incident face, including non-manifold edges.
    for (int child=subedg0; child<subedg0+mainfold+2; ++child) {
        setDelSurEdg(child);
        BndEdg.erase(SurEdgs[child].iStart,SurEdgs[child].iEnd);
    }
    SurEdgs[lostedg].info=1;
    BndEdg.add(p1,p2,lostedg);

	for (i = 0; i < SurEdgs[subedg0].face.size(); i++) {
		SurTri* s = &SurTris[SurEdgs[subedg0].face[i]];   //sub Tri
		setDelSurTri(SurEdgs[subedg0].face[i]);
		//SurTris[s->parent].info = 1;                      //set parent bnd Tri is recovered, Comment for findTriParent
		BndTri.erase(s->form[0], s->form[1], s->form[2]); //update BndTri Hash
	}
	for (i = 0; i < SurEdgs[subedg1].face.size(); i++) {
		SurTri* s = &SurTris[SurEdgs[subedg1].face[i]];   //sub Tri
		setDelSurTri(SurEdgs[subedg1].face[i]);
		//SurTris[s->parent].info = 1;                      //set parent bnd Tri is recovered
		BndTri.erase(s->form[0], s->form[1], s->form[2]); //update BndTri Hash
	}

	//checkMeshError();

	for (i = 0; i < newN.size(); i++) {
		if (!removePnt(newN[i])) {
			smooth_volume(newN[i], true); 
		}
	}

	//checkMeshError();
	//Recover BndTri
	for (i = 0; i < mainfold; i++) {
		int Fidx = SurEdgs[lostedg].face[i];
		BndTri.add(SurTris[Fidx].form[0], SurTris[Fidx].form[1], SurTris[Fidx].form[2], Fidx);
		if (!isMeshFace(SurTris[Fidx].form[0], SurTris[Fidx].form[1], SurTris[Fidx].form[2])) {
			SurTris[Fidx].info = -1;
			int ret = recoverFace(Fidx, 0);
			if (ret < 0) meshLogger->warn("May Error happen when remove steiner point in edge {}", idx);
		}
		SurTris[Fidx].info = 1;
		for (int j = 0; j < 3; j++) {
			int tempp1 = SurTris[Fidx].form[j];
			int tempp2 = SurTris[Fidx].form[(j + 1) % 3];
			int tempEdg = BndEdg.get(tempp1, tempp2);
			if (tempEdg == -1)
				continue;
			for (int k = 0; k < SurEdgs[tempEdg].face.size(); k++) {
				int subsubsubf = SurEdgs[tempEdg].face[k];
				for (int t = 0; t < 3; t++) {
					if (SurTris[subsubsubf].form[t] == iNod) {
						SurEdgs[tempEdg].face[k] = Fidx;
						break;
					}
				}
			}
		}
	}

	return 1;
}
//#pragma optimize("",on)

int DT::removeTriStiner(const int idx) {
	int iNod = TriSteiner[idx].first;				//steiner point
	int lostTriidx = TriSteiner[idx].second;		//lost edges idx
	SurTri* lostTri = &SurTris[lostTriidx];			//lost edges
	int p1 = lostTri->form[0];
	int p2 = lostTri->form[1];
	int p3 = lostTri->form[2];
	int a, b, c, d, i, j, k, m;

	//If this face has been recovered due to accuracy issues
	if (!isMeshFace(p1, p2, p3)) {
		/************************** Build subTri Hash ************************/
		TriHasher<int> subTri;
		for (i = lostTri->info; i < lostTri->info + 3; i++) {
			SurTri* s = &SurTris[i];
			subTri.add(s->form[0], s->form[1], s->form[2], i);
		}
		/**************************** color sphere **************************/
		std::vector<int> sph;
		findSphere(iNod, sph);
		std::vector<int> sphColor(sph.size(), -1);
		std::map<int, int> mp2sph;
		for (i = 0; i < sph.size(); i++)
			mp2sph[sph[i]] = i;
		int color = 0;
		for (i = 0; i < sph.size(); i++) {
			if (sphColor[i] != -1)
				continue;
			std::queue<int> q;//for BFS
			q.push(i);
			sphColor[i] = color;
			while (!q.empty()) {
				int src = sph[q.front()];
				q.pop();
				for (m = 0; m < 4; m++) {
					int neig = getNeig(src, m);
					if (mp2sph.find(neig) == mp2sph.end())
						continue;//neig don't in sph
					if (sphColor[mp2sph[neig]] != -1)
						continue;
					DNC(m, a, b, c, d);
					b = Elems[src].form[b];
					c = Elems[src].form[c];
					d = Elems[src].form[d];
					if (subTri.find(b, c, d))
						continue;
					q.push(mp2sph[neig]);
					sphColor[mp2sph[neig]] = color;
				}
			}
			color++;
		}
		if (color != 2) {
			if (infolevel > 0) meshLogger->debug("Remove Face steiner point {} error,classify sph fail!", idx);
			return 0;
		}
		/*********************** determine norm and new pt ***********************/
		std::vector<std::vector<int>> sphcla(2);
		double normal[2][3] = { 0 };
		int bndElem[2][3] = { 0 };
		for (i = 0; i < 2; i++) {
			//continue hull tet,ignore volume
			bool hullflag = false;
			//------------------ Classifiy sphere ---------------
			for (j = 0; j < sph.size(); j++) {
				if (sphColor[j] == i)
					sphcla[i].push_back(sph[j]);
				if (ishulltet(sph[j]))
					hullflag = true;
			}
			//--------------- get sub Tri and bnd ---------------
			int nsubsph = sphcla[i].size();
			std::vector<std::vector<int>> bndpt(nsubsph, std::vector<int>(3, -1));
			std::map<int, int> PTV;//parent bnd tri have been visited?
			int SubTriord = 0;
			for (j = 0; j < nsubsph; j++) {
				int t0 = sphcla[i][j];
				for (m = 0; m < 4; m++) {
					int neig = getNeig(t0, m);
					DNC(m, a, b, c, d);
					b = Elems[t0].form[b];
					c = Elems[t0].form[c];
					d = Elems[t0].form[d];
					if (mp2sph.find(neig) == mp2sph.end()) {
						bndpt[j][0] = b;//store bnd point
						bndpt[j][1] = c;
						bndpt[j][2] = d;
						continue;//neig don't in sph
					}
					if (sphColor[mp2sph[neig]] == i)
						continue;
					//normal [b,c,d]->a
					double f[3];
					calnormal(b, d, c, f);
					double flen = lenvec(f);
					for (k = 0; k < 3; k++)
						normal[i][k] += f[k] / flen;
					//determine ord of last two tri
				}
			}
			//----------------- determine the new point -----------------
			double f[3];
			calnormal(p1, p2, p3, f);
			if (dot(f, normal[i]) > 0) {
				bndElem[i][0] = p1; bndElem[i][1] = p3; bndElem[i][2] = p2;
			}
			else {
				bndElem[i][0] = p1; bndElem[i][1] = p2; bndElem[i][2] = p3;
			}
			//Normalization
			for (k = 0; k < 3; k++)
				normal[i][k] /= 3;
			double edgelen = std::min(std::min(distance(Nodes[p1].pt, Nodes[iNod].pt), distance(Nodes[p2].pt, Nodes[iNod].pt)), distance(Nodes[p3].pt, Nodes[iNod].pt));
			double oldp[3] = { Nodes[iNod].pt[0],Nodes[iNod].pt[1],Nodes[iNod].pt[2] };
			//try to find a position,let all tet's volume is positive
			double ShortestDistance = 1e-10;
			while (edgelen > ShortestDistance) {
				bool allpositive = true;
				double newp[3] = { oldp[0] + edgelen * normal[i][0] ,oldp[1] + edgelen * normal[i][1] ,oldp[2] + edgelen * normal[i][2] };
				for (j = 0; j < bndpt.size(); j++) {
					if (bndpt[j][0] == ghost || bndpt[j][1] == ghost || bndpt[j][2] == ghost)
						continue;
					double ori = dt::GEOM_FUNC::orient3d(newp, Nodes[bndpt[j][0]].pt, Nodes[bndpt[j][2]].pt, Nodes[bndpt[j][1]].pt);
					if (ori <= 0) {
						edgelen /= 2;
						allpositive = false;
						break;
					}
				}
				if (allpositive) {
					for (k = 0; k < 3; k++)
						normal[i][k] = newp[k];//stor new position in normal for easily
					break;
				}
			}
			//can't find a position let all tet's volume is positive
			if (edgelen < ShortestDistance) {
				for (k = 0; k < 3; k++)
					normal[i][k] = oldp[k];//stor new position in normal for easily
			}
		}
		/************************ creat new tet and neig info **************************/
		int newElemIdx[2] = { 0 };
		int newN[2] = { 0 };
		for (i = 0; i < 2; i++) {
			if (i == 0) {
				newN[i] = iNod;
				for (k = 0; k < 3; k++)
					Nodes[newN[i]].pt[k] = normal[i][k];//obtain new position
			}//iNod used for first class
			else {
				newN[i] = addNode(normal[i][0], normal[i][1], normal[i][2], Nodes[iNod].space);
				//---------- update elems connect to iNod -----------
				for (j = 0; j < sphcla[i].size(); j++) {
					int t0 = sphcla[i][j];
					for (m = 0; m < 4; m++) {
						if (Elems[t0].form[m] == iNod) {
							Elems[t0].form[m] = newN[i];
							break;
						}
					}
				}
			}
			clearbndpnt(newN[i]);
			//reset point to tet
			setP2T(newN[i], sphcla[i][0]);
			//create new tet
			newElemIdx[i] = addElem(newN[i], bndElem[i][0], bndElem[i][1], bndElem[i][2]);
			//connect 3 neig information
			int t0 = newElemIdx[i];
			for (k = 1; k < 4; k++) {//neig[0] connect to another class
				if (getNeig(t0, k) != -1)
					continue;
				DNC(k, a, b, c, d);
				b = Elems[t0].form[b];
				c = Elems[t0].form[c];
				d = Elems[t0].form[d];

				for (m = 0; m < sphcla[i].size(); m++) {
					int t1 = sphcla[i][m];
					int* form = Elems[t1].form;
					int ia = -1, ib = -1, ic = -1, id = -1;
					for (int k1 = 0; k1 < 4; k1++) {
						if (form[k1] == b)ib = k1;
						else if (form[k1] == c)ic = k1;
						else if (form[k1] == d)id = k1;
						else ia = k1;
					}
					if (ib == -1 || ic == -1 || id == -1)
						continue;
					//find neig
					bond(t0, k, t1, ia);
					break;
				}
			}
		}
		//conect last one neig
		bond(newElemIdx[0], 0, newElemIdx[1], 0);

		for (i = 0; i < 2; i++) {
			if (!removePnt(newN[i]))
				smooth_volume(newN[i], true); 
		}
	}
	else {
		clearbndpnt(iNod);
		//Recover BndTri
		BndTri.add(p1, p2, p3, lostTriidx);
		if (!removePnt(iNod))
			smooth_volume(iNod, true);
	}
	/************************ update SurTris and SurEdgs **************************/
	//Recover BndTri
		BndTri.add(p1, p2, p3, lostTriidx);

	//SurEdgs
	for (int i = 0; i < 3; i++) {
		int tempEdg = BndEdg.get(iNod, lostTri->form[i]);
		setDelSurEdg(tempEdg);								 //set edge destroyed
		BndEdg.erase(iNod, lostTri->form[i]);				 //update BndEdg Hash
	}
	for (i = lostTri->info; i < lostTri->info + 3; i++) {
		SurTri* s = &SurTris[i];
		setDelSurTri(i);									 //set Tri destroyed
		BndTri.erase(s->form[0], s->form[1], s->form[2]);	 //update BndTri Hash
	}
	for (int i = 0; i < 3; i++) {
		p1 = lostTri->form[i];
		p2 = lostTri->form[(i + 1) % 3];
		int edgid = BndEdg.get(p1, p2);
		for (int j = 0; j < SurEdgs[edgid].face.size(); j++) {
			if (lostTri->info <= SurEdgs[edgid].face[j] && SurEdgs[edgid].face[j] < lostTri->info + 3) {
				SurEdgs[edgid].face[j] = lostTriidx;
				break;
			}
		}
	}
	lostTri->info = 1;										 //set old tri recovered, comment for findTriParent

	return 1;
}

int DT::removePnt(const int iNod, int tryTime) {
	if (isbndpnt(iNod) || isDelNod(iNod))
		return 0;

	std::vector<int> sph;
	findSphere(iNod, sph);

	int minE_i, minE_j0, minE_j1;
	double minedge = DBL_MAX;
	std::map<double, std::array<int, 3>> mp;
	for (int j = 0; j < sph.size(); j++) {
		for (int k = 0; k < 6; k++) {
			if (isDelEle(sph[j]))
				continue;
			int ta = Elems[sph[j]].form[Egid[k][0]];
			int tb = Elems[sph[j]].form[Egid[k][1]];
			if (ta == iNod || tb == iNod) {
				double dis;
				if (ta == ghost || tb == ghost)  continue;//dis = DBL_MAX;
				else dis = distance2(Nodes[ta].pt, Nodes[tb].pt);
				if (mp.find(dis) == mp.end()) {
					minE_i = sph[j];
					minE_j0 = ta == iNod ? Egid[k][0] : Egid[k][1];
					minE_j1 = Egid[k][0] + Egid[k][1] - minE_j0;
					mp[dis] = { minE_i, minE_j0, minE_j1 };
				}
			}
		}
	}

	int i = 0, ret = 0;
	for (auto it : mp) {
		ret = destroyShortEdge(it.second[0], it.second[1], it.second[2]);
		if (i++ > tryTime || ret == 1) {
			break;
		}
	}

	if (ret == 1) {
		return 1;
	}

	if (sph.size() == 4) {
		int newe = flip41(sph, iNod);
		return 1;
	}

	return 0;
}

int DT::removeface(std::vector<int>& oldtet, int ia, int info, int thread_n) {
	int  a, b, c, d, e, ib = 0, ic = 0, id = 0, flipflag = 0, rmvEdg = -1, success = 0, ret = 0;
	double ori, * pa, * pb, * pc, * pd, * pe;
	std::vector<int> flat;
	oldtet.resize(2);

	DFC(ia, a, ib, ic, id);
	a = Elems[oldtet[0]].form[a];
	b = Elems[oldtet[0]].form[ib];
	c = Elems[oldtet[0]].form[ic];
	d = Elems[oldtet[0]].form[id];

	if (isBndTri(b, c, d)) {
		return 0;
	}
	oldtet[1] = getNeig(oldtet[0], ia);
    if (oldtet[1] < 0 || oldtet[1] >= static_cast<int>(Elems.size()) || isDelEle(oldtet[1])) {
        oldtet.resize(1);
        return 0;
    }
	e = getNeigOrd(oldtet[0], ia);
	e = Elems[oldtet[1]].form[e];

	if (ishulltet(oldtet[0])) {
		//for non-convex Tetrahedronlization
		std::vector<int> sph;
		findSphere(a, sph);
		for (int i = 0; i < sph.size(); i++)
			if ((isNod_in_Tet(a, sph[i]) != -1 && isNod_in_Tet(e, sph[i]) != -1 && isNod_in_Tet(b, sph[i]) != -1 && isNod_in_Tet(c, sph[i]) != -1) ||
				(isNod_in_Tet(a, sph[i]) != -1 && isNod_in_Tet(e, sph[i]) != -1 && isNod_in_Tet(c, sph[i]) != -1 && isNod_in_Tet(d, sph[i]) != -1) ||
				(isNod_in_Tet(a, sph[i]) != -1 && isNod_in_Tet(e, sph[i]) != -1 && isNod_in_Tet(d, sph[i]) != -1 && isNod_in_Tet(b, sph[i]) != -1)
				) {
				return 0;
			}

		DFC(3, d, a, b, c);
		pa = Nodes[Elems[oldtet[0]].form[a]].pt;
		pb = Nodes[Elems[oldtet[0]].form[b]].pt;
		pc = Nodes[Elems[oldtet[0]].form[c]].pt;
		pe = Nodes[e].pt;
		ori = dt::GEOM_FUNC::orient3d(pa, pb, pc, pe);
		if (ori >= 0) {
			return 0;
		}
		
		flip23(oldtet, ia, thread_n);
		return 1;
	}

	if (thread_n != -1) {
		if (!tryOccupying(Nodes[e].occupying, thread_n)) {
			return -1;
		}
	}

	pa = Nodes[a].pt;
	pb = Nodes[b].pt;
	pc = Nodes[c].pt;
	pd = Nodes[d].pt;
	pe = Nodes[e].pt;

	ori = dt::GEOM_FUNC::orient3d(pa, pb, pd, pe);
	rmvEdg = 0;
	if (ori < 0) {
		ori = dt::GEOM_FUNC::orient3d(pa, pd, pc, pe);
		rmvEdg = 1;
		if (ori < 0) {
			ori = dt::GEOM_FUNC::orient3d(pa, pc, pb, pe);
			rmvEdg = 2;
			if (ori < 0) {
				// Found a 2-to-3 flip.
				flipflag = 1;
			}
		}
	}

	if (flipflag == 1) {
		// A 2-to-3 flip is found.
		if (improve_step) {
			//in improvement step,we should update quality
			double minq = std::min(Elems[oldtet[0]].q, Elems[oldtet[1]].q);
			double q1 = 0, q2 = 0, q3 = 0;
			if (improve_Metric != 8) {
				double AniMetric1[6] = { 0 }, AniMetric2[6] = { 0 }, AniMetric3[6] = { 0 };
				if (improve_step && AniSol.size() != 0) {
					getmm(a, b, d, e, AniMetric1);
					getmm(a, d, c, e, AniMetric2);
					getmm(a, c, b, e, AniMetric3);
				}
				q1 = tetquality(pa, pb, pd, pe, AniMetric1, improve_Metric);
				q2 = tetquality(pa, pd, pc, pe, AniMetric2, improve_Metric);
				q3 = tetquality(pa, pc, pb, pe, AniMetric3, improve_Metric);

			}
			else {
				q1 = orthogonal(pa, pb, pd, pe,
					getoppoP(oldtet[1], isNod_in_Tet(c, oldtet[1])) == ghost ? NULL : Nodes[getoppoP(oldtet[1], isNod_in_Tet(c, oldtet[1]))].pt, pc, pc,
					getoppoP(oldtet[0], ic) == ghost ? NULL : Nodes[getoppoP(oldtet[0], ic)].pt);
				q2 = orthogonal(pa, pd, pc, pe,
					getoppoP(oldtet[1], isNod_in_Tet(b, oldtet[1])) == ghost ? NULL : Nodes[getoppoP(oldtet[1], isNod_in_Tet(b, oldtet[1]))].pt, pb, pb,
					getoppoP(oldtet[0], ib) == ghost ? NULL : Nodes[getoppoP(oldtet[0], ib)].pt);
				q3 = orthogonal(pa, pc, pb, pe,
					getoppoP(oldtet[1], isNod_in_Tet(d, oldtet[1])) == ghost ? NULL : Nodes[getoppoP(oldtet[1], isNod_in_Tet(d, oldtet[1]))].pt, pd, pd,
					getoppoP(oldtet[0], id) == ghost ? NULL : Nodes[getoppoP(oldtet[0], id)].pt);
			}
			if (q1 < minq || q2 < minq || q3 < minq) {
				if (thread_n != -1) clearOccupying(Nodes[e].occupying, thread_n);
				return 0;//if create worse tet
			}
		}
		flip23(oldtet, ia, thread_n);
		if (thread_n != -1) {
			clearOccupying(Nodes[a].occupying, thread_n);
			clearOccupying(Nodes[b].occupying, thread_n);
			clearOccupying(Nodes[c].occupying, thread_n);
			clearOccupying(Nodes[d].occupying, thread_n);
			clearOccupying(Nodes[e].occupying, thread_n);
		}
		return 1;
	}

	if (thread_n != -1)
		clearOccupying(Nodes[e].occupying, thread_n);

	if (flipflag != 1 && !improve_step) {
		if (rmvEdg == 0) {
			ret = removeEdge(oldtet, ib, id, info, thread_n);
		}
		else if (rmvEdg == 1) {
			ret = removeEdge(oldtet, id, ic, info, thread_n);
		}
		else {
			ret = removeEdge(oldtet, ic, ib, info, thread_n);
		}
		if (ret == 1) {
			return 1;
		}
	}
	return -(rmvEdg + 1);
}

int DT::removeEdge(std::vector<int>& oldtet, int ia, int ib, int info, int thread_n) {
    auto& flipHistory = thread_n < 0 ? flipnmRecll : parallelFlipHistory[thread_n];
    auto& flipCount = thread_n < 0 ? tempfliptime : parallelFlipCount[thread_n];

	int ie = oldtet[0];
	if (info < 0 || isDelEle(ie))
		return 0;
	double minq = 1;
	int pa = Elems[ie].form[ia];
	int pb = Elems[ie].form[ib];

	if (isBndEdg(pa, pb)) {//can't remove Boundary Edge
		return 0;
	}
	std::vector<int> shell_point;
	//find shell around edge pa_pb,store at old tet
	if (findShell(ie, ia, ib, oldtet, shell_point, thread_n) == -1)
		return -1;

    struct ShellRelease {
        DT* mesh; const std::vector<int>& points; int owner;
        ~ShellRelease() {
            if (owner >= 0) for (size_t i = 2; i < points.size(); ++i)
                mesh->clearOccupying(mesh->Nodes[points[i]].occupying, owner);
        }
    } shellRelease{this, shell_point, thread_n};

	// Mesh may contains inverted tetrahedra.
	if (oldtet.size() < 3) {
		//checkMeshError();
		//printf("%d %d\n", pa, pb);
		meshLogger->error("Opt remove edge {} {} connect to two tet.", pa, pb);
		throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
	}

	int maxlevel = info;
	if (improve_step) {
		for (int i = 0; i < oldtet.size(); i++) {
			if (Elems[oldtet[i]].q == -1)
				updateQuality(oldtet[i]);
			minq = std::min(minq, Elems[oldtet[i]].q);
		}
		minq = std::max(minq, 0.0);
		flipHistory.clear();
	}

	//set Elem's info mean this tet is in star.
	for (int i = 0; i < oldtet.size(); i++)
		Elems[oldtet[i]].info++;

	int ret = flipnm(oldtet, ia, ib, 0, maxlevel, minq, thread_n);

	if (ret == 2) {
		if (thread_n != -1) {
			clearOccupying(Nodes[pa].occupying, thread_n);
			clearOccupying(Nodes[pb].occupying, thread_n);
			for (int i = 0; i < shell_point.size(); i++)
				clearOccupying(Nodes[shell_point[i]].occupying, thread_n);
		}
		flipCount = 0;
		return 1;//remove success
	}

	//clear Elem's info mean this tet is in star.
	for (int i = 0; i < oldtet.size(); i++) {
		Elems[oldtet[i]].info--;
	}
	
	if (thread_n != -1)
		for (int i = 2; i < shell_point.size(); i++)
			clearOccupying(Nodes[shell_point[i]].occupying, thread_n);

	flipCount = 0;
	return 0;
}

/*
* oldtet:all tet around edge
* ia,ib:Edge's point of oldtet[0]'s order
* level: deep of flipnm
*/
int DT::flipnm(std::vector<int>& oldtet, const int ia, const int ib, int level, int maxlevel, double minq, int thread_n) {
    auto& flipHistory = thread_n < 0 ? flipnmRecll : parallelFlipHistory[thread_n];
    auto& flipCount = thread_n < 0 ? tempfliptime : parallelFlipCount[thread_n];

	int pa, pb, pc, pd, pe, pf, a, b, c = 0, d = 0, e = 0, n, nn, i, j, k, hullflag, flipflag, newtet, flatflag = 0, flatedge;
	double ori, qa, qb, qc;
	int oldtet0P[4] = { Elems[oldtet[0]].form[0] ,Elems[oldtet[0]].form[1] ,Elems[oldtet[0]].form[2] ,Elems[oldtet[0]].form[3] };
	n = oldtet.size();

	if (flipCount > maxfliptime)
		return n;

	if (n == 3) {//flip32
		hullflag = flipflag = 0;
		pa = Elems[oldtet[0]].form[ia];
		pb = Elems[oldtet[0]].form[ib];
		//if have hull tet,store at second and third oldtet
		if (pa != ghost && pb != ghost) {
			while (ishulltet(oldtet[0])) {
				std::swap(oldtet[0], oldtet[1]);
				std::swap(oldtet[0], oldtet[2]);
				hullflag = true;
			}
			if (ishulltet(oldtet[1])) {
				hullflag = true;
				for (j = 0; j < 4; j++) {//update a,b
					if (Elems[oldtet[0]].form[j] == pa)a = j;
					else if (Elems[oldtet[0]].form[j] == pb)b = j;
				}
			}
			else {
				a = ia;
				b = ib;
			}
		}
		else {
			a = ia;
			b = ib;
			hullflag = true;
		}

		//quickly found pa,pb,pc,pd,pe
		DDNC(c, d, a, b);

		pc = Elems[oldtet[0]].form[c];
		pd = Elems[oldtet[0]].form[d];
		pe = Elems[oldtet[1]].form[getNeigOrd(oldtet[0], c)];

		if (!hullflag) {
			ori = dt::GEOM_FUNC::orient3d(Nodes[pd].pt, Nodes[pc].pt, Nodes[pe].pt, Nodes[pa].pt);
			if (ori < 0) {
				ori = dt::GEOM_FUNC::orient3d(Nodes[pc].pt, Nodes[pd].pt, Nodes[pe].pt, Nodes[pb].pt);
				if (ori < 0) flipflag = true;
			}
		}
		else {
			flipflag = true;
		}

		if (flipflag) {
			if (flipintersectcheck(2, pc, pd, pe, pb, pa)) {
				flipflag = false;
			}
		}
		if (flipflag) {
			if (improve_step) {
				if (!hullflag) {

					double qa = 0, qb = 0;
					if (improve_Metric != 8) {
						double AniMetric1[6] = { 0 }, AniMetric2[6] = { 0 };
						if (AniSol.size() != 0) {
							getmm(pd, pc, pe, pa, AniMetric1);
							getmm(pc, pd, pe, pb, AniMetric2);
						}

						qa = tetquality(Nodes[pd].pt, Nodes[pc].pt, Nodes[pe].pt, Nodes[pa].pt, AniMetric1, improve_Metric);
						qb = tetquality(Nodes[pc].pt, Nodes[pd].pt, Nodes[pe].pt, Nodes[pb].pt, AniMetric2, improve_Metric);
					}
					else {
						int tetde = -1;
						if (isNod_in_Tet(pd, oldtet[1]) != -1) tetde = oldtet[1];
						else  tetde = oldtet[2];
						int tetce = oldtet[1] + oldtet[2] - tetde;

						qa = orthogonal(Nodes[pd].pt, Nodes[pc].pt, Nodes[pe].pt, Nodes[pa].pt,
							getoppoP(tetce, isNod_in_Tet(pb, tetce)) == ghost ? NULL : Nodes[getoppoP(tetce, isNod_in_Tet(pb, tetce))].pt,
							getoppoP(tetde, isNod_in_Tet(pb, tetde)) == ghost ? NULL : Nodes[getoppoP(tetde, isNod_in_Tet(pb, tetde))].pt,
							getoppoP(oldtet[0], isNod_in_Tet(pb, oldtet[0])) == ghost ? NULL : Nodes[getoppoP(oldtet[0], isNod_in_Tet(pb, oldtet[0]))].pt, Nodes[pb].pt);
						qb = orthogonal(Nodes[pc].pt, Nodes[pd].pt, Nodes[pe].pt, Nodes[pb].pt, 
							getoppoP(tetde, isNod_in_Tet(pa, tetde)) == ghost ? NULL : Nodes[getoppoP(tetde, isNod_in_Tet(pa, tetde))].pt,
							getoppoP(tetce, isNod_in_Tet(pa, tetce)) == ghost ? NULL : Nodes[getoppoP(tetce, isNod_in_Tet(pa, tetce))].pt,
							getoppoP(oldtet[0], isNod_in_Tet(pa, oldtet[0])) == ghost ? NULL : Nodes[getoppoP(oldtet[0], isNod_in_Tet(pa, oldtet[0]))].pt, Nodes[pa].pt);
					}

					if (qa <= minq + 1e-15 || qb <= minq + 1e-15) {
						//backflag = true;
						return 3;//if create worse tet
					}
				}
				flipHistory.push_back({ pd,pc,pe });
			}

			flip32(oldtet, a, b, thread_n);

			if (hullflag) {
				for (int j = 0; j < 2; j++) {
					for (int i = 0; i < 3; i++) {
						if (Elems[oldtet[j]].form[i] == ghost) {
							int p[4] = { Elems[oldtet[j]].form[(i + 1) % 4] ,Elems[oldtet[j]].form[(i + 2) % 4] ,Elems[oldtet[j]].form[(i + 3) % 4] ,Elems[oldtet[j]].form[i] };
							matchtet(p, oldtet[j]);
							break;
						}
					}
				}
			}
			flipCount++;

			oldtet.resize(2);
			return 2;//success
		}
	}//oldtet.size() == 3
	else {//oldtet.size() >3
		pa = Elems[oldtet[0]].form[ia];
		pb = Elems[oldtet[0]].form[ib];
		for (i = n - 1; i >= 0; i--) {
			a = b = -1;

			for (j = 0; j < 4; j++) {
				if (Elems[oldtet[i]].form[j] == pa) a = j;
				else if (Elems[oldtet[i]].form[j] == pb) b = j;
			}
			DDNC(c, d, a, b);
			int lefttet = getNeig(oldtet[i], d);
			if (Elems[oldtet[i]].info > 1 || Elems[lefttet].info > 1) {//if there are two Stars involved.
				continue;
			}
			e = getNeigOrd(oldtet[i], d);

			pc = Elems[oldtet[i]].form[c];
			pd = Elems[oldtet[i]].form[d];
			pe = Elems[lefttet].form[e];

			if (pd == ghost || pe == ghost)//[a,b,d] is hull face
				continue;
			//target face is [a,b,c]. e left d right
			if (isBndTri(pa, pb, pc))//can't remove bnd Tri
				continue;
			flipflag = 0;
			hullflag = (pc == ghost);

			if (hullflag == 0) {
				ori = dt::GEOM_FUNC::orient3d(Nodes[pa].pt, Nodes[pc].pt, Nodes[pd].pt, Nodes[pe].pt);//t[a,c,d,e]
				if (ori < 0) {
					ori = dt::GEOM_FUNC::orient3d(Nodes[pc].pt, Nodes[pb].pt, Nodes[pd].pt, Nodes[pe].pt);//t[c,b,d,e]
					if (ori < 0) {
						ori = dt::GEOM_FUNC::orient3d(Nodes[pb].pt, Nodes[pa].pt, Nodes[pd].pt, Nodes[pe].pt);//t[b,a,d,e]
						if (ori < 0) {
							flipflag = 1;
						}
						else if (ori == 0) {
							// [a,b] is flat.
							if (n == 4) {
								// The "flat" tet can be removed immediately by a 3-to-2 flip.
								flipflag = 1;
							}
						}
					}
				}
				if (!flipflag) {
					flatflag = 1;
				}
			}//hullflag == 0
			else {// 'c' is dummypoint.
				if (n == 4) {
					// Let the vertex opposite to 'c' is 'f'.
					// A 4-to-4 flip is possible if the two tets [d,e,f,a] and [e,d,f,b] are valid tets.
					pf = getoppoP(oldtet[i], c);
					ori = dt::GEOM_FUNC::orient3d(Nodes[pd].pt, Nodes[pe].pt, Nodes[pf].pt, Nodes[pa].pt);//t[d,e,f,a]
					if (ori < 0) {
						ori = dt::GEOM_FUNC::orient3d(Nodes[pe].pt, Nodes[pd].pt, Nodes[pf].pt, Nodes[pb].pt);//t[d,e,f,b]
						if (ori < 0) {
							flipflag = 1;
							ori = 0; // Signal as a 4-to-4 flip (like a co-planar case).
						}
					}
				}
			}

			if (flipflag) {
				if (flipintersectcheck(1, pa, pb, pc, pd, pe))
					flipflag = false;
			}

			if (flipflag && improve_step) {

				double qa = 0, qb = 0;
				if (improve_Metric != 8) {
					double AniMetric1[6] = { 0 }, AniMetric2[6] = { 0 };
					if (AniSol.size() != 0) {
						getmm(pc, pb, pd, pe, AniMetric1);
						getmm(pa, pc, pd, pe, AniMetric2);
					}

					qa = tetquality(Nodes[pc].pt, Nodes[pb].pt, Nodes[pd].pt, Nodes[pe].pt, AniMetric1, improve_Metric);
					qb = tetquality(Nodes[pa].pt, Nodes[pc].pt, Nodes[pd].pt, Nodes[pe].pt, AniMetric2, improve_Metric);
				}
				else {
					qa = orthogonal(Nodes[pc].pt, Nodes[pb].pt, Nodes[pd].pt, Nodes[pe].pt,
						Nodes[pa].pt, Nodes[pa].pt,
						getoppoP(lefttet, isNod_in_Tet(pa, lefttet)) == ghost ? NULL : Nodes[getoppoP(lefttet, isNod_in_Tet(pa, lefttet))].pt,
						getoppoP(oldtet[i], isNod_in_Tet(pa, oldtet[i])) == ghost ? NULL : Nodes[getoppoP(oldtet[i], isNod_in_Tet(pa, oldtet[i]))].pt);
					qb = orthogonal(Nodes[pa].pt, Nodes[pc].pt, Nodes[pd].pt, Nodes[pe].pt,
						Nodes[pb].pt, Nodes[pb].pt,
						getoppoP(lefttet, isNod_in_Tet(pb, lefttet)) == ghost ? NULL : Nodes[getoppoP(lefttet, isNod_in_Tet(pb, lefttet))].pt,
						getoppoP(oldtet[i], isNod_in_Tet(pb, oldtet[i])) == ghost ? NULL : Nodes[getoppoP(oldtet[i], isNod_in_Tet(pb, oldtet[i]))].pt);
				}

				if (qa <= minq + 1e-15 || qb <= minq + 1e-15) {
					flipflag = false;
				}
				else {
					flipHistory.push_back({ pd,pe,-1 });
				}	
			}
			if (flipflag) {
				//could flip
				nn = INT16_MAX;//reset nn
				std::vector<int> flipvec(3);
				flipvec[0] = lefttet;//find letf tet
				flipvec[1] = oldtet[i];

				flip23(flipvec, e, thread_n);
				flipCount++;

				//find new tet belong to new tet
				int tet_in_star = -1;
				for (j = 0; j < 3; j++) {
					for (k = 0; k < 4; k++) {
						if (Elems[flipvec[j]].form[k] == pc)
							break;
					}
					if (k == 4) {
						tet_in_star = j;
						break;
					}
				}
				Elems[flipvec[tet_in_star]].info++;//set this tet in star;
				//rebuild oldtet
				k = 0;
				for (j = 0; j < n; j++) {
					if (j == i) { oldtet[k++] = flipvec[tet_in_star]; }
					else if (oldtet[j] == lefttet) { continue; }
					else { oldtet[k++] = oldtet[j]; }
				}

				oldtet.resize(n - 1);
				//find pa,pb's idx in oldtet[0].form
				int newia = -1, newib = -1;
				for (j = 0; j < 4; j++) {
					if (Elems[oldtet[0]].form[j] == pa)newia = j;
					if (Elems[oldtet[0]].form[j] == pb)newib = j;
				}
				//next flip

				nn = flipnm(oldtet, newia, newib, level, maxlevel, minq, thread_n);

				if (nn == 2) {
					//flip success
					return nn;
				}
				else if (ori == 0) {
					//in mesh improvement,we will backtrack if the quality worse
					int havpepd = -1;
					std::vector<int> backtrack;
					//find a tet don't have pc
					for (j = 0; j < 3; j++) {
						if (j == tet_in_star)
							continue;//old tet_in_star must destroy
						d = -1, e = -1;
						for (k = 0; k < 4; k++) {
							if (Elems[flipvec[j]].form[k] == pd) d = k;
							else if (Elems[flipvec[j]].form[k] == pe) e = k;
							if (d != -1 && e != -1) {
								havpepd = flipvec[j];
								break;
							}
						}
						if (havpepd != -1)
							break;
					}
					if (havpepd == -1) {//old flipvec is useless
						if (!isMeshEdge(pd, pe, &havpepd)) {
							meshLogger->error("Pd and Pe is missing in backtrack flipnm!");
							throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
						}
						d = -1, e = -1;
						for (k = 0; k < 4; k++) {
							if (Elems[havpepd].form[k] == pd) d = k;
							else if (Elems[havpepd].form[k] == pe) e = k;
						}
					}
					//find a tet e pd and pe success,and find all backtrack
					std::vector<int> shell_point;
					findShell(havpepd, d, e, backtrack, shell_point);
					if (backtrack.size() != 3) {
						meshLogger->error("Backtrack flipnm error!");
						throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
					}

					//if have hull tet,store at second and third oldtet
					while (ishulltet(backtrack[0])) {
						std::swap(backtrack[0], backtrack[1]);
						std::swap(backtrack[0], backtrack[2]);
					}
					if (ishulltet(backtrack[1])) {
						for (j = 0; j < 4; j++) {//update a,b
							if (Elems[backtrack[0]].form[j] == pd) d = j;
							else if (Elems[backtrack[0]].form[j] == pe) e = j;
						}
					}

					flip32(backtrack, d, e, thread_n);//backtrack

					Elems[backtrack[0]].info++;//set this tet in star;
					Elems[backtrack[1]].info++;//set this tet in star;
					//re find oldtet from backtrack[0]
					if (i <= 1) {
						int oldtet0 = findtet(oldtet0P, backtrack);//find and correct direction of oldtet[0]
						findShell(oldtet0, ia, ib, oldtet, shell_point);
						return oldtet.size();
					}
					else {
						findShell(oldtet[0], ia, ib, oldtet, shell_point);
					}
				}
				else if (improve_step) {
					while (!flipHistory.empty())
					{
						auto it = flipHistory.back();
						flipHistory.pop_back();
						int te = -1;
						if (it[2] == -1) {
							if (isMeshEdge(it[0], it[1], &te)) {
								int ia = isNod_in_Tet(it[0], te);
								int ib = isNod_in_Tet(it[1], te);
								std::vector<int> shell, shellp;
								findShell(te, ia, ib, shell, shellp);
								flip32(shell, ia, ib, thread_n);
							}
							if (it[0] == pd && it[1] == pe) {

								break;
							}
						}
						else {
							if (isMeshFace(it[0], it[1],it[2], & te)) {
								int ia = isNod_in_Tet(it[0], te);
								int ib = isNod_in_Tet(it[1], te);
								int ic = isNod_in_Tet(it[2], te);
								int id = 0 + 1 + 2 + 3 - ia - ib - ic;
								int neig = getNeig(te, id);
								std::vector<int> shell = { te,neig };
								flip23(shell, id, thread_n);
							}
						}
					}
					if (isMeshEdge(pa, pb, &newtet)) {//ab has exists
						a = b = -1;
						for (k = 0; k < 4; k++) {
							if (Elems[newtet].form[k] == pa) a = k;
							else if (Elems[newtet].form[k] == pb) b = k;
						}
						std::vector<int> shellp;
						findShell(newtet, b, a, oldtet, shellp);//find new shell ab
						for (k = 0; k < oldtet.size(); k++) {
							if (Elems[oldtet[k]].info == 0) {//set thit tet in star
								Elems[oldtet[k]].info++;
							}
						}
					}
					return oldtet.size();
				}
				else {
					return oldtet.size();
				}
			}
		}//for (i = 0; i < n; i++)

		//The size of oldtet don't reduce
		if (flatflag == 1) {
			//remove ac/bc
			if (level < maxlevel) {
				for (i = 0; i < n; i++) {
					int falttet = oldtet[i];
					a = b = -1;
					for (j = 0; j < 4; j++) {
						if (Elems[falttet].form[j] == pa) a = j;
						else if (Elems[falttet].form[j] == pb) b = j;
					}
					DDNC(c, d, a, b);
					int lefttet = getNeig(falttet, d);
					if (Elems[falttet].info > 1 || Elems[lefttet].info > 1) {
						//if there are two Stars involved.
						continue;
					}
					e = getNeigOrd(oldtet[i], d);
					pc = Elems[falttet].form[c];
					pd = Elems[falttet].form[d];
					pe = Elems[lefttet].form[e];

					if (pc == ghost || pd == ghost || pe == ghost)//[a,b,d] is hull face
						continue;

					flatedge = -1;

					double ori1 = dt::GEOM_FUNC::orient3d(Nodes[pc].pt, Nodes[pb].pt, Nodes[pd].pt, Nodes[pe].pt);//t[c,b,d,e];
					double ori2 = dt::GEOM_FUNC::orient3d(Nodes[pa].pt, Nodes[pc].pt, Nodes[pd].pt, Nodes[pe].pt);//t[a,c,d,e]

					if (ori1 < 0 && ori2 >= 0) {
						flatedge = a;
					}
					else if (ori2 < 0 && ori1 >= 0) {
						flatedge = b;
					}
					else if (ori2 >= 0 && ori1 >= 0) {
						//printf("%.16lf %.16lf \n", ori1, ori2);
						if (isBndEdg(pc, Elems[falttet].form[b]) || ori2 > ori1) {
							flatedge = a;
						}
						else if (isBndEdg(pc, Elems[falttet].form[a]) || ori1 > ori2) {
							flatedge = b;
						}
					}

					if (flatedge == -1)
						continue;

					//Comment those code for test boundary recover
					if (isBndEdg(pc, Elems[falttet].form[flatedge]))
						continue;

					std::vector<int> flatvec;
					std::vector<int> shell_point;
					findShell(falttet, flatedge, c, flatvec, shell_point);
					int overlaps = 0;
					for (k = 0; k < flatvec.size(); k++)
						overlaps += Elems[flatvec[k]].info;//set thit tet in star
					if (overlaps > 2)
						continue;
					for (k = 0; k < flatvec.size(); k++) {
						Elems[flatvec[k]].info++;//set thit tet in star
					}
					
					double tempminq = 1;
					if (improve_step) {
						for (int j = 0; j < flatvec.size(); j++) {
							tempminq = std::min(tempminq, Elems[flatvec[j]].q);
						}
					}
					nn = flipnm(flatvec, flatedge, c, level + 1, maxlevel, tempminq);
					if (nn == 2) {
						if (isMeshEdge(pa, pb, &newtet)) {//ab has exists
							a = b = -1;
							for (k = 0; k < 4; k++) {
								if (Elems[newtet].form[k] == pa) a = k;
								else if (Elems[newtet].form[k] == pb) b = k;
							}
							std::vector<int> shell_point;
							findShell(newtet, b, a, oldtet, shell_point);//find new shell ab
							for (k = 0; k < oldtet.size(); k++) {
								if (Elems[oldtet[k]].info == 0) {//set thit tet in star
									Elems[oldtet[k]].info++;
								}
							}

							return flipnm(oldtet, b, a, level, maxlevel, minq);
						}//if (isMeshEdge(pa, pb, &newtet)) {//ab has exists
					}
					else {
						for (k = 0; k < nn; k++) {
							Elems[flatvec[k]].info--;//set this tet in star
						}
					}
				}
			}
		}
	}
	return oldtet.size();
}

/*
* oldtet:three tet store in
* ia,ib:Edge of oldtet[0]
* s.t. only flip,if it could be fliped should be deterine in last step
*/
int DT::flip32(std::vector<int>& oldtet, int ia, int ib, int thread_n) {
	int a[3], b[3], c[2], d[2], e[2], i, pa, pb, pc, pd, pe, adj1, adj2, newa, newb;

	pa = Elems[oldtet[0]].form[ia];
	pb = Elems[oldtet[0]].form[ib];
	//determine a,b,c,d,e order in three tet
	a[0] = ia;
	b[0] = ib;
	DDNC(c[0], d[0], a[0], b[0]);
	pc = Elems[oldtet[0]].form[c[0]];
	pd = Elems[oldtet[0]].form[d[0]];

	//check order and swap
	adj1 = getNeig(oldtet[0], c[0]);
	adj2 = getNeig(oldtet[0], d[0]);
	if (adj1 == oldtet[2] && adj2 == oldtet[1])
		std::swap(oldtet[1], oldtet[2]);

	for (int i = 0; i < 4; i++) {
		if (Elems[oldtet[1]].form[i] == pa) a[1] = i;
		else if (Elems[oldtet[1]].form[i] == pb) b[1] = i;
		else if (Elems[oldtet[1]].form[i] == pd) d[1] = i;
		else {
			pe = Elems[oldtet[1]].form[i];
			e[0] = i;
		}
	}

	for (int i = 0; i < 4; i++) {
		if (Elems[oldtet[2]].form[i] == pa) a[2] = i;
		else if (Elems[oldtet[2]].form[i] == pb) b[2] = i;
		else if (Elems[oldtet[2]].form[i] == pc) c[1] = i;
		else {
			e[1] = i;
		}
	}

	//create new tet,pe should at forth pos,it might be ghost
	if (thread_n == -1) {
		newa = addElem(pa, pc, pd, pe);
		newb = addElem(pb, pd, pc, pe);
	}
	else {
		newa = addElem(pa, pc, pd, pe, thread_n);
		newb = addElem(pb, pd, pc, pe, thread_n);
	}
 
	//bond adjacent
	bond(newb, 0, newa, 0);
	bond(newa, 1, getNeig(oldtet[1], b[1]), getNeigOrd(oldtet[1], b[1]));
	bond(newa, 2, getNeig(oldtet[2], b[2]), getNeigOrd(oldtet[2], b[2]));
	bond(newa, 3, getNeig(oldtet[0], b[0]), getNeigOrd(oldtet[0], b[0]));

	bond(newb, 1, getNeig(oldtet[2], a[2]), getNeigOrd(oldtet[2], a[2]));
	bond(newb, 2, getNeig(oldtet[1], a[1]), getNeigOrd(oldtet[1], a[1]));
	bond(newb, 3, getNeig(oldtet[0], a[0]), getNeigOrd(oldtet[0], a[0]));
	//set point to tet
	setP2T(pa, newa);
	setP2T(pb, newb);
	setP2T(pc, newb);
	setP2T(pd, newb);
	setP2T(pe, newb);
	//update geo

	Elems[newa].geo = Elems[oldtet[0]].geo;
	Elems[newb].geo = Elems[oldtet[0]].geo;

	if (improve_step) {
		updateQuality(newa);
		updateQuality(newb);
	}

	//delet oldtet
	if (thread_n == -1) {
		DelEle(oldtet[0]);
		DelEle(oldtet[1]);
		DelEle(oldtet[2]);
	}
	else {
		DelEle(oldtet[0], thread_n);
		DelEle(oldtet[1], thread_n);
		DelEle(oldtet[2], thread_n);
	}
	//store new tet, newa at 0, newb at 1
	oldtet[0] = newa;
	oldtet[1] = newb;
	return 1;
}

/*
* oldtet:three tet store in
* ia:point of oldtet[0]
* s.t. only flip,if it could be fliped should be deterine in last step
*     b
*  c  d  e
*	  a
*/
int DT::flip23(std::vector<int>& oldtet, int a, int thread_n) {
	std::vector<int> sph;
	int  b, c[2], d[2], e[2], i, pa, pb, pc, pd, pe, ib, newc, newd, newe;
	oldtet[1] = getNeig(oldtet[0], a);
	oldtet.resize(3);
	b = getNeigOrd(oldtet[0], a);

	//determine a,b,c,d,e order in three tet
	pa = Elems[oldtet[0]].form[a];
	pb = Elems[oldtet[1]].form[b];

	DFC(b, ib, c[1], d[1], e[1]);

	if (Elems[oldtet[1]].form[c[1]] == ghost) {//to deal with hulltet,set ghost at e
		std::swap(c[1], e[1]);
		std::swap(d[1], c[1]);
	}
	else if (Elems[oldtet[1]].form[d[1]] == ghost) {//to deal with hulltet,set ghost at e
		std::swap(d[1], e[1]);
		std::swap(d[1], c[1]);
	}
	for (i = 0; i < 4; i++) {
		if (Elems[oldtet[0]].form[i] == Elems[oldtet[1]].form[c[1]]) { c[0] = i; continue; }
		if (Elems[oldtet[0]].form[i] == Elems[oldtet[1]].form[d[1]]) { d[0] = i; continue; }
		if (Elems[oldtet[0]].form[i] == Elems[oldtet[1]].form[e[1]]) { e[0] = i; continue; }
	}
	pc = Elems[oldtet[0]].form[c[0]];
	pd = Elems[oldtet[0]].form[d[0]];
	pe = Elems[oldtet[0]].form[e[0]];

	//create new tet
	if (thread_n == -1) {
		newc = addElem(pb, pd, pa, pe);
		newd = addElem(pb, pa, pc, pe);
		newe = addElem(pb, pa, pd, pc);
	}
	else {
		newc = addElem(pb, pd, pa, pe, thread_n);
		newd = addElem(pb, pa, pc, pe, thread_n);
		newe = addElem(pb, pa, pd, pc, thread_n);
	}
	// A 2-to-3 flip is found.
	//bond adjacent
	bond(newc, 1, newd, 2);
	bond(newd, 3, newe, 2);
	bond(newc, 3, newe, 3);
	bond(newc, 0, getNeig(oldtet[0], c[0]), getNeigOrd(oldtet[0], c[0]));
	bond(newc, 2, getNeig(oldtet[1], c[1]), getNeigOrd(oldtet[1], c[1]));
	bond(newd, 0, getNeig(oldtet[0], d[0]), getNeigOrd(oldtet[0], d[0]));
	bond(newd, 1, getNeig(oldtet[1], d[1]), getNeigOrd(oldtet[1], d[1]));
	bond(newe, 0, getNeig(oldtet[0], e[0]), getNeigOrd(oldtet[0], e[0]));
	bond(newe, 1, getNeig(oldtet[1], e[1]), getNeigOrd(oldtet[1], e[1]));

	//update geo
	Elems[newc].geo = Elems[oldtet[0]].geo;
	Elems[newd].geo = Elems[oldtet[0]].geo;
	Elems[newe].geo = Elems[oldtet[0]].geo;

	//set point to tet
	setP2T(pa, newe);
	setP2T(pb, newe);
	setP2T(pc, newe);
	setP2T(pd, newe);
	setP2T(pe, newd);

	if (improve_step) {
		//in improvement step,we should update quality
		updateQuality(newc);
		updateQuality(newd);
		updateQuality(newe);
	}

	//delet oldtet
	if (thread_n == -1) {
		DelEle(oldtet[0]);
		DelEle(oldtet[1]);
	}
	else {
		DelEle(oldtet[0], thread_n);
		DelEle(oldtet[1], thread_n);
	}
	//store new tet, newa at 0, newb at 1
	oldtet[0] = newc;
	oldtet[1] = newd;
	oldtet[2] = newe;

	return 1;
}

/* flip41 to remove a pnt */
int DT::flip41(std::vector<int>& oldtet, int iNod) {
	int form[4] = { 0 };
	std::unordered_map<int, int> mp;

	for (int i = 0; i < 4; i++) {
		if (Elems[oldtet[0]].form[i] != iNod) {
			form[i] = Elems[oldtet[0]].form[i];
			mp[form[i]] = oldtet[0] + oldtet[1] + oldtet[2] + oldtet[3];
		}
		else {
			for (int j = 0; j < 4; j++) {
				if (Elems[oldtet[1]].form[j] == Elems[oldtet[0]].form[0] ||
					Elems[oldtet[1]].form[j] == Elems[oldtet[0]].form[1] ||
					Elems[oldtet[1]].form[j] == Elems[oldtet[0]].form[2] ||
					Elems[oldtet[1]].form[j] == Elems[oldtet[0]].form[3])
					continue;
				form[i] = Elems[oldtet[1]].form[j];
				mp[form[i]] = oldtet[0] + oldtet[1] + oldtet[2] + oldtet[3];
				break;
			}
		}
	}

	for (int i = 0; i < 4; i++) {
		for (int k = 0; k < 4; k++) {
			if (Elems[oldtet[i]].form[k] == iNod) {
				continue;
			}
			mp[Elems[oldtet[i]].form[k]] -= oldtet[i];
		}
	}

	int newe = -1;

	newe = addElem(form[0], form[1], form[2], form[3]);
	for (int i = 0; i < 4; i++) {
		int neigtet = getNeig(mp[form[i]], isNod_in_Tet(iNod, mp[form[i]]));
		int neigord = getNeigOrd(mp[form[i]], isNod_in_Tet(iNod, mp[form[i]]));
		bond(newe, i, neigtet, neigord);
	}

	Elems[newe].geo = Elems[oldtet[0]].geo;

	DelNod(iNod);
	for (int i = 0; i < 4; i++)
		DelEle(oldtet[i]);
	for (auto i : form) {
		setP2T(i, newe);
	}

	for (int i = 0; i < 4; i++) {
		if (Elems[newe].form[i] == ghost) {
			int p[4] = { Elems[newe].form[(i + 1) % 4] ,Elems[newe].form[(i + 2) % 4] ,Elems[newe].form[(i + 3) % 4] ,Elems[newe].form[i] };
			matchtet(p, newe);
			break;
		}
	}
	return newe;
}

/*
* fliptype 1 : 2-to-3 and 2 : 3-to-2, respectively.
* [a,b,c] is the flip face, and [d,e] is the flip edge.
* NOTE: 'pc' may be 'ghost'
*/
bool DT::flipintersectcheck(int fliptype, int a, int b, int c, int d, int e) {
	int  intTyp, rejflag = 0, i, k, intCod;
	double linep[2][3], facept[3][3], intPnt[3];
	double* pa = Nodes[a].pt;
	double* pb = Nodes[b].pt;
	double* pc = Nodes[c].pt;
	double* pd = Nodes[d].pt;
	double* pe = Nodes[e].pt;

	if (seg[0] != -1) {
		for (k = 0; k < 3; k++) {
			linep[0][k] = Nodes[seg[0]].pt[k];
			linep[1][k] = Nodes[seg[1]].pt[k];
		}
		// A constraining edge is given (e.g., for edge recovery).
		if (fliptype == 1) {
			if ((d == seg[0] || d == seg[1]) && (e == seg[0] || e == seg[1]))
				return false;
			// A 2-to-3 flip: [a,b,c] => [e,d,a], [e,d,b], [e,d,c].
			for (k = 0; k < 3; k++) {
				facept[0][k] = pd[k];
				facept[1][k] = pe[k];
			}

			for (i = 0; i < 3 && !rejflag; i++) {
				if (i == 0) {
					for (k = 0; k < 3; k++)
						facept[2][k] = pa[k];
				}
				else if (i == 1) {
					for (k = 0; k < 3; k++)
						facept[2][k] = pb[k];
				}
				else if (i == 2) {
					if (c == ghost)break;
					for (k = 0; k < 3; k++)
						facept[2][k] = pc[k];
				}
				// Test if the face [e,d,#] intersects the edge.
				dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, nullptr);
				if (intTyp == dt::GEOM_FUNC::LTI_INTERSECT_FAC ||
					intTyp == dt::GEOM_FUNC::LTI_INTERSECT_INS) {
					rejflag = 1;
				}
				else if (intTyp == dt::GEOM_FUNC::LTI_INTERSECT_EDG) {
					if (intCod == 0)//intersect with [e,d]
						rejflag = 1;
				}
			} // i
		}
		else if (fliptype == 2) {
			// A 3-to-2 flip: [e,d,a], [e,d,b], [e,d,c] => [a,b,c]
			if (c != ghost) {
				if ((a == seg[0] || a == seg[1]) + (b == seg[0] || b == seg[1]) + (c == seg[0] || c == seg[1]) == 2)
					return false;

				for (k = 0; k < 3; k++) {
					facept[0][k] = pa[k];
					facept[1][k] = pb[k];
					facept[2][k] = pc[k];
				}
				dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, intPnt);

				if (intTyp == dt::GEOM_FUNC::LTI_INTERSECT_FAC ||
					intTyp == dt::GEOM_FUNC::LTI_INTERSECT_INS) {
					rejflag = 1;
				}
			}
		}
	} // if (fc->seg[0] != NULL)

	if (fac[0] != -1 && !rejflag) {
		for (k = 0; k < 3; k++) {
			facept[0][k] = Nodes[fac[0]].pt[k];
			facept[1][k] = Nodes[fac[1]].pt[k];
			facept[2][k] = Nodes[fac[2]].pt[k];
		}
		// A constraining face is given (e.g., for face recovery).
		if (fliptype == 1) {
			// A 2-to-3 flip.
			// Test if the new edge [e,d] intersects the face.
			for (k = 0; k < 3; k++) {
				linep[0][k] = pd[k];
				linep[1][k] = pe[k];
			}
			dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, nullptr);

			if (/*intTyp == dt::GEOM_FUNC::LTI_INTERSECT_NOD ||*/
				intTyp == dt::GEOM_FUNC::LTI_INTERSECT_EDG ||
				intTyp == dt::GEOM_FUNC::LTI_INTERSECT_FAC ||
				intTyp == dt::GEOM_FUNC::LTI_INTERSECT_INS) {
				rejflag = 1;
			}
		} // if (fliptype == 1)
	} // if (fc->fac[0] != NULL)

	return rejflag;
}

/* Transport dt's Nodes and Elems to mesh*/
int DT::outMesh(Mesh& mesh, Args& args) {
    MeshStageLog stageLog(*this, "Export mesh", 1);
	int k = 0;

	//waiting for add clean delete Nodes
	clearNodesElems();
	if (periodic_P.size() != 0) {
		args.periodic_P.clear();
		for (const auto& it : periodic_P) {
			int old_key = it.first;
			int new_key = Nodes[old_key].info;

			if (new_key < 0) continue;

			for (int old_v : it.second) {
				int new_v = Nodes[old_v].info;

				if (new_v < 0) continue;

				args.periodic_P.push_back(new_key);
				args.periodic_P.push_back(new_v);
			}
		}
	}

	//update mesh's segments' info need before getMeshEdg
	mesh.S.clear();
	mesh.S.resize(SurEdgs.size());
	for (int i = 0; i < SurEdgs.size(); i++) {
		if (isDelSurEdg(i))
			continue;
		if (SurEdgs[i].constrain == 2) {
			int p1 = SurEdgs[i].iStart;
			int p2 = SurEdgs[i].iEnd;

			mesh.S[k][0] = Nodes[p1].info;
			mesh.S[k][1] = Nodes[p2].info;
			mesh.S[k][2] = SurEdgs[i].geo;
			k++;
		}
	}
	mesh.S.resize(k);

	if (args.getmeshedg) {
		//if you want to get Mesh Bnd edge based on SurTris geo information
		getMeshEdgebyGeo(mesh);
	}
	mesh.V.clear();
	mesh.V.resize(Nodes.size());
	k = 0;
	for (int i = 0; i < Nodes.size(); i++) {//have ghost
		if (isDelNod(i))
			continue;
		for (int j = 0; j < 3; j++)
			mesh.V[k][j] = Nodes[i].pt[j];
		k++;
	}
	mesh.V.resize(k);
	mesh.F.clear();
	mesh.F.resize(SurTris.size());

	k = 0;
	for (int i = 0; i < SurTris.size(); i++) {
		//SurTris[i].info > 1 ,this SurTri is divided
		if (isDelSurTri(i) || SurTris[i].info > 1)
			continue;
		mesh.F[k][3] = SurTris[i].parent;
		for (int j = 0; j < 3; j++)
			mesh.F[k][j] = SurTris[i].form[j];
		k++;
	}
	mesh.F.resize(k);
	SurTris.clear();
	SurEdgs.clear();
	mesh.T.clear();
	mesh.T.resize(Elems.size());

	k = 0;
	for (int i = 0; i < Elems.size(); i++) {
		if (ishulltet(i) || isDelEle(i))
			continue;
		mesh.T[k][4] = Elems[i].geo;
		for (int j = 0; j < 4; j++) {
			mesh.T[k][j] = Elems[i].form[j];
		}
		k++;
	}
	mesh.T.resize(k);
	return 1;
}

void DT::buildBndInfo(Mesh& mesh, Args& args, bool buildSize) {
    MeshStageLog stageLog(*this, "Build boundary topology", 2);
	int i, j, k, m;
	if (infolevel > 0)
		meshLogger->debug("Build faces&edges' TOPO.");
	//cal Area and destory it
	SurMeshClean(mesh, args);
	nSurTris = mesh.F.size(); //update num of surTris, because of surMeshClean

	SurTris.reserve(nSurTris * 1.5);
	SurEdgs.reserve(nSurTris * 2);//nSurTris*3/2


	//std::map<int, int> mp;
	//int fff = 100;

	int EdgNum = 0;
	for (i = 0; i < nSurTris; i++) {//input SurTri
		SurTris.emplace_back(SurTri());

		//Duplicate Tri
		const bool remapped = mesh.F[i][0] >= 0 && mesh.F[i][1] >= 0 && mesh.F[i][2] >= 0 &&
			(Nodes[mesh.F[i][0]].info != 0 || Nodes[mesh.F[i][1]].info != 0 || Nodes[mesh.F[i][2]].info != 0);
		if ((mesh.F[i][0] == -1 && mesh.F[i][1] == -1 && mesh.F[i][2] == -1) ||
			(remapped && BndTri.find(mesh.F[i][0], mesh.F[i][1], mesh.F[i][2]))) {
			setDelSurTri(i);
			continue;
		}

		//Degenerate triangle, because of repeated points.
		for (m = 0; m < 3; m++) {
			int idx = mesh.F[i][m];
			if (Nodes[idx].info == 0)
				SurTris[i].form[m] = idx;
			else //deal with duplicate points
				SurTris[i].form[m] = Nodes[idx].info - 1;//from 1
		}
		if (SurTris[i].form[0] == SurTris[i].form[1] ||
			SurTris[i].form[0] == SurTris[i].form[2] ||
			SurTris[i].form[1] == SurTris[i].form[2]) {
			setDelSurTri(i);
			continue;
		}
		//store hash
		auto faceEntry = BndTri.try_emplace(SurTris[i].form[0], SurTris[i].form[1], SurTris[i].form[2], i);
		if (!faceEntry.second) {
			if (!remapped) {
				setDelSurTri(i);
				continue;
			}
			*faceEntry.first = i; // Preserve replacement semantics for remapped faces.
		}
		//Set Face ID

		SurTris[i].parent = mesh.F[i][3];

		for (m = 0; m < 3; m++) {
			j = SurTris[i].form[(m + 1) % 3];
			k = SurTris[i].form[(m + 2) % 3];
			setbndpnt(j);

			AvgEdgLen += distance(Nodes[j].pt, Nodes[k].pt);
			EdgNum++;
			auto edgeEntry = BndEdg.try_emplace(j, k, static_cast<int>(SurEdgs.size()));
			if (edgeEntry.second) {//new edge
				//point_point mp edge,merge int_int to uint64_t
				int newE = SurEdgs.size();
				SurEdgs.emplace_back(SurEdg());
				SurEdgs[newE].iStart = j;
				SurEdgs[newE].iEnd = k;
				SurEdgs[newE].face.reserve(2);
				SurEdgs[newE].face.emplace_back(i);//first face
			}
			else {
				int edgid = *edgeEntry.first;
				SurEdgs[edgid].face.emplace_back(i);
			}
		}
	}//for (i = 0; i < nSurTris; i++)

	AvgEdgLen /= EdgNum;
	if (infolevel > 0)
		meshLogger->debug("AveEdgLen: {}", AvgEdgLen);

	setAllP2T();//set all point to tet

	if (buildSize) {
		buildspace(mesh, args);//build space
		//if some node don't get space from bnd
		for (int i = 0; i < Nodes.size(); i++) {
			if (isDelNod(i) || i == ghost)
				continue;
			if (Nodes[i].space != 0)
				continue;
			int neigNum = 0;
			std::unordered_set<int> Sphere_pnt;
			findSphere_pnt(i, Sphere_pnt);
			for (auto it : Sphere_pnt) {
				if (Nodes[it].space != 0) {
					Nodes[i].space += Nodes[it].space;
					neigNum++;
				}
			}
			Nodes[i].space /= neigNum;
		}
	}

	int outputnum = 0;
	if (/*args.watertightcheck|| args.constrain == 0*/ mesh.T.size()==0) {
		for (int ie = 0; ie < SurEdgs.size(); ie++) {
			if (isDelSurEdg(ie) || SurEdgs[ie].face.size() != 1)
				continue;
				if(outputnum++ < 3)
					meshLogger->debug("Single edge: {} {}", SurEdgs[ie].iStart, SurEdgs[ie].iEnd);
				AttachSeg2Pnt(ie);
		}

		for (int ie = 0; ie < SurEdgs.size(); ie++) {
			if (isDelSurEdg(ie) || SurEdgs[ie].face.size() != 1)
				continue;
			DelSingleEdge(ie);
		}
	}

	return;
}

int DT::AttachPnt2Seg(int iNod,int targetE) {
	//Try to attach in segment
	int nSe = SurEdgs.size();
	for (int ie = targetE; ie < nSe; ie++) {
		if (isDelSurEdg(ie))
			continue;
		int p1 = SurEdgs[ie].iStart;
		int p2 = SurEdgs[ie].iEnd;
		if (p1 == iNod || p2 == iNod)
			continue;
		if (P_in_Line(Nodes[iNod].pt, Nodes[p1].pt, Nodes[p2].pt))
		{
			if (infolevel > 0) meshLogger->debug("Attach {} in {},{}", iNod, p1, p2);

			int  e1, e2, manifold;
			std::vector<int> f;
			std::vector<int> e;

			for (auto itf : SurEdgs[ie].face) {
				if (std::find(SurTris[itf].form, SurTris[itf].form + 3, iNod) != SurTris[itf].form + 3) {
					auto& vec = SurEdgs[ie].face;
					setDelSurTri(itf);
					vec.erase(std::remove(vec.begin(), vec.end(), itf), vec.end());
					break;
				}
			}

			manifold = SurEdgs[ie].face.size();
			f.resize(manifold * 3);
			e.resize(manifold);
			for (int i = 0; i < manifold; i++) {
				f[i] = SurEdgs[ie].face[i];
			}

			for (int i = manifold; i < manifold * 3; i++) {
				f[i] = SurTris.size();
				SurTris.emplace_back(SurTri());//add new surTri
			}

			//add [p1,newp]
			auto e1Entry = BndEdg.try_emplace(iNod, p1, static_cast<int>(SurEdgs.size()));
			if (!e1Entry.second) {
				e1 = *e1Entry.first;
				SurEdgs[e1].face.resize(manifold + 1);
				for (int i = 0; i < manifold; i++)
					SurEdgs[e1].face[i + 1] = f[manifold + i * 2];//first face
			}
			else {
				e1 = SurEdgs.size();
				SurEdgs.emplace_back(SurEdg());
				SurEdgs[e1].iStart = p1;
				SurEdgs[e1].iEnd = iNod;
				SurEdgs[e1].info = 0;//don't recover
				SurEdgs[e1].face.resize(manifold);
				for (int i = 0; i < manifold; i++)
					SurEdgs[e1].face[i] = f[manifold + i * 2];//first face
			}

			//add [p2,newp]
			auto e2Entry = BndEdg.try_emplace(iNod, p2, static_cast<int>(SurEdgs.size()));
			if (!e2Entry.second) {
				e2 = *e2Entry.first;
				SurEdgs[e2].face.resize(manifold + 1);
				for (int i = 0; i < manifold; i++)
					SurEdgs[e2].face[i + 1] = f[manifold + i * 2 + 1];//first face
			}
			else {
				e2 = SurEdgs.size();
				SurEdgs.emplace_back(SurEdg());
				SurEdgs[e2].iStart = p2;
				SurEdgs[e2].iEnd = iNod;
				SurEdgs[e2].info = 0;//don't recover
				SurEdgs[e2].face.resize(manifold);
				for (int i = 0; i < manifold; i++)
					SurEdgs[e2].face[i] = f[manifold + i * 2 + 1];//first face
			}

			//update lostE's info,it's son is e1
			setDelSurEdg(ie);
			BndEdg.erase(p1, p2);
			//add new edge [newp,f[i].form[x]!=p1,p2]
			for (int i = 0; i < manifold; i++) {
				int p3 = -1;
				for (int j = 0; j < 3; j++) {
					if (SurTris[f[i]].form[j] != p1 && SurTris[f[i]].form[j] != p2) {
						p3 = SurTris[f[i]].form[j];
						break;
					}
				}
				e[i] = SurEdgs.size();
				BndEdg.add(p3, iNod, e[i]);
				SurEdgs.emplace_back(SurEdg());
				SurEdgs[e[i]].iStart = p3;
				SurEdgs[e[i]].iEnd = iNod;
				//SurEdgs[e[i]].parent = lostE;
				SurEdgs[e[i]].info = 0;//don't recover
				SurEdgs[e[i]].face.resize(2);
				SurEdgs[e[i]].face[0] = f[manifold + i * 2];//first face
				SurEdgs[e[i]].face[1] = f[manifold + i * 2 + 1];//first face
			}
			//updata old surTri
			for (int i = 0; i < manifold; i++) {
				setDelSurTri(f[i]);
				BndTri.erase(SurTris[f[i]].form[0], SurTris[f[i]].form[1], SurTris[f[i]].form[2]);
			}
			//updata surTri connect p1
			for (int i = manifold; i < manifold * 3; i = i + 2) {
				for (int j = 0; j < 3; j++) {
					int parentform = SurTris[f[(i - manifold) / 2]].form[j];
					SurTris[f[i]].form[j] = parentform == p2 ? iNod : parentform;
				}
				if (!BndTri.try_emplace(SurTris[f[i]].form[0], SurTris[f[i]].form[1], SurTris[f[i]].form[2], f[i]).second) {
					//have exist
					setDelSurTri(f[i]);
					continue;
				}
				SurTris[f[i]].parent = SurTris[f[(i - manifold) / 2]].parent;
				for (int j = 0; j < 3; j++) {
					int tempE = BndEdg.get(SurTris[f[i]].form[(j + 1) % 3], SurTris[f[i]].form[(j + 2) % 3]);
					//SurTris[f[i]].edgs[j] = tempE;
					if (SurTris[f[i]].form[(j + 1) % 3] != iNod && SurTris[f[i]].form[(j + 2) % 3] != iNod) {
						for (int k = 0; k < SurEdgs[tempE].face.size(); k++) {
							if (SurEdgs[tempE].face[k] == f[(i - manifold) / 2]) {
								SurEdgs[tempE].face[k] = f[i];
								break;
							}
						}
					}
				}
				SurTris[f[i]].info = 0;
			}
			//updata surTri connect p2
			for (int i = manifold + 1; i < manifold * 3; i = i + 2) {
				for (int j = 0; j < 3; j++) {
					int parentform = SurTris[f[(i - manifold - 1) / 2]].form[j];
					SurTris[f[i]].form[j] = parentform == p1 ? iNod : parentform;
				}
				if (!BndTri.try_emplace(SurTris[f[i]].form[0], SurTris[f[i]].form[1], SurTris[f[i]].form[2], f[i]).second) {
					//have exist
					setDelSurTri(f[i]);
					continue;
				}
				SurTris[f[i]].parent = SurTris[f[(i - manifold - 1) / 2]].parent;
				for (int j = 0; j < 3; j++) {
					int tempE = BndEdg.get(SurTris[f[i]].form[(j + 1) % 3], SurTris[f[i]].form[(j + 2) % 3]);
					//SurTris[f[i]].edgs[j] = tempE;
					if (SurTris[f[i]].form[(j + 1) % 3] != iNod && SurTris[f[i]].form[(j + 2) % 3] != iNod) {
						for (int k = 0; k < SurEdgs[tempE].face.size(); k++) {
							if (SurEdgs[tempE].face[k] == f[(i - manifold - 1) / 2]) {
								SurEdgs[tempE].face[k] = f[i];
								break;
							}
						}
					}
				}
				SurTris[f[i]].info = 0;
			}
			return nSe;
		}
	}

	//if (infolevel > 0)
	//	meshLogger->debug("{} Attach fail,need try to do Attach to Facet, connect to WYF~",iNod);
	return nSe;
}

void  DT::AttachSeg2Pnt(int ie) {
	//Try to attach in segment
	int p1 = SurEdgs[ie].iStart;
	int p2 = SurEdgs[ie].iEnd;

	for (int iNod = 0; iNod < Nodes.size(); iNod++) {
		if (isDelNod(iNod) || iNod ==ghost)
			continue;
		if (p1 == iNod || p2 == iNod)
			continue;
		if (P_in_Line(Nodes[iNod].pt, Nodes[p1].pt, Nodes[p2].pt))
		{
			if (infolevel > 0) meshLogger->debug("Attach {} in {},{}", iNod, p1, p2);

			if (auto* boundaryEntry = BndTri.find(iNod, p1, p2)) {
				const int boundaryIndex = *boundaryEntry;
				int delf=boundaryIndex;
				setDelSurTri(delf);
				int e1 = BndEdg.get(iNod, p1);
				int e2 = BndEdg.get(iNod, p2);
				auto& vec1 = SurEdgs[e1].face;
				vec1.erase(std::remove(vec1.begin(), vec1.end(), delf), vec1.end());
				auto& vec2 = SurEdgs[e2].face;
				vec2.erase(std::remove(vec2.begin(), vec2.end(), delf), vec2.end());
				return;
			}
			int  e1, e2, manifold;
			std::vector<int> f;
			std::vector<int> e;

			manifold = SurEdgs[ie].face.size();
			f.resize(manifold * 3);
			e.resize(manifold);
			for (int i = 0; i < manifold; i++) {
				f[i] = SurEdgs[ie].face[i];
			}
			for (int i = manifold; i < manifold * 3; i++) {
				f[i] = SurTris.size();
				SurTris.emplace_back(SurTri());//add new surTri
			}

			//add [p1,newp]
			auto e1Entry = BndEdg.try_emplace(iNod, p1, static_cast<int>(SurEdgs.size()));
			if (!e1Entry.second) {
				e1 = *e1Entry.first;
				SurEdgs[e1].face.resize(manifold + 1);
				for (int i = 0; i < manifold; i++)
					SurEdgs[e1].face[i + 1] = f[manifold + i * 2];//first face
			}
			else {
				e1 = SurEdgs.size();
				SurEdgs.emplace_back(SurEdg());
				SurEdgs[e1].iStart = p1;
				SurEdgs[e1].iEnd = iNod;
				SurEdgs[e1].info = 0;//don't recover
				SurEdgs[e1].face.resize(manifold);
				for (int i = 0; i < manifold; i++)
					SurEdgs[e1].face[i] = f[manifold + i * 2];//first face
			}

			//add [p2,newp]
			auto e2Entry = BndEdg.try_emplace(iNod, p2, static_cast<int>(SurEdgs.size()));
			if (!e2Entry.second) {
				e2 = *e2Entry.first;
				SurEdgs[e2].face.resize(manifold + 1);
				for (int i = 0; i < manifold; i++)
					SurEdgs[e2].face[i + 1] = f[manifold + i * 2 + 1];//first face
			}
			else {
				e2 = SurEdgs.size();
				SurEdgs.emplace_back(SurEdg());
				SurEdgs[e2].iStart = p2;
				SurEdgs[e2].iEnd = iNod;
				SurEdgs[e2].info = 0;//don't recover
				SurEdgs[e2].face.resize(manifold);
				for (int i = 0; i < manifold; i++)
					SurEdgs[e2].face[i] = f[manifold + i * 2 + 1];//first face
			}

			//update lostE's info,it's son is e1
			setDelSurEdg(ie);
			BndEdg.erase(p1, p2);
			//printf("%d %d %d %d\n", p1, p2, ie, SurEdgs[ie].info);
			//add new edge [newp,f[i].form[x]!=p1,p2]
			for (int i = 0; i < manifold; i++) {
				int p3 = -1;
				for (int j = 0; j < 3; j++) {
					if (SurTris[f[i]].form[j] != p1 && SurTris[f[i]].form[j] != p2) {
						p3 = SurTris[f[i]].form[j];
						break;
					}
				}
				e[i] = SurEdgs.size();
				BndEdg.add(p3, iNod, e[i]);
				SurEdgs.emplace_back(SurEdg());
				SurEdgs[e[i]].iStart = p3;
				SurEdgs[e[i]].iEnd = iNod;
				//SurEdgs[e[i]].parent = lostE;
				SurEdgs[e[i]].info = 0;//don't recover
				SurEdgs[e[i]].face.resize(2);
				SurEdgs[e[i]].face[0] = f[manifold + i * 2];//first face
				SurEdgs[e[i]].face[1] = f[manifold + i * 2 + 1];//first face
			}
			//updata old surTri
			for (int i = 0; i < manifold; i++) {
				setDelSurTri(f[i]);
				BndTri.erase(SurTris[f[i]].form[0], SurTris[f[i]].form[1], SurTris[f[i]].form[2]);
			}
			//updata surTri connect p1
			for (int i = manifold; i < manifold * 3; i = i + 2) {
				for (int j = 0; j < 3; j++) {
					int parentform = SurTris[f[(i - manifold) / 2]].form[j];
					SurTris[f[i]].form[j] = parentform == p2 ? iNod : parentform;
				}
				SurTris[f[i]].parent = SurTris[f[(i - manifold) / 2]].parent;
				for (int j = 0; j < 3; j++) {
					int tempE = BndEdg.get(SurTris[f[i]].form[(j + 1) % 3], SurTris[f[i]].form[(j + 2) % 3]);
					//SurTris[f[i]].edgs[j] = tempE;
					if (SurTris[f[i]].form[(j + 1) % 3] != iNod && SurTris[f[i]].form[(j + 2) % 3] != iNod) {
						for (int k = 0; k < SurEdgs[tempE].face.size(); k++) {
							if (SurEdgs[tempE].face[k] == f[(i - manifold) / 2]) {
								SurEdgs[tempE].face[k] = f[i];
								break;
							}
						}
					}
				}
				SurTris[f[i]].info = 0;
				BndTri.add(SurTris[f[i]].form[0], SurTris[f[i]].form[1], SurTris[f[i]].form[2], f[i]);//add this sub bnd Tris
			}
			//updata surTri connect p2
			for (int i = manifold + 1; i < manifold * 3; i = i + 2) {
				for (int j = 0; j < 3; j++) {
					int parentform = SurTris[f[(i - manifold - 1) / 2]].form[j];
					SurTris[f[i]].form[j] = parentform == p1 ? iNod : parentform;
				}
				SurTris[f[i]].parent = SurTris[f[(i - manifold - 1) / 2]].parent;
				for (int j = 0; j < 3; j++) {
					int tempE = BndEdg.get(SurTris[f[i]].form[(j + 1) % 3], SurTris[f[i]].form[(j + 2) % 3]);
					//SurTris[f[i]].edgs[j] = tempE;
					if (SurTris[f[i]].form[(j + 1) % 3] != iNod && SurTris[f[i]].form[(j + 2) % 3] != iNod) {
						for (int k = 0; k < SurEdgs[tempE].face.size(); k++) {
							if (SurEdgs[tempE].face[k] == f[(i - manifold - 1) / 2]) {
								SurEdgs[tempE].face[k] = f[i];
								break;
							}
						}
					}
				}
				SurTris[f[i]].info = 0;
				BndTri.add(SurTris[f[i]].form[0], SurTris[f[i]].form[1], SurTris[f[i]].form[2], f[i]);//add this sub bnd Tris
			}
			return;
		}
	}

	//if (infolevel > 0)
	//	meshLogger->debug("{} Attach fail,need try to do Attach to Facet, connect to WYF~",iNod);
	return;
}

void DT::DelSingleEdge(int ie) {
	const double COP_EPS = 1e-10;
	const double DIS_EPS = 1e-10;

	int p1 = SurEdgs[ie].iStart;
	int p2 = SurEdgs[ie].iEnd;

	if (SurEdgs[ie].face.empty())
		return;

	int iFace = SurEdgs[ie].face[0];
	int p3 = SurTris[iFace].form[0] + SurTris[iFace].form[1] + SurTris[iFace].form[2] - p1 - p2;

	for (int i = 0; i < 2; i++) {
		int px = (i == 0 ? p1 : p2);
		int ie2 = BndEdg.get(p3, px);
		if (ie2 < 0) continue;

		for (auto it : SurEdgs[ie2].face) {
			if (it == iFace) continue;

			// 其中不为 px, p3 的点叫 py
			int py = -1;
			for (int k = 0; k < 3; k++) {
				int pp = SurTris[it].form[k];
				if (pp != px && pp != p3) {
					py = pp;
					break;
				}
			}
			if (py < 0) continue;

			double* P1 = Nodes[p1].pt;
			double* P2 = Nodes[p2].pt;
			double* P3 = Nodes[p3].pt;
			double* PY = Nodes[py].pt;

			double ori = std::fabs(dt::GEOM_FUNC::orient3d(PY, P3, P1, P2));

			bool cross = false;
			// 近共面时，再检查线段是否相交/近相交

			if (ori < COP_EPS) {
				double dis = segmentSegmentDistance(PY, P3, P1, P2);
				if (dis < DIS_EPS) {
					cross = true;
				}
			}

			if (cross) {
				if (infolevel > 0)
					meshLogger->debug("Ignore single Edge {},{}", p1, p2);
				ignoreE(ie);
				return;
			}
		}
	}
	return;
}

//O(N^2),slow, but seldom
void DT::SurMeshClean(Mesh& mesh, Args& args){
	for(int targetF=0; targetF< mesh.F.size(); targetF++){
	
		int p1 = mesh.F[targetF][0];//one of old Tri point
		int p2 = mesh.F[targetF][1];//one of old Tri point
		int p3 = mesh.F[targetF][2];//one of old Tri point

		if (p1 < 0 || p2 < 0 || p3 < 0) continue;
		double A = calArea(Nodes[p1].pt, Nodes[p2].pt, Nodes[p3].pt);

		//std::vector<double> angle;
		//double min1 = calculateTriangleAngles(Nodes[p1].pt, Nodes[p2].pt, Nodes[p3].pt, angle);

		if (A >= 1e-15) continue;

		if (infolevel > 0)
			meshLogger->warn("Faceid:{} is degenerate unit, area: {}", targetF, A);

		//continue;
		//if (args.constrain == 1) continue;

		double L1 = distance2(Nodes[p2].pt, Nodes[p3].pt);
		double L2 = distance2(Nodes[p1].pt, Nodes[p3].pt);
		double L3 = distance2(Nodes[p1].pt, Nodes[p2].pt);

		double shortE = 1e-16;
		if (L1 < shortE || L2 < shortE || L3 < shortE) {
			//because of short edge, destroy it
			if (L2 < shortE) {
				swap(p1, p2);
			}else if (L3 < shortE) {
				swap(p1, p3);
			}

			//contract p2,p3 (destory p3)
			for (int i = 0; i < mesh.F.size(); i++) {
				int f[3] = { mesh.F[i][0],mesh.F[i][1],mesh.F[i][2] };
				for (int j = 0; j < 3; j++) {
					if (f[j] == p3) {
						if (f[(j + 1) % 3] == p2 || f[(j + 2) % 3] == p2) {
							mesh.F[i][0] = mesh.F[i][1] = mesh.F[i][2] = -1;//destory tri
						}
						else {
							mesh.F[i][j] = p2;//p3->p2
						}
						break;
					}
				}
			}
		}
		else {//three point at one edge
			// Long Edge <p2,p3>
			if (L2 > L1 && L2 > L3) swap(p1, p2);
			else if (L3 > L1 && L3 > L2) swap(p1, p3);

			mesh.F[targetF][0] = mesh.F[targetF][1] = mesh.F[targetF][2] = -1;//destory tri

			for (int i = 0; i < mesh.F.size(); i++) {
				int f[3] = { mesh.F[i][0],mesh.F[i][1],mesh.F[i][2] };
				for (int j = 0; j < 3; j++) {
					if (f[j] == p3) {
						if (f[(j + 1) % 3] == p2 || f[(j + 2) % 3] == p2) {
							int p4 = f[0] + f[1] + f[2] - p3 - p2;
							//Have p3,p2
							mesh.F[i][0] = p1;
							mesh.F[i][1] = p4;
							mesh.F[i][2] = p2;

							mesh.F.push_back({ p1,p4,p3,mesh.F[i][3] });
						}
						break;
					}
				}
			}
		}
	}
	return;
}

int DT::ColorTets() {
	int layer = 0, iElem, iSrch, ia, ib, ic, id, pb, pc, pd;

    std::vector<unsigned> marks(Elems.size(), 0);
    std::vector<int> queue;
    queue.reserve(256);
    unsigned stamp = 0;
    auto colorComponent = [&](int seed, int color) {
        return colorTetComponent(seed, color, marks, ++stamp, queue);
    };
	//find hulltet
	if (ghost != -1) {
		for (int i = 0; i < Elems.size(); i++) {
			if (isDelEle(i))
				continue;
			if (ishulltet(i)) {
				colorComponent(i, virtualID);
                ++layer; // Keep physical region numbering starting at 1.
				break;
			}
		}
	}

	std::vector<int> waitColor;

	for (int i = 0; i < Elems.size(); i++) {
		if (isDelEle(i)) //this step have no hull tet must
			continue;
		if (Elems[i].geo == -1) {//try to BFS
			for (int j = 0; j < 4; j++) {
				int neig = getNeig(i, j);
				if (neig < 0 || neig >= static_cast<int>(Elems.size()) || Elems[neig].geo == -1) {
					continue;
				}
				int Ncolor = colorComponent(i, layer++);
				//printf("%d  ", Ncolor);
				break;
			}
			if (Elems[i].geo == -1) {
				waitColor.push_back(i);
			}
		}
	}


	for (auto it:waitColor) {
		if (Elems[it].geo == -1) {//try to BFS
			int Ncolor = colorComponent(it, layer++);
		}
	}

	return layer - 1;
}

int DT::ColorTetNeig(int iElm, int color) {
    std::vector<unsigned> marks(Elems.size(), 0);
    std::vector<int> queue;
    queue.reserve(256);
    return colorTetComponent(iElm, color, marks, 1, queue);
}

int DT::colorTetComponent(int iElm, int color, std::vector<unsigned>& marks,
    unsigned stamp, std::vector<int>& q) {
    if (iElm < 0 || iElm >= static_cast<int>(Elems.size()) || isDelEle(iElm)) return 0;
    q.clear();
    marks[iElm] = stamp;
    size_t head = 0;
    q.push_back(iElm);
	int Ncolor = 0;

	while (head < q.size())
	{
		int tempE = q[head++];

		//if (Elems[tempE].geo >= 0) {
		//	printf("%d %d\n", tempE, Elems[tempE].geo);
		//	std::vector<int> ooo{tempE};			
		//	printSph_VTK(ooo, "./holetet.vtk");
		//	throw(1);
		//}

		Elems[tempE].geo = color;
		++Ncolor;

		for (int m = 0; m < 4; ++m)
		{
			int iSrch = getNeig(tempE, m);
			if (iSrch < 0 || iSrch >= static_cast<int>(Elems.size()) || isDelEle(iSrch) || marks[iSrch] == stamp)
				continue;

			int ia, ib, ic, id;
			DFC(m, ia, ib, ic, id);

			int pb = Elems[tempE].form[ib];
			int pc = Elems[tempE].form[ic];
			int pd = Elems[tempE].form[id];

			if (!isBndTri(pb, pc, pd))
			{
				marks[iSrch] = stamp;
				q.push_back(iSrch);
			}
		}
	}

	return Ncolor;
}


int DT::shellextra() {
	int ia, ib, ic, id, pb, pc, pd, nid = 0;
	std::set<int> Shell;
	for (int i = 0; i < Elems.size(); i++) {
		if (isvirtualtet(i)) {
			for (int j = 0; j < 4; j++) {
				DNC(j, ia, ib, ic, id);
				pb = Elems[i].form[ib];
				pc = Elems[i].form[ic];
				pd = Elems[i].form[id];
				if (auto* boundaryEntry = BndTri.find(pb, pc, pd)) {
					const int boundaryIndex = *boundaryEntry;
					Shell.insert(boundaryIndex);
				}
			}
		}
	}
	Mesh mesh;
	mesh.V.reserve(Nodes.size());
	mesh.F.reserve(Shell.size());
	std::vector<int> Nidx(Nodes.size());
	for (int i = 0; i < Nodes.size(); i++) {//have ghost
		if (i == ghost)
			continue;
		mesh.V.push_back({ Nodes[i].pt[0],Nodes[i].pt[1],Nodes[i].pt[2] });
		Nidx[i] = nid++;
	}
	for (auto it : Shell) {
		mesh.F.push_back({ Nidx[SurTris[it].form[0]], Nidx[SurTris[it].form[1]], Nidx[SurTris[it].form[2]] });
	}
	std::string shellname = "./Shell.vtk";
	dt::writeVTK(shellname, mesh, true);

	return 1;
}
void DT::updateSize(int iNod, Args& args) {
	Node* tempN = &Nodes[iNod];
	if (args.sizingFunc) {
		tempN->space = std::min(tempN->space, (args.sizingFunc)(tempN->pt[0], tempN->pt[1], tempN->pt[2]));
	}
	return;
}
//build init node spcace
void DT::buildspace(Mesh& mesh, Args& args) {
	int i, p1, p2, p3;
	double d1, d2, d3;

	//build space from boundary
	for (i = 0; i < SurTris.size(); i++) {
		if (SurTris[i].info < 0 || SurTris[i].info >1) {//don't recovered or split surTri
			continue;
		}
		p1 = SurTris[i].form[0];
		p2 = SurTris[i].form[1];
		p3 = SurTris[i].form[2];

		d1 = distance(Nodes[p1].pt, Nodes[p2].pt);
		d2 = distance(Nodes[p2].pt, Nodes[p3].pt);
		d3 = distance(Nodes[p3].pt, Nodes[p1].pt);

		Nodes[p1].space += d1 + d3;
		Nodes[p2].space += d2 + d1;
		Nodes[p3].space += d3 + d2;

		Nodes[p1].info += 2;
		Nodes[p2].info += 2;
		Nodes[p3].info += 2;
	}
	for (i = 0; i < Nodes.size(); i++) {//from 0
		if (isDelNod(i))
			continue;
		if ((Nodes[i].info & 0xEFFFFFF) == 0)//inner point
		{
			//if inner point,don't give size = 0, lead to refine error
			Nodes[i].space = AvgEdgLen;
			continue;
		}
		Nodes[i].space /= (Nodes[i].info & 0xEFFFFFF);
		Nodes[i].info &= ~0xEFFFFFF;
	}


	//if input bndpnt size, don't do global smooth now, we don't have anisotropic opt now
	if (mesh.pointSize.size() > 0) {
		for (i = 0; i < Nodes.size(); i++) {//from 0
			if (isDelNod(i))
				continue;
			if (i >= mesh.pointSize.size())
				continue;
			double s = mesh.pointSize[i];
			if (s <= 0) {
				meshLogger->warn("Input Size <= 0");
				continue;
			}
			Nodes[i].space = std::min(Nodes[i].space, s);
		}
	}

	if (args.sizingFunc) {
		for (i = 0; i < Nodes.size(); i++) {//from 0, because of ghost
			if (isDelNod(i))
				continue;
			//querySizeNum++;
			double s = (args.sizingFunc)(Nodes[i].pt[0], Nodes[i].pt[1], Nodes[i].pt[2]);
			//querySizeTime += getTime(t_1, t_2);
			if (s <= 0) {
				meshLogger->warn("Input Size <= 0");
				continue;
			}
			Nodes[i].space = std::min(Nodes[i].space, s);
		}
	}

	return;
}

/*
* Version gravity
* input i : Elems id
* output -1:don't need add point
* output else :new point id
*/
// Compatibility helper: direct callers retain the original element-bit contract.
// MeshRefine uses an explicit call-local rejection map instead.
int DT::createRefineCandidate(int i, Args& args, bool& rejected) {
	int  j, k;
	double pnt[3] = { 0 }, dp = 0, dt[4] = { 0 }, SizeAlpha = args.size;

	if (isDelEle(i))
		return -1;
	if (rejected)
		return -1;

	//get barycenter of this tet
	calBarycenter(i, pnt);

	for (j = 0; j < 4; j++) {
		dt[j] = distance(pnt, Nodes[Elems[i].form[j]].pt);
		dp += Nodes[Elems[i].form[j]].space + dt[j] * (growsize - 1);
	}
	dp /= 4.0;

	//dp=dbl_+max;
	//for (j = 0; j < 4; j++) {
	//	dt[j] = distance(pnt, Nodes[Elems[i].form[j]].pt);
	//	dp = std::min(dp, Nodes[Elems[i].form[j]].space + dt[j] * (growsize - 1));
	//}

	if (maxEdge != -1) {
		dp = std::min(maxEdge, dp);
	}

	if (minEdge != -1) {
		dp = std::max(minEdge, dp);
	}

	//Determine whether to insert this point
	for (j = 0; j < 4; j++) {
		if (dt[j] < SizeAlpha * dp && dt[j] < SizeAlpha * Nodes[Elems[i].form[j]].space)
			break;
	}

	if (j == 4) {//this point allow to insert
		int newN = addNode(pnt[0], pnt[1], pnt[2], dp);
		return newN;
	}
	else {
		//set this tet has don't need check
		rejected = true;
	}
	return -1;//create new point fail
}
