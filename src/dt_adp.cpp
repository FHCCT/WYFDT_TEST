#include "./dt.h"

int DT::adaptation_by_pError(Mesh& mesh, Args& args, std::vector<std::array<double, 4>>& addVertex, double GrowRatio) {
	int p1, p2;

	/************************ Prepare ******************************/
	//Rebuild TOPO
	if (!dt_init(mesh, args)) return 1;
    MeshStageLog stageLog(*this, "adaptation_by_pError");
	buildTetInfo(mesh, args);

	//Get min volume
	updateminVolume();
	//determine init tet num
	int initHullTet = Elems.size() - mesh.T.size();	//dt's hull tet
	int initBoxtet = 0;								//HFSS' box tet
	for (int i = 0; i < Elems.size(); i++) {
		if (Elems[i].geo == 0) {
			initBoxtet++;
		}
	}
	int inTetNum = mesh.T.size() - initBoxtet;		//wait split tet

	//sort Vertex base on Error
	std::sort(addVertex.begin(), addVertex.end(),
		[](const std::array<double, 4>& a, const std::array<double, 4>& b) {
			return a[3] > b[3];
		});

	modifyBnd = true;	//set Bnd Modification
	Type_Vertex_Edg(args.adpangle, mesh); //determine the constrained edge by geoinfo

	/************************ Where to split ******************************/
	EdgeHasher<int> EdgCount;
	std::vector<std::pair<int, int>> EdgOrd;
	for (int i = 0; i < addVertex.size(); i++) {
		int te = 0;
		while (1) {
			if (!isDelEle(te) && !ishulltet(te))
				break;
			te++;
		}

		int tempV = addNode(addVertex[i][0], addVertex[i][1], addVertex[i][2], 0);
		int loc = locate_pnt(tempV, te);
		DelNod(tempV);

		double maxEdgLen = 0;
		int maxEdgIdx = -1;
		for (int j = 0; j < 6; j++) {
			int p1 = Elems[te].form[Egid[j][0]];
			int p2 = Elems[te].form[Egid[j][1]];
			////Old Edg idndex may be destroy
			double dis = distance(Nodes[p1].pt, Nodes[p2].pt);

			if (auto* boundaryEntry = BndEdg.find(p1, p2)) {
				const int boundaryIndex = *boundaryEntry;
				//dis *= 1.1;
				int Edgid = boundaryIndex;
				if(SurEdgs[Edgid].constrain > 0)
					dis *= 1.25;
			}
			if (dis > maxEdgLen) {
				maxEdgLen = dis;
				maxEdgIdx = j;
			}
		}
		p1 = Elems[te].form[Egid[maxEdgIdx][0]];
		p2 = Elems[te].form[Egid[maxEdgIdx][1]];

		EdgOrd.push_back({ p1,p2 });
		EdgCount.try_emplace(p1, p2, 1);
	}

	/************************ Start Refine ******************************/
	for (auto it : EdgOrd) {
		if (Elems.size() - initHullTet - initBoxtet >= inTetNum * GrowRatio) {
			break;
		}

		p1 = it.first;
		p2 = it.second;

		if (!EdgCount.erase(p1, p2))
			continue;

		if (auto* boundaryEntry = BndEdg.find(p1, p2)) {
			const int boundaryIndex = *boundaryEntry;
			int Edgid = boundaryIndex;
			if (isDelSurEdg(Edgid))
				continue;
			int ret = splitEdg(Edgid);
		}
		else {
			double addV[3] = { 0.5 * (Nodes[p1].pt[0] + Nodes[p2].pt[0]),
			 0.5 * (Nodes[p1].pt[1] + Nodes[p2].pt[1]),
			 0.5 * (Nodes[p1].pt[2] + Nodes[p2].pt[2]) 
			};

			int iNod = addNode(addV[0], addV[1], addV[2], 0);

			int searchtet = getP2T(p1);
			locate_pnt(iNod, searchtet);
			if (isvirtualtet(searchtet)) {
				DelNod(iNod);
				continue;
			}

			//add inner point
			std::vector<int> nearTets = { getP2T(p1) };
			int BW_recall = BW_insert_vertex(iNod, nearTets, 1);
			if (BW_recall == 1) {//success
				continue;
			}
			else {//failed
				DelNod(iNod);
			}
		}
	}

	int ne = 0;
	for (int i = 0; i < Elems.size(); i++) {
		if (isDelEle(i) || ishulltet(i))
			continue;
		ne++;
	}
	meshLogger->debug("Init tet: {}  Refine tet: {}  Ratio: {:.3f}", inTetNum, ne, 1.0 * ne / inTetNum);
	/************************ Mesh Improvement ******************************/
	args.optlevel = 7;
	MeshImprove(args);

	RemoveTet(args);
	outMesh(mesh, args);

	meshLogger->debug("Init tet: {}  Adaptation tet: {}  Ratio: {:.3f}", inTetNum, mesh.T.size() - initBoxtet, 1.0 * (mesh.T.size() - initBoxtet) / inTetNum);
	return 1;
}

int DT::adaptation_by_Tetid (Mesh& mesh, Args& args, std::vector<int> refine_tri_id, std::vector<int> refine_tet_id, double GrowRatio) {
	int p1, p2;

	/************************ Prepare ******************************/
	//Rebuild TOPO
	if (!dt_init(mesh, args)) return 1;
    MeshStageLog stageLog(*this, "adaptation_by_Tetid ");
	buildTetInfo(mesh, args);

	int inTetNum = mesh.T.size();
	int initHullTet = Elems.size() - mesh.T.size();	//dt's hull tet

	//Split Long Edge
	modifyBnd = true;	//set Bnd Modification
	//determine the constrained edge by geoinfo
	Type_Vertex_Edg(args.adpangle, mesh);

	//Get min volume
	updateminVolume();
	/************************ Where to split ******************************/
	EdgeHasher<int> EdgCount;
	std::vector<std::array<int, 3>> EdgOrd;

	// Refine by tri
	for (int i = 0; i < refine_tri_id.size(); i++) {
		int te = refine_tri_id[i];
		double maxEdgLen = 0;
		int maxEdgIdx = -1;
		for (int j = 0; j < 3; j++) {
			p1 = SurTris[te].form[j];
			p2 = SurTris[te].form[(j + 1) % 3];
			////Old Edg idndex may be destroy
			double dis = distance(Nodes[p1].pt, Nodes[p2].pt);

			int Edgid = BndEdg.get(p1, p2);
			if (SurEdgs[Edgid].constrain > 0)
				dis *= 1.25;

			if (dis > maxEdgLen) {
				maxEdgLen = dis;
				maxEdgIdx = j;
			}
		}
		p1 = SurTris[te].form[maxEdgIdx];
		p2 = SurTris[te].form[(maxEdgIdx + 1) % 3];

		EdgOrd.push_back({ p1,p2,i });
		EdgCount.try_emplace(p1, p2, 1);
	}

	// Refine by tet
	for (int i = 0; i < refine_tet_id.size(); i++) {
		int te = refine_tet_id[i];
		double maxEdgLen = 0;
		int maxEdgIdx = -1;
		for (int j = 0; j < 6; j++) {
			p1 = Elems[te].form[Egid[j][0]];
			p2 = Elems[te].form[Egid[j][1]];
			////Old Edg idndex may be destroy
			double dis = distance(Nodes[p1].pt, Nodes[p2].pt);

			if (auto* boundaryEntry = BndEdg.find(p1, p2)) {
				const int boundaryIndex = *boundaryEntry;
				//dis *= 1.1;
				int Edgid = boundaryIndex;
				if (SurEdgs[Edgid].constrain > 0)
					dis *= 1.12;//1.118 1.25
			}

			if (dis > maxEdgLen) {
				maxEdgLen = dis;
				maxEdgIdx = j;
			}
		}
		p1 = Elems[te].form[Egid[maxEdgIdx][0]];
		p2 = Elems[te].form[Egid[maxEdgIdx][1]];

		EdgOrd.push_back({ p1,p2,i });
		EdgCount.try_emplace(p1, p2, 1);
	}

	/************************ Start Refine ******************************/
	for (auto it : EdgOrd) {
		p1 = it[0];
		p2 = it[1];

		if (Elems.size() - initHullTet >= inTetNum * GrowRatio) {
//#ifdef _WIN32
//			Mesh tempmesh;
//			tempmesh.V = mesh.V;
//
//			for (int i = 0; i < it[2]; i++) {
//				tempmesh.T.push_back(mesh.T[refine_tet_id[i]]);
//			}
//
//			// Current time: hour_minute_second
//			std::time_t now = std::time(nullptr);
//			std::tm localTime;
//			localtime_s(&localTime, &now);
//
//			std::ostringstream oss;
//			oss << "./UsedRefine_"
//				<< std::put_time(&localTime, "%H_%M_%S")
//				<< ".vtk";
//
//			std::string tempfilename = oss.str();
//
//			writeVTK(tempfilename, tempmesh, false);
//#endif
			break;
		}

		if (!EdgCount.erase(p1, p2))
			continue;

		if (auto* boundaryEntry = BndEdg.find(p1, p2)) {
			const int boundaryIndex = *boundaryEntry;
			int Edgid = boundaryIndex;
			if (isDelSurEdg(Edgid))
				continue;
			int ret = splitEdg(Edgid);
		}
		else {
			double addV[3] = { 0.5 * (Nodes[p1].pt[0] + Nodes[p2].pt[0]),
			 0.5 * (Nodes[p1].pt[1] + Nodes[p2].pt[1]),
			 0.5 * (Nodes[p1].pt[2] + Nodes[p2].pt[2]) };

			int iNod = addNode(addV[0], addV[1], addV[2], 0);

			int searchtet = getP2T(p1);
			locate_pnt(iNod, searchtet);
			if (isvirtualtet(searchtet)) {
				DelNod(iNod);
				continue;
			}

			//add inner point
			std::vector<int> nearTets = { searchtet };
			int BW_recall = BW_insert_vertex(iNod, nearTets, 1);
			if (BW_recall == 1) {//success
				continue;
			}
			else {//failed
				DelNod(iNod);
			}
		}
	}

	int ne = 0;
	for (int i = 0; i < Elems.size(); i++) {
		if (isDelEle(i) || ishulltet(i))
			continue;
		ne++;
	}
	meshLogger->info("Init tet: {}  Refine tet: {}  Ratio: {:.3f}", inTetNum, ne, 1.0 * ne / inTetNum);
	/************************ Mesh Improvement ******************************/
	args.optlevel = 7;
	MeshImprove(args);

	RemoveTet(args);
	outMesh(mesh, args);


	meshLogger->info("Init tet: {}  Adaptation tet: {}  Ratio: {:.3f}", inTetNum, mesh.T.size(), 1.0 * mesh.T.size() / inTetNum);
	return 1;
}

int DT::adaptation_by_YunBoSzieControl(
	Mesh& mesh, Args& args, double lamasize,
	std::unordered_map<int, double>& facetSize,
	std::unordered_map<int, int>& facetNum,
	std::unordered_map<int, double>& elementSize,
	std::unordered_map<int, int>& elementNum) {

	//Rebuild TOPO
	if (!dt_init(mesh, args)) return 1;
    MeshStageLog stageLog(*this, "adaptation_by_YunBoSzieControl", 1, MeshStageSummary::MeshCount);
	buildTetInfo(mesh, args);

	//Split Long Edge
	modifyBnd = true;	//set Bnd Modification
	fliplevel = 0;		//set a big flip level, if modify bnd

	//determine the constrained edge by geoinfo
	Type_Vertex_Edg(args.adpangle, mesh);
	
	//Mesh Improvement
	args.optlevel = 20;
	args.optloop = std::max(args.optloop, 4);
	improve_init(args);
	OptSizeControl_yunbo(lamasize, facetSize, facetNum, elementSize, elementNum, args);

	RemoveTet(args);
	outMesh(mesh, args);
    stageLog.finish(mesh.T.size());

	return 1;
}

//try it's best coarse
int DT::adaptation_Coarse(Mesh& mesh, Args& args) {
	int inTetNum = mesh.T.size();

	//Rebuild TOPO
	if (!dt_init(mesh, args)) return 1;
    MeshStageLog stageLog(*this, "adaptation_Coarse");
	buildTetInfo(mesh, args);

	//Split Long Edge
	modifyBnd = true;	//set Bnd Modification
	fliplevel = 0;		//set a big flip level, if modify bnd

	//determine the constrained edge by geoinfo
	Type_Vertex_Edg(args.adpangle, mesh);

	//Mesh Improvement
	args.optlevel = 8;

    // Coarsening relaxes this bound only for its own call.
    struct RestoreUpper {
        double& value;
        double saved;
        ~RestoreUpper() { value = saved; }
    } restoreUpper{ani_upper, ani_upper};
	ani_upper = 1e6;

	MeshImprove(args);
	RemoveTet(args);
	outMesh(mesh, args);

	meshLogger->debug("Init tet: {}  Adaptation tet: {}  Reduce to: {:.3f}%", inTetNum, mesh.T.size(), 100.0 * mesh.T.size() / inTetNum);
	return 1;
}

int DT::adaptation_by_pSize(Mesh& mesh, Args& args) {
	if (!dt_init(mesh, args)) return 1;
    MeshStageLog stageLog(*this, "adaptation_by_pSize");
	//Rebuild TOPO
	buildTetInfo(mesh, args);
	modifyBnd = true;	//set Bnd Modification
	fliplevel = 0;		//set a big flip level, if modify bnd
	Type_Vertex_Edg(args.adpangle, mesh); //determine the constrained edge by geoinfo

	args.optlevel = 6;
	MeshImprove(args);

	RemoveTet(args);
	outMesh(mesh, args);

	return 1;
}

int DT::adaptation_by_ani(Mesh& mesh, Args& args, std::vector<std::array<double, 6>>& anisol, std::vector<int> lockFactes, std::vector<int> lockVertex) {

	if (!dt_init(mesh, args)) return 1;
    MeshStageLog stageLog(*this, "adaptation_by_ani");
	AniSol = anisol;
	//Smooth_size_ani(mesh, anisol, lockFactes, lockVertex);
	//Rebuild TOPO
	buildTetInfo(mesh, args);
	modifyBnd = true;	//set Bnd Modification
	fliplevel = 0;		//set a big flip level, if modify bnd
	Type_Vertex_Edg(args.adpangle, mesh); //determine the constrained edge by geoinfo
	lockingPass(lockFactes, lockVertex);

	args.optlevel = 9;
	MeshImprove(args);

	RemoveTet(args);
	outMesh(mesh, args);

	meshLogger->debug("Final Tet: {}", mesh.T.size());
	return 1;
}

int DT::lockingPass(std::vector<int> lockFactes, std::vector<int> lockVertex) {
	for (int i = 0; i < lockFactes.size(); i++) {
		lockF.insert(lockFactes[i]);
		for (int j = 0; j < 3; j++) {
			int p1 = SurTris[lockFactes[i]].form[j];
			int p2 = SurTris[lockFactes[i]].form[(j + 1) % 3];
			lockV.insert(p1);
			lockE.insert(BndEdg.get(p1, p2));
		}
	}
	for (int i = 0; i < lockVertex.size(); i++) {
		lockV.insert(lockVertex[i]);
	}
	return 0;
}

int DT::optimization(Mesh& mesh, Args& args) {
	if (!dt_init(mesh, args)) return 1;
    MeshStageLog stageLog(*this, "optimization");
	//Rebuild TOPO
	buildTetInfo(mesh, args);
	MeshImprove(args);

	RemoveTet(args);
	outMesh(mesh, args);

	return 1;
}

// for adaptive mesh
void DT::buildTetInfo(Mesh& mesh, Args& args) {
    MeshStageLog stageLog(*this, "Rebuild topology", 1);
	// add Nodes
	buildPntInfo(mesh);

	//add Elem
	for (int i = 0; i < mesh.T.size(); i++) {
		int newE;
		double ori = dt::GEOM_FUNC::orient3d(
			Nodes[mesh.T[i][0]].pt,
			Nodes[mesh.T[i][1]].pt,
			Nodes[mesh.T[i][3]].pt,
			Nodes[mesh.T[i][2]].pt);
		if (ori >= 0) {
			newE = addElem(mesh.T[i][0], mesh.T[i][1], mesh.T[i][2], mesh.T[i][3]);
		}
		else {
			meshLogger->error("Tet: {} is invert.", i);
			newE = addElem(mesh.T[i][0], mesh.T[i][1], mesh.T[i][3], mesh.T[i][2]);
		}
		Elems[newE].geo = mesh.T[i][4];
	}

	//build neig info
	TriHasher<int64_t> Tri;
	int ia, ib, ic, id, p1, p2, p3;
	for (int i = 0; i < Elems.size(); i++) {
		for (int j = 0; j < 4; j++) {
			DNC(j, ia, ib, ic, id);
			p1 = Elems[i].form[ib];
			p2 = Elems[i].form[ic];
			p3 = Elems[i].form[id];
			auto faceEntry = Tri.try_emplace(p1, p2, p3, (((int64_t)i << 2) | j));
			if (!faceEntry.second) {
				int64_t neiginfo = *faceEntry.first;
				bond(i, j, neiginfo >> 2, neiginfo & 3);
			}
		}
	}

	//build hull tet, and don't build space
	args.constrain = 1;
	buildBndInfo(mesh, args, true);
	ghost = addNode();
	int firstHulltet = Elems.size();
	for (int i = 0; i < firstHulltet; i++) {
		for (int j = 0; j < 4; j++) {
			if (getNeig(i, j) == -1) {
				DFC(j, ia, ib, ic, id);
				p1 = Elems[i].form[ib];
				p2 = Elems[i].form[ic];
				p3 = Elems[i].form[id];
				int newE = addElem(p1, p3, p2, ghost);
				bond(i, j, newE, 3);
				setvirtualtet(newE);
			}
		}
	}

    if (firstHulltet < static_cast<int>(Elems.size())) setP2T(ghost, firstHulltet);

	EdgeHasher<int64_t> EdgeAdj;
	EdgeHasher<bool> EdgeNum;
	//build hull tet topo
	for (int i = firstHulltet; i < Elems.size(); i++) {
		for (int j = 0; j < 3; j++) {
			DNC(j, ia, ib, ic, id);
			p1 = Elems[i].form[ib == 3 ? id : ib];
			p2 = Elems[i].form[ic == 3 ? id : ic];
			auto edgeEntry = EdgeAdj.try_emplace(p1, p2, (((int64_t)i << 2) | j));
			if (!edgeEntry.second) {
				if (EdgeNum.try_emplace(p1, p2, true).second) {
					int64_t neiginfo = *edgeEntry.first;
					bond(i, j, neiginfo >> 2, neiginfo & 3);
				}
				else {
					//try lock it
					int eid = BndEdg.get(p1, p2);
					lockE.insert(eid);
					//if (meshLogger->level() != spdlog::level::off) printf("Non-main-fold: %d,%d\n", p1, p2);
					meshLogger->warn("Non-main-fold: {},{}", p1, p2);
					//throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
				}

			}
		}
	}
    // Local star traversals require one connected incident fan. A vertex shared
    // by disconnected bodies must stay fixed; otherwise deleting it updates only
    // one fan and leaves other tetrahedra referencing a deleted vertex.
    std::vector<int> incidentCount(Nodes.size(), 0);
    for (const auto& tet : Elems)
        for (int point : tet.form) ++incidentCount[point];
    std::unordered_set<int> disconnectedVertices;
    std::vector<int> sphere;
    for (int point = 0; point < static_cast<int>(Nodes.size()); ++point) {
        if (point == ghost || incidentCount[point] == 0) continue;
        if (findSphere(point, sphere) != incidentCount[point]) {
            lockV.insert(point);
            disconnectedVertices.insert(point);
        }
    }
    for (int edge = 0; edge < static_cast<int>(SurEdgs.size()); ++edge)
        if (disconnectedVertices.count(SurEdgs[edge].iStart) ||
            disconnectedVertices.count(SurEdgs[edge].iEnd)) lockE.insert(edge);
	return;
}