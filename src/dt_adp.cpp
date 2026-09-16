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
        if (te < 0 || te >= static_cast<int>(SurTris.size()) || isDelSurTri(te)) continue;
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
    lockF.insert(lockFactes.begin(), lockFactes.end());
    lockV.insert(lockVertex.begin(), lockVertex.end());
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
        if (lockFactes[i] < 0 || lockFactes[i] >= static_cast<int>(SurTris.size()) || isDelSurTri(lockFactes[i])) continue;
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

void DT::mergeThinLayers(Mesh& mesh, TriHasher<int64_t>& tetFaces,
    const std::vector<double>& tetVolumes) {
    const int nt = static_cast<int>(mesh.T.size());
    if (nt == 0 || mesh.F.empty()) return;
    // Dimensionless tolerances: small component, small thickness, small angle.
    constexpr double maxVolumeFraction = 1e-6;
    constexpr double maxVolumeAreaRatio = 1e-3;
    constexpr double maxDihedral = PI / 180.0;
    const double parallelCos = std::cos(5.0 * PI / 180.0);
    std::vector<std::array<int, 4>> boundary(nt, std::array<int, 4>{{-1,-1,-1,-1}});
    for (int f = 0; f < static_cast<int>(mesh.F.size()); ++f) {
        const auto& face = mesh.F[f];
        if (face[0] < 0 || face[1] < 0 || face[2] < 0) continue;
        const int64_t* found = tetFaces.find(face[0], face[1], face[2]);
        if (!found) return; // Incomplete input constraints: leave cleanup to the caller.
        const int t = static_cast<int>(*found >> 2), j = static_cast<int>(*found & 3);
        if (boundary[t][j] >= 0) return; // Ambiguous duplicate constraint.
        boundary[t][j] = f;
        const int other = getNeig(t, j);
        if (other >= 0) boundary[other][getNeigOrd(t, j)] = f;
    }
    const bool hasProtectedPoints = !periodic_P.empty() || !lockV.empty() || !mesh.S.empty();
    std::vector<unsigned char> protectedPoint(hasProtectedPoints ? Nodes.size() : 0, 0);
    for (const auto& point : periodic_P) protectedPoint[point.first] = 1;
    for (int point : lockV) if (point >= 0 && point < static_cast<int>(Nodes.size())) protectedPoint[point] = 1;
    for (const auto& edge : mesh.S)
        for (int j = 0; j < 2; ++j)
            if (edge[j] >= 0 && edge[j] < static_cast<int>(Nodes.size())) protectedPoint[edge[j]] = 1;
    struct Component {
        size_t begin, end;
        double volume;
        int geo;
        bool closed;
        Component(size_t first, int region) : begin(first), end(first), volume(0), geo(region), closed(true) {}
    };
    std::vector<Component> components;
    std::vector<int> label(nt, -1), order;
    order.reserve(nt);
    double totalVolume = 0;
    // Constraints split connected components even when both sides share a geo ID.
    for (int seed = 0; seed < nt; ++seed) {
        if (label[seed] >= 0) continue;
        const int c = static_cast<int>(components.size());
        components.emplace_back(order.size(), Elems[seed].geo);
        Component& component = components.back();
        label[seed] = c;
        order.push_back(seed);
        for (size_t pos = component.begin; pos < order.size(); ++pos) {
            const int t = order[pos];
            const auto& tet = Elems[t];
            const double volume = tetVolumes[t];
            if (!(volume > 0) || !std::isfinite(volume)) component.closed = false;
            else component.volume += volume;
            if (component.closed && hasProtectedPoints)
                for (int j = 0; j < 4; ++j)
                    if (protectedPoint[tet.form[j]]) component.closed = false;
            for (int j = 0; j < 4; ++j) {
                const int other = getNeig(t, j), f = boundary[t][j];
                if (other < 0) { component.closed = false; continue; }
                if (f >= 0) {
                    if (lockF.count(f)) component.closed = false;
                    continue;
                }
                // An already visited cell in this component has the same region.
                // Avoid fetching its full element again just to confirm that fact.
                const int otherLabel = label[other];
                if (otherLabel == c) continue;
                if (Elems[other].geo != component.geo) { component.closed = false; continue; }
                if (otherLabel < 0) { label[other] = c; order.push_back(other); }
            }
        }
        component.end = order.size();
        totalVolume += component.volume;
    }
    if (!(totalVolume > 0) || !std::isfinite(totalVolume)) return;
    const auto faceNormal = [&](int t, int j) -> std::array<double, 3> {
        const auto& face = mesh.F[boundary[t][j]];
        const double* a = Nodes[face[0]].pt;
        const double* b = Nodes[face[1]].pt;
        const double* d = Nodes[face[2]].pt;
        std::array<double,3> n = {{(b[1]-a[1])*(d[2]-a[2])-(b[2]-a[2])*(d[1]-a[1]),
            (b[2]-a[2])*(d[0]-a[0])-(b[0]-a[0])*(d[2]-a[2]),
            (b[0]-a[0])*(d[1]-a[1])-(b[1]-a[1])*(d[0]-a[0])}};
        const double* opposite = Nodes[Elems[t].form[j]].pt;
        if (n[0]*(opposite[0]-a[0])+n[1]*(opposite[1]-a[1])+n[2]*(opposite[2]-a[2]) > 0)
            for (double& value : n) value = -value;
        return n;
    };
    const auto norm = [](const std::array<double,3>& n) { return std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]); };
    std::vector<int> target(components.size());
    for (int c = 0; c < static_cast<int>(components.size()); ++c) target[c] = c;
    std::vector<unsigned char> recipient(components.size(), 0);
    std::set<int> mergedFaceIDs;
    std::vector<double> angles;
    for (int c = 0; c < static_cast<int>(components.size()); ++c) {
        const Component& component = components[c];
        if (recipient[c] || !component.closed || component.geo < 0 || component.volume > totalVolume * maxVolumeFraction) continue;
        std::map<int,double> areaByID;
        double totalArea = 0, minAngle = PI;
        bool valid = true;
        for (size_t pos = component.begin; pos < component.end; ++pos) {
            const int t = order[pos];
            const auto& tet = Elems[t];
            double lo, hi;
            if (!CalDihedral(Nodes[tet.form[0]].pt, Nodes[tet.form[1]].pt,
                Nodes[tet.form[3]].pt, Nodes[tet.form[2]].pt, lo, hi, angles)) { valid = false; break; }
            minAngle = std::min(minAngle, lo);
            for (int j = 0; j < 4; ++j) if (boundary[t][j] >= 0) {
                if (label[getNeig(t,j)] == c) { valid = false; break; }
                const double area = 0.5 * norm(faceNormal(t,j));
                if (!(area > 0) || !std::isfinite(area)) { valid = false; break; }
                totalArea += area;
                areaByID[mesh.F[boundary[t][j]][3]] += area;
            }
            if (!valid) break;
        }
        if (!valid || minAngle > maxDihedral || !(totalArea > 0) ||
            component.volume > maxVolumeAreaRatio * totalArea * std::sqrt(totalArea)) continue;
        int faceID = -1;
        double pairedArea = 0;
        for (const auto& entry : areaByID)
            if (entry.second > pairedArea) { faceID = entry.first; pairedArea = entry.second; }
        // Permit small side caps, but require the two main sheets to share an ID.
        if (pairedArea < 0.8 * totalArea) continue;
        std::array<double,3> reference = {{0,0,0}};
        double largest = 0;
        for (size_t pos = component.begin; pos < component.end; ++pos)
            for (int j = 0; j < 4; ++j) {
                const int t = order[pos], f = boundary[t][j];
                if (f < 0 || mesh.F[f][3] != faceID) continue;
                const auto n = faceNormal(t,j);
                const double length = norm(n);
                if (length > largest) { largest = length; reference = n; }
            }
        if (!(largest > 0)) continue;
        for (double& value : reference) value /= largest;
        double sideArea[2] = {0,0};
        int sideGeo[2] = {-1,-1}, neighbor[2] = {-1,-1};
        for (size_t pos = component.begin; pos < component.end; ++pos)
            for (int j = 0; j < 4; ++j) {
                const int t = order[pos], f = boundary[t][j];
                if (f < 0 || mesh.F[f][3] != faceID) continue;
                const auto n = faceNormal(t,j);
                const double length = norm(n);
                const double cosine = (n[0]*reference[0]+n[1]*reference[1]+n[2]*reference[2])/length;
                if (std::abs(cosine) < parallelCos) continue;
                const int side = cosine > 0 ? 0 : 1;
                const int other = label[getNeig(t,j)], geo = components[other].geo;
                if (geo < 0 || (sideGeo[side] >= 0 && sideGeo[side] != geo)) valid = false;
                sideGeo[side] = geo;
                sideArea[side] += length * 0.5;
                if (neighbor[side] < 0 || components[other].volume > components[neighbor[side]].volume)
                    neighbor[side] = other;
            }
        if (!valid || sideArea[0] < pairedArea*0.2 || sideArea[1] < pairedArea*0.2 ||
            sideArea[0]+sideArea[1] < pairedArea*0.95) continue;
        int chosen = -1;
        for (int side = 0; side < 2; ++side) {
            const int other = neighbor[side];
            // Prefer a matching region even when its adjacent component is smaller.
            // Keep each recipient intact for this pass; no cycles or cascading relabels.
            if (other < 0 || target[other] != other) continue;
            if (chosen < 0 || (sideGeo[side] == component.geo && sideGeo[chosen] != component.geo) ||
                ((sideGeo[side] == component.geo) == (sideGeo[chosen] == component.geo) && sideArea[side] > sideArea[chosen])) chosen = side;
        }
        if (chosen < 0) continue;
        target[c] = neighbor[chosen];
        recipient[target[c]] = 1;
        mergedFaceIDs.insert(faceID);
    }
    if (mergedFaceIDs.empty()) return;
    // Only merged components change. Their boundary lists also contain every
    // constraint that can become internal; untouched components need no writes.
    for (int c = 0; c < static_cast<int>(components.size()); ++c) {
        if (target[c] == c) continue;
        const int geo = components[target[c]].geo;
        for (size_t pos = components[c].begin; pos < components[c].end; ++pos) {
            const int t = order[pos];
            Elems[t].geo = mesh.T[t][4] = geo;
            for (int j = 0; j < 4; ++j) {
                const int f = boundary[t][j];
                if (f < 0) continue;
                const int other = getNeig(t,j);
                if (other < 0 || components[target[label[other]]].geo != geo) continue;
                // Keep input face indices until the existing export compacts them.
                mesh.F[f][0] = mesh.F[f][1] = mesh.F[f][2] = -1;
            }
        }
    }
    if (infolevel > 0)
        for (int faceID : mergedFaceIDs) meshLogger->warn("Merge thin layers : {}", faceID);
}

// for adaptive mesh
void DT::buildTetInfo(Mesh& mesh, Args& args) {
    MeshStageLog stageLog(*this, "Rebuild topology", 1);
	// add Nodes
	buildPntInfo(mesh);

    // Retain the exact volume calculation for thin-layer detection.
    std::vector<double> tetVolumes;
    tetVolumes.reserve(mesh.T.size());
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
        if (ori < 0) {
            // Match the predicate evaluation on the corrected vertex order exactly.
            const auto& tet = Elems[newE];
            ori = dt::GEOM_FUNC::orient3d(Nodes[tet.form[0]].pt, Nodes[tet.form[1]].pt,
                Nodes[tet.form[3]].pt, Nodes[tet.form[2]].pt);
        }
        tetVolumes.push_back(ori / 6.0);
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

    mergeThinLayers(mesh, Tri, tetVolumes);
    std::vector<double>().swap(tetVolumes);

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