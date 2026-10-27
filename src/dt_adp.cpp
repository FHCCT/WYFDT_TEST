#include "./dt.h"
#include <limits>

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
			splitEdg(Edgid);
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
			int BW_recall = BW_insert_vertex(iNod, nearTets, 1, -1, {{p1, p2}});
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
    const int firstRefinementNode = static_cast<int>(Nodes.size());
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
			splitEdg(Edgid);
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
			int BW_recall = BW_insert_vertex(iNod, nearTets, 1, -1, {{p1, p2}});
			if (BW_recall == 1) {//success
				continue;
			}
			else {//failed
				DelNod(iNod);
			}
		}
	}

    // Failed insertion slots can be reused. Register live points only after
    // all front-end splits, including redirected and matched boundary splits.
    for (int n = firstRefinementNode; n < static_cast<int>(Nodes.size()); ++n)
        if (!isDelNod(n)) refinementProtectedPoints.emplace(n,true);

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

    //mergeThinLayers(mesh, Tri, tetVolumes);
    //std::vector<double>().swap(tetVolumes);

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
    buildRectangularFaces();
	return;
}

namespace dt {
namespace {
using Vec = Eigen::Vector3d;
using Edge = std::pair<int, int>;
constexpr double rectangularSynchronizationWidthRatio = 2.0;
constexpr double rectangularRailMatchWidthRatio = 0.25;
Vec point(const DT& mesh, int node) { return Vec(mesh.Nodes[node].pt); }
Vec vector(const std::array<double, 3>& v) { return Vec(v.data()); }
std::array<double, 3> array(const Vec& v) { return {{v[0], v[1], v[2]}}; }

struct RailLocation {
    int node = -1, edge = -1;
    double fraction = 0;
};

bool plane(DT& mesh, const std::vector<int>& faces, DT::PlanarFace& result) {
    std::set<int> nodes;
    Vec origin = Vec::Zero(), normal = Vec::Zero();
    double maximumAreaSquared = 0;
    for (int f : faces) {
        const auto& v = mesh.SurTris[f].form;
        for (int n : v) nodes.insert(n);
        const Vec a = point(mesh, v[0]);
        const Vec n = (point(mesh, v[1]) - a).cross(point(mesh, v[2]) - a);
        if (n.squaredNorm() > maximumAreaSquared) {
            maximumAreaSquared = n.squaredNorm(); origin = a; normal = n;
        }
    }
    if (!(maximumAreaSquared > 0) || !std::isfinite(maximumAreaSquared)) return false;
    normal.normalize();
    Vec lo = point(mesh, *nodes.begin()), hi = lo;
    double magnitude = 0;
    for (int n : nodes) {
        const Vec p = point(mesh, n);
        if (!p.allFinite()) return false;
        lo = lo.cwiseMin(p); hi = hi.cwiseMax(p);
        magnitude = std::max(magnitude, p.cwiseAbs().maxCoeff());
    }
    const double tolerance = std::max(1e-8 * (hi - lo).norm(),
        128 * std::numeric_limits<double>::epsilon() * magnitude);
    // Fit the plane offset as well as its normal. Input rounding can place
    // vertices on either side: anchoring at an extreme vertex doubles error.
    double minimumOffset = DBL_MAX, maximumOffset = -DBL_MAX;
    for (int n : nodes) {
        const double offset = (point(mesh,n)-origin).dot(normal);
        minimumOffset = std::min(minimumOffset,offset); maximumOffset = std::max(maximumOffset,offset);
    }
    origin += .5*(minimumOffset+maximumOffset)*normal;
    for (int n : nodes)
        if (std::abs((point(mesh, n) - origin).dot(normal)) > tolerance) return false;
    const Vec along = normal.unitOrthogonal(), across = normal.cross(along);
    result.origin = array(origin); result.normal = array(normal);
    result.along = array(along); result.across = array(across); result.tolerance = tolerance;
    return true;
}

RailLocation locateRail(DT& mesh, const DT::RectangularFace& face, int geo, int side, double position) {
    RailLocation result;
    const Vec origin = vector(face.origin), along = vector(face.along), across = vector(face.across);
    int before = -1, after = -1; double lower = -DBL_MAX, upper = DBL_MAX;
    // A transverse station need not agree to floating-point precision. Reuse
    // the closest live station within a quarter width before splitting its rail.
    const double matchTolerance = std::max(face.tolerance, rectangularRailMatchWidthRatio * face.width);
    double nearest = DBL_MAX;
    for (int n : face.longSideNodes[side]) {
        if (mesh.isDelNod(n) || !mesh.isbndpnt(n)) continue;
        const Vec p = point(mesh, n) - origin;
        if (std::abs(p.dot(across) - side * face.width) > face.tolerance ||
            std::abs(p.dot(along.cross(across))) > face.tolerance) continue;
        const double x = p.dot(along);
        if (std::abs(x - position) <= matchTolerance && std::abs(x - position) < nearest) {
            // Recycled slots belonging to another geometry cannot be reused.
            std::unordered_set<int> current;
            mesh.findSphere_tri(n, current);
            for (int f : current) if (!mesh.isDelSurTri(f) && mesh.SurTris[f].parent == geo) {
                result.node = n; nearest = std::abs(x - position); break;
            }
            continue;
        }
        if (x < position && x > lower) { lower = x; before = n; }
        if (x > position && x < upper) { upper = x; after = n; }
    }
    if (result.node >= 0) return result;
    const int* found = before >= 0 && after >= 0 ? mesh.BndEdg.find(before, after) : nullptr;
    if (!found || mesh.isDelSurEdg(*found) || mesh.lockE.count(*found)) return result;
    int uses = 0;
    for (int f : mesh.SurEdgs[*found].face) {
        if (mesh.lockF.count(f)) return result;
        if (!mesh.isDelSurTri(f) && mesh.SurTris[f].parent == geo) ++uses;
    }
    if (uses != 1) return result;
    result.edge = *found;
    result.fraction = (position - lower) / (upper - lower);
    if (mesh.SurEdgs[*found].iStart != before) result.fraction = 1 - result.fraction;
    return result;
}

bool rectangle(DT& mesh, const std::vector<int>& faces, DT::RectangularFace& result) {
    struct Use { int face = -1, count = 0, direction = 0; };
    std::map<Edge, Use> edges;
    std::set<int> nodes;
    std::vector<std::vector<std::pair<int, int>>> adjacent(faces.size());
    for (int i = 0; i < static_cast<int>(faces.size()); ++i) {
        const auto& f = mesh.SurTris[faces[i]];
        for (int k = 0; k < 3; ++k) {
            const int a = f.form[k], b = f.form[(k + 1) % 3];
            if (a == b) return false;
            nodes.insert(a);
            auto& use = edges[std::minmax(a, b)];
            const int direction = a < b ? 1 : -1;
            if (++use.count > 2) return false;
            if (use.face >= 0) {
                const int sign = direction == use.direction ? -1 : 1;
                adjacent[i].push_back({use.face, sign}); adjacent[use.face].push_back({i, sign});
            } else { use.face = i; use.direction = direction; }
        }
    }
    if (faces.empty() || static_cast<int64_t>(nodes.size()) - static_cast<int64_t>(edges.size()) +
        static_cast<int64_t>(faces.size()) != 1)
        return false;
    std::vector<int> signs(faces.size(), 0);
    std::vector<int> connected(1, 0); signs[0] = 1;
    for (size_t i = 0; i < connected.size(); ++i)
        for (const auto& entry : adjacent[connected[i]]) {
            const int n = entry.first, sign = signs[connected[i]] * entry.second;
            if (!signs[n]) { signs[n] = sign; connected.push_back(n); }
            else if (signs[n] != sign) return false;
        }
    if (connected.size() != faces.size()) return false;

    std::map<int, std::vector<int>> boundary;
    size_t boundaryEdges = 0;
    for (const auto& edge : edges) if (edge.second.count == 1) {
        boundary[edge.first.first].push_back(edge.first.second);
        boundary[edge.first.second].push_back(edge.first.first); ++boundaryEdges;
    }
    if (boundary.size() < 4) return false;
    for (const auto& entry : boundary) if (entry.second.size() != 2) return false;
    std::vector<int> loop;
    int previous = -1, current = boundary.begin()->first;
    do {
        if (loop.size() >= boundary.size()) return false;
        loop.push_back(current);
        const auto& next = boundary.at(current);
        const int n = next[0] == previous ? next[1] : next[0];
        previous = current; current = n;
    } while (current != loop.front());
    // Exactly one boundary loop: holes and disconnected components are excluded.
    if (loop.size() != boundary.size() || loop.size() != boundaryEdges) return false;

    Vec lo = point(mesh, loop.front()), hi = lo;
    double magnitude = 0;
    for (int n : nodes) {
        const Vec p = point(mesh, n); lo = lo.cwiseMin(p); hi = hi.cwiseMax(p);
        magnitude = std::max(magnitude, p.cwiseAbs().maxCoeff());
    }
    const double tolerance = std::max((hi - lo).norm() * 1e-8,
        128 * std::numeric_limits<double>::epsilon() * magnitude);
    std::vector<int> corners;
    for (size_t i = 0; i < loop.size(); ++i) {
        const Vec a = point(mesh, loop[(i + loop.size() - 1) % loop.size()]);
        const Vec b = point(mesh, loop[i]), c = point(mesh, loop[(i + 1) % loop.size()]);
        const Vec u = b - a, v = c - b;
        if (u.norm() <= tolerance || v.norm() <= tolerance) return false;
        if (u.dot(v) <= 0 || u.cross(v).norm() > tolerance * (u.norm() + v.norm()))
            corners.push_back(loop[i]);
    }
    if (corners.size() != 4) return false;
    const Vec origin = point(mesh, corners[0]);
    Vec along = point(mesh, corners[1]) - origin;
    Vec across = point(mesh, corners[3]) - origin;
    if ((point(mesh, corners[2]) - origin - along - across).norm() > 4 * tolerance) return false;
    if (std::abs(along.dot(across)) > 1e-7 * along.norm() * across.norm()) return false;
    if (along.squaredNorm() < across.squaredNorm()) std::swap(along, across);
    const double length = along.norm(), width = across.norm();
    if (!(width > 16 * tolerance && length > mesh.rectangularFaceAspectRatio * width)) return false;
    along /= length; across /= width;
    const Vec normal = along.cross(across);
    for (int n : nodes) {
        const Vec p = point(mesh, n) - origin;
        const double x = p.dot(along), y = p.dot(across);
        if (std::abs(p.dot(normal)) > tolerance || x < -tolerance || x > length + tolerance ||
            y < -tolerance || y > width + tolerance) return false;
    }
    for (const auto& edge : edges) if (edge.second.count == 1) {
        const Vec a = point(mesh, edge.first.first) - origin, b = point(mesh, edge.first.second) - origin;
        const double x1 = a.dot(along), x2 = b.dot(along), y1 = a.dot(across), y2 = b.dot(across);
        const bool side = (std::abs(x1) <= tolerance && std::abs(x2) <= tolerance) ||
            (std::abs(x1 - length) <= tolerance && std::abs(x2 - length) <= tolerance) ||
            (std::abs(y1) <= tolerance && std::abs(y2) <= tolerance) ||
            (std::abs(y1 - width) <= tolerance && std::abs(y2 - width) <= tolerance);
        if (!side) return false;
    }
    double area = 0, orientation = 0;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& v = mesh.SurTris[faces[i]].form;
        const double a = (point(mesh, v[1]) - point(mesh, v[0])).cross(
            point(mesh, v[2]) - point(mesh, v[0])).dot(normal) * signs[i];
        if (a == 0 || (orientation != 0 && a * orientation <= 0)) return false;
        orientation = a; area += std::abs(a) * .5;
    }
    if (std::abs(area - length * width) > tolerance * (length + width) * 4) return false;
    result.origin = array(origin); result.along = array(along); result.across = array(across);
    result.length = length; result.width = width; result.tolerance = tolerance;
    for (const auto& entry : boundary) {
        const double y = (point(mesh, entry.first) - origin).dot(across);
        if (std::abs(y) <= tolerance) result.longSideNodes[0].insert(entry.first);
        if (std::abs(y - width) <= tolerance) result.longSideNodes[1].insert(entry.first);
    }
    return true;
}
}

void DT::buildRectangularFaces() {
    rectangularFaces.clear();
    planarFaces.clear();
    std::map<int, std::vector<int>> groups;
    for (int f = 0; f < static_cast<int>(SurTris.size()); ++f)
        if (!isDelSurTri(f)) groups[SurTris[f].parent].push_back(f);
    for (const auto& group : groups) {
        RectangularFace face;
        if (rectangle(*this, group.second, face)) rectangularFaces.emplace(group.first, std::move(face));
        else {
            PlanarFace planar;
            if (plane(*this, group.second, planar)) planarFaces.emplace(group.first, std::move(planar));
        }
    }
    if (infolevel > 1) meshLogger->debug("Cached boundary geometry: rectangles={} otherPlanes={}",
        rectangularFaces.size(), planarFaces.size());
}

int DT::cachedFaceFlipGeometry(int a, int b, int c, int d) {
    const int* edge = BndEdg.find(a, b);
    if (!edge || SurEdgs[*edge].constrain > 0 || SurEdgs[*edge].face.size() != 2) return -1;
    const int f = SurEdgs[*edge].face[0], g = SurEdgs[*edge].face[1];
    if (isDelSurTri(f) || isDelSurTri(g) || SurTris[f].parent != SurTris[g].parent) return -1;
    const auto cached = rectangularFaces.find(SurTris[f].parent);
    const RectangularFace* rectangle = cached == rectangularFaces.end() ? nullptr : &cached->second;
    Vec origin, along, across, normal;
    double tolerance, length = 1, width = 1;
    if (rectangle) {
        origin = vector(rectangle->origin); along = vector(rectangle->along); across = vector(rectangle->across);
        normal = along.cross(across); tolerance = rectangle->tolerance;
        length = rectangle->length; width = rectangle->width;
    } else {
        const auto planar = planarFaces.find(SurTris[f].parent);
        if (planar == planarFaces.end()) return -1;
        const auto& face = planar->second;
        origin = vector(face.origin); along = vector(face.along); across = vector(face.across);
        normal = vector(face.normal); tolerance = face.tolerance;
    }
    const int ids[] = {a,b,c,d};
    long double xy[4][2];
    for (int i = 0; i < 4; ++i) {
        if (ids[i] < 0 || ids[i] >= static_cast<int>(Nodes.size()) || isDelNod(ids[i])) return 0;
        const Vec p = point(*this, ids[i]) - origin;
        const double x = p.dot(along), y = p.dot(across);
        if (!p.allFinite() || std::abs(p.dot(normal)) > tolerance) return 0;
        if (rectangle && (x < -tolerance || x > length + tolerance ||
            y < -tolerance || y > width + tolerance)) return 0;
        xy[i][0] = x / length; xy[i][1] = y / width;
    }
    auto orientation = [&](int i, int j, int k) {
        const long double x = (xy[j][0]-xy[i][0])*(xy[k][1]-xy[i][1]);
        const long double y = (xy[j][1]-xy[i][1])*(xy[k][0]-xy[i][0]);
        const long double value = x-y;
        const long double error = 64*std::numeric_limits<double>::epsilon()*(std::abs(x)+std::abs(y));
        return value > error ? 1 : value < -error ? -1 : 0;
    };
    // Both diagonals must lie inside the quadrilateral. This also excludes
    // folds/overlaps and degenerate replacements in the original plane frame.
    return orientation(0,1,2)*orientation(0,1,3) == -1 &&
        orientation(2,3,0)*orientation(2,3,1) == -1 ? 1 : 0;
}

int DT::redirectRectangularSplit(int edge, int deep, bool retryVolume, double fraction) {
    // -1 means an ordinary split. Once a long interior edge is recognized,
    // failure/no-op must not fall back to inserting another point on it.
    if (!modifyBnd || rectangularPairingInProgress || rectangularFaces.empty() ||
        SurEdgs[edge].constrain > 0 || SurEdgs[edge].face.size() != 2) return -1;
    const int first = SurEdgs[edge].face[0], second = SurEdgs[edge].face[1];
    if (isDelSurTri(first) || isDelSurTri(second)) return -1;
    const int geo = SurTris[first].parent;
    if (SurTris[second].parent != geo) return -1;
    const auto cached = rectangularFaces.find(geo);
    if (cached == rectangularFaces.end()) return -1;
    const auto& face = cached->second;
    const Vec origin = vector(face.origin), along = vector(face.along), across = vector(face.across);
    const Vec a = point(*this, SurEdgs[edge].iStart) - origin;
    const Vec b = point(*this, SurEdgs[edge].iEnd) - origin;
    // A local diagonal must be splittable to replace its fan connectivity.
    // Redirect only diagonals that are still long relative to the strip width.
    if ((b - a).norm() <= rectangularInteriorSplitWidthRatio * face.width + face.tolerance)
        return -1;
    if (std::abs(a.dot(along.cross(across))) > face.tolerance ||
        std::abs(b.dot(along.cross(across))) > face.tolerance) return -1;
    const double position = ((1 - fraction) * a + fraction * b).dot(along);
    if (position <= face.tolerance || position >= face.length - face.tolerance) return -1;
    const RailLocation sides[] = {locateRail(*this, face, geo, 0, position),
        locateRail(*this, face, geo, 1, position)};
    if (sides[0].node >= 0 && sides[1].node >= 0) return 0;
    for (const auto& side : sides) if (side.edge >= 0) {
        const int inserted = splitEdgImpl(side.edge, deep, retryVolume, side.fraction, true);
        if (inserted > 0) {
            ++rectangularRedirectedSplits;
            if (infolevel > 1) meshLogger->debug(
                "Rectangle split redirected: geometry={} diagonal={} railNode={} position={:.12g}",
                geo, edge, inserted, position);
            return inserted;
        }
    }
    ++rectangularPairFailures;
    if (infolevel > 1) meshLogger->debug("Rectangle redirection rejected: geometry={} diagonal={}", geo, edge);
    return 0;
}

void DT::synchronizeRectangularSplit(int node, double sourceLength) {
    if (node < 0 || isDelNod(node)) return;
    std::unordered_set<int> incident;
    findSphere_tri(node, incident);
    std::set<int> geometries;
    for (int f : incident) if (!isDelSurTri(f) && rectangularFaces.count(SurTris[f].parent))
        geometries.insert(SurTris[f].parent);
    // Register child nodes even during paired/periodic splits. Edge IDs are not
    // cached: later queries resolve current edges after splitting or collapse.
    for (int geo : geometries) {
        auto& face = rectangularFaces.at(geo);
        const double y = (point(*this, node) - vector(face.origin)).dot(vector(face.across));
        if (std::abs(y) <= face.tolerance) face.longSideNodes[0].insert(node);
        if (std::abs(y - face.width) <= face.tolerance) face.longSideNodes[1].insert(node);
    }
    if (!modifyBnd) return;
    pendingBoundaryPairs.push_back(node);
    if (rectangularPairingInProgress) return;
    struct Scope {
        bool& active;
        explicit Scope(bool& value) : active(value) { active = true; }
        ~Scope() { active = false; }
    } scope(rectangularPairingInProgress);
    for (int geo : geometries) {
        auto& face = rectangularFaces.at(geo);
        // A large rectangle can already contain a fine local triangulation.
        // Synchronize only an actual split that is long relative to its width;
        // local diagonals and short rail segments must not refine remote sides.
        if (!(sourceLength > rectangularSynchronizationWidthRatio * face.width + face.tolerance)) continue;
        const Vec origin = vector(face.origin), along = vector(face.along);
        const double position = (point(*this, node) - origin).dot(along);
        if (position <= face.tolerance || position >= face.length - face.tolerance) continue;
        for (int side = 0; side < 2; ++side) {
            const auto target = locateRail(*this, face, geo, side, position);
            if (target.node >= 0) continue;
            if (target.edge < 0) {
                ++rectangularPairFailures;
                if (infolevel > 1) meshLogger->debug("Rectangle pair unavailable: geometry={} node={} side={}", geo, node, side);
                continue;
            }
            struct Acceptance {
                BWAcceptance& value; BWAcceptance saved;
                explicit Acceptance(BWAcceptance& v) : value(v), saved(v) { value = BWAcceptance::VolumeOnly; }
                ~Acceptance() { value = saved; }
            } acceptance(bwAcceptance);
            const int paired = splitEdgImpl(target.edge, 0, false, target.fraction, true);
            if (paired > 0) {
                ++rectangularPairedSplits;
                if (infolevel > 1) meshLogger->debug("Rectangle pair inserted: geometry={} node={} paired={} fraction={:.12g}", geo, node, paired, target.fraction);
            } else {
                ++rectangularPairFailures;
                if (infolevel > 1) meshLogger->debug("Rectangle pair rejected by volume/insertion checks: geometry={} node={}", geo, node);
            }
        }
    }
    // Paired points can border another, nonrectangular geometry. Process that
    // local cross-section too, without recursive splitting around a face chain.
    std::set<int> pairedGeometries;
    const size_t originalPairs = pendingBoundaryPairs.size();
    for (size_t i = 0; i < originalPairs; ++i) {
        const int current = pendingBoundaryPairs[i];
        if (!qualityBoundarySplitInProgress || isDelNod(current)) continue;
        std::unordered_set<int> star;
        findSphere_tri(current,star);
        std::set<int> local;
        for (int f : star) if (!isDelSurTri(f) && !rectangularFaces.count(SurTris[f].parent))
            local.insert(SurTris[f].parent);
        for (int geo : local) if (pairedGeometries.insert(geo).second) {
            const int paired = pairLocalStrip(current,geo);
            if (paired > 0) ++localStripPairedSplits;
        }
    }
    pendingBoundaryPairs.clear();
}

int DT::pairLocalStrip(int node, int geo) {
    // The source must lie inside a straight feature rail, not at a corner.
    std::unordered_set<int> star;
    findSphere_tri(node,star);
    std::set<int> rail;
    Vec normal = Vec::Zero(); double area = 0;
    for (int f : star) if (!isDelSurTri(f) && SurTris[f].parent == geo) {
        const auto& v = SurTris[f].form;
        const Vec n = (point(*this,v[1])-point(*this,v[0])).cross(point(*this,v[2])-point(*this,v[0]));
        if (n.squaredNorm() > area) { normal = n; area = n.squaredNorm(); }
        for (int other : v) if (other != node) {
            const int* e = BndEdg.find(node,other);
            if (e && SurEdgs[*e].constrain > 0) rail.insert(other);
        }
    }
    if (rail.size() != 2 || !(area > 0)) return 0;
    const Vec origin = point(*this,node);
    auto it = rail.begin(); const Vec left = point(*this,*it++)-origin, right = point(*this,*it)-origin;
    const double span = left.norm()+right.norm();
    if (!(left.norm()>0 && right.norm()>0) || left.normalized().dot(right.normalized()) > -.999999) return 0;
    const Vec along = left.normalized(); normal.normalize();
    Vec across = normal.cross(along).normalized();
    // Use a local stable plane: a remote part of the same geometry may bend.
    const double tolerance = std::max(1e-8*span,128*std::numeric_limits<double>::epsilon()*origin.cwiseAbs().maxCoeff());
    bool inward = false;
    for (int f : star) if (!isDelSurTri(f) && SurTris[f].parent == geo)
        for (int n : SurTris[f].form) if (std::abs((point(*this,n)-origin).dot(across))>tolerance) {
            if ((point(*this,n)-origin).dot(across)<0) across = -across;
            inward = true; break;
        }
    if (!inward) return 0;
    std::vector<int> queue;
    for (int f : star) if (!isDelSurTri(f) && SurTris[f].parent == geo) queue.push_back(f);
    std::set<int> visited(queue.begin(),queue.end());
    int target = -1; double fraction = 0, width = DBL_MAX;
    // Walk only triangles intersected by this inward ray. Bound the geometric
    // search by the rail's local support, not by the size of the whole face.
    for (size_t i=0;i<queue.size();++i) {
        const auto& v = SurTris[queue[i]].form;
        for (int n : v) if (std::abs((point(*this,n)-origin).dot(normal))>tolerance) return 0;
        for (int k=0;k<3;++k) {
            const Vec a = point(*this,v[k])-origin, b = point(*this,v[(k+1)%3])-origin;
            const double x = a.dot(along), dx = (b-a).dot(along);
            if (std::abs(dx)<=tolerance) continue;
            const double t = -x/dx;
            if (t < -1e-10 || t > 1+1e-10) continue;
            const double y = (a+t*(b-a)).dot(across);
            if (y<=tolerance || y>span/rectangularFaceAspectRatio || y>width+tolerance) continue;
            const int* e = BndEdg.find(v[k],v[(k+1)%3]);
            if (!e) return 0;
            if (SurEdgs[*e].constrain > 0) {
                if (std::abs((b-a).normalized().dot(along)) < std::cos(5*PI/180)) return 0;
                if (y<width) { target=*e; width=y; fraction=SurEdgs[*e].iStart==v[k]?t:1-t; }
            } else {
                for (int f : SurEdgs[*e].face) if (!isDelSurTri(f) && SurTris[f].parent==geo && visited.insert(f).second)
                    queue.push_back(f);
            }
        }
    }
    if (target<0 || lockE.count(target)) return 0;
    for (int f : SurEdgs[target].face) if (lockF.count(f)) return 0;
    const double length = distance(Nodes[SurEdgs[target].iStart].pt,Nodes[SurEdgs[target].iEnd].pt);
    // Reuse an existing transverse match; do not create an almost coincident point.
    if (std::min(fraction,1-fraction)*length<=std::max(tolerance,.5*width)) return 0;
    struct Acceptance {
        BWAcceptance& value; BWAcceptance saved;
        ~Acceptance() { value=saved; }
    } acceptance{bwAcceptance,bwAcceptance};
    bwAcceptance=BWAcceptance::VolumeOnly;
    const int paired=splitEdgImpl(target,0,false,fraction,true);
    if (paired>0 && infolevel>1) meshLogger->debug("Local strip pair: geometry={} node={} paired={} width={:.12g}",geo,node,paired,width);
    return paired;
}
}
