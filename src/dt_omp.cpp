#include "dt.h"

// Select candidates first; retain incidence only for nodes that will be colored.
// The second scan still includes GOOD cells, including disconnected vertex stars.
void DT::colorBadQualityNodes(std::vector<std::vector<int>>& colors, double improve_goal) {
    colors.clear();
    std::vector<int> candidateIndex(Nodes.size(), -1), candidates;
    const int nElems = static_cast<int>(Elems.size());
    for (int t = 0; t < nElems; ++t) {
        if (isDelEle(t) || isvirtualtet(t) || (ghost != -1 && ishulltet(t))) continue;
        if (!(improve_goal >= 1 || Elems[t].q < improve_goal)) continue;
        for (int n : Elems[t].form)
            if (n >= 0 && n < static_cast<int>(candidateIndex.size())) candidateIndex[n] = 0;
    }
    for (int n = 0; n < static_cast<int>(Nodes.size()); ++n) {
        if (candidateIndex[n] < 0) continue;
        candidateIndex[n] = -1;
        if (n == ghost || isDelNod(n) || isbndpnt(n) || isCornerpnt(n) || lockV.count(n) ||
            periodic_P.count(n) || (addBoxFlag && n > ghost && n < ghost + 9)) continue;
        candidateIndex[n] = static_cast<int>(candidates.size());
        candidates.push_back(n);
    }
    if (candidates.empty()) return;
    std::vector<std::vector<int>> incident(candidates.size());
    for (int t = 0; t < nElems; ++t) {
        if (isDelEle(t) || isvirtualtet(t) || (ghost != -1 && ishulltet(t))) continue;
        for (int n : Elems[t].form) {
            if (n < 0 || n >= static_cast<int>(candidateIndex.size())) continue;
            const int index = candidateIndex[n];
            if (index >= 0) incident[index].push_back(t);
        }
    }
    std::vector<int> assigned(candidates.size(), -1), used;
    for (int index = 0; index < static_cast<int>(candidates.size()); ++index) {
        used.resize(colors.size(), -1);
        for (int t : incident[index]) for (int neighbor : Elems[t].form) {
            if (neighbor < 0 || neighbor >= static_cast<int>(candidateIndex.size())) continue;
            const int other = candidateIndex[neighbor];
            if (other >= 0 && assigned[other] >= 0) used[assigned[other]] = index;
        }
        int color = 0;
        while (color < static_cast<int>(used.size()) && used[color] == index) ++color;
        if (color == static_cast<int>(colors.size())) colors.emplace_back();
        assigned[index] = color;
        colors[color].push_back(candidates[index]);
    }
}

//obtain nodes to smooth
void DT::evalNodesToSmooth(std::vector<int>& nodesToSmooth, double improve_goal) {
	if (improve_goal >= 1) {//global
		nodesToSmooth.resize(Nodes.size());
		int nNodesToSmooth = 0;
		for (int i = 0; i < Nodes.size(); i++) {//from  nSurNodes
			if (isDelNod(i)) continue;
			//if (isbndpnt(i)) continue; //boundary modification
			if (i == ghost) continue;
			if (addBoxFlag && i > ghost && i < ghost + 9) continue;
			nodesToSmooth[nNodesToSmooth++] = i;
		}
		nodesToSmooth.resize(nNodesToSmooth);
	}
	else {
		const int nNodes = Nodes.size();
		const int nElems = Elems.size();
		std::vector<bool> isIllNode(nNodes);
		for (int i = 0; i < nElems; i++) {
			if (isDelEle(i) || isvirtualtet(i))
				continue;
			if (Elems[i].q < improve_goal) {
				for (int j = 0; j < 4; j++)
					isIllNode[Elems[i].form[j]] = true;
			}
		}

		nodesToSmooth.reserve(nNodes);
		for (int i = 0; i < nNodes; i++) {
			if (isIllNode[i]) {
				nodesToSmooth.emplace_back(i);
			}
		}
	}
	return;
}

int DT::addElem(int thread_n, bool UseVacancy) {
	int newE = -1;

	if (UseVacancy && !Evacancy_thread[thread_n].empty()) {
		newE = Evacancy_thread[thread_n].front();
		Evacancy_thread[thread_n].pop();
		Elems[newE].form[0] = 0;
		Elems[newE].form[1] = 0;
		Elems[newE].form[2] = 0;
		Elems[newE].form[3] = 0;
		Elems[newE].neig[0] = -1;
		Elems[newE].neig[1] = -1;
		Elems[newE].neig[2] = -1;
		Elems[newE].neig[3] = -1;
		Elems[newE].info = 0;
		Elems[newE].q = -1;
	}
	else {
        if (parallelTopologyBatch)
            throw std::runtime_error("Topology batch exhausted preallocated element slots");
#pragma omp critical
		{
			newE = Elems.size();
			Elems.push_back(Elem());
		}
	}
	return newE;
}

int DT::addElem(int pa, int pb, int pc, int pd, int thread_n, bool UseVacancy) {
	int newE = -1;

	if (UseVacancy && !Evacancy_thread[thread_n].empty()) {
		newE = Evacancy_thread[thread_n].front();
		Evacancy_thread[thread_n].pop();
		Elems[newE].form[0] = pa;
		Elems[newE].form[1] = pb;
		Elems[newE].form[2] = pc;
		Elems[newE].form[3] = pd;
		Elems[newE].neig[0] = -1;
		Elems[newE].neig[1] = -1;
		Elems[newE].neig[2] = -1;
		Elems[newE].neig[3] = -1;
		Elems[newE].info = 0;
		Elems[newE].q = -1;
	}
	else {
        if (parallelTopologyBatch)
            throw std::runtime_error("Topology batch exhausted preallocated element slots");
#pragma omp critical
		{
			newE = Elems.size();
			Elems.push_back(Elem(pa, pb, pc, pd));
		}
	}
	return newE;
}

int DT::addNode(double x, double y, double z, double space, int thread_n, bool UseVacancy) {
	int newN = -1;
#pragma omp critical
	{
		if (UseVacancy && !Nvacancy.empty()) {
			newN = Nvacancy.front();
			Nvacancy.pop();
			Nodes[newN].pt[0] = x;
			Nodes[newN].pt[1] = y;
			Nodes[newN].pt[2] = z;
			Nodes[newN].space = space;
			Nodes[newN].info = 0;
			Nodes[newN].tet = 0;
			Nodes[newN].type = 0;
		}
		else {
			newN = Nodes.size();
			Nodes.push_back(Node(x, y, z, space));
		}
	}
	Nodes[newN].occupying.store(thread_n);

	return newN;
}

void DT::DelEle(int E, int thread_n) {
	if (isDelEle(E))return;
	Elems[E].form[0] = Elems[E].form[1] = Elems[E].form[2] = Elems[E].form[3] = -1;
	Elems[E].neig[0] = Elems[E].neig[1] = Elems[E].neig[2] = Elems[E].neig[3] = -1;
	Elems[E].info = -1;	Elems[E].geo = -1; Elems[E].q = -1;
	Evacancy_thread[thread_n].push(E);//Collect delete id
	return;
}
