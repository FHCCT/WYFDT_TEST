#ifndef mesh_io_h
#define mesh_io_h

#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iostream>
#include "dt.h"

namespace dt {
	int readMesh(std::string& filename, dt::Mesh& mesh);
	int readVTK(std::string& filename, dt::Mesh& mesh);
	int readOBJ(std::string& filename, dt::Mesh& mesh);
	int writeVTK(std::string& filename, dt::Mesh& mesh, bool addSurTri = false);
	int readVertex(std::string& filename, std::vector<std::array<double, 4>>& addVertex);
	int readRefineT(std::string& filename, std::vector<int>& refine_tet_id);
	int readNodesSize(std::string& filename, std::vector<double>& NodesPoint);
	int readPeriodicP(std::string& filename, std::vector<int>& PeriodicP);
	int readAniSol(std::string& filename, std::vector<std::array<double, 6>>& anisol);
}//namespace MeshIO

#endif
