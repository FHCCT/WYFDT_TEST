#include "mesh_io.h"
#include <array>
#include <cstring>

// Seperate string origin by given a set of patterns.
std::vector<std::string> seperate_string(std::string origin) {
	std::vector<std::string> result;
	std::stringstream ss(origin);
	while (ss >> origin) result.push_back(origin);
	return result;
}
int dt::readVTK(std::string& filename, Mesh& mesh)
{
	auto& V = mesh.V;
	auto& S = mesh.S;
	auto& F = mesh.F;
	auto& T = mesh.T;

	int nPoints = 0;
	int nE = 0;
	double vtkVersion = 0.0;
	std::vector<int> order;
	std::ifstream vtk_file(filename);

	if (!vtk_file.is_open()) {
		spdlog::error("No such file: {}", filename);
		return 0;
	}

	std::string line;
	std::string vtk_type_str = "POLYDATA ";

	// 读取第一行，例如：# vtk DataFile Version 5.1
	if (!std::getline(vtk_file, line)) return 0;

	const std::string& versionLine = line;
	std::size_t versionPos = versionLine.find("Version");

	if (versionPos != std::string::npos) {
		std::stringstream ss(versionLine.substr(versionPos + 7));
		ss >> vtkVersion;
	}

	const bool vtk51 = vtkVersion > 5.0;

	while (std::getline(vtk_file, line)) {

		if (line.length() < 2 || line[0] == '#') continue;

		if (line.find("DATASET") != std::string::npos) {
			std::vector<std::string> words = seperate_string(line);

			if (words[1] == "POLYDATA") vtk_type_str = "POLYGONS ";
			else if (words[1] == "UNSTRUCTURED_GRID") vtk_type_str = "CELLS ";
			else {
				spdlog::error("The format of VTK file is illegal, No clear DATASET name. - {}", filename);
				return 0;
			}
		}

		if (line.find("POINTS ") != std::string::npos) {
			std::vector<std::string> words = seperate_string(line);
			nPoints = std::stoi(words[1]);
			V.resize(nPoints);

			for (int i = 0; i < nPoints; ++i) {
				if (!(vtk_file >> V[i][0] >> V[i][1] >> V[i][2])) return 0;
			}
		}

		if (line.find(vtk_type_str) != std::string::npos) {
			std::vector<std::string> words = seperate_string(line);

			// VTK 5.1：OFFSETS + CONNECTIVITY
			if (vtk51) {
				int nOffsets = std::stoi(words[1]);
				int nConnectivity = std::stoi(words[2]);
				std::string keyword;
				std::string dataType;

				if (!(vtk_file >> keyword >> dataType) || keyword != "OFFSETS") {
					spdlog::error("Invalid VTK 5.1 OFFSETS section.");
					return 0;
				}

				std::vector<long long> offsets(nOffsets);

				for (int i = 0; i < nOffsets; ++i) {
					if (!(vtk_file >> offsets[i])) return 0;
				}

				if (!(vtk_file >> keyword >> dataType) || keyword != "CONNECTIVITY") {
					spdlog::error("Invalid VTK 5.1 CONNECTIVITY section.");
					return 0;
				}

				std::vector<long long> connectivity(nConnectivity);

				for (int i = 0; i < nConnectivity; ++i) {
					if (!(vtk_file >> connectivity[i])) return 0;
				}

				for (int i = 0; i < nOffsets - 1; ++i) {
					long long begin = offsets[i];
					int nindices = static_cast<int>(offsets[i + 1] - offsets[i]);
					order.push_back(nindices);

					if (nindices == 2) {
						std::array<int, 3> tmp = { 0 };
						tmp[0] = static_cast<int>(connectivity[begin]);
						tmp[1] = static_cast<int>(connectivity[begin + 1]);
						S.push_back(tmp);
						++nE;
					}
					else if (nindices == 3) {
						std::array<int, 4> tmp = { 0 };
						tmp[0] = static_cast<int>(connectivity[begin]);
						tmp[1] = static_cast<int>(connectivity[begin + 1]);
						tmp[2] = static_cast<int>(connectivity[begin + 2]);
						F.push_back(tmp);
						++nE;
					}
					else if (nindices == 4) {
						std::array<int, 5> tmp = { 0 };
						tmp[0] = static_cast<int>(connectivity[begin]);
						tmp[1] = static_cast<int>(connectivity[begin + 1]);
						tmp[2] = static_cast<int>(connectivity[begin + 2]);
						tmp[3] = static_cast<int>(connectivity[begin + 3]);
						T.push_back(tmp);
						++nE;
					}
				}
			}

			// VTK 5.0 及以前：原读取逻辑
			else {
				int a = std::stoi(words[1]);

				for (int i = 0; i < a; ++i) {
					int nindices;
					vtk_file >> nindices;
					order.push_back(nindices);

					if (nindices == 2) {
						std::array<int, 3> tmp = { 0 };
						vtk_file >> tmp[0] >> tmp[1];
						S.push_back(tmp);
						++nE;
					}
					else if (nindices == 3) {
						std::array<int, 4> tmp = { 0 };
						vtk_file >> tmp[0] >> tmp[1] >> tmp[2];
						F.push_back(tmp);
						++nE;
					}
					else if (nindices == 4) {
						std::array<int, 5> tmp = { 0 };
						vtk_file >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3];
						T.push_back(tmp);
						++nE;
					}
				}
			}
		}

		if (line.find("CELL_DATA ") != std::string::npos) {
			std::getline(vtk_file, line);
			std::getline(vtk_file, line);

			int nS = 0;
			int nF = 0;
			int nT = 0;

			for (int i = 0; i < nE; ++i) {
				if (order[i] == 2) vtk_file >> S[nS++][2];
				else if (order[i] == 3) vtk_file >> F[nF++][3];
				else if (order[i] == 4) vtk_file >> T[nT++][4];
			}

			break;
		}
	}

	return 1;
}
//int dt::readVTK(std::string& filename, Mesh& mesh) //ignore now
//{
//	auto& V = mesh.V;
//	auto& S = mesh.S;
//	auto& F = mesh.F;
//	auto& T = mesh.T;
//	int nPoints = 0;
//	int nE = 0;
//	std::vector<int> order;
//	std::ifstream vtk_file;
//	vtk_file.open(filename);
//	if (!vtk_file.is_open()) {
//		std::cout << "No such file: " << filename << std::endl;
//		return 0;
//	}
//	std::string vtk_type_str = "POLYDATA ";
//	char buffer[BUFFER_LENGTH];
//	while (!vtk_file.eof()) {
//		vtk_file.getline(buffer, BUFFER_LENGTH);
//		std::string line = (std::string)buffer;
//		if (line.length() < 2 || buffer[0] == '#') continue;
//		if (line.find("DATASET") != std::string::npos) {
//			std::vector<std::string> words = seperate_string(line);
//			if (words[1] == "POLYDATA")
//				vtk_type_str = "POLYGONS ";
//			else if (words[1] == "UNSTRUCTURED_GRID")
//				vtk_type_str = "CELLS ";
//			else {
//				std::cout
//					<< "The format of VTK file is illegal, No clear DATASET name. - "
//					<< filename << std::endl;
//			}
//		}
//		if (line.find("POINTS ") != std::string::npos) {
//			std::vector<std::string> words = seperate_string(line);
//			nPoints = stoi(words[1]);
//			V.resize(nPoints);
//			for (int i = 0; i < nPoints; i++) {
//				if (!(vtk_file >> V[i][0] >> V[i][1] >> V[i][2]))
//					return 0;
//			}
//		}
//		if (line.find(vtk_type_str) != std::string::npos) {
//			std::vector<std::string> words = seperate_string(line);
//			int a = stoi(words[1]);
//			for (int i = 0; i < a; i++) {
//				int nindices;
//				vtk_file >> nindices;
//				order.push_back(nindices);
//				if (nindices == 2) {
//					std::array<int, 3> tmp = { 0 };
//					vtk_file >> tmp[0] >> tmp[1];
//					S.push_back(tmp);
//					nE++;
//				}
//				else if (nindices == 3) {
//					std::array<int, 4> tmp = { 0 };
//					vtk_file >> tmp[0] >> tmp[1] >> tmp[2];
//					F.push_back(tmp);
//					nE++;
//				}
//				else if (nindices == 4) {
//					std::array<int, 5> tmp;
//					vtk_file >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3];
//					T.push_back(tmp);
//					nE++;
//				}
//			}
//		}
//		if (line.find("CELL_DATA ") != std::string::npos) {
//			vtk_file.getline(buffer, BUFFER_LENGTH);
//			vtk_file.getline(buffer, BUFFER_LENGTH);
//			int nS = 0, nF = 0, nT = 0;
//			for (int i = 0; i < nE; i++) {
//				if (order[i] == 2) vtk_file >> S[nS++][2];
//				else if (order[i] == 3) vtk_file >> F[nF++][3];
//				else if (order[i] == 4) vtk_file >> T[nT++][4];
//			}
//			break;
//		}
//	}
//	vtk_file.close();
//	return 1;
//}
// Implement readOBJ
int dt::readOBJ(std::string& filename, Mesh& mesh) {
	std::ifstream obj_file(filename);
	if (!obj_file.is_open()) {
		spdlog::error("Could not open file: {}", filename);
		return 0;
	}

	std::string line;
	while (std::getline(obj_file, line)) {
		// Ignore empty lines and comment lines
		if (line.empty() || line[0] == '#')
			continue;

		std::vector<std::string> tokens = seperate_string(line);
		if (tokens.empty())
			continue;

		// Vertex definition
		if (tokens[0] == "v") {
			if (tokens.size() < 4) {
				spdlog::error("Incomplete vertex definition: {}", line);
				return 0;
			}
			std::array<double, 3> vertex = { 0 };
			try {
				vertex[0] = std::stod(tokens[1]);
				vertex[1] = std::stod(tokens[2]);
				vertex[2] = std::stod(tokens[3]);
			}
			catch (const std::invalid_argument&) {
				spdlog::error("Invalid vertex coordinate: {}", line);
				return 0;
			}
			mesh.V.push_back(vertex);
		}
		// Face definition
		else if (tokens[0] == "f") {
			// Only process triangles here
			// If there's a polygon with more than 3 vertices, you'd need to handle that
			if (tokens.size() < 4) {
				spdlog::error("Incomplete face definition: {}", line);
				return 0;
			}
			std::array<int, 4> face = { 0 };
			try {
				// In OBJ files, indices start from 1, so we convert to 0-based.
				for (int i = 0; i < 3; ++i) {
					std::string vertex_str = tokens[i + 1];
					// Handle potential texture/normal data (format "f v/t/n")
					size_t pos = vertex_str.find('/');
					if (pos != std::string::npos) {
						vertex_str = vertex_str.substr(0, pos);
					}
					face[i] = std::stoi(vertex_str) - 1;
				}
			}
			catch (const std::invalid_argument&) {
				spdlog::error("Invalid face index: {}", line);
				return 0;
			}
			mesh.F.push_back(face);
		}
		// You can handle other types of data here if needed
	}

	obj_file.close();
	return 1;
}

int dt::readMesh(std::string& filename, Mesh& mesh) {
	// Find the last '.' in the filename
	size_t dot_pos = filename.find_last_of('.');

	if (dot_pos == std::string::npos) {
		spdlog::error("Could not determine the file extension for: {}", filename);
		return 0;
	}

	// Extract extension
	std::string extension = filename.substr(dot_pos);

	// Convert extension to lowercase
	std::transform(
		extension.begin(),
		extension.end(),
		extension.begin(),
		[](unsigned char c) { return std::tolower(c); }
	);

	if (extension == ".vtk") {

		// readVTK currently takes a non-const std::string&
		std::string filename_non_const = filename;
		return readVTK(filename_non_const, mesh);
	}
	else if (extension == ".obj") {

		return readOBJ(filename, mesh);
	}
	else {
		spdlog::error(
			"[readMesh] Unsupported file extension: {} for file: {}",
			extension,
			filename
		);

		return 0;
	}
}

namespace {
// Keep the existing ASCII format; batch integer conversion and stream writes.
// Each invocation owns its stream and buffers, so there is no shared state.
class VtkIntegerBuffer {
    FILE* file;
    std::array<char, 65536> buffer;
    size_t used = 0;
public:
    explicit VtkIntegerBuffer(FILE* out) : file(out) {}
    ~VtkIntegerBuffer() { flush(); }
    void flush() {
        if (used) { fwrite(buffer.data(), 1, used, file); used = 0; }
    }
    void number(int value, char separator) {
        if (buffer.size() - used < 32) flush();
        // Unsigned subtraction also handles INT_MIN without signed overflow.
        unsigned magnitude = static_cast<unsigned>(value);
        if (value < 0) { buffer[used++] = '-'; magnitude = 0u - magnitude; }
        char digits[32];
        char* end = digits + sizeof(digits), *begin = end;
        do { *--begin = static_cast<char>('0' + magnitude % 10); magnitude /= 10; } while (magnitude);
        const size_t length = static_cast<size_t>(end - begin);
        std::memcpy(buffer.data() + used, begin, length); used += length;
        buffer[used++] = separator;
    }
    void newline() { buffer[used++] = '\n'; }
};

template<size_t N, bool trailingSpace = false>
void writeVtkCells(FILE* out, const std::vector<std::array<int, N>>& cells) {
    VtkIntegerBuffer buffer(out);
    for (const auto& cell : cells) {
        buffer.number(static_cast<int>(N - 1), ' ');
        for (size_t j = 0; j < N - 1; ++j)
            buffer.number(cell[j], j == N - 2 && !trailingSpace ? '\n' : ' ');
        if (trailingSpace) buffer.newline();
    }
}

template<size_t N>
void writeVtkGeo(FILE* out, const std::vector<std::array<int, N>>& cells) {
    VtkIntegerBuffer buffer(out);
    for (const auto& cell : cells) buffer.number(cell[N - 1], '\n');
}

void writeVtkTypes(FILE* out, int type, size_t count) {
    std::array<char, 12288> block;
    const size_t width = static_cast<size_t>(snprintf(block.data(), block.size(), "%d\n", type));
    const size_t records = block.size() / width;
    for (size_t i = 1; i < records; ++i) std::memcpy(block.data() + i * width, block.data(), width);
    while (count) {
        const size_t n = std::min(count, records);
        fwrite(block.data(), width, n, out); count -= n;
    }
}
} // namespace

int dt::writeVTK(std::string& filename, Mesh& mesh, bool addSurTri)
{
	auto& V = mesh.V;
	auto& S = mesh.S;
	auto& F = mesh.F;
	auto& T = mesh.T;
	auto& PS = mesh.pointSize;
	std::array<char, 65536> fileBuffer;
	FILE* outFile = fopen(filename.c_str(), "w");
	if (outFile == nullptr) {
		spdlog::error("Write VTK file failed. - {}", filename);
		return -1;
	}

	setvbuf(outFile, fileBuffer.data(), _IOFBF, fileBuffer.size());

	fprintf(outFile, "# vtk DataFile Version 2.0\n");
	fprintf(outFile, "TetWild Mesh\n");
	fprintf(outFile, "ASCII\n");
	fprintf(outFile, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(outFile, "POINTS %d double\n", (int)V.size());
	for (int i = 0; i < V.size(); i++)
		fprintf(outFile, "%.16lf %.16lf %.16lf\n", V[i][0], V[i][1], V[i][2]);

	if (addSurTri) {
		fprintf(outFile, "CELLS %ld %ld\n", S.size() + F.size() + T.size(), S.size() * 3 + F.size() * 4 + T.size() * 5);
		writeVtkCells<3, true>(outFile, S);
		writeVtkCells(outFile, F);
	}
	else {
		fprintf(outFile, "CELLS %ld %ld\n", T.size(), T.size() * 5);
	}
	writeVtkCells(outFile, T);

	if (addSurTri) {
		fprintf(outFile, "CELL_TYPES  %ld\n", S.size() + F.size() + T.size());
		writeVtkTypes(outFile, 3, S.size());
		writeVtkTypes(outFile, 5, F.size());
	}
	else {
		fprintf(outFile, "CELL_TYPES  %ld\n", T.size());
	}
	writeVtkTypes(outFile, 10, T.size());

	if (PS.size() != 0) {
		fprintf(outFile, "POINT_DATA  %ld\n", PS.size());
		fprintf(outFile, "SCALARS PointSize float 1\n");
		fprintf(outFile, "LOOKUP_TABLE default\n");
		for (int i = 0; i < PS.size(); i++)
			fprintf(outFile, "%lf\n", PS[i]);
	}

	//out geominfo
	if (addSurTri) {
		fprintf(outFile, "CELL_DATA  %ld\n",  S.size() + F.size() + T.size());
	}
	else {
		fprintf(outFile, "CELL_DATA  %ld\n",   T.size());
	}
	fprintf(outFile, "SCALARS GeoId float 1\n");
	fprintf(outFile, "LOOKUP_TABLE default\n");

	if (addSurTri) {
		writeVtkGeo(outFile, S);
		writeVtkGeo(outFile, F);
	}
	writeVtkGeo(outFile, T);

	fclose(outFile);
	return 1;
}

int dt::readVertex(std::string& filename, std::vector<std::array<double, 4>>& addVertex) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		spdlog::error("无法打开文件: {}", filename);
		return 0;
	}

	size_t pointCount;
	file >> pointCount; // 读取点的个数

	addVertex.resize(pointCount); // 调整 vector 大小

	for (size_t i = 0; i < pointCount; ++i) {
		for (size_t j = 0; j < 4; ++j) {
			file >> addVertex[i][j]; // 读取每个点的4个 double 值
		}
	}

	file.close(); // 关闭文件

	return 1;
}

int dt::readRefineT(std::string& filename, std::vector<int>& refine_tet_id) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		spdlog::error("无法打开文件: {}", filename);
		return 0;
	}

	size_t pointCount;
	file >> pointCount; // 读取点的个数

	refine_tet_id.resize(pointCount); // 调整 vector 大小

	for (size_t i = 0; i < pointCount; ++i) {
		file >> refine_tet_id[i];
	}

	file.close(); // 关闭文件

	return 1;
}
int dt::readNodesSize(std::string& filename, std::vector<double>& NodeSize) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		spdlog::error("无法打开文件: {}", filename);
		return 0;
	}

	size_t pointCount;
	file >> pointCount; // 读取点的个数

	NodeSize.resize(pointCount); // 调整 vector 大小

	for (size_t i = 0; i < pointCount; ++i) {
		file >> NodeSize[i];
	}

	file.close(); // 关闭文件

	return 1;
}
int dt::readPeriodicP(std::string& filename, std::vector<int>& PeriodicP) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		spdlog::error("无法打开文件: {}", filename);
		return 0;
	}

	size_t pointCount;
	file >> pointCount; // 读取点的个数

	for (size_t i = 0; i < pointCount; ++i) {
		int a, b;
		file >> a >>b;
		PeriodicP.push_back(a);
		PeriodicP.push_back(b);
	}

	file.close(); // 关闭文件

	return 1;
}
int dt::readAniSol(std::string& filename, std::vector<std::array<double, 6>>& anisol) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		spdlog::error("无法打开文件: {}", filename);
		return 0;
	}

	size_t pointCount;
	file >> pointCount; // 读取点的个数

	anisol.resize(pointCount); // 调整 vector 大小

	for (size_t i = 0; i < pointCount; ++i) {
		file >> anisol[i][0] >> anisol[i][1] >> anisol[i][2] >> anisol[i][3] >> anisol[i][4] >> anisol[i][5];
	}

	file.close(); // 关闭文件

	return 1;
}