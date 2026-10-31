#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace dt {
class DT;

enum class BWStatus { Ready, Inserted, Duplicate, Rejected, RetryConflict, Stale };

// Volume-only insertion skips shape thresholds; callers may add local volume floors.
enum class BWAcceptance { Default, VolumeOnly };

// A request owns candidate data; node == -1 denotes an unallocated point.
struct BWRequest {
    int node = -1;
    std::array<double, 3> point{};
    double space = 0;
    std::array<double, 6> metric{};
    bool anisotropic = false;
    bool boundary = false;
    bool trackAccess = true;
    // The initial serial Delaunay pass has no threaded occupation check.
    bool collectLockNodes = true;
    // Refinement only: commit writes the new node; finishBW serially updates old-node P2T.
    bool deferP2T = false;
    int info = 0;
    BWAcceptance acceptance = BWAcceptance::Default;
    int initShell = 0;
    std::vector<int> seeds;
    // Source segment identity belongs to this insertion, including frozen workers.
    // Initially a midpoint; bounded searches may move it inside the same segment.
    // A face containing both endpoints contains every such point exactly.
    std::array<int, 2> midpointEdge{{-1, -1}};
    // Matched rectangle/periodic points must retain their requested position.
    bool fixedEdgePoint = false;
    // A point constructed inside a source face is exactly on that face,
    // independent of cancellation in the orientation predicate.
    std::array<int, 3> sourceFace{{-1, -1, -1}};
    bool midpointOnFace(int a, int b, int c) const {
        if (sourceFace[0] >= 0 &&
            std::find(sourceFace.begin(), sourceFace.end(), a) != sourceFace.end() &&
            std::find(sourceFace.begin(), sourceFace.end(), b) != sourceFace.end() &&
            std::find(sourceFace.begin(), sourceFace.end(), c) != sourceFace.end()) return true;
        if (midpointEdge[0] < 0 || midpointEdge[1] < 0 || midpointEdge[0] == midpointEdge[1]) return false;
        return (a == midpointEdge[0] || b == midpointEdge[0] || c == midpointEdge[0]) &&
            (a == midpointEdge[1] || b == midpointEdge[1] || c == midpointEdge[1]);
    }
    // Optional borrowed membership from a frozen BndTri table, for info == 2.
    // The owner must outlive all plans; never retain it in asynchronous work.
    const std::vector<uint8_t>* boundaryNodes = nullptr;
};

struct BWScratchTable {
    struct Slot { int key = -1; int64_t value = 0; };
    std::array<Slot, 128> small;
    std::vector<Slot> large;
    size_t count = 0;
    int64_t& operator[](int key);
    void clear();
    // Dense entity IDs in the serial input pass avoid repeated hashing.
    // Resize before planning, while no references to entries are held.
    void index(size_t size);
    // Initial serial insertion only: borrow mesh flags during planning and
    // restore them with clear() before allocating or committing any cells.
    void borrowInitialFlags(DT& mesh);
private:
    DT* initialMesh = nullptr;
    std::vector<std::pair<int, int64_t>> savedFlags;
    std::vector<int64_t> indexedValues;
    std::vector<uint32_t> indexedStamps;
    uint32_t indexedEpoch = 0;
    void grow();
};

struct BWEdgeTable {
    struct Slot {
        uint64_t key = UINT64_MAX;
        int face = -1;
        uint8_t side = 0;
        bool paired = false;
    };
    std::array<Slot, 64> small;
    std::vector<Slot> large;
    size_t count = 0;
    void clear();
    Slot& insert(uint64_t key, int face, int side, bool& inserted);
};

struct BWFace {
    std::array<int, 3> vertices{};
    int outside = -1;
    int outsideFace = -1;
    int geo = -1;
    // Faces 1..3 refer to other new tetrahedra using local face indices.
    std::array<int, 3> adjacent{{-1, -1, -1}};
    std::array<int, 3> adjacentFace{{-1, -1, -1}};
};

struct BWPlan {
    BWRequest request;
    BWStatus status = BWStatus::Rejected;
    int duplicateNode = -1;
    std::vector<int> readTets;
    std::vector<int> writeTets;
    std::vector<int> writeNodes;
    std::vector<int> cavity;
    std::vector<BWFace> faces;
    std::vector<int> newElements;
    int ghostTet = -1;
    // Scratch is local to one plan and is never stored on mesh entities.
    std::vector<int> workingCavity;
    std::vector<int> boundaryFaces;
    std::vector<int> lockNodes;
    std::vector<uint8_t> adjustQueued;
    std::vector<int> adjustQueue;
    std::vector<int> shell;
    BWScratchTable tetInfo, nodeInfo;
    BWEdgeTable edges;
    std::array<uint32_t, 1024> denseEdges{};
    uint32_t denseEpoch = 0;
    std::vector<std::array<uint8_t, 3>> denseVertices;
};

struct BWSerialWorkspace {
    BWPlan plan;
    std::vector<int> slots;
};

BWRequest makeBWRequest(DT& mesh, int node, const std::vector<int>& seeds, int info);
BWStatus findBWCavity(DT& mesh, BWPlan& plan);
BWStatus adjustBWCavity(DT& mesh, BWPlan& plan);
BWStatus prepareBWFill(DT& mesh, BWPlan& plan, bool checkInitialVisibility = false);
// Internal adjacency builder, also exercised independently by regression tests.
bool pairBWBoundary(BWPlan& plan);
BWStatus planBW(DT& mesh, const BWRequest& request, BWPlan& plan);
// Isotropic edge points only: retry poor/failed candidates at bounded positions.
// Planning stays read-only; commitBW installs the chosen position and size.
BWStatus planBWOnEdge(DT& mesh, const BWRequest& request, BWPlan& plan);

// Caller freezes topology throughout planning, selects independent plans, and
// allocates every node/element slot before any worker enters commitBW.
// commitBW does not grow containers, recycle slots, or update global counters.
void commitBW(DT& mesh, BWPlan& plan, int node, const std::vector<int>& slots);
// Run after all commits have completed. This updates recycling, counters and P2T.
// With deferP2T, serialize finishes and allow no planner until all finishes complete.
void finishBW(DT& mesh, BWPlan& plan, int thread = -1);
}
