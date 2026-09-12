#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace dt {
class DT;

enum class BWStatus { Ready, Inserted, Duplicate, Rejected, RetryConflict, Stale };

// A request owns candidate data; node == -1 denotes an unallocated point.
struct BWRequest {
    int node = -1;
    std::array<double, 3> point{};
    double space = 0;
    std::array<double, 6> metric{};
    bool anisotropic = false;
    bool boundary = false;
    bool trackAccess = true;
    // Refinement only: commit writes the new node; finishBW serially updates old-node P2T.
    bool deferP2T = false;
    int info = 0;
    int initShell = 0;
    std::vector<int> seeds;
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
private:
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
BWStatus prepareBWFill(DT& mesh, BWPlan& plan);
// Internal adjacency builder, also exercised independently by regression tests.
bool pairBWBoundary(BWPlan& plan);
BWStatus planBW(DT& mesh, const BWRequest& request, BWPlan& plan);

// Caller freezes topology throughout planning, selects independent plans, and
// allocates every node/element slot before any worker enters commitBW.
// commitBW does not grow containers, recycle slots, or update global counters.
void commitBW(DT& mesh, BWPlan& plan, int node, const std::vector<int>& slots);
// Run after all commits have completed. This updates recycling, counters and P2T.
// With deferP2T, serialize finishes and allow no planner until all finishes complete.
void finishBW(DT& mesh, BWPlan& plan, int thread = -1);
}
