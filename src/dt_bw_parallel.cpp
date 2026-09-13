#include "dt.h"
#include <exception>
#include <limits>

namespace dt {
double refineBWSizeLowerBound(DT& mesh, int tet) {
    double lower = 0;
    if (mesh.growsize >= 1 && mesh.growsize <= 1e100) {
        lower = std::numeric_limits<double>::max();
        for (int node : mesh.Elems[tet].form) {
            const auto& n = mesh.Nodes[node];
            // Bound distances and products, including the inf*0 case when grow=1.
            // Outside this safe finite range, use only the existing minEdge bound.
            if (!(n.space >= 1e-100 && n.space <= 1e100) ||
                !(std::abs(n.pt[0]) <= 1e100 && std::abs(n.pt[1]) <= 1e100 && std::abs(n.pt[2]) <= 1e100)) {
                lower = 0; break;
            }
            lower = std::min(lower, n.space);
        }
        // Every interpolation term is >= its source size. Allow more than the
        // four-term summation's rounding error before applying the same clamps.
        lower *= 1 - 8 * std::numeric_limits<double>::epsilon();
        if (mesh.maxEdge > 0) lower = std::min(mesh.maxEdge, lower);
    }
    if (mesh.minEdge > 0) lower = std::max(mesh.minEdge, lower);
    return lower;
}
namespace {

struct SpatialEntry { uint32_t key; int id; };

int regionCount(int workers, size_t points) {
    int regions = std::min(workers * BW_MIN_REGIONS_PER_THREAD, BW_MAX_REGIONS);
    const int limit = std::min(workers * BW_MAX_REGIONS_PER_THREAD, BW_MAX_REGIONS);
    while (regions < limit && points / static_cast<size_t>(regions) > BW_TARGET_POINTS_PER_REGION)
        regions = std::min(regions * 2, limit);
    return static_cast<int>(std::min(static_cast<size_t>(regions),
        std::max(size_t(1), points / BW_MIN_POINTS_PER_REGION)));
}

uint32_t spreadBits(uint32_t value) {
    value = (value | (value << 16)) & 0x030000ffu;
    value = (value | (value << 8)) & 0x0300f00fu;
    value = (value | (value << 4)) & 0x030c30c3u;
    return (value | (value << 2)) & 0x09249249u;
}

struct SpatialFrame {
    std::array<double, 3> lower, inverse;
    explicit SpatialFrame(const DT& mesh) {
        for (int axis = 0; axis < 3; ++axis) {
            lower[axis] = mesh.minW[axis];
            const double width = mesh.maxW[axis] - mesh.minW[axis];
            inverse[axis] = width > 0 ? 1023.0 / width : 0;
        }
    }
    uint32_t key(const double* point) const {
        uint32_t result = 0;
        for (int axis = 0; axis < 3; ++axis) {
            const double scaled = (point[axis] - lower[axis]) * inverse[axis];
            const uint32_t coordinate = static_cast<uint32_t>(std::max(0.0, std::min(1023.0, scaled)));
            result |= spreadBits(coordinate) << axis;
        }
        return result;
    }
};

// Three linear passes, stable within a voxel. The existing multiscale Hilbert
// order remains useful locally, but its separate levels are not global regions.
void spatialSort(std::vector<SpatialEntry>& entries, std::vector<SpatialEntry>& scratch, int workers = 1) {
    DTParallelScope parallelScope;
    scratch.resize(entries.size());
    if (workers > 1) {
        std::vector<std::array<size_t, 1024>> offsets(workers);
#pragma omp parallel num_threads(workers)
        {
            const int thread = omp_get_thread_num(), team = omp_get_num_threads();
            const size_t begin = entries.size() * thread / team;
            const size_t end = entries.size() * (thread + 1) / team;
            auto& offset = offsets[thread];
            for (int shift = 0; shift < 30; shift += 10) {
                offset.fill(0);
                for (size_t i = begin; i < end; ++i) ++offset[(entries[i].key >> shift) & 1023];
#pragma omp barrier
#pragma omp single
                {
                    size_t total = 0;
                    for (int bin = 0; bin < 1024; ++bin)
                        for (int t = 0; t < team; ++t) {
                            const size_t count = offsets[t][bin];
                            offsets[t][bin] = total; total += count;
                        }
                }
                // Contiguous input chunks and ordered offsets preserve serial tie order.
                for (size_t i = begin; i < end; ++i)
                    scratch[offset[(entries[i].key >> shift) & 1023]++] = entries[i];
#pragma omp barrier
#pragma omp single
                entries.swap(scratch);
            }
        }
        return;
    }
    for (int shift = 0; shift < 30; shift += 10) {
        std::array<size_t, 1024> offset{};
        for (const auto& entry : entries) ++offset[(entry.key >> shift) & 1023];
        size_t total = 0;
        for (auto& count : offset) { const size_t old = count; count = total; total += old; }
        for (const auto& entry : entries) scratch[offset[(entry.key >> shift) & 1023]++] = entry;
        entries.swap(scratch);
    }
}

// A batch reads a stable mesh. Membership stamps are touched only by its
// coordinator, and all storage growth happens outside worker regions.
struct BWBatchWorkspace {
    std::vector<uint32_t> readers, writers, nodeWriters;
    uint32_t epoch = 0;
    std::vector<int> selected, order;
    std::vector<std::pair<double, int>> spacingOrder;
    bool spacingIndexValid = true;
    int spacingAxis = 0;
    std::vector<std::vector<int>> slots;
    std::vector<int> nodes;
    std::vector<std::exception_ptr> errors;

    void begin(DT& mesh, size_t jobs, bool nodeAccess) {
        readers.resize(mesh.Elems.size(), 0);
        writers.resize(mesh.Elems.size(), 0);
        if (nodeAccess) nodeWriters.resize(mesh.Nodes.size(), 0);
        if (++epoch == 0) {
            std::fill(readers.begin(), readers.end(), 0);
            std::fill(writers.begin(), writers.end(), 0);
            std::fill(nodeWriters.begin(), nodeWriters.end(), 0);
            epoch = 1;
        }
        selected.clear(); order.clear(); spacingOrder.clear(); spacingIndexValid = true;
        spacingAxis = 0;
        for (int axis = 1; axis < 3; ++axis)
            if (mesh.maxW[axis] - mesh.minW[axis] > mesh.maxW[spacingAxis] - mesh.minW[spacingAxis]) spacingAxis = axis;
        slots.resize(jobs);
        nodes.resize(jobs, -1);
        errors.assign(jobs, nullptr);
    }

    bool overlaps(const BWPlan& p) const {
        for (int t : p.writeTets)
            if (writers[t] == epoch || readers[t] == epoch) return true;
        for (int t : p.readTets)
            if (writers[t] == epoch) return true;
        for (int n : p.writeNodes)
            if (nodeWriters[n] == epoch) return true;
        return false;
    }

    void include(const BWPlan& p) {
        for (int t : p.readTets) readers[t] = epoch;
        for (int t : p.writeTets) writers[t] = epoch;
        for (int n : p.writeNodes) nodeWriters[n] = epoch;
    }
};

bool candidatesTooClose(DT& mesh, BWRequest& a, BWRequest& b) {
    if (a.point == b.point) return true;
    if (a.info != 2 || b.info != 2 || a.boundary || b.boundary) return false;
    const bool anisotropic = a.anisotropic && b.anisotropic;
    if (!anisotropic) {
        // The spacing predicate can only succeed inside this axis-aligned box.
        // Most candidates belong to distant regions, so avoid their square root.
        const double radius = std::max(mesh.minEdge,
            std::min(mesh.meshSize * a.space, mesh.meshSize * b.space));
        // Preserve the original floating-point predicate when squaring this
        // scale can underflow (or overflow), instead of using a stronger bound.
        if (std::isnormal(radius * radius))
            for (int axis = 0; axis < 3; ++axis)
                if (std::abs(a.point[axis] - b.point[axis]) >= radius) return false;
    }
    double squared = 0;
    for (int j = 0; j < 3; ++j) {
        const double v = a.point[j] - b.point[j];
        squared += v * v;
    }
    if (squared == 0) return true;
    const double length = std::sqrt(squared);
    if (anisotropic) {
        if (mesh.minEdge > 0 && length < mesh.minEdge * mesh.ani_lower) return true;
        return mesh.cal_ani_length(a.point.data(), b.point.data(),
                                  a.metric.data(), b.metric.data()) < mesh.ani_lower;
    }
    return (mesh.minEdge > 0 && length < mesh.minEdge) ||
        (length < mesh.meshSize * a.space && length < mesh.meshSize * b.space);
}

uint64_t planWork(const BWPlan& p) {
    return static_cast<uint64_t>(p.cavity.size()) + p.faces.size() + p.readTets.size() + p.writeTets.size();
}

// Sorted-coordinate broad phase; the existing exact spacing predicate remains
// authoritative. Unusual scales and anisotropy use the exhaustive fallback.
bool spacingConflict(DT& mesh, std::vector<BWPlan>& plans, int index,
                     BWBatchWorkspace& work, bool refinement) {
    BWRequest& a = plans[index].request;
    auto close = [&](int other) { return candidatesTooClose(mesh, a, plans[other].request); };
    const double radius = std::max(mesh.minEdge, mesh.meshSize * a.space);
    if (refinement && work.spacingIndexValid && !a.anisotropic && a.info == 2 && !a.boundary &&
        std::isfinite(a.point[work.spacingAxis]) && radius > 0 && std::isnormal(radius * radius)) {
        const double lower = std::nextafter(a.point[work.spacingAxis] - radius, -std::numeric_limits<double>::infinity());
        const double upper = std::nextafter(a.point[work.spacingAxis] + radius, std::numeric_limits<double>::infinity());
        auto first = std::lower_bound(work.spacingOrder.begin(), work.spacingOrder.end(), lower,
            [](const std::pair<double,int>& value, double bound) { return value.first < bound; });
        auto last = std::upper_bound(first, work.spacingOrder.end(), upper,
            [](double bound, const std::pair<double,int>& value) { return bound < value.first; });
        for (; first != last; ++first) if (close(first->second)) return true;
    } else for (int other : work.selected) if (close(other)) return true;
    return false;
}

struct BWBatchResult {
    uint64_t inserted = 0, deferred = 0, work = 0, wastedWork = 0;
};

BWBatchResult commitBatch(DT& mesh, std::vector<BWPlan>& plans, BWBatchWorkspace& work,
                 Args* refineArgs, int workers) {
    DTParallelScope parallelScope;
    BWBatchResult result;
    work.begin(mesh, plans.size(), refineArgs == nullptr);
    for (int i = 0; i < static_cast<int>(plans.size()); ++i)
        if (plans[i].status == BWStatus::Ready) work.order.push_back(i);
    if (refineArgs) {
        // Expensive cavities win before cheap overlapping candidates. Ties use
        // stable lane IDs, never worker completion order.
        std::sort(work.order.begin(), work.order.end(), [&](int a, int b) {
            const uint64_t ca = planWork(plans[a]), cb = planWork(plans[b]);
            return ca != cb ? ca > cb : a < b;
        });
    }
    for (int i : work.order) {
        BWPlan& plan = plans[i];
        if (refineArgs) result.work += planWork(plan);
        bool conflict = work.overlaps(plan);
        if (!conflict && spacingConflict(mesh, plans, i, work, refineArgs != nullptr)) {
            conflict = true;
        }
        if (conflict) {
            plan.status = BWStatus::RetryConflict;
            ++result.deferred;
            if (refineArgs) result.wastedWork += planWork(plan);
            continue;
        }
        work.include(plan);
        work.selected.push_back(i);
        if (refineArgs && work.spacingIndexValid) {
            const double key = plan.request.point[work.spacingAxis];
            if (plan.request.anisotropic || !std::isfinite(key)) work.spacingIndexValid = false;
            else {
                const std::pair<double,int> entry(key, i);
                work.spacingOrder.insert(std::lower_bound(work.spacingOrder.begin(), work.spacingOrder.end(), entry), entry);
            }
        }
    }
    if (work.selected.empty()) return result;

    // Reserve actual IDs, not just vector capacity: SmallVector::operator[]
    // reads vec_size, so even capacity-preserving push_back is forbidden below.
    for (int index : work.selected) {
        BWPlan& plan = plans[index];
        int node = plan.request.node;
        if (node < 0) {
            node = mesh.addNode(plan.request.point[0], plan.request.point[1],
                                plan.request.point[2], plan.request.space);
            if (plan.request.anisotropic) {
                mesh.AniSol.resize(mesh.Nodes.size());
                mesh.AniSol[node] = plan.request.metric;
            }
        }
        work.nodes[index] = node;
        auto& slots = work.slots[index];
        slots.clear();
        slots.reserve(plan.faces.size());
        for (size_t f = 0; f < plan.faces.size(); ++f) slots.push_back(mesh.addElem());
    }

    const int count = static_cast<int>(work.selected.size());
    // prepareBWFill has performed all validation/allocation. Only independent
    // mesh fields are written here; no traversal, pool growth or callback occurs.
#pragma omp parallel for num_threads(workers) schedule(static) if(count > 1 && workers > 1)
    for (int j = 0; j < count; ++j) {
        const int index = work.selected[j];
        try {
            commitBW(mesh, plans[index], work.nodes[index], work.slots[index]);
        } catch (...) { work.errors[index] = std::current_exception(); }
    }
    for (const auto& error : work.errors) if (error) std::rethrow_exception(error);
    // Deletion/recycling, the ghost representative, and callbacks are finalized
    // after every worker has finished and before the next planning phase.
    for (int index : work.selected) {
        finishBW(mesh, plans[index]);
        if (refineArgs) mesh.updateSize(work.nodes[index], *refineArgs);
        ++result.inserted;
    }
    return result;
}

int liveHint(DT& mesh, int anchor) {
    if (anchor >= 0 && anchor < static_cast<int>(mesh.Nodes.size())) {
        const int tet = mesh.getP2T(anchor);
        if (tet >= 0 && tet < static_cast<int>(mesh.Elems.size()) &&
            !mesh.isDelEle(tet) && mesh.isNod_in_Tet(anchor, tet) >= 0) return tet;
    }
    for (int t = 0; t < static_cast<int>(mesh.Elems.size()); ++t)
        if (!mesh.isDelEle(t) && !mesh.ishulltet(t)) return t;
    throw std::runtime_error("No live tetrahedron for BW point location");
}

void recordInitialResult(DT& mesh, BWPlan& plan, int& anchor) {
    if (plan.status == BWStatus::Inserted) {
        anchor = plan.request.node;
    } else if (plan.status == BWStatus::Duplicate) {
        mesh.Nodes[plan.request.node].info = plan.duplicateNode + 1;
        anchor = plan.duplicateNode;
        if (mesh.infolevel > 0)
            mesh.meshLogger->warn("Duplicate point {} and {}", plan.request.node, plan.duplicateNode);
    } else if (plan.status == BWStatus::Stale) {
        // The mesh is frozen during planning and liveHint has been checked.
        // Stale here therefore indicates broken adjacency, not contention.
        throw std::runtime_error("Invalid mesh adjacency while planning BW input point " + std::to_string(plan.request.node));
    } else if (plan.status != BWStatus::RetryConflict) {
        throw std::runtime_error("BW could not insert input point " + std::to_string(plan.request.node));
    }
}

// Candidate construction never changes Nodes, element flags or allocation
// queues. Rejection flags are applied only after all readers have finished.
// Both callers have already checked eligibility on the frozen mesh.
bool makeRefineRequest(DT& mesh, int tet, const Args& args, BWRequest& request) {
    mesh.calBarycenter(tet, request.point.data());
    request.space = 0;
    const int* vertices = mesh.Elems[tet].form;
    const double lower = args.size * refineBWSizeLowerBound(mesh, tet);
    double length[4];
    for (int j = 0; j < 4; ++j) {
        const int node = vertices[j];
        length[j] = mesh.distance(request.point.data(), mesh.Nodes[node].pt);
        if (length[j] < lower && length[j] < args.size * mesh.Nodes[node].space) return false;
        request.space += mesh.Nodes[node].space + length[j] * (mesh.growsize - 1);
    }
    request.space /= 4.;
    if (mesh.maxEdge > 0) request.space = std::min(mesh.maxEdge, request.space);
    if (mesh.minEdge > 0) request.space = std::max(mesh.minEdge, request.space);
    for (int j = 0; j < 4; ++j)
        if (length[j] < args.size * request.space &&
            length[j] < args.size * mesh.Nodes[vertices[j]].space) return false;

    // Most inspected cells fail the cheap size test. Initialize planning-only
    // state after it passes; calBarycenter already overwrites all coordinates.
    request.node = -1;
    request.boundary = false;
    request.trackAccess = true;
    request.deferP2T = false;
    request.info = 2;
    request.seeds.resize(1);
    request.seeds[0] = tet;
    request.anisotropic = !mesh.AniSol.empty();
    if (request.anisotropic) {
        request.metric.fill(0);
        for (int n : mesh.Elems[tet].form) for (int j = 0; j < 6; ++j)
            request.metric[j] += mesh.AniSol[n][j] * .25;
    }
    return true;
}

struct RefineLane {
    int next = 0, end = 0, candidate = -1;
    std::array<int, 4> candidateVertices{};
    std::vector<int> rejected;
};

bool currentCandidate(DT& mesh, const RefineLane& lane) {
    const int t = lane.candidate;
    if (t < 0 || t >= static_cast<int>(mesh.Elems.size()) || mesh.isDelEle(t)) return false;
    for (int j = 0; j < 4; ++j)
        if (mesh.Elems[t].form[j] != lane.candidateVertices[j]) return false;
    return true;
}
}

void insertDelaunayPoints(DT& mesh, const std::vector<int>& order) {
    MeshStageLog stageLog(mesh, "Delaunay insertion", 2);
    DTParallelScope parallelScope;
    if (order.size() <= 4) return;
    if (mesh.infolevel > 0) mesh.meshLogger->debug("Incremental insert points.");
    const int workers = activeDTThreads(mesh, order.size());
    const int remainingCount = static_cast<int>(order.size()) - 4;
    BWPlan serialPlan;
    std::vector<int> serialSlots;
    int anchor = order[0];
    auto insertSerial = [&](int node, int& hint) {
        BWRequest request = makeBWRequest(mesh, node, {liveHint(mesh, hint)}, 0);
        request.trackAccess = false;
        planBW(mesh, request, serialPlan);
        if (serialPlan.status == BWStatus::Ready) {
            serialSlots.clear();
            for (size_t f = 0; f < serialPlan.faces.size(); ++f) serialSlots.push_back(mesh.addElem());
            commitBW(mesh, serialPlan, node, serialSlots);
            finishBW(mesh, serialPlan);
        }
        recordInitialResult(mesh, serialPlan, hint);
    };
    if (workers == 1) {
        for (int i = 4; i < static_cast<int>(order.size()); ++i) insertSerial(order[i], anchor);
        return;
    }

    const int laneCount = regionCount(workers, remainingCount);
    const SpatialFrame frame(mesh);
    std::vector<SpatialEntry> spatial, scratch;
    spatial.reserve(remainingCount);
    for (size_t i = 4; i < order.size(); ++i)
        spatial.push_back({frame.key(mesh.Nodes[order[i]].pt), order[i]});
    spatialSort(spatial, scratch);
    // Partition by global position, then keep the original multiscale ordering
    // inside each region. A pure spatial sweep grows a long insertion front;
    // the early levels supply a progressively denser local scaffold instead.
    std::vector<int> regionOfNode(mesh.Nodes.size(), -1), offset(laneCount);
    for (int lane = 0; lane < laneCount; ++lane) {
        const int begin = static_cast<int>(int64_t(remainingCount) * lane / laneCount);
        const int end = static_cast<int>(int64_t(remainingCount) * (lane + 1) / laneCount);
        offset[lane] = begin;
        for (int i = begin; i < end; ++i) regionOfNode[spatial[i].id] = lane;
    }
    for (size_t i = 4; i < order.size(); ++i) {
        const int node = order[i];
        spatial[offset[regionOfNode[node]]++] = {0, node};
    }
    // Select real input points near each region's six extremal faces, plus a
    // small evenly spread sample. These scaffold points reduce large cavities;
    // they are not constraints and never replace exact conflict checks.
    std::vector<unsigned char> warmed(remainingCount, 0);
    for (int lane = 0; lane < laneCount; ++lane) {
        const int begin = static_cast<int>(int64_t(remainingCount) * lane / laneCount);
        const int end = static_cast<int>(int64_t(remainingCount) * (lane + 1) / laneCount);
        std::array<int, 6> extreme; extreme.fill(begin);
        for (int i = begin + 1; i < end; ++i) {
            const double* point = mesh.Nodes[spatial[i].id].pt;
            for (int axis = 0; axis < 3; ++axis) {
                if (point[axis] < mesh.Nodes[spatial[extreme[2 * axis]].id].pt[axis]) extreme[2 * axis] = i;
                if (point[axis] > mesh.Nodes[spatial[extreme[2 * axis + 1]].id].pt[axis]) extreme[2 * axis + 1] = i;
            }
        }
        for (int i : extreme) warmed[i] = 1;
        const int samples = std::min(BW_SCAFFOLD_SAMPLES, end - begin);
        for (int k = 0; k < samples; ++k) warmed[begin + (end - begin) * k / samples] = 1;
    }
    for (int i = 0; i < remainingCount; ++i) if (warmed[i]) {
        insertSerial(spatial[i].id, anchor);
    }
    std::vector<int> pending, next(laneCount), ends(laneCount), anchors(laneCount, anchor);
    pending.reserve(remainingCount);
    for (int lane = 0; lane < laneCount; ++lane) {
        const int begin = static_cast<int>(int64_t(remainingCount) * lane / laneCount);
        const int end = static_cast<int>(int64_t(remainingCount) * (lane + 1) / laneCount);
        next[lane] = static_cast<int>(pending.size());
        for (int i = begin; i < end; ++i) {
            if (!warmed[i]) pending.push_back(spatial[i].id);
            else anchors[lane] = mesh.Nodes[spatial[i].id].info > 0 ?
                mesh.Nodes[spatial[i].id].info - 1 : spatial[i].id;
        }
        ends[lane] = static_cast<int>(pending.size());
    }
    BWBatchWorkspace work;
    std::vector<BWPlan> plans(laneCount);
    std::vector<std::exception_ptr> errors(laneCount);
    for (;;) {
        int remaining = 0, activeLanes = 0;
        for (int lane = 0; lane < laneCount; ++lane) {
            plans[lane].status = BWStatus::Rejected;
            errors[lane] = nullptr;
            remaining += ends[lane] - next[lane];
            activeLanes += next[lane] < ends[lane];
        }
        if (!remaining) break;
        const int team = static_cast<size_t>(remaining) < mesh.parallelMinPointsPerThread / BW_SERIAL_TAIL_DIVISOR ? 1 :
            activeDTThreads(mesh, order.size(), activeLanes);
        // The small tail no longer benefits from conflict tracking, batches or
        // barriers. Keep each region's nearby anchor while draining it serially.
        if (team == 1) {
            for (int lane = 0; lane < laneCount; ++lane)
                while (next[lane] < ends[lane]) insertSerial(pending[next[lane]++], anchors[lane]);
            break;
        }
#pragma omp parallel for num_threads(team) schedule(dynamic, 1)
        for (int lane = 0; lane < laneCount; ++lane) {
            if (next[lane] >= ends[lane]) continue;
            try {
                BWRequest& request = plans[lane].request;
                request.node = pending[next[lane]];
                request.info = 0;
                request.trackAccess = true;
                request.seeds.resize(1);
                request.seeds[0] = liveHint(mesh, anchors[lane]);
                std::copy(mesh.Nodes[request.node].pt, mesh.Nodes[request.node].pt + 3, request.point.begin());
                request.space = mesh.Nodes[request.node].space;
                request.boundary = mesh.isbndpnt(request.node);
                request.anisotropic = !mesh.AniSol.empty();
                if (request.anisotropic) request.metric = mesh.AniSol[request.node];
                planBW(mesh, request, plans[lane]);
            } catch (...) { errors[lane] = std::current_exception(); }
        }
        for (const auto& error : errors) if (error) std::rethrow_exception(error);
        commitBatch(mesh, plans, work, nullptr, team);
        for (int lane = 0; lane < laneCount; ++lane) {
            if (next[lane] >= ends[lane]) continue;
            recordInitialResult(mesh, plans[lane], anchors[lane]);
            if (plans[lane].status == BWStatus::Inserted || plans[lane].status == BWStatus::Duplicate)
                ++next[lane];
        }
    }
}

void refineBWParallel(DT& mesh, Args& args) {
    MeshStageLog stageLog(mesh, "Parallel refinement", 2);
    DTParallelScope parallelScope;
    uint64_t inserted = 0, serialPlans = 0;
    std::vector<uint8_t> rejected(mesh.Elems.size(), 0);
    // Volume refinement never changes BndTri. Derive membership from its actual
    // keys, not Node::type, which need not mirror the boundary table.
    const std::vector<uint8_t> boundaryNodes = [&]() {
        std::vector<uint8_t> nodes(mesh.Nodes.size(), 0);
        const auto end = mesh.BndTri.end();
        for (auto face = mesh.BndTri.begin(); face != end; face = mesh.BndTri.next(face)) {
            for (int node : {face.n1, face.n2, face.n3}) {
                if (node < 0) continue;
                if (static_cast<size_t>(node) >= nodes.size()) nodes.resize(static_cast<size_t>(node) + 1, 0);
                nodes[node] = 1;
            }
        }
        return nodes;
    }();
    // Declared after the bitmap so every borrowed request is destroyed first.
    // No plan or worker survives this function, including exceptional exits.
    std::vector<RefineLane> lanes;
    std::vector<BWPlan> plans;
    std::vector<std::exception_ptr> errors;
    BWBatchWorkspace work;
    BWRefineTeamPolicy policy;
    const BWRefineTeamPolicy workloadPolicy;
    BWPlan serialPlan;
    std::vector<int> serialSlots;
    std::vector<int> pending(mesh.Elems.size()), nextPending;
    for (int t = 0; t < static_cast<int>(pending.size()); ++t) pending[t] = t;
    std::vector<uint32_t> queued;
    uint32_t queueEpoch = 1;
    bool serialFinish = false, serialAuditStarted = false;
    auto enqueueNext = [&](int t) {
        if (queued[t] != queueEpoch) {
            queued[t] = queueEpoch; nextPending.push_back(t);
        }
    };
    // Queue accepted replacements and untried conflicting source cells. IDs
    // may be recycled repeatedly; one queued ID denotes its latest geometry
    auto enqueueNew = [&](const BWPlan& plan) {
        for (int t : plan.newElements) { rejected[t] = 0; enqueueNext(t); }
    };
    auto eligible = [&](int tet) {
        return !mesh.isDelEle(tet) && !mesh.isvirtualtet(tet) && !mesh.ishulltet(tet) &&
            !rejected[tet];
    };
    auto insertSerial = [&](int tet) {
        if (!eligible(tet)) return;
        BWRequest& request = serialPlan.request;
        if (!makeRefineRequest(mesh, tet, args, request)) {
            rejected[tet] = 1;
            return;
        }
        request.boundaryNodes = &boundaryNodes;
        request.trackAccess = false; // No other planner/committer is running.
        planBW(mesh, request, serialPlan);
        ++serialPlans;
        if (serialPlan.status == BWStatus::Stale)
            throw std::runtime_error("Invalid adjacency during serial BW refinement");
        if (serialPlan.status != BWStatus::Ready) {
            rejected[tet] = 1;
            return;
        }
        const int node = mesh.addNode(request.point[0], request.point[1], request.point[2], request.space);
        if (request.anisotropic) {
            mesh.AniSol.resize(mesh.Nodes.size()); mesh.AniSol[node] = request.metric;
        }
        serialSlots.clear();
        for (size_t face = 0; face < serialPlan.faces.size(); ++face) serialSlots.push_back(mesh.addElem());
        commitBW(mesh, serialPlan, node, serialSlots);
        finishBW(mesh, serialPlan);
        mesh.updateSize(node, args); // Original calling thread, outside OpenMP.
        ++inserted;
        queued.resize(mesh.Elems.size(), 0);
        rejected.resize(mesh.Elems.size(), 0);
        enqueueNew(serialPlan);
    };
    const size_t workPerThread = mesh.parallelMinPointsPerThread <= 1 ? 1 :
        std::max(size_t(256), mesh.parallelMinPointsPerThread / 4);
    // Coarse sweep grain amortizes BW access tracking, partitioning and allocation.
    // The much smaller workPerThread gate still handles the within-sweep tail.
    const size_t sweepWorkPerThread = mesh.parallelMinPointsPerThread <= 1 ? 1 : size_t(262144);
    int loop = 0;
    std::vector<SpatialEntry> spatial, spatialScratch;
    const SpatialFrame frame(mesh);
    for (;;) {
        // Enter a permanent serial finish when the previous sweep made little
        // progress or the workload gate no longer offers useful parallelism.
        if (!serialFinish && (pending.empty() ||
            workloadPolicy.choose(activeDTThreads(mesh, 0, pending.size()), pending.size(), BW_MAX_REGIONS, sweepWorkPerThread) == 1))
            serialFinish = true;
        if (serialFinish && !serialAuditStarted) {
            serialAuditStarted = true;
            // One full serial recheck includes every deferred/recycled cell and
            // re-evaluates previous rejections. Subsequent sweeps use the queue.
            pending.clear(); pending.reserve(mesh.Elems.size());
            for (int t = 0; t < static_cast<int>(mesh.Elems.size()); ++t)
                if (!mesh.isDelEle(t) && !mesh.isvirtualtet(t) && !mesh.ishulltet(t)) {
                    rejected[t] = 0; pending.push_back(t);
                }
        }

        const uint64_t before = inserted;
        spatial.clear(); nextPending.clear();
        // Sparse sweeps skip Morton construction/sorting entirely.
        const int cap = serialFinish ? 1 : activeDTThreads(mesh, 0, pending.size());
        const int sortThreads = workloadPolicy.choose(cap, pending.size(), BW_MAX_REGIONS, sweepWorkPerThread);
        if (sortThreads > 1) {
            spatial.resize(pending.size());
#pragma omp parallel for num_threads(sortThreads) schedule(static)
            for (int i = 0; i < static_cast<int>(pending.size()); ++i) {
                const int tet = pending[i];
                SpatialEntry entry{0, -1};
                if (eligible(tet)) {
                    double center[3]; mesh.calBarycenter(tet, center);
                    entry = {frame.key(center), tet};
                }
                spatial[i] = entry;
            }
            size_t count = 0;
            for (size_t i = 0; i < spatial.size(); ++i)
                if (spatial[i].id >= 0) spatial[count++] = spatial[i];
            spatial.resize(count);
            spatialSort(spatial, spatialScratch, sortThreads);
        } else {
            for (int tet : pending) if (eligible(tet)) spatial.push_back({0, tet});
        }
        const int n = static_cast<int>(spatial.size());
        if (!n) {
            if (serialFinish) break;
            serialFinish = true; continue;
        }
        const int workers = workloadPolicy.choose(cap, spatial.size(), BW_MAX_REGIONS, sweepWorkPerThread);
        // Bound quadratic candidate-spacing selection when the caller requests many cores.
        const int laneCount = std::min(128, regionCount(workers, spatial.size()));
        lanes.resize(laneCount); plans.resize(laneCount); errors.resize(laneCount);
        for (int lane = 0; lane < laneCount; ++lane) {
            lanes[lane].next = static_cast<int>(int64_t(n) * lane / laneCount);
            lanes[lane].end = static_cast<int>(int64_t(n) * (lane + 1) / laneCount);
            lanes[lane].rejected.reserve(512);
            plans[lane].status = BWStatus::Rejected;
        }
        for (;;) {
            int activeLanes = 0;
            size_t remaining = 0;
            for (int lane = 0; lane < laneCount; ++lane) {
                activeLanes += lanes[lane].next < lanes[lane].end;
                remaining += lanes[lane].end - lanes[lane].next;
            }
            if (!activeLanes) break;
            const int team = policy.choose(workers, remaining, activeLanes, workPerThread);
            if (team == 1) {
                // Actual sequential insertion removes access-set construction,
                // artificial intra-batch conflicts, and batch allocation scans.
                size_t processed = 0;
                const uint64_t serialPlansBefore = serialPlans;
                for (int lane = 0; lane < laneCount && processed < 512; ++lane) {
                    while (lanes[lane].next < lanes[lane].end && processed < 512) {
                        const int tet = spatial[lanes[lane].next++].id;
                        insertSerial(tet); ++processed;
                    }
                    plans[lane].status = BWStatus::Rejected;
                }
                policy.serialProgress(static_cast<size_t>(serialPlans - serialPlansBefore));
                continue;
            }
            for (int lane = 0; lane < laneCount; ++lane) {
                plans[lane].status = BWStatus::Rejected;
                errors[lane] = nullptr;
                lanes[lane].candidate = -1; lanes[lane].rejected.clear();
            }
#pragma omp parallel for num_threads(team) schedule(dynamic, 1)
            for (int index = 0; index < laneCount; ++index) {
                RefineLane& lane = lanes[index];
                try {
                    for (int scanned = 0; lane.next < lane.end && scanned < 512; ++scanned) {
                        const int tet = spatial[lane.next].id;
                        if (!eligible(tet)) { ++lane.next; continue; }
                        BWRequest& request = plans[index].request;
                        if (!makeRefineRequest(mesh, tet, args, request)) {
                            lane.rejected.push_back(tet); ++lane.next; continue;
                        }
                        request.boundaryNodes = &boundaryNodes;
                        request.deferP2T = true;
                        lane.candidate = tet;
                        std::copy(mesh.Elems[tet].form, mesh.Elems[tet].form + 4, lane.candidateVertices.begin());
                        planBW(mesh, request, plans[index]);
                        break;
                    }
                } catch (...) { errors[index] = std::current_exception(); }
            }
            for (const auto& error : errors) if (error) std::rethrow_exception(error);
            uint64_t ready = 0;
            for (int index = 0; index < laneCount; ++index) {
                RefineLane& lane = lanes[index];
                for (int t : lane.rejected) rejected[t] = 1;
                if (lane.candidate < 0) continue;
                if (plans[index].status == BWStatus::Rejected || plans[index].status == BWStatus::Duplicate) {
                    rejected[lane.candidate] = 1; ++lane.next;
                } else if (plans[index].status == BWStatus::Stale) {
                    // Topology is frozen while planning; this is an error, not contention.
                    throw std::runtime_error("Invalid mesh adjacency while planning BW refinement at tetrahedron " +
                                             std::to_string(lane.candidate));
                } else if (plans[index].status == BWStatus::Ready) ++ready;
            }
            const BWBatchResult batch = commitBatch(mesh, plans, work, &args, team);
            inserted += batch.inserted;
            policy.observe(team, ready, batch.inserted, batch.deferred, batch.work, batch.wastedWork);
            // Growth and queue stamps are coordinator-only, after every commit.
            queued.resize(mesh.Elems.size(), 0);
            rejected.resize(mesh.Elems.size(), 0);
            for (int index = 0; index < laneCount; ++index) {
                if (plans[index].status == BWStatus::Inserted) {
                    ++lanes[index].next; enqueueNew(plans[index]);
                } else if (plans[index].status == BWStatus::RetryConflict) {
                    if (currentCandidate(mesh, lanes[index])) enqueueNext(lanes[index].candidate);
                    // Replacement geometry is queued by its successful plan.
                    ++lanes[index].next; // Defer this source until the next sweep.
                }
            }
        }
        const uint64_t success = inserted - before;
        if (mesh.infolevel > 0)
            mesh.meshLogger->debug("Loop:{: <2} add:{: <5} Elem:{: <8} Node:{: <8}",
                         loop++, success, mesh.Elems.size(), mesh.Nodes.size());
        if (!success && serialFinish) break;
        // Low progress makes the following sweep permanently serial.
        if (!serialFinish && success < std::max(uint64_t(1024), uint64_t(mesh.Nodes.size() / 100)))
            serialFinish = true;
        pending.swap(nextPending);
        if (++queueEpoch == 0) { std::fill(queued.begin(), queued.end(), 0); queueEpoch = 1; }
    }
}
}
