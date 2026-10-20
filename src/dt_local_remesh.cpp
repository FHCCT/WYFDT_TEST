#include "dt.h"
#include "spdlog/sinks/null_sink.h"

namespace dt {
namespace {
using Face = std::array<int, 3>;
Face faceKey(const int* v, int opposite) {
    Face key; int k = 0;
    for (int i = 0; i < 4; ++i) if (i != opposite) key[k++] = v[i];
    std::sort(key.begin(), key.end());
    return key;
}
struct PatchBoundary { int cell, face, outside, outsideFace; };
struct PatchFace { int cell = -1, face = -1, count = 0; };

// All allocation and validation precede the first edit to the host mesh.
bool replacePatch(DT& host, const std::vector<int>& cells, Mesh& mesh,
    const std::vector<int>& toGlobal, const std::map<Face, PatchBoundary>& boundary,
    const OptimizationQuality& before, const DT::TopologyCandidate& seed,
    OptimizationQuality* globalQuality) {
    if (mesh.V.size() != toGlobal.size() || mesh.T.empty() || mesh.T.size() > cells.size() * 2)
        return false;
    // outMesh may compact/reorder points. Match coordinates, never assume that
    // a returned point number is the original local number.
    std::map<std::array<double, 3>, int> coordinates;
    for (int n : toGlobal)
        if (!coordinates.emplace(std::array<double,3>{{host.Nodes[n].pt[0],
            host.Nodes[n].pt[1], host.Nodes[n].pt[2]}}, n).second) return false;
    std::vector<int> returned;
    std::set<int> matched, used;
    for (const auto& point : mesh.V) {
        const auto it = coordinates.find(point);
        if (it == coordinates.end() || it->first != point || !matched.insert(it->second).second) return false;
        returned.push_back(it->second);
    }
    std::vector<Elem> replacement;
    replacement.reserve(mesh.T.size());
    std::map<Face, PatchFace> faces;
    std::set<std::array<int,4>> uniqueCells;
    auto seedKey = seed.form; std::sort(seedKey.begin(), seedKey.end());
    OptimizationQuality after;
    long double oldVolume = 0, newVolume = 0;
    for (int t : cells) oldVolume += host.calVolume(t);
    for (const auto& cell : mesh.T) {
        std::array<int,4> v;
        for (int i = 0; i < 4; ++i) {
            if (cell[i] < 0 || cell[i] >= static_cast<int>(returned.size())) return false;
            v[i] = returned[cell[i]]; used.insert(v[i]);
        }
        auto key = v; std::sort(key.begin(), key.end());
        if (std::adjacent_find(key.begin(),key.end()) != key.end() || key == seedKey ||
            !uniqueCells.insert(key).second) return false;
        double *a = host.Nodes[v[0]].pt, *b = host.Nodes[v[1]].pt,
            *c = host.Nodes[v[2]].pt, *d = host.Nodes[v[3]].pt;
        const double volume = GEOM_FUNC::orient3d(a,b,d,c);
        if (!std::isfinite(volume) || !(volume > std::max(0.0,host.minVolume_bw))) return false;
        newVolume += volume;
        addOptimizationQuality(host,after,a,b,c,d);
        Elem element(v[0],v[1],v[2],v[3]); element.geo = host.Elems[seed.tet].geo;
        element.q = host.tetquality(a,b,c,d,nullptr,host.improve_Metric);
        if (!(element.q > 0) || !std::isfinite(element.q)) return false;
        const int local = static_cast<int>(replacement.size());
        for (int f = 0; f < 4; ++f) {
            auto& entry = faces[faceKey(v.data(),f)];
            if (++entry.count > 2) return false;
            if (entry.count == 1) { entry.cell = local; entry.face = f; }
            else {
                // Opposite vertices must lie on opposite sides of this face.
                const auto& other = replacement[entry.cell];
                const Face triangle = faceKey(v.data(),f);
                const double x = GEOM_FUNC::orient3d(host.Nodes[triangle[0]].pt,
                    host.Nodes[triangle[1]].pt,host.Nodes[triangle[2]].pt,host.Nodes[v[f]].pt);
                const double y = GEOM_FUNC::orient3d(host.Nodes[triangle[0]].pt,
                    host.Nodes[triangle[1]].pt,host.Nodes[triangle[2]].pt,host.Nodes[other.form[entry.face]].pt);
                if (!((x < 0 && y > 0) || (x > 0 && y < 0))) return false;
                element.neig[f] = (int64_t(entry.cell) << 2) | entry.face;
                replacement[entry.cell].neig[entry.face] = (int64_t(local) << 2) | f;
            }
        }
        replacement.push_back(element);
    }
    if (used != matched || !after.valid || !after.count ||
        std::abs(newVolume-oldVolume) > 1e-10L * std::max(oldVolume,1e-30L)) return false;
    if (host.infolevel > 1) host.meshLogger->debug(
        "Local remesh quality: minA={:.6g}->{:.6g} minSJ={:.6g}->{:.6g} meanA={:.6g}->{:.6g} meanSJ={:.6g}->{:.6g}",
        before.minAngle,after.minAngle,before.minJacobian,after.minJacobian,
        before.meanAngle(),after.meanAngle(),before.meanJacobian(),after.meanJacobian());
    // Preserve both bad tails and local means. SUS is only the generator's
    // objective; accept a block by its actual angle and scaled Jacobian.
    if (after.minAngle + 1e-10 < before.minAngle ||
        after.minJacobian + 1e-12 < before.minJacobian ||
        after.tinyAngle > before.tinyAngle || after.tinyJacobian > before.tinyJacobian ||
        after.meanAngle() + 1e-10 < before.meanAngle() ||
        after.meanJacobian() + 1e-12 < before.meanJacobian() ||
        !(after.minAngle > before.minAngle * 1.01 + 1e-10 ||
          after.minJacobian > before.minJacobian * 1.01 + 1e-12)) return false;
    if (globalQuality) {
        const double delta = double(after.count)-double(before.count);
        if (!globalQuality->valid ||
            after.angleSum-before.angleSum + 1e-10L < globalQuality->meanAngle()*delta ||
            after.jacobianSum-before.jacobianSum + 1e-12L < globalQuality->meanJacobian()*delta)
            return false;
    }
    size_t exterior = 0;
    for (const auto& item : faces) if (item.second.count == 1) {
        ++exterior;
        const auto it = boundary.find(item.first);
        if (it == boundary.end()) return false;
        const auto& b = it->second;
        if (b.outside < 0 || b.outside >= static_cast<int>(host.Elems.size()) ||
            host.isDelEle(b.outside) || host.getNeig(b.outside,b.outsideFace) != b.cell ||
            host.getNeigOrd(b.outside,b.outsideFace) != b.face) return false;
        const Face& triangle = item.first;
        const double oldSide = GEOM_FUNC::orient3d(host.Nodes[triangle[0]].pt,
            host.Nodes[triangle[1]].pt,host.Nodes[triangle[2]].pt,
            host.Nodes[host.Elems[b.cell].form[b.face]].pt);
        const double newSide = GEOM_FUNC::orient3d(host.Nodes[triangle[0]].pt,
            host.Nodes[triangle[1]].pt,host.Nodes[triangle[2]].pt,
            host.Nodes[replacement[item.second.cell].form[item.second.face]].pt);
        if (!((oldSide < 0 && newSide < 0) || (oldSide > 0 && newSide > 0))) return false;
    }
    if (exterior != boundary.size() || !host.topologyCandidateCurrent(seed)) return false;

    const size_t count = replacement.size(), originalSize = host.Elems.size();
    const size_t extra = count > cells.size() ? count-cells.size() : 0;
    std::vector<int> slots(count);
    for (size_t i = 0; i < count; ++i)
        slots[i] = i < cells.size() ? cells[i] : static_cast<int>(originalSize+i-cells.size());
    // Prepare the free list too: queue growth must not throw mid-commit.
    auto freeCells = host.Evacancy;
    for (size_t i = count; i < cells.size(); ++i) freeCells.push(cells[i]);
    const int region = host.Elems[seed.tet].geo;
    auto regionCount = host.GeoNum.find(region);
    if (host.QuantityControl && regionCount == host.GeoNum.end()) return false;
    host.Elems.reserve(originalSize+extra);
    for (size_t i = 0; i < extra; ++i) host.addElem(false);
    for (auto& element : replacement) for (int f = 0; f < 4; ++f)
        if (element.neig[f] >= 0)
            element.neig[f] = (int64_t(slots[element.neig[f] >> 2]) << 2) | (element.neig[f] & 3);
    for (const auto& item : faces) if (item.second.count == 1) {
        const auto& b = boundary.at(item.first);
        auto& element = replacement[item.second.cell];
        element.neig[item.second.face] = (int64_t(b.outside) << 2) | b.outsideFace;
    }
    for (size_t i = 0; i < count; ++i) {
        host.Elems[slots[i]] = replacement[i];
        for (int n : replacement[i].form) host.Nodes[n].tet = slots[i];
    }
    for (const auto& item : faces) if (item.second.count == 1) {
        const auto& b = boundary.at(item.first);
        host.Elems[b.outside].neig[b.outsideFace] =
            (int64_t(slots[item.second.cell]) << 2) | item.second.face;
    }
    for (size_t i = count; i < cells.size(); ++i) {
        auto& e = host.Elems[cells[i]];
        std::fill(e.form,e.form+4,-1); std::fill(e.neig,e.neig+4,int64_t(-1));
        e.info = -1; e.geo = -1; e.q = -1;
    }
    host.Evacancy.swap(freeCells);
    if (host.QuantityControl) regionCount->second += static_cast<int>(count)-static_cast<int>(cells.size());
    if (globalQuality) {
        globalQuality->count = globalQuality->count-before.count+after.count;
        globalQuality->angleSum += after.angleSum-before.angleSum;
        globalQuality->jacobianSum += after.jacobianSum-before.jacobianSum;
        globalQuality->tinyAngle = globalQuality->tinyAngle-before.tinyAngle+after.tinyAngle;
        globalQuality->tinyJacobian = globalQuality->tinyJacobian-before.tinyJacobian+after.tinyJacobian;
        // Local minima cannot decrease, so these remain conservative floors.
    }
    for (int n : toGlobal) if (n < static_cast<int>(host.susIdleStates.size()))
        host.susIdleStates[n] = DT::SusIdleState{};
    if (host.infolevel > 1) host.meshLogger->debug(
        "Local remesh: cells={}->{} minAngle={:.6g}->{:.6g} minSJ={:.6g}->{:.6g}",
        cells.size(),count,before.minAngle,after.minAngle,before.minJacobian,after.minJacobian);
    return true;
}
}

int DT::remeshTopologyPatch(int tet, OptimizationQuality* globalQuality) {
    if (parallelTopologyBatch || omp_in_parallel())
        throw std::logic_error("Local remeshing must run outside a parallel region");
    if (!improve_step || !AniSol.empty() || topologyRemeshTrialBudget == 0 || tet < 0 ||
        tet >= static_cast<int>(Elems.size()) || isDelEle(tet) || ishulltet(tet) ||
        isvirtualtet(tet)) return 0;
    --topologyRemeshTrialBudget;
    // Standalone failure repairs also need count-corrected global acceptance.
    // The final batch supplies a shared, incrementally updated measurement.
    OptimizationQuality measuredGlobal;
    if (!globalQuality) {
        measuredGlobal = measureOptimizationQuality(*this);
        if (!measuredGlobal.valid) return 0;
        globalQuality = &measuredGlobal;
    }
    const auto seed = topologyCandidate(tet);
    const int region = Elems[tet].geo;
    constexpr size_t maxCells = 96, maxNodes = 64;
    std::vector<int> cells(1,tet), frontier(1,tet);
    std::set<int> included{tet};
    std::unique_ptr<DT> trial;
    for (int ring = 1; ring <= 4; ++ring) {
        std::vector<int> next;
        for (int t : frontier) for (int f = 0; f < 4; ++f) {
            const Face face = faceKey(Elems[t].form,f);
            if (isBndTri(face[0],face[1],face[2])) continue;
            const int n = getNeig(t,f);
            if (n < 0 || n >= static_cast<int>(Elems.size()) || included.count(n) ||
                isDelEle(n) || ishulltet(n) || isvirtualtet(n) || Elems[n].geo != region) continue;
            if (cells.size() >= maxCells) continue;
            included.insert(n); cells.push_back(n); next.push_back(n);
        }
        frontier.swap(next);
        if (ring == 1) continue;
        Mesh mesh;
        std::map<int,int> local;
        std::vector<int> toGlobal;
        for (int t : cells) for (int n : Elems[t].form)
            if (local.emplace(n,static_cast<int>(toGlobal.size())).second) {
                toGlobal.push_back(n);
                mesh.V.push_back({{Nodes[n].pt[0],Nodes[n].pt[1],Nodes[n].pt[2]}});
            }
        if (toGlobal.size() > maxNodes) break;
        std::map<Face,PatchBoundary> boundary;
        for (int t : cells) for (int f = 0; f < 4; ++f) {
            const Face face = faceKey(Elems[t].form,f);
            if (included.count(getNeig(t,f))) {
                // An internal constrained sheet can be reached around an edge.
                // This first implementation only regenerates one unsplit block.
                if (isBndTri(face[0],face[1],face[2])) return 0;
                continue;
            }
            if (!boundary.emplace(face,PatchBoundary{t,f,getNeig(t,f),getNeigOrd(t,f)}).second) return 0;
            int ia,ib,ic,id; DNC(f,ia,ib,ic,id);
            mesh.F.push_back({{local[Elems[t].form[ib]],local[Elems[t].form[ic]],
                local[Elems[t].form[id]],static_cast<int>(mesh.F.size())}});
        }
        const auto before = measureOptimizationQuality(*this,&cells);
        if (!before.valid || boundary.empty()) return 0;
        std::map<std::array<int,2>,int> edgeUses;
        for (const auto& entry : boundary) for (int k = 0; k < 3; ++k) {
            std::array<int,2> edge{{entry.first[k],entry.first[(k+1)%3]}};
            std::sort(edge.begin(),edge.end()); ++edgeUses[edge];
        }
        bool closed = true;
        for (const auto& edge : edgeUses) if (edge.second != 2) { closed = false; break; }
        if (!closed) {
            if (frontier.empty()) break;
            continue;
        }
        if (!trial) {
            trial.reset(new DT());
            trial->meshLogger = std::make_shared<spdlog::logger>("local_remesh_trial",
                std::make_shared<spdlog::sinks::null_sink_mt>());
        }
        Args args; args.infolevel = 0; args.refine = 0; args.optlevel = 0;
        args.nthread = 1; args.constrain = 1; args.outlogfile = 0;
        bool generated = false;
        try { generated = trial->tetrahedralize(mesh,args) != 0; }
        catch (const std::bad_alloc&) { throw; }
        catch (...) { trial.reset(); }
        if (infolevel > 1) meshLogger->debug(
            "Local remesh trial: ring={} cells={}->{} points={}->{} generated={}",
            ring,cells.size(),mesh.T.size(),toGlobal.size(),mesh.V.size(),generated);
        if (generated && replacePatch(*this,cells,mesh,toGlobal,boundary,before,seed,globalQuality)) return 1;
        if (frontier.empty() || cells.size() >= maxCells) break;
    }
    return 0;
}

int DT::repairBadTopologyPass(double angleDegrees) {
    if (parallelTopologyBatch || omp_in_parallel())
        throw std::logic_error("Local remeshing must run outside a parallel region");
    if (!topologyRepairEnabled || topologyRepairDeferred) return 0;
    MeshStageLog stageLog(*this,"Topology repair",2);
    std::vector<TopologyCandidate> candidates;
    for (int t = 0; t < static_cast<int>(Elems.size()); ++t)
        if (needsTopologyRepair(t,angleDegrees))
            candidates.push_back(topologyCandidate(t));
    std::stable_sort(candidates.begin(),candidates.end(),[&](const TopologyCandidate& a,
        const TopologyCandidate& b) { return Elems[a.tet].q < Elems[b.tet].q; });
    if (candidates.empty()) return 0;
    auto global = measureOptimizationQuality(*this);
    if (!global.valid) return 0;
    const size_t budget = topologyRemeshTrialBudget;
    int accepted = 0;
    for (const auto& candidate : candidates) {
        if (topologyCandidateCurrent(candidate))
            accepted += tryTopologyRepair(candidate,angleDegrees,&global);
    }
    if (infolevel > 1) meshLogger->debug("Topology repair pass: accepted={} remeshTrials={}",accepted,budget-topologyRemeshTrialBudget);
    return accepted;
}
}
