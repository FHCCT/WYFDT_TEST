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
    const OptimizationQuality& before, const DT::TopologyCandidate& seed, detail::QualityBudget* budget) {
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
    detail::QualityTotals budgetBefore, budgetAfter;
    if (budget) budgetBefore = budget->cells(host,cells);
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
        if (budget) budget->add(host,budgetAfter,a,b,c,d);
        Elem element(v[0],v[1],v[2],v[3]); element.geo = host.Elems[seed.tet].geo;
        element.q = host.tetquality(a,b,c,d,nullptr,host.improve_Metric,false);
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
    if (!budget && host.infolevel > 1) host.meshLogger->debug(
        "Local remesh quality: minA={:.6g}->{:.6g} minSJ={:.6g}->{:.6g} meanA={:.6g}->{:.6g} meanSJ={:.6g}->{:.6g}",
        before.minAngle,after.minAngle,before.minJacobian,after.minJacobian,
        before.meanAngle(),after.meanAngle(),before.meanJacobian(),after.meanJacobian());
    // Accept by the cavity's actual worst angle and scaled Jacobian: neither
    // may decrease and at least one must improve by more than one percent.
    // Local means/counts are diagnostics, not a veto on worst-cell repair.
    // Near zero, an absolute tolerance must not mask a large relative loss.
    const double angleTolerance = std::min(1e-10,before.minAngle*1e-10);
    const double jacobianTolerance = std::min(1e-12,before.minJacobian*1e-10);
    if (after.minAngle + angleTolerance < before.minAngle ||
        after.minJacobian + jacobianTolerance < before.minJacobian ||
        !(after.minAngle > before.minAngle * 1.01 + 1e-10 ||
          after.minJacobian > before.minJacobian * 1.01 + 1e-12)) return false;
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

    if (budget && !budget->accepts(budgetBefore,budgetAfter)) return false;

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
    if (budget) budget->commit(budgetBefore,budgetAfter);
    if (host.QuantityControl) regionCount->second += static_cast<int>(count)-static_cast<int>(cells.size());
    for (int n : toGlobal) if (n < static_cast<int>(host.susIdleStates.size()))
        host.susIdleStates[n] = DT::SusIdleState{};
    if (!budget && host.infolevel > 1) host.meshLogger->debug(
        "Local remesh: cells={}->{} minAngle={:.6g}->{:.6g} minSJ={:.6g}->{:.6g}",
        cells.size(),count,before.minAngle,after.minAngle,before.minJacobian,after.minJacobian);
    return true;
}
}

int DT::remeshTopologyPatch(int tet, detail::QualityBudget* budget) {
    if (omp_in_parallel())
        throw std::logic_error("Local remeshing must run outside a parallel region");
    if (!improve_step || !AniSol.empty() || topologyRemeshTrialBudget == 0 || tet < 0 ||
        tet >= static_cast<int>(Elems.size()) || isDelEle(tet) || ishulltet(tet) ||
        isvirtualtet(tet)) return 0;
    --topologyRemeshTrialBudget;
    const auto seed = topologyCandidate(tet);
    const int region = Elems[tet].geo;
    constexpr size_t maxCells = 192, maxNodes = 112;
    std::unique_ptr<DT> trial;
    std::set<std::vector<int>> tried;
    auto physical = [&](int t) {
        return t >= 0 && t < static_cast<int>(Elems.size()) && !isDelEle(t) &&
            !ishulltet(t) && !isvirtualtet(t) && Elems[t].geo == region;
    };
    auto canAdd = [&](int t, const std::set<int>& included) {
        if (!physical(t) || included.count(t)) return false;
        for (int f = 0; f < 4; ++f) {
            const Face face = faceKey(Elems[t].form,f);
            if (included.count(getNeig(t,f)) && isBndTri(face[0],face[1],face[2])) return false;
        }
        return true;
    };
    auto expand = [&](std::vector<int>& cells, std::set<int>& included,
                      const std::set<int>* allowed) {
        const size_t original = cells.size();
        for (size_t i = 0; i < original; ++i) {
            const int t = cells[i];
            for (int f = 0; f < 4; ++f) {
                const Face face = faceKey(Elems[t].form,f);
                if (isBndTri(face[0],face[1],face[2])) continue;
                const int n = getNeig(t,f);
                if (cells.size() >= maxCells || (allowed && !allowed->count(n)) || !canAdd(n,included)) continue;
                included.insert(n); cells.push_back(n);
            }
        }
        return cells.size() > original;
    };
    auto attempt = [&](const std::vector<int>& cells, const std::set<int>& included) {
        auto signature = cells; std::sort(signature.begin(),signature.end());
        if (!tried.insert(signature).second) return false;
        Mesh mesh;
        std::map<int,int> local;
        std::vector<int> toGlobal;
        for (int t : cells) for (int n : Elems[t].form)
            if (local.emplace(n,static_cast<int>(toGlobal.size())).second) {
                toGlobal.push_back(n);
                mesh.V.push_back({{Nodes[n].pt[0],Nodes[n].pt[1],Nodes[n].pt[2]}});
            }
        if (toGlobal.size() > maxNodes) return false;
        std::map<Face,PatchBoundary> boundary;
        for (int t : cells) for (int f = 0; f < 4; ++f) {
            const Face face = faceKey(Elems[t].form,f);
            if (included.count(getNeig(t,f))) {
                if (isBndTri(face[0],face[1],face[2])) return false;
                continue;
            }
            if (!boundary.emplace(face,PatchBoundary{t,f,getNeig(t,f),getNeigOrd(t,f)}).second) return false;
            int ia,ib,ic,id; DNC(f,ia,ib,ic,id);
            mesh.F.push_back({{local[Elems[t].form[ib]],local[Elems[t].form[ic]],
                local[Elems[t].form[id]],static_cast<int>(mesh.F.size())}});
        }
        const auto before = measureOptimizationQuality(*this,&cells);
        if (!before.valid || boundary.empty()) return false;
        std::map<std::array<int,2>,int> edgeUses;
        for (const auto& entry : boundary) for (int k = 0; k < 3; ++k) {
            std::array<int,2> edge{{entry.first[k],entry.first[(k+1)%3]}};
            std::sort(edge.begin(),edge.end()); ++edgeUses[edge];
        }
        for (const auto& edge : edgeUses) if (edge.second != 2) return false;
        if (!trial) {
            trial.reset(new DT());
            trial->meshLogger = std::make_shared<spdlog::logger>("local_remesh_trial",
                std::make_shared<spdlog::sinks::null_sink_mt>());
        }
        // Under a quality ledger, try deep flips privately before regenerating
        // this same constrained patch. Export may compact points; replacePatch
        // always matches actual host coordinates and validates the whole cavity.
        if (budget) {
            try {
                Mesh flipMesh = mesh;
                for (int cell : cells) {
                    const auto& v = Elems[cell].form;
                    flipMesh.T.push_back({{local[v[0]],local[v[1]],local[v[2]],local[v[3]],1}});
                }
                Args flipArgs; flipArgs.infolevel=0; flipArgs.outlogfile=0; flipArgs.nthread=1; flipArgs.optlevel=1;
                if (trial->dt_init(flipMesh,flipArgs)) {
                    trial->buildTetInfo(flipMesh,flipArgs); trial->improve_init(flipArgs);
                    trial->modifyBnd=false; trial->improve_Metric=improve_Metric;
                    trial->optflipdeep=2; trial->prepareQuality();
                    int localSeed=0;
                    if (trial->removebadtet(localSeed,-1)==1) {
                        trial->RemoveTet(flipArgs); trial->outMesh(flipMesh,flipArgs);
                        if (replacePatch(*this,cells,flipMesh,toGlobal,boundary,before,seed,budget)) return true;
                    }
                }
            } catch (const std::bad_alloc&) { throw; }
            catch (...) { trial.reset(); }
            if (!trial) {
                trial.reset(new DT());
                trial->meshLogger = std::make_shared<spdlog::logger>("local_remesh_trial",
                    std::make_shared<spdlog::sinks::null_sink_mt>());
            }
        }
        Args args; args.infolevel = 0; args.refine = 0; args.optlevel = 0;
        args.nthread = 1; args.constrain = 1; args.outlogfile = 0;
        bool generated = false;
        try { generated = trial->tetrahedralize(mesh,args) != 0; }
        catch (const std::bad_alloc&) { throw; }
        catch (...) { trial.reset(); }
        if (!budget && infolevel > 1) meshLogger->debug(
            "Local remesh trial: cells={}->{} points={}->{} generated={}",
            cells.size(),mesh.T.size(),toGlobal.size(),mesh.V.size(),generated);
        if (generated && replacePatch(*this,cells,mesh,toGlobal,boundary,before,seed,budget)) return true;
        return false;
    };
    // Start with small edge stars; a whole face ball can dilute a repair
    // with unrelated cells. Each star stays in the seed's constrained sector.
    for (int edge = 0; edge < 6; ++edge) {
        const int a = seed.form[Egid[edge][0]], b = seed.form[Egid[edge][1]];
        std::vector<int> star; findSphere(a,star);
        std::set<int> allowed;
        for (int t : star) if (physical(t) &&
            std::find(Elems[t].form,Elems[t].form+4,b) != Elems[t].form+4)
            allowed.insert(t);
        std::vector<int> cells(1,tet); std::set<int> included{tet};
        while (expand(cells,included,&allowed)) {}
        if (attempt(cells,included)) return 1;
    }
    // If small stars fail, try face rings two through six with fixed limits.
    std::vector<int> cells(1,tet); std::set<int> included{tet};
    for (int ring = 1; ring <= 6; ++ring) {
        const bool grew = expand(cells,included,nullptr);
        if (ring > 1 && attempt(cells,included)) return 1;
        if (!grew || cells.size() >= maxCells) break;
    }
    return 0;
}

}
