#include "dt.h"
#include <iostream>
#include <stdexcept>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static int physicalSeed(dt::DT& d) {
    for (int t = 0; t < static_cast<int>(d.Elems.size()); ++t)
        if (!d.isDelEle(t) && !d.ishulltet(t) && !d.isvirtualtet(t)) return t;
    return -1;
}

static void checkPolicy(double scale) {
    dt::Mesh mesh;
    mesh.V = {{{0, 0, 0}}, {{scale, 0, 0}}, {{0, scale, 0}}, {{0, 0, scale}}};
    mesh.F = {{{0, 1, 2, 0}}, {{0, 3, 1, 0}}, {{0, 2, 3, 0}}, {{1, 3, 2, 0}}};
    dt::Args args;
    args.refine = 0; args.optlevel = 0; args.nthread = 1; args.infolevel = 0;
    dt::DT d;
    require(d.dt_init(mesh, args) && d.BndPntInst(mesh, args), "tetrahedron setup failed");
    d.BoudaryRecover(mesh, args); d.ColorVirtualTet(args);
    const int seed = physicalSeed(d);
    require(seed >= 0, "missing physical tetrahedron");
    d.improve_step = true;
    d.improve_Metric = dt::DT::SUS_METRIC;
    d.prepareQuality();
    d.minVolume_bw = scale * scale * scale * 1e-12;
    const int node = d.addNode(scale * .25, scale * .25, scale * 1e-9, 1);
    const std::vector<int> seeds(1, seed);
    const size_t elementCount = d.Elems.size();
    for (int info : {1, 3}) {
        dt::BWPlan plan;
        d.bwAcceptance = dt::BWAcceptance::Default;
        auto request = dt::makeBWRequest(d, node, seeds, info);
        request.trackAccess = false;
        require(dt::planBW(d, request, plan) == dt::BWStatus::Rejected,
            "optimization accepted the thin candidate");
        d.bwAcceptance = dt::BWAcceptance::VolumeOnly;
        request = dt::makeBWRequest(d, node, seeds, info);
        request.trackAccess = false;
        // Acceptance belongs to the request, not a mutable flag at planning time.
        d.bwAcceptance = dt::BWAcceptance::Default;
        require(dt::planBW(d, request, plan) == dt::BWStatus::Ready,
            "volume-only refinement inherited the shape-quality threshold");
        double minimum = DBL_MAX;
        for (const auto& face : plan.faces) {
            if (face.geo == dt::DT::virtualID ||
                std::find(face.vertices.begin(), face.vertices.end(), d.ghost) != face.vertices.end()) continue;
            minimum = std::min(minimum, d.quality_sus(request.point.data(),
                d.Nodes[face.vertices[0]].pt, d.Nodes[face.vertices[1]].pt,
                d.Nodes[face.vertices[2]].pt));
        }
        require(minimum > 0 && minimum < .01, "fixture does not exercise low-quality refinement");
        d.minVolume_bw = scale * scale * scale * 1e-6;
        require(dt::planBW(d, request, plan) == dt::BWStatus::Rejected,
            "volume-only refinement ignored its volume bound");
        d.minVolume_bw = scale * scale * scale * 1e-12;
        require(d.Elems.size() == elementCount, "planning modified the mesh");
    }
    d.bwAcceptance = dt::BWAcceptance::VolumeOnly;
    std::vector<int> insertionSeeds = seeds;
    require(d.BW_insert_vertex(node, insertionSeeds, 1) == 1, "volume-only commit failed");
    require(d.improve_step, "refinement disabled quality-cache updates");
    std::vector<int> star;
    d.findSphere(node, star);
    for (int t : star) {
        if (d.ishulltet(t) || d.isvirtualtet(t)) continue;
        const auto& v = d.Elems[t].form;
        const double quality = d.quality_sus(d.Nodes[v[0]].pt, d.Nodes[v[1]].pt,
            d.Nodes[v[2]].pt, d.Nodes[v[3]].pt);
        require(d.Elems[t].q == quality, "refinement left a stale quality cache");
        require(d.calVolume(t) > d.minVolume_bw, "refinement committed a tiny or inverted cell");
        for (int f = 0; f < 4; ++f)
            require(d.getNeig(d.getNeig(t, f), d.getNeigOrd(t, f)) == t, "invalid adjacency");
    }
    d.bwAcceptance = dt::BWAcceptance::Default;
    d.Type_Vertex_Edg(170, mesh);
    d.improve_init(args);
    args.optTh = -1; args.optanglestrict = 0; args.optloop = 3;
    std::unordered_map<int, double> facetSize, elementSize;
    std::unordered_map<int, int> facetNum, elementNum;
    require(d.OptSizeControl_yunbo(0, facetSize, facetNum, elementSize, elementNum, args) == 1,
        "size-control early-exit check failed");
    require(d.bwAcceptance == dt::BWAcceptance::Default && d.improve_step,
        "size-control early exit leaked its acceptance policy");
    d.bwAcceptance = dt::BWAcceptance::VolumeOnly;
    require(d.dt_init(mesh, args) == 1 && d.bwAcceptance == dt::BWAcceptance::Default,
        "mesh reinitialization retained the refinement policy");
}

static void checkExistingFeature(std::string filename) {
    dt::Mesh mesh;
    dt::Args args;
    args.nthread = 1; args.infolevel = 0;
    dt::DT d;
    require(dt::readMesh(filename, mesh) && d.dt_init(mesh, args), "read existing mesh failed");
    d.buildTetInfo(mesh, args);
    d.modifyBnd = true;
    d.Type_Vertex_Edg(170, mesh);
    d.improve_init(args);
    d.improve_Metric = dt::DT::SUS_METRIC;
    d.prepareQuality(); d.updateminVolume();
    int longest = -1;
    double length = 0;
    for (int e = 0; e < static_cast<int>(d.SurEdgs.size()); ++e) {
        if (d.isDelSurEdg(e) || d.SurEdgs[e].constrain <= 0) continue;
        const auto& edge = d.SurEdgs[e];
        const double candidate = d.distance(d.Nodes[edge.iStart].pt, d.Nodes[edge.iEnd].pt);
        if (candidate > length) { length = candidate; longest = e; }
    }
    require(longest >= 0 && length > 800, "missing the original long feature edge");
    require(d.splitEdg(longest) == 0, "optimization baseline unexpectedly accepted the original edge");
    d.bwAcceptance = dt::BWAcceptance::VolumeOnly;
    const int node = d.splitEdg(longest);
    require(node > 0 && d.isDelSurEdg(longest), "volume-only refinement did not split the original edge");
    require(d.isSegmentpnt(node), "new feature point lost its constraint");
    std::cout << "Original feature edge split: length=" << length << " node=" << node << '\n';
}

static void inspectOutput(std::string filename) {
    dt::Mesh mesh;
    dt::Args args;
    args.nthread = 1; args.infolevel = 0;
    dt::DT d;
    require(dt::readMesh(filename, mesh) && d.dt_init(mesh, args), "read result failed");
    d.buildTetInfo(mesh, args);
    d.Type_Vertex_Edg(170, mesh);
    size_t featureCount = 0, oversized = 0, physical = 0, invalid = 0;
    double longest = 0, minimumVolume = DBL_MAX, minimumQuality = DBL_MAX;
    for (int e = 0; e < static_cast<int>(d.SurEdgs.size()); ++e) {
        if (d.isDelSurEdg(e) || d.SurEdgs[e].constrain <= 0) continue;
        const auto& edge = d.SurEdgs[e];
        const double length = d.distance(d.Nodes[edge.iStart].pt, d.Nodes[edge.iEnd].pt);
        ++featureCount; oversized += length > 33.3069 * 4 / 3;
        longest = std::max(longest, length);
    }
    for (int t = 0; t < static_cast<int>(d.Elems.size()); ++t) {
        if (d.isDelEle(t) || d.ishulltet(t) || d.isvirtualtet(t)) continue;
        ++physical;
        const auto& v = d.Elems[t].form;
        const double volume = d.calVolume(t);
        invalid += !std::isfinite(volume) || !(volume > 0);
        minimumVolume = std::min(minimumVolume, volume);
        minimumQuality = std::min(minimumQuality, d.quality_sus(d.Nodes[v[0]].pt,
            d.Nodes[v[1]].pt, d.Nodes[v[2]].pt, d.Nodes[v[3]].pt));
    }
    std::cout << std::setprecision(17) << "Final features: count=" << featureCount
        << " max_length=" << longest << " over_44.4092=" << oversized << '\n'
        << "Final cells: count=" << physical << " invalid_volume=" << invalid
        << " minimum_volume6=" << minimumVolume << " minimum_quality=" << minimumQuality << '\n';
    require(oversized == 0, "result retains an oversized feature edge");
    require(invalid == 0, "result has non-positive or non-finite cells");
}

int main(int argc, char** argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--inspect") {
            inspectOutput(argv[2]);
            return 0;
        }
        for (double scale : {1.0, 1e3, 1e-3}) checkPolicy(scale);
        if (argc > 1) checkExistingFeature(argv[1]);
        std::cout << "Refinement policy passed: optimization guard, volume-only acceptance, "
            "volume rejection, request snapshot, quality cache, adjacency, early-exit restore, reset.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
