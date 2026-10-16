#include "dt_refinement_guard.h"
#include <iostream>
#include <stdexcept>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static void check(double scale, bool feature) {
    dt::Mesh mesh;
    mesh.V = {{{0, 0, 0}}, {{scale, 0, 0}}, {{0, scale, 0}}, {{0, 0, scale}}};
    mesh.T = {{{0, 1, 2, 3, 42}}};
    mesh.F = {{{0, 1, 2, 7}}, {{0, 3, 1, 7}}, {{0, 2, 3, 7}}, {{1, 3, 2, 7}}};
    dt::Args args;
    args.nthread = 1; args.infolevel = 0;
    dt::DT d;
    require(d.dt_init(mesh, args) == 1, "initialize fixture");
    d.buildTetInfo(mesh, args);
    d.modifyBnd = true;
    d.Type_Vertex_Edg(170, mesh);
    d.improve_init(args);
    d.improve_Metric = dt::DT::SUS_METRIC;
    d.prepareQuality();
    d.bwAcceptance = dt::BWAcceptance::VolumeOnly;
    d.minVolume_bw = 1e-12 * scale * scale * scale;
    const int* entry = d.BndEdg.find(0, 1);
    require(entry != nullptr, "missing fixture edge");
    const int edge = *entry;
    d.SurEdgs[edge].constrain = feature ? 1 : 0;
    dt::RefinementBoundaryEdges index(d);
    double nearPoint[] = {scale * .25, scale * .00002, scale * .00003};
    require(index.nearby(nearPoint, scale, 42) == edge, "near-edge candidate was not redirected");
    require(index.nearby(nearPoint, scale, 99) == -1, "candidate redirected to a different region");
    require(index.nearby(nearPoint, scale * .01, 42) == -1, "local insertion scale was ignored");
    double formerlyNear[] = {scale * .25, scale * .002, scale * .003};
    require(index.nearby(formerlyNear, scale, 42) == -1,
        "moderately close candidate still redirected under the strict threshold");
    double justInside[] = {scale * .25, scale * .0009, 0};
    double justOutside[] = {scale * .25, scale * .0011, 0};
    require(index.nearby(justInside, scale, 42) == edge, "strict near threshold missed");
    require(index.nearby(justOutside, scale, 42) == -1, "strict threshold exceeded");
    double farPoint[] = {scale * .25, scale * .25, scale * .25};
    require(index.nearby(farPoint, scale, 42) == -1, "far candidate was redirected");
    double endpointOnly[] = {-scale * .00002, -scale * .00002, -scale * .00002};
    require(index.nearby(endpointOnly, scale, 42) == -1, "endpoint-only proximity split an edge");
    d.lockE.insert(edge);
    require(index.protectedEdge(edge), "locked edge was not protected");
    d.lockE.erase(edge);
    const int face = d.SurEdgs[edge].face[0];
    d.lockF.insert(face);
    require(index.protectedEdge(edge), "locked face was not protected");
    d.lockF.erase(face);
    const int midpoint = d.splitEdg(edge);
    require(midpoint > 0 && d.isDelSurEdg(edge), "boundary redirection failed to split the edge");
    require(d.isbndpnt(midpoint), "redirection inserted an interior point");
    const int* child = d.BndEdg.find(0, midpoint);
    require(child && d.SurEdgs[*child].constrain == (feature ? 1 : 0), "child lost edge classification");
    double nearChild[] = {scale * .125, scale * .00002, scale * .00003};
    require(index.nearby(nearChild, scale, 42) == *child,
        "index returned a deleted edge or missed an appended child");
    std::vector<int> star;
    d.findSphere(midpoint, star);
    for (int t : star) {
        if (d.isvirtualtet(t) || d.ishulltet(t)) continue;
        require(d.calVolume(t) > 0, "redirection produced a non-positive cell");
        for (int f = 0; f < 4; ++f)
            require(d.getNeig(d.getNeig(t, f), d.getNeigOrd(t, f)) == t, "redirection broke adjacency");
    }
}

static void checkExample(std::string filename) {
    dt::Mesh mesh;
    dt::Args args;
    args.nthread = 1; args.infolevel = 0;
    dt::DT d;
    require(dt::readMesh(filename, mesh) && d.dt_init(mesh, args), "read example");
    require(mesh.T.size() > 781165, "example index out of range");
    const auto example = mesh.T[781165];
    require(example[0] == 192842 && example[1] == 248856 && example[2] == 249810 && example[3] == 256660,
        "example mesh changed; tetrahedron IDs no longer match the reported case");
    d.buildTetInfo(mesh, args);
    d.Type_Vertex_Edg(170, mesh);
    require(d.isbndpnt(192842) && d.isbndpnt(248856) && d.isbndpnt(249810)
        && !d.isbndpnt(256660), "example should have three boundary vertices and one interior vertex");
    dt::RefinementBoundaryEdges index(d);
    const double length = d.distance(d.Nodes[192842].pt, d.Nodes[248856].pt);
    const int found = index.nearby(d.Nodes[256660].pt, length, example[4]);
    const int* expected = d.BndEdg.find(192842, 248856);
    require(expected && found == *expected, "reported thin-cell point did not trigger redirection");
    require(d.SurEdgs[found].constrain == 0, "reported edge should be an ordinary boundary edge");
    double* a = d.Nodes[example[0]].pt;
    double* b = d.Nodes[example[1]].pt;
    double* c = d.Nodes[example[2]].pt;
    double* p = d.Nodes[example[3]].pt;
    std::cout << std::setprecision(17) << "Example tetrahedron (zero-based) 781165: quality="
        << d.quality_sus(a, b, c, p) << " volume6=" << d.calVolume(781165)
        << " redirected_edge=" << found << " ordinary_boundary=true\n";
}

int main(int argc, char** argv) {
    try {
        for (double scale : {1.0, 1e3, 1e-3}) for (bool feature : {false, true}) check(scale, feature);
        if (argc > 1) checkExample(argv[1]);
        std::cout << "Near-boundary refinement checks passed: 3 scales, ordinary/feature edges, "
            "region/endpoint/lock checks, dynamic children, positive volume and adjacency.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
