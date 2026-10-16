#include "dt_smoothing_guard.h"
#include <iostream>
#include <iomanip>

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    dt::Mesh mesh;
    std::string input = argv[1];
    if (!dt::readVTK(input, mesh)) return 3;
    dt::Args args;
    args.nthread = argc > 3 ? std::atoi(argv[3]) : 1;
    args.infolevel = 1;
    dt::DT d;
    if (!d.dt_init(mesh, args)) return 4;
    d.buildTetInfo(mesh, args);
    d.modifyBnd = true;
    d.Type_Vertex_Edg(170, mesh);
    d.improve_init(args);
    d.improve_Metric = argc > 2 ? std::atoi(argv[2]) : dt::DT::SUS_METRIC;
    d.prepareQuality();
    std::vector<int> cells, nodes;
    for (int t = 0; t < static_cast<int>(d.Elems.size()); ++t)
        if (!d.isDelEle(t) && !d.ishulltet(t) && !d.isvirtualtet(t)) cells.push_back(t);
    for (int n = 0; n < static_cast<int>(d.Nodes.size()); ++n)
        if (n != d.ghost && !d.isDelNod(n)) nodes.push_back(n);
    size_t improved = 0, unchanged = 0;
    for (int round = 0; round < 3; ++round) {
        std::vector<std::array<double, 3>> positions(nodes.size());
        std::vector<double> qualities;
        for (size_t i = 0; i < nodes.size(); ++i)
            std::copy(d.Nodes[nodes[i]].pt, d.Nodes[nodes[i]].pt + 3, positions[i].begin());
        for (int t : cells) qualities.push_back(d.Elems[t].q);
        const auto before = dt::measureSmoothingQuality(d, cells);
        d.SmoothPass(1, 0.5);
        const auto after = dt::measureSmoothingQuality(d, cells);
        bool moved = false, cacheChanged = false;
        for (size_t i = 0; i < nodes.size(); ++i)
            moved = moved || !std::equal(positions[i].begin(), positions[i].end(), d.Nodes[nodes[i]].pt);
        for (size_t i = 0; i < cells.size(); ++i)
            cacheChanged = cacheChanged || qualities[i] != d.Elems[cells[i]].q;
        if (moved) {
            if (!dt::smoothingImprovesQuality(before, after)) return 5;
            ++improved;
        } else {
            if (cacheChanged || !dt::smoothingPreservesQuality(before, after)) return 6;
            ++unchanged;
        }
        std::cout << std::setprecision(15) << "SMOOTH_PASS round=" << round
            << " mean_shape=" << before.meanShape() << "->" << after.meanShape()
            << " mean_angle_deg=" << RADIO2ANGLE(before.meanAngle())
            << "->" << RADIO2ANGLE(after.meanAngle()) << " unchanged=" << !moved << '\n';
    }
    // This optional case is the recorded angle-only regression on the supplied mesh.
    if (argc > 4 && std::atoi(argv[4]) != 0 && unchanged != 3) return 7;
    std::cout << "SMOOTH_ACCEPTANCE metric=" << d.improve_Metric << " improved=" << improved
        << " unchanged=" << unchanged << " regressions=0 dirty_rejections=0\n";
    return 0;
}
