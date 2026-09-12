#include "dt.h"
#include "dt_parallel.h"

namespace dt {
int activeDTThreads(const DT& mesh, size_t pointCount, size_t jobs) {
#ifdef _OPENMP
    if (omp_in_parallel()) return 1;
    size_t limit = static_cast<size_t>(std::max(1, std::min(128, mesh.num_threads)));
    limit = std::min(limit, std::max(size_t(1), jobs));
    // Zero/one disables the scale gate for small correctness fixtures.
    if (mesh.parallelMinPointsPerThread > 1) {
        if (pointCount == 0) pointCount = mesh.Nodes.size();
        const size_t bySize = std::max(size_t(1), pointCount / mesh.parallelMinPointsPerThread);
        limit = std::min(limit, bySize);
    }
    int workers = 1;
    while (static_cast<size_t>(workers) <= limit / 2) workers *= 2;
    return workers;
#else
    return 1;
#endif
}

void initializeDTThreads(DT& mesh, const Args& args, size_t pointCount) {
    if (mesh.threadsInitialized) return;
    int available = 1;
#ifdef _OPENMP
    available = std::max(1, omp_get_num_procs());
    mesh.num_threads = std::min(128, std::min(available,
        args.nthread > 0 ? args.nthread : available));
#else
    mesh.num_threads = 1;
#endif
    mesh.threadsInitialized = true;
    if (args.infolevel > 0) {
        mesh.meshLogger->info("DT Thread : {}", mesh.num_threads);
    }
}
}
