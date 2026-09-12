#pragma once

#include <cstddef>
#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace dt {
// OpenMP control variables belong to the calling task. Bound the temporary
// dynamic-team override to the algorithm call and restore it on exceptions too.
class DTParallelScope {
#ifdef _OPENMP
    const int dynamic = omp_get_dynamic();
public:
    DTParallelScope() { if (dynamic) omp_set_dynamic(0); }
    ~DTParallelScope() { if (dynamic) omp_set_dynamic(dynamic); }
#else
public:
    DTParallelScope() = default;
#endif
    DTParallelScope(const DTParallelScope&) = delete;
    DTParallelScope& operator=(const DTParallelScope&) = delete;
};

class DT;
struct Args;

// num_threads is the only user-configured upper bound. Each stage chooses a
// smaller power-of-two team from its point count and available independent jobs.
int activeDTThreads(const DT& mesh, size_t pointCount = 0, size_t jobs = SIZE_MAX);
void initializeDTThreads(DT& mesh, const Args& args, size_t pointCount = 0);
}
