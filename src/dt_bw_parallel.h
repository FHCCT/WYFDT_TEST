#pragma once

#include <array>
#include <algorithm>
#include <cstdint>
#include <vector>

namespace dt {
class DT;
struct Args;

// Scheduler policy shared by all BW entry points.
constexpr int BW_MIN_REGIONS_PER_THREAD = 8;
constexpr int BW_MAX_REGIONS_PER_THREAD = 16;
constexpr int BW_MAX_REGIONS = 256;
constexpr int BW_TARGET_POINTS_PER_REGION = 4096;
constexpr int BW_MIN_POINTS_PER_REGION = 64;
constexpr int BW_SCAFFOLD_SAMPLES = 8;
constexpr int BW_SERIAL_TAIL_DIVISOR = 4;

// A high requested thread count is a ceiling, not a target. Decisions use
// remaining work and measured independent insertions, never just node slots.
class BWRefineTeamPolicy {
    int ceiling = 128;
    size_t cooldown = 0;
    uint64_t ready = 0, inserted = 0, conflicts = 0, usefulWorkTarget = 0;
    uint64_t work = 0, wastedWork = 0;
    int samples = 0;
public:
    int choose(int cap, size_t remaining, int lanes, size_t workPerThread) const {
        if (cooldown) return 1;
        size_t limit = std::min(static_cast<size_t>(std::max(1, std::min(cap, ceiling))),
            std::max(size_t(1), remaining / std::max(size_t(1), workPerThread)));
        limit = std::min(limit, static_cast<size_t>(std::max(1, lanes)));
        int team = 1;
        while (static_cast<size_t>(team) <= limit / 2) team *= 2;
        return team;
    }
    bool observe(int team, uint64_t candidates, uint64_t accepted, uint64_t deferred,
                 uint64_t estimatedWork = 0, uint64_t deferredWork = 0) {
        if (team <= 1) return false;
        ready += candidates; inserted += accepted; conflicts += deferred; ++samples;
        work += estimatedWork; wastedWork += deferredWork;
        usefulWorkTarget += static_cast<uint64_t>(team * 2);
        // An isolated conflict does not change the team. Sustained contention
        // or too few independent jobs gets a real serial progress interval.
        if (samples < 8) return false;
        const bool congested = ready >= 16 &&
            (conflicts * 2 >= ready || (work && wastedWork >= work / 4));
        const bool sparse = inserted < usefulWorkTarget;
        ready = inserted = conflicts = usefulWorkTarget = work = wastedWork = 0; samples = 0;
        if (congested || sparse) { ceiling = 1; cooldown = 512; return true; }
        ceiling = std::min(128, std::max(2, ceiling * 2));
        return false;
    }
    void serialProgress(size_t attempts) {
        if (!cooldown) return;
        if (attempts >= cooldown) { cooldown = 0; ceiling = 2; }
        else cooldown -= attempts;
    }
};

void insertDelaunayPoints(DT& mesh, const std::vector<int>& order);
void refineBWParallel(DT& mesh, Args& args);
// Conservative lower bound used only for early candidate rejection.
double refineBWSizeLowerBound(DT& mesh, int tet);
}
