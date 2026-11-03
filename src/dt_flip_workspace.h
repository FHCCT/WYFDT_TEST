#pragma once

#include <memory>
#include <vector>

namespace dt { namespace detail {

// A nested flip or a private boundary trial must never borrow its caller's
// buffers. Keep each leased workspace stable when the pool grows.
struct FlipWorkspace {
    std::vector<int> first, second, third;
};

class FlipWorkspaceLease {
    struct Pool {
        size_t depth = 0;
        std::vector<std::unique_ptr<FlipWorkspace>> slots;
    };
    static Pool& pool() {
        static thread_local Pool value;
        return value;
    }
    Pool& owner;
    FlipWorkspace* buffers;
public:
    FlipWorkspaceLease() : owner(pool()) {
        if (owner.depth == owner.slots.size())
            owner.slots.push_back(std::unique_ptr<FlipWorkspace>(new FlipWorkspace()));
        buffers = owner.slots[owner.depth].get();
        ++owner.depth;
        buffers->first.clear();
        buffers->second.clear();
        buffers->third.clear();
    }
    ~FlipWorkspaceLease() { --owner.depth; }
    FlipWorkspace& get() { return *buffers; }
    FlipWorkspaceLease(const FlipWorkspaceLease&) = delete;
    FlipWorkspaceLease& operator=(const FlipWorkspaceLease&) = delete;
};

} }
