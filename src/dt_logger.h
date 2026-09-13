#pragma once
#include "spdlog/spdlog.h"
#include <chrono>

namespace dt {
class DT;
// Copies of a DT may share synchronized output sinks, but never mutable log
// levels. The host owns its default logger; a mesh never reconfigures it.
class InstanceLogger {
    std::shared_ptr<spdlog::logger> value;
public:
    InstanceLogger() { reset(); }
    InstanceLogger(const InstanceLogger& other)
        : value(other.value->clone(other.value->name())) {}
    InstanceLogger& operator=(const InstanceLogger& other) {
        if (this != &other) value = other.value->clone(other.value->name());
        return *this;
    }
    InstanceLogger& operator=(std::shared_ptr<spdlog::logger> logger) {
        value = std::move(logger);
        return *this;
    }
    spdlog::logger* operator->() const { return value.get(); }
    void reset() {
        const auto host = spdlog::default_logger();
        value = host->clone(host->name());
    }
};

// One timer per algorithm stage, outside worker loops. Quiet runs do not read the clock.
class MeshStageLog {
    DT& mesh;
    const char* name;
    spdlog::level::level_enum level;
    bool active;
    std::chrono::steady_clock::time_point started;
public:
    MeshStageLog(DT& mesh, const char* name, int verbosity = 1);
    ~MeshStageLog() noexcept;
    MeshStageLog(const MeshStageLog&) = delete;
    MeshStageLog& operator=(const MeshStageLog&) = delete;
};
}
