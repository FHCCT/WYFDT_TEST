#pragma once
#include "spdlog/spdlog.h"

namespace dt {
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
}
