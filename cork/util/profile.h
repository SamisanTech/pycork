// Lightweight phase profiler. Enabled at runtime with env CORK_PROFILE=1.
#pragma once
#include <chrono>
#include <cstdio>
#include <cstdlib>

namespace cork_prof {

inline bool enabled() {
    static int e = -1;
    if (e < 0) e = std::getenv("CORK_PROFILE") ? 1 : 0;
    return e == 1;
}

struct Scope {
    const char *name;
    std::chrono::steady_clock::time_point t0;
    explicit Scope(const char *n) : name(n), t0(std::chrono::steady_clock::now()) {}
    ~Scope() {
        if (!enabled()) return;
        double ms = std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - t0).count();
        std::fprintf(stderr, "[cork] %-36s %9.2f ms\n", name, ms);
        std::fflush(stderr);
    }
};

inline void note(const char *name, double value) {
    if (!enabled()) return;
    std::fprintf(stderr, "[cork] %-36s %9.0f\n", name, value);
    std::fflush(stderr);
}

} // namespace cork_prof

#define CORK_PROF_CAT2(a, b) a##b
#define CORK_PROF_CAT(a, b) CORK_PROF_CAT2(a, b)
#define CORK_PROF(name) cork_prof::Scope CORK_PROF_CAT(_cork_prof_, __LINE__)(name)
