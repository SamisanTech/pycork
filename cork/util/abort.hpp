#pragma once

#include <stdexcept>

namespace cork {

// Shewchuk Triangle calls triexit() -> exit() on a hostile PSLG.
// Throw instead so the repair entry can catch it inside this DLL.
// Do not let this exception cross the memimporter boundary.
struct TriangleAbort : std::runtime_error {
    int status;
    explicit TriangleAbort(int status)
        : std::runtime_error("triangle aborted"), status(status ? status : 1) {}
};

}  // namespace cork
