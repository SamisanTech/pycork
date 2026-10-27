#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include <cstdio>

#include "seh_guard.h"

int cork_seh_invoke(void (*fn)(void *), void *arg) {
    __try {
        fn(arg);
        return 0;
    } __except (EXCEPTION_EXECUTE_HANDLER) {
        const unsigned code = (unsigned)GetExceptionCode();
        std::fprintf(stderr, "fast_cork: caught native exception 0x%08X\n", code);
        std::fflush(stderr);
        return code ? (int)code : 1;
    }
}
