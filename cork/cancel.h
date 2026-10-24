// Cooperative cancel for long-running Cork SI / boolean ops.
#ifndef CORK_CANCEL_H_HEADER_HAS_BEEN_INCLUDED
#define CORK_CANCEL_H_HEADER_HAS_BEEN_INCLUDED

#include <atomic>
#include <stdexcept>

#include "CORK_Export.h"

struct CorkCancelled : public std::runtime_error {
	CorkCancelled() : std::runtime_error("cancelled") {}
};

CORK_EXPORT void cork_request_cancel();
CORK_EXPORT void cork_reset_cancel();
CORK_EXPORT bool cork_is_cancelled();

inline void cork_throw_if_cancelled()
{
	if (cork_is_cancelled()) {
		throw CorkCancelled();
	}
}

// Hot-path helper: local counter + power-of-two stride so we almost never
// touch the atomic. Default 4096 ≈ sub-ms abort latency on SI resolve.
inline void cork_throw_if_cancelled_strided(unsigned &counter,
                                           unsigned stride_mask = 4095u)
{
	if ((++counter & stride_mask) == 0u) {
		cork_throw_if_cancelled();
	}
}

#endif
