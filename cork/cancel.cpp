#include "cancel.h"

static std::atomic<bool> g_cork_cancel{false};

void cork_request_cancel()
{
	g_cork_cancel.store(true, std::memory_order_relaxed);
}

void cork_reset_cancel()
{
	g_cork_cancel.store(false, std::memory_order_relaxed);
}

bool cork_is_cancelled()
{
	return g_cork_cancel.load(std::memory_order_relaxed);
}
