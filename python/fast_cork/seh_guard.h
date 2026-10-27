#pragma once

// Run fn(arg). Returns 0 when fn returns. Returns the SEH exception code
// when fn hits an access violation or similar native fault.
// fn must catch C++ exceptions itself; this function has no C++ objects.
int cork_seh_invoke(void (*fn)(void *), void *arg);
