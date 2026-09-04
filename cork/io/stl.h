#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include <cork/core/surface.h>
#include <cork/util/parallel.h>
#include <cork/util/profile.h>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

namespace cork {
namespace io {

#pragma pack(push, 1)
struct StlTri {
    float n[3];
    float v[3][3];
    std::uint16_t attr;
};
#pragma pack(pop)
static_assert(sizeof(StlTri) == 50, "binary STL triangle is 50 bytes");

#ifdef _WIN32
inline std::wstring utf8_to_wide(const std::string &s) {
    if (s.empty()) return std::wstring();
    int n = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), (int)s.size(), nullptr, 0);
    std::wstring w((size_t)n, L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.c_str(), (int)s.size(), &w[0], n);
    return w;
}
#endif

inline const char *map_read(const std::string &path, size_t &n, void **keep) {
#ifdef _WIN32
    std::wstring w = utf8_to_wide(path);
    HANDLE f = CreateFileW(w.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                           OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, nullptr);
    if (f == INVALID_HANDLE_VALUE)
        throw std::runtime_error("readSTL: cannot open " + path);
    LARGE_INTEGER sz;
    GetFileSizeEx(f, &sz);
    n = (size_t)sz.QuadPart;
    HANDLE m = CreateFileMappingW(f, nullptr, PAGE_READONLY, 0, 0, nullptr);
    const char *p = m ? (const char *)MapViewOfFile(m, FILE_MAP_READ, 0, 0, 0) : nullptr;
    CloseHandle(f);
    if (!p) {
        if (m) CloseHandle(m);
        throw std::runtime_error("readSTL: map failed " + path);
    }
    keep[0] = m;
    keep[1] = (void *)p;
    return p;
#else
    (void)keep;
    throw std::runtime_error("readSTL: mmap not implemented");
#endif
}

inline char *map_write(const std::string &path, size_t n, void **keep) {
#ifdef _WIN32
    std::wstring w = utf8_to_wide(path);
    HANDLE f = CreateFileW(w.c_str(), GENERIC_READ | GENERIC_WRITE, 0, nullptr,
                           CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (f == INVALID_HANDLE_VALUE)
        throw std::runtime_error("writeSTL: cannot open " + path);
    LARGE_INTEGER sz;
    sz.QuadPart = (LONGLONG)n;
    SetFilePointerEx(f, sz, nullptr, FILE_BEGIN);
    SetEndOfFile(f);
    HANDLE m = CreateFileMappingW(f, nullptr, PAGE_READWRITE, sz.HighPart, sz.LowPart, nullptr);
    char *p = m ? (char *)MapViewOfFile(m, FILE_MAP_WRITE, 0, 0, n) : nullptr;
    CloseHandle(f);
    if (!p) {
        if (m) CloseHandle(m);
        throw std::runtime_error("writeSTL: map failed " + path);
    }
    keep[0] = m;
    keep[1] = p;
    return p;
#else
    (void)path;
    (void)n;
    (void)keep;
    throw std::runtime_error("writeSTL: mmap not implemented");
#endif
}

inline void unmap(void **keep) {
#ifdef _WIN32
    if (keep[1]) UnmapViewOfFile(keep[1]);
    if (keep[0]) CloseHandle((HANDLE)keep[0]);
#endif
    keep[0] = keep[1] = nullptr;
}

inline std::uint32_t hash3(std::uint32_t a, std::uint32_t b, std::uint32_t c) {
    std::uint32_t h = a * 0x9e3779b1u;
    h ^= b + 0x9e3779b1u + (h << 6) + (h >> 2);
    h ^= c + 0x9e3779b1u + (h << 6) + (h >> 2);
    return h;
}

// Native STL is float32. Weld on exact bits. Caller casts if needed.
inline Surface<float> read_stl(const std::string &path) {
    CORK_PROF("io.read_stl");
    size_t nbytes = 0;
    void *keep[2] = {nullptr, nullptr};
    const char *buf = map_read(path, nbytes, keep);
    if (nbytes < 84) {
        unmap(keep);
        throw std::runtime_error("readSTL: file too small " + path);
    }
    std::uint32_t ntri = 0;
    std::memcpy(&ntri, buf + 80, 4);
    if (nbytes < 84ull + (size_t)ntri * 50ull) {
        unmap(keep);
        throw std::runtime_error("readSTL: truncated binary STL " + path);
    }

    const StlTri *rec = reinterpret_cast<const StlTri *>(buf + 84);
    const size_t ncorner = (size_t)ntri * 3;
    size_t cap = 1;
    while (cap < ((size_t)ntri < 8 ? 8 : (size_t)ntri)) cap <<= 1;
    if (cap < 1024) cap = 1024;
    struct Slot {
        std::uint32_t h, id;
    };
    std::vector<Slot> tab(cap);
    std::memset(tab.data(), 0xff, cap * sizeof(Slot));
    std::uint32_t mask = (std::uint32_t)cap - 1;

    std::vector<std::uint32_t> bits;
    bits.reserve(((ncorner >> 1) + 16) * 3);
    Surface<float> out;
    out.f.resize(ncorner);

    std::uint32_t next = 0;
    auto grow = [&]() {
        size_t ncap = cap << 1;
        std::vector<Slot> ntab(ncap);
        std::memset(ntab.data(), 0xff, ncap * sizeof(Slot));
        const std::uint32_t nmask = (std::uint32_t)ncap - 1;
        for (size_t i = 0; i < cap; ++i) {
            if (tab[i].id == 0xffffffffu) continue;
            std::uint32_t j = tab[i].h & nmask;
            while (ntab[j].id != 0xffffffffu) j = (j + 1) & nmask;
            ntab[j] = tab[i];
        }
        tab.swap(ntab);
        cap = ncap;
        mask = nmask;
    };

    for (std::uint32_t t = 0; t < ntri; ++t) {
        for (int k = 0; k < 3; ++k) {
            std::uint32_t b0, b1, b2;
            std::memcpy(&b0, &rec[t].v[k][0], 4);
            std::memcpy(&b1, &rec[t].v[k][1], 4);
            std::memcpy(&b2, &rec[t].v[k][2], 4);
            const std::uint32_t h = hash3(b0, b1, b2);
            std::uint32_t i = h & mask;
            for (;;) {
                Slot &s = tab[i];
                if (s.id == 0xffffffffu) {
                    if (next * 4 > (std::uint32_t)cap * 3) {
                        grow();
                        i = h & mask;
                        continue;
                    }
                    s.h = h;
                    s.id = next;
                    bits.push_back(b0);
                    bits.push_back(b1);
                    bits.push_back(b2);
                    out.f[t * 3 + k] = next++;
                    break;
                }
                if (s.h == h) {
                    const std::uint32_t *p = &bits[s.id * 3];
                    if (p[0] == b0 && p[1] == b1 && p[2] == b2) {
                        out.f[t * 3 + k] = s.id;
                        break;
                    }
                }
                i = (i + 1) & mask;
            }
        }
    }
    unmap(keep);

    out.xyz.resize((size_t)next * 3);
    const float *fp = reinterpret_cast<const float *>(bits.data());
    for (size_t i = 0, n = out.xyz.size(); i < n; ++i) out.xyz[i] = fp[i];
    return out;
}

template <class S, class I>
inline void write_stl(const std::string &path, const S *vp, const I *fp,
                      std::uint32_t ntri) {
    CORK_PROF("io.write_stl");
    const size_t nbytes = 84ull + (size_t)ntri * 50ull;
    void *keep[2] = {nullptr, nullptr};
    char *buf = map_write(path, nbytes, keep);
    std::memset(buf, 0, 84);
    std::memcpy(buf + 80, &ntri, 4);
    StlTri *rec = reinterpret_cast<StlTri *>(buf + 84);
    cork_par::for_range((size_t)ntri, 4096, [&](size_t b, size_t e) {
        for (size_t t = b; t < e; ++t) {
            const size_t a = (size_t)fp[t * 3 + 0];
            const size_t bb = (size_t)fp[t * 3 + 1];
            const size_t c = (size_t)fp[t * 3 + 2];
            rec[t].n[0] = rec[t].n[1] = rec[t].n[2] = 0.f;
            rec[t].attr = 0;
            rec[t].v[0][0] = (float)vp[a * 3 + 0];
            rec[t].v[0][1] = (float)vp[a * 3 + 1];
            rec[t].v[0][2] = (float)vp[a * 3 + 2];
            rec[t].v[1][0] = (float)vp[bb * 3 + 0];
            rec[t].v[1][1] = (float)vp[bb * 3 + 1];
            rec[t].v[1][2] = (float)vp[bb * 3 + 2];
            rec[t].v[2][0] = (float)vp[c * 3 + 0];
            rec[t].v[2][1] = (float)vp[c * 3 + 1];
            rec[t].v[2][2] = (float)vp[c * 3 + 2];
        }
    });
    unmap(keep);
}

}  // namespace io
}  // namespace cork
