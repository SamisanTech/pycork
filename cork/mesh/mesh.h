// +-------------------------------------------------------------------------
// | mesh.h
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2013
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the Cork library.
// |
// |    Cork is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    Cork is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with Cork.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#ifndef CORK_MESH_H_HEADER_HAS_BEEN_INCLUDED
#define CORK_MESH_H_HEADER_HAS_BEEN_INCLUDED

// SIMPLE USAGE:
//  In order to get standard template inclusion behavior/usage
//  just include this file.  Then the entire template code
//  will be included in the compilation unit.
// ADVANCED USAGE:
//  Only include "mesh.decl.h" where-ever you would normally include a
//  header file.  This will avoid including the implementation code in
//  the current compilation unit.
//  Then, create a seperate cpp file which includes "mesh.h" and
//  explicitly instantiates the template with the desired template
//  parameters.
//  By following this scheme, you can prevent re-compiling the entire
//  template implementation in every usage compilation unit and every
//  time those compilation units are recompiled during development.

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>
#include <set>
#include <sstream>

#include <map>

#include <cork/accel/aabvh.h>
#include <cork/accel/lbvh.h>

#include <cork/isct/empty3d.h>
#include <cork/isct/quantization.h>
#include <cork/isct/unsafeRayTriIsct.h>
#include <cork/isct/triangle.h>

#include <cork/math/bbox.h>
#include <cork/math/vec.h>
#include <cork/math/ray.h>

#include <cork/rawmesh/rawMesh.h>

#include <cork/util/iterPool.h>
#include <cork/util/memPool.h>
#include <cork/util/prelude.h>
#include <cork/util/profile.h>
#include <cork/util/shortVec.h>

#include <cork/util/unionFind.h>

#define REAL double

extern "C" {
    #include <cork/isct/triangle.h>
}

namespace cork_hull { struct HullStats; }


struct BoolVertexData {
};

struct BoolTriangleData {
    byte bool_alg_data; // internal use by algorithm
    // please copy value when the triangle is subdivided
};

template<class VertData, class TriData>
struct IsctVertEdgeTriInput
{
    VertData*   e[2];
    VertData*   t[3];
};

template<class VertData, class TriData>
struct IsctVertTriTriTriInput
{
    VertData*   t[3][3];
};

template<class VertData, class TriData>
struct SubdivideTriInput
{
    TriData*    pt;
    VertData*   pv[3];
    VertData*   v[3];
};

// in order to perform intersections, VertData and TriData must support
struct IsctVertexData {
    // specify how to compute new data for vertices formed by intersections
    /*
    // vertices on edge and triangle forming an intersection...
    void isct(IsctVertEdgeTriInput input);
    void isct(IsctVertTriTriTriInput input);
    void isctInterpolate(const VertData &v0, const VertData &v1);
    */
};

struct IsctTriangleData {
    // specify how to compute new data for a triangle in the event
    // that it is merged with another triangle (merge)
    // split into two triangles (split)
    // or that the triangle is moved (move)
    /*
    void subdivide(SubdivideTriInput input);
    */
};

// in order to perform remeshing, VertData and TriData must support
struct RemeshVertexData {
    bool manifold; // whether this point is manifold.
    // useful for modifying interpolation behavior
    // specify how to compute new data for a vertex in the event of
    // either an edge collapse (via merge) or edge split (via interpolate)
    /*
    void merge(const VertData &v0, const VertData &v1);
    void interpolate(const VertData &v0, const VertData &v1);
    */
};

struct RemeshTriangleData {
    // specify how to compute new data for a triangle in the event
    // that it is merged with another triangle (merge)
    // split into two triangles (split)
    // or that the triangle is moved (move)
    /*
    void merge(const TriData &t0, const TriData &t1);
    static void split(TriData &t0, TriData &t1, const TriData &t_orig);
    void move(const TriData &t_old);
    */
};

struct RemeshOptions
{
    double maxEdgeLength;
    double minEdgeLength;
    double minAngle;
    double maxAngle;

    RemeshOptions() :
            maxEdgeLength(1.0),
            minEdgeLength(0.3),
            minAngle(5.0),
            maxAngle(170.0)
    {}
};

// only for internal use, please do not use as client
struct TopoVert;
struct TopoEdge;
struct TopoTri;

typedef TopoVert* Vptr;
typedef TopoEdge* Eptr;
typedef TopoTri*  Tptr;
//using Vptr = TopoVert*;
//using Eptr = TopoEdge*;
//using Tptr = TopoTri*;
// end internal items

template<class VertData, class TriData>
class Mesh
{
public:
    Mesh();
    Mesh(Mesh &&src);
    Mesh(const RawMesh<VertData,TriData> &raw);
    Mesh(RawMesh<VertData,TriData> &&raw);
    virtual ~Mesh();

    void operator=(Mesh &&src);

    // validity check:
    //  - all numbers are well-defined and finite
    //  - all triangle vertex indices are in the right range
    bool valid() const;

    RawMesh<VertData,TriData> raw() const;

    // Parallel, copy-free export: fv(i, vert), ft(i, a, b, c)
    template<class FV, class FT>
    void export_parallel(FV fv, FT ft) const {
        cork_par::for_each_idx(verts.size(), 8192, [&](size_t i) { fv(i, verts[i]); });
        cork_par::for_each_idx(tris.size(), 8192, [&](size_t i) {
            ft(i, tris[i].a, tris[i].b, tris[i].c);
        });
    }

    inline int numVerts() const { return verts.size(); }
    inline int numTris() const { return tris.size(); }

    inline void for_verts(std::function<void(VertData &)> func);
    inline void for_tris(std::function<void(TriData &, VertData &, VertData &, VertData &)> func);
    inline void for_edges(std::function<void(VertData &, VertData &)> start,
                          std::function<void(TriData &t,VertData &, VertData &, VertData &)> each_tri);

    // form the disjoint union of two meshes
    void disjointUnion(const Mesh &cp);

    struct Isct {
        Ray3d   ray;
        bool    exists;

        uint    tri_id;
        Vec3d   isct;
        Vec3d   bary;
    };

    Isct pick(Ray3d ray);
    inline void accessIsct(const Isct &isct,
                           std::function<void(TriData &,
                                              VertData &, VertData &, VertData &)> func);

    // checks if the mesh is closed
    bool isClosed();

public: // REMESHING module
    // REQUIRES:
    //  - MinimalData
    //  - RemeshData
    void remesh();
    RemeshOptions remesh_options;

public: // ISCT (intersections) module
    void resolveIntersections(); // makes all intersections explicit
    bool isSelfIntersecting(); // is the mesh self-intersecting?

    void testingComputeStaticIsctPoints(std::vector<Vec3d> *points);
    void testingComputeStaticIsct(std::vector<Vec3d> *points,
                                  std::vector< std::pair<Vec3d,Vec3d> > *edges);

public: // OUTER HULL module (single mesh; call after resolveIntersections)
    // keeps only faces whose one side has generalized winding number 0,
    // flipping faces whose *back* side is exterior.  See mesh.hull.tpp.
    void outerHull(int raysPerPatch = 5, cork_hull::HullStats *stats = nullptr,
                   double leftoverAreaFrac = 0.0, bool exact = false);
    // Same patch hull as outerHull, but Manifold-style axis-aligned winding
    // (+Z, then +X/+Y).  Opt-in; default path stays random-ray outerHull().
    void outerHullExact(cork_hull::HullStats *stats = nullptr,
                        double leftoverAreaFrac = 0.0);
    // True if a busy vertex has an edge used 8+ times (stacked duplicate
    // faces).  Cheap; no mesh copy.  Used to decide pre-SI dup drop.
    bool hasStackedDuplicateFaces() const;
    // Sheet-pile soups: LBVH SI instead of packed-grid 336M cell pairs.
    bool preferLbvhIsct = false;

public: // BOOLean operation module
    // all of the form
    //      this = this OP rhs
    void boolUnion(Mesh &rhs);
    void boolDiff(Mesh &rhs);
    void boolIsct(Mesh &rhs);
    void boolXor(Mesh &rhs);

private:    // Internal Formats
    struct Tri {
        TriData data;
        union {
            struct {
                uint a, b, c; // vertex ids
            };
            uint v[3];
        };

        inline Tri() {}
    };

    inline void merge_tris(uint tid_result, uint tid0, uint tid1);
    inline void split_tris(uint t0ref, uint t1ref, uint t_orig_ref);
    inline void move_tri(Tri &t_new, Tri &t_old);
    inline void subdivide_tri(uint t_piece_ref, uint t_parent_ref);

private:    // DATA
    std::vector<Tri>        tris;
    std::vector<VertData>   verts;

private:    // caches
    struct NeighborEntry {
        uint vid;
        ShortVec<uint, 2> tids;
        inline NeighborEntry() {}
        inline NeighborEntry(uint vid_) : vid(vid_) {}
    };

    struct NeighborCache {
        std::vector< ShortVec<NeighborEntry, 8> > skeleton;
        inline NeighborEntry& operator()(uint i, uint j) {
            uint N = skeleton[i].size();
            for(uint k = 0; k < N; k++) {
                if(skeleton[i][k].vid == j)
                    return skeleton[i][k];
            }
            skeleton[i].push_back(NeighborEntry(j));
            return skeleton[i][N];
        }
    };

    NeighborCache createNeighborCache();

    // parallel to vertex array
    std::vector<uint> getComponentIds();

    // like the neighbor cache, but more customizable
    template<class Edata>
    struct EGraphEntry {
        uint                vid;
        ShortVec<uint, 2>   tids;
        Edata               data;
        inline EGraphEntry() {}
        inline EGraphEntry(uint vid_) : vid(vid_) {}
    };
    template<class Edata>
    struct EGraphCache {
        std::vector< ShortVec<EGraphEntry<Edata>, 8> > skeleton;
        inline EGraphEntry<Edata> & operator()(uint i, uint j) {
            uint N = skeleton[i].size();
            for(uint k = 0; k < N; k++) {
                if(skeleton[i][k].vid == j)
                    return skeleton[i][k];
            }
            skeleton[i].push_back(EGraphEntry<Edata>(j));
            return skeleton[i][N];
        }
        inline void for_each(std::function<void(
                uint i, uint j, EGraphEntry<Edata> &entry
        )> action
        ) {
            for(uint i=0; i<skeleton.size(); i++) {
                for(auto &entry : skeleton[i]) {
                    action(i, entry.vid, entry);
                }
            }
        }
    };
    template<class Edata>
    EGraphCache<Edata> createEGraphCache();


private:    // TopoCache Support
    struct TopoCache;
private:    // Isct Support
    class  IsctProblem; // implements intersection functionality
    class TriangleProblem; // support type for IsctProblem
    typedef TriangleProblem* Tprob;
    //using Tprob = TriangleProblem*;
private:    // Bool Support
    class BoolProblem;

private:    // Remeshing Support
    struct RemeshScratchpad;

    Eptr allocateRemeshEdge(RemeshScratchpad &);
    void deallocateRemeshEdge(RemeshScratchpad &, Eptr);

    void edgeSplit(RemeshScratchpad &,
                   Eptr e_split);
    void edgeCollapse(RemeshScratchpad &,
                      Eptr e_collapse,
                      bool collapsing_tetrahedra_disappear);

    // Need edge scoring routines...
    void scoreAndEnqueue(std::set< std::pair<double, Eptr> > &queue, Eptr edge);
    void dequeue(std::set< std::pair<double, Eptr> > &queue, Eptr edge);
    double computeEdgeScore(Eptr edge);

    // support functions
    void populateTriFromTopoTri(Tptr t);
    // calls the first function once, then the second once for each triangle
    inline void edgeNeighborhood(
            Eptr edge,
            std::function<void(VertData &v0, VertData &v1)> once,
            std::function<void(VertData &v0, VertData &v1,
                               VertData &vopp, TriData &t)> each_tri
    );
};

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::for_verts(std::function<void(VertData &v)> func) {
    for(auto &v : verts)
        func(v);
}

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::for_tris(std::function<void(TriData &, VertData &, VertData &, VertData &)> func) {
    for(auto &tri : tris) {
        auto &a = verts[tri.a];
        auto &b = verts[tri.b];
        auto &c = verts[tri.c];
        func(tri.data, a, b, c);
    }
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::for_edges(std::function<void(VertData &, VertData &)> start,
                                       std::function<void(TriData &t, VertData &, VertData &, VertData &)> each_tri
) {
    NeighborCache cache = createNeighborCache();
    for(uint i=0; i<cache.skeleton.size(); i++) {
        for(auto &entry : cache.skeleton[i]) {
            uint j = entry.vid;
            start(verts[i], verts[j]);
            for(uint tid : entry.tids) {
                Tri &tri = tris[tid];
                each_tri(tri.data, verts[tri.a], verts[tri.b], verts[tri.c]);
            }
        }
    }
}

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::accessIsct(const Isct &isct,
                                               std::function<void(TriData &, VertData &, VertData &, VertData &)> func) {
    Tri &tri = tris[isct.tri_id];
    auto &a = verts[tri.a];
    auto &b = verts[tri.b];
    auto &c = verts[tri.c];
    func(tri.data, a, b, c);
}


/*
 * Implementation of mesh
 */
// constructors
template<class VertData, class TriData>
Mesh<VertData,TriData>::Mesh() {}
template<class VertData, class TriData>
Mesh<VertData,TriData>::Mesh(Mesh &&cp)
        : tris(cp.tris), verts(cp.verts)
{}
template<class VertData, class TriData>
Mesh<VertData,TriData>::Mesh(const RawMesh<VertData,TriData> &raw) :
        tris(raw.triangles.size()), verts(raw.vertices)
{
    const size_t nt = raw.triangles.size();
    cork_par::for_each_idx(nt, 8192, [&](size_t i) {
        tris[i].data = raw.triangles[i];
        tris[i].a = raw.triangles[i].a;
        tris[i].b = raw.triangles[i].b;
        tris[i].c = raw.triangles[i].c;
    });
}
template<class VertData, class TriData>
Mesh<VertData,TriData>::Mesh(RawMesh<VertData,TriData> &&raw) :
        tris(raw.triangles.size()), verts(std::move(raw.vertices))
{
    const size_t nt = raw.triangles.size();
    cork_par::for_each_idx(nt, 8192, [&](size_t i) {
        tris[i].data = raw.triangles[i];
        tris[i].a = raw.triangles[i].a;
        tris[i].b = raw.triangles[i].b;
        tris[i].c = raw.triangles[i].c;
    });
}
template<class VertData, class TriData>
Mesh<VertData,TriData>::~Mesh()
{

}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::operator=(Mesh &&src)
{
    tris = src.tris;
    verts = src.verts;
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::valid() const
{
    for(uint i=0; i<verts.size(); i++) {
        if(!std::isfinite(verts[i].pos.x) ||
           !std::isfinite(verts[i].pos.y) ||
           !std::isfinite(verts[i].pos.z)) {
            std::ostringstream message;
            message << "vertex #" << i << " has non-finite coordinates: "
                    << verts[i].pos;
            CORK_ERROR(message.str());
            return false;
        }
    }

    for(uint i=0; i<tris.size(); i++) {
        if(tris[i].a >= verts.size() ||
           tris[i].b >= verts.size() ||
           tris[i].c >= verts.size()) {
            std::ostringstream message;
            message << "triangle #" << i << " should have indices in "
                    << "the range 0 to " << (verts.size()-1)
                    << ", but it has invalid indices: "
                    << tris[i].a << ", " << tris[i].b << ", " << tris[i].c;
            CORK_ERROR(message.str());
            return false;
        }
    }

    return true;
}

template<class VertData, class TriData>
RawMesh<VertData,TriData> Mesh<VertData,TriData>::raw() const
{
    RawMesh<VertData,TriData> result;
    result.vertices = verts;
    result.triangles.resize(tris.size());
    for(uint i=0; i<tris.size(); i++) {
        result.triangles[i]   = tris[i].data;
        result.triangles[i].a = tris[i].a;
        result.triangles[i].b = tris[i].b;
        result.triangles[i].c = tris[i].c;
    }
    return result;
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::disjointUnion(const Mesh &cp)
{
    uint oldVsize = verts.size();
    uint oldTsize = tris.size();
    uint cpVsize  = cp.verts.size();
    uint cpTsize  = cp.tris.size();
    uint newVsize = oldVsize + cpVsize;
    uint newTsize = oldTsize + cpTsize;

    std::vector<int> v_remap(cpVsize); // oh this is obvious...
    verts.resize(newVsize);
    tris.resize(newTsize);

    for(uint i=0; i<cpVsize; i++)
        verts[oldVsize + i] = cp.verts[i];

    for(uint i=0; i<cpTsize; i++) {
        auto &tri = tris[oldTsize + i];
        tri = cp.tris[i];
        tri.a += oldVsize;
        tri.b += oldVsize;
        tri.c += oldVsize;
    }
}

// Picking.
// Dumb Implementation just passes over all triangles w/o any precomputed
// acceleration structure
template<class VertData, class TriData>
typename Mesh<VertData,TriData>::Isct
Mesh<VertData,TriData>::pick(Ray3d ray)
{
    Isct result;
    result.ray = ray;
    result.exists = false;

    double mint = DBL_MAX;

    // pass all triangles over ray
    for(uint i=0; i<tris.size(); i++) {
        const Tri  &tri = tris[i];

        uint   a = tri.a;
        uint   b = tri.b;
        uint   c = tri.c;
        Vec3d va = verts[a].pos;
        Vec3d vb = verts[b].pos;
        Vec3d vc = verts[c].pos;
        // normalize vertex order (to prevent leaks)
        if(a > b) { std::swap(a, b); std::swap(va, vb); }
        if(b > c) { std::swap(b, c); std::swap(vb, vc); }
        if(a > b) { std::swap(a, b); std::swap(va, vb); }

        double t;
        Vec3d  bary;
        if(isct_ray_triangle(ray, va, vb, vc, &t, &bary)) {
            if(t > 0 && t < mint) {
                result.exists = true;
                mint = t;
                result.tri_id = i;
                result.isct = ray.p + t * ray.r;
                result.bary = bary;
            }
        }
    }

    return result;
}




template<class VertData, class TriData>
bool Mesh<VertData,TriData>::isClosed()
{
    EGraphCache<int> chains = createEGraphCache<int>();
    chains.for_each([&](uint i, uint j, EGraphEntry<int> &entry) {
        entry.data = 0;
    });
    // count up how many times each edge is encountered in one
    // orientation vs. the other
    for(Tri &tri : tris) {
        chains(tri.a, tri.b).data ++;
        chains(tri.b, tri.a).data --;

        chains(tri.b, tri.c).data ++;
        chains(tri.c, tri.b).data --;

        chains(tri.c, tri.a).data ++;
        chains(tri.a, tri.c).data --;
    }
    // now go through and see if any of these are non-zero
    bool closed = true;
    chains.for_each([&](uint i, uint j, EGraphEntry<int> &entry) {
        if(entry.data != 0)
            closed = false;
    });
    return closed;
}




static inline
bool contains(const ShortVec<uint, 8> &list, uint item)
{
    for(uint k : list)
        if(k == item)
            return true;
    return false;
}

template<class VertData, class TriData>
typename Mesh<VertData,TriData>::NeighborCache
Mesh<VertData,TriData>::createNeighborCache()
{
    NeighborCache result;
    result.skeleton.resize(verts.size());

    for(uint tid = 0; tid < tris.size(); tid++) {
        const Tri &tri = tris[tid];

        result(tri.a, tri.b).tids.push_back(tid);
        result(tri.b, tri.a).tids.push_back(tid);

        result(tri.a, tri.c).tids.push_back(tid);
        result(tri.c, tri.a).tids.push_back(tid);

        result(tri.b, tri.c).tids.push_back(tid);
        result(tri.c, tri.b).tids.push_back(tid);
    }

    return result;
}

// This function signature is an amazing disaster...
#ifdef _WIN32
template<class VertData, class TriData>
template<class Edata>
typename Mesh<VertData,TriData>::EGraphCache<Edata>
#else
template<class VertData, class TriData>
template<class Edata>
typename Mesh<VertData,TriData>::template EGraphCache<Edata>
#endif
Mesh<VertData,TriData>::createEGraphCache()
{
    EGraphCache<Edata> result;
    result.skeleton.resize(verts.size());

    for(uint tid = 0; tid < tris.size(); tid++) {
        const Tri &tri = tris[tid];

        result(tri.a, tri.b).tids.push_back(tid);
        result(tri.b, tri.a).tids.push_back(tid);

        result(tri.a, tri.c).tids.push_back(tid);
        result(tri.c, tri.a).tids.push_back(tid);

        result(tri.b, tri.c).tids.push_back(tid);
        result(tri.c, tri.b).tids.push_back(tid);
    }

    return result;
}


template<class VertData, class TriData>
std::vector<uint> Mesh<VertData,TriData>::getComponentIds()
{
    UnionFind uf(verts.size());
    for(const Tri &tri : tris) {
        uf.unionIds(tri.a, tri.b);
        uf.unionIds(tri.a, tri.c);
    }

    return uf.dump();
}

#include <cork/util/iterPool.h>

/*
 *  Allows for topological algorithms to manipulate
 *  a more familiar pointer data structure based on a simplicial complex.
 *  This structure can be regenerated from the more basic
 *  vertex/triangle arrays using
 *      createTopoCache()
 *  Once manipulations have been satisfactorily performed,
 *  the underlying vertex/triangle arrays can be cleaned up for
 *  further use by topologically insensitive algorithms by
 *      commitTopoCache()
 */

#define INVALID_ID uint(-1)

struct TopoVert {
    uint                    ref;        // index to actual data
    void*                   data;       // algorithm specific handle

    ShortVec<Tptr, 8>       tris;       // triangles this vertex is incident on
    ShortVec<Eptr, 8>       edges;      // edges this vertex is incident on
};

struct TopoEdge {
    void*                   data;       // algorithm specific handle

    Vptr                    verts[2];   // endpoint vertices
    ShortVec<Tptr, 2>       tris;       // incident triangles
};

struct TopoTri {
    uint                    ref;        // index to actual data
    void*                   data;       // algorithm specific handle

    Vptr                    verts[3];   // vertices of this triangle
    Eptr                    edges[3];   // edges of this triangle
                                        // opposite to the given vertex
};


template<class VertData, class TriData>
struct Mesh<VertData, TriData>::TopoCache {
    IterPool<TopoVert>    verts;
    IterPool<TopoEdge>    edges;
    IterPool<TopoTri>     tris;

    Mesh *mesh;
    TopoCache(Mesh *owner);
    // vertEdges=false skips vertex->edge lists; vertTris=false skips
    // vertex->triangle lists (IsctProblem never reads either after init).
    TopoCache(Mesh *owner, bool vertEdges);
    TopoCache(Mesh *owner, bool vertEdges, bool vertTris);
    virtual ~TopoCache() {}

    // until commit() is called, the Mesh::verts and Mesh::tris
    // arrays will still contain garbage entries
    void commit();

    bool isValid();
    void print();

    // helpers to create bits and pieces
    inline Vptr newVert();
    inline Eptr newEdge();
    inline Tptr newTri();

    // helpers to release bits and pieces
    inline void freeVert(Vptr);
    inline void freeEdge(Eptr);
    inline void freeTri(Tptr);

    // helper to delete geometry in a structured way
    inline void deleteTri(Tptr);

    // helper to flip triangle orientation
    inline void flipTri(Tptr);

    // Contiguous first-wave allocations (valid until individual free()).
    IterPool<TopoVert>::Bulk vbulk;
    IterPool<TopoEdge>::Bulk ebulk;
    IterPool<TopoTri>::Bulk  tbulk;

private:
    void init(bool vertEdges = true, bool vertTris = true);
};


template<class VertData, class TriData> inline
Vptr Mesh<VertData, TriData>::TopoCache::newVert()
{
    uint        ref         = mesh->verts.size();
                mesh->verts.push_back(VertData());
    Vptr        v           = verts.alloc(); // cache.verts
                v->ref      = ref;
                return v;
}
template<class VertData, class TriData> inline
Eptr Mesh<VertData, TriData>::TopoCache::newEdge()
{
    Eptr        e           = edges.alloc(); // cache.edges
                return e;
}
template<class VertData, class TriData> inline
Tptr Mesh<VertData, TriData>::TopoCache::newTri()
{
    uint        ref         = mesh->tris.size();
                mesh->tris.push_back(Tri());
    Tptr        t           = tris.alloc(); // cache.tris
                t->ref      = ref;
                return t;
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::freeVert(Vptr v)
{
    verts.free(v);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::freeEdge(Eptr e)
{
    edges.free(e);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::freeTri(Tptr t)
{
    tris.free(t);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::deleteTri(Tptr tri)
{
    // first, unhook the triangle from its faces
    for(uint k=0; k<3; k++) {
        Vptr            v                   = tri->verts[k];
                        v->tris.erase(tri);
        Eptr            e                   = tri->edges[k];
                        e->tris.erase(tri);
    }
    // now, let's check for any edges which no longer border triangles
    for(uint k=0; k<3; k++) {
        Eptr            e                   = tri->edges[k];
        if(e->tris.size() == 0) {
            // delete edge
            // unhook from vertices
            Vptr        v0                  = e->verts[0];
                        v0->edges.erase(e);
            Vptr        v1                  = e->verts[1];
                        v1->edges.erase(e);
            freeEdge(e);
        }
    }
    // now, let's check for any vertices which no longer border triangles
    for(uint k=0; k<3; k++) {

        Vptr v = tri->verts[k];

        if(v->tris.size() == 0) {
            freeVert(v);
        }

    }

    // finally, release the triangle
    freeTri(tri);
}

template<class VertData, class TriData> inline
void Mesh<VertData, TriData>::TopoCache::flipTri(Tptr t)
{
    std::swap(t->verts[0], t->verts[1]);
    std::swap(t->edges[0], t->edges[1]);
    std::swap(mesh->tris[t->ref].v[0], mesh->tris[t->ref].v[1]);
}

template<class VertData, class TriData>
Mesh<VertData, TriData>::TopoCache::TopoCache(Mesh *owner) : mesh(owner)
{
    init(true);
}
template<class VertData, class TriData>
Mesh<VertData, TriData>::TopoCache::TopoCache(Mesh *owner, bool vertEdges) : mesh(owner)
{
    init(vertEdges, true);
}
template<class VertData, class TriData>
Mesh<VertData, TriData>::TopoCache::TopoCache(Mesh *owner, bool vertEdges, bool vertTris) : mesh(owner)
{
    init(vertEdges, vertTris);
}


// support structure for cache construction
struct TopoEdgePrototype {
    uint vid;
    ShortVec<Tptr, 2> tris;
    TopoEdgePrototype() {}
    TopoEdgePrototype(uint v) : vid(v) {}
};

inline TopoEdgePrototype& getTopoEdgePrototype(uint a, uint b,
                                               std::vector< ShortVec<TopoEdgePrototype, 8> > &prototypes) {
    uint N = prototypes[a].size();

    for(uint i=0; i<N; i++) {
        if(prototypes[a][i].vid == b)
            return prototypes[a][i];
    }
    prototypes[a].push_back(TopoEdgePrototype(b));
    return prototypes[a][N];
}

// Parallel topology construction.
//  * vertices and triangles are bulk allocated and initialised in parallel
//  * edges are found by sorting the 3 (min,max) vertex-id keys of every
//    triangle; each run of equal keys is one edge
//  * vertex incidence lists are filled from CSR buckets, one vertex per task
// The resulting structure is identical to the original serial build
// (same verts[0]/verts[1] ordering, same tri->edges[k] correspondence, same
//  ascending-triangle order in incidence lists).
template<class VertData, class TriData>
void Mesh<VertData, TriData>::TopoCache::init(bool vertEdges, bool vertTris)
{
    CORK_PROF("    TopoCache::init");
    const size_t nv = mesh->verts.size();
    const size_t nt = mesh->tris.size();

    // ---- vertices ----
    vbulk = verts.alloc_bulk((uint)nv);
    auto vb = vbulk;
    cork_par::for_each_idx(nv, 8192, [&](size_t i) {
        Vptr v  = vb[i];
        v->ref  = (uint)i;
        v->data = nullptr;
    });

    // ---- triangles ----
    tbulk = tris.alloc_bulk((uint)nt);
    auto tb = tbulk;
    cork_par::for_each_idx(nt, 8192, [&](size_t i) {
        Tptr t  = tb[i];
        t->ref  = (uint)i;
        t->data = nullptr;
        const Tri &rt = mesh->tris[i];
        for(uint k=0; k<3; k++) t->verts[k] = vb[rt.v[k]];
    });

    // ---- vertex -> triangle (skipped by IsctProblem: never read) ----
    if (vertTris) {
        std::vector<unsigned> off, ord;
        cork_par::build_csr(nv, 3*nt,
            [&](size_t c) { return (size_t)mesh->tris[c/3].v[c%3]; }, off, ord);
        cork_par::for_each_idx(nv, 4096, [&](size_t v) {
            Vptr vp = vb[v];
            for(unsigned k = off[v]; k < off[v+1]; ++k)
                vp->tris.push_back(tb[ord[k]/3]);
        });
    }

    // ---- edges ----
    struct EKey { uint64_t key; uint32_t tri; uint32_t k; };
    std::vector<EKey> ek(3*nt);
    cork_par::for_each_idx(nt, 8192, [&](size_t i) {
        const Tri &rt = mesh->tris[i];
        for(uint k=0; k<3; k++) {
            uint a = rt.v[(k+1)%3], b = rt.v[(k+2)%3];
            uint lo = std::min(a,b), hi = std::max(a,b);
            EKey &e = ek[3*i+k];
            e.key = ((uint64_t)lo << 32) | (uint64_t)hi;
            e.tri = (uint32_t)i;
            e.k   = k;
        }
    });
    cork_par::sort(ek.begin(), ek.end(), [](const EKey &a, const EKey &b) {
        if(a.key != b.key) return a.key < b.key;
        if(a.tri != b.tri) return a.tri < b.tri;
        return a.k < b.k;
    });

    std::vector<uint32_t> runStart;
    runStart.reserve(3*nt/2 + 1);
    for(size_t i=0; i<ek.size(); i++)
        if(i == 0 || ek[i].key != ek[i-1].key) runStart.push_back((uint32_t)i);
    runStart.push_back((uint32_t)ek.size());
    const size_t ne = runStart.size() - 1;

    ebulk = edges.alloc_bulk((uint)ne);
    auto eb = ebulk;
    cork_par::for_each_idx(ne, 4096, [&](size_t r) {
        Eptr e   = eb[r];
        e->data  = nullptr;
        uint64_t key = ek[runStart[r]].key;
        e->verts[0] = vb[(size_t)(key >> 32)];
        e->verts[1] = vb[(size_t)(key & 0xffffffffu)];
        for(uint32_t idx = runStart[r]; idx < runStart[r+1]; ++idx) {
            Tptr t = tb[ek[idx].tri];
            e->tris.push_back(t);
            t->edges[ek[idx].k] = e;
        }
    });

    // ---- vertex -> edge ----
    if(vertEdges) {
        std::vector<unsigned> off, ord;
        cork_par::build_csr(nv, 2*ne, [&](size_t j) {
            uint64_t key = ek[runStart[j/2]].key;
            return (size_t)((j & 1) ? (key & 0xffffffffu) : (key >> 32));
        }, off, ord);
        cork_par::for_each_idx(nv, 4096, [&](size_t v) {
            Vptr vp = vb[v];
            for(unsigned k = off[v]; k < off[v+1]; ++k)
                vp->edges.push_back(eb[ord[k]/2]);
        });
    }

    //ENSURE(isValid());
    //print();
}




template<class VertData, class TriData>
void Mesh<VertData, TriData>::TopoCache::commit()
{
    //ENSURE(isValid());
    std::vector<Vptr> vv;   verts.collect(vv);
    std::vector<Tptr> tv;   tris.collect(tv);
    const size_t NV = mesh->verts.size();
    const size_t NT = mesh->tris.size();

    // record which vertices are live
    std::vector<uint8_t> live_verts(NV, 0);
    cork_par::for_each_idx(vv.size(), 8192, [&](size_t i) {
        live_verts[vv[i]->ref] = 1;
    });

    // record which triangles are live, and record connectivity
    std::vector<uint8_t> live_tris(NT, 0);
    cork_par::for_each_idx(tv.size(), 8192, [&](size_t i) {
        Tptr tri = tv[i];
        live_tris[tri->ref] = 1;
        for(uint k=0; k<3; k++)
            mesh->tris[tri->ref].v[k] = tri->verts[k]->ref;
    });

    // compact the vertices and build a remapping function
    std::vector<uint> vmap(NV);
    std::vector<uint32_t> vps(NV);
    cork_par::for_each_idx(NV, 8192, [&](size_t i) { vps[i] = live_verts[i]; });
    cork_par::prefix_sum_inplace(vps.data(), NV);
    uint write = NV ? vps[NV - 1] : 0;
    cork_par::for_each_idx(NV, 8192, [&](size_t read) {
        if (live_verts[read]) vmap[read] = vps[read] - 1;
        else                  vmap[read] = INVALID_ID;
    });
    {
        std::vector<VertData> nverts(write);
        cork_par::for_each_idx(NV, 8192, [&](size_t read) {
            if(live_verts[read]) nverts[vmap[read]] = mesh->verts[read];
        });
        mesh->verts.swap(nverts);
    }

    // rewrite the vertex reference ids
    cork_par::for_each_idx(vv.size(), 8192, [&](size_t i) {
        vv[i]->ref = vmap[vv[i]->ref];
    });

    std::vector<uint> tmap(NT);
    std::vector<uint32_t> tps(NT);
    cork_par::for_each_idx(NT, 8192, [&](size_t i) { tps[i] = live_tris[i]; });
    cork_par::prefix_sum_inplace(tps.data(), NT);
    write = NT ? tps[NT - 1] : 0;
    cork_par::for_each_idx(NT, 8192, [&](size_t read) {
        if (live_tris[read]) tmap[read] = tps[read] - 1;
        else                 tmap[read] = INVALID_ID;
    });
    {
        std::vector<Tri> ntris(write);
        cork_par::for_each_idx(NT, 8192, [&](size_t read) {
            if(!live_tris[read]) return;
            Tri &dst = ntris[tmap[read]];
            dst = mesh->tris[read];
            for(uint k=0; k<3; k++)
                dst.v[k] = vmap[dst.v[k]];
        });
        mesh->tris.swap(ntris);
    }

    // rewrite the triangle reference ids
    cork_par::for_each_idx(tv.size(), 8192, [&](size_t i) {
        tv[i]->ref = tmap[tv[i]->ref];
    });
}



// support functions for validity check
template<class T, class Container> inline
bool count(const Container &contain, const T &val) {
    uint c=0;
    for(const T &t : contain)
        if(t == val)    c++;
    return c;
}
template<class T> inline
bool count2(const T arr[], const T &val) {
    return ((arr[0] == val)? 1 : 0) + ((arr[1] == val)? 1 : 0);
}
template<class T> inline
bool count3(const T arr[], const T &val) {
    return ((arr[0] == val)? 1 : 0) + ((arr[1] == val)? 1 : 0)
                                    + ((arr[2] == val)? 1 : 0);
}

template<class VertData, class TriData>
bool Mesh<VertData, TriData>::TopoCache::isValid()
{
    //print();
    std::set<Vptr> vaddr;
    std::set<Eptr> eaddr;
    std::set<Tptr> taddr;
    verts.for_each([&vaddr](Vptr v) { vaddr.insert(v); });
    edges.for_each([&eaddr](Eptr e) { eaddr.insert(e); });
    tris.for_each( [&taddr](Tptr t) { taddr.insert(t); });

    // check verts
    verts.for_each([&](Vptr v) {
        ENSURE(v->ref < mesh->verts.size());
        // make sure each edge pointer goes somewhere and that
        // the pointed-to site also points back correctly
        for(Eptr e : v->edges) {
            ENSURE(eaddr.count(e) > 0); // pointer is good
            ENSURE(count2(e->verts, v) == 1); // back-pointer is good
        }
        for(Tptr t : v->tris) {
            ENSURE(taddr.count(t) > 0);
            ENSURE(count3(t->verts, v) == 1);
        }
    });

    // check edges
    edges.for_each([&](Eptr e) {
        // check for non-degeneracy
        ENSURE(e->verts[0] != e->verts[1]);
        for(uint k=0; k<2; k++) {
            Vptr v = e->verts[k];
            ENSURE(vaddr.count(v) > 0);
            ENSURE(count(v->edges, e) == 1);
        }
        for(Tptr t : e->tris) {
            ENSURE(taddr.count(t) > 0);
            ENSURE(count3(t->edges, e) == 1);
        }
    });

    // check triangles
    tris.for_each([&](Tptr t) {
        // check for non-degeneracy
        ENSURE(t->verts[0] != t->verts[1] && t->verts[1] != t->verts[2]
                                          && t->verts[0] != t->verts[2]);
        for(uint k=0; k<3; k++) {
            Vptr v = t->verts[k];
            ENSURE(vaddr.count(v) > 0);
            ENSURE(count(v->tris, t) == 1);

            Eptr e = t->edges[k];
            ENSURE(eaddr.count(e) == 1);
            ENSURE(count(e->tris, t) == 1);

            // also need to ensure that the edges are opposite the
            // vertices as expected
            Vptr v0 = e->verts[0];
            Vptr v1 = e->verts[1];
            ENSURE((v0 == t->verts[(k+1)%3] && v1 == t->verts[(k+2)%3])
                || (v0 == t->verts[(k+2)%3] && v1 == t->verts[(k+1)%3]));
        }
    });

    return true;
}


CORK_EXPORT std::ostream& operator<<(std::ostream &out, const TopoVert& vert);
CORK_EXPORT std::ostream& operator<<(std::ostream &out, const TopoEdge& edge);
CORK_EXPORT std::ostream& operator<<(std::ostream &out, const TopoTri& tri);

template<class VertData, class TriData>
void Mesh<VertData, TriData>::TopoCache::print()
{
    using std::cout;
    using std::endl;

    cout << "dumping remeshing cache for debug..." << endl;
    cout << "TRIS" << endl;
    int tri_count = 0;
    tris.for_each([&](Tptr t) {
        cout << " " << t << ": " << *t << endl;
        tri_count++;
    });
    cout << "There were " << tri_count << " TRIS" << endl;
    cout << "EDGES" << endl;
    int edge_count = 0;
    edges.for_each([&](Eptr e) {
        cout << " " << e << ": " << endl;
        cout << "  v " << e->verts[0] << "; "
                       << e->verts[1] << endl;
        cout << "  t (" << e->tris.size() << ")" << endl;
        for(Tptr t : e->tris)
        cout << "    " << t << endl;
        edge_count++;
    });
    cout << "There were " << edge_count << " EDGES" << endl;
    cout << "VERTS" << endl;
    int vert_count = 0;
    verts.for_each([&](Vptr v) {
        cout << " " << v << ": ref(" << v->ref << ")" << endl;
        cout << "  e (" << v->edges.size() << ")" << endl;
        for(Eptr e : v->edges)
        cout << "    " << e << endl;
        cout << "  t (" << v->tris.size() << ")" << endl;
        for(Tptr t : v->tris)
        cout << "    " << t << endl;
        vert_count++;
    });
    cout << "There were " << vert_count << " VERTS" << endl;
}


// Implementation

#include "mesh.remesh.tpp"
#include "mesh.isct.tpp"
#include "mesh.bool.tpp"
#include "mesh.hull.tpp"
#include "mesh.hull_exact.tpp"


#endif

