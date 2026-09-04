// +-------------------------------------------------------------------------
// | mesh.isct.tpp
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
#pragma once

#include <utility>
#include <vector>
#include <array>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <cork/accel/bbox_avx.h>
#include <cork/util/parallel.h>
#include <cork/util/profile.h>

struct GenericVertType;
    struct IsctVertType;
    struct OrigVertType;
struct GenericEdgeType;
    struct IsctEdgeType;
    struct OrigEdgeType;
    struct SplitEdgeType;
struct GenericTriType;

struct GluePointMarker;

//using GVptr             = GenericVertType*;
//    using IVptr         = IsctVertType*;
//    using OVptr         = OrigVertType*;
//using GEptr             = GenericEdgeType*;
//    using IEptr         = IsctEdgeType*;
//    using OEptr         = OrigEdgeType*;
//    using SEptr         = SplitEdgeType*;
//using GTptr             = GenericTriType*;

//using GluePt = GluePointMarker*;

typedef GenericVertType*    GVptr;
typedef IsctVertType*           IVptr;
typedef OrigVertType*	        OVptr;
typedef GenericEdgeType*    GEptr;
typedef IsctEdgeType*           IEptr;
typedef OrigEdgeType*	        OEptr;
typedef SplitEdgeType*	        SEptr;
typedef GenericTriType*     GTptr;

typedef GluePointMarker*    GluePt;

struct GenericVertType
{
    virtual ~GenericVertType() {}
    Vptr                    concrete;
    Vec3d                   coord;
    
    bool                    boundary;
    uint                    idx; // temporary for triangulation marshalling
    
    ShortVec<GEptr,2>       edges;
};
struct IsctVertType : public GenericVertType
{
    GluePt                  glue_marker;
};
struct OrigVertType : public GenericVertType {};

struct GenericEdgeType
{
    virtual ~GenericEdgeType() {}
    Eptr                    concrete;
    
    bool                    boundary;
    uint                    idx; // temporary for triangulation marshalling
    
    GVptr                   ends[2];
    ShortVec<IVptr, 1>      interior;       
};
struct IsctEdgeType : public GenericEdgeType
{
public:
    // use to detect duplicate instances within a triangle
    Tptr                    other_tri_key;
};
struct OrigEdgeType : public GenericEdgeType {};
struct SplitEdgeType : public GenericEdgeType {};

struct GenericTriType
{
    Tptr                    concrete;
    
    GVptr                   verts[3];
};

struct GluePointMarker
{
    // list of all the vertices to be glued...
    ShortVec<IVptr, 3>      copies;
    bool                    split_type; // splits are introduced
                                        // manually, not via intersection
                                        // and therefore use only e pointer
    bool                    edge_tri_type; // true if edge-tri intersection
                                           // false if tri-tri-tri
    Eptr                    e;
    Tptr                    t[3];
};



template<uint LEN> inline
IEptr find_edge(ShortVec<IEptr,LEN> &vec, Tptr key)
{
    for(IEptr ie : vec) {
        if(ie->other_tri_key == key)
            return ie;
    }
    return nullptr;
}

inline Vptr commonVert(Tptr t0, Tptr t1)
{
    for(uint i=0; i<3; i++) {
      for(uint j=0; j<3; j++) {
        if(t0->verts[i] == t1->verts[j])
            return t0->verts[i];
      }
    }
    return nullptr;
}

inline bool hasCommonVert(Tptr t0, Tptr t1)
{
    return (t0->verts[0] == t1->verts[0] ||
            t0->verts[0] == t1->verts[1] ||
            t0->verts[0] == t1->verts[2] ||
            t0->verts[1] == t1->verts[0] ||
            t0->verts[1] == t1->verts[1] ||
            t0->verts[1] == t1->verts[2] ||
            t0->verts[2] == t1->verts[0] ||
            t0->verts[2] == t1->verts[1] ||
            t0->verts[2] == t1->verts[2]);
}

inline bool hasCommonVert(Eptr e, Tptr t)
{
    return (e->verts[0] == t->verts[0] ||
            e->verts[0] == t->verts[1] ||
            e->verts[0] == t->verts[2] ||
            e->verts[1] == t->verts[0] ||
            e->verts[1] == t->verts[1] ||
            e->verts[1] == t->verts[2]);
}

inline void disconnectGE(GEptr ge)
{
    ge->ends[0]->edges.erase(ge);
    ge->ends[1]->edges.erase(ge);
    for(IVptr iv : ge->interior)
        iv->edges.erase(ge);
}

// should deal with via pointers
template<class VertData, class TriData>
class Mesh<VertData, TriData>::TriangleProblem
{
public:
    TriangleProblem() {}
    ~TriangleProblem() {}
    
    inline void init(IsctProblem *iprob, Tptr t) {
        the_tri             = t;
        // extract original edges/verts
        for(uint k=0; k<3; k++)
            overts[k]       = iprob->newOrigVert(the_tri->verts[k]);
        for(uint k=0; k<3; k++) {
            oedges[k]       = iprob->newOrigEdge(the_tri->edges[k],
                                                 overts[(k+1)%3],
                                                 overts[(k+2)%3]);
        }
    }
    
private: // may actually not add edge, but instead just hook up endpoint
    inline void addEdge(
        IsctProblem *iprob, IVptr iv, Tptr tri_key
    ) {
        IEptr       ie              = find_edge(iedges, tri_key);
        if(ie) { // if the edge is already present
                    ie->ends[1]     = iv;
                    iv->edges.push_back(ie);
        } else { // if the edge is being added
                    ie              = iprob->newIsctEdge(iv, tri_key);
                    iedges.push_back(ie);
        }
    }
    void addBoundaryHelper(
        Eptr edge, IVptr iv
    ) {
                    iv->boundary    = true;
                    iverts.push_back(iv);
        // hook up point to boundary edge interior!
        for(uint k=0; k<3; k++) {
            OEptr   oe              = oedges[k];
            if(oe->concrete == edge) {
                    oe->interior.push_back(iv);
                    iv->edges.push_back(oe);
                    break;
            }
        }
    }
public:
    // specify reference glue point and edge piercing this triangle.
    IVptr addInteriorEndpoint(
        IsctProblem *iprob, Eptr edge, GluePt glue
    ) {
        IVptr       iv              = iprob->newIsctVert(edge, the_tri, glue);
                    iv->boundary    = false;
                    iverts.push_back(iv);
        for(Tptr tri_key : edge->tris) {
                    addEdge(iprob, iv, tri_key);
        }
        return iv;
    }
    // same, with the (exact) intersection coordinate already computed
    IVptr addInteriorEndpoint(
        IsctProblem *iprob, Eptr edge, GluePt glue, const Vec3d &coord
    ) {
        IVptr       iv              = iprob->newSplitIsctVert(coord, glue);
                    iv->boundary    = false;
                    iverts.push_back(iv);
        for(Tptr tri_key : edge->tris) {
                    addEdge(iprob, iv, tri_key);
        }
        return iv;
    }
    // specify the other triangle cutting this one, the edge cut,
    // and the resulting point of intersection
    void addBoundaryEndpoint(
        IsctProblem *iprob, Tptr tri_key, Eptr edge, IVptr iv
    ) {
                    iv              = iprob->copyIsctVert(iv);
                    addBoundaryHelper(edge, iv);
        // handle edge extending into interior
                    addEdge(iprob, iv, tri_key);
    }
    IVptr addBoundaryEndpoint(
        IsctProblem *iprob, Tptr tri_key, Eptr edge, Vec3d coord, GluePt glue
    ) {
        IVptr       iv              = iprob->newSplitIsctVert(coord, glue);
                    addBoundaryHelper(edge, iv);
        // handle edge extending into interior
                    addEdge(iprob, iv, tri_key);
        return iv;
    }
    // Should only happen for manually inserted split points on
    // edges, not for points computed via intersection...
    IVptr addBoundaryPointAlone(
        IsctProblem *iprob, Eptr edge, Vec3d coord, GluePt glue
    ) {
        IVptr       iv              = iprob->newSplitIsctVert(coord, glue);
                    addBoundaryHelper(edge, iv);
        return iv;
    }
    void addInteriorPoint(
        IsctProblem *iprob, Tptr t0, Tptr t1, GluePt glue
    ) {
        addInteriorPoint(iprob, t0, t1, glue,
                         iprob->computeCoords(the_tri, t0, t1));
    }
    void addInteriorPoint(
        IsctProblem *iprob, Tptr t0, Tptr t1, GluePt glue, const Vec3d &coord
    ) {
        IVptr       iv              = iprob->newSplitIsctVert(coord, glue);
                    iv->boundary    = false;
                    iverts.push_back(iv);
        // find the 2 interior edges
        for(IEptr ie : iedges) {
            if(ie->other_tri_key == t0 ||
               ie->other_tri_key == t1) {
                    ie->interior.push_back(iv);
                    iv->edges.push_back(ie);
            }
        }
    }
    
    // run after we've accumulated all the elements
    void consolidate(IsctProblem *iprob) {
        // identify all intersection edges missing endpoints
        // and check to see if we can assign an original vertex
        // as the appropriate endpoint.
        for(IEptr ie : iedges) {
            if(ie->ends[1] == nullptr) {
                // try to figure out which vertex must be the endpoint...
                Vptr vert = commonVert(the_tri, ie->other_tri_key);
                if(!vert) {
                    std::cout << "the  edge is "
                              << ie->ends[0] << ",  "
                              << ie->ends[1] << std::endl;
                    IVptr iv = dynamic_cast<IVptr>(ie->ends[0]);
                    std::cout << "   "
                              << iv->glue_marker->edge_tri_type
                              << std::endl;
                    std::cout << "the   tri is " << the_tri << ": "
                              << *the_tri << std::endl;
                    std::cout << "other tri is " << ie->other_tri_key << ": "
                              << *(ie->other_tri_key) << std::endl;
                    std::cout << "coordinates for triangles" << std::endl;
                    std::cout << "the tri" << std::endl;
                    for(uint k=0; k<3; k++)
                        std::cout << iprob->vPos(the_tri->verts[k])
                                  << std::endl;
                    for(uint k=0; k<3; k++)
                        std::cout << iprob->vPos(ie->other_tri_key->verts[k])
                                  << std::endl;
                    std::cout << "degen count:"
                              << empty3d::degeneracy_count << std::endl;
                    std::cout << "exact count: "
                              << empty3d::exact_count << std::endl;
                }
                ENSURE(vert); // bad if we can't find a common vertex
                // then, find the corresponding OVptr, and connect
                for(uint k=0; k<3; k++) {
                    if(overts[k]->concrete == vert) {
                        ie->ends[1] = overts[k];
                        overts[k]->edges.push_back(ie);
                        break;
                    }
                }
            }
        }
        
        ENSURE(isValid());
    }
    
    bool isValid() const {
        ENSURE(the_tri);
        
        return true;
    }
    
    // Subdivision is split in three phases so that the expensive middle
    // one (constrained Delaunay triangulation) can run in parallel:
    //   prepare  (serial)   -- gathers points, splits edges (pool allocs)
    //   triangulate (parallel) -- pure function of SubdivData, no shared state
    //   finish   (serial)   -- creates generic triangles (pool allocs)
    struct SubdivData {
        ShortVec<GVptr, 7>  points;
        ShortVec<GEptr, 8>  edges;
        Vec3d               origin;
        Vec3d               u, v;           // orthonormal basis of the face
        std::vector<int>    tris;           // 3 indices per output triangle
    };

    static Vec2d proj2(const SubdivData &d, const Vec3d &p) {
        Vec3d w = p - d.origin;
        return Vec2d(dot(w, d.u), dot(w, d.v));
    }

    static void choose_face_basis(SubdivData &d) {
        const Vec3d p0 = d.points[0]->coord;
        const Vec3d e0 = d.points[1]->coord - p0;
        const Vec3d e1 = d.points[2]->coord - p0;
        const Vec3d n  = cross(e0, e1);
        d.origin = p0;
        const double n2 = len2(n);
        const double e00 = len2(e0);
        if (n2 > 1e-30 * (e00 + 1.0) && e00 > 1e-30) {
            d.u = e0 / sqrt(e00);
            d.v = normalized(cross(n, d.u));
            return;
        }
        // Degenerate face: axis-aligned drop, oriented CCW when possible.
        uint normdim = maxDim(abs(n) + Vec3d(1e-30, 1e-31, 1e-32));
        d.u = Vec3d(0, 0, 0);
        d.v = Vec3d(0, 0, 0);
        d.u.v[(normdim + 1) % 3] = 1.0;
        d.v.v[(normdim + 2) % 3] = (n.v[normdim] < 0.0) ? -1.0 : 1.0;
    }

    // Make the 2D PSLG a real planar arrangement: unique points, no
    // zero-length or duplicate segments, T-junctions split, crossings
    // split at an existing or newly interpolated vertex. Triangle's
    // insertsegment/locate die on the uncleaned axis-drop graph.
    void sanitize_pslg(IsctProblem *iprob, SubdivData &d) {
        std::vector<GVptr> pts;
        pts.reserve(d.points.size() + 4);
        for (GVptr p : d.points) pts.push_back(p);
        std::vector<GEptr> eds;
        eds.reserve(d.edges.size() + 8);
        for (GEptr e : d.edges) eds.push_back(e);

        auto xy = [&](GVptr p) { return proj2(d, p->coord); };

        double extent = 0.0;
        for (GVptr p : pts) {
            Vec2d q = xy(p);
            extent = std::max(extent, std::max(std::fabs(q.x), std::fabs(q.y)));
        }
        const double eps  = std::max(1e-14, 1e-10 * std::max(extent, 1.0));
        const double eps2 = eps * eps;

        auto on_seg = [&](GVptr a, GVptr b, GVptr c) -> bool {
            if (c == a || c == b) return false;
            Vec2d A = xy(a), B = xy(b), C = xy(c);
            Vec2d ab = B - A, ac = C - A;
            double ab2 = len2(ab);
            if (ab2 <= eps2) return false;
            double dist = std::fabs(ab.x * ac.y - ab.y * ac.x) / std::sqrt(ab2);
            if (dist > eps) return false;
            double t = dot(ac, ab) / ab2;
            return t > 1e-8 && t < 1.0 - 1e-8;
        };

        auto proper_isct = [&](GVptr a, GVptr b, GVptr c, GVptr d,
                               double &t, double &s) -> bool {
            if (a == c || a == d || b == c || b == d) return false;
            Vec2d A = xy(a), B = xy(b), C = xy(c), D = xy(d);
            Vec2d ab = B - A, cd = D - C, ac = C - A;
            double den = ab.x * cd.y - ab.y * cd.x;
            if (std::fabs(den) <= eps2) return false;
            t = (ac.x * cd.y - ac.y * cd.x) / den;
            s = (ac.x * ab.y - ac.y * ab.x) / den;
            return t > 1e-8 && t < 1.0 - 1e-8 && s > 1e-8 && s < 1.0 - 1e-8;
        };

        auto drop_bad_edges = [&]() {
            std::vector<GEptr> keep;
            keep.reserve(eds.size());
            for (GEptr e : eds) {
                if (!e->ends[0] || !e->ends[1] || e->ends[0] == e->ends[1])
                    continue;
                if (len2(xy(e->ends[1]) - xy(e->ends[0])) <= eps2)
                    continue;
                bool dup = false;
                for (GEptr k : keep) {
                    if ((k->ends[0] == e->ends[0] && k->ends[1] == e->ends[1]) ||
                        (k->ends[0] == e->ends[1] && k->ends[1] == e->ends[0])) {
                        dup = true;
                        break;
                    }
                }
                if (!dup) keep.push_back(e);
            }
            eds.swap(keep);
        };

        bool changed = true;
        for (int pass = 0; changed && pass < 24; ++pass) {
            changed = false;

            // Merge 2D-coincident vertices (keep the earlier one).
            for (size_t i = 0; i < pts.size(); ++i) {
                Vec2d A = xy(pts[i]);
                for (size_t j = i + 1; j < pts.size(); ++j) {
                    if (len2(xy(pts[j]) - A) > eps2) continue;
                    GVptr keep = pts[i], drop = pts[j];
                    for (GEptr e : eds) {
                        if (e->ends[0] == drop) e->ends[0] = keep;
                        if (e->ends[1] == drop) e->ends[1] = keep;
                    }
                    pts.erase(pts.begin() + (std::ptrdiff_t)j);
                    changed = true;
                    --j;
                }
            }
            drop_bad_edges();

            // T-junctions: an existing vertex sits on a segment.
            for (size_t ei = 0; ei < eds.size(); ++ei) {
                GEptr e = eds[ei];
                GVptr a = e->ends[0], b = e->ends[1];
                bool split = false;
                for (GVptr c : pts) {
                    if (!on_seg(a, b, c)) continue;
                    SEptr s0 = iprob->newSplitEdge(a, c, e->boundary);
                    SEptr s1 = iprob->newSplitEdge(c, b, e->boundary);
                    eds[ei] = s0;
                    eds.push_back(s1);
                    changed = true;
                    split = true;
                    break;
                }
                if (split) break; // restart scans after topology change
            }
            if (changed) continue;

            // Crossings: two segments meet in the face without a vertex.
            bool did_cross = false;
            for (size_t i = 0; i < eds.size() && !did_cross; ++i) {
                for (size_t j = i + 1; j < eds.size(); ++j) {
                    double t = 0, s = 0;
                    GVptr a = eds[i]->ends[0], b = eds[i]->ends[1];
                    GVptr c = eds[j]->ends[0], dd = eds[j]->ends[1];
                    if (!proper_isct(a, b, c, dd, t, s)) continue;

                    Vec3d p3 = a->coord + t * (b->coord - a->coord);
                    GVptr hit = nullptr;
                    Vec2d Wp = proj2(d, p3);
                    for (GVptr p : pts) {
                        if (len2(xy(p) - Wp) <= eps2) { hit = p; break; }
                    }
                    if (!hit) {
                        GluePt g = iprob->newGluePt();
                        g->split_type = true;
                        g->e = the_tri->edges[0];
                        IVptr iv = iprob->newSplitIsctVert(p3, g);
                        iv->boundary = false;
                        hit = iv;
                        pts.push_back(hit);
                    }
                    bool bi = eds[i]->boundary, bj = eds[j]->boundary;
                    eds[i] = iprob->newSplitEdge(a, hit, bi);
                    eds.push_back(iprob->newSplitEdge(hit, b, bi));
                    eds[j] = iprob->newSplitEdge(c, hit, bj);
                    eds.push_back(iprob->newSplitEdge(hit, dd, bj));
                    changed = true;
                    did_cross = true;
                    break;
                }
            }
        }
        drop_bad_edges();

        d.points.resize(0);
        for (GVptr p : pts) d.points.push_back(p);
        d.edges.resize(0);
        for (GEptr e : eds) d.edges.push_back(e);
        for (uint i = 0; i < d.points.size(); i++) d.points[i]->idx = i;
        for (uint i = 0; i < d.edges.size(); i++)  d.edges[i]->idx  = i;
    }

    void subdivide_prepare(IsctProblem *iprob, SubdivData &d) {
        ShortVec<GVptr, 7> &points = d.points;
        for (uint k = 0; k < 3; k++)
            points.push_back(overts[k]);
        for (IVptr iv : iverts)
            points.push_back(iv);

        ShortVec<GEptr, 8> &edges = d.edges;
        for (uint k = 0; k < 3; k++) {
            subdivideEdge(iprob, oedges[k], edges);
            oedges[k] = nullptr;
        }
        for (IEptr &ie : iedges) {
            subdivideEdge(iprob, ie, edges);
            ie = nullptr;
        }

        choose_face_basis(d);
        sanitize_pslg(iprob, d);
    }

    static void subdivide_triangulate(SubdivData &d) {
        const ShortVec<GVptr, 7> &points = d.points;
        const ShortVec<GEptr, 8> &edges  = d.edges;
        struct triangulateio in, out;

        in.numberofpoints           = (int)points.size();
        in.numberofpointattributes  = 0;
        in.pointlist                = (double*)malloc(sizeof(double) * in.numberofpoints * 2);
        in.pointattributelist       = nullptr;
        in.pointmarkerlist          = (int*)malloc(sizeof(int) * in.numberofpoints);
        for (int k = 0; k < in.numberofpoints; k++) {
            Vec2d q = proj2(d, points[k]->coord);
            in.pointlist[k*2 + 0] = q.x;
            in.pointlist[k*2 + 1] = q.y;
            in.pointmarkerlist[k] = (points[k]->boundary)? 1 : 0;
            points[k]->idx = (uint)k;
        }

        in.numberofsegments = (int)edges.size();
        in.numberofholes = 0;
        in.numberofregions = 0;
        in.segmentlist = (int*)malloc(sizeof(int) * in.numberofsegments * 2);
        in.segmentmarkerlist = (int*)malloc(sizeof(int) * in.numberofsegments);
        for (int k = 0; k < in.numberofsegments; k++) {
            in.segmentlist[k*2 + 0] = (int)edges[k]->ends[0]->idx;
            in.segmentlist[k*2 + 1] = (int)edges[k]->ends[1]->idx;
            in.segmentmarkerlist[k] = (edges[k]->boundary)? 1 : 0;
        }

        in.numberoftriangles = 0;
        in.numberoftriangleattributes = 0;

        out.pointlist = nullptr;
        out.pointattributelist = nullptr;
        out.pointmarkerlist = nullptr;
        out.trianglelist = nullptr;
        out.segmentlist = nullptr;
        out.segmentmarkerlist = nullptr;

        char *params = (char*)("pzQYY");
        triangulate(params, &in, &out, nullptr);

        ENSURE(out.numberofpoints == in.numberofpoints);
        d.tris.assign(out.trianglelist, out.trianglelist + 3*out.numberoftriangles);

        free(in.pointlist);
        free(in.pointmarkerlist);
        free(in.segmentlist);
        free(in.segmentmarkerlist);
        free(out.pointlist);
        free(out.pointmarkerlist);
        free(out.trianglelist);
        free(out.segmentlist);
        free(out.segmentmarkerlist);
    }

    void subdivide_finish(IsctProblem *iprob, SubdivData &d) {
        const ShortVec<GVptr, 7> &points = d.points;
        const uint ntri = (uint)(d.tris.size() / 3);
        gtris.resize(ntri);
        for(uint k=0; k<ntri; k++) {
            GVptr       gv0         = points[d.tris[(k*3)+0]];
            GVptr       gv1         = points[d.tris[(k*3)+1]];
            GVptr       gv2         = points[d.tris[(k*3)+2]];
                        gtris[k]    = iprob->newGenericTri(gv0, gv1, gv2);
        }
    }

    void subdivide(IsctProblem *iprob) {
        SubdivData d;
        subdivide_prepare(iprob, d);
        subdivide_triangulate(d);
        subdivide_finish(iprob, d);
    }

private:
    void subdivideEdge(IsctProblem *iprob, GEptr ge, ShortVec<GEptr, 8> &edges)
    {
        if(ge->interior.size() == 0) {
            //if(typeid(ge) == typeid(IEptr)) { // generate new edge
            //    iprob->buildConcreteEdge(ge);
            //}
            edges.push_back(ge);
        } else if(ge->interior.size() == 1) { // common case
            SEptr       se0     = iprob->newSplitEdge(ge->ends[0],
                                                      ge->interior[0],
                                                      ge->boundary);
            SEptr       se1     = iprob->newSplitEdge(ge->interior[0],
                                                      ge->ends[1],
                                                      ge->boundary);
                        //iprob->buildConcreteEdge(se0);
                        //iprob->buildConcreteEdge(se1);
                        edges.push_back(se0);
                        edges.push_back(se1);
            
            // get rid of old edge
                        iprob->releaseEdge(ge);
        } else { // sorting is the uncommon case
            // Parameterize along the actual 3D edge. Sorting on a single
            // axis (old code) reverses or ties when the edge is diagonal
            // or nearly perpendicular to its longest component, which
            // then feeds Triangle a self-overlapping PSLG.
            Vec3d       dir     = ge->ends[1]->coord - ge->ends[0]->coord;
            Vec3d       origin  = ge->ends[0]->coord;
            
            std::vector< std::pair<double,IVptr> > verts;
            for(IVptr iv : ge->interior) {
                        verts.push_back(std::make_pair(
                            dot(iv->coord - origin, dir),
                            iv
                        ));
            }
            // ... and sort the vector
                        std::sort(verts.begin(), verts.end());
            // then, write the verts into a new container with the endpoints
            std::vector<GVptr>  allv(verts.size()+2);
                        allv[0]             = ge->ends[0];
                        allv[allv.size()-1] = ge->ends[1];
            for(uint k=0; k<verts.size(); k++)
                        allv[k+1]           = verts[k].second;
            
            // now create and accumulate new split edges
            for(uint i=1; i<allv.size(); i++) {
                SEptr   se      = iprob->newSplitEdge(allv[i-1],
                                                      allv[i],
                                                      ge->boundary);
                        edges.push_back(se);
            }
            // get rid of old edge
                        iprob->releaseEdge(ge);
        }
    }
    
public: // data
    ShortVec<IVptr, 4>      iverts;
    ShortVec<IEptr, 2>      iedges;
    // original triangle elements
    OVptr                   overts[3];
    OEptr                   oedges[3];
    
    ShortVec<GTptr, 8>      gtris;
    
    Tptr                    the_tri;
};


template<class VertData, class TriData>
class Mesh<VertData,TriData>::IsctProblem : public TopoCache
{
public:
    IsctProblem(Mesh *owner) : TopoCache(owner, /*vertEdges=*/false, /*vertTris=*/false)
    {
        // (TopoCache::init already leaves every t->data == nullptr)
        
        // Callibrate the quantization unit...
        double maxMag = 0.0;
        {
            cork_par::Local<double> acc([] { return 0.0; });
            const auto &vs = TopoCache::mesh->verts;
            cork_par::for_range(vs.size(), 8192, [&](size_t b, size_t e) {
                double m = 0.0;
                for (size_t i = b; i < e; ++i) m = std::max(m, max(abs(vs[i].pos)));
                acc.local() = std::max(acc.local(), m);
            });
            acc.combine_each([&](double m) { maxMag = std::max(maxMag, m); });
        }
        quantization::calibrate(maxMag);
        
        // and use vertex auxiliary data to store quantized vertex coordinates
        const size_t nv = (size_t)TopoCache::vbulk.n;
        quantized_coords.resize(nv);
        Mesh *m = TopoCache::mesh;
        auto vb = TopoCache::vbulk;
        cork_par::for_each_idx(nv, 8192, [&](size_t i) {
            Vptr v = vb[i];
            Vec3d raw = m->verts[v->ref].pos;
            quantized_coords[i].x = quantization::quantize(raw.x);
            quantized_coords[i].y = quantization::quantize(raw.y);
            quantized_coords[i].z = quantization::quantize(raw.z);
            v->data = &(quantized_coords[i]);
        });
    }
    
    virtual ~IsctProblem() {
        // IsctProblem never filled vert incidence lists (inline-empty ShortVecs)
        // and never free()'d individual edges, so those pools can drop memory
        // after a parallel dtor sweep of the original bulk (edges) or none (verts).
        cork_par::invoke(
            [&] { glue_pts.release(); },
            [&] { tprobs.release(); },
            [&] { ivpool.release(); },
            [&] { ovpool.release(); },
            [&] { iepool.release(); },
            [&] { oepool.release(); },
            [&] { sepool.release(); },
            [&] { gtpool.release(); },
            [&] { TopoCache::verts.release_memory(); },
            [&] { TopoCache::edges.release_bulk(TopoCache::ebulk); },
            [&] { TopoCache::tris.release_memory(); });
    }
    
    // access auxiliary quantized coordinates
    inline Vec3d vPos(Vptr v) const {
        return *(reinterpret_cast<Vec3d*>(v->data));
    }
    
    Tprob getTprob(Tptr t) {
        Tprob prob = reinterpret_cast<Tprob>(t->data);
        if(!prob) {
            t->data = prob = tprobs.alloc();
            prob->init(this, t);
        }
        return prob;
    }
    GluePt newGluePt() {
        GluePt glue = glue_pts.alloc();
        glue->split_type = false;
        return glue;
    }
    
    inline IVptr newIsctVert(Eptr e, Tptr t, GluePt glue) {
        IVptr       iv                  = ivpool.alloc();
                    iv->concrete        = nullptr;
                    iv->coord           = computeCoords(e, t);
                    iv->glue_marker     = glue;
                    glue->copies.push_back(iv);
        return      iv;
    }
    inline IVptr newIsctVert(Tptr t0, Tptr t1, Tptr t2, GluePt glue) {
        IVptr       iv                  = ivpool.alloc();
                    iv->concrete        = nullptr;
                    iv->coord           = computeCoords(t0, t1, t2);
                    iv->glue_marker     = glue;
                    glue->copies.push_back(iv);
        return      iv;
    }
    inline IVptr newSplitIsctVert(Vec3d coords, GluePt glue) {
        IVptr       iv                  = ivpool.alloc();
                    iv->concrete        = nullptr;
                    iv->coord           = coords;
                    iv->glue_marker     = glue;
                    glue->copies.push_back(iv);
        return      iv;
    }
    inline IVptr copyIsctVert(IVptr orig) {
        IVptr       iv                  = ivpool.alloc();
                    iv->concrete        = nullptr;
                    iv->coord           = orig->coord;
                    iv->glue_marker     = orig->glue_marker;
                    orig->glue_marker->copies.push_back(iv);
        return      iv;
    }
    inline IEptr newIsctEdge(IVptr endpoint, Tptr tri_key) {
        IEptr       ie                  = iepool.alloc();
                    ie->concrete        = nullptr;
                    ie->boundary        = false;
                    ie->ends[0]         = endpoint;
                    endpoint->edges.push_back(ie);
                    ie->ends[1]         = nullptr; // other end null
                    ie->other_tri_key   = tri_key;
        return      ie;
    }
    
    inline OVptr newOrigVert(Vptr v) {
        OVptr       ov                  = ovpool.alloc();
                    ov->concrete        = v;
                    ov->coord           = vPos(v);
                    ov->boundary        = true;
        return      ov;
    }
    inline OEptr newOrigEdge(Eptr e, OVptr v0, OVptr v1) {
        OEptr       oe                  = oepool.alloc();
                    oe->concrete        = e;
                    oe->boundary        = true;
                    oe->ends[0]         = v0;
                    oe->ends[1]         = v1;
                    v0->edges.push_back(oe);
                    v1->edges.push_back(oe);
        return      oe;
    }
    inline SEptr newSplitEdge(GVptr v0, GVptr v1, bool boundary) {
        SEptr       se                  = sepool.alloc();
                    se->concrete        = nullptr;
                    se->boundary        = boundary;
                    se->ends[0]         = v0;
                    se->ends[1]         = v1;
                    v0->edges.push_back(se);
                    v1->edges.push_back(se);
        return      se;
    }
    
    inline GTptr newGenericTri(GVptr v0, GVptr v1, GVptr v2) {
        GTptr       gt                  = gtpool.alloc();
                    gt->verts[0]        = v0;
                    gt->verts[1]        = v1;
                    gt->verts[2]        = v2;
                    gt->concrete        = nullptr;
        return      gt;
    }
    
    inline void releaseEdge(GEptr ge) {
        disconnectGE(ge);
        IEptr       ie      = dynamic_cast<IEptr>(ge);
        if(ie) {
                    iepool.free(ie);
        } else {
            OEptr   oe      = dynamic_cast<OEptr>(ge);
                    ENSURE(oe);
                    oepool.free(oe);
        }
    }
    
    inline void killIsctVert(IVptr iv) {
        iv->glue_marker->copies.erase(iv);
        if(iv->glue_marker->copies.size() == 0)
            glue_pts.free(iv->glue_marker);
        
        for(GEptr ge : iv->edges) {
            // disconnect
            ge->interior.erase(iv);
            if(ge->ends[0] == iv)   ge->ends[0] = nullptr;
            if(ge->ends[1] == iv)   ge->ends[1] = nullptr;
        }
        
        ivpool.free(iv);
    }
    
    inline void killIsctEdge(IEptr ie) {
        // an endpoint may be an original vertex
        if(ie->ends[1])
            ie->ends[1]->edges.erase(ie);
        iepool.free(ie);
    }
    
    inline void killOrigVert(OVptr ov) {
        ovpool.free(ov);
    }
    inline void killOrigEdge(OEptr oe) {
        oepool.free(oe);
    }
    
    /*
    inline void buildConcreteVert(GVptr gv) {
                    gv->concrete        = newVert();
                    // NEED DATA SOMEHOW??
    }
    // make sure the endpoints have concrete versions first!
    inline void buildConcreteEdge(GEptr ge) {
        Eptr        e   = ge->concrete  = newEdge();
        Vptr        v0  = e->verts[0]   = ge->ends[0]->concrete;
        Vptr        v1  = e->verts[1]   = ge->ends[1]->concrete;
                    v0->edges.push_back(e);
                    v1->edges.push_back(e);
    }
    inline void buildConcreteTri(GVptr gv0, GVptr gv1, GVptr gv2) {
        // create edges as necessary...
    }*/
    
    bool hasIntersections(); // test for iscts, exit if one is found
    
    void findIntersections();
    void resolveAllIntersections();
private:
    // if we encounter ambiguous degeneracies, then this
    // routine returns false, indicating that the computation aborted.
    bool tryToFindIntersections();
    // In that case, we can perturb the positions of points
    void perturbPositions();
    // in order to give things another try, discard partial work
    void reset();
public:
    
    void dumpIsctPoints(std::vector<Vec3d> *points);
    void dumpIsctEdges(std::vector< std::pair<Vec3d,Vec3d> > *edges);
    
protected: // DATA
    IterPool<GluePointMarker>   glue_pts;
    IterPool<TriangleProblem>   tprobs;
    
    IterPool<IsctVertType>      ivpool;
    IterPool<OrigVertType>      ovpool;
    IterPool<IsctEdgeType>      iepool;
    IterPool<OrigEdgeType>      oepool;
    IterPool<SplitEdgeType>     sepool;
    IterPool<GenericTriType>    gtpool;
private:
    std::vector<Vec3d>          quantized_coords;
private:
    inline void for_edge_tri(std::function<bool(Eptr e, Tptr t)>);
    inline void bvh_edge_tri(std::function<bool(Eptr e, Tptr t)>);

    inline GeomBlob<Eptr> edge_blob(Eptr e);
    inline BBox3d bboxFromTptr(Tptr t);
    
    inline BBox3d buildBox(Eptr e) const;
    inline BBox3d buildBox(Tptr t) const;
    
    inline void marshallArithmeticInput(empty3d::TriIn &input, Tptr t) const;
    inline void marshallArithmeticInput(empty3d::EdgeIn &input, Eptr e) const;
    inline void marshallArithmeticInput(
            empty3d::TriEdgeIn &input, Eptr e, Tptr t) const;
    inline void marshallArithmeticInput(
            empty3d::TriTriTriIn &input, Tptr t0, Tptr t1, Tptr t2) const;
    
    bool checkIsct(Eptr e, Tptr t) const;
    bool checkIsct(Tptr t0, Tptr t1, Tptr t2) const;
    
public:
    Vec3d computeCoords(Eptr e, Tptr t) const;
    Vec3d computeCoords(Tptr t0, Tptr t1, Tptr t2) const;
private:
    
    void fillOutVertData(GluePt glue, VertData &data);
    void fillOutTriData(Tptr tri, Tptr parent);
private:
    class EdgeCache;
    
private: // functions here to get around a GCC bug...
    void createRealPtFromGluePt(GluePt glue);
    void createRealTriangles(Tprob tprob, EdgeCache &ecache);
    long long crt_ns[4] = {0,0,0,0};
};

template<class T, uint LEN> inline
void for_pairs(
    ShortVec<T,LEN> &vec,
    std::function<void(T&,T&)> func
) {
    for(uint i=0; i<vec.size(); i++)
        for(uint j=i+1; j<vec.size(); j++)
            func(vec[i], vec[j]);
}

struct TriTripleTemp
{
    Tptr t0, t1, t2;
    TriTripleTemp(Tptr tp0, Tptr tp1, Tptr tp2) :
        t0(tp0), t1(tp1), t2(tp2)
    {}
};

template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::IsctProblem::for_edge_tri(
    std::function<bool(Eptr e, Tptr t)> func
) {
    bool aborted = false;
    TopoCache::tris.for_each([&](Tptr t) {
        TopoCache::edges.for_each([&](Eptr e) {
            if(!aborted) {
                if(!func(e, t))
                    aborted = true;
            }
        });
    });
}

template<class VertData, class TriData> inline
GeomBlob<Eptr> Mesh<VertData,TriData>::IsctProblem::edge_blob(
    Eptr e
) {
    GeomBlob<Eptr>  blob;
    blob.bbox = buildBox(e);
    blob.point = (blob.bbox.minp + blob.bbox.maxp) / 2.0;
    blob.id = e;
    return blob;
}

template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::IsctProblem::bvh_edge_tri(
    std::function<bool(Eptr e, Tptr t)> func
) {
    std::vector< GeomBlob<Eptr> > edge_geoms;
    edge_geoms.reserve(TopoCache::edges.size());
    TopoCache::edges.for_each([&](Eptr e) {
        edge_geoms.push_back(edge_blob(e));
    });
    AABVH<Eptr> edgeBVH(edge_geoms);

    bool aborted = false;
    TopoCache::tris.for_each([&](Tptr t) {
        if(aborted) return;
        BBox3d bbox = buildBox(t);
        edgeBVH.for_each_in_box(bbox, [&](Eptr e) {
            if(aborted) return;
            if(!func(e, t))
                aborted = true;
        });
    });
}

// ---------------------------------------------------------------------------
// Broad phase: fine uniform cells, sorted entries, merge join.
//
// A dense 3D grid is a poor fit for a triangle *surface* (2D data in 3D):
// almost all cells are empty while the few occupied ones are packed.  Here
// we choose a cell size on the order of an edge, bin every edge and every
// triangle into the cells its bbox touches, sort the (cell, prim) entries and
// merge-join the two sorted streams.  A pair is tested in exactly one cell
// (the cell containing the lower corner of the two boxes' intersection), so
// no de-duplication is needed.  The set of pairs that reach the exact test
// is exactly the set of edge/tri pairs with overlapping boxes -- identical to
// what the original BVH traversal produced.
// ---------------------------------------------------------------------------
namespace cork_si {

// Scratch array without value-initialisation: std::vector<T>(n) zero-fills
// (or default-constructs) hundreds of MB serially here; every element is
// written by a parallel loop before it is read, so skip that.
template<class T>
struct RawArray {
    T *p = nullptr; size_t n = 0;
    RawArray() {}
    explicit RawArray(size_t n_) { alloc(n_); }
    RawArray(const RawArray&) = delete;
    RawArray &operator=(const RawArray&) = delete;
    ~RawArray() { release(); }
    void alloc(size_t n_) {
        release(); n = n_;
        p = static_cast<T*>(::operator new(sizeof(T) * (n_ ? n_ : 1)));
    }
    void release() { if (p) ::operator delete(p); p = nullptr; n = 0; }
    inline T &operator[](size_t i)             { return p[i]; }
    inline const T &operator[](size_t i) const { return p[i]; }
    size_t size() const { return n; }
    T *data() { return p; }
    const T *data() const { return p; }
    T *begin() { return p; }  T *end() { return p + n; }
    const T *begin() const { return p; }  const T *end() const { return p + n; }
};

struct CellGrid {
    Vec3d   org;
    Vec3d   inv;
    int     nx = 1, ny = 1, nz = 1;

    inline void cellOf(const Vec3d &p, int &x, int &y, int &z) const {
        x = (int)std::floor((p.x - org.x) * inv.x);
        y = (int)std::floor((p.y - org.y) * inv.y);
        z = (int)std::floor((p.z - org.z) * inv.z);
        x = x < 0 ? 0 : (x >= nx ? nx - 1 : x);
        y = y < 0 ? 0 : (y >= ny ? ny - 1 : y);
        z = z < 0 ? 0 : (z >= nz ? nz - 1 : z);
    }
    inline uint32_t key(int x, int y, int z) const {
        return (uint32_t)(((uint64_t)z * (uint64_t)ny + (uint64_t)y) * (uint64_t)nx + (uint64_t)x);
    }
    inline uint32_t ownerKey(const BBox3d &a, const BBox3d &b) const {
        Vec3d L = max(a.minp, b.minp);
        int x, y, z;
        cellOf(L, x, y, z);
        return key(x, y, z);
    }
};

inline uint32_t entryKey(uint64_t e) { return (uint32_t)(e >> 32); }
inline uint32_t entryIdx(uint64_t e) { return (uint32_t)(e & 0xffffffffu); }

// sorted (cell<<32 | prim) entries for a set of boxes
template<class BoxArray>
inline void buildEntries(const CellGrid &g, const BoxArray &boxes,
                         RawArray<uint64_t> &out)
{
    const size_t n = boxes.size();
    RawArray<uint32_t> cnt(n + 1);
    cnt[0] = 0;
    cork_par::for_each_idx(n, 4096, [&](size_t i) {
        int a[3], b[3];
        g.cellOf(boxes[i].minp, a[0], a[1], a[2]);
        g.cellOf(boxes[i].maxp, b[0], b[1], b[2]);
        cnt[i + 1] = (uint32_t)((b[0]-a[0]+1) * (b[1]-a[1]+1) * (b[2]-a[2]+1));
    });
    cork_par::prefix_sum_inplace(cnt.data(), n + 1);
    out.alloc(cnt[n]);
    cork_par::for_each_idx(n, 4096, [&](size_t i) {
        int a[3], b[3];
        g.cellOf(boxes[i].minp, a[0], a[1], a[2]);
        g.cellOf(boxes[i].maxp, b[0], b[1], b[2]);
        size_t w = cnt[i];
        for (int z = a[2]; z <= b[2]; ++z)
            for (int y = a[1]; y <= b[1]; ++y)
                for (int x = a[0]; x <= b[0]; ++x)
                    out[w++] = ((uint64_t)g.key(x, y, z) << 32) | (uint64_t)i;
    });
    cork_par::sort(out.begin(), out.end()); // tbb parallel_sort; 16-byte-key radix was slower here
}

// chunk the sorted tri entries without splitting a cell run
inline void chunkRuns(const RawArray<uint64_t> &ent, size_t chunk,
                      std::vector<size_t> &starts)
{
    starts.clear();
    starts.push_back(0);
    const size_t n = ent.size();
    for (size_t s = chunk; s < n; s += chunk) {
        size_t i = s;
        while (i < n && entryKey(ent[i]) == entryKey(ent[i - 1])) ++i;
        if (i < n && i > starts.back()) starts.push_back(i);
    }
    starts.push_back(n);
}

struct EdgeTriHit {
    uint32_t ti, ei;
    Vec3d    coord;
};

// float box rounded outward: overlap test on these is conservative w.r.t.
// the exact double test, which is still applied to every pair that passes.
struct BoxF { float mnx, mny, mnz, mxx, mxy, mxz; };
inline float fdown(double v) { float f = (float)v; return ((double)f > v) ? std::nextafter(f, -INFINITY) : f; }
inline float fup(double v)   { float f = (float)v; return ((double)f < v) ? std::nextafter(f,  INFINITY) : f; }
inline BoxF toBoxF(const BBox3d &b) {
    BoxF f;
    f.mnx = fdown(b.minp.x); f.mny = fdown(b.minp.y); f.mnz = fdown(b.minp.z);
    f.mxx = fup(b.maxp.x);   f.mxy = fup(b.maxp.y);   f.mxz = fup(b.maxp.z);
    return f;
}
inline bool overlapF(const BoxF &a, const BoxF &b) {
    return a.mnx <= b.mxx && b.mnx <= a.mxx &&
           a.mny <= b.mxy && b.mny <= a.mxy &&
           a.mnz <= b.mxz && b.mnz <= a.mxz;
}

// Structure-of-arrays float boxes with an AVX2 range query.
struct BoxSoA {
    RawArray<float> mnx, mny, mnz, mxx, mxy, mxz;
    void resize(size_t n) {
        // +8 padding so the vector tail can be loaded unconditionally;
        // only the padding is initialised here (to never-matching boxes)
        RawArray<float> *all[6] = { &mnx, &mny, &mnz, &mxx, &mxy, &mxz };
        for (int k = 0; k < 6; ++k) {
            all[k]->alloc(n + 8);
            for (size_t i = n; i < n + 8; ++i) (*all[k])[i] = (k < 3) ? INFINITY : -INFINITY;
        }
    }
    inline void set(size_t i, const BoxF &b) {
        mnx[i] = b.mnx; mny[i] = b.mny; mnz[i] = b.mnz;
        mxx[i] = b.mxx; mxy[i] = b.mxy; mxz[i] = b.mxz;
    }
    // indices in [b0,b1) whose box overlaps q -> out ; returns count
    inline int overlaps(const BoxF &q, size_t b0, size_t b1, std::vector<size_t> &out) const {
        if (out.size() < (b1 - b0) + 8) out.resize((b1 - b0) + 8);
        int w = 0;
#if defined(__AVX2__) || defined(_MSC_VER)
        const __m256 qmnx = _mm256_set1_ps(q.mnx), qmny = _mm256_set1_ps(q.mny), qmnz = _mm256_set1_ps(q.mnz);
        const __m256 qmxx = _mm256_set1_ps(q.mxx), qmxy = _mm256_set1_ps(q.mxy), qmxz = _mm256_set1_ps(q.mxz);
        size_t b = b0;
        for (; b < b1; b += 8) {
            __m256 ok = _mm256_cmp_ps(qmnx, _mm256_loadu_ps(&mxx[b]), _CMP_LE_OQ);
            ok = _mm256_and_ps(ok, _mm256_cmp_ps(_mm256_loadu_ps(&mnx[b]), qmxx, _CMP_LE_OQ));
            ok = _mm256_and_ps(ok, _mm256_cmp_ps(qmny, _mm256_loadu_ps(&mxy[b]), _CMP_LE_OQ));
            ok = _mm256_and_ps(ok, _mm256_cmp_ps(_mm256_loadu_ps(&mny[b]), qmxy, _CMP_LE_OQ));
            ok = _mm256_and_ps(ok, _mm256_cmp_ps(qmnz, _mm256_loadu_ps(&mxz[b]), _CMP_LE_OQ));
            ok = _mm256_and_ps(ok, _mm256_cmp_ps(_mm256_loadu_ps(&mnz[b]), qmxz, _CMP_LE_OQ));
            unsigned mask = (unsigned)_mm256_movemask_ps(ok);
            // lanes past b1 read padding (inf boxes) and never match
            while (mask) {
                unsigned long bit;
#if defined(_MSC_VER)
                _BitScanForward(&bit, mask);
#else
                bit = (unsigned long)__builtin_ctz(mask);
#endif
                size_t idx = b + bit;
                if (idx < b1) out[w++] = idx;
                mask &= mask - 1;
            }
        }
#else
        for (size_t b = b0; b < b1; ++b) {
            if (q.mnx <= mxx[b] && mnx[b] <= q.mxx &&
                q.mny <= mxy[b] && mny[b] <= q.mxy &&
                q.mnz <= mxz[b] && mnz[b] <= q.mxz) out[w++] = b;
        }
#endif
        return w;
    }
};

} // namespace cork_si

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::IsctProblem::tryToFindIntersections()
{
    using namespace cork_si;
    empty3d::degeneracy_count = 0;
    empty3d::exact_count = 0;

    std::vector<Eptr> edges;
    std::vector<Tptr> tris;
    {
        CORK_PROF("      collect edges/tris");
        edges.resize(TopoCache::ebulk.n);
        tris.resize(TopoCache::tbulk.n);
        cork_par::for_each_idx((size_t)TopoCache::ebulk.n, 8192, [&](size_t i) {
            edges[i] = TopoCache::ebulk[i];
        });
        cork_par::for_each_idx((size_t)TopoCache::tbulk.n, 8192, [&](size_t i) {
            tris[i] = TopoCache::tbulk[i];
        });
    }
    const size_t ne = edges.size();
    const size_t nt = tris.size();
    if (ne == 0 || nt == 0) return true;

    // ---- boxes (parallel) + world bbox / mean edge extent (reduce) ----
    RawArray<BBox3d> ebb(ne), tbb_(nt);
    // flat copies of what the narrow phase needs (vertex ids + quantized
    // positions) so the hot loop never chases Eptr/Tptr/Vptr pointers
    struct EdgeRec { Vec3d p[2]; uint32_t v[2]; };
    struct TriRec  { Vec3d p[3]; uint32_t v[3]; };
    RawArray<EdgeRec> erec(ne);
    RawArray<TriRec>  trec(nt);
    const Vec3d *qbase = quantized_coords.data();
    BBox3d world;
    double meanExt = 0.0;
    {
        CORK_PROF("      boxes");
        struct Acc { BBox3d box; double sum = 0.0; };
        cork_par::Local<Acc> acc;
        cork_par::for_range(ne, 8192, [&](size_t b, size_t e) {
            Acc &a = acc.local();
            for (size_t i = b; i < e; ++i) {
                ebb[i] = buildBox(edges[i]);
                a.box = convex(a.box, ebb[i]);
                Vec3d d = ebb[i].maxp - ebb[i].minp;
                a.sum += std::max(d.x, std::max(d.y, d.z));
                for (int k = 0; k < 2; ++k) {
                    const Vec3d *q = (const Vec3d*)edges[i]->verts[k]->data;
                    erec[i].p[k] = *q; erec[i].v[k] = (uint32_t)(q - qbase);
                }
            }
        });
        acc.combine_each([&](const Acc &a) { world = convex(world, a.box); meanExt += a.sum; });
        meanExt /= (double)ne;
        cork_par::for_each_idx(nt, 8192, [&](size_t i) {
            tbb_[i] = buildBox(tris[i]);
            for (int k = 0; k < 3; ++k) {
                const Vec3d *q = (const Vec3d*)tris[i]->verts[k]->data;
                trec[i].p[k] = *q; trec[i].v[k] = (uint32_t)(q - qbase);
            }
        });
    }

    const bool useLBVH = std::getenv("CORK_LBVH") != nullptr;
    std::vector<EdgeTriHit> hits;
    int any_degen = 0;

    if (useLBVH) {
        const bool queryTris = std::getenv("CORK_LBVH_TRI") != nullptr;
        cork_lbvh::LBVH tree;
        {
            CORK_PROF("      LBVH build");
            if (queryTris) tree.build(ebb.data(), ne, world);
            else           tree.build(tbb_.data(), nt, world);
            cork_prof::note("      LBVH leaves", (double)tree.n);
        }
        {
            CORK_PROF("      LBVH query + narrow");
            cork_par::Local<std::vector<EdgeTriHit>> tls_hits;
            cork_par::Local<int> tls_degen([] { return 0; });
            cork_par::Local<std::array<long long, 3>> tls_cnt([] { return std::array<long long, 3>{0, 0, 0}; });
            auto handle_pair = [&](uint32_t ei, uint32_t ti, std::vector<EdgeTriHit> &local_hits, std::array<long long, 3> &cnt) {
                const EdgeRec &er = erec[ei];
                const TriRec &tr = trec[ti];
                if (er.v[0] == tr.v[0] || er.v[0] == tr.v[1] || er.v[0] == tr.v[2] ||
                    er.v[1] == tr.v[0] || er.v[1] == tr.v[1] || er.v[1] == tr.v[2]) return;
                if (!hasIsct(tbb_[ti], ebb[ei])) return;
                cnt[0]++;
                empty3d::TriEdgeIn input;
                input.edge.p[0] = er.p[0]; input.edge.p[1] = er.p[1];
                input.tri.p[0] = tr.p[0]; input.tri.p[1] = tr.p[1]; input.tri.p[2] = tr.p[2];
                if (!empty3d::emptyExact(input)) {
                    EdgeTriHit h;
                    h.ti = ti; h.ei = ei;
                    h.coord = empty3d::coordsExact(input);
                    local_hits.push_back(h);
                }
            };
            if (queryTris) {
            cork_par::for_each_idx(nt, 64, [&](size_t ti) {
                auto &local_hits = tls_hits.local();
                auto &cnt = tls_cnt.local();
                empty3d::degeneracy_count = 0;
                cork_lbvh::BoxF q = cork_lbvh::toBoxF(tbb_[ti]);
                tree.query(q, [&](uint32_t ei) { handle_pair(ei, (uint32_t)ti, local_hits, cnt); });
                tls_degen.local() += empty3d::degeneracy_count;
                empty3d::degeneracy_count = 0;
            });
            } else {
            cork_par::for_each_idx(ne, 64, [&](size_t ei) {
                auto &local_hits = tls_hits.local();
                auto &cnt = tls_cnt.local();
                empty3d::degeneracy_count = 0;
                cork_lbvh::BoxF q = cork_lbvh::toBoxF(ebb[ei]);
                tree.query(q, [&](uint32_t ti) { handle_pair((uint32_t)ei, ti, local_hits, cnt); });
                tls_degen.local() += empty3d::degeneracy_count;
                empty3d::degeneracy_count = 0;
            });
            }
            tls_hits.combine_each([&](std::vector<EdgeTriHit> &v) { hits.insert(hits.end(), v.begin(), v.end()); });
            tls_degen.combine_each([&](int d) { any_degen += d; });
            long long tot0 = 0;
            tls_cnt.combine_each([&](const std::array<long long, 3> &c) { tot0 += c[0]; });
            cork_prof::note("      #filter calls (bbox pass)", (double)tot0);
            cork_prof::note("      #edge-tri hits", (double)hits.size());
            cork_par::sort(hits.begin(), hits.end(), [](const EdgeTriHit &a, const EdgeTriHit &b) {
                return a.ti < b.ti || (a.ti == b.ti && a.ei < b.ei);
            });
        }
    } else {

    // ---- cell grid ----
    CellGrid grid;
    {
        Vec3d ext = world.maxp - world.minp;
        double cellFactor = 3.0;   // cells ~3x mean edge extent (measured best)
        if (const char *cf = std::getenv("CORK_CELL_FACTOR")) cellFactor = std::atof(cf);
        double h = cellFactor * meanExt;
        double maxExt = std::max(ext.x, std::max(ext.y, ext.z));
        if (!(h > 0.0)) h = (maxExt > 0.0) ? maxExt : 1.0;
        for (;;) {
            long long nx = std::max(1LL, std::min(1LL << 20, (long long)std::ceil(ext.x / h)));
            long long ny = std::max(1LL, std::min(1LL << 20, (long long)std::ceil(ext.y / h)));
            long long nz = std::max(1LL, std::min(1LL << 20, (long long)std::ceil(ext.z / h)));
            if (nx * ny * nz < (1LL << 32)) { grid.nx = (int)nx; grid.ny = (int)ny; grid.nz = (int)nz; break; }
            h *= 1.5;
        }
        double pad = 1e-9 * std::max(maxExt, 1.0);
        grid.org = world.minp - Vec3d(pad, pad, pad);
        grid.inv = Vec3d(1.0 / h, 1.0 / h, 1.0 / h);
        cork_prof::note("      grid nx", (double)grid.nx);
        cork_prof::note("      grid ny", (double)grid.ny);
        cork_prof::note("      grid nz", (double)grid.nz);
    }

    RawArray<uint64_t> eent, tent;
    {
        CORK_PROF("      bin + sort entries");
        cork_par::invoke([&] { buildEntries(grid, ebb, eent); },
                         [&] { buildEntries(grid, tbb_, tent); });
        cork_prof::note("      #edge entries", (double)eent.size());
        cork_prof::note("      #tri entries", (double)tent.size());
    }
    // edge boxes in sorted-entry order: the inner join loop then streams
    // through contiguous memory instead of gathering random 48 byte boxes
    // edge boxes in sorted-entry order as SoA floats (rounded outward), so the
    // join streams contiguous memory and tests 8 edges per AVX2 instruction
    BoxSoA esoa;
    RawArray<uint32_t> evid0(eent.size()), evid1(eent.size());
    {
        CORK_PROF("      stream boxes (SoA)");
        esoa.resize(eent.size());
        cork_par::for_each_idx(eent.size(), 16384, [&](size_t b) {
            uint32_t ei = entryIdx(eent[b]);
            esoa.set(b, toBoxF(ebb[ei]));
            evid0[b] = erec[ei].v[0];
            evid1[b] = erec[ei].v[1];
        });
    }

    // ---- merge join + narrow phase + exact coordinates (parallel) ----
    {
        CORK_PROF("      parallel edge-tri SI");
        std::vector<size_t> starts;
        chunkRuns(tent, 2048, starts);
        const size_t nchunks = starts.size() - 1;
        const size_t NE = eent.size();

        cork_par::Local<std::vector<EdgeTriHit>> tls_hits;
        cork_par::Local<int> tls_degen([] { return 0; });
        cork_par::Local<std::array<long long, 9>> tls_cnt([] { return std::array<long long, 9>{0, 0, 0, 0, 0, 0, 0, 0, 0}; });
        const bool verify_coords = std::getenv("CORK_VERIFY_COORDS") != nullptr;

        cork_par::Local<std::vector<size_t>> tls_cand;
        cork_par::for_range(nchunks, 1, [&](size_t cb, size_t ce) {
            auto &local_hits = tls_hits.local();
            auto &cnt = tls_cnt.local();
            auto &cand = tls_cand.local();
            auto chunk_t0 = std::chrono::steady_clock::now();
            empty3d::degeneracy_count = 0;
            empty3d::exact_count = 0;
            for (size_t c = cb; c < ce; ++c) {
                size_t i = starts[c], end = starts[c + 1];
                if (i >= end) continue;
                uint32_t k0 = entryKey(tent[i]);
                size_t ep = (size_t)(std::lower_bound(eent.begin(), eent.end(), (uint64_t)k0 << 32) - eent.begin());
                while (i < end) {
                    uint32_t k = entryKey(tent[i]);
                    size_t j = i + 1;
                    while (j < end && entryKey(tent[j]) == k) ++j;
                    while (ep < NE && entryKey(eent[ep]) < k) ++ep;
                    if (ep < NE && entryKey(eent[ep]) == k) {
                        size_t eq = ep + 1;
                        while (eq < NE && entryKey(eent[eq]) == k) ++eq;
                        for (size_t a = i; a < j; ++a) {
                            uint32_t ti = entryIdx(tent[a]);
                            const BBox3d &tb = tbb_[ti];
                            const BoxF tbf = toBoxF(tb);
                            const TriRec &tr = trec[ti];
                            empty3d::TriEdgeIn input;
                            input.tri.p[0] = tr.p[0]; input.tri.p[1] = tr.p[1]; input.tri.p[2] = tr.p[2];
                            cnt[5] += (long long)(eq - ep);
                            // AVX2: candidate bitmask over the cell's edge entries
                            int ncand = esoa.overlaps(tbf, ep, eq, cand);
                            cnt[7] += ncand;
                            for (int ci = 0; ci < ncand; ++ci) {
                                size_t b = cand[ci];
                                // most candidates are edges touching the triangle:
                                // reject them from the streamed vertex ids first
                                const uint32_t ev0 = evid0[b], ev1 = evid1[b];
                                if (ev0 == tr.v[0] || ev0 == tr.v[1] || ev0 == tr.v[2] ||
                                    ev1 == tr.v[0] || ev1 == tr.v[1] || ev1 == tr.v[2]) continue;
                                uint32_t ei = entryIdx(eent[b]);
                                const BBox3d &eb = ebb[ei];
                                if (!hasIsct(tb, eb)) continue;
                                if (grid.ownerKey(tb, eb) != k) continue;
                                const EdgeRec &er = erec[ei];
                                cnt[0]++;
                                input.edge.p[0] = er.p[0]; input.edge.p[1] = er.p[1];
                                if (!empty3d::emptyExact(input)) {
                                    EdgeTriHit h;
                                    h.ti = ti; h.ei = ei;
                                    auto c0 = std::chrono::steady_clock::now();
                                    h.coord = empty3d::coordsExact(input);
                                    cnt[2] += std::chrono::duration_cast<std::chrono::nanoseconds>(
                                        std::chrono::steady_clock::now() - c0).count();
                                    if (verify_coords) {
                                        Vec3d g = empty3d::coordsExactGmp(input);
                                        if (g.x != h.coord.x || g.y != h.coord.y || g.z != h.coord.z)
                                            cnt[3]++;
                                    }
                                    local_hits.push_back(h);
                                }
                            }
                        }
                        ep = eq;
                    }
                    i = j;
                }
            }
            tls_degen.local() += empty3d::degeneracy_count;
            cnt[1] += empty3d::exact_count;
            cnt[4] += std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - chunk_t0).count();
            empty3d::degeneracy_count = 0;
            empty3d::exact_count = 0;
        });

        tls_hits.combine_each([&](std::vector<EdgeTriHit> &v) { hits.insert(hits.end(), v.begin(), v.end()); });
        tls_degen.combine_each([&](int d) { any_degen += d; });
        long long tot[9] = {0,0,0,0,0,0,0,0,0};
        tls_cnt.combine_each([&](const std::array<long long, 9> &c) { for (int i=0;i<9;i++) tot[i] += c[i]; });
        cork_prof::note("      #cell cross pairs", (double)tot[5]);
        cork_prof::note("      #float-box candidates", (double)tot[7]);
        cork_prof::note("      #filter calls (bbox pass)", (double)tot[0]);
        cork_prof::note("      #exact fallbacks", (double)tot[1]);
        cork_prof::note("      coords CPU ms (sum threads)", (double)tot[2] * 1e-6);
        cork_prof::note("      SI busy CPU ms (sum threads)", (double)tot[4] * 1e-6);
        if (verify_coords) cork_prof::note("      #coords != GMP reference", (double)tot[3]);
        cork_prof::note("      #edge-tri hits", (double)hits.size());

        // deterministic order for a given perturbation
        cork_par::sort(hits.begin(), hits.end(), [](const EdgeTriHit &a, const EdgeTriHit &b) {
            return a.ti < b.ti || (a.ti == b.ti && a.ei < b.ei);
        });
    }

    } // !useLBVH

    if (any_degen > 0) {
        empty3d::degeneracy_count = any_degen;
        return false;
    }
    empty3d::degeneracy_count = 0;

    {
        CORK_PROF("      apply hits (glue pts)");
        for (const EdgeTriHit &h : hits) {
            Eptr eisct = edges[h.ei];
            Tptr tisct = tris[h.ti];
            GluePt glue = newGluePt();
            glue->edge_tri_type = true;
            glue->e = eisct;
            glue->t[0] = tisct;
            IVptr iv = getTprob(tisct)->addInteriorEndpoint(this, eisct, glue, h.coord);
            for (Tptr tri : eisct->tris) {
                getTprob(tri)->addBoundaryEndpoint(this, tisct, eisct, iv);
            }
        }
    }

    // ---- tri-tri-tri: collect (serial), test + coords (parallel), apply (serial) ----
    CORK_PROF("      tri-tri-tri");
    std::vector<TriTripleTemp> triples;
    tprobs.for_each([&](Tprob tprob) {
        Tptr t0 = tprob->the_tri;
        for_pairs<IEptr,2>(tprob->iedges, [&](IEptr &ie1, IEptr &ie2){
            Tptr t1 = ie1->other_tri_key;
            Tptr t2 = ie2->other_tri_key;
            if(t0 < t1 && t0 < t2) {
                Tprob prob1 = reinterpret_cast<Tprob>(t1->data);
                for(IEptr ie : prob1->iedges) {
                    if(ie->other_tri_key == t2) {
                        triples.push_back(TriTripleTemp(t0, t1, t2));
                    }
                }
            }
        });
    });
    struct TripleRes { bool ok; Vec3d c0, c1, c2; };
    std::vector<TripleRes> tres(triples.size());
    int tri_degen = 0;
    {
        cork_par::Local<int> tls_degen([] { return 0; });
        cork_par::for_range(triples.size(), 64, [&](size_t b, size_t e) {
            empty3d::degeneracy_count = 0;
            for (size_t i = b; i < e; ++i) {
                const TriTripleTemp &t = triples[i];
                TripleRes &r = tres[i];
                r.ok = checkIsct(t.t0, t.t1, t.t2);
                if (r.ok) {
                    // same orderings as the three addInteriorPoint calls below
                    r.c0 = computeCoords(t.t0, t.t1, t.t2);
                    r.c1 = computeCoords(t.t1, t.t0, t.t2);
                    r.c2 = computeCoords(t.t2, t.t0, t.t1);
                }
            }
            tls_degen.local() += empty3d::degeneracy_count;
            empty3d::degeneracy_count = 0;
        });
        tls_degen.combine_each([&](int d) { tri_degen += d; });
    }
    if (tri_degen > 0) {
        empty3d::degeneracy_count = tri_degen;
        return false;
    }
    for (size_t i = 0; i < triples.size(); ++i) {
        if (!tres[i].ok) continue;
        const TriTripleTemp &t = triples[i];
        GluePt      glue                    = newGluePt();
                    glue->edge_tri_type     = false;
                    glue->t[0]              = t.t0;
                    glue->t[1]              = t.t1;
                    glue->t[2]              = t.t2;
        getTprob(t.t0)->addInteriorPoint(this, t.t1, t.t2, glue, tres[i].c0);
        getTprob(t.t1)->addInteriorPoint(this, t.t0, t.t2, glue, tres[i].c1);
        getTprob(t.t2)->addInteriorPoint(this, t.t0, t.t1, glue, tres[i].c2);
    }

    return true;
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::perturbPositions()
{
    const double EPSILON = 1.0e-5; // perturbation epsilon
    // Per-vertex splitmix: thread-safe, same magnitude as drand(-E,E)^3.
    // (std::rand is not safe to call from many threads.)
    const uint64_t seed = (uint64_t)std::rand() ^ 0xA5A5A5A5A5A5A5A5ULL;
    cork_par::for_each_idx(quantized_coords.size(), 4096, [&](size_t i) {
        uint64_t x = seed + 0x9e3779b97f4a7c15ULL * (i + 1);
        auto u01 = [](uint64_t &s) {
            s += 0x9e3779b97f4a7c15ULL;
            uint64_t z = s;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            z ^= z >> 31;
            return (double)(z >> 11) * (1.0 / 9007199254740992.0);
        };
        Vec3d &coord = quantized_coords[i];
        Vec3d perturbation(quantization::quantize((u01(x) * 2.0 - 1.0) * EPSILON),
                           quantization::quantize((u01(x) * 2.0 - 1.0) * EPSILON),
                           quantization::quantize((u01(x) * 2.0 - 1.0) * EPSILON));
        coord += perturbation;
    });
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::reset()
{
    // the data pointer in the triangles points to tproblems
    // that we're about to destroy,
    // so zero out all those pointers first!
    tprobs.for_each([](Tprob tprob) {
        Tptr t = tprob->the_tri;
        t->data = nullptr;
    });
    
    glue_pts.clear();
    tprobs.clear();
    
    ivpool.clear();
    ovpool.clear();
    iepool.clear();
    oepool.clear();
    sepool.clear();
    gtpool.clear();
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::findIntersections()
{
    int nTrys = 5;
    {
        CORK_PROF("    perturbPositions");
        perturbPositions(); // always perturb for safety...
    }
    while(nTrys > 0) {
        bool ok;
        {
            CORK_PROF("    tryToFindIntersections");
            ok = tryToFindIntersections();
        }
        if(!ok) {
            CORK_PROF("    reset+perturb (retry)");
            reset();
            perturbPositions();
            nTrys--;
        } else {
            break;
        }
    }
    if(nTrys <= 0) {
        CORK_ERROR("Ran out of tries to perturb the mesh");
        exit(1);
    }
    
    // ok all points put together,
    // all triangle problems assembled.
    // Some intersection edges may have original vertices as endpoints
    // we consolidate the problems to check for cases like these.
    CORK_PROF("    consolidate");
    cork_prof::note("    #tprobs", (double)tprobs.size());
    tprobs.for_each([&](Tprob tprob) {
        tprob->consolidate(this);
    });
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::IsctProblem::hasIntersections()
{
    bool foundIsct = false;
    empty3d::degeneracy_count = 0;
    // Find some edge-triangle intersection point...
    bvh_edge_tri([&](Eptr eisct, Tptr tisct)->bool{
      if(checkIsct(eisct,tisct)) {
        foundIsct = true;
        return false; // break;
      }
      if(empty3d::degeneracy_count > 0) {
        return false; // break;
      }
      return true; // continue
    });
    
    if(empty3d::degeneracy_count > 0 || foundIsct) {
        std::cout << "This self-intersection might be spurious. "
                     "Degeneracies were detected." << std::endl;
        return true;
    } else {
        return false;
    }
}


template<class VertData, class TriData> inline
BBox3d Mesh<VertData,TriData>::IsctProblem::buildBox(Eptr e) const
{
    Vec3d p0 = vPos(e->verts[0]);
    Vec3d p1 = vPos(e->verts[1]);
    return BBox3d(min(p0, p1), max(p0, p1));
}
template<class VertData, class TriData> inline
BBox3d Mesh<VertData,TriData>::IsctProblem::buildBox(Tptr t) const
{
    Vec3d p0 = vPos(t->verts[0]);
    Vec3d p1 = vPos(t->verts[1]);
    Vec3d p2 = vPos(t->verts[2]);
    return BBox3d(min(p0, min(p1, p2)), max(p0, max(p1, p2)));
}

template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::IsctProblem::marshallArithmeticInput(
        empty3d::EdgeIn &input, Eptr e
) const {
    input.p[0] = vPos(e->verts[0]);
    input.p[1] = vPos(e->verts[1]);
}
template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::IsctProblem::marshallArithmeticInput(
        empty3d::TriIn &input, Tptr t
) const {
    input.p[0] = vPos(t->verts[0]);
    input.p[1] = vPos(t->verts[1]);
    input.p[2] = vPos(t->verts[2]);
}
template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::IsctProblem::marshallArithmeticInput(
        empty3d::TriEdgeIn &input,
        Eptr e, Tptr t
) const {
    marshallArithmeticInput(input.edge, e);
    marshallArithmeticInput(input.tri, t);
}
template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::IsctProblem::marshallArithmeticInput(
        empty3d::TriTriTriIn &input,
        Tptr t0, Tptr t1, Tptr t2
) const {
    marshallArithmeticInput(input.tri[0], t0);
    marshallArithmeticInput(input.tri[1], t1);
    marshallArithmeticInput(input.tri[2], t2);
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::IsctProblem::checkIsct(Eptr e, Tptr t) const
{
    // BVH already tested boxes for most callers; keep a cheap reject for serial
    BBox3d      ebox        = buildBox(e);
    BBox3d      tbox        = buildBox(t);
    if(!hasIsct(ebox, tbox))
                return      false;

    if(hasCommonVert(e, t))
                return      false;

    empty3d::TriEdgeIn input;
    marshallArithmeticInput(input, e, t);
    bool empty = empty3d::emptyExact(input);
    return !empty;
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::IsctProblem::checkIsct(
    Tptr t0, Tptr t1, Tptr t2
) const {
    // This function should only be called if we've already
    // identified that the intersection edges
    //      (t0,t1), (t0,t2), (t1,t2)
    // exist.
    // From this, we can conclude that each pair of triangles
    // shares no more than a single vertex in common.
    //  If each of these shared vertices is different from each other,
    // then we could legitimately have a triple intersection point,
    // but if all three pairs share the same vertex in common, then
    // the intersection of the three triangles must be that vertex.
    // So, we must check for such a single vertex in common amongst
    // the three triangles
    Vptr common = commonVert(t0, t1);
    if(common) {
        for(uint i=0; i<3; i++)
            if(common == t2->verts[i])
                return      false;
    }
    
    empty3d::TriTriTriIn input;
    marshallArithmeticInput(input, t0, t1, t2);
    //bool empty = empty3d::isEmpty(input);
    bool empty = empty3d::emptyExact(input);
    return !empty;
}

template<class VertData, class TriData>
Vec3d Mesh<VertData,TriData>::IsctProblem::computeCoords(Eptr e, Tptr t) const
{
    empty3d::TriEdgeIn input;
    marshallArithmeticInput(input, e, t);
    Vec3d coords = empty3d::coordsExact(input);
    return coords;
}

template<class VertData, class TriData>
Vec3d Mesh<VertData,TriData>::IsctProblem::computeCoords(
    Tptr t0, Tptr t1, Tptr t2
) const {
    empty3d::TriTriTriIn input;
    marshallArithmeticInput(input, t0, t1, t2);
    Vec3d coords = empty3d::coordsExact(input);
    return coords;
}


template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::fillOutVertData(
    GluePt glue, VertData &data
) {
    if(glue->split_type) { // manually inserted split point
        uint v0i = glue->e->verts[0]->ref;
        uint v1i = glue->e->verts[1]->ref;
        data.isctInterpolate(TopoCache::mesh->verts[v0i],
                             TopoCache::mesh->verts[v1i]);
    } else
    if(glue->edge_tri_type) { // edge-tri type
        IsctVertEdgeTriInput<VertData,TriData>      input;
        for(uint k=0; k<2; k++) {
            uint    vid                 = glue->e->verts[k]->ref;
                    input.e[k]          = &(TopoCache::mesh->verts[vid]);
        }
        for(uint k=0; k<3; k++) {
            uint    vid                 = glue->t[0]->verts[k]->ref;
                    input.t[k]          = &(TopoCache::mesh->verts[vid]);
        }
        data.isct(input);
    } else { // tri-tri-tri type
        IsctVertTriTriTriInput<VertData,TriData>    input;
        for(uint i=0; i<3; i++) {
          for(uint j=0; j<3; j++) {
            uint    vid                 = glue->t[i]->verts[j]->ref;
                    input.t[i][j]       = &(TopoCache::mesh->verts[vid]);
        }}
        data.isct(input);
    }
}

template<class VertData, class TriData> inline
void Mesh<VertData,TriData>::subdivide_tri(
    uint t_piece_ref, uint t_parent_ref
) {
    SubdivideTriInput<VertData,TriData>     input;
                input.pt        = &(tris[t_parent_ref].data);
    for(uint k=0; k<3; k++) {
                input.pv[k]     = &(verts[tris[t_parent_ref].v[k]]);
                input.v[k]      = &(verts[tris[t_piece_ref].v[k]]);
    }
    tris[t_piece_ref].data.subdivide(input);
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::fillOutTriData(
    Tptr piece, Tptr parent
) {
    TopoCache::mesh->subdivide_tri(piece->ref, parent->ref);
}


// Open-addressing hash from (min ref, max ref) -> Eptr.
// The original kept a ShortVec<_,8> per mesh vertex (~150 bytes each, so
// ~100MB of zero-initialised memory for a 700k vertex mesh) even though only
// the handful of vertices touched by intersections ever hold an entry.
template<class VertData, class TriData>
class Mesh<VertData,TriData>::IsctProblem::EdgeCache
{
public:
    EdgeCache(IsctProblem *ip, uint expected = 1024) : iprob(ip) {
        size_t cap = 16;
        while(cap < (size_t)expected * 2) cap <<= 1;
        keys.assign(cap, (uint64_t)EMPTY);
        vals.assign(cap, nullptr);
        mask = cap - 1;
    }
    
    Eptr operator()(Vptr v0, Vptr v1) {
        uint i = v0->ref;
        uint j = v1->ref;
        if(i > j) std::swap(i,j);
        uint64_t key = ((uint64_t)i << 32) | (uint64_t)j;
        size_t slot = find(key);
        if(keys[slot] == key) return vals[slot];
        // if not existing, create it
        Eptr e = iprob->newEdge();
        e->verts[0] = v0;
        e->verts[1] = v1;
        v0->edges.push_back(e);
        v1->edges.push_back(e);
        insertAt(slot, key, e);
        return e;
    }
    
    // k = 0, 1, or 2
    Eptr getTriangleEdge(GTptr gt, uint k, Tptr big_tri)
    {
        GVptr   gv0             = gt->verts[(k+1)%3];
        GVptr   gv1             = gt->verts[(k+2)%3];
        Vptr    v0              = gv0->concrete;
        Vptr    v1              = gv1->concrete;
            // if neither of these are intersection points,
            // then this is a pre-existing edge...
        // NOTE: upstream tested `typeid(gv0) == typeid(OVptr)`, which compares
        // the *static* pointer types and is therefore always false, so every
        // boundary edge was re-created through the cache (a duplicate edge
        // between the same two vertices).  Output verts/tris are unaffected
        // either way; we simply reuse the original edge when both endpoints
        // are corners of the big triangle.  Falls through to the cache if
        // the corner pair is somehow not one of its edges.
        for(uint c=0; c<3; c++) {
            Vptr corner0 = big_tri->verts[(c+1)%3];
            Vptr corner1 = big_tri->verts[(c+2)%3];
            if((corner0 == v0 && corner1 == v1) ||
               (corner0 == v1 && corner1 == v0)) {
                return big_tri->edges[c];
            }
        }
        (void)gv0; (void)gv1;
        return operator()(v0, v1);
    }
    
    Eptr maybeEdge(GEptr ge)
    {
        uint i = ge->ends[0]->concrete->ref;
        uint j = ge->ends[1]->concrete->ref;
        if(i > j) std::swap(i,j);
        uint64_t key = ((uint64_t)i << 32) | (uint64_t)j;
        size_t slot = find(key);
        return (keys[slot] == key) ? vals[slot] : nullptr;
    }
    
private:
    enum : uint64_t { EMPTY = ~(uint64_t)0 };
    static inline size_t hash(uint64_t k) {
        k ^= k >> 33; k *= 0xff51afd7ed558ccdULL;
        k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL;
        k ^= k >> 33;
        return (size_t)k;
    }
    inline size_t find(uint64_t key) const {
        size_t s = hash(key) & mask;
        while(keys[s] != EMPTY && keys[s] != key) s = (s + 1) & mask;
        return s;
    }
    void insertAt(size_t slot, uint64_t key, Eptr e) {
        keys[slot] = key;
        vals[slot] = e;
        if(++count * 2 > keys.size()) grow();
    }
    void grow() {
        std::vector<uint64_t> ok(std::move(keys));
        std::vector<Eptr>     ov(std::move(vals));
        size_t cap = ok.size() * 2;
        keys.assign(cap, (uint64_t)EMPTY);
        vals.assign(cap, nullptr);
        mask = cap - 1;
        for(size_t s=0; s<ok.size(); s++) {
            if(ok[s] == EMPTY) continue;
            size_t t = find(ok[s]);
            keys[t] = ok[s];
            vals[t] = ov[s];
        }
    }
    
    IsctProblem *iprob;
    std::vector<uint64_t>   keys;
    std::vector<Eptr>       vals;
    size_t                  mask  = 0;
    size_t                  count = 0;
};

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::createRealPtFromGluePt(GluePt glue) {
    ENSURE(glue->copies.size() > 0);
    Vptr        v               = TopoCache::newVert();
    VertData    &data           = TopoCache::mesh->verts[v->ref];
                data.pos        = glue->copies[0]->coord;
                fillOutVertData(glue, data);
    for(IVptr iv : glue->copies)
                iv->concrete    = v;
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::createRealTriangles(
    Tprob tprob, EdgeCache &ecache
) {
    using clk = std::chrono::steady_clock;
    const bool prof = cork_prof::enabled();
    for(GTptr gt : tprob->gtris) {
        clk::time_point a0; if(prof) a0 = clk::now();
        Tptr        t               = TopoCache::newTri();
                    gt->concrete    = t;
        Tri         &tri            = TopoCache::mesh->tris[t->ref];
        for(uint k=0; k<3; k++) {
            Vptr    v               = gt->verts[k]->concrete;
                    t->verts[k]     = v;
                    v->tris.push_back(t);
                    tri.v[k]        = v->ref;
        }
        clk::time_point a1; if(prof) { a1 = clk::now(); crt_ns[0] += (a1-a0).count(); }
        // Edge topology for the new pieces is intentionally NOT built.
        // The IsctProblem's TopoCache is discarded right after commit()
        // (which only reads verts/tris), and boolean operations rebuild
        // their own edge cache from the committed mesh, so the ~1.2M hash
        // lookups + edge allocations here never influenced any output.
        for(uint k=0; k<3; k++) t->edges[k] = nullptr;
        (void)ecache;
        clk::time_point a2; if(prof) { a2 = clk::now(); crt_ns[1] += (a2-a1).count(); }
                    fillOutTriData(t, tprob->the_tri);
        if(prof) { crt_ns[2] += (clk::now()-a2).count(); }
    }
    // Once all the pieces are hooked up, let's kill the old triangle!
    clk::time_point d0; if(prof) d0 = clk::now();
    TopoCache::deleteTri(tprob->the_tri);
    if(prof) crt_ns[3] += (clk::now()-d0).count();
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::resolveAllIntersections()
{
    // solve a subdivision problem in each triangle
    std::vector<Tprob> probs;
    tprobs.collect(probs);
    size_t n_new_tris = 0;
    {
        std::vector<typename TriangleProblem::SubdivData> sd(probs.size());
        {
            CORK_PROF("    subdivide: prepare");
            for(size_t i=0; i<probs.size(); i++)
                probs[i]->subdivide_prepare(this, sd[i]);
        }
        {
            CORK_PROF("    subdivide: triangulate (par)");
            cork_par::for_each_idx(probs.size(), 16, [&](size_t i) {
                TriangleProblem::subdivide_triangulate(sd[i]);
            });
        }
        {
            CORK_PROF("    subdivide: finish");
            for(size_t i=0; i<probs.size(); i++) {
                probs[i]->subdivide_finish(this, sd[i]);
                n_new_tris += sd[i].tris.size() / 3;
            }
        }
    }
    
    // now we have diced up triangles inside each triangle problem
    
    // Let's go through the glue points and create a new concrete
    // vertex object for each of these.
    {
        CORK_PROF("    createRealPtFromGluePt");
        TopoCache::mesh->verts.reserve(TopoCache::mesh->verts.size() + glue_pts.size());
        glue_pts.for_each([&](GluePt glue) {
            createRealPtFromGluePt(glue);
        });
    }
    
    EdgeCache ecache(this, 16);   // unused now, see createRealTriangles
    
    // Now that we have concrete vertices plugged in, we can
    // go through the diced triangle pieces and create concrete triangles
    // for each of those.
    // Along the way, let's go ahead and hook up edges as appropriate
    {
        CORK_PROF("    createRealTriangles");
        // One bulk allocation + one mesh.tris resize; old triangles are
        // simply unlinked (commit only cares about live pool objects).
        const size_t oldN = TopoCache::mesh->tris.size();
        TopoCache::mesh->tris.resize(oldN + n_new_tris);
        auto nb = TopoCache::tris.alloc_bulk((uint)n_new_tris);
        size_t w = 0;
        for(Tprob tprob : probs) {
            for(GTptr gt : tprob->gtris) {
                Tptr t = nb[w];
                t->ref = (uint)(oldN + w);
                t->data = nullptr;
                gt->concrete = t;
                Tri &tri = TopoCache::mesh->tris[t->ref];
                for(uint k=0; k<3; k++) {
                    Vptr v = gt->verts[k]->concrete;
                    t->verts[k] = v;
                    t->edges[k] = nullptr;
                    tri.v[k] = v->ref;
                }
                // resolveIntersection does not consume bool_alg_data
                ++w;
            }
            TopoCache::freeTri(tprob->the_tri);
        }
        (void)ecache;
        cork_prof::note("    crt: ns newTri+hook", (double)crt_ns[0]);
        cork_prof::note("    crt: ns edges", (double)crt_ns[1]);
        cork_prof::note("    crt: ns fillOutTriData", (double)crt_ns[2]);
        cork_prof::note("    crt: ns deleteTri", (double)crt_ns[3]);
    }
    
    // "mark isct edges" (e->data = 1 on intersection edges) used to run here.
    // Nothing ever read those marks: the IsctProblem (and its TopoCache) is
    // destroyed immediately after commit(), and BoolProblem builds its own
    // edge cache from the committed mesh.  Dropped together with the edge
    // hookup in createRealTriangles.
    
    // This basically takes care of everything EXCEPT one detail
    // *) The base mesh data structures still need to be compacted
    
    // This detail should be handled by the calling code...
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::dumpIsctPoints(
    std::vector<Vec3d> *points
) {
    points->resize(glue_pts.size());
    uint write = 0;
    glue_pts.for_each([&](GluePt glue) {
        ENSURE(glue->copies.size() > 0);
        IVptr       iv                  = glue->copies[0];
                    (*points)[write]    = iv->coord;
                    write++;
    });
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::IsctProblem::dumpIsctEdges(
    std::vector< std::pair<Vec3d,Vec3d> > *edges
) {
    edges->clear();
    tprobs.for_each([&](Tprob tprob) {
        for(IEptr ie : tprob->iedges) {
            GVptr gv0 = ie->ends[0];
            GVptr gv1 = ie->ends[1];
            edges->push_back(std::make_pair(gv0->coord, gv1->coord));
        }
    });
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::testingComputeStaticIsctPoints(
    std::vector<Vec3d> *points
) {
    IsctProblem iproblem(this);
    
    iproblem.findIntersections();
    
    iproblem.dumpIsctPoints(points);
}


template<class VertData, class TriData>
void Mesh<VertData,TriData>::testingComputeStaticIsct(
    std::vector<Vec3d> *points,
    std::vector< std::pair<Vec3d,Vec3d> > *edges
) {
    IsctProblem iproblem(this);
    
    iproblem.findIntersections();
    
    iproblem.dumpIsctPoints(points);
    iproblem.dumpIsctEdges(edges);
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::resolveIntersections()
{
    CORK_PROF("resolveIntersections total");
    IsctProblem *ip = nullptr;
    {
        CORK_PROF("  IsctProblem ctor (topo+quant)");
        ip = new IsctProblem(this);
    }
    IsctProblem &iproblem = *ip;
    {
        CORK_PROF("  findIntersections");
        iproblem.findIntersections();
    }
    {
        CORK_PROF("  resolveAllIntersections");
        iproblem.resolveAllIntersections();
    }
    {
        CORK_PROF("  commit");
        iproblem.commit();
    }
    {
        CORK_PROF("  ~IsctProblem");
        delete ip;
    }
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::isSelfIntersecting()
{
    IsctProblem iproblem(this);
    
    return iproblem.hasIntersections();
}



