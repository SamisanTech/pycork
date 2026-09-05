// +-------------------------------------------------------------------------
// | mesh.hull_exact.tpp
// |
// | Opt-in Manifold-style outer hull.  Default outerHull() is unchanged.
// |
// | Manifold's robust hull is not a separate extractor: after SI / self-union
// | it keeps the regularized |winding| >= 1 solid, with winding from
// | PointWinding / Kernel02 — axis-aligned +Z signed face crossings
// | (src/boolean3.cpp, src/manifold.cpp).  Random rays are the fragile part
// | of our hull; this path uses the same patch / seed / leftover machinery
// | with +Z (then +X/+Y) only.
// +-------------------------------------------------------------------------
#pragma once

template<class VertData, class TriData>
void Mesh<VertData,TriData>::outerHullExact(cork_hull::HullStats *stats,
                                            double leftoverAreaFrac)
{
    outerHull(5, stats, leftoverAreaFrac, /*exact=*/true);
}
