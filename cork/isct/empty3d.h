// +-------------------------------------------------------------------------
// | empty3d.h
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

#ifndef CORK_EMPTY3D_H_HEADER_HAS_BEEN_INCLUDED
#define CORK_EMPTY3D_H_HEADER_HAS_BEEN_INCLUDED

#include <cork/math/vec.h>
#include <CORK_EXPORT.h>

namespace empty3d {

CORK_EXPORT struct TriIn
{
    Vec3d p[3];
};

CORK_EXPORT struct EdgeIn
{
    Vec3d p[2];
};

CORK_EXPORT struct TriEdgeIn
{
    TriIn   tri;
    EdgeIn  edge;
};

CORK_EXPORT bool isEmpty(const TriEdgeIn &input);
CORK_EXPORT Vec3d coords(const TriEdgeIn &input);
CORK_EXPORT bool emptyExact(const TriEdgeIn &input);
CORK_EXPORT Vec3d coordsExact(const TriEdgeIn &input);

CORK_EXPORT struct TriTriTriIn
{
    TriIn tri[3];
};

CORK_EXPORT bool isEmpty(const TriTriTriIn &input);
CORK_EXPORT Vec3d coords(const TriTriTriIn &input);
CORK_EXPORT bool emptyExact(const TriTriTriIn &input);
CORK_EXPORT Vec3d coordsExact(const TriTriTriIn &input);

CORK_EXPORT extern int degeneracy_count; // count degeneracies encountered
CORK_EXPORT extern int exact_count; // count of filter calls failed
CORK_EXPORT extern int callcount; // total call count

/*
// exact versions


bool emptyExact(const Cell3d0 &c0,
                const Complex3d2 &complex,
                const Metric3d2 &metric);

void cell3d0toPointExact(SmVector3 &point,
                         const Cell3d0 &c0,
                         const Complex3d2 &complex,
                         const Metric3d2 &metric);
*/
} // end namespace empty3d

#endif

