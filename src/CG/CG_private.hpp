/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

# ifndef __BITPIT_CG_PRIVATE_HPP__
# define __BITPIT_CG_PRIVATE_HPP__

# include <array>
# include <vector>

namespace bitpit{

namespace CGElem{

typedef std::array<double,3> array3D ;

void _projectPointsTriangle( int, array3D const *, array3D const &, array3D const &, array3D const &, array3D *, double *);
void _projectPointsPlane( int, array3D const *, array3D const &, array3D const &, array3D const &, array3D *, double *);
bool _intersectSegmentBox( array3D const &, array3D const &, array3D const &, array3D const &, bool, bool, std::vector<array3D> *, std::vector<int> *, int  dim = 3 ) ;
bool _intersectPlaneBox( array3D const &, array3D const &, array3D const &, array3D const &, std::vector<array3D> *, int dim=3);
bool _intersectBoxTriangle( array3D const &, array3D const &, array3D const &, array3D const &, array3D const &, bool, bool, bool, std::vector<array3D> *, std::vector<int> *, int dim=3 ) ;
bool _intersectBoxPolygon( array3D const &, array3D const &, std::size_t, array3D const *, bool, bool, bool, std::vector<array3D> *, std::vector<int> *, int dim=3 );
BITPIT_DEPRECATED( bool _intersectBoxSimplex( array3D const &, array3D const &, std::vector<array3D> const &, bool, bool, bool, std::vector<array3D> *, std::vector<int> *, int dim=3 ) );

}

}

# endif
