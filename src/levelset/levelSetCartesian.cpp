/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

# include "levelSet.hpp"

# include "bitpit_SA.hpp"
# include "bitpit_CG.hpp"

namespace bitpit {

/*!
	@ingroup    levelset
	@class      LevelSetCartesian
	@brief      Implements LevelSetKernel for cartesian meshes
*/

/*!
 * Destructor
 */
LevelSetCartesian::~LevelSetCartesian( ){
    m_cartesian = NULL ;
};

/*!
 * Constructor
 */
LevelSetCartesian::LevelSetCartesian(VolCartesian &patch ): LevelSetKernel( (static_cast<VolumeKernel*>(&patch)) ){
    m_cartesian = &patch ;
};

/*!
 * Calculate size of narrow band in order to guarantee one element on each side of geometry
 * @param[in] visitor reference object 
 */
double LevelSetCartesian::computeSizeNarrowBand( LevelSetObject *visitor ){

    BITPIT_UNUSED(visitor) ;

    double RSearch = -1.;

    for( int d=0; d<m_cartesian->getDimension(); ++d){
        RSearch = std::max( RSearch, m_cartesian->getSpacing(d) ) ;
    };

    return RSearch;
};

/*!
 * Update the size of the narrow band after an adaptation of the cartesian mesh
 * @param[in] mapper mesh modifications
 */
double LevelSetCartesian::updateSizeNarrowBand( const std::vector<adaption::Info> &mapper ){

    BITPIT_UNUSED(mapper) ;

    double newRSearch = -1. ;

    for( int d=0; d<m_cartesian->getDimension(); ++d){
        newRSearch = std::max( newRSearch, m_cartesian->getSpacing(d) ) ;
    };

    return newRSearch;

};

/*! 
 * Update scalar field value at mesh vertex on a by locally solving the 3D Eikonal equation.
 * @param[in] s Flag for inwards/outwards propagation (s = -+1).
 * @param[in] g Propagation speed for the 3D Eikonal equation.
 * @param[in] I index of the cartesian cell to be updated.
 * @return Updated value at mesh vertex
 */
double LevelSetCartesian::updateEikonal( double s, double g, const long &I ){

    int                d;
    long               J;
    double             h2, delta, value, a(0), b(0), c(0);

    Cell&   cell = m_cartesian->getCell(I) ;

    for( d=0; d<m_cartesian->getDimension(); ++d){ // COMPUTE QUADRATIC FORM COEFFICIENTS FROM UPWIND STENCIL

        value   = levelSetDefaults::VALUE ;

        // Left neighbor
        J   = cell.getAdjacency( 2*d, 0) ;

        LSInfo  &lsInfo = m_ls[J] ;

        if( J >= 0 && lsInfo.active == 0){
            value = std::min(s*lsInfo.value, value);
        };


        // Right neighbor
        J   = cell.getAdjacency( 2*d+1, 0) ;

        if( J >= 0 && lsInfo.active == 0){
            value = std::min(s*lsInfo.value, value);
        };


        // Update coeffs in the quadratic form
        if (value < levelSetDefaults::VALUE) {
            h2 = pow(m_cartesian->getSpacing(d), 2);

            a += 1.0/h2;
            b += -2.0 * value/h2;
            c += std::pow(value, 2)/h2;
        }

    };


    { // SOLVE THE QUADRATIC FORM
        // Quadratic form determinant
        delta = pow(b, 2) - 4.0*a*(c - pow(g, 2));

        // Solution
        value = -(b - sqrt(delta))/(2.0*a);

    }

    return(value); 

};

}
