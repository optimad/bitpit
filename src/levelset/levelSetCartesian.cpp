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

VolCartesian* LevelSetCartesian::getCartesianMesh() const{
    return m_cartesian ;
}
/*!
 * Calculate size of narrow band in order to guarantee one element on each side of geometry
 * @param[in] visitor reference object 
 * @return size of narrow band
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
 * @param[in] visitor reference object 
 * @return size of narrow band
 */
double LevelSetCartesian::updateSizeNarrowBand( const std::vector<adaption::Info> &mapper, LevelSetObject *visitor ){

    BITPIT_UNUSED(mapper) ;
    BITPIT_UNUSED(visitor) ;

    return computeRSearchFromCell( 0 ) ;

};

/*!
 * Compute size of narrow band given a cell.
 * @param[in] id is the id of the cell
 * @return size of narrow band
 */
double LevelSetCartesian::computeRSearchFromCell( long id ){

    BITPIT_UNUSED(id) ;

    double newRSearch = 0. ;

    for( int d=0; d<m_cartesian->getDimension(); ++d){
        newRSearch = std::max( newRSearch, m_cartesian->getSpacing(d) ) ;
    };

    return newRSearch;

};

}
