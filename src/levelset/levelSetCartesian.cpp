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

# include "bitpit_volcartesian.hpp"

# include "levelSetKernel.hpp"
# include "levelSetCartesian.hpp"


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
}

/*!
 * Constructor
 */
LevelSetCartesian::LevelSetCartesian(VolCartesian &patch ): LevelSetKernel( (static_cast<VolumeKernel*>(&patch)) ){
    m_cartesian = &patch ;
}

/*!
 * Returns a pointer to VolCartesian
 * @return pointer to VolCartesian
 */
VolCartesian* LevelSetCartesian::getCartesianMesh() const{
    return m_cartesian ;
}

/*!
 * Computes the radius of the incircle of the specfified cell.
 * @param[in] id is the index of cell
 * @return radius of incircle
 */
double LevelSetCartesian::computeCellIncircle(long id) {

    BITPIT_UNUSED(id);

    int dim = m_cartesian->getDimension();
    std::array<double,3> spacing = m_cartesian->getSpacing();
    double minSpacing = std::numeric_limits<double>::max() ;

    for(int i=0; i<dim; ++i){
        minSpacing = std::min( minSpacing, spacing[i] );
    }

    return 0.5*minSpacing;
}

/*!
 * Computes the radius of the circumcircle of the specfified cell.
 * @param[in] id is the index of cell
 * @return radius of incircle
 */
double LevelSetCartesian::computeCellCircumcircle( long id ) {

    BITPIT_UNUSED(id);

    int dim = m_cartesian->getDimension();
    std::array<double,3> spacing = m_cartesian->getSpacing();
    double maxSpacing = -std::numeric_limits<double>::max() ;

    for(int i=0; i<dim; ++i){
        maxSpacing = std::max( maxSpacing, spacing[i] );
    }

    return 0.5*sqrt((float) dim)*maxSpacing;
}

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
    }

    return newRSearch;

}

}
