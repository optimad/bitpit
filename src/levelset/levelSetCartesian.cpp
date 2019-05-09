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

# include "bitpit_volcartesian.hpp"
# include "bitpit_CG.hpp"

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
    double diagonalSquare = 0.0;

    for(int i=0; i<dim; ++i){
        diagonalSquare += spacing[i]*spacing[i];
    }

    return 0.5*sqrt(diagonalSquare);
}

/*!
 * Checks if a plane intersects the cell
 * @param[in] id is the index of cell
 * @param[in] root is a point on the plane
 * @param[in] normal is the normal of the plane
 * @return true if intersect
 */
bool LevelSetCartesian::intersectCellPlane( long id, const std::array<double,3> &root, const std::array<double,3> &normal ) {

    std::array<double,3> centroid( computeCellCentroid(id) );
    std::array<double,3> spacing( m_cartesian->getSpacing() );
    std::array<double,3> minPoint( centroid -0.5*spacing );
    std::array<double,3> maxPoint( centroid +0.5*spacing );

    int dim = m_cartesian->getDimension();
    return CGElem::intersectPlaneBox( root, normal, minPoint, maxPoint, dim);
}

}
