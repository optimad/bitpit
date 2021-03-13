/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
 * Constructor
 */
LevelSetCartesian::LevelSetCartesian(VolCartesian &patch ): LevelSetKernel( (static_cast<VolumeKernel*>(&patch)) ){
    m_cartesian = &patch ;

    clearCellCirclesCache();
    updateCellCirclesCache();
}

/*!
 * Returns a pointer to VolCartesian
 * @return pointer to VolCartesian
 */
VolCartesian* LevelSetCartesian::getCartesianMesh() const{
    return m_cartesian ;
}

/*!
 * Get the radius of the incircle of the cells.
 * @return radius of incircle
 */
double LevelSetCartesian::getCellIncircle() {

    return m_cellIncircle;
}

/*!
 * Computes the radius of the incircle of the specfified cell.
 * @param[in] id is the index of cell
 * @return radius of incircle
 */
double LevelSetCartesian::computeCellIncircle(long id) {

    BITPIT_UNUSED(id);

    return getCellIncircle();
}

/*!
 * Get the radius of the circumcircle of the cells.
 * @return radius of circumcircle
 */
double LevelSetCartesian::getCellCircumcircle() {

    return m_cellCircumcircle;
}

/*!
 * Computes the radius of the circumcircle of the specfified cell.
 * @param[in] id is the index of cell
 * @return radius of incircle
 */
double LevelSetCartesian::computeCellCircumcircle( long id ) {

    BITPIT_UNUSED(id);

    return getCellCircumcircle();
}

/*!
 * Checks if a plane intersects the cell
 * @param[in] id is the index of cell
 * @param[in] root is a point on the plane
 * @param[in] normal is the normal of the plane
 * @param[in] tolerance is the tolerance used for distance comparisons
 * @return true if intersect
 */
bool LevelSetCartesian::intersectCellPlane( long id, const std::array<double,3> &root, const std::array<double,3> &normal, double tolerance ) {

    std::array<double,3> minPoint;
    std::array<double,3> maxPoint;
    m_cartesian->evalCellBoundingBox(id, &minPoint, &maxPoint);

    int dim = m_cartesian->getDimension();
    return CGElem::intersectPlaneBox( root, normal, minPoint, maxPoint, dim, tolerance);
}

/*!
 * Clears the geometry cache.
 */
void LevelSetCartesian::clearGeometryCache(  ) {

    LevelSetKernel::clearGeometryCache();

    clearCellCirclesCache();

}

/*!
 * Updates the geometry cache after an adaption.
 */
void LevelSetCartesian::updateGeometryCache( const std::vector<adaption::Info> &mapper ) {

    LevelSetKernel::updateGeometryCache(mapper);

    updateCellCirclesCache();

}

/*!
 * Clears the cache that hold information about cell incircle and circumcircle.
 */
void LevelSetCartesian::clearCellCirclesCache(  ) {

    m_cellIncircle     = 0.;
    m_cellCircumcircle = 0.;

}

/*!
 * Updates the cache that hold information about cell incircle and circumcircle.
 */
void LevelSetCartesian::updateCellCirclesCache(  ) {

    int dimension = m_cartesian->getDimension();
    std::array<double,3> spacing = m_cartesian->getSpacing();

    m_cellIncircle     = std::numeric_limits<double>::max();
    m_cellCircumcircle = 0.;
    for(int i=0; i<dimension; ++i){
        m_cellIncircle      = std::min(m_cellIncircle, spacing[i]);
        m_cellCircumcircle += spacing[i] * spacing[i];
    }
    m_cellIncircle     = 0.5 * m_cellIncircle;
    m_cellCircumcircle = 0.5 * std::sqrt(m_cellCircumcircle);

}
}
