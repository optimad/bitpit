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
# include "levelSetCartesianKernel.hpp"

namespace bitpit {

/*!
	@ingroup    levelset
	@class      LevelSetCartesianKernel
	@brief      Implements LevelSetKernel for cartesian meshes
*/

/*!
 * Constructor
 */
LevelSetCartesianKernel::LevelSetCartesianKernel(VolCartesian &patch ): LevelSetKernel( (static_cast<VolumeKernel*>(&patch)) ){

    clearCellCirclesCache();
    updateCellCirclesCache();
}

/*!
 * Returns a pointer to VolCartesian
 * @return pointer to VolCartesian
 */
VolCartesian * LevelSetCartesianKernel::getMesh() const{
    return static_cast<VolCartesian *>(LevelSetKernel::getMesh()) ;
}

/*!
 * Get the radius of the incircle of the cells.
 * @return radius of incircle
 */
double LevelSetCartesianKernel::getCellIncircle() const {

    return m_cellIncircle;
}

/*!
 * Computes the radius of the incircle of the specfified cell.
 * @param[in] id is the index of cell
 * @return radius of incircle
 */
double LevelSetCartesianKernel::computeCellIncircle(long id) const {

    BITPIT_UNUSED(id);

    return getCellIncircle();
}

/*!
 * Get the radius of the circumcircle of the cells.
 * @return radius of circumcircle
 */
double LevelSetCartesianKernel::getCellCircumcircle() const {

    return m_cellCircumcircle;
}

/*!
 * Computes the radius of the circumcircle of the specfified cell.
 * @param[in] id is the index of cell
 * @return radius of incircle
 */
double LevelSetCartesianKernel::computeCellCircumcircle( long id ) const {

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
bool LevelSetCartesianKernel::intersectCellPlane( long id, const std::array<double,3> &root, const std::array<double,3> &normal, double tolerance ) {

    std::array<double,3> minPoint;
    std::array<double,3> maxPoint;
    m_cartesian->evalCellBoundingBox(id, &minPoint, &maxPoint);

    int dim = m_cartesian->getDimension();
    return CGElem::intersectPlaneBox( root, normal, minPoint, maxPoint, dim, tolerance);
}

/*!
 * Clears the geometry cache.
 */
void LevelSetCartesianKernel::clearGeometryCache(  ) {

    LevelSetKernel::clearGeometryCache();

    clearCellCirclesCache();

}

/*!
 * Updates the geometry cache after an adaption.
 */
void LevelSetCartesianKernel::updateGeometryCache( const std::vector<adaption::Info> &adaptionData ) {

    LevelSetKernel::updateGeometryCache(adaptionData);

    updateCellCirclesCache();

}

/*!
 * Clears the cache that hold information about cell incircle and circumcircle.
 */
void LevelSetCartesianKernel::clearCellCirclesCache(  ) {

    m_cellIncircle     = 0.;
    m_cellCircumcircle = 0.;

}

/*!
 * Updates the cache that hold information about cell incircle and circumcircle.
 */
void LevelSetCartesianKernel::updateCellCirclesCache(  ) {

    const VolCartesian *mesh = getMesh();
    int dimension = mesh->getDimension();
    std::array<double,3> spacing = mesh->getSpacing();

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
