/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

# include "bitpit_voloctree.hpp"
# include "bitpit_CG.hpp"
# include "levelSetKernel.hpp"
# include "levelSetOctree.hpp"

namespace bitpit {

/*!
 *  @ingroup    levelset
 *  @class      LevelSetOctree
 *  @brief      Implements LevelSetKernel for octree meshes
 */

/*!
 * Destructor
 */
LevelSetOctree::~LevelSetOctree( ){
    m_octree = NULL ;
}

/*!
 * Constructor
 */
LevelSetOctree::LevelSetOctree(VolOctree & patch ): LevelSetKernel( (static_cast<VolumeKernel*>(&patch)) ){
    m_octree = &patch ;
}

/*!
 * Returns a pointer to VolOctree
 * @return pointer to VolOctree
 */
VolOctree* LevelSetOctree::getOctreeMesh() const{
    return m_octree ;
}

/*!
 * Computes the radius of the incircle of the specfified cell.
 * @param[in] id is the index of cell
 * @return radius of incircle
 */
double LevelSetOctree::computeCellIncircle(long id) {
    return 0.5*m_octree->evalCellSize(id);
}

/*!
 * Computes the radius of the circumcircle of the specfified cell.
 * @param[in] id is the index of cell
 * @return radius of incircle
 */
double LevelSetOctree::computeCellCircumcircle( long id ) {
    int dim = m_octree->getDimension();
    return 0.5*sqrt((float) dim)*m_octree->evalCellSize(id);
}

/*!
 * Checks if a plane intersects the cell
 * @param[in] id is the index of cell
 * @param[in] root is a point on the plane
 * @param[in] normal is the normal of the plane
 * @return true if intersect
 */
bool LevelSetOctree::intersectCellPlane( long id, const std::array<double,3> &root, const std::array<double,3> &normal ) {

    std::array<double,3> centroid( computeCellCentroid(id) );
    double spacing( m_octree->evalCellSize(id) );
    std::array<double,3> minPoint( centroid -0.5*spacing );
    std::array<double,3> maxPoint( centroid +0.5*spacing );

    int dim = m_octree->getDimension();
    return CGElem::intersectPlaneBox( root, normal, minPoint, maxPoint, dim);
}
}
