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

# include "bitpit_voloctree.hpp"
# include "bitpit_CG.hpp"
# include "levelSetKernel.hpp"
# include "levelSetOctreeKernel.hpp"

namespace bitpit {

/*!
 *  @ingroup    levelset
 *  @class      LevelSetOctreeKernel
 *  @brief      Implements LevelSetKernel for octree meshes
 */

/*!
 * Constructor
 */
LevelSetOctreeKernel::LevelSetOctreeKernel(VolOctree & patch ): LevelSetKernel( &patch ){
    clearCellCirclesCache();
    updateCellCirclesCache();
}

/*!
 * Returns a pointer to VolOctree
 * @return pointer to VolOctree
 */
VolOctree * LevelSetOctreeKernel::getMesh() const{
    return static_cast<VolOctree *>(LevelSetKernel::getMesh()) ;
}

/*!
 * Computes the radius of the incircle of the specfified cell.
 * @param[in] id is the index of cell
 * @return radius of incircle
 */
double LevelSetOctreeKernel::computeCellIncircle(long id) const {
    const VolOctree *mesh = getMesh();
    int cellLevel = mesh->getCellLevel(id);
    return m_levelToCellIncircle[cellLevel];
}

/*!
 * Computes the radius of the circumcircle of the specfified cell.
 * @param[in] id is the index of cell
 * @return radius of incircle
 */
double LevelSetOctreeKernel::computeCellCircumcircle( long id ) const {
    const VolOctree *mesh = getMesh();
    int cellLevel = mesh->getCellLevel(id);
    return m_levelToCellCircumcircle[cellLevel];
}

/*!
 * Clears the geometry cache.
 */
void LevelSetOctreeKernel::clearGeometryCache(  ) {

    LevelSetKernel::clearGeometryCache();

    clearCellCirclesCache();

}

/*!
 * Updates the geometry cache after an adaption.
 */
void LevelSetOctreeKernel::updateGeometryCache( const std::vector<adaption::Info> &adaptionData ) {

    LevelSetKernel::updateGeometryCache(adaptionData);

    updateCellCirclesCache();

}

/*!
 * Clears the cache that hold information about cell incircle and circumcircle.
 */
void LevelSetOctreeKernel::clearCellCirclesCache(  ) {

    m_levelToCellIncircle.clear();
    m_levelToCellCircumcircle.clear();

}

/*!
 * Updates the cache that hold information about cell incircle and circumcircle.
 */
void LevelSetOctreeKernel::updateCellCirclesCache(  ) {

    VolOctree *mesh = getMesh();
    int dimension = mesh->getDimension();
    int maxLevel  = mesh->getTree().getMaxLevel();

    m_levelToCellIncircle.resize(maxLevel + 1);
    m_levelToCellCircumcircle.resize(maxLevel + 1);
    for (int level = 0; level <= maxLevel; ++level) {
        double levelSize = mesh->getTree().levelToSize(level);

        m_levelToCellIncircle[level]     = 0.5 * levelSize;
        m_levelToCellCircumcircle[level] = 0.5 * std::sqrt(dimension) * levelSize;
    }

}
}
