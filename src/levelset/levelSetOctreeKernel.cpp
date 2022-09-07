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
 *  @class      LevelSetOctreeCellCacheEntry
 *  @brief      Cache entry associated with cells of octree patches, the entry will
 *              store cell information needed to evaluate the levelset.
 */

/*!
 * Constructor
 *
 * @param[in] patch is the patch the cell belongs to
 * @param[in] cellId is the id of the cell
 */
LevelSetOctreeCellCacheEntry::LevelSetOctreeCellCacheEntry(const VolumeKernel &patch, long cellId) {

    // Get mesh information
    assert(dynamic_cast<const VolOctree *>(&patch)) ;
    const VolOctree &octreePatch = static_cast<const VolOctree &>(patch) ;

    // Evaluate cell centroid
    centroid = octreePatch.evalCellCentroid(cellId) ;

}

/*!
 *  @ingroup    levelset
 *  @class      LevelSetOctreeKernel
 *  @brief      Implements LevelSetKernel for octree meshes
 */

/*!
 * Constructor
 */
LevelSetOctreeKernel::LevelSetOctreeKernel(VolOctree & patch ): LevelSetCachedKernel<LevelSetOctreeCellCacheEntry>( &patch ){

    // Mesh information
    VolOctree *mesh = getMesh();
    int dimension = mesh->getDimension();
    int maxLevel  = TreeConstants::instance(dimension).maxLevel;

    // Initialize bounding and tangent radii
    m_octantTangentRadii.resize(maxLevel + 1);
    m_octantBoundingRadii.resize(maxLevel + 1);
    for (int level = 0; level <= maxLevel; ++level) {
        double levelSize = mesh->getTree().levelToSize(level);

        m_octantTangentRadii[level]  = 0.5 * levelSize;
        m_octantBoundingRadii[level] = 0.5 * std::sqrt(dimension) * levelSize;
    }

}

/*!
 * Returns a pointer to VolOctree
 * @return pointer to VolOctree
 */
VolOctree * LevelSetOctreeKernel::getMesh() const{
    return static_cast<VolOctree *>(LevelSetKernel::getMesh()) ;
}

/*!
 * Computes the radius of the tangent sphere associated with the specified level.
 *
 * The tangent sphere is a sphere having the center in the level centroid and tangent
 * to the cell.
 *
 * @param[in] level is the level of the octant
 * @return The radius of the tangent sphere.
 */
double LevelSetOctreeKernel::getOctantTangentRadius( int level ) const {

    return m_octantTangentRadii[level];

}

/*!
 * Computes the radius of the bounding sphere associated with the specified level.
 *
 * The bounding sphere is the sphere with the minimum radius that contains all the
 * level vertices and has the center in the level centroid.
 *
 * @param[in] level is the level of the octant
 * @return The radius of the bounding sphere.
 */
double LevelSetOctreeKernel::getOctantBoundingRadius( int level ) const {

    return m_octantBoundingRadii[level];

}

/*!
 * Computes the centroid of the specfified cell.
 *
 * @param[in] id is the id of cell
 * @return The centroid of the cell.
 */
std::array<double, 3> LevelSetOctreeKernel::computeCellCentroid( long id ) const {

    const LevelSetOctreeCellCacheEntry &cellCacheEntry = computeCellCacheEntry( id );

    return cellCacheEntry.centroid;

}

/*!
 * Computes the radius of the tangent sphere associated with the specified cell.
 *
 * The tangent sphere is a sphere having the center in the cell centroid and tangent
 * to the cell.
 *
 * @param[in] id is the id of cell
 * @return The radius of the tangent sphere.
 */
double LevelSetOctreeKernel::computeCellTangentRadius( long id ) const {

    const VolOctree *mesh = getMesh();
    int cellLevel = mesh->getCellLevel(id);

    return getOctantTangentRadius(cellLevel);

}

/*!
 * Computes the radius of the bounding sphere associated with the specified cell.
 *
 * The bounding sphere is the sphere with the minimum radius that contains all the
 * cell vertices and has the center in the cell centroid.
 *
 * @param[in] id is the id of cell
 * @return The radius of the bounding sphere.
 */
double LevelSetOctreeKernel::computeCellBoundingRadius( long id ) const {

    const VolOctree *mesh = getMesh();
    int cellLevel = mesh->getCellLevel(id);

    return getOctantBoundingRadius(cellLevel);

}

}
