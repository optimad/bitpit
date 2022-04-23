/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2022 OPTIMAD engineering Srl
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

# include "levelSetUnstructuredKernel.hpp"

# include "bitpit_CG.hpp"

namespace bitpit {

/*!
 *  @ingroup    levelset
 *  @class      LevelSetUnstructuredCellCacheEntry
 *  @brief      Cache entry associated with cells of unstructured patches, the entry will
 *              store cell information needed to evaluate the levelset.
 */

/*!
 * Constructor
 *
 * @param[in] patch is the patch the cell belongs to
 * @param[in] cellId is the id of the cell
 */
LevelSetUnstructuredCellCacheEntry::LevelSetUnstructuredCellCacheEntry(const VolumeKernel &patch, long cellId) {

    // Get mesh information
    assert(dynamic_cast<const VolUnstructured *>(&patch)) ;
    const VolUnstructured &unstructuredPatch = static_cast<const VolUnstructured &>(patch) ;

    // Cell information
    const Cell &cell = unstructuredPatch.getCell(cellId);
    ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
    int nCellVertices = cellVertexIds.size();

    BITPIT_CREATE_WORKSPACE(vertexCoordinates, std::array<double BITPIT_COMMA 3>, nCellVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    unstructuredPatch.getVertexCoords(nCellVertices, cellVertexIds.data(), vertexCoordinates);

    // Evaluate centroid
    centroid = cell.evalCentroid(vertexCoordinates);

    // Evaluate the radius of the tangent sphere
    //
    // Since the center of the tangent sphere is the cell centroid, its
    // radius can be evaluated as distance between the cell centroid and
    // the cell.
    tangentRadius = cell.evalPointDistance(centroid, vertexCoordinates);

    // Evaluate the radius of the bounding sphere
    //
    // Since the center of the bounding sphere is the cell centroid, its
    // radius can be evaluated as distance between the cell centroid and
    // the farthest vertex.
    boundingRadius = 0.;
    for (int k = 0; k < nCellVertices; ++k) {
        double vertexDistance = 0.;
        for (int d = 0; d < 3; ++d) {
            vertexDistance += uipow(centroid[d] - vertexCoordinates[k][d], 2);
        }

        boundingRadius = std::max(vertexDistance, boundingRadius);
    }
    boundingRadius = std::sqrt(boundingRadius);

}

/*!
 *  @ingroup    levelset
 *  @class      LevelSetUnstructuredKernel
 *  @brief      Implements LevelSetKernel for unstructured meshes
 */

/*!
 * Constructor
 */
LevelSetUnstructuredKernel::LevelSetUnstructuredKernel(VolUnstructured & patch ): LevelSetCachedKernel<LevelSetUnstructuredCellCacheEntry>( &patch ){

}

/*!
 * Returns a pointer to VolUnstructured
 * @return pointer to VolUnstructured
 */
VolUnstructured * LevelSetUnstructuredKernel::getMesh() const{
    return static_cast<VolUnstructured *>(LevelSetKernel::getMesh()) ;
}

/*!
 * Computes the centroid of the specfified cell.
 *
 * @param[in] id is the id of cell
 * @return The centroid of the cell.
 */
std::array<double, 3> LevelSetUnstructuredKernel::computeCellCentroid(long id) const {

    const LevelSetUnstructuredCellCacheEntry &cacheEntry = computeCellCacheEntry( id );

    return cacheEntry.centroid;

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
double LevelSetUnstructuredKernel::computeCellTangentRadius(long id) const {

    const LevelSetUnstructuredCellCacheEntry &cacheEntry = computeCellCacheEntry( id );

    return cacheEntry.tangentRadius;

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
double LevelSetUnstructuredKernel::computeCellBoundingRadius( long id ) const {

    const LevelSetUnstructuredCellCacheEntry &cacheEntry = computeCellCacheEntry( id );

    return cacheEntry.boundingRadius;

}

}
