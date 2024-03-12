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
 *  @class      LevelSetUnstructuredKernel
 *  @brief      Implements LevelSetKernel for unstructured meshes
 */

/*!
 * Constructor
 *
 * \param patch is the underlying mesh
 * \param fillIn expected kernel fill-in
 */
LevelSetUnstructuredKernel::LevelSetUnstructuredKernel(VolUnstructured &patch, LevelSetFillIn fillIn ): LevelSetCachedKernel( &patch, fillIn ){

    CellCacheCollection &cacheCollection = getCellCacheCollection();
    if (fillIn == LevelSetFillIn::SPARSE) {
        m_cellCentroidCacheId       = cacheCollection.insert<CellSparseCacheContainer<std::array<double, 3>>>(CellCacheCollection::NULL_CACHE_ID);
        m_cellTangentRadiusCacheId  = cacheCollection.insert<CellSparseCacheContainer<double>>(CellCacheCollection::NULL_CACHE_ID);
        m_cellBoundingRadiusCacheId = cacheCollection.insert<CellSparseCacheContainer<double>>(CellCacheCollection::NULL_CACHE_ID);
    } else if (fillIn == LevelSetFillIn::DENSE) {
        m_cellCentroidCacheId       = cacheCollection.insert<CellDenseCacheContainer<std::array<double, 3>>>(CellCacheCollection::NULL_CACHE_ID);
        m_cellTangentRadiusCacheId  = cacheCollection.insert<CellDenseCacheContainer<double>>(CellCacheCollection::NULL_CACHE_ID);
        m_cellBoundingRadiusCacheId = cacheCollection.insert<CellDenseCacheContainer<double>>(CellCacheCollection::NULL_CACHE_ID);
    } else {
        m_cellCentroidCacheId       = CellCacheCollection::NULL_CACHE_ID;
        m_cellTangentRadiusCacheId  = CellCacheCollection::NULL_CACHE_ID;
        m_cellBoundingRadiusCacheId = CellCacheCollection::NULL_CACHE_ID;
    }
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
std::array<double, 3> LevelSetUnstructuredKernel::computeCellCentroid( long id ) const {

    // Try fetching the value from the cache
    CellCacheCollection::ValueCache<std::array<double, 3>> *cache = (*m_cellCacheCollection)[m_cellCentroidCacheId].getCache<std::array<double, 3>>();
    if (cache) {
        typename CellCacheCollection::ValueCache<std::array<double, 3>>::Entry cacheEntry = cache->findEntry(id);
        if (cacheEntry.isValid()) {
            return *cacheEntry;
        }
    }

    // Evaluate the centroid
    std::array<double, 3> centroid = getMesh()->evalCellCentroid(id);

    // Update the cache
    if (cache) {
        cache->insertEntry(id, centroid);
    }

    return centroid;

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

    // Try fetching the value from the cache
    CellCacheCollection::ValueCache<double> *cache = (*m_cellCacheCollection)[m_cellTangentRadiusCacheId].getCache<double>();
    if (cache) {
        typename CellCacheCollection::ValueCache<double>::Entry cacheEntry = cache->findEntry(id);
        if (cacheEntry.isValid()) {
            return *cacheEntry;
        }
    }

    // Get mesh information
    const VolUnstructured &unstructuredPatch = *(this->getMesh()) ;

    // Cell information
    const Cell &cell = unstructuredPatch.getCell(id);
    ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
    int nCellVertices = cellVertexIds.size();
    std::array<double, 3> centroid = computeCellCentroid(id);

    BITPIT_CREATE_WORKSPACE(vertexCoordinates, std::array<double BITPIT_COMMA 3>, nCellVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    unstructuredPatch.getVertexCoords(nCellVertices, cellVertexIds.data(), vertexCoordinates);

    // Evaluate the radius of the tangent sphere
    //
    // Since the center of the tangent sphere is the cell centroid, its
    // radius can be evaluated as distance between the cell centroid and
    // the cell.
    double tangentRadius = cell.evalPointDistance(centroid, vertexCoordinates);

    // Update the cache
    if (cache) {
        cache->insertEntry(id, tangentRadius);
    }

    return tangentRadius;

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

    // Try fetching the value from the cache
    CellCacheCollection::ValueCache<double> *cache = (*m_cellCacheCollection)[m_cellBoundingRadiusCacheId].getCache<double>();
    if (cache) {
        typename CellCacheCollection::ValueCache<double>::Entry cacheEntry = cache->findEntry(id);
        if (cacheEntry.isValid()) {
            return *cacheEntry;
        }
    }

    // Get mesh information
    const VolUnstructured &unstructuredPatch = *(this->getMesh()) ;

    // Cell information
    const Cell &cell = unstructuredPatch.getCell(id);
    ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
    int nCellVertices = cellVertexIds.size();
    std::array<double, 3> centroid = computeCellCentroid(id);

    BITPIT_CREATE_WORKSPACE(vertexCoordinates, std::array<double BITPIT_COMMA 3>, nCellVertices, ReferenceElementInfo::MAX_ELEM_VERTICES);
    unstructuredPatch.getVertexCoords(nCellVertices, cellVertexIds.data(), vertexCoordinates);

    // Evaluate the radius of the bounding sphere
    //
    // Since the center of the bounding sphere is the cell centroid, its
    // radius can be evaluated as distance between the cell centroid and
    // the farthest vertex.
    double boundingRadius = 0.;
    for (int k = 0; k < nCellVertices; ++k) {
        double vertexDistance = 0.;
        for (int d = 0; d < 3; ++d) {
            vertexDistance += uipow(centroid[d] - vertexCoordinates[k][d], 2);
        }

        boundingRadius = std::max(vertexDistance, boundingRadius);
    }
    boundingRadius = std::sqrt(boundingRadius);

    // Update the cache
    if (cache) {
        cache->insertEntry(id, boundingRadius);
    }

    return boundingRadius;

}

}
