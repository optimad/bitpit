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

# ifndef __BITPIT_LEVELSET_OBJECT_TPP__
# define __BITPIT_LEVELSET_OBJECT_TPP__

namespace bitpit {

/*!
 * Register the cache that will be used for storing cell information of the specified field.
 *
 * \param field is the field for which the caches will be added
 * \param cacheMode is the type of cache that will be used for caching cell information
 * \result The id associated with the registered cache.
 */
template<typename value_t>
std::size_t LevelSetObject::registerFieldCellCache(LevelSetField field, LevelSetCacheMode cacheMode)
{
    // Update the collection
    std::size_t cacheId = registerCellCache<value_t>(cacheMode);

    // Set cache information
    std::size_t fieldIndex = static_cast<std::size_t>(field);
    if (cacheId != CellCacheCollection::NULL_ID && cacheMode != LevelSetCacheMode::NONE) {
        m_cellFieldCacheIds[fieldIndex]   = cacheId;
        m_cellFieldCacheModes[fieldIndex] = cacheMode;
    } else {
        m_cellFieldCacheIds[fieldIndex]   = CellCacheCollection::NULL_ID;
        m_cellFieldCacheModes[fieldIndex] = LevelSetCacheMode::NONE;
    }

    return cacheId;
}

/*!
 * Register the cache that will be used for storing cell information of the specified field.
 *
 * \param cacheMode is the type of cache that will be used for caching cell information
 * \result The id associated with the registered cache.
 */
template<typename value_t>
std::size_t LevelSetObject::registerCellCache(LevelSetCacheMode cacheMode)
{
    LevelSetFillIn fillIn = m_kernel->getExpectedFillIn();

    PiercedVector<Cell, long> *piercedCacheKernel = &(m_kernel->getMesh()->getCells());
    PiercedSyncMaster::SyncMode piercedCacheSyncMode = PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED;

    std::size_t cacheId = CellCacheCollection::NULL_ID;
    if (dynamic_cast<const LevelSetCartesianKernel *>(m_kernel)){
        if (fillIn == LevelSetFillIn::DENSE || cacheMode == LevelSetCacheMode::FULL) {
            cacheId = m_cellCacheCollection->insert<LevelSetCartesianKernel::CellDenseCacheContainer<value_t>>();
        } else if (fillIn == LevelSetFillIn::SPARSE) {
            cacheId = m_cellCacheCollection->insert<LevelSetCartesianKernel::CellSparseCacheContainer<value_t>>();
        }
    } else if (dynamic_cast<const LevelSetOctreeKernel *>(m_kernel)){
        if (fillIn == LevelSetFillIn::DENSE || cacheMode == LevelSetCacheMode::FULL) {
            cacheId = m_cellCacheCollection->insert<LevelSetOctreeKernel::CellDenseCacheContainer<value_t>>(piercedCacheKernel, piercedCacheSyncMode);
        } else if (fillIn == LevelSetFillIn::SPARSE) {
            cacheId = m_cellCacheCollection->insert<LevelSetOctreeKernel::CellSparseCacheContainer<value_t>>();
        }
    } else if (dynamic_cast<const LevelSetUnstructuredKernel *>(m_kernel)){
        if (fillIn == LevelSetFillIn::DENSE || cacheMode == LevelSetCacheMode::FULL) {
            cacheId = m_cellCacheCollection->insert<LevelSetUnstructuredKernel::CellDenseCacheContainer<value_t>>(piercedCacheKernel, piercedCacheSyncMode);
        } else if (fillIn == LevelSetFillIn::SPARSE) {
            cacheId = m_cellCacheCollection->insert<LevelSetUnstructuredKernel::CellSparseCacheContainer<value_t>>();
        }
    }

    return cacheId;
}

/*!
 * Get a pointer to the cell cache for the specified field.
 *
 * If no cache was registered for the specified field, a null pointer is returned.
 *
 * \param fieldset is the fieldset for which the caches will be added
 * \result A pointer to the cell cache for the specified field.
 */
template<typename value_t>
typename LevelSetObject::CellValueCache<value_t> * LevelSetObject::getFieldCellCache(LevelSetField field) const
{
    std::size_t fieldIndex = static_cast<std::size_t>(field);
    std::size_t cacheId = m_cellFieldCacheIds[fieldIndex];
    if (cacheId == CellCacheCollection::NULL_ID) {
        return nullptr;
    }

    return getCellCache<value_t>(cacheId);
}

/*!
 * Get a pointer to the specified cell cache.
 *
 * If specified cell cache was registered, a null pointer is returned.
 *
 * \param cacheId the id of the cell that will be unregistered
 * \result A pointer to the specified cell cache.
 */
template<typename value_t>
typename LevelSetObject::CellValueCache<value_t> * LevelSetObject::getCellCache(std::size_t cacheId) const
{
    return m_cellCacheCollection->at<value_t>(cacheId);
}

}

#endif
