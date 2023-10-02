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
 * Get a pointer to the cell cache for the specified field.
 *
 * If a cache was not enabled for the specified field, a null pointer is returned. If a cache
 * was enabled for the field, but it was not yet created, it will be created now.
 *
 * \param field is the field for which the caches is requested
 * \result A pointer to the cell cache for the specified field.
 */
template<typename value_t>
typename LevelSetObject::CellCacheCollection::ValueCache<value_t> * LevelSetObject::getFieldCellCache(LevelSetField field) const
{
    LevelSetCacheMode cacheMode = getFieldCellCacheMode(field);
    if (cacheMode == LevelSetCacheMode::NONE) {
        return nullptr;
    }

    std::size_t cacheId = getFieldCellCacheId(field);

    return getCellCache<value_t>(cacheId);
}

/*!
 * Get a pointer to the specified cell cache.
 *
 * If the specified cell cache was not registered, a null pointer is returned.
 *
 * \param cacheId the id of the cache that will be unregistered
 * \result A pointer to the specified cell cache or a null pointer if the cache was not registered.
 */
template<typename value_t>
typename LevelSetObject::CellCacheCollection::ValueCache<value_t> * LevelSetObject::getCellCache(std::size_t cacheId) const
{
    if (cacheId == CellCacheCollection::NULL_CACHE_ID) {
        return nullptr;
    }

    return (*m_cellCacheCollection)[cacheId].getCache<value_t>();
}

/*!
 * Create the cache that will be used for storing cell information of the specified field.
 *
 * \param field is the field for which the caches will be added
 * \result The id associated with the registered cache.
 */
template<typename value_t>
std::size_t LevelSetObject::createFieldCellCache(LevelSetField field)
{
    // Create cache
    LevelSetCacheMode cacheMode = getFieldCellCacheMode(field);

    std::size_t cacheId;
    if (cacheMode != LevelSetCacheMode::NONE) {
        LevelSetFillIn expectedFillIn;
        if (m_kernel->getExpectedFillIn() == LevelSetFillIn::DENSE || cacheMode == LevelSetCacheMode::NARROW_BAND) {
            expectedFillIn = LevelSetFillIn::DENSE;
        } else {
            expectedFillIn = LevelSetFillIn::SPARSE;
        }

        cacheId = createCellCache<value_t>(expectedFillIn);
    } else {
        cacheId = CellCacheCollection::NULL_CACHE_ID;
    }

    // Update field properties
    std::size_t fieldIndex = static_cast<std::size_t>(field);
    m_cellFieldCacheIds[fieldIndex] = cacheId;

    return cacheId;
}

/*!
 * Create the cache that will be used for storing cell information of the specified field.
 *
 * \param expectedFillIn is the expected fill-in of the cache
 * \result The id associated with the registered cache.
 */
template<typename value_t>
std::size_t LevelSetObject::createCellCache(LevelSetFillIn expectedFillIn)
{
    // Create the cache
    std::size_t cacheId = CellCacheCollection::NULL_CACHE_ID;
    if (dynamic_cast<const LevelSetCartesianKernel *>(m_kernel)){
        if (expectedFillIn == LevelSetFillIn::DENSE) {
            cacheId = m_cellCacheCollection->insert<LevelSetCartesianKernel::CellDenseCacheContainer<value_t>>();
        } else if (expectedFillIn == LevelSetFillIn::SPARSE) {
            cacheId = m_cellCacheCollection->insert<LevelSetCartesianKernel::CellSparseCacheContainer<value_t>>();
        }
    } else if (dynamic_cast<const LevelSetOctreeKernel *>(m_kernel)){
        if (expectedFillIn == LevelSetFillIn::DENSE) {
            cacheId = m_cellCacheCollection->insert<LevelSetOctreeKernel::CellDenseCacheContainer<value_t>>();
        } else if (expectedFillIn == LevelSetFillIn::SPARSE) {
            cacheId = m_cellCacheCollection->insert<LevelSetOctreeKernel::CellSparseCacheContainer<value_t>>();
        }
    } else if (dynamic_cast<const LevelSetUnstructuredKernel *>(m_kernel)){
        if (expectedFillIn == LevelSetFillIn::DENSE) {
            cacheId = m_cellCacheCollection->insert<LevelSetUnstructuredKernel::CellDenseCacheContainer<value_t>>();
        } else if (expectedFillIn == LevelSetFillIn::SPARSE) {
            cacheId = m_cellCacheCollection->insert<LevelSetUnstructuredKernel::CellSparseCacheContainer<value_t>>();
        }
    }

    return cacheId;
}

/*!
 * Write the VTK data associated with the specified field to the given stream.
 *
 * Only data currently stored in the cache will be written, no new field data will be evaluated
 * by the function.
 *
 * @param[in] stream is the output stream
 * @param[in] format is the format which must be used. Supported options
 * are "ascii" or "appended". For "appended" type an unformatted binary
 * stream must be used
 * @param[in] field is the field
 * @param[in] evaluator is the functor that should be used to evaluate the field
 * @param[in] fallback is the functor that should be used to evaluate the field fallback value
 */
template<typename value_t, typename evaluator_t, typename fallback_t>
void LevelSetObject::flushVTKOutputData(std::fstream &stream, VTKFormat format, LevelSetField field,
                                        const evaluator_t evaluator, const fallback_t fallback) const
{
    LevelSetCacheMode fieldCacheMode = getFieldCellCacheMode(field);
    switch (fieldCacheMode) {

    case LevelSetCacheMode::ON_DEMAND:
    {
        CellCacheCollection::Cache *cache = getFieldCellCache(field);
        for (const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange()) {
            long cellId = cell.getId();
            if (cache->contains(cellId)) {
                flushValue(stream, format, evaluator(cellId));
            } else {
                flushValue(stream, format, fallback(cellId));
            }
        }

        break;
    }

    case LevelSetCacheMode::NARROW_BAND:
    {
        for (const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange()) {
            long cellId = cell.getId();
            if (isCellInNarrowBand(cellId)) {
                flushValue(stream, format, evaluator(cellId));
            } else {
                flushValue(stream, format, fallback(cellId));
            }
        }

        break;
    }

    case LevelSetCacheMode::FULL:
    {
        for (const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange()) {
            long cellId = cell.getId();
            flushValue(stream, format, evaluator(cellId));
        }

        break;
    }

    default:
    {
        for (const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange()) {
            long cellId = cell.getId();
            flushValue(stream, format, fallback(cellId));
        }

        break;
    }

    }
}

/*!
 * Evaluate the specified field for the given cell.
 *
 * The value is first searched in the cache. If the cache doesn't contain an entry for the
 * specified cell, the value is evaluated from scratch and the cache is updated.
 *
 * @param[in] field is the field that should be evaluated
 * @param[in] id is the id of the cell where the field should be evaluated
 * @param[in] evaluator is the function that should be used to evaluate the field when an exact
 * value is requested
 * @param[in] fallback is the function that should be used to evaluate the field when a dummy
 * value is requested
 */
template<typename value_t, typename evaluator_t, typename fallback_t>
value_t LevelSetObject::evalCellFieldCached(LevelSetField field, long id, const evaluator_t &evaluator, const fallback_t &fallback) const
{
    // Try fetching the field from the cache
    CellCacheCollection::ValueCache<value_t> *cache = getFieldCellCache<value_t>(field);
    if (cache) {
        typename CellCacheCollection::ValueCache<value_t>::Entry cacheEntry = cache->findEntry(id);
        if (cacheEntry.isValid()) {
            return *cacheEntry;
        }
    }

    // Evaluate the field
    value_t value = evalCellField<value_t>(field, id, evaluator, fallback);

    // Update the cache
    fillFieldCellCache(field, id, value);

    // Return the value
    return value;
}

/*!
 * Evaluate the specified field for the given cell.
 * @param[in] field is the field that should be evaluated
 * @param[in] id is the id of the cell where the field should be evaluated
 * @param[in] evaluator is the function that should be used to evaluate the field when an exact
 * value is requested
 * @param[in] fallback is the function that should be used to evaluate the field when a dummy
 * value is requested
 */
template<typename value_t, typename evaluator_t, typename fallback_t>
value_t LevelSetObject::evalCellField(LevelSetField field, long id, const evaluator_t &evaluator, const fallback_t &fallback) const
{
    BITPIT_UNUSED(field);

    // Get cell location
    LevelSetCellLocation cellLocation = getCellLocation(id);

    // Early return if the zone has not been detected
    if (cellLocation == LevelSetCellLocation::UNKNOWN) {
        return evaluator(id);
    }

    // Early return if the fallback value should be returned
    if (cellLocation == LevelSetCellLocation::BULK) {
        if (m_cellBulkEvaluationMode != LevelSetBulkEvaluationMode::EXACT) {
            return fallback(id);
        }
    }

    // Return the value
    return evaluator(id);
}

/*!
 * Fill the cache value associated with the given cell ids for the specified field.
 *
 * Depending on the cache mode and on the bulk evaluation mode the cached may not need to be
 * filled for the specified cell. If the cache already contains a value for the specified cell,
 * that value will be replaced with the given one.
 *
 * @param[in] field is the field that should be evaluated
 * @param[in] id is the id of the cell where the field should be evaluated
 * @param[in] value is the value that will be added to the cell
 */
template<typename value_t>
void LevelSetObject::fillFieldCellCache(LevelSetField field, long id, const value_t &value) const
{
    // Early return if no cache is associated with the field
    LevelSetCacheMode levelsetCacheMode = getFieldCellCacheMode(field);
    if (levelsetCacheMode == LevelSetCacheMode::NONE) {
        return;
    }

    // Early return if the cache doesn't need to be updated
    //
    // There are some cases were the cache doesn't need to be updated:
    //  - cells outside the narrow band should not be added to caches in "narrow band" mode.
    //
    // If the zone of the cell has not been detected, it's not possible to check the aforementioned
    // conditions. In this case the cache can be updated only if it is operating in "full" mode and
    // the bulk is evaluated using an exact method.
    LevelSetCellLocation cellLocation = getCellLocation(id);
    if (cellLocation == LevelSetCellLocation::UNKNOWN) {
        if (levelsetCacheMode != LevelSetCacheMode::FULL) {
            return;
        }

        if (getCellBulkEvaluationMode() != LevelSetBulkEvaluationMode::EXACT) {
            return;
        }
    } else if (cellLocation == LevelSetCellLocation::BULK) {
        if (levelsetCacheMode == LevelSetCacheMode::NARROW_BAND) {
            return;
        }
    }

    // Update the cache
    CellCacheCollection::ValueCache<value_t> *cache = getFieldCellCache<value_t>(field);
    cache->insertEntry(id, value);
}

}

#endif
