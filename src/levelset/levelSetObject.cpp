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

#include "levelSetObject.hpp"
#include "levelSetCommon.hpp"

#if BITPIT_ENABLE_MPI
#    include <mpi.h>
#endif

#include <algorithm>
#include <numeric>
#include <vector>

namespace bitpit {

/*!
    @interface LevelSetObject
    @ingroup levelset
    @brief Interface class for all objects with respect to whom the levelset function may be computed.

    Evaluation of the fields inside the narrow band is always performed using an exact algorithm,
    on the other hand evaluation of the fields in the bulk can be performed choosing one of the
    following modes:
        - NONE, no data is evaluated (a dummy value will be returned when requesting any
        - data);
        - SIGN_PROPAGATION, sign is propagated from the narrow band, no other data will
          be evaluated (a dummy value will be returned when requesting data other than
          the sign, whenever possible the evaluated sign will be used for assigning the
          correct sign to the dummy value);
        - EXACT, exact data is evaluated.

    Fields evaluated by the object on cell centroids can be cached. Cache can be enabled for each
    individual field separately and one of the following modes can be chosen:
        - NONE, no caching will be performed;
        - ON_DEMAND, data are cached only where explicitly evaluated;
        - NARROW_BAND, Data are cached only inside the narrow band;
        - FULL, data are cached in the whole domain.
*/

/*!
 * Controls if the cache stores signed values.
 *
 * The implemented algorithms requires that the cache stores unsigned values.
 */
const bool LevelSetObject::CELL_CACHE_IS_SIGNED = false;

/*!
 * Intersection mode that should be used when detecting the location of a cell.
 *
 * This intersection mode needs to be FAST_GUARANTEE_FALSE because we can accept cells tagged
 * as intersecting the surface but actually not intersecting it, but we cannot accept cells
 * tagged as not intersecting the surface but actually intersecting it. There are places where
 * we make use of the fact that the identification of the location uses the FAST_GUARANTEE_FALSE
 * mode (for example when checking if a cell intersects the zero-levelset iso-surface).
 */
const LevelSetIntersectionMode LevelSetObject::CELL_LOCATION_INTERSECTION_MODE = LevelSetIntersectionMode::FAST_GUARANTEE_FALSE;

/*!
 * Constructor
 * @param[in] id id assigned to object
 */
LevelSetObject::LevelSetObject(int id)
    : m_kernel(nullptr),
      m_defaultSignedLevelSet(false),
      m_narrowBandSize(levelSetDefaults::NARROW_BAND_SIZE),
      m_cellLocationCacheId(CellCacheCollection::NULL_CACHE_ID),
      m_cellPropagatedSignCacheId(CellCacheCollection::NULL_CACHE_ID),
      m_nReferences(0),
      m_cellBulkEvaluationMode(LevelSetBulkEvaluationMode::SIGN_PROPAGATION),
      m_cellFieldCacheModes(static_cast<std::size_t>(LevelSetField::COUNT), LevelSetCacheMode::NONE),
      m_cellFieldCacheIds(static_cast<std::size_t>(LevelSetField::COUNT), CellCacheCollection::NULL_CACHE_ID)
{
    setId(id);
}

/*!
 * Copy constructor
 * @param[in] other is another object whose content is copied in this object
 */
LevelSetObject::LevelSetObject(const LevelSetObject &other)
    : m_kernel(other.m_kernel),
      m_defaultSignedLevelSet(other.m_defaultSignedLevelSet),
      m_enabledOutputFields(other.m_enabledOutputFields),
      m_narrowBandSize(other.m_narrowBandSize),
      m_cellLocationCacheId(other.m_cellLocationCacheId),
      m_cellPropagatedSignCacheId(other.m_cellPropagatedSignCacheId),
      m_id(other.m_id),
      m_nReferences(other.m_nReferences),
      m_cellBulkEvaluationMode(other.m_cellBulkEvaluationMode),
      m_cellCacheCollection(new CellCacheCollection(*(other.m_cellCacheCollection))),
      m_cellFieldCacheModes(other.m_cellFieldCacheModes),
      m_cellFieldCacheIds(other.m_cellFieldCacheIds)

{
    for ( const auto &fieldEntry : m_enabledOutputFields ) {
        enableVTKOutput(fieldEntry.first, true);
    }
}

/*!
 * Move constructor
 * @param[in] other is another object whose content is moved in this object
 */
LevelSetObject::LevelSetObject(LevelSetObject &&other)
    : m_kernel(std::move(other.m_kernel)),
      m_defaultSignedLevelSet(std::move(other.m_defaultSignedLevelSet)),
      m_narrowBandSize(std::move(other.m_narrowBandSize)),
      m_cellLocationCacheId(std::move(other.m_cellLocationCacheId)),
      m_cellPropagatedSignCacheId(std::move(other.m_cellPropagatedSignCacheId)),
      m_id(std::move(other.m_id)),
      m_nReferences(std::move(other.m_nReferences)),
      m_cellBulkEvaluationMode(std::move(other.m_cellBulkEvaluationMode)),
      m_cellCacheCollection(std::move(other.m_cellCacheCollection)),
      m_cellFieldCacheModes(std::move(other.m_cellFieldCacheModes)),
      m_cellFieldCacheIds(std::move(other.m_cellFieldCacheIds))
{
    for ( const auto &fieldEntry : other.m_enabledOutputFields ) {
        enableVTKOutput(fieldEntry.first, true);
    }

    for ( const auto &fieldEntry : other.m_enabledOutputFields ) {
        other.enableVTKOutput(fieldEntry.first, false);
    }
}

/*!
 * Destructor.
 */
LevelSetObject::~LevelSetObject() {
    // Disable all output for the object
    if (m_kernel) {
        try {
            // Disable output
            LevelSetFieldset enabledOutputFieldset;
            enabledOutputFieldset.reserve(m_enabledOutputFields.size());
            for ( const auto &fieldEntry : m_enabledOutputFields ) {
                enabledOutputFieldset.push_back(fieldEntry.first);
            }

            enableVTKOutput(enabledOutputFieldset, false);
        } catch (const std::exception &exception) {
            BITPIT_UNUSED(exception);

            // Nothing to do
        }
    }
}

/*!
 * Get the list of supported field.
 * @result The list of supported field.
 */
LevelSetFieldset LevelSetObject::getSupportedFields() const {

    LevelSetFieldset supportedFields;
    supportedFields.push_back(LevelSetField::SIGN);
    supportedFields.push_back(LevelSetField::VALUE);
    supportedFields.push_back(LevelSetField::GRADIENT);

    return supportedFields;

}

/*!
 * Sets the identifier of object.
 * @param[in] id is the identifier of the object
 */
void LevelSetObject::setId(int id) {
    m_id = id;
}

/*!
 * Increment reference count.
 */
std::size_t LevelSetObject::incrementReferenceCount() {

    ++m_nReferences;

    return m_nReferences;

}

/*!
 * Decrement reference count.
 */
std::size_t LevelSetObject::decrementReferenceCount() {

    assert(m_nReferences > 0);
    --m_nReferences;

    return m_nReferences;

}

/*!
 * Count how many times the object is referenced by other objects.
 * @result The number of times the object is referenced by other objects.
 */
std::size_t LevelSetObject::getReferenceCount() const {

    return m_nReferences ;

}

/*!
 * Set whether a signed or unsigned levelset is used when signdness is not explicitly specified.
 * This function is only needed for guarantee backwards compatibility with older versions. In
 * such versions, functions for evaluating levelset information were not taking in input the
 * signdness.
 * @param signedLevelSet controls if signed levelset function will be used
 */
void LevelSetObject::setDefaultLevelSetSigndness(bool signedLevelSet) {
    m_defaultSignedLevelSet = signedLevelSet;
}

/*!
 * Sets the kernel for the object.
 * @param[in] kernel is the kernel that will be associated with the object
 */
void LevelSetObject::setKernel(LevelSetKernel *kernel) {
    // Set kernel
    m_kernel = kernel;

    // Enable output
    for ( const auto &fieldEntry : m_enabledOutputFields ) {
        enableVTKOutput( fieldEntry.first, true ) ;
    }

    // Create cell cache collection
    PiercedVector<Cell, long> *cacheKernel = &(m_kernel->getMesh()->getCells());
    m_cellCacheCollection = std::unique_ptr<CellCacheCollection>(new CellCacheCollection(cacheKernel));

    // Evaluate data
    evaluate();
}

/*!
 * Gets a pointer to the kernel of the object.
 * @return A pointer to the kernel of the object.
 */
LevelSetKernel * LevelSetObject::getKernel() {
    return m_kernel;
}

/*!
 * Gets a constant pointer to the kernel of the object.
 * @return A constant pointer to the kernel of the object.
 */
const LevelSetKernel * LevelSetObject::getKernel() const {
    return m_kernel;
}

/*!
 * Get the id of the object.
 * @return The id of the object.
 */
int LevelSetObject::getId( ) const {
    return m_id ;
}

/*!
 * Check if the levelset is a primary object (e.g. of a surface triangulation) or not (e.g. derived
 * by boolean operations between two levelsets)
 * @return Returns true if object is a primary object, false otherwise.
 */
bool LevelSetObject::isPrimary( ) const {
    return true;
}

/*!
 * Evaluate object data.
 */
void LevelSetObject::evaluate()
{
    // Get fields with empty caches
    std::vector<LevelSetField> fillableFields;
    for (std::size_t fieldIndex = 0; fieldIndex < static_cast<std::size_t>(LevelSetField::COUNT); ++fieldIndex) {
        LevelSetField field = static_cast<LevelSetField>(fieldIndex);

        LevelSetCacheMode fieldCacheMode = getFieldCellCacheMode(field);
        if (fieldCacheMode == LevelSetCacheMode::NONE) {
            continue;
        }

        std::size_t fieldCacheId = getFieldCellCacheId(field);
        if (fieldCacheId != CellCacheCollection::NULL_CACHE_ID) {
            continue;
        }

        fillableFields.push_back(field);
    }

    // Create field cell caches
    for (LevelSetField field : fillableFields) {
        createFieldCellCache(field);
    }

    // Evaluate narrow band data
    evaluateCellNarrowBandData();

    // Fill field cell caches inside the narrow band
    fillFieldCellCaches(LevelSetZone::NARROW_BAND, fillableFields);

    // Evaluate bulk data
    evaluateCellBulkData();

    // Fill field cell caches inside the bulk
    fillFieldCellCaches(LevelSetZone::BULK, fillableFields);
}

/*!
 * Updates object data after a mesh update.
 * @param[in] adaptionData are the information about the mesh update
 */
void LevelSetObject::update( const std::vector<adaption::Info> &adaptionData )
{
    // Get fields with non-empty caches
    std::vector<LevelSetField> fillableFields;
    for (std::size_t fieldIndex = 0; fieldIndex < static_cast<std::size_t>(LevelSetField::COUNT); ++fieldIndex) {
        LevelSetField field = static_cast<LevelSetField>(fieldIndex);

        LevelSetCacheMode fieldCacheMode = getFieldCellCacheMode(field);
        if (fieldCacheMode == LevelSetCacheMode::NONE) {
            continue;
        }

        std::size_t fieldCacheId = getFieldCellCacheId(field);
        if (fieldCacheId == CellCacheCollection::NULL_CACHE_ID) {
            continue;
        }

        fillableFields.push_back(field);
    }

    // Update cell caches
    adaptCellCaches(adaptionData);

    // Evaluate narrow band data
    updateCellNarrowBandData(adaptionData);

    // Fill field cell caches inside the narrow band
    fillFieldCellCaches(LevelSetZone::NARROW_BAND, fillableFields, adaptionData);

    // Evaluate bulk data
    updateCellBulkData(adaptionData);

    // Fill field cell caches inside the bulk
    fillFieldCellCaches(LevelSetZone::BULK, fillableFields, adaptionData);
}

/*!
 * Get the physical size of the narrow band.
 *
 * The size of the narrow band is the absolute distance from the zero-levelset iso surface below
 * which a point is considered belonging to the narrow band. Setting the size of the narrow band
 * to LEVELSET_NARROW_BAND_UNLIMITED means that the whole domain belongs to the narrow band.
 * Regardless of the specified size, the narrow band will always contain the intersected cells
 * and their neighbours.
 *
 * \result The physical size of the narrow band.
 */
double LevelSetObject::getNarrowBandSize() const
{
    return m_narrowBandSize;
}

/*!
 * Set the physical size of the narrow band.
 *
 * The size of the narrow band is the absolute distance from the zero-levelset iso surface below
 * which a point is considered belonging to the narrow band. Setting the size of the narrow band
 * to LEVELSET_NARROW_BAND_UNLIMITED means that the whole domain belongs to the narrow band.
 * Regardless of the specified size, the narrow band will always contain the intersected cells
 * and their neighbours.
 *
 * \param size is the physical size of the narrow band
 */
void LevelSetObject::setNarrowBandSize(double size)
{
    // Early return if size doesn't need to be updated
    if (m_narrowBandSize == size) {
        return;
    }

    // Destroy current data
    //
    // Current information need to be destroyed (not just cleared) because the updated narrow
    // band size may require different storage types.
    if (getKernel()) {
        destroyCellNarrowBandData();
        destroyCellBulkData();

        for (std::size_t fieldIndex = 0; fieldIndex < static_cast<std::size_t>(LevelSetField::COUNT); ++fieldIndex) {
            LevelSetField field = static_cast<LevelSetField>(fieldIndex);
            LevelSetCacheMode fieldCacheMode = getFieldCellCacheMode(field);

            bool destroyFieldCache = true;
            if (fieldCacheMode == LevelSetCacheMode::NONE) {
                destroyFieldCache = false;
            } else if (fieldCacheMode == LevelSetCacheMode::FULL && getCellBulkEvaluationMode() == LevelSetBulkEvaluationMode::EXACT) {
                destroyFieldCache = false;
            }

            if (destroyFieldCache) {
                destroyFieldCellCache(field);
            }
        }
    }

    // Set the narrow band size
    m_narrowBandSize = size;

    // Evaluate updated information
    if (getKernel()) {
        evaluate();
    }
}

/*!
 * Evaluate cell data inside the narrow band.
 */
void LevelSetObject::evaluateCellNarrowBandData()
{
    // Identify cells locations
    if (m_narrowBandSize != LEVELSET_NARROW_BAND_UNLIMITED) {
        createCellLocationCache();
        fillCellLocationCache();
    }
}

/*!
 * Update cell data inside the narrow band after a mesh update.
 *
 * \param adaptionData are the information about the mesh update
 */
void LevelSetObject::updateCellNarrowBandData( const std::vector<adaption::Info> &adaptionData )
{
    // Update cells locations
    if (m_narrowBandSize != LEVELSET_NARROW_BAND_UNLIMITED) {
        fillCellLocationCache(adaptionData);
    }
}

/*!
 * Destroy cell data inside the narrow band.
 */
void LevelSetObject::destroyCellNarrowBandData( )
{
    // Identify cells inside the narrow band
    if (m_narrowBandSize != LEVELSET_NARROW_BAND_UNLIMITED) {
        destroyCellLocationCache();
    }
}

/**
 * Get the location associated with the specified cell.
 *
 * @param[in] id is the cell id
 * @return The location associated with the specified cell.
 */
LevelSetCellLocation LevelSetObject::getCellLocation(long id) const
{
    // Early return if the object is empty
    if (empty()) {
        return LevelSetCellLocation::BULK;
    }

    // Early return if the narrow band is unlimited
    if (m_narrowBandSize == LEVELSET_NARROW_BAND_UNLIMITED) {
        return LevelSetCellLocation::NARROW_BAND_UNDEFINED;
    }

    // Early return if the narrow band cells has not been identified yet
    if (m_cellLocationCacheId == CellCacheCollection::NULL_CACHE_ID) {
        return LevelSetCellLocation::UNKNOWN;
    }

    // Get the location from the cache
    CellCacheCollection::ValueCache<char> *locationCache = getCellCache<char>(m_cellLocationCacheId);
    CellCacheCollection::ValueCache<char>::Entry locationCacheEntry = locationCache->findEntry(id);

    if (locationCacheEntry.isValid()) {
        return static_cast<LevelSetCellLocation>(*locationCacheEntry);
    } else {
        return LevelSetCellLocation::UNKNOWN;
    }
}

/*!
 * Get the zone of the specified cell.
 *
 * @param[in] id is the cell id
 * \result Return the zone of the specified cell.
 */
LevelSetZone LevelSetObject::getCellZone(long id) const
{
    LevelSetCellLocation cellLocation = getCellLocation(id);
    switch (cellLocation) {

    case LevelSetCellLocation::BULK:
        return LevelSetZone::BULK;

    case LevelSetCellLocation::NARROW_BAND_DISTANCE:
    case LevelSetCellLocation::NARROW_BAND_INTERSECTED:
    case LevelSetCellLocation::NARROW_BAND_NEIGHBOUR:
    case LevelSetCellLocation::NARROW_BAND_UNDEFINED:
        return LevelSetZone::NARROW_BAND;

    case LevelSetCellLocation::UNKNOWN:
        throw std::runtime_error("Cell location identification has not been performed!");

    default:
        throw std::runtime_error("Unable to identify cell zone!");

    }
}

/*!
 * Fill the cache that contains the location associated to the cells.
 *
 * A cell can be either in the narrow band or in the bulk. It will be considered inside the narrow
 * band if one of the following conditions holds:
 *  - its distance from the surface is less than the narrow band size;
 *  - it intersects the zero-levelset iso-surface (intersections are checked using the
 *    FAST_GUARANTEE_FALSE criterium);
 *  - one of its neighbors intersects the zero-levelset iso-surface.
 */
void LevelSetObject::fillCellLocationCache()
{
    // Mesh information
    const VolumeKernel &mesh = *(m_kernel->getMesh()) ;

    // Get the ids associated with internal cells
    std::vector<long> internalCellIds;
    if (const VolCartesian *cartesianMesh = dynamic_cast<const VolCartesian *>(&mesh)) {
        long nCells = cartesianMesh->getCellCount();
        internalCellIds.resize(nCells);
        std::iota(internalCellIds.begin(), internalCellIds.end(), 0);
    } else {
        long nCells = mesh.getCellCount();
        VolumeKernel::CellConstIterator internalCellsBegin = mesh.internalCellConstBegin();
        VolumeKernel::CellConstIterator internalCellsEnd = mesh.internalCellConstEnd();

        internalCellIds.reserve(nCells);
        for (VolumeKernel::CellConstIterator cellItr = internalCellsBegin; cellItr != internalCellsEnd; ++cellItr) {
            long cellId = cellItr.getId();
            internalCellIds.push_back(cellId);
        }
    }

    // Get cell location cache
    CellCacheCollection::ValueCache<char> *locationCache = getCellCache<char>(m_cellLocationCacheId);

    // Clear cell location cache
    clearCellCache(m_cellLocationCacheId, false);

    // Assign an unknown location to the internal cells.
    for (long cellId : internalCellIds) {
        locationCache->insertEntry(cellId, static_cast<char>(LevelSetCellLocation::UNKNOWN));
    }

    // Identify internal cells that are geometrically inside the narrow band
    //
    // A cell is geometrically inside the narrow band if its distance from the surface is smaller
    // than the narrow band side or if it intersects the surface.
    //
    // Cells that intersect the zero-levelset iso-surface need to be identified because they will
    // be further processed for finding which of their neighbour is inside the narrow band.
    std::vector<long> intersectedCellIds;
    if (!empty()) {
        for (long cellId : internalCellIds) {
            // Fill location cache for cells geometrically inside the narrow band
            LevelSetCellLocation cellLocation = fillCellGeometricNarrowBandLocationCache(cellId);

            // Track intersected cells
            if (cellLocation == LevelSetCellLocation::NARROW_BAND_INTERSECTED) {
                intersectedCellIds.push_back(cellId);
            }
        }
    }

    // Identify face neighbours of cells that intersect the surface
    for (long cellId : intersectedCellIds) {
        auto neighProcessor = [this, &locationCache](long neighId, int layer) {
            BITPIT_UNUSED(layer);

            // Skip neighbours whose region have already been identified
            if (getCellLocation(neighId) != LevelSetCellLocation::UNKNOWN) {
                return false;
            }

            // The neighbour is inside the narrow band
            locationCache->insertEntry(neighId, static_cast<char>(LevelSetCellLocation::NARROW_BAND_NEIGHBOUR));

            // Continue processing the other neighbours
            return false;
        };

        mesh.processCellFaceNeighbours(cellId, 1, neighProcessor);
    }

    // Cells whose location is still unknown are in the bulk
    for (long cellId : internalCellIds) {
        CellCacheCollection::ValueCache<char>::Entry locationCacheEntry = locationCache->findEntry(cellId);
        assert(locationCacheEntry.isValid());
        if (static_cast<LevelSetCellLocation>(*locationCacheEntry) == LevelSetCellLocation::UNKNOWN) {
            locationCache->insertEntry(cellId, static_cast<char>(LevelSetCellLocation::BULK));
        }
    }

#if BITPIT_ENABLE_MPI==1
    // Exchange ghost data
    if (mesh.isPartitioned()) {
        std::unique_ptr<DataCommunicator> dataCommunicator = m_kernel->createDataCommunicator();
        startCellCacheExchange(mesh.getGhostCellExchangeSources(), m_cellLocationCacheId, dataCommunicator.get());
        completeCellCacheExchange(mesh.getGhostCellExchangeTargets(), m_cellLocationCacheId, dataCommunicator.get());
    }
#endif
}

/*!
 * Fill the cache that contains the location associated to the cells.
 *
 * A cell can be either in the narrow band or in the bulk. It will be considered inside the narrow
 * band if one of the following conditions holds:
 *  - its distance from the surface is less than the narrow band size;
 *  - it intersects the zero-levelset iso-surface (intersections are checked using the
 *    FAST_GUARANTEE_FALSE criterium);
 *  - one of its neighbors intersects the zero-levelset iso-surface.
 *
 * \param adaptionData are the information about the mesh update
 */
void LevelSetObject::fillCellLocationCache(const std::vector<adaption::Info> &adaptionData)
{
    // Mesh information
    const VolumeKernel &mesh = *(m_kernel->getMesh()) ;

    // Get cell location cache
    CellCacheCollection::ValueCache<char> *locationCache = getCellCache<char>(m_cellLocationCacheId);

    // Initialize cells updated by the mesh adaption
    std::unordered_set<long> unknownZoneCellIds;
    for (const adaption::Info &adaptionInfo : adaptionData) {
        if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
            continue;
        }

        // Skip received cells when the cache is non-volatile
        //
        // If the cache is non-volatile, data on exchanged cells has been communicated during
        // cache adaption.
        if (adaptionInfo.type == adaption::Type::TYPE_PARTITION_RECV && !locationCache->isVolatile()) {
            continue;
        }

        // Add current cells to the process list
        for (long cellId : adaptionInfo.current) {

            // Assign an unknown location to the newly created cells
            locationCache->insertEntry(cellId, static_cast<char>(LevelSetCellLocation::UNKNOWN));

            // Add internal cells to the process list
            bool isCellProcessable = true;
#if BITPIT_ENABLE_MPI==1
            if (mesh.isPartitioned()) {
                VolumeKernel::CellConstIterator cellItr = mesh.getCellConstIterator(cellId);
                isCellProcessable = cellItr->isInterior();
            }
#endif

            if (isCellProcessable) {
                unknownZoneCellIds.insert(cellId);
            }
        }
    }

    // Identify internal cells that are geometrically inside the narrow band
    //
    // A cell is geometrically inside the narrow band if its distance from the surface is smaller
    // than the narrow band side or if it intersects the surface.
    //
    // Cells that intersect the surface are identified, because they will be further processed for
    // finding which of their neighbours is inside the narrow band.
    std::unordered_set<long> intersectedCellIds;
    if (!empty()) {
        auto cellIdItr = unknownZoneCellIds.begin();
        while (cellIdItr != unknownZoneCellIds.end()) {
            long cellId = *cellIdItr;

            // Fill location cache for cells geometrically inside the narrow band
            LevelSetCellLocation cellLocation = fillCellGeometricNarrowBandLocationCache(cellId);

            // Track intersected cells
            if (cellLocation == LevelSetCellLocation::NARROW_BAND_INTERSECTED) {
                intersectedCellIds.insert(cellId);
            }

            // Remove from the list cell whose region has been identified
            if (cellLocation != LevelSetCellLocation::UNKNOWN) {
                cellIdItr = unknownZoneCellIds.erase(cellIdItr);
            } else {
                ++cellIdItr;
            }
        }
    }

#if BITPIT_ENABLE_MPI==1
    // Exchange ghost data
    if (mesh.isPartitioned()) {
        std::unique_ptr<DataCommunicator> dataCommunicator = m_kernel->createDataCommunicator();
        startCellCacheExchange(mesh.getGhostCellExchangeSources(), m_cellLocationCacheId, dataCommunicator.get());
        completeCellCacheExchange(mesh.getGhostCellExchangeTargets(), m_cellLocationCacheId, dataCommunicator.get());
    }
#endif

    // Track intersected neighbours of cells whose location was not yet identified
    //
    // We need to identified cells that are the narrow band because of a neighbour that was
    // not affected by the mesh update.
    if (!empty()) {
        for (long cellId : unknownZoneCellIds) {
            // Process cells whose location has not been identified
            //
            // For cells whose location has not been identified, the neighbours that intersect
            // the surface need to be tracked.
            auto neighProcessor = [this, &intersectedCellIds](long neighId, int layer) {
                BITPIT_UNUSED(layer);

                // Skip neighbours that are aready identified as intersected
                if (intersectedCellIds.count(neighId) > 0) {
                    return false;
                }

                // Skip neighbours that doesn't intersect the surface
                LevelSetCellLocation neighLocation = getCellLocation(neighId);
                if (neighLocation != LevelSetCellLocation::NARROW_BAND_INTERSECTED) {
                    return false;
                }

                // Track neighbours that intersect the surface
                intersectedCellIds.insert(neighId);

                // Continue processing the other neighbours
                return false;
            };

            mesh.processCellFaceNeighbours(cellId, 1, neighProcessor);
        }
    }

    // Identify neighbours of cells that intersect the surface
    //
    // We may need the distance of cells that are not processed,
    for (long cellId : intersectedCellIds) {
        // Process face neighbours
        auto neighProcessor = [this, &locationCache, &unknownZoneCellIds](long neighId, int layer) {
            BITPIT_UNUSED(layer);

            // Skip neighbours whose region has already been identified
            if (getCellLocation(neighId) != LevelSetCellLocation::UNKNOWN) {
                return false;
            }

            // The neighbour is inside the narrow band
            locationCache->insertEntry(neighId, static_cast<char>(LevelSetCellLocation::NARROW_BAND_NEIGHBOUR));
            unknownZoneCellIds.erase(neighId);

            // Continue processing the other neighbours
            return false;
        };

        mesh.processCellFaceNeighbours(cellId, 1, neighProcessor);
    }

    // Cells whose location is still unknown are in the bulk
    for (long cellId : unknownZoneCellIds) {
        locationCache->insertEntry(cellId, static_cast<char>(LevelSetCellLocation::BULK));
    }

#if BITPIT_ENABLE_MPI==1
    // Exchange ghost data
    if (mesh.isPartitioned()) {
        std::unique_ptr<DataCommunicator> dataCommunicator = m_kernel->createDataCommunicator();
        startCellCacheExchange(mesh.getGhostCellExchangeSources(), m_cellLocationCacheId, dataCommunicator.get());
        completeCellCacheExchange(mesh.getGhostCellExchangeTargets(), m_cellLocationCacheId, dataCommunicator.get());
    }
#endif
}

/*!
 * Fill location cache for the specified cell if it is geometrically inside the narrow band
 *
 * A cell is geometrically inside the narrow band if its distance from the surface is smaller
 * than the narrow band side or if it intersects the surface.
 *
 * This function may require the evaluation of some levelset fields. To improve performance,
 * it is important to attempt filling the cache of the evaluated fields. It is then up to the
 * caches to decide if the fields can be cached.
 *
 * \param[in] id is the cell id
 * \return The cell location of the cache or LevelSetCellLocation::UNKNOWN if the cell
 * location was not identified.
 */
LevelSetCellLocation LevelSetObject::fillCellGeometricNarrowBandLocationCache(long id)
{
    // Get cell information
    double cellCacheValue    = evalCellValue(id, CELL_CACHE_IS_SIGNED);
    double cellUnsignedValue = std::abs(cellCacheValue);

    // Identify cells that are geometrically inside the narrow band
    //
    // First we need to check if the cell intersectes the surface, and only if it
    // deosn't we should check if its distance is lower than the narrow band size.
    LevelSetCellLocation cellLocation = LevelSetCellLocation::UNKNOWN;
    if (_intersectSurface(id, cellUnsignedValue, CELL_LOCATION_INTERSECTION_MODE) == LevelSetIntersectionStatus::TRUE) {
        cellLocation = LevelSetCellLocation::NARROW_BAND_INTERSECTED;
    } else if (cellUnsignedValue <= m_narrowBandSize) {
        cellLocation = LevelSetCellLocation::NARROW_BAND_DISTANCE;
    }

    if (cellLocation != LevelSetCellLocation::UNKNOWN) {
        CellCacheCollection::ValueCache<char> *locationCache = getCellCache<char>(m_cellLocationCacheId);
        locationCache->insertEntry(id, static_cast<char>(cellLocation));
    }

    // Fill the cache of the evaluated fields
    //
    // Now that the cell location has been identified, we can fill the cache of the
    // evaluated fields.
    fillFieldCellCache(LevelSetField::VALUE, id, cellCacheValue);

    // Return the location
    return cellLocation;
}

/**
 * Create the cache that will be used for storing cell location information.
 *
 * \param cacheId is the id that will be associated with the cache, if a NULL_ID is specified
 * the cache id will be assigned automatically
 * \result The id associated with the cache.
 */
std::size_t LevelSetObject::createCellLocationCache(std::size_t cacheId)
{
    m_cellLocationCacheId = createCellCache<char>(LevelSetFillIn::DENSE, cacheId);

    return m_cellLocationCacheId;
}

/**
 * Destroy identification of cell location.
 */
void LevelSetObject::destroyCellLocationCache()
{
    destroyCellCache(m_cellLocationCacheId);
    m_cellLocationCacheId = CellCacheCollection::NULL_CACHE_ID;
}

/*!
 * Check if the specified cell lies inside the narrow band.
 *
 * A cell is considered inside the narrow band if one of the following conditions hold:
 *  - its distance from the surface is less than the narrow band size;
 *  - it intersects the zero-levelset iso-surface (intersections are checked using the
 *    FAST_GUARANTEE_FALSE criterium);
 *  - one of its neighbors intersects the zero-levelset iso-surface.
 *
 * If no caches with "narrow band" mode have been filled, the function may return wrong
 * results if the cell is on the last layer of ghosts.
 *
 * \param[in] id is the cell id
 * \result Return true if the cell is in the narrow band, false otherwise.
 */
bool LevelSetObject::isCellInNarrowBand(long id)const
{
    LevelSetZone zone = getCellZone(id);

    return (zone == LevelSetZone::NARROW_BAND);
}

/*!
 * Check if the specified point lies inside the narrow band.
 *
 * The value of the levelset is evaluated and compared with the specified narrow band size.
 *
 * If a point is inside a cell that belongs to the narrow band because it is a neighbour of a
 * cell with a different levelset sign, its function may identify the point as outside the
 * narrow band. That's because it will only compare the levelset of the point with the narrow
 * band size. The extreme case is when the narrow band size is set to zero and being a neighbour
 * of a cell with a different levelset sign is the only criterion to identify cells inside the
 * narrow band. In this situation, this function will identify as inside the narrow band only
 * the points that lie on the levelset-zero iso surface.
 *
 * \param point are the coordinates of the point
 * \result Return true if the cell is in the narrow band, false otherwise.
 */
bool LevelSetObject::isInNarrowBand(const std::array<double,3> &point)const
{
    // Early return if the object is empty
    if (empty()) {
        return false;
    }

    // Early return if no narrow band was set
    if (m_narrowBandSize == LEVELSET_NARROW_BAND_UNLIMITED) {
        return true;
    }

    // Explicitly check if the cell is inside the narrow band
    double value = evalValue(point, false);
    if (value <= m_narrowBandSize) {
        return true;
    }

    return false;
}

/**
 * Get the mode that will be used to evaluate cell data in the bulk.
 *
 * \result The mode that will be used to evaluate cell data in the bulk.
 */
LevelSetBulkEvaluationMode LevelSetObject::getCellBulkEvaluationMode() const
{
    return m_cellBulkEvaluationMode;
}

/**
 * Set the mode that will be used to evaluate cell data in the bulk.
 *
 * \param evaluationMode is the mode that will be used to evaluate cell data in the bulk.
 */
void LevelSetObject::setCellBulkEvaluationMode(LevelSetBulkEvaluationMode evaluationMode)
{
    // Early return if the bulk evaluation mode doesn't need to be updated
    if (getCellBulkEvaluationMode() == evaluationMode) {
        return;
    }

    // Destroy bulk data
    //
    // Current data need to be destroyed (not just cleared) because the updated bulk evaluation
    // mode may require different storage types.
    if (getKernel()) {
        destroyCellBulkData();

        for (std::size_t fieldIndex = 0; fieldIndex < static_cast<std::size_t>(LevelSetField::COUNT); ++fieldIndex) {
            LevelSetField field = static_cast<LevelSetField>(fieldIndex);
            LevelSetCacheMode fieldCacheMode = getFieldCellCacheMode(field);

            bool destroyFieldCache = false;
            if (fieldCacheMode == LevelSetCacheMode::ON_DEMAND) {
                destroyFieldCache = true;
            } else if (fieldCacheMode == LevelSetCacheMode::FULL) {
                destroyFieldCache = true;
            }

            if (destroyFieldCache) {
                destroyFieldCellCache(field);
            }
        }
    }

    // Set bulk evaluation mode
    m_cellBulkEvaluationMode = evaluationMode;

    // Evaluate data
    if (getKernel()) {
        evaluate();
    }
}

/*!
 * Evaluate cell data inside the bulk.
 */
void LevelSetObject::evaluateCellBulkData()
{
    if (m_cellBulkEvaluationMode == LevelSetBulkEvaluationMode::SIGN_PROPAGATION) {
        createCellPropagatedSignCache();
        fillCellPropagatedSignCache();
    }
}

/*!
 * Update cell data inside the bulk.
 */
void LevelSetObject::updateCellBulkData( const std::vector<adaption::Info> &adaptionData )
{
    BITPIT_UNUSED(adaptionData);

    if (m_cellBulkEvaluationMode == LevelSetBulkEvaluationMode::SIGN_PROPAGATION) {
        fillCellPropagatedSignCache();
    }
}

/**
 * Destroy cell data in the bulk.
 */
void LevelSetObject::destroyCellBulkData()
{
    if (m_cellBulkEvaluationMode == LevelSetBulkEvaluationMode::SIGN_PROPAGATION) {
        destroyCellPropagatedSignCache();
    }
}

/*!
 * Fill the cache that contains the propagated cell sign.
 */
void LevelSetObject::fillCellPropagatedSignCache()
{
    const VolumeKernel &mesh = *(m_kernel->getMesh());

#if BITPIT_ENABLE_MPI==1
    // Check if the communications are needed
    bool ghostExchangeNeeded = mesh.isPartitioned();

    // Initialize ghost communications
    std::unique_ptr<DataCommunicator> ghostCommunicator;
    if (ghostExchangeNeeded) {
        ghostCommunicator = std::unique_ptr<DataCommunicator>(new DataCommunicator(mesh.getCommunicator()));

        for (auto &entry : mesh.getGhostCellExchangeSources()) {
            const int rank = entry.first;
            const std::vector<long> &sources = entry.second;
            ghostCommunicator->setSend(rank, sources.size() * (sizeof(unsigned char) + sizeof(char)));
        }

        for (auto &entry : mesh.getGhostCellExchangeTargets()) {
            const int rank = entry.first;
            const std::vector<long> &targets = entry.second;
            ghostCommunicator->setRecv(rank, targets.size() * (sizeof(unsigned char) + sizeof(char)));
        }
    }
#endif

    // Get the ids associated with internal cells
    std::vector<long> internalCellIds;
    if (const VolCartesian *cartesianMesh = dynamic_cast<const VolCartesian *>(&mesh)) {
        long nCells = cartesianMesh->getCellCount();
        internalCellIds.resize(nCells);
        std::iota(internalCellIds.begin(), internalCellIds.end(), 0);
    } else {
        long nCells = mesh.getCellCount();
        VolumeKernel::CellConstIterator internalCellsBegin = mesh.internalCellConstBegin();
        VolumeKernel::CellConstIterator internalCellsEnd = mesh.internalCellConstEnd();

        internalCellIds.reserve(nCells);
        for (VolumeKernel::CellConstIterator cellItr = internalCellsBegin; cellItr != internalCellsEnd; ++cellItr) {
            long cellId = cellItr.getId();
            internalCellIds.push_back(cellId);
        }
    }

    // Get cache for sign propagation
    CellCacheCollection::ValueCache<char> *propagatedSignCache = getCellCache<char>(m_cellPropagatedSignCacheId);

    // Clear sign propagation cache
    clearCellCache(m_cellPropagatedSignCacheId, false);

    // Initialize sign propagation
    //
    // Sign propagation will start from the internal cells within the narrow band.
    // It is not possible to popagate a null sign.
    std::unordered_set<long> seedIds;
    for (long cellId : internalCellIds) {
        if (isCellInNarrowBand(cellId)) {
            char cellSign = static_cast<char>(evalCellSign(cellId));
            propagatedSignCache->insertEntry(cellId, cellSign);
            if (cellSign != 0) {
                seedIds.insert(cellId);
            }
        }
    }

    // Propagate the sign
    std::vector<long> processList;
    while (true) {
        bool signMismatch = false;
        while (!seedIds.empty()) {
            // Get seed
            long seedId = *(seedIds.begin());
            seedIds.erase(seedIds.begin());

            // Get seed sign
            CellCacheCollection::ValueCache<char>::Entry seedSignEntry = propagatedSignCache->findEntry(seedId);
            char seedSign = *(seedSignEntry);

            // Propagate seed sign
            //
            // We need to continue processing until the process list is empty, even
            // if we have already processed the internal cells. The process list may
            // contain ghost cells that carry information received by other processes.
            processList.assign(1, seedId);
            while (!processList.empty()) {
                // Get cell to process
                long cellId = processList.back();
                processList.pop_back();

                // Set the sign of the cell
                //
                // If the sign of the cell has already been set, we check if it matches
                // the seed sign and then we skip the cell, otherwise we set the sign of
                // the cell and we process its neighoburs. Once the sign of a cell is
                // set, we can remove it form the seed list (there is no need to start
                // again the porpagation from that cell, beause it will not reach any
                // new portions of the domain).
                if (cellId != seedId) {
                    CellCacheCollection::ValueCache<char>::Entry cellSignEntry = propagatedSignCache->findEntry(cellId);
                    if (cellSignEntry.isValid()) {
                        char cellSign = *(cellSignEntry);
                        if (cellSign != seedSign) {
                            signMismatch = true;
                            log::error() << "Sign mismatch on cell " << cellId << "!" << std::endl;
                            break;
                        } else {
                            continue;
                        }
                    } else {
                        propagatedSignCache->insertEntry(cellId, seedSign);
                        seedIds.erase(cellId);
                    }
                }

                // Add face neighbours to the process list
                auto neighProcessor = [&mesh, &propagatedSignCache, &processList](long neighId, int layer) {
                    BITPIT_UNUSED(layer);
#if BITPIT_ENABLE_MPI==0
                    BITPIT_UNUSED(mesh);
#endif

                    // Skip neighbours whose sign has already been set
                    if (propagatedSignCache->contains(neighId)) {
                        return false;
                    }

#if BITPIT_ENABLE_MPI==1
                    // Skip ghost neighbours
                    if (mesh.isPartitioned()) {
                        VolumeKernel::CellConstIterator neighItr = mesh.getCellConstIterator(neighId);
                        const Cell &neigh = *neighItr;
                        if (!neigh.isInterior()) {
                            return false;
                        }
                    }
#endif

                    // Add neighbour to the process list
                    processList.push_back(neighId);

                    // Continue processing the other neighbours
                    return false;
                };

                mesh.processCellFaceNeighbours(cellId, 1, neighProcessor);
            }

            // Stop processing if a sign mismatch was detected
            if (signMismatch) {
                break;
            }
        }

        // Raise an error if a sign mismatch was detected
#if BITPIT_ENABLE_MPI==1
        if (ghostExchangeNeeded) {
            MPI_Allreduce(MPI_IN_PLACE, &signMismatch, 1, MPI_C_BOOL, MPI_LOR, mesh.getCommunicator());
        }
#endif

        if (signMismatch) {
            throw std::runtime_error("Sign mismatch in sign propagation!");
        }

#if BITPIT_ENABLE_MPI==1
        // Add ghost to the process list
        //
        // A ghost should be added to the process list if it has been processed by the owner
        // but is hasn't yet processed by the current process.
        if (ghostExchangeNeeded) {
            // Exchange ghost data
            ghostCommunicator->startAllRecvs();

            for(auto &entry : mesh.getGhostCellExchangeSources()) {
                int rank = entry.first;
                SendBuffer &buffer = ghostCommunicator->getSendBuffer(rank);
                for (long cellId : entry.second) {
                    CellCacheCollection::ValueCache<char>::Entry cellSignEntry = propagatedSignCache->findEntry(cellId);

                    char cellHasSign = cellSignEntry.isValid() ? 1 : 0;
                    buffer << cellHasSign;

                    if (cellHasSign) {
                        char cellSign = *(cellSignEntry);
                        buffer << cellSign;
                    }
                }
                ghostCommunicator->startSend(rank);
            }

            bool ghostSignMismatch = false;
            int nPendingRecvs = ghostCommunicator->getRecvCount();
            while (nPendingRecvs != 0) {
                int rank = ghostCommunicator->waitAnyRecv();
                RecvBuffer buffer = ghostCommunicator->getRecvBuffer(rank);

                for (long cellId : mesh.getGhostCellExchangeTargets(rank)) {
                    char sourceHasSign;
                    buffer >> sourceHasSign;

                    if (sourceHasSign == 1) {
                        char sourceSign;
                        buffer >> sourceSign;

                        if (!propagatedSignCache->contains(cellId)) {
                            seedIds.insert(cellId);
                            propagatedSignCache->insertEntry(cellId, sourceSign);
                        } else {
                            char cellCachedSign = *(propagatedSignCache->findEntry(cellId));
                            if (cellCachedSign != sourceSign) {
                                ghostSignMismatch = true;
                                log::error() << "Sign mismatch on ghost cell " << cellId << "!" << getId() << std::endl;
                                break;
                            }
                        }
                    }
                }

                --nPendingRecvs;

                // Stop processing if a sign mismatch was detected
                if (ghostSignMismatch) {
                    break;
                }
            }

            ghostCommunicator->waitAllSends();

            // Raise an error if there are ghosts with a sign mismatch
            if (ghostExchangeNeeded) {
                MPI_Allreduce(MPI_IN_PLACE, &ghostSignMismatch, 1, MPI_C_BOOL, MPI_LOR, mesh.getCommunicator());
            }

            if (ghostSignMismatch) {
                throw std::runtime_error("Sign mismatch in sign propagation!");
            }
        }
#endif

        // Detect if the propagation is completed
        bool completed = seedIds.empty();
#if BITPIT_ENABLE_MPI==1
        if (ghostExchangeNeeded) {
            MPI_Allreduce(MPI_IN_PLACE, &completed, 1, MPI_C_BOOL, MPI_LAND, mesh.getCommunicator());
        }
#endif

        if (completed) {
            break;
        }
    }
}

/**
 * Create the cache that will be used for storing cell propagated sign.
 *
 * \param cacheId is the id that will be associated with the cache, if a NULL_ID is specified
 * the cache id will be assigned automatically
 * \result The id associated with the cache.
 */
std::size_t LevelSetObject::createCellPropagatedSignCache(std::size_t cacheId)
{
    m_cellPropagatedSignCacheId = createCellCache<char>(LevelSetFillIn::DENSE, cacheId);

    return m_cellPropagatedSignCacheId;
}

/**
 * Destroy data related to sign propagation in the bulk.
 */
void LevelSetObject::destroyCellPropagatedSignCache()
{
    destroyCellCache(m_cellPropagatedSignCacheId);
    m_cellPropagatedSignCacheId = CellCacheCollection::NULL_CACHE_ID;
}

/*!
 * Function for checking if the specified cell intersects the zero-levelset iso-surface.
 *
 * If mode==LevelSetIntersectionMode::FAST_FUZZY the method will compare the levelset
 * value to tangent and bounding radius of a cell. If the value is smaller than the
 * tangent radius LevelSetIntersectionStatus::TRUE is returned, if it is larger than the
 * bounding radius LevelSetIntersectionStatus::FALSE is returned. If it is in-between
 * LevelSetIntersectionStatus::CLOSE is returned.
 *
 * If mode==LevelSetIntersectionMode::FAST_GUARANTEE_TRUE and the levelset value is
 * smaller than the rangent radius LevelSetIntersectionStatus::TRUE is returned,
 * otherwise LevelSetIntersectionStatus::FALSE.
 *
 * If mode==LevelSetIntersectionMode::FAST_GURANTEE_FALSE and the levelset value is
 * larger than the bounding radius LevelSetIntersectionStatus::FALSE is returned,
 * otherwise LevelSetIntersectionStatus::TRUE.
 *
 * If mode==LevelSetIntersectionMode::ACCURATE, the same checks of fuzzy mode are
 * performed, however, in the cases where fuzzy mode would return CLOSE, an additional
 * check on the intersection between the tangent plane at the projection point and the
 * cell is performed. Errors of the method are related to the ratio of surface curvature
 * over cell size.
 *
 * The bounding sphere is the sphere with the minimum radius that contains all the
 * cell vertices and has the center in the cell centroid.
 *
 * The tangent sphere is a sphere having the center in the level centroid and tangent
 * to the cell.
 *
 * @param[in] id cell id
 * @param[in] mode describes the types of check that should be performed
 * @return indicator regarding intersection
 */
LevelSetIntersectionStatus LevelSetObject::intersectSurface(long id, LevelSetIntersectionMode mode) const
{
    // Try evaluating intersection information from the cell location
    //
    // Regardless from the requested mode, a cell can intersect the zero-levelset iso-surface
    // only if it is inside the narrow band.
    //
    // If the requested mode is the same mode used for identify the cell location, we can get the
    // intersection information directly from the location.
    LevelSetCellLocation cellLocation = getCellLocation(id);
    if (cellLocation == LevelSetCellLocation::BULK) {
        return LevelSetIntersectionStatus::FALSE;
    } else if (mode == CELL_LOCATION_INTERSECTION_MODE) {
        if (cellLocation == LevelSetCellLocation::NARROW_BAND_INTERSECTED) {
            return LevelSetIntersectionStatus::TRUE;
        } else if (cellLocation != LevelSetCellLocation::NARROW_BAND_UNDEFINED) {
            return LevelSetIntersectionStatus::FALSE;
        }
    }

    // Check for intersection with zero-levelset iso-surface
    double distance = evalCellValue(id, false);

    return _intersectSurface(id, distance, mode);
}

/*!
 * Check if the specified cell intersects the zero-levelset iso-surface.
 *
 * If mode==LevelSetIntersectionMode::FAST_FUZZY the method will compare the levelset
 * value to tangent and bounding radius of a cell. If the value is smaller than the
 * tangent radius LevelSetIntersectionStatus::TRUE is returned, if it is larger than the
 * bounding radius LevelSetIntersectionStatus::FALSE is returned. If it is in-between
 * LevelSetIntersectionStatus::CLOSE is returned.
 *
 * If mode==LevelSetIntersectionMode::FAST_GUARANTEE_TRUE and the levelset value is
 * smaller than the rangent radius LevelSetIntersectionStatus::TRUE is returned,
 * otherwise LevelSetIntersectionStatus::FALSE.
 *
 * If mode==LevelSetIntersectionMode::FAST_GURANTEE_FALSE and the levelset value is
 * larger than the bounding radius LevelSetIntersectionStatus::FALSE is returned,
 * otherwise LevelSetIntersectionStatus::TRUE.
 *
 * If mode==LevelSetIntersectionMode::ACCURATE, the same checks of fuzzy mode are
 * performed, however, in the cases where fuzzy mode would return CLOSE, an additional
 * check on the intersection between the tangent plane at the projection point and the
 * cell is performed. Errors of the method are related to the ratio of surface curvature
 * over cell size.
 *
 * The bounding sphere is the sphere with the minimum radius that contains all the
 * cell vertices and has the center in the cell centroid.
 *
 * The tangent sphere is a sphere having the center in the level centroid and tangent
 * to the cell.
 *
 * @param[in] id cell id
 * @param[in] distance is the unsigned distance of the cell centroid from the zero-levelset
 * iso-surface
 * @param[in] mode describes the types of check that should be performed
 * @return indicator regarding intersection
 */
LevelSetIntersectionStatus LevelSetObject::_intersectSurface(long id, double distance, LevelSetIntersectionMode mode) const
{
    double distanceTolerance = m_kernel->getDistanceTolerance();

    switch(mode){
        case LevelSetIntersectionMode::FAST_GUARANTEE_TRUE:
        {
            double tangentSphere = m_kernel->computeCellTangentRadius(id) ;
            if(utils::DoubleFloatingLessEqual()(distance, tangentSphere, distanceTolerance, distanceTolerance)){
                return LevelSetIntersectionStatus::TRUE;
            } else {
                return LevelSetIntersectionStatus::FALSE;
            }

            break;
        }

        case LevelSetIntersectionMode::FAST_GUARANTEE_FALSE:
        {
            double boundingSphere = m_kernel->computeCellBoundingRadius(id) ;
            if(utils::DoubleFloatingGreater()(distance, boundingSphere, distanceTolerance, distanceTolerance)){
                return LevelSetIntersectionStatus::FALSE;
            } else {
                return LevelSetIntersectionStatus::TRUE;
            }

            break;
        }

        case LevelSetIntersectionMode::FAST_FUZZY:
        {
            double boundingSphere = m_kernel->computeCellBoundingRadius(id) ;
            if(utils::DoubleFloatingGreater()(distance, boundingSphere, distanceTolerance, distanceTolerance)){
                return LevelSetIntersectionStatus::FALSE;
            }

            double tangentSphere = m_kernel->computeCellTangentRadius(id) ;
            if(utils::DoubleFloatingLessEqual()(distance, tangentSphere, distanceTolerance, distanceTolerance)){
                return LevelSetIntersectionStatus::TRUE;
            }

            return LevelSetIntersectionStatus::CLOSE;

            break;
        }

        case LevelSetIntersectionMode::ACCURATE:
        {
            double boundingSphere = m_kernel->computeCellBoundingRadius(id) ;
            if(utils::DoubleFloatingGreater()(distance, boundingSphere, distanceTolerance, distanceTolerance)){
                return LevelSetIntersectionStatus::FALSE;
            }

            double tangentSphere = m_kernel->computeCellTangentRadius(id) ;
            if(utils::DoubleFloatingLessEqual()(distance, tangentSphere, distanceTolerance, distanceTolerance)){
                return LevelSetIntersectionStatus::TRUE;
            }

            std::array<double,3> root = evalCellProjectionPoint(id);
            std::array<double,3> normal = evalCellGradient(id, true);
            if( m_kernel->intersectCellPlane(id,root,normal, distanceTolerance) ){
                return LevelSetIntersectionStatus::TRUE;
            } else {
                return LevelSetIntersectionStatus::FALSE;
            }

            break;
        }
    }

    BITPIT_UNREACHABLE("cannot reach");

}

/*!
 * Evaluate levelset sign at the specified cell.
 * @param[in] id cell id
 * @result The sign of the levelset at the specified cell.
 */
short LevelSetObject::evalCellSign(long id) const {

    // Define sign evaluators
    auto evaluator = [this] (long id)
        {
            // Try fetching the sign from sign propagation
            CellCacheCollection::ValueCache<char> *propagatedSignCache = getCellCache<char>(m_cellPropagatedSignCacheId);
            if (propagatedSignCache) {
                CellCacheCollection::ValueCache<char>::Entry propagatedSignCacheEntry = propagatedSignCache->findEntry(id);
                if (propagatedSignCacheEntry.isValid()) {
                    return static_cast<short>(*propagatedSignCacheEntry);
                }
            }

            // Evaluate the sign
            return _evalCellSign(id);
        };

    auto fallback = [this] (long id)
        {
            // Try fetching the sign from sign propagation
            CellCacheCollection::ValueCache<char> *propagatedSignCache = getCellCache<char>(m_cellPropagatedSignCacheId);
            if (propagatedSignCache) {

                CellCacheCollection::ValueCache<char>::Entry propagatedSignCacheEntry = propagatedSignCache->findEntry(id);
                if (propagatedSignCacheEntry.isValid()) {
                    return static_cast<short>(*propagatedSignCacheEntry);
                }
            }

            // Return a dummy value
            return levelSetDefaults::SIGN;
        };

    LevelSetField field = LevelSetField::SIGN;
    short sign = evalCellFieldCached<short>(field, id, evaluator, fallback);

    return sign;
}

/*!
 * Evaluate levelset value at the specified cell.
 * @param[in] id cell id
 * @param[in] signedLevelSet controls if signed levelset function will be used
 * @result The value of the levelset at the specified cell.
 */
double LevelSetObject::evalCellValue(long id, bool signedLevelSet) const {

    // Evaluate signed value
    //
    // The value stored in the cache is unsigned.
    auto evaluator = [this] (long id)
        {
            return _evalCellValue(id, false);
        };

    auto fallback = [] (long id)
        {
            BITPIT_UNUSED(id);

            return levelSetDefaults::VALUE;
        };

    LevelSetField field = LevelSetField::VALUE;
    double value = evalCellFieldCached<double>(field, id, evaluator, fallback);

    // Evaluate the value with the correct signdness
    if (signedLevelSet) {
        value *= evalCellSign(id);
    }

    return value;
}

/*!
 * Evaluate levelset gradient at the specified cell.
 * @param[in] id cell id
 * @param[in] signedLevelSet controls if signed levelset function will be used
 * @result The gradient of the levelset at the specified cell.
 */
std::array<double,3> LevelSetObject::evalCellGradient(long id, bool signedLevelSet) const {
    // Evaluate signed gradient
    //
    // The gradient stored in the cache is unsigned.
    auto evaluator = [this] (long id)
        {
            return _evalCellGradient(id, false);
        };

    auto fallback = [] (long id)
        {
            BITPIT_UNUSED(id);

            return levelSetDefaults::GRADIENT;
        };

    LevelSetField field = LevelSetField::GRADIENT;
    std::array<double, 3> gradient = evalCellFieldCached<std::array<double, 3>>(field, id, evaluator, fallback);

    // Evaluate the gradient with the correct signdness
    if (signedLevelSet) {
        short cellSign = evalCellSign(id);
        if (cellSign <= 0) {
            gradient *= static_cast<double>(cellSign);
        }
    }

    return gradient;
}

/*!
 * Computes the projection point of the cell center, i.e. the closest
 * point to the cell center on the zero level set
 * @param[in] id cell id
 * @return the projection point
 */
std::array<double,3> LevelSetObject::evalCellProjectionPoint(long id) const {
    return m_kernel->computeCellCentroid(id) - evalCellValue(id, true) * evalCellGradient(id, true);
}

/*!
 * Evaluate levelset sign at the specified point.
 * @param point are the coordinates of the point
 * @result The sign of the levelset at the specified point.
 */
short LevelSetObject::evalSign(const std::array<double,3> &point) const {
    return _evalSign(point);
}

/*!
 * Evaluate levelset value at the specified point.
 * @param point are the coordinates of the point
 * @param[in] signedLevelSet controls if signed levelset function will be used
 * @result The value of the levelset at the specified point.
 */
double LevelSetObject::evalValue(const std::array<double,3> &point, bool signedLevelSet) const {
    return _evalValue(point, signedLevelSet);
}

/*!
 * Evaluate levelset gradient at the specified point.
 * @param point are the coordinates of the point
 * @param[in] signedLevelSet controls if signed levelset function will be used
 * @result The gradient of the levelset at the specified point.
 */
std::array<double,3> LevelSetObject::evalGradient(const std::array<double,3> &point, bool signedLevelSet) const {
    return _evalGradient(point, signedLevelSet);
}

/*!
 * Evaluate levelset sign at the specified point.
 *
 * \param point are the coordinates of the point
 * \result The sign of the levelset at the specified point.
 */
short LevelSetObject::_evalSign(const std::array<double,3> &point) const {
    return evalValueSign(this->evalValue(point, true));
}

/*!
 * Projects a vertex on the zero levelset
 * @param[in] point point coordinates
 * @return the projected point
 */
std::array<double,3> LevelSetObject::evalProjectionPoint(const std::array<double,3> &point) const{
    return point - evalValue(point, true) * evalGradient(point, true);
}

/*!
 * Evaluates the sign of the given levelset value.
 * @param[in] value is the levleset value
 * @return The sign of the given levelset value.
 */
short LevelSetObject::evalValueSign(double value) const{
    return static_cast<short>(sign(value));
}

/*!
 * Writes LevelSetObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::dump( std::ostream &stream ){
    // Identifier
    utils::binary::write(stream, m_id) ;
    utils::binary::write(stream, m_nReferences) ;

    // Object properties
    utils::binary::write(stream, m_defaultSignedLevelSet) ;
    utils::binary::write(stream, m_narrowBandSize) ;
    utils::binary::write(stream, m_cellBulkEvaluationMode) ;

    // Restore cell caches
    utils::binary::write(stream, m_cellLocationCacheId);
    utils::binary::write(stream, m_cellPropagatedSignCacheId);

    std::size_t nCaches = m_cellCacheCollection->size();
    utils::binary::write(stream, nCaches);
    for (std::size_t cacheId = 0; cacheId < nCaches; ++cacheId) {
        CellCacheCollection::Item &cacheItem = (*m_cellCacheCollection)[cacheId];

        // Skip items without a factory
        bool hasFactory = cacheItem.hasFactory();
        utils::binary::write(stream, hasFactory);
        if (!hasFactory) {
            continue;
        }

        // Field information
        LevelSetField field = LevelSetField::UNDEFINED;
        for (std::size_t fieldIndex = 0; fieldIndex < static_cast<std::size_t>(LevelSetField::COUNT); ++fieldIndex) {
            field = static_cast<LevelSetField>(fieldIndex);
            if (cacheId == getFieldCellCacheId(field)) {
                break;
            } else {
                field = LevelSetField::UNDEFINED;
            }
        }

        utils::binary::write(stream, field);
        if (field != LevelSetField::UNDEFINED) {
            utils::binary::write(stream, getFieldCellCacheMode(field));
        }

        // Dump the cache
        bool hasCache = cacheItem.hasCache();
        utils::binary::write(stream, hasCache);
        if (hasCache) {
            CellCacheCollection::Cache *cache = cacheItem.getCache();
            cache->dump(stream);
        }
    }

    // Output fields
    std::size_t nEnabledOutputFields = m_enabledOutputFields.size() ;
    utils::binary::write(stream, nEnabledOutputFields) ;
    for (const auto &fieldEntry : m_enabledOutputFields) {
        utils::binary::write(stream, fieldEntry.first) ;
        utils::binary::write(stream, fieldEntry.second) ;
    }
}

/*!
 * Reads LevelSetObject from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::restore( std::istream &stream ){
    // Disable output
    LevelSetFieldset enabledOutputFieldset;
    enabledOutputFieldset.reserve(m_enabledOutputFields.size());
    for ( const auto &fieldEntry : m_enabledOutputFields ) {
        enabledOutputFieldset.push_back(fieldEntry.first);
    }

    enableVTKOutput(enabledOutputFieldset, false);

    // Reset object information
    for (std::size_t fieldIndex = 0; fieldIndex < static_cast<std::size_t>(LevelSetField::COUNT); ++fieldIndex) {
        m_cellFieldCacheModes[fieldIndex] = LevelSetCacheMode::NONE;
        m_cellFieldCacheIds[fieldIndex]   = CellCacheCollection::NULL_CACHE_ID;
    }

    m_cellPropagatedSignCacheId = CellCacheCollection::NULL_CACHE_ID;
    m_cellLocationCacheId       = CellCacheCollection::NULL_CACHE_ID;

    m_cellCacheCollection->clear();

    // Identifier
    utils::binary::read(stream, m_id) ;
    utils::binary::read(stream, m_nReferences) ;

    // Object properties
    utils::binary::read(stream, m_defaultSignedLevelSet) ;
    utils::binary::read(stream, m_narrowBandSize) ;
    utils::binary::read(stream, m_cellBulkEvaluationMode) ;

    // Cell caches
    std::size_t expectedCellLocationCacheId;
    utils::binary::read(stream, expectedCellLocationCacheId);

    std::size_t expectedCellPropagatedSignCacheId;
    utils::binary::read(stream, expectedCellPropagatedSignCacheId);

    std::size_t nCaches;
    utils::binary::read(stream, nCaches);
    for (std::size_t cacheId = 0; cacheId < nCaches; ++cacheId) {
        // Skip items without a factory
        bool hasFactory;
        utils::binary::read(stream, hasFactory);
        if (!hasFactory) {
            continue;
        }

        // Field information
        LevelSetField field;
        utils::binary::read(stream, field);
        if (field != LevelSetField::UNDEFINED) {
            std::size_t fieldIndex = static_cast<std::size_t>(field);
            utils::binary::read(stream, m_cellFieldCacheModes[fieldIndex]);
        }

        // Restore the cache
        bool hasCache;
        utils::binary::read(stream, hasCache);
        if (hasCache) {
            // Create the cache
            if (cacheId == expectedCellLocationCacheId) {
                createCellLocationCache(cacheId);
            } else if (cacheId == expectedCellPropagatedSignCacheId) {
                createCellPropagatedSignCache(cacheId);
            } else if (field != LevelSetField::UNDEFINED) {
                createFieldCellCache(field, cacheId);
            } else {
                throw std::runtime_error("Unable to restore levelset object " + std::to_string(getId()) + "!");
            }

            // Restore cache contents
            CellCacheCollection::Cache *cache = getCellCache(cacheId);
            cache->restore(stream);
        } else {
        }
    }

    // Enable output
    std::size_t nEnabledVTKOutputs ;
    utils::binary::read(stream, nEnabledVTKOutputs) ;
    for (std::size_t i = 0; i < nEnabledVTKOutputs; ++i) {
        LevelSetField field ;
        std::string fieldName;
        utils::binary::read(stream, field) ;
        utils::binary::read(stream, fieldName) ;

        enableVTKOutput(field, fieldName);
    }
}

/*!
 * Enables or disables the VTK output
 * @param[in] fieldset is the fieldset that that should be enabled/disabled
 * @param[in] enable true for enabling, false for disabling
 */
void LevelSetObject::enableVTKOutput( const LevelSetFieldset &fieldset, bool enable) {

    std::stringstream objectNameStream;
    objectNameStream << getId();

    enableVTKOutput(fieldset, objectNameStream.str(), enable);

}

/*!
 * Enables or disables the VTK output
 * The output will be enabled only if the object supports it.
 * @param[in] fieldset is the fieldset that that should be enabled/disabled
 * @param[in] objectName is the name that will be associated with the object
 * @param[in] enable true for enabling, false for disabling
 */
void LevelSetObject::enableVTKOutput( const LevelSetFieldset &fieldset, const std::string &objectName, bool enable) {

    for (LevelSetField field : fieldset) {
        enableVTKOutput(field, objectName, enable);
    }

}

/*!
 * Enables or disables the VTK output
 * @param[in] field is the field that that should be enabled/disabled
 * @param[in] enable true for enabling, false for disabling
 */
void LevelSetObject::enableVTKOutput( LevelSetField field, bool enable) {

    std::stringstream objectNameStream;
    objectNameStream << getId();

    enableVTKOutput(field, objectNameStream.str(), enable);

}

/*!
 * Enables or disables the VTK output
 * The output will be enabled only if the object supports it.
 * @param[in] field is the field that that should be enabled/disabled
 * @param[in] objectName is the name that will be associated with the object
 * @param[in] enable true for enabling, false for disabling
 */
void LevelSetObject::enableVTKOutput( LevelSetField field, const std::string &objectName, bool enable) {

    // Discard fields that are not supported
    LevelSetFieldset supportedFields = getSupportedFields();
    if (std::find(supportedFields.begin(), supportedFields.end(), field) == supportedFields.end()) {
        return;
    }

    // Check if the state of the filed is already the requested one
    if (enable == hasVTKOutputData(field, objectName)) {
        return;
    }

    // Process the field
    if (!enable) {
        removeVTKOutputData(field, objectName) ;
        m_enabledOutputFields.erase(field) ;
    } else {
        addVTKOutputData(field, objectName) ;
        m_enabledOutputFields.insert({field, getVTKOutputDataName(field, objectName)}) ;
    }

}

/*!
 * Enables or disables the VTK output
 * @param[in] field is the field that that should be enabled/disabled
 * @param[in] enable true for enabling, false for disabling
 */
void LevelSetObject::enableVTKOutput( LevelSetWriteField field, bool enable) {

    std::stringstream objectNameStream;
    objectNameStream << getId();

    enableVTKOutput(field, objectNameStream.str(), enable);

}

/*!
 * Enables or disables the VTK output
 * The output will be enabled only if the object supports it.
 * @param[in] writeField is the write field that that should be enabled/disabled
 * @param[in] objectName is the name that will be associated with the object
 * @param[in] enable true for enabling, false for disabling
 */
void LevelSetObject::enableVTKOutput( LevelSetWriteField writeField, const std::string &objectName, bool enable) {

    LevelSetFieldset fieldset;
    if( writeField==LevelSetWriteField::ALL){
        fieldset = getSupportedFields();

    } else if ( writeField==LevelSetWriteField::DEFAULT){
        fieldset.push_back(LevelSetField::VALUE);
        fieldset.push_back(LevelSetField::GRADIENT);

    } else {
        LevelSetField field = static_cast<LevelSetField>(writeField);

        LevelSetFieldset supportedFields = getSupportedFields();
        if (std::find(supportedFields.begin(), supportedFields.end(), field) == supportedFields.end()) {
            log::warning() << "The specified field is not supported by the levelset object" << std::endl;
            return;
        }

        fieldset.push_back(field);
    }

    enableVTKOutput( fieldset, objectName, enable);

}

/*!
 * Check if the VTK writer has data associated with the specified field.
 *
 * @param[in] field is the field
 * @param[in] objectName is the name that will be associated with the object
 * @result True if the VTK writer has data associated with the specified field,
 * false otherwise.
 */
bool LevelSetObject::hasVTKOutputData( LevelSetField field, const std::string &objectName) const {

    VTK &vtkWriter = m_kernel->getMesh()->getVTK() ;
    std::string name = getVTKOutputDataName(field, objectName);

    return vtkWriter.hasData(name);

}

/*!
 * Remove the VTK data associated with the specified field.
 *
 * @param[in] field is the field
 * @param[in] objectName is the name that will be associated with the object
 */
void LevelSetObject::removeVTKOutputData( LevelSetField field, const std::string &objectName) {

    VTK &vtkWriter = m_kernel->getMesh()->getVTK() ;
    std::string name = getVTKOutputDataName(field, objectName);

    vtkWriter.removeData(name);

}

/*!
 * Add the VTK data associated with the specified field.
 *
 * @param[in] field is the field
 * @param[in] objectName is the name that will be associated with the object
 */
void LevelSetObject::addVTKOutputData( LevelSetField field, const std::string &objectName) {

    VTK &vtkWriter = m_kernel->getMesh()->getVTK() ;
    std::string name = getVTKOutputDataName(field, objectName);

    switch(field){

        case LevelSetField::VALUE:
            vtkWriter.addData<double>( name, VTKFieldType::SCALAR, VTKLocation::CELL, this);
            break;

        case LevelSetField::SIGN:
            vtkWriter.addData<short>( name, VTKFieldType::SCALAR, VTKLocation::CELL, this);
            break;

        case LevelSetField::GRADIENT:
            vtkWriter.addData<double>( name, VTKFieldType::VECTOR, VTKLocation::CELL, this);
            break;

        default:
            throw std::runtime_error("Unsupported value of field in LevelSetObject::addDataToVTK() ");
            break;

    }

}

/*!
 * Get the name that will be used by the VTK writer for the specifed data.
 *
 * @param[in] field is the field
 * @param[in] objectName is the name that will be associated with the object
 * @result The name that will be used by the VTK writer for the specifed data.
 */
std::string LevelSetObject::getVTKOutputDataName( LevelSetField field, const std::string &objectName) const {

    std::stringstream nameStream;
    nameStream << "levelset" << getVTKOutputFieldName(field) << "_" << objectName;
    std::string name = nameStream.str();

    return name;

}

/*!
 * Get the name that will be used by the VTK writer for the specifed field.
 *
 * @param[in] field is the field
 * @result The name that will be used by the VTK writer for the specifed field.
 */
std::string LevelSetObject::getVTKOutputFieldName( LevelSetField field) const {

    switch(field){

        case LevelSetField::VALUE:
            return "Value";

        case LevelSetField::SIGN:
            return "Sign";

        case LevelSetField::GRADIENT:
            return "Gradient";

        default:
            throw std::runtime_error("Unsupported value of field in LevelSetObject::addDataToVTK() ");
            break;

    }

}

/*!
 * Interface for writing data to the VTK stream.
 *
 * @param[in] stream output stream
 * @param[in] name is the name of the data to be written. Either user
 * data or patch data
 * @param[in] format is the format which must be used. Supported options
 * are "ascii" or "appended". For "appended" type an unformatted binary
 * stream must be used
 */
void LevelSetObject::flushData( std::fstream &stream, const std::string &name, VTKFormat format){

    for ( const auto &fieldEntry : m_enabledOutputFields ) {
        const std::string &fieldName = fieldEntry.second;
        if (utils::string::keywordInString(name, fieldName)) {
            LevelSetField field = fieldEntry.first;
            flushVTKOutputData(stream, format, field);
        }
    }

}

/*!
 * Write the VTK data associated with the specified field to the given stream.
 *
 * Only data currently stored in the cache will be written, no new field data will be evaluated
 * by the function.
 *
 * @param[in] field is the field
 * @param[in] stream is the output stream
 * @param[in] format is the format which must be used. Supported options
 * are "ascii" or "appended". For "appended" type an unformatted binary
 * stream must be used
 */
void LevelSetObject::flushVTKOutputData(std::fstream &stream, VTKFormat format, LevelSetField field) const {

    switch(field) {

    case LevelSetField::VALUE:
    {
        auto evaluator = [this] (long id) { return evalCellValue(id, true); };
        auto fallback = [] (long id) { BITPIT_UNUSED(id); return levelSetDefaults::VALUE; };
        flushVTKOutputData<double>(stream, format, field, evaluator, fallback);
        break;
    }

    case LevelSetField::SIGN:
    {
        auto evaluator = [this] (long id) { return (short) evalCellSign(id); };
        auto fallback = [] (long id) { BITPIT_UNUSED(id); return levelSetDefaults::SIGN; };
        flushVTKOutputData<short>(stream, format, field, evaluator, fallback);
        break;
    }

    case LevelSetField::GRADIENT:
    {
        auto evaluator = [this] (long id) { return evalCellGradient(id, true); };
        auto fallback = [] (long id) { BITPIT_UNUSED(id); return levelSetDefaults::GRADIENT; };
        flushVTKOutputData<std::array<double, 3>>(stream, format, field, evaluator, fallback);
        break;
    }

    default:
    {
        throw std::runtime_error("Unable to write the field.");
    }

    }
}

#if BITPIT_ENABLE_MPI
/*!
 * Start data exchange for the specified cell cache.
 * @param[in] sendCellIds is the list of cell ids to send
 * @param[in] cacheId is the id of the caches whose data will be exchanged
 * @param[in,out] dataCommunicator is the data communicator that will be used
 * for the data exchange
 */
void LevelSetObject::startCellCacheExchange( const std::unordered_map<int, std::vector<long>> &sendCellIds,
                                            std::size_t cacheId, DataCommunicator *dataCommunicator) const
{
    startCellCachesExchange(sendCellIds, std::vector<std::size_t>(1, cacheId), dataCommunicator);
}

/*!
 * Start data exchange for the specified cell caches.
 * @param[in] sendCellIds is the list of cell ids to send
 * @param[in] cacheIds are the ids of the caches whose data will be exchanged
 * @param[in,out] dataCommunicator is the data communicator that will be used
 * for the data exchange
 */
void LevelSetObject::startCellCachesExchange( const std::unordered_map<int, std::vector<long>> &sendCellIds,
                                              const std::vector<std::size_t> &cacheIds,
                                              DataCommunicator *dataCommunicator) const {

    // Fill the exchange buffer with the content of the cache
    dataCommunicator->clearAllSends(false);
    for (const auto &entry : sendCellIds) {
        // Create an empty send
        int rank = entry.first;
        dataCommunicator->setSend(rank, 0);

        // Write data in the buffer
        const std::vector<long> &rankSendCellIds = entry.second;
        SendBuffer &buffer = dataCommunicator->getSendBuffer(rank);
        for (std::size_t cacheId : cacheIds) {
            getCellCache(cacheId)->writeBuffer(rankSendCellIds, buffer);
        }
        buffer.squeeze( ) ;
    }

    // Discover the receives
    dataCommunicator->discoverRecvs();

    // Start the sends
    dataCommunicator->startAllRecvs();

    // Start the sends
    dataCommunicator->startAllSends();

}

/*!
 * Complete exchange of cell cache data.
 * @param[in] recvCellIds is the list of cell ids to receive
 * @param[in] cacheId is the id of the cache whose data will be exchanged
 * @param[in,out] dataCommunicator is the data communicator that will be used
 * for the data exchange
 */
void LevelSetObject::completeCellCacheExchange( const std::unordered_map<int, std::vector<long>> &recvCellIds,
                                                std::size_t cacheId, DataCommunicator *dataCommunicator)
{
    completeCellCachesExchange(recvCellIds, std::vector<std::size_t>(1, cacheId), dataCommunicator);
}

/*!
 * Complete exchange of cell cache data.
 * @param[in] recvCellIds is the list of cell ids to receive
 * @param[in] cacheId is the id of the cache whose data will be exchanged
 * @param[in,out] dataCommunicator is the data communicator that will be used
 * for the data exchange
 */
void LevelSetObject::completeCellCachesExchange( const std::unordered_map<int, std::vector<long>> &recvCellIds,
                                                 const std::vector<std::size_t> &cacheIds,
                                                 DataCommunicator *dataCommunicator){

    // Read cache data from the exchange buffer
    int nCompletedRecvs = 0;
    while (nCompletedRecvs < dataCommunicator->getRecvCount()) {
        int rank = dataCommunicator->waitAnyRecv();

        const std::vector<long> &rankRecvCellIds = recvCellIds.at(rank);
        RecvBuffer &buffer = dataCommunicator->getRecvBuffer(rank);
        for (std::size_t cacheId : cacheIds) {
            getCellCache(cacheId)->readBuffer(rankRecvCellIds, buffer);
        }

        ++nCompletedRecvs;
    }

    // Wait completion of all sends
    dataCommunicator->waitAllSends();

}
#endif

/*!
 * Get the cache mode associated with the specified field.
 *
 * \param field is the specified field
 * \result The cache mode associated with the specified field.
 */
LevelSetCacheMode LevelSetObject::getFieldCellCacheMode(LevelSetField field) const
{
    std::size_t fieldIndex = static_cast<std::size_t>(field);

    return m_cellFieldCacheModes[fieldIndex];
}

/*!
 * Enable the cache for the specified specified field.
 *
 * If a cache with the same mode is already defined for the specified field, the function will
 * exit without performing any action. If a cache with a different mode is already defined for
 * the specified field, the existing cache will be destroyed and a new cache with the requested
 * mode will be created from scratch.
 *
 * \param field is the field for which cache will be enabled
 * \param cacheMode is the cache mode that will be used for field
 */
void LevelSetObject::enableFieldCellCache(LevelSetField field, LevelSetCacheMode cacheMode)
{
    // Early return if the cache mode doesn't need to be updated
    if (getFieldCellCacheMode(field) == cacheMode) {
        return;
    }

    // Early return if the cache should be disabled
    if (cacheMode == LevelSetCacheMode::NONE) {
        disableFieldCellCache(field);
        return;
    }

    // Destroy existing cache
    destroyFieldCellCache(field);

    // Update field properties
    std::size_t fieldIndex = static_cast<std::size_t>(field);
    m_cellFieldCacheModes[fieldIndex] = cacheMode;

    // Update data
    if (getKernel()) {
        // Create field cell caches
        createFieldCellCache(field);

        // Fill field cell caches inside the narrow band
        fillFieldCellCaches(LevelSetZone::NARROW_BAND, std::vector<LevelSetField>(1, field));

        // Fill field cell caches inside the bulk
        fillFieldCellCaches(LevelSetZone::BULK, std::vector<LevelSetField>(1, field));
    }
}

/*!
 * Disable the cache for the specified specified field.
 *
 * \param field is the field for which cache will be disabled
 */
void LevelSetObject::disableFieldCellCache(LevelSetField field)
{
    // Early return if no cache is associated with the field
    if (getFieldCellCacheMode(field) == LevelSetCacheMode::NONE) {
        return;
    }

    // Destroy the cache
    destroyFieldCellCache(field);

    // Update field properties
    std::size_t fieldIndex = static_cast<std::size_t>(field);
    m_cellFieldCacheModes[fieldIndex] = LevelSetCacheMode::NONE;
}

/*!
 * Adapt cell cache after a mesh update.
 *
 * Only the transformation listed in the adaption data will be applied, entries associated
 * with new cell will not be filled.
 *
 * \param adaptionData are the information about the mesh update
 */
void LevelSetObject::adaptCellCaches( const std::vector<adaption::Info> &adaptionData )
{
#if BITPIT_ENABLE_MPI==1
    // Get mesh information
    const VolumeKernel &mesh = *(m_kernel->getMesh());
#endif

    // Early return if all the caches are empty
    bool allCachesEmpty = true;
    for (std::size_t cacheId = 0; cacheId < m_cellCacheCollection->size(); ++cacheId) {
        CellCacheCollection::Item &cacheItem = (*m_cellCacheCollection)[cacheId];
        if (cacheItem.hasCache()) {
            allCachesEmpty = false;
            break;
        }
    }

#if BITPIT_ENABLE_MPI==1
    if (mesh.isPartitioned()) {
        MPI_Allreduce(MPI_IN_PLACE, &allCachesEmpty, 1, MPI_C_BOOL, MPI_LAND, m_kernel->getCommunicator()) ;
    }
#endif

    if (allCachesEmpty) {
        return;
    }

#if BITPIT_ENABLE_MPI
    // Identify the ids of the non-volatile caches
    std::vector<std::size_t> nonVolatileCellCacheIds;

    CellCacheCollection::Cache *locationCache = getCellCache(m_cellLocationCacheId);
    if (locationCache && !locationCache->isVolatile()) {
        nonVolatileCellCacheIds.push_back(m_cellLocationCacheId);
    }

    for (std::size_t fieldIndex = 0; fieldIndex < static_cast<std::size_t>(LevelSetField::COUNT); ++fieldIndex) {
        LevelSetField field = static_cast<LevelSetField>(fieldIndex);
        CellCacheCollection::Cache *fieldCache = getFieldCellCache(field);
        if (fieldCache && !fieldCache->isVolatile()) {
            std::size_t fieldCacheId = getFieldCellCacheId(field);
            nonVolatileCellCacheIds.push_back(fieldCacheId);
        }
    }

    // Initialize data exchange
    //
    // Data exchange can be performed only for non-volatile caches.
    std::unordered_map<int, std::vector<long>> exchangeSendList ;
    std::unordered_map<int, std::vector<long>> exchangeRecvList ;
    if (!nonVolatileCellCacheIds.empty()) {
        if (m_kernel->getMesh()->isPartitioned()) {
            for( const adaption::Info &adaptionInfo : adaptionData){
                if( adaptionInfo.entity != adaption::Entity::ENTITY_CELL){
                    continue;
                }

                switch (adaptionInfo.type) {

                case adaption::Type::TYPE_PARTITION_SEND:
                    exchangeSendList.insert({{adaptionInfo.rank,adaptionInfo.previous}}) ;
                    break;

                case adaption::Type::TYPE_PARTITION_RECV:
                    exchangeRecvList.insert({{adaptionInfo.rank,adaptionInfo.current}}) ;
                    break;

                default:
                    break;

                }
            }
        }
    }

    bool exchangeCacheData = (!exchangeSendList.empty() || !exchangeRecvList.empty());
    if (m_kernel->getMesh()->isPartitioned()) {
        MPI_Allreduce(MPI_IN_PLACE, &exchangeCacheData, 1, MPI_C_BOOL, MPI_LOR, m_kernel->getCommunicator()) ;
    }

    // Create data communicator
    std::unique_ptr<DataCommunicator> dataCommunicator;
    if (exchangeCacheData) {
        dataCommunicator = m_kernel->createDataCommunicator() ;
    }
#endif

#if BITPIT_ENABLE_MPI
    // Start data exchange
    if (exchangeCacheData) {
        startCellCachesExchange(exchangeSendList, nonVolatileCellCacheIds, dataCommunicator.get());
    }
#endif

    // Remove stale entries from the caches
    std::vector<long> staleCellIds = evalCellCacheStaleIds(adaptionData);

    for (std::size_t cacheId = 0; cacheId < m_cellCacheCollection->size(); ++cacheId) {
        CellCacheCollection::Item &cacheItem = (*m_cellCacheCollection)[cacheId];
        if (!cacheItem.hasCache()) {
            continue;
        }

        pruneCellCache(cacheId, staleCellIds);
    }

#if BITPIT_ENABLE_MPI
    // Complete data exchange
    if (exchangeCacheData) {
        completeCellCachesExchange(exchangeRecvList, nonVolatileCellCacheIds, dataCommunicator.get());
    }
#endif
}

/*!
 * Clear the specified cell cache.
 *
 * \param cacheId is the id of the cached that will be cleared
 * \param release if set to true the memory hold by the caches will be released, otherwise
 * the caches will be cleared but its memory may not be released
 */
void LevelSetObject::clearCellCache( std::size_t cacheId, bool release )
{
    CellCacheCollection::Item &cacheItem = (*m_cellCacheCollection)[cacheId];

    if (release) {
        cacheItem.destroyCache();
    } else if (cacheItem.hasCache()) {
        CellCacheCollection::Cache *cache = cacheItem.getCache();
        cache->clear();
    }
}

/*!
 * Remove the values associated with the given cell ids from the specified cache.
 *
 * \param cacheId is the id of the cache
 * \param cellIds are the ids of the cells whose values will be removed from the cache
 */
void LevelSetObject::pruneCellCache(std::size_t cacheId, const std::vector<long> &cellIds)
{
    // Get cache
    CellCacheCollection::Item &cacheItem = (*m_cellCacheCollection)[cacheId];
    if (!cacheItem.hasCache()) {
        return;
    }

    CellCacheCollection::Cache *cache = cacheItem.getCache();

    // Remove entries
    cache->erase(cellIds);

    // Reclaim unused memory
    if (m_kernel->getExpectedFillIn() == LevelSetFillIn::SPARSE) {
        cache->shrink_to_fit();
    }
}

/*!
 * Identify the cells that should be inserted in a cache operating in the specified mode.
 *
 * \param zone is the zone for which the fillable cells are requested
 * \param cacheMode is the cache mode for which the fillable cells are requested
 * \result The ids of the cells that should be inserted in a cache operating in the specified mode.
 */
std::vector<long> LevelSetObject::evalCellCacheFillIds(LevelSetZone zone, LevelSetCacheMode cacheMode) const
{
    switch (cacheMode) {

    case LevelSetCacheMode::ON_DEMAND:
        return evalCellOnDemandCacheFillIds(zone);

    case LevelSetCacheMode::NARROW_BAND:
        return evalCellNarrowBandCacheFillIds(zone);

    case LevelSetCacheMode::FULL:
        return evalCellFullCacheFillIds(zone);

    default:
        throw std::runtime_error("Unsupported cache mode!");

    }
}

/*!
 * Identify the newly added cells that should be inserted after a mesh update in a cache operating
 * in the specified mode.
 *
 * \param zone is the zone for which the fillable cells are requested
 * \param cacheMode is the cache mode for which the fillable cells are requested
 * \param adaptionData are the information about the mesh update
 * \result The ids of newly added cells that should be inserted after a mesh update in a cache
 * operating in the specified mode.
 */
std::vector<long> LevelSetObject::evalCellCacheFillIds(LevelSetZone zone, LevelSetCacheMode cacheMode,
                                                       const std::vector<adaption::Info> &adaptionData) const
{
    switch (cacheMode) {

    case LevelSetCacheMode::ON_DEMAND:
        return evalCellOnDemandCacheFillIds(zone, adaptionData);

    case LevelSetCacheMode::NARROW_BAND:
        return evalCellNarrowBandCacheFillIds(zone, adaptionData);

    case LevelSetCacheMode::FULL:
        return evalCellFullCacheFillIds(zone, adaptionData);

    default:
        throw std::runtime_error("Unsupported cache mode!");

    }
}

/*!
 * Identify the cells that should be inserted in a cache operating in "on-demand" mode.
 *
 * No cells should be automatically added to caches operating in "on-demand" mode.
 *
 * \param zone is the zone for which the fillable cells are requested
 * \result The ids of cells that should be inserted in a cache operating in "on-demand" mode.
 */
std::vector<long> LevelSetObject::evalCellOnDemandCacheFillIds(LevelSetZone zone) const
{
    BITPIT_UNUSED(zone);

    return std::vector<long>();
}

/*!
 * Identify the newly added cells that should be inserted after a mesh update in a cache operating
 * in "on-demand" mode.
 *
 * No cells should be automatically added to caches operating in "on-demand" mode.
 *
 * \param zone is the zone for which the fillable cells are requested
 * \param adaptionData are the information about the mesh update
 * \result The ids of the cells that should be inserted after a mesh update in a cache operating
 * in "on-demand" mode.
 */
std::vector<long> LevelSetObject::evalCellOnDemandCacheFillIds(LevelSetZone zone, const std::vector<adaption::Info> &adaptionData) const
{
    BITPIT_UNUSED(zone);
    BITPIT_UNUSED(adaptionData);

    return std::vector<long>();
}

/*!
 * Identify the cells that should be inserted in a cache operating in "narrow band" mode.
 *
 * Cache operating in "narrow band" mode should contain the cells inside the narrow band.
 *
 * \param zone is the zone for which the fillable cells are requested
 * \result The ids of cells that should be inserted in a cache operating in "narrow band" mode.
 */
std::vector<long> LevelSetObject::evalCellNarrowBandCacheFillIds(LevelSetZone zone) const
{
    // Early return if the object is empty
    if (empty()) {
        return std::vector<long>();
    }

    // There are no narrow band cells in the bulk
    if (zone == LevelSetZone::BULK) {
        return std::vector<long>();
    }

    // Get mesh information
    const VolumeKernel &mesh = *(m_kernel->getMesh()) ;

    // Detect the cells to be filled
    //
    // If no narrow band size was set, all cells should be filled, otherwise only
    // the cells inside the narrow band should be filled.
    if (m_narrowBandSize == LEVELSET_NARROW_BAND_UNLIMITED) {
        if (const VolCartesian *cartesianMesh = dynamic_cast<const VolCartesian *>(&mesh)) {
            long nCells = cartesianMesh->getCellCount();
            std::vector<long> cellFillIds(nCells);
            std::iota(cellFillIds.begin(), cellFillIds.end(), 0);

            return cellFillIds;
        } else {
            return mesh.getCells().getIds(false);
        }
    } else {
        std::vector<long> cellFillIds;
        if (const VolCartesian *cartesianMesh = dynamic_cast<const VolCartesian *>(&mesh)) {
            long nCells = cartesianMesh->getCellCount();
            for (long cellId = 0; cellId < nCells; ++cellId) {
                if (isCellInNarrowBand(cellId)) {
                    cellFillIds.push_back(cellId);
                }
            }
        } else {
            for (const Cell &cell : mesh.getCells()) {
                long cellId = cell.getId();
                if (isCellInNarrowBand(cellId)) {
                    cellFillIds.push_back(cellId);
                }
            }
        }

        return cellFillIds;
    }
}

/*!
 * Identify the newly added cells that should be inserted after a mesh update in a cache operating
 * in "narrow band" mode.
 *
 * Cache operating in "narrow band" mode should contain the cells inside the narrow band.
 *
 * \param zone is the zone for which the fillable cells are requested
 * \param adaptionData are the information about the mesh update
 * \result The ids of the cells that should be inserted after a mesh update in a cache operating
 * in "narrow band" mode.
 */
std::vector<long> LevelSetObject::evalCellNarrowBandCacheFillIds(LevelSetZone zone, const std::vector<adaption::Info> &adaptionData) const
{
    // Early return if the object is empty
    if (empty()) {
        return std::vector<long>();
    }

    // There are no narrow band cells in the bulk
    if (zone == LevelSetZone::BULK) {
        return std::vector<long>();
    }

    // Detect the cells to be filled
    //
    // If no narrow band size was set, all new cells should be filled, otherwise only
    // the cells inside the narrow band should be filled.
    std::vector<long> cellFillIds;
    for (const adaption::Info &adaptionInfo : adaptionData) {
        if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
            continue;
        }

        if (m_narrowBandSize == LEVELSET_NARROW_BAND_UNLIMITED) {
            cellFillIds.insert(cellFillIds.end(), adaptionInfo.current.begin(), adaptionInfo.current.end());
        } else {
            for (long cellId : adaptionInfo.current) {
                if (isCellInNarrowBand(cellId)) {
                    cellFillIds.push_back(cellId);
                }
            }
        }
    }

    return cellFillIds;
}

/*!
 * Identify the cells that should be inserted in a cache operating in "full" mode.
 *
 * Cache operating in "full" mode should contain all the cells inside the requested zone.
 *
 * \param zone is the zone for which the fillable cells are requested
 * \result The ids of cells that should be inserted in a cache operating in "full" mode.
 */
std::vector<long> LevelSetObject::evalCellFullCacheFillIds(LevelSetZone zone) const
{
    // Early return if the object is empty
    if (empty()) {
        return std::vector<long>();
    }

    // Detect the cells to be filled
    //
    // If the requested zone is the narrow band, it is possible to call the same function used
    // for the caches in narrow band mode.
    //
    // If the requested zone is the bulk, only the functions outside the narrow band should be
    // filled. This implies that, if the narrow band size was not not set, no cells should be
    // filled.
    if (zone == LevelSetZone::NARROW_BAND) {
        return evalCellNarrowBandCacheFillIds(zone);
    } else {
        if (m_narrowBandSize == LEVELSET_NARROW_BAND_UNLIMITED) {
            return std::vector<long>();
        } else {
            // Get mesh information
            const VolumeKernel &mesh = *(m_kernel->getMesh()) ;

            // Identify cells outside the narrow band
            std::vector<long> cellFillIds;
            if (const VolCartesian *cartesianMesh = dynamic_cast<const VolCartesian *>(&mesh)) {
                long nCells = cartesianMesh->getCellCount();
                for (long cellId = 0; cellId < nCells; ++cellId) {
                    if (!isCellInNarrowBand(cellId)) {
                        cellFillIds.push_back(cellId);
                    }
                }
            } else {
                for (const Cell &cell : mesh.getCells()) {
                    long cellId = cell.getId();
                    if (!isCellInNarrowBand(cellId)) {
                        cellFillIds.push_back(cellId);
                    }
                }
            }

            return cellFillIds;
        }
    }
}

/*!
 * Identify the newly added cells that should be inserted after a mesh update in a cache operating
 * in "full" mode.
 *
 * Cache operating in "full" mode should contain all the cells inside the requested zone.
 *
 * \param zone is the zone for which the fillable cells are requested
 * \param adaptionData are the information about the mesh update
 * \result The ids of the cells that should be inserted after a mesh update in a cache operating
 * in "full" mode.
 */
std::vector<long> LevelSetObject::evalCellFullCacheFillIds(LevelSetZone zone, const std::vector<adaption::Info> &adaptionData) const
{
    // Early return if the object is empty
    if (empty()) {
        return std::vector<long>();
    }

    // Detect the cells to be filled
    //
    // If the requested zone is the narrow band, it is possible to call the same function used
    // for the caches in narrow band mode.
    //
    // If the requested zone is the bulk, only the functions outside the narrow band should be
    // filled. This implies that, if the narrow band size was not not set, no cells should be
    // filled.
    if (zone == LevelSetZone::NARROW_BAND) {
        return evalCellNarrowBandCacheFillIds(zone);
    } else {
        if (m_narrowBandSize == LEVELSET_NARROW_BAND_UNLIMITED) {
            return std::vector<long>();
        } else {
            // Identify cells outside the narrow band
            std::vector<long> cellFillIds;
            for (const adaption::Info &adaptionInfo : adaptionData) {
                if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
                    continue;
                }

                for (long cellId : adaptionInfo.current) {
                    if (!isCellInNarrowBand(cellId)) {
                        cellFillIds.push_back(cellId);
                    }
                }
            }

            return cellFillIds;
        }
    }
}

/*!
 * Identify the stale cell ids that should be removed from the cache after a mesh update.
 *
 * \param adaptionData are the information about the mesh update
 * \result The stale cell ids that should be removed from the cache after a mesh update.
 */
std::vector<long> LevelSetObject::evalCellCacheStaleIds(const std::vector<adaption::Info> &adaptionData) const
{
    // Identify cells whose entries should be pruned
    //
    // We need to prune entries for both cells that have been removed and newly added cells.
    // When a new cell is added to the cache, synchronized caches will automatically create
    // an entry for the new cell, however the information associated to that entry is not
    // initialized and may contain random data. We need to explicitly state that the newly
    // added cells are not in the cached yet.
    std::vector<long> cellPruneIds;
    for (const adaption::Info &adaptionInfo : adaptionData) {
        if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
            continue;
        }

        cellPruneIds.insert(cellPruneIds.end(), adaptionInfo.previous.begin(), adaptionInfo.previous.end());
        cellPruneIds.insert(cellPruneIds.end(), adaptionInfo.current.begin(), adaptionInfo.current.end());
    }

    return cellPruneIds;
}

/*!
 * Get a pointer to the cell cache for the specified field.
 *
 * If no cache was registered for the specified field, a null pointer is returned.
 *
 * \param field is the field whose cache will be retrieved
 * \result A pointer to the cell cache for the specified field.
 */
typename LevelSetObject::CellCacheCollection::Cache * LevelSetObject::getFieldCellCache(LevelSetField field) const
{
    // Early return if no cache is associated with the field
    if (getFieldCellCacheMode(field) == LevelSetCacheMode::NONE) {
        return nullptr;
    }

    // Get field cache
    std::size_t cacheId = getFieldCellCacheId(field);

    return getCellCache(cacheId);
}

/*!
 * Get the id of the cache associated with the specified field.
 *
 * \param field is the specified field
 * \result The id of the cache associated with the specified field.
 */
std::size_t LevelSetObject::getFieldCellCacheId(LevelSetField field) const
{
    std::size_t fieldIndex = static_cast<std::size_t>(field);

    return m_cellFieldCacheIds[fieldIndex];
}

/*!
 * Create the cache that will be used for storing cell information of the specified field.
 *
 * \param field is the field for which the caches will be added
 * \param cacheId is the id that will be associated with the cache, if a NULL_ID is specified
 * the cache id will be assigned automatically
 * \result The id associated with the cache.
 */
std::size_t LevelSetObject::createFieldCellCache(LevelSetField field, std::size_t cacheId)
{
    switch(field) {

    case LevelSetField::VALUE:
        return createFieldCellCache<double>(field, cacheId);

    case LevelSetField::SIGN:
        return createFieldCellCache<short>(field, cacheId);

    case LevelSetField::GRADIENT:
        return createFieldCellCache<std::array<double, 3>>(field, cacheId);

    default:
        throw std::runtime_error("The requested field is not supported!");

    }
}

/*!
 * Unregister the specified field cache.
 *
 * \param field is the field for which the caches will be added
 */
void LevelSetObject::destroyFieldCellCache(LevelSetField field)
{
    // Update field properties
    std::size_t fieldIndex = static_cast<std::size_t>(field);;
    m_cellFieldCacheIds[fieldIndex] = CellCacheCollection::NULL_CACHE_ID;

    // Destroy cache
    std::size_t cacheId = getFieldCellCacheId(field);
    destroyCellCache(cacheId);
}

/*!
 * Fill the cell caches of the specified fields.
 *
 * \param zone is the zone where the cell caches will be filled
 * \param fields are the fields whose caches will be filled
 */
void LevelSetObject::fillFieldCellCaches(LevelSetZone zone, const std::vector<LevelSetField> &fields)
{
    // Early return if there are no caches to update
    if (fields.empty()) {
        return;
    }

    // It's more efficient to process the fields grouped by their cache mode
    for (std::size_t cacheModeIndex = 0; cacheModeIndex < static_cast<std::size_t>(LevelSetCacheMode::COUNT); ++cacheModeIndex) {
        LevelSetCacheMode cacheMode = static_cast<LevelSetCacheMode>(cacheModeIndex);

        // Get fields that operate in the processed cache mode
        std::vector<LevelSetField> cacheModeFields;
        for (LevelSetField field : fields) {
            if (getFieldCellCacheMode(field) != cacheMode) {
                continue;
            }

            cacheModeFields.push_back(field);
        }

        if (cacheModeFields.empty()) {
            continue;
        }

        // Identify the cell ids that should be added to the cache
        std::vector<long> cellCacheFillIds = evalCellCacheFillIds(zone, cacheMode);

        // Fill the caches
        for (LevelSetField field : cacheModeFields) {
            fillFieldCellCache(field, cellCacheFillIds);
        }
    }
}

/*!
 * Fill the cell caches of the specified fields after a mesh update.
 *
 * \param zone is the zone where the cell caches will be filled
 * \param fields are the fields whose caches will be filled
 * \param adaptionData are the information about the mesh update
 */
void LevelSetObject::fillFieldCellCaches(LevelSetZone zone, const std::vector<LevelSetField> &fields,
                                         const std::vector<adaption::Info> &adaptionData)
{
    // Early return if there are no caches to update
    if (fields.empty()) {
        return;
    }

#if BITPIT_ENABLE_MPI==1
    // Get mesh information
    const VolumeKernel &mesh = *(m_kernel->getMesh());

    // Check if there are non-volatile caches to process
    bool areNonVolatileCacheProcessed = false;
    for (LevelSetField field : fields) {
        CellCacheCollection::Cache *fieldCache = getFieldCellCache(field);
        if (fieldCache && !fieldCache->isVolatile()) {
            areNonVolatileCacheProcessed = true;
            break;
        }
    }

    // Identify the cells that have been received.
    std::unordered_set<long> recievedCellIds;
    if (areNonVolatileCacheProcessed) {
        if (m_kernel->getMesh()->isPartitioned()) {
            for( const adaption::Info &adaptionInfo : adaptionData){
                if( adaptionInfo.entity != adaption::Entity::ENTITY_CELL){
                    continue;
                }

                switch (adaptionInfo.type) {

                case adaption::Type::TYPE_PARTITION_RECV:
                    recievedCellIds.insert(adaptionInfo.current.begin(), adaptionInfo.current.end());
                    break;

                default:
                    break;

                }
            }
        }
    }
#endif

    // It's more efficient to process the fields grouped by their cache mode
    for (std::size_t cacheModeIndex = 0; cacheModeIndex < static_cast<std::size_t>(LevelSetCacheMode::COUNT); ++cacheModeIndex) {
        LevelSetCacheMode cacheMode = static_cast<LevelSetCacheMode>(cacheModeIndex);

        // Get fields with caches that operate in the processed mode
        std::vector<LevelSetField> cacheModeFields;
        for (LevelSetField field : fields) {
            if (getFieldCellCacheMode(field) != cacheMode) {
                continue;
            }

            cacheModeFields.push_back(field);
        }

        if (cacheModeFields.empty()) {
            continue;
        }

        // Identify the cell ids that should be filled for volatile caches
        std::vector<long> cellVolatileCacheFillIds = evalCellCacheFillIds(zone, cacheMode, adaptionData);

        bool emptyCellVolatileCacheFillIds = cellVolatileCacheFillIds.empty();
#if BITPIT_ENABLE_MPI==1
        if (mesh.isPartitioned()) {
            MPI_Allreduce(MPI_IN_PLACE, &emptyCellVolatileCacheFillIds, 1, MPI_C_BOOL, MPI_LAND, m_kernel->getCommunicator()) ;
        }
#endif

        if (emptyCellVolatileCacheFillIds) {
            continue;
        }

        // Identify the cell ids that should be filled for non-volatile caches
        //
        // For non-volatile caches, entries associated with cell received from other processes
        // are exchanged when the caches are adapted and there is no need to evaluate them.
        std::vector<long> cellNonVolatileCacheFillIds;
#if BITPIT_ENABLE_MPI
        if (areNonVolatileCacheProcessed && !recievedCellIds.empty()) {
            cellNonVolatileCacheFillIds.reserve(cellVolatileCacheFillIds.size());
            for (long cellId : cellVolatileCacheFillIds) {
                if (recievedCellIds.count(cellId)) {
                    continue;
                }

                cellNonVolatileCacheFillIds.push_back(cellId);
            }
        } else {
            cellNonVolatileCacheFillIds = cellVolatileCacheFillIds;
        }
#else
        cellNonVolatileCacheFillIds = cellVolatileCacheFillIds;
#endif

        // Fill field caches
        for (LevelSetField field : cacheModeFields) {
            CellCacheCollection::Cache *cache = getFieldCellCache(field);
            if (cache->isVolatile()) {
                fillFieldCellCache(field, cellVolatileCacheFillIds);
            } else {
                fillFieldCellCache(field, cellNonVolatileCacheFillIds);
            }
        }
    }
}

/*!
 * Fill the cache values associated with the given cell ids for the specified field.
 *
 * \param field is the field whose cache will be filled
 * \param cellIds are the ids of the cells whose values will be filled
 */
void LevelSetObject::fillFieldCellCache(LevelSetField field, const std::vector<long> &cellIds)
{
    // Early return if no cache is associated with the field
    if (getFieldCellCacheMode(field) == LevelSetCacheMode::NONE) {
        return;
    }

#if BITPIT_ENABLE_MPI==1
    // Get mesh information
    const VolumeKernel &mesh = *(m_kernel->getMesh());
#endif

    // Early return if there are no cells to fill
    //
    // Even if there are no cells to fill, the cache will now be considered non-empty, in this
    // way it will be processed by the functions that update the caches.
    bool emptyCellIds = cellIds.empty();
#if BITPIT_ENABLE_MPI==1
    if (mesh.isPartitioned()) {
        MPI_Allreduce(MPI_IN_PLACE, &emptyCellIds, 1, MPI_C_BOOL, MPI_LAND, m_kernel->getCommunicator()) ;
    }
#endif

    if (emptyCellIds) {
        return;
    }

    // Get list of ghost exchange sources
    std::unordered_set<long> ghostExchangeSources;
    std::unordered_set<long> ghostExchangeTargets;
#if BITPIT_ENABLE_MPI
    for (const auto &entry : mesh.getGhostCellExchangeSources()) {
        const std::vector<long> &rankSourceList = entry.second;
        ghostExchangeSources.insert(rankSourceList.begin(), rankSourceList.end());
    }

    for (const auto &entry : mesh.getGhostCellExchangeTargets()) {
        const std::vector<long> &rankTargetList = entry.second;
        ghostExchangeTargets.insert(rankTargetList.begin(), rankTargetList.end());
    }
#endif

#if BITPIT_ENABLE_MPI
    // Get field cache
    std::size_t cacheId = getFieldCellCacheId(field);

    // Create data communicator
    std::unique_ptr<DataCommunicator> dataCommunicator;
    if (mesh.isPartitioned()) {
        dataCommunicator = m_kernel->createDataCommunicator();
    }

    // Fill cache values of cells that are sources for ghost exchange
    for (long cellId : cellIds) {
        if (ghostExchangeSources.count(cellId) == 0) {
            continue;
        }

        fillFieldCellCache(field, cellId);
    }

    // Start ghost exchange
    if (dataCommunicator) {
        startCellCacheExchange(mesh.getGhostCellExchangeSources(), cacheId, dataCommunicator.get());
    }
#endif

    // Fill cache values of internal cells that are not sources for ghost exchange
    //
    // Values on ghost cells are not evaluated, because those values will be communicated
    // by the ghost data exchange.
    for (long cellId : cellIds) {
        if (ghostExchangeSources.count(cellId) > 0) {
            continue;
        } else if (ghostExchangeTargets.count(cellId) > 0) {
            continue;
        }

        fillFieldCellCache(field, cellId);
    }

#if BITPIT_ENABLE_MPI
    // Complete ghost exchange
    if (dataCommunicator) {
        completeCellCacheExchange(mesh.getGhostCellExchangeTargets(), cacheId, dataCommunicator.get());
    }
#endif
}

/*!
 * Fill the specified field cache of the given cell.
 *
 * \param field is the field whose cache will be filled
 * \param id is the id of the cell whose cache will be filled
 */
void LevelSetObject::fillFieldCellCache(LevelSetField field, long id)
{
    // Early return if no cache is associated with the field
    if (getFieldCellCacheMode(field) == LevelSetCacheMode::NONE) {
        return;
    }

    // Fill cache
    switch (field) {

    case LevelSetField::SIGN:
        evalCellSign(id);
        break;

    case LevelSetField::VALUE:
        evalCellValue(id, true);
        break;

    case LevelSetField::GRADIENT:
        evalCellGradient(id, true);
        break;

    default:
        throw std::runtime_error("Unsupported field!");

    }
}

/*!
 * Get a pointer to the specified cell cache.
 *
 * If specified cell cache was not registered or if an invalid cache id is specified, a null
 * pointer is returned.
 *
 * \param cacheId the id of the cell that will be unregistered
 * \result A pointer to the specified cell cache.
 */
typename LevelSetObject::CellCacheCollection::Cache * LevelSetObject::getCellCache(std::size_t cacheId) const
{
    if (cacheId == CellCacheCollection::NULL_CACHE_ID) {
        return nullptr;
    }

    return (*m_cellCacheCollection)[cacheId].getCache();
}

/*!
 * Destroy the specified cache.
 *
 * \param cacheId the id of the cell that will be unregistered
 */
void LevelSetObject::destroyCellCache(std::size_t cacheId)
{
    m_cellCacheCollection->erase(cacheId);
}

/*!
 * Computes the projection point of the cell center, i.e. the closest
 * point to the cell center on the zero level set
 * @param[in] cellId cell id
 * @return the projection point
 */
std::array<double,3> LevelSetObject::computeProjectionPoint(long cellId) const {
    return evalCellProjectionPoint(cellId);
}

/*!
 * Computes the projection point of the vertex, i.e. the closest point to the
 * vertex on the zero level set
 * @param[in] vertexId vertex id
 * @return the projection point
 */
std::array<double,3> LevelSetObject::computeVertexProjectionPoint(long vertexId) const {
    const std::array<double,3> &vertexCoords = m_kernel->getMesh()->getVertexCoords(vertexId);
    return evalProjectionPoint(vertexCoords);
}

/*!
 * Projects a vertex on the zero levelset
 * @param[in] point point coordinates
 * @return the projected point
 */
std::array<double,3> LevelSetObject::computeProjectionPoint(const std::array<double,3> &point) const{
    return evalProjectionPoint(point);
}

/*!
 * Get the sign of the levelset function
 * @param[in] cellId cell id
 * @return sign of levelset
 */
short LevelSetObject::getSign(long cellId) const {

    return evalCellSign(cellId);

}

/*!
 * Get the levelset value of cell
 * @param[in] id cell id
 * @return levelset value in cell
 */
double LevelSetObject::getValue(long cellId) const {

    return evalCellValue(cellId, m_defaultSignedLevelSet);

}

/*!
 * Get the levelset gradient of cell
 * @param[in] cellId cell id
 * @return levelset gradient in cell
 */
std::array<double,3> LevelSetObject::getGradient(long cellId) const {

    return evalCellGradient(cellId, m_defaultSignedLevelSet);

}

/*!
 * Get LevelSetInfo of cell
 * @param[in] cellId cell idex
 * @return LevelSetInfo of cell
*/
LevelSetInfo LevelSetObject::getLevelSetInfo(long cellId) const {

    return LevelSetInfo(evalCellValue(cellId, m_defaultSignedLevelSet), evalCellGradient(cellId, m_defaultSignedLevelSet));

}

/*!
 * Get the levelset value of cell
 * @param[in] cellId cell id
 * @return levelset value in cell
 */
double LevelSetObject::getLS(long cellId) const {

    return evalCellValue(cellId, m_defaultSignedLevelSet);

}

/*!
 * Get the current size of the narrow band.
 * The function will always return an "infinite" distance.
 * @return size of the current narrow band
 */
double LevelSetObject::getSizeNarrowBand()const{
    return getNarrowBandSize();
}

/*!
 * Manually set the size of the narrow band.
 * The function is a no-op.
 * @param[in] r size of the narrow band.
 */
void LevelSetObject::setSizeNarrowBand(double r){
    setNarrowBandSize(r);
}

}
