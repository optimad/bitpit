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

# include <vector>

# include "levelSetObject.hpp"

namespace bitpit {

/*!
	@interface LevelSetObject
	@ingroup levelset
	@brief Interface class for all objects with respect to whom the levelset function may be computed.
*/

/*!
 * Constructor
 * @param[in] id id assigned to object
 */
LevelSetObject::LevelSetObject(int id)
    : m_kernel(nullptr),
      m_narrowBandSize(levelSetDefaults::NARROWBAND_SIZE),
      m_cellNarrowBandCacheId(CellCacheCollection::NULL_ID),
      m_defaultSignedLevelSet(false),
      m_nReferences(0),
      m_cellFieldCacheIds(static_cast<std::size_t>(LevelSetField::COUNT), CellCacheCollection::NULL_ID),
      m_cellFieldCacheModes(static_cast<std::size_t>(LevelSetField::COUNT), LevelSetCacheMode::NONE)
{
    setId(id);
}

/*!
 * Copy constructor
 * @param[in] other is another object whose content is copied in this object
 */
LevelSetObject::LevelSetObject(const LevelSetObject &other)
    : m_kernel(other.m_kernel),
      m_narrowBandSize(other.m_narrowBandSize),
      m_cellNarrowBandCacheId(other.m_cellNarrowBandCacheId),
      m_defaultSignedLevelSet(other.m_defaultSignedLevelSet),
      m_enabledOutputFields(other.m_enabledOutputFields),
      m_id(other.m_id),
      m_nReferences(other.m_nReferences),
      m_cellCacheCollection(new CellCacheCollection(*(other.m_cellCacheCollection))),
      m_cellFieldCacheIds(other.m_cellFieldCacheIds),
      m_cellFieldCacheModes(other.m_cellFieldCacheModes)

{
    for ( const auto &fieldEntry : m_enabledOutputFields ) {
        enableVTKOutput(fieldEntry.first, true);
    }
}

/*!
 * Move constructor
 * @param[in] other is another object whose content is copied in this object
 */
LevelSetObject::LevelSetObject(LevelSetObject &&other)
    : m_kernel(std::move(other.m_kernel)),
      m_narrowBandSize(std::move(other.m_narrowBandSize)),
      m_cellNarrowBandCacheId(std::move(other.m_cellNarrowBandCacheId)),
      m_defaultSignedLevelSet(std::move(other.m_defaultSignedLevelSet)),
      m_id(std::move(other.m_id)),
      m_nReferences(std::move(other.m_nReferences)),
      m_cellCacheCollection(std::move(other.m_cellCacheCollection)),
      m_cellFieldCacheIds(std::move(other.m_cellFieldCacheIds)),
      m_cellFieldCacheModes(std::move(other.m_cellFieldCacheModes))
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
    // Disable output
    LevelSetFieldset enabledOutputFieldset;
    for ( const auto &fieldEntry : m_enabledOutputFields ) {
        enabledOutputFieldset.insert(fieldEntry.first);
    }

    enableVTKOutput(enabledOutputFieldset, false);
}

/*!
 * Get the list of supported field.
 * @result The list of supported field.
 */
LevelSetFieldset LevelSetObject::getSupportedFields() const {

    LevelSetFieldset supportedFields;
    supportedFields.insert(LevelSetField::SIGN);
    supportedFields.insert(LevelSetField::VALUE);
    supportedFields.insert(LevelSetField::GRADIENT);

    return supportedFields;

}

/*!
 * Sets the identifier of object
 * @param[in] id is the identifier
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
 * Sets the kernel for the object
 * @param[in] kernel is the LevelSetKernel
 */
void LevelSetObject::setKernel(LevelSetKernel *kernel) {
    // Set kernel
    m_kernel = kernel;

    // Enable output
    for ( const auto &fieldEntry : m_enabledOutputFields ) {
        enableVTKOutput( fieldEntry.first, true ) ;
    }

    // Create cell cache collection
    m_cellCacheCollection = std::unique_ptr<CellCacheCollection>(new CellCacheCollection());
}

/*!
 * Gets a pointer to the kernel for the object
 * @return A pointer to the kernel for the object
 */
LevelSetKernel * LevelSetObject::getKernel() {
    return m_kernel;
}

/*!
 * Gets a constant pointer to the kernel for the object
 * @return A constant pointer to the kernel for the object
 */
const LevelSetKernel * LevelSetObject::getKernel() const {
    return m_kernel;
}

/*!
 * Set whether a signed or unsigned levelset is used when signdness is not explicitly specified.
 * This function is only needed for guarantee backwards compatibility with older versions. In
 * such versions, functions for evaluating levelset information were not taking in input the
 * signdness.
 * @param signedLevelSet controls if signed levelset function will be used
 */
void LevelSetObject::setDefaultLevelSetSigned(bool signedLevelSet) {
    m_defaultSignedLevelSet = signedLevelSet;
}

/*!
 * Get the id 
 * @return id of the object
 */
int LevelSetObject::getId( ) const {
    return m_id ;
}

/*!
 * If the levelset is primary (e.g. of a surface triangulation) or not (e.g. derived by boolean operations between two levelsets)
 * @return if object is primary
 */
bool LevelSetObject::isPrimary( ) const {
    return true;
}

/*!
 * Check if cell intersects the surface
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
LevelSetIntersectionStatus LevelSetObject::intersectSurface(long id, LevelSetIntersectionMode mode) const{

    double absoluteDistance  = evalCellValue(id, false);
    double distanceTolerance = m_kernel->getMesh()->getTol();

    switch(mode){
        case LevelSetIntersectionMode::FAST_GUARANTEE_TRUE:
        {
            double tangentSphere = m_kernel->computeCellTangentRadius(id) ;
            if(utils::DoubleFloatingLessEqual()(absoluteDistance, tangentSphere, distanceTolerance, distanceTolerance)){
                return LevelSetIntersectionStatus::TRUE;
            } else {
                return LevelSetIntersectionStatus::FALSE;
            }

            break;
        }

        case LevelSetIntersectionMode::FAST_GUARANTEE_FALSE:
        {
            double boundingSphere = m_kernel->computeCellBoundingRadius(id) ;
            if(utils::DoubleFloatingGreater()(absoluteDistance, boundingSphere, distanceTolerance, distanceTolerance)){
                return LevelSetIntersectionStatus::FALSE;
            } else {
                return LevelSetIntersectionStatus::TRUE;
            }

            break;
        }

        case LevelSetIntersectionMode::FAST_FUZZY:
        {
            double boundingSphere = m_kernel->computeCellBoundingRadius(id) ;
            if(utils::DoubleFloatingGreater()(absoluteDistance, boundingSphere, distanceTolerance, distanceTolerance)){
                return LevelSetIntersectionStatus::FALSE;
            }

            double tangentSphere = m_kernel->computeCellTangentRadius(id) ;
            if(utils::DoubleFloatingLessEqual()(absoluteDistance, tangentSphere, distanceTolerance, distanceTolerance)){
                return LevelSetIntersectionStatus::TRUE;
            }

            return LevelSetIntersectionStatus::CLOSE;

            break;
        }

        case LevelSetIntersectionMode::ACCURATE:
        {
            double boundingSphere = m_kernel->computeCellBoundingRadius(id) ;
            if(utils::DoubleFloatingGreater()(absoluteDistance, boundingSphere, distanceTolerance, distanceTolerance)){
                return LevelSetIntersectionStatus::FALSE;
            }

            double tangentSphere = m_kernel->computeCellTangentRadius(id) ;
            if(utils::DoubleFloatingLessEqual()(absoluteDistance, tangentSphere, distanceTolerance, distanceTolerance)){
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
 * Get the physical size of the narrow band.
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
 * All caches that use the "narrow band" mode will be destroyed and re-generated.
 *
 * \param size is the physical size of the narrow band
 */
void LevelSetObject::setNarrowBandSize(double size)
{
    // Set narrow band size
    m_narrowBandSize = size;

    // Update caches
    fillNarrowBandCellCaches();
}

/*!
 * Check if the specified cell lies within the narrow band.
 *
 * A cell is considered within the narrow band if one of the following conditions hold:
 *  - its distance from the surface is less than the narrow band size;
 *  - it intersects the zero-levelset iso-surface (intersections are checked using the
 *    FAST_GUARANTEE_FALSE criterium);
 *  - it has at least one negihbour that intersects the zero-levelset iso-surface and its
 *    sign differs form the sign of the neighbour that intersects the surface.
 *
 * If no caches with "narrow band" mode have been filled, the function may return wrong
 * results if the cell is on the last layer of ghosts.
 *
 * \param[in] id is the cell id
 * \result Return true if the cell is in the narrow band, false otherwise.
 */
bool LevelSetObject::isCellInNarrowBand(long id)const
{
    if (m_cellNarrowBandCacheId != CellCacheCollection::NULL_ID) {
        CellValueCache<bool> *narrowBandCache = getCellCache<bool>(m_cellNarrowBandCacheId);
        CellValueCache<bool>::Entry narrowBandCacheEntry = narrowBandCache->findEntry(id);
        if (narrowBandCacheEntry.isValid()) {
            return *narrowBandCacheEntry;
        } else {
            return false;
        }
    }

    return _isCellInNarrowBand(id, true);
}

/*!
 * Internal function to check if the specified cell lies within the narrow band.
 *
 * A cell is considered within the narrow band if one of the following conditions hold:
 *  - its distance from the surface is less than the narrow band size;
 *  - it intersects the zero-levelset iso-surface (intersections are checked using the
 *    FAST_GUARANTEE_FALSE criterium);
 *  - it has at least one negihbour that intersects the zero-levelset iso-surface and its
 *    sign differs form the sign of the neighbour that intersects the surface.
 *
 * Neighbour check is not reliable if the cell is on the last layer of ghosts.
 *
 * \param[in] id is the cell id
 * \param[in] checkNeighbours is set to true, neighbours are check to detect if the cell
 * should be added to the levelset because it has at least one negihbour that intersects
 * the zero-levelset iso-surface and its sign differs form the sign of the neighbour that
 * intersects the surface
 * \param[out] maximumDistance if a valid pointer is provided and the cell is inside the
 * narrow band, on output will contain a conservative estimate for the distance of the
 * cell from the surface, the distance of the cell from the surface will always be less
 * or equal than the provided estimate
 * \result Return true if the cell is in the narrow band, false otherwise.
 */
bool LevelSetObject::_isCellInNarrowBand(long id, bool checkNeighbours, double *maximumDistance) const
{
    // Check if the distance from the surface is less than the narrow band size
    double value = evalCellValue(id, false);
    if (utils::DoubleFloatingLessEqual()(value, m_narrowBandSize)) {
        if (maximumDistance) {
            *maximumDistance = value;
        }

        return true;
    }

    // Check if the cell intersects the surface
    //
    // Cells that intersect the surface should always be included in the narrow band, even
    // if their distance from the surface is greater than then narrow band size.
    if (intersectSurface(id, LevelSetIntersectionMode::FAST_GUARANTEE_FALSE) == LevelSetIntersectionStatus::TRUE) {
        if (maximumDistance) {
            *maximumDistance = value;
        }

        return true;
    }

    // Process cells with neighbours that intersect the zero-levelset iso-surface
    //
    // Cells with at least a negihbour that intersect the zero-levelset iso-surface need to be
    // added to the narrow band if their sign differs form the sign of the neighbour that
    // intersects the surface.
    if (checkNeighbours) {
        const VolumeKernel *mesh = m_kernel->getMesh();
        bool cellAdjacenciesAvailable = (mesh->getAdjacenciesBuildStrategy() != VolumeKernel::ADJACENCIES_NONE);

        const long *cellNeighs;
        int nCellNeighs;
        if (cellAdjacenciesAvailable) {
            const Cell &cell = mesh->getCell(id);
            cellNeighs = cell.getAdjacencies();
            nCellNeighs = cell.getAdjacencyCount();
        } else {
            std::vector<long> cellNeighStorage;
            mesh->findCellFaceNeighs(id, &cellNeighStorage);
            cellNeighs = cellNeighStorage.data();
            nCellNeighs = cellNeighStorage.size();
        }

        for(int n = 0; n < nCellNeighs; ++n){
            long neighId = cellNeighs[n];

            // Skip neighbours with the same sign
            long neighSign = evalCellSign(neighId);
            if (neighId == neighSign) {
                continue;
            }

            // Skip neighbours that doesn't intersect the surface
            if (intersectSurface(neighId, LevelSetIntersectionMode::FAST_GUARANTEE_FALSE) != LevelSetIntersectionStatus::TRUE){
                continue;
            }

            // Cell is inside the narrow band
            //
            // The cell has a neighbour with opposite sign the intersects the zero-levelset
            // iso-surface.
            return true;
        }
    }

    // The cell is not in the narrow band
    if (maximumDistance) {
        *maximumDistance = -1.;
    }

    return false;
}

/*!
 * Check if the specified point lies within the narrow band.
 *
 * The value of the levelset is evaluated and compared with the specified narrow band size.
 *
 * \param[in] id is the cell id
 * \result Return true if the cell is in the narrow band, false otherwise.
 */
bool LevelSetObject::isInNarrowBand(const std::array<double,3> &point)const
{
    double value = evalValue(point, false);
    if (value <= m_narrowBandSize) {
        return true;
    }

    return false;
}

/*!
 * Evaluate levelset sign at the specified cell.
 * @param[in] id cell id
 * @result The sign of the levelset at the specified cell.
 */
short LevelSetObject::evalCellSign(long id) const {
    // Try fetching the value from the sign cache
    CellValueCache<short> *signCache = getFieldCellCache<short>(LevelSetField::SIGN);
    if (signCache) {
        CellValueCache<short>::Entry signCacheEntry = signCache->findEntry(id);
        if (signCacheEntry.isValid()) {
            return *signCacheEntry;
        }
    }

    // Try evaluating the sign using the cached value
    CellValueCache<double> *valueCache = getFieldCellCache<double>(LevelSetField::VALUE);
    if (valueCache) {
        CellValueCache<double>::Entry valueCacheEntry = valueCache->findEntry(id);
        if (valueCacheEntry.isValid()) {
            return static_cast<short>(sign(*valueCacheEntry));
        }
    }

    return _evalCellSign(id);
}

/*!
 * Evaluate levelset value at the specified cell.
 * @param[in] id cell id
 * @param[in] signedLevelSet controls if signed levelset function will be used
 * @result The value of the levelset at the specified cell.
 */
double LevelSetObject::evalCellValue(long id, bool signedLevelSet) const {
    // Early return if not cache is registered for the field
    CellValueCache<double> *cache = getFieldCellCache<double>(LevelSetField::VALUE);
    if (!cache) {
        return _evalCellValue(id, signedLevelSet);
    }

    // Try fetching the value from the cache
    CellValueCache<double>::Entry cacheEntry = cache->findEntry(id);
    if (cacheEntry.isValid()) {
        if (signedLevelSet) {
            return *cacheEntry;
        } else {
            return *cacheEntry / evalCellSign(id);
        }
    }

    // Evaluate and store the value in the cache
    //
    // The value stored in the cache is always signed.
    double value = _evalCellValue(id, true);
    cache->insertEntry(id, value);

    // Evaluate the value with the correct signdness
    if (signedLevelSet) {
        return value;
    } else {
        return value / evalCellSign(id);
    }
}

/*!
 * Evaluate levelset gradient at the specified cell.
 * @param[in] id cell id
 * @param[in] signedLevelSet controls if signed levelset function will be used
 * @result The gradient of the levelset at the specified cell.
 */
std::array<double,3> LevelSetObject::evalCellGradient(long id, bool signedLevelSet) const {
    // Early return if not cache is registered for the field
    CellValueCache<std::array<double, 3>> *cache = getFieldCellCache<std::array<double, 3>>(LevelSetField::GRADIENT);
    if (!cache) {
        return _evalCellGradient(id, signedLevelSet);
    }

    // Try fetching the gradient from the cache
    CellValueCache<std::array<double, 3>>::Entry cacheEntry = cache->findEntry(id);
    if (cacheEntry.isValid()) {
        if (signedLevelSet) {
            return *cacheEntry;
        } else {
            return (1. / evalCellSign(id)) * (*cacheEntry);
        }
    }

    // Evaluate and store the gradient in the cache
    //
    // The gradient stored in the cache is always signed.
    std::array<double, 3> gradient = _evalCellGradient(id, true);
    cache->insertEntry(id, gradient);

    // Evaluate the gradient with the correct signdness
    if (signedLevelSet) {
        return gradient;
    } else {
        return (1. / evalCellSign(id)) * gradient;
    }
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
 * @param[in] id cell id
 * @result The sign of the levelset at the specified point.
 */
short LevelSetObject::evalSign(const std::array<double,3> &point) const {
    return _evalSign(point);
}

/*!
 * Evaluate levelset value at the specified point.
 * @param[in] id cell id
 * @param[in] signedLevelSet controls if signed levelset function will be used
 * @result The value of the levelset at the specified point.
 */
double LevelSetObject::evalValue(const std::array<double,3> &point, bool signedLevelSet) const {
    return _evalValue(point, signedLevelSet);
}

/*!
 * Evaluate levelset gradient at the specified point.
 * @param[in] id cell id
 * @param[in] signedLevelSet controls if signed levelset function will be used
 * @result The gradient of the levelset at the specified point.
 */
std::array<double,3> LevelSetObject::evalGradient(const std::array<double,3> &point, bool signedLevelSet) const {
    return _evalGradient(point, signedLevelSet);
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
 * Writes LevelSetObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::dump( std::ostream &stream ){
    // Identifier
    utils::binary::write(stream, m_id) ;

    // Write fields
    std::size_t nEnabledOutputFields = m_enabledOutputFields.size() ;
    utils::binary::write(stream, nEnabledOutputFields) ;
    for (const auto &fieldEntry : m_enabledOutputFields) {
        utils::binary::write(stream, fieldEntry.first) ;
        utils::binary::write(stream, fieldEntry.second) ;
    }

    // Additional information
    _dump(stream) ;
}

/*!
 * Writes LevelSetObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::_dump( std::ostream &stream ){

    // Dump cache
    for (std::unique_ptr<CellCache> &cache : *m_cellCacheCollection) {
        if (!cache) {
            continue;
        }

        cache->dump(stream);
    }
}

/*!
 * Reads LevelSetObject from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::restore( std::istream &stream ){
    // Identifier
    utils::binary::read(stream, m_id) ;

    // Write fields
    std::size_t nEnabledVTKOutputs ;
    utils::binary::read(stream, nEnabledVTKOutputs) ;
    for (std::size_t i = 0; i < nEnabledVTKOutputs; ++i) {
        LevelSetField field ;
        std::string fieldName;
        utils::binary::read(stream, field) ;
        utils::binary::read(stream, fieldName) ;
        m_enabledOutputFields.insert({field, fieldName}) ;
    }

    // Additional information
    _restore(stream) ;
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
 * @param[in] field is the field that that should be enabled/disabled
 * @param[in] enable true for enabling, false for disabling
 */
void LevelSetObject::enableVTKOutput( const LevelSetFieldset &fieldset, bool enable) {

    std::stringstream objectNameStream;
    objectNameStream << getId();

    enableVTKOutput(fieldset, objectNameStream.str(), enable);

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
 * @param[in] writeField is the write field that that should be enabled/disabled
 * @param[in] objectName is the name that will be associated with the object
 * @param[in] enable true for enabling, false for disabling
 */
void LevelSetObject::enableVTKOutput( LevelSetWriteField writeField, const std::string &objectName, bool enable) {

    LevelSetFieldset fieldset;
    if( writeField==LevelSetWriteField::ALL){
        fieldset = getSupportedFields();

    } else if ( writeField==LevelSetWriteField::DEFAULT){
        fieldset.insert(LevelSetField::VALUE);
        fieldset.insert(LevelSetField::GRADIENT);

    } else {
        LevelSetField field = static_cast<LevelSetField>(writeField);
        if (getSupportedFields().count(field) == 0) {
            log::warning() << "The specified field is not supported by the levelset object" << std::endl;
            return;
        }

        fieldset.insert(field);
    }

    enableVTKOutput( fieldset, objectName, enable);

}

/*!
 * Enables or disables the VTK output
 * The output will be enabled only if the object supports it.
 * @param[in] field is the field that that should be enabled/disabled
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
 * The output will be enabled only if the object supports it.
 * @param[in] field is the field that that should be enabled/disabled
 * @param[in] objectName is the name that will be associated with the object
 * @param[in] enable true for enabling, false for disabling
 */
void LevelSetObject::enableVTKOutput( LevelSetField field, const std::string &objectName, bool enable) {

    // Discard fields that are not supported
    if (getSupportedFields().count(field) == 0) {
        return;
    }

    // Get name of the field
    std::stringstream nameStream;
    nameStream << "levelset";
    switch(field){

        case LevelSetField::VALUE:
            nameStream << "Value_" << objectName;
            break;

        case LevelSetField::GRADIENT:
            nameStream << "Gradient_" << objectName;
            break;

        case LevelSetField::NORMAL:
            nameStream << "Normal_" << objectName;
            break;

        case LevelSetField::PART:
            nameStream << "PartId_" << objectName;
            break;

        default:
            throw std::runtime_error ("Unsupported value of field in LevelSetObject::addDataToVTK() ");
            break;
    }

    std::string name = nameStream.str();

    // Check if the state of the filed is already the requested one
    VTK &vtkWriter = m_kernel->getMesh()->getVTK() ;
    if (enable == vtkWriter.hasData(name)) {
        return;
    }

    // Process the field
    if (!enable) {
        vtkWriter.removeData( name);
        m_enabledOutputFields.erase(field) ;
    } else {
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

            case LevelSetField::NORMAL:
                vtkWriter.addData<double>( name, VTKFieldType::VECTOR, VTKLocation::CELL, this);
                break;

            case LevelSetField::PART:
                vtkWriter.addData<int>( name, VTKFieldType::SCALAR, VTKLocation::CELL, this);
                break;

            default:
                throw std::runtime_error ("Unsupported value of field in LevelSetObject::addDataToVTK() ");
                break;

        }
        m_enabledOutputFields.insert({field, name}) ;
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
            flushField(field, stream, format);
        }
    }

}

/*!
 * Write the specified field to the given stream.
 *
 * @param[in] field is the field that will be written
 * @param[in] stream output stream
 * @param[in] format is the format which must be used. Supported options
 * are "ascii" or "appended". For "appended" type an unformatted binary
 * stream must be used
 */
void LevelSetObject::flushField(LevelSetField field, std::fstream &stream, VTKFormat format) const {

    switch(field) {

    case LevelSetField::VALUE:
    {
        void (*writeFunctionPtr)(std::fstream &, const double &) = nullptr;

        if(format==VTKFormat::APPENDED){
            writeFunctionPtr = genericIO::flushBINARY<double>;
        } else if(format==VTKFormat::ASCII){
            writeFunctionPtr = genericIO::flushASCII<double>;
        } else {
            BITPIT_UNREACHABLE("Non-existent VTK format.");
        }

        CellValueCache<double> *cache = getFieldCellCache<double>(LevelSetField::VALUE);
        for( const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange() ){
            if (cache) {
                long cellId = cell.getId();
                CellValueCache<double>::Entry cacheEntry = cache->findEntry(cellId);
                if (cacheEntry.isValid()) {
                    (*writeFunctionPtr)(stream, *cacheEntry);
                } else {
                    (*writeFunctionPtr)(stream, levelSetDefaults::VALUE);
                }
            } else {
                (*writeFunctionPtr)(stream, levelSetDefaults::VALUE);
            }
        }

        break;
    }

    case LevelSetField::SIGN:
    {
        void (*writeFunctionPtr)(std::fstream &, const short &) = nullptr;

        if(format==VTKFormat::APPENDED){
            writeFunctionPtr = genericIO::flushBINARY<short>;
        } else if(format==VTKFormat::ASCII){
            writeFunctionPtr = genericIO::flushASCII<short>;
        } else {
            BITPIT_UNREACHABLE("Non-existent VTK format.");
        }

        CellValueCache<short> *cache = getFieldCellCache<short>(LevelSetField::SIGN);
        for( const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange() ){
            if (cache) {
                long cellId = cell.getId();
                CellValueCache<short>::Entry cacheEntry = cache->findEntry(cellId);
                if (cacheEntry.isValid()) {
                    (*writeFunctionPtr)(stream, *cacheEntry);
                } else {
                    (*writeFunctionPtr)(stream, levelSetDefaults::SIGN);
                }
            } else {
                (*writeFunctionPtr)(stream, levelSetDefaults::SIGN);
            }
        }

        break;
    }

    case LevelSetField::GRADIENT:
    {
        void (*writeFunctionPtr)(std::fstream &, const std::array<double,3> &) = nullptr;

        if(format==VTKFormat::APPENDED){
            writeFunctionPtr = genericIO::flushBINARY<std::array<double,3>>;
        } else if(format==VTKFormat::ASCII){
            writeFunctionPtr = genericIO::flushASCII<std::array<double,3>>;
        } else {
            BITPIT_UNREACHABLE("Non-existent VTK format.");
        }

        CellValueCache<std::array<double, 3>> *cache = getFieldCellCache<std::array<double, 3>>(LevelSetField::GRADIENT);
        for( const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange() ){
            if (cache) {
                long cellId = cell.getId();
                CellValueCache<std::array<double, 3>>::Entry cacheEntry = cache->findEntry(cellId);
                if (cacheEntry.isValid()) {
                    (*writeFunctionPtr)(stream, *cacheEntry);
                } else {
                    (*writeFunctionPtr)(stream, levelSetDefaults::GRADIENT);
                }
            } else {
                (*writeFunctionPtr)(stream, levelSetDefaults::GRADIENT);
            }
        }

        break;
    }

    default:
    {
        throw std::runtime_error("Unable to write the field.");
    }

    }
}

/*!
 * Reads LevelSetObject from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::_restore( std::istream &stream ){

    // Restore cache
    for (std::unique_ptr<CellCache> &cache : *m_cellCacheCollection) {
        if (!cache) {
            continue;
        }

        cache->restore(stream);
    }
}

#if BITPIT_ENABLE_MPI
/*!
 * Exchange of data structures of kernel and objects on ghost cells.
 */
void LevelSetObject::exchangeGhosts(){

    if (!m_kernel->getMesh()->isPartitioned()) {
        return;
    }

    std::unique_ptr<DataCommunicator> dataCommunicator = m_kernel->createDataCommunicator();
    startExchange(m_kernel->getMesh()->getGhostCellExchangeSources(), dataCommunicator.get());
    completeExchange(m_kernel->getMesh()->getGhostCellExchangeTargets(), dataCommunicator.get());
}

/*!
 * Start exchange of data structures of kernel and objects.
 * @param[in] sendList list of elements to be send
 * @param[in,out] dataCommunicator is the data communicator that will be used
 * for the data exchange
 */
void LevelSetObject::startExchange( const std::unordered_map<int,std::vector<long>> &sendList,
                                    DataCommunicator *dataCommunicator){

    // Fill the send buffer with the  content from the LevelSetObject base
    // class and the specific derived class.
    for (const auto &entry : sendList) {
        // Create an empty send
        int rank = entry.first;
        dataCommunicator->setSend(rank, 0);

        // Write data in the buffer
        SendBuffer &buffer = dataCommunicator->getSendBuffer(rank);
        writeCommunicationBuffer(entry.second, buffer);
    }

    // Discover the receives
    dataCommunicator->discoverRecvs();

    // Start the sends
    dataCommunicator->startAllRecvs();

    // Start the sends
    dataCommunicator->startAllSends();

}

/*!
 * Complete exchange of data structures of kernel and objects.
 * @param[in] recvList list of elements to be received
 * @param[in,out] dataCommunicator is the data communicator that will be used
 * for the data exchange
 */
void LevelSetObject::completeExchange( const std::unordered_map<int,std::vector<long>> &recvList,
                                       DataCommunicator *dataCommunicator){

    // Check which data communications have arrived. For those which are
    // available start reading the databuffer into the data structure of
    // LevelSetObject and its derived classes.
    int nCompletedRecvs = 0;
    while (nCompletedRecvs < dataCommunicator->getRecvCount()) {
        int rank = dataCommunicator->waitAnyRecv();

        RecvBuffer &dataBuffer = dataCommunicator->getRecvBuffer(rank);
        readCommunicationBuffer(recvList.at(rank), dataBuffer);

        ++nCompletedRecvs;
    }

    dataCommunicator->waitAllSends();
    dataCommunicator->finalize();
}

/*!
 * Flushing of data to communication buffers for partitioning
 * @param[in] sendList list of cells to be sent
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
void LevelSetObject::writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &dataBuffer ){
    _writeCommunicationBuffer( sendList, dataBuffer) ;
    dataBuffer.squeeze( ) ;
}

/*!
 * Flushing of data to communication buffers for partitioning
 * @param[in] sendList list of cells to be sent
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
void LevelSetObject::_writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &dataBuffer ){

    // Write cache data
    for (std::unique_ptr<CellCache> &cache : *m_cellCacheCollection) {
        if (!cache) {
            continue;
        }

        cache->writeBuffer(sendList, dataBuffer);
    }
}

/*!
 * Processing of communication buffer into data structure
 * @param[in] recvList list of cells to be received
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetObject::readCommunicationBuffer( const std::vector<long> &recvList, RecvBuffer &dataBuffer ){
    _readCommunicationBuffer(recvList, dataBuffer) ;
}

/*!
 * Processing of communication buffer into data structure
 * @param[in] recvList list of cells to be received
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetObject::_readCommunicationBuffer( const std::vector<long> &recvList, RecvBuffer &dataBuffer ){

    // Read cache data
    for (std::unique_ptr<CellCache> &cache : *m_cellCacheCollection) {
        if (!cache) {
            continue;
        }

        cache->readBuffer(recvList, dataBuffer);
    }
}
#endif 

/*!
 * Set the cache to be used for the specified fieldset.
 *
 * If an existing cache is already defined for one of the the specified fields, it will be
 * destroyed and recreated from scratch.
 *
 * \param fieldset are the fields for which cache mode will be set
 * \param cacheMode is the type of cache that will be used for caching field information
 */
void LevelSetObject::setFieldCache(const LevelSetFieldset &fieldset, LevelSetCacheMode cacheMode)
{
    for (LevelSetField field : fieldset) {
        setFieldCache(field, cacheMode);
    }
}

/*!
 * Set the cache to be used for the specified field.
 *
 * If an existing cache is already defined for the specified field, it will be destroyed and
 * recreated from scratch.
 *
 * \param field is the field for which cache mode will be set
 * \param cacheMode is the type of cache that will be used for caching field information
 */
void LevelSetObject::setFieldCache(LevelSetField field, LevelSetCacheMode cacheMode)
{
    switch(field) {

    case LevelSetField::VALUE:
        unregisterFieldCellCache(field);
        registerFieldCellCache<double>(field, cacheMode);
        break;

    case LevelSetField::SIGN:
        unregisterFieldCellCache(field);
        registerFieldCellCache<short>(field, cacheMode);
        break;

    case LevelSetField::GRADIENT:
        unregisterFieldCellCache(field);
        registerFieldCellCache<std::array<double, 3>>(field, cacheMode);
        break;

    default:
        throw std::runtime_error ("The requested field is not supported!");

    }
}

/*!
 * Fill object's cache.
 */
void LevelSetObject::fillCache()
{
    // Initialize narrow band cache
    bool isNarrowBandCacheNeeded = false;
    for (std::size_t fieldIndex = 0; fieldIndex < static_cast<std::size_t>(LevelSetField::COUNT); ++fieldIndex) {
        LevelSetField field = static_cast<LevelSetField>(fieldIndex);
        if (getFieldCellCacheMode(field) == LevelSetCacheMode::NARROW_BAND) {
            isNarrowBandCacheNeeded = true;
            break;
        }
    }

    if (m_cellNarrowBandCacheId == CellCacheCollection::NULL_ID) {
        if (isNarrowBandCacheNeeded) {
            m_cellNarrowBandCacheId = registerCellCache<bool>(LevelSetCacheMode::ON_DEMAND);
        }
    } else {
        if (!isNarrowBandCacheNeeded) {
            unregisterCellCache(m_cellNarrowBandCacheId);
            m_cellNarrowBandCacheId = CellCacheCollection::NULL_ID;
        }
    }

    // Fill the caches
    //
    // Caches in narrow band mode should be filled before cahces in full mode, because the
    // latter may need information evaluated by the former.
    fillNarrowBandCellCaches();
    fillFullCellCaches();
}

/*!
 * Updates object's cache after a mesh update.
 * @param[in] adaptionData are the information about the adaption
 */
void LevelSetObject::updateCache( const std::vector<adaption::Info> &adaptionData ){

    // Remove stale objects from the caches
    for (const adaption::Info &adaptionInfo : adaptionData) {
        if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
            continue;
        }

        for (std::unique_ptr<CellCache> &cache : *m_cellCacheCollection) {
            if (!cache) {
                continue;
            }
            cache->erase(adaptionInfo.previous);
        }
    }

    // Fill the caches
    //
    // Caches in narrow band mode should be filled before cahces in full mode, because the
    // latter may need information evaluated by the former.
    fillNarrowBandCellCaches(adaptionData);
    fillFullCellCaches(adaptionData);

    // Reclaim unused cache space
    for (std::unique_ptr<CellCache> &cache : *m_cellCacheCollection) {
        if (m_kernel->getExpectedFillIn() == LevelSetFillIn::SPARSE) {
            cache->shrink_to_fit();
        }
    }
}

/*!
 * Clear object's cache.
 *
 * \param release if set to true the memory hold by the cache will be released, otherwise
 * the cache will be cleared but its memory may not be released
 */
void LevelSetObject::clearCache( bool release ){

    // Clear cache
    for (std::unique_ptr<CellCache> &cache : *m_cellCacheCollection) {
        if (!cache) {
            continue;
        }

        if (release) {
            cache.reset();
        } else {
            cache->clear();
        }
    }
}

/*!
 * Fill cell caches in "full" mode.
 */
void LevelSetObject::fillFullCellCaches()
{
    LevelSetFieldset fullCacheFieldset = getCachedFields(LevelSetCacheMode::FULL);
    if (fullCacheFieldset.empty()) {
        return;
    }

    // Fill caches
    //
    // Sign cache should be handled separately.
    std::vector<LevelSetField> fieldProcessList;
    for (LevelSetField field : fullCacheFieldset) {
        if (field == LevelSetField::SIGN) {
            continue;
        }

        fieldProcessList.emplace_back(field);
    }

    for (LevelSetField field : fieldProcessList) {
        CellCache *cache = getFieldCellCache(field);
        cache->reserve(m_kernel->getMesh()->getCells().size());
    }

    for (const Cell &cell : m_kernel->getMesh()->getCells()) {
        for (LevelSetField field : fieldProcessList) {
            fillFieldCellCache(cell.getId(), field);
        }
    }

    // Fill sign cache
    if (fullCacheFieldset.count(LevelSetField::SIGN) != 0) {
        CellCache *cache = getFieldCellCache(LevelSetField::SIGN);
        cache->reserve(m_kernel->getMesh()->getCells().size());

        fillFullCellSignCache();
    }
}

/*!
 * Fill cell caches in "full" mode after a mesh update.
 *
 * \param adaptionData are the information about the adaption
 */
void LevelSetObject::fillFullCellCaches(const std::vector<adaption::Info> &adaptionData)
{
    LevelSetFieldset fullCacheFieldset = getCachedFields(LevelSetCacheMode::FULL);
    if (fullCacheFieldset.empty()) {
        return;
    }

    // Fill caches
    //
    // Sign cache should be handled separately.
    std::vector<LevelSetField> fieldProcessList;
    for (LevelSetField field : fullCacheFieldset) {
        if (field == LevelSetField::SIGN) {
            continue;
        }

        fieldProcessList.emplace_back(field);
    }

    for (LevelSetField field : fieldProcessList) {
        CellCache *cache = getFieldCellCache(field);
        cache->reserve(m_kernel->getMesh()->getCells().size());
    }

    for (const adaption::Info &adaptionInfo : adaptionData) {
        if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
            continue;
        }

        for (long cellId : adaptionInfo.current) {
            for (LevelSetField field : fieldProcessList) {
                fillFieldCellCache(cellId, field);
            }
        }
    }

    // Fill sign cache
    if (fullCacheFieldset.count(LevelSetField::SIGN) != 0) {
        CellCache *cache = getFieldCellCache(LevelSetField::SIGN);
        cache->reserve(m_kernel->getMesh()->getCells().size());

        fillFullCellSignCache(adaptionData);
    }
}

/*!
 * Fill the sign cache on the whole domain.
 *
 * If value cache is in narrow band mode, sign can be propagated using narrow band information
 * as seeds. Otherwise, sign should be evaluated from scratch.
 */
void LevelSetObject::fillFullCellSignCache()
{
    if (getFieldCellCacheMode(LevelSetField::SIGN) != LevelSetCacheMode::FULL) {
        return;
    }

    const VolumeKernel *mesh = m_kernel->getMesh();
    if (getFieldCellCacheMode(LevelSetField::VALUE) == LevelSetCacheMode::NARROW_BAND) {
        CellValueCache<short> *signCache = getFieldCellCache<short>(LevelSetField::SIGN);

        VolumeKernel::CellConstIterator cellBegin = mesh->cellConstBegin();
        VolumeKernel::CellConstIterator cellEnd   = mesh->cellConstEnd();
        std::size_t nUnaccountedCells = mesh->getInternalCellCount();

        // Initialize sign using narrow band information
        for (VolumeKernel::CellConstIterator cellItr = cellBegin; cellItr != cellEnd; ++cellItr) {
            long cellId = cellItr.getId();
            if (isCellInNarrowBand(cellId)) {
                signCache->insertEntry(cellId, evalCellSign(cellId));
                --nUnaccountedCells;
            }
        }

        // Propagate sign
        //
        // Sign propagation is done in two stages:
        //  - in the first stage, the seeds for the propagation are the cells inside the
        //    narrow band;
        //  - in the second stage, the seeds for the propagation are the cells whose sign
        //    is still unknown.
        bool cellAdjacenciesAvailable = (mesh->getAdjacenciesBuildStrategy() != VolumeKernel::ADJACENCIES_NONE);

        std::vector<long> processList;
        std::vector<long> cellNeighStorage;
        for (int stage = 0; stage < 2; ++stage) {
            for (VolumeKernel::CellConstIterator seedItr = cellBegin; seedItr != cellEnd; ++seedItr) {
                // Select seeds
                //
                // Seed selection criterium depend on the current stage:
                //  - in the first stage, the seeds for the propagation are the cells inside the
                //    narrow band;
                //  - in the second stage, the seeds for the propagation are the cells whose sign
                //    is still unknown.
                long seedId = seedItr.getId();
                CellValueCache<short>::Entry seedSignCacheEntry = signCache->findEntry(seedId);
                if (stage == 0) {
                    if (!seedSignCacheEntry.isValid()) {
                        continue;
                    }
                } else {
                    if (seedSignCacheEntry.isValid()) {
                        continue;
                    }

                    seedSignCacheEntry = signCache->insertEntry(seedId, evalCellSign(seedId));
                    --nUnaccountedCells;
                    if (nUnaccountedCells == 0) {
                        break;
                    }
                }

                short seedSign = *seedSignCacheEntry;

                // Propagate sign seed to neighbour cells
                processList.assign(1, seedId);
                while (!processList.empty()) {
                    long cellId = processList.back();
                    processList.resize(processList.size() - 1);

                    // Set the sign of the cell
                    if (cellId != seedId) {
                        signCache->insertEntry(cellId, seedSign);
                        --nUnaccountedCells;
                        if (nUnaccountedCells == 0) {
                            break;
                        }
                    }

                    // Add cell neighbours with no sign to the process list
                    const long *cellNeighs;
                    int nCellNeighs;
                    if (cellAdjacenciesAvailable) {
                        const Cell &cell = *(mesh->getCells().rawFind(seedItr.getRawIndex()));
                        cellNeighs = cell.getAdjacencies();
                        nCellNeighs = cell.getAdjacencyCount();
                    } else {
                        cellNeighStorage.clear();
                        mesh->findCellFaceNeighs(cellId, &cellNeighStorage);
                        cellNeighs = cellNeighStorage.data();
                        nCellNeighs = cellNeighStorage.size();
                    }

                    for(int n = 0; n < nCellNeighs; ++n){
                        long neighId = cellNeighs[n];
                        if (!signCache->findEntry(neighId).isValid()) {
                            processList.push_back(neighId);
                        }
                    }
                }

                // Check if propagation has already reached all the cells
                if (nUnaccountedCells == 0) {
                    break;
                }
            }

            // Check if propagation has already reached all the cells
            if (nUnaccountedCells == 0) {
                break;
            }
        }
    } else {
        // Evaluate sign from scratch
        for (const Cell &cell : mesh->getCells()) {
            fillFieldCellCache(cell.getId(), LevelSetField::SIGN);
        }
    }
}

/*!
 * Fill the sign cache on the whole domain.
 *
 * If value cache is in narrow band mode, sign can be propagated using narrow band information
 * as seeds. Otherwise, sign should be evaluated from scratch.
 *
 * \param[in] adaptionData are the information about the adaption
 */
void LevelSetObject::fillFullCellSignCache(const std::vector<adaption::Info> &adaptionData)
{
    if (getFieldCellCacheMode(LevelSetField::SIGN) != LevelSetCacheMode::FULL) {
        return;
    }

    if (getFieldCellCacheMode(LevelSetField::VALUE) == LevelSetCacheMode::NARROW_BAND) {
        // Propagate sign using narrow band information
        fillFullCellSignCache();
    } else {
        // Evaluate sign from scratch
        for (const adaption::Info &adaptionInfo : adaptionData) {
            if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
                continue;
            }

            for (long cellId : adaptionInfo.current) {
                fillFieldCellCache(cellId, LevelSetField::SIGN);
            }
        }
    }
}

/*!
 * Fill cell caches in "narrow band" mode.
 *
 * If the size of the narrow band has been set, the method will fill the caches on the cells
 * that intersect the surface, on all their first neighbours and on the cells with a distance
 * from the surface less than the defined narrow band size.
 *
 * In case the size of the narrow band has not been set, the method will fill the caches on
 * the cells that intersect the surface and on all their first neighbours.
 */
void LevelSetObject::fillNarrowBandCellCaches()
{
    // Get field to process
    LevelSetFieldset narrowBandCacheFieldset = getCachedFields(LevelSetCacheMode::NARROW_BAND);
    std::vector<LevelSetField> fieldProcessList(narrowBandCacheFieldset.begin(), narrowBandCacheFieldset.end());
    if (fieldProcessList.empty()) {
        return;
    }

    // Get narrow band cache
    CellValueCache<bool> *narrowBandCache = getCellCache<bool>(m_cellNarrowBandCacheId);
    narrowBandCache->clear();

    // Get mesh information
    const VolumeKernel &mesh = *(m_kernel->getMesh()) ;

    // Identify cells that are within the narrow band
    //
    // A cell is within the narrow band if its distance from the surface is smaller than
    // the narrow band side or if it intersects the surface. Cells that intersect the
    // surface need to be identified because they will be further processed for finding
    // which of their neighbour is within the narrow band.
    VolumeKernel::CellConstIterator cellBegin = mesh.cellConstBegin();
    VolumeKernel::CellConstIterator cellEnd   = mesh.cellConstEnd();

    std::unordered_set<std::size_t> intersectedRawCellIds;
    for (VolumeKernel::CellConstIterator cellItr = cellBegin; cellItr != cellEnd; ++cellItr) {
        long cellId = cellItr.getId();

        // Check if the cell is within the narrow band
        //
        // No neighbour check is performed, cells with neighbours that intersect the
        // zero-levelset iso-surface the surface will be processed later.
        //
        // The check should be performed ignoring the narrow band cache, because we are
        // in the process of building that cache.
        double maximumDistance;
        bool cellInsideNarrowBand = _isCellInNarrowBand(cellId, false, &maximumDistance);
        if (!cellInsideNarrowBand) {
            continue;
        }

        // Fill cell caches
        narrowBandCache->insertEntry(cellId, true);
        for (LevelSetField field : fieldProcessList) {
            fillFieldCellCache(cellId, field, maximumDistance);
        }

        // Update the list of cells that intersects the surface
        //
        //
        // When the narrow band size is not explicitly set, the cell will always
        // intersects the surface because only cells that intersect the surface
        // are considered, otherwise we need to explicitly check if the cell
        // intersects the surface.
        if (m_narrowBandSize < 0 || intersectSurface(cellId, LevelSetIntersectionMode::FAST_GUARANTEE_FALSE) == LevelSetIntersectionStatus::TRUE) {
            std::size_t cellRawId = cellItr.getRawIndex();
            intersectedRawCellIds.insert(cellRawId);
        }
    }

    // Process cells that intersect the zero-levelset iso-surface
    //
    // Neighbours of cells that intersect the zero-levelset iso-surface need to be added to the
    // narrow band if their sign differs from the sign of the cell that intersect the surface.
    for (std::size_t cellRawId : intersectedRawCellIds) {
        // Compute cell projection point
        VolumeKernel::CellConstIterator cellItr = mesh.getCells().rawFind(cellRawId);
        std::array<double,3> cellProjectionPoint = evalCellProjectionPoint(cellItr.getId());
        short cellSign = evalCellSign(cellItr.getId());

        // Process cell adjacencies
        const long *neighbours = cellItr->getAdjacencies() ;
        int nNeighbours = cellItr->getAdjacencyCount() ;
        for (int n = 0; n < nNeighbours; ++n) {
            // Skip neighbours that have already been processed
            long neighId = neighbours[n];
            if (isCellInNarrowBand(neighId)) {
                continue;
            }

            // Skip neighbours with the same sign
            long neighSign = evalCellSign(neighId);
            if (neighSign == cellSign) {
                continue;
            }

            // Fill cell caches
            //
            // We use the information about the intersected cell to provide a distance hint to
            // the function that will fill the cache: we know that the levelset value will be
            // smaller that the distance between the neighbour centroid and the projection point
            // of the intersected cell (we increase that value by an arbitrary 5% to avoid any
            // approximation error).
            std::array<double,3> neighCentroid = m_kernel->computeCellCentroid(neighId);
            double searchRadius = 1.05 * norm2(neighCentroid - cellProjectionPoint);

            narrowBandCache->insertEntry(neighId, true);
            for (LevelSetField field : fieldProcessList) {
                fillFieldCellCache(neighId, field, searchRadius);
            }
        }
    }

#if BITPIT_ENABLE_MPI
    // Exchange data among processes
    //
    // It's not possible to evaluate if the cells on the last layer of ghost cells is inside the
    // narrow band only using local informaiont (we lack information on their negihbours).
    if (mesh.isPartitioned()) {
    }
#endif
}

/*!
 * Fill cell caches in "full" mode after a mesh update.
 *
 * If the size of the narrow band has been set, the method will fill the caches on the cells
 * that intersect the surface, on all their first neighbours and on the cells with a distance
 * from the surface less than the defined narrow band size.
 *
 * In case the size of the narrow band has not been set, the method will fill the caches on
 * the cells that intersect the surface and on all their first neighbours.
 *
 * \param adaptionData are the information about the adaption
 */
void LevelSetObject::fillNarrowBandCellCaches(const std::vector<adaption::Info> &adaptionData)
{
    // Get field to process
    LevelSetFieldset narrowBandCacheFieldset = getCachedFields(LevelSetCacheMode::NARROW_BAND);
    std::vector<LevelSetField> fieldProcessList(narrowBandCacheFieldset.begin(), narrowBandCacheFieldset.end());
    if (fieldProcessList.empty()) {
        return;
    }

    // Get narrow band cache
    CellValueCache<bool> *narrowBandCache = getCellCache<bool>(m_cellNarrowBandCacheId);

    // Identify updated cells that are within the narrow band
    //
    // A cell is within the narrow band if its levelset value is smaller than the narrow band
    // size or if it intersects the zero-levelset iso-surface. New cells that are outside the
    // narrow band need to be identified because if one of their neighbours intersects the
    // zero-levelset iso-surface, they need to be added to the narrow band.
    std::vector<long> cellsOutsideNarrowband;
    for (const adaption::Info &adaptionInfo : adaptionData) {
        if (adaptionInfo.entity != adaption::Entity::ENTITY_CELL) {
            continue;
        }

        if (adaptionInfo.type == adaption::Type::TYPE_PARTITION_SEND){
            continue;
        } else if (adaptionInfo.type == adaption::Type::TYPE_PARTITION_RECV){
            continue;
        }

        for (long cellId : adaptionInfo.current) {
            // Check if the cell is within the narrow band
            //
            // No neighbour check is performed, cells with neighbours that intersect the
            // zero-levelset iso-surface the surface will be processed later.
            //
            // The check should be performed ignoring the narrow band cache, because we are
            // in the process of building that cache.
            double maximumDistance;
            bool cellInsideNarrowBand = _isCellInNarrowBand(cellId, false, &maximumDistance);
            if (!cellInsideNarrowBand) {
                cellsOutsideNarrowband.push_back(cellId);
                continue;
            }

            // Fill field caches
            narrowBandCache->insertEntry(cellId, true);
            for (LevelSetField field : fieldProcessList) {
                fillFieldCellCache(cellId, field, maximumDistance);
            }
        }
    }

    // Process cells with neighbours that intersect the zero-levelset iso-surface
    //
    // Cells with at least a negihbour that intersect the zero-levelset iso-surface need to be
    // added to the narrow band if their sign differs form the sign of the neighbour that
    // intersects the surface.
    const VolumeKernel &mesh = *(m_kernel->getMesh()) ;
    for (long cellId : cellsOutsideNarrowband){
        const Cell &cell = mesh.getCell(cellId);
        short cellSign = evalCellSign(cellId);

        // Identify cells with a neighbour that intersects the zero-levelset iso-surface
        // and has a different sign.
        const long *neighbours = cell.getAdjacencies() ;
        int nNeighbours = cell.getAdjacencyCount() ;

        long intersectedNeighId = Cell::NULL_ID;
        for (int n = 0; n < nNeighbours; ++n) {
            long neighId = neighbours[n];

            // Skip neighoburs outside the narrow band
            if (!isCellInNarrowBand(neighId)) {
                continue;
            }

            // Skip neighbours with the same sign of the cell
            long neighSign = evalCellSign(neighId);
            if (neighSign == cellSign) {
                continue;
            }

            // Skip neighbours that doens't intersect the surface
            if (intersectSurface(neighId, LevelSetIntersectionMode::FAST_GUARANTEE_FALSE) != LevelSetIntersectionStatus::TRUE){
                continue;
            }

            intersectedNeighId = neighId;
            break;
        }

        if (intersectedNeighId == Cell::NULL_ID) {
            continue;
        }

        // Fill cell caches
        //
        // We use the information about the intersected neighbour to provide a distance hint
        // to the function that will fill the cache: we know that the levelset value will be
        // smaller that the distance between the cell centroid and the projection point of
        // the intersected neighbour (we increase that value by an arbitrary 5% to avoid any
        // approximation error).
        std::array<double,3> cellCentroid = m_kernel->computeCellCentroid(cellId);
        std::array<double,3> neighProjectionPoint = evalCellProjectionPoint(intersectedNeighId);
        double searchRadius = 1.05 * norm2(cellCentroid - neighProjectionPoint);

        narrowBandCache->insertEntry(cellId, true);
        for (LevelSetField field : fieldProcessList) {
            fillFieldCellCache(cellId, field, searchRadius);
        }
    }
}

/*!
 * Get the fields cached using the specified cache mode.
 *
 * \result The fields cached using the specified cache mode.
 */
LevelSetFieldset LevelSetObject::getCachedFields(LevelSetCacheMode cacheMode) const
{
    LevelSetFieldset fieldset;
    for (std::size_t fieldIndex = 0; fieldIndex < static_cast<std::size_t>(LevelSetField::COUNT); ++fieldIndex) {
        LevelSetField field = static_cast<LevelSetField>(fieldIndex);
        if (getFieldCellCacheMode(field) == cacheMode) {
            fieldset.insert(field);
        }
    }

    return fieldset;
}

/*!
 * Get a pointer to the cell cache for the specified field.
 *
 * If no cache was registered for the specified field, a null pointer is returned.
 *
 * \param field is the field whose cache will be retrieved
 * \result A pointer to the cell cache for the specified field.
 */
typename LevelSetObject::CellCache * LevelSetObject::getFieldCellCache(LevelSetField field) const
{
    std::size_t fieldIndex = static_cast<std::size_t>(field);

    return getCellCache(fieldIndex);
}

/*!
 * Get the cache mode associated with the specified field.
 *
 * \param field is the field whose cache mode will be retrieved
 * \result The cache mode associated with the specified field.
 */
LevelSetCacheMode LevelSetObject::getFieldCellCacheMode(LevelSetField field) const
{
    std::size_t fieldIndex = static_cast<std::size_t>(field);

    return m_cellFieldCacheModes[fieldIndex];
}

/*!
 * Unregister the specified field cache.
 *
 * \param fieldset is the fieldset for which the caches will be added
 */
void LevelSetObject::unregisterFieldCellCache(LevelSetField field)
{
    std::size_t fieldIndex = static_cast<std::size_t>(field);
    unregisterCellCache(fieldIndex);
}

/*!
 * Fill the specified field cache of the given cell.
 *
 * \param id is the id of the cell whose cache will be filled
 * \param field is the field whose cache will be filled
 * \param searchRadius all the portions of the surface with a distance greater than the search
 * radius will not be considered when evaluating the levelset. Trying to fill the cache cache
 * of a cell whose distance is greater than the search radius results in undefined behaviour.
 * Reducing the search radius could speedup the evaluation of levelset information.
 */
void LevelSetObject::fillFieldCellCache(long id, LevelSetField field, double searchRadius)
{
    BITPIT_UNUSED(searchRadius);

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
        throw std::runtime_error ("Unsupported field!");

    }
}

/*!
 * Get a pointer to the specified cell cache.
 *
 * If specified cell cache was registered, a null pointer is returned.
 *
 * \param cacheId the id of the cell that will be unregistered
 * \result A pointer to the specified cell cache.
 */
typename LevelSetObject::CellCache * LevelSetObject::getCellCache(std::size_t cacheId) const
{
    return m_cellCacheCollection->at(m_cellFieldCacheIds[cacheId]);
}

/*!
 * Unregister the specified cache.
 *
 * \param cacheId the id of the cell that will be unregistered
 */
void LevelSetObject::unregisterCellCache(std::size_t cacheId)
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
