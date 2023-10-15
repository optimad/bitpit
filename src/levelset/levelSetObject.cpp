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

# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"
# include "bitpit_patchkernel.hpp"

# include "levelSetKernel.hpp"
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
LevelSetObject::LevelSetObject(int id) : m_nReferences(0), m_kernel(nullptr), m_narrowBandSize(levelSetDefaults::NARROWBAND_SIZE) {
    setId(id);
}

/*!
 * Copy constructor
 * @param[in] other is another object whose content is copied in this object
 */
LevelSetObject::LevelSetObject(const LevelSetObject &other)
    : m_id(other.m_id),
      m_nReferences(other.m_nReferences),
      m_enabledOutputFields(other.m_enabledOutputFields),
      m_kernel(other.m_kernel),
      m_narrowBandSize(other.m_narrowBandSize)
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
    : m_id(other.m_id),
      m_nReferences(other.m_nReferences),
      m_enabledOutputFields(other.m_enabledOutputFields),
      m_kernel(other.m_kernel),
      m_narrowBandSize(other.m_narrowBandSize)
{
    for ( const auto &fieldEntry : other.m_enabledOutputFields ) {
        other.enableVTKOutput(fieldEntry.first, false);
    }

    for ( const auto &fieldEntry : m_enabledOutputFields ) {
        enableVTKOutput(fieldEntry.first, true);
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
            for ( const auto &fieldEntry : m_enabledOutputFields ) {
                enabledOutputFieldset.insert(fieldEntry.first);
            }

            enableVTKOutput(enabledOutputFieldset, false);
        } catch (const std::exception &exception) {
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
    m_kernel = kernel;

    for ( const auto &fieldEntry : m_enabledOutputFields ) {
        enableVTKOutput( fieldEntry.first, true ) ;
    }
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
 * Computes the projection point of the cell center, i.e. the closest
 * point to the cell center on the zero level set
 * @param[in] id cell id
 * @return the projection point
 */
std::array<double,3> LevelSetObject::computeProjectionPoint(long id) const{
    double value = getValue(id);
    if(utils::DoubleFloatingEqual()(value,levelSetDefaults::VALUE)){
        return levelSetDefaults::POINT;
    }

    return m_kernel->computeCellCentroid(id) -value *getGradient(id);
}

/*!
 * Projects a vertex on the zero levelset
 * @param[in] coords point coordinates
 * @return the projected point
 */
std::array<double,3> LevelSetObject::computeProjectionPoint(const std::array<double,3> &coords) const{

    LevelSetInfo info = computeLevelSetInfo(coords);
    return coords -info.value *info.gradient;
}

/*!
 * Projects a vertex on the zero levelset
 * @param[in] vertexId index of the vertex
 * @return the projected point
 */
std::array<double,3> LevelSetObject::computeVertexProjectionPoint(long vertexId) const{

    const std::array<double,3> &coords = m_kernel->getMesh()->getVertexCoords(vertexId);
    return computeProjectionPoint(coords);
}

/*!
 * Get LevelSetInfo of cell
 * @param[in] cellId cell idex
 * @return LevelSetInfo of cell
*/
LevelSetInfo LevelSetObject::getLevelSetInfo(long cellId) const {

    return LevelSetInfo(getValue(cellId), getGradient(cellId));

}

/*!
 * Get the levelset value of cell
 * @param[in] cellId cell id
 * @return levelset value in cell
 */
double LevelSetObject::getLS(long cellId) const {

    return getValue(cellId);

}

/*!
 * Get the sign of the levelset function
 * @param[in] id cell id
 * @return sign of levelset
 */
short LevelSetObject::getSign(long id)const{
    return evalValueSign(getValue(id));
}

/*!
 * Eval the sign of the specified levelset value
 * @param[in] value is the levelset value
 * @return sign of levelset
 */
short LevelSetObject::evalValueSign(double value)const{
    return static_cast<short>(sign(value));
}

/*!
 * Get the current size of the narrow band.
 * A size equal or less than zero means that the levelset will be evaluated
 * only on cells that intersect the surface.
 * @return size of the current narrow band
 */
double LevelSetObject::getSizeNarrowBand()const{
    return m_narrowBandSize;
}

/*!
 * Manually set the size of the narrow band.
 * Setting a size equal or less than zero, levelset will be evaluated only on
 * the cells that intersect the surface and on all their first neighbours.
 * After setting the size of the narrowband, the levelset is not automatically
 * updated. It's up to the caller to make sure the levelset will be properly
 * updated if the size of the narrowband changes.
 * @param[in] r size of the narrow band.
 */
void LevelSetObject::setSizeNarrowBand(double r){
    m_narrowBandSize = r;
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

    double absoluteDistance  = std::abs(getValue(id));
    double distanceTolerance = m_kernel->getDistanceTolerance();

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

            std::array<double,3> root = computeProjectionPoint(id);
            std::array<double,3> normal = getGradient(id);
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
 * Updates the object after an adaption.
 *
 * @param[in] adaptionData are the information about the adaption
 * @param[in] signedDistance controls if signed- or unsigned- distance function should be calculated
 */
void LevelSetObject::update( const std::vector<adaption::Info> &adaptionData, bool signedDistance ) {

#if BITPIT_ENABLE_MPI
    UpdateStrategy partitioningUpdateStrategy = getPartitioningUpdateStrategy();
#endif

    std::vector<long> pruneList ;
    std::vector<long> updateList ;
#if BITPIT_ENABLE_MPI
    std::unordered_map<int, std::vector<long>> exchangeSendList ;
    std::unordered_map<int, std::vector<long>> exchangeRecvList ;
#endif
    for( const adaption::Info &adaptionInfo : adaptionData){
        if( adaptionInfo.entity != adaption::Entity::ENTITY_CELL){
            continue;
        }

        switch (adaptionInfo.type) {

#if BITPIT_ENABLE_MPI
        case adaption::Type::TYPE_PARTITION_SEND:
            if (partitioningUpdateStrategy == UPDATE_STRATEGY_EXCHANGE) {
                exchangeSendList.insert({{adaptionInfo.rank,adaptionInfo.previous}}) ;
            }
            break;

        case adaption::Type::TYPE_PARTITION_RECV:
            if (partitioningUpdateStrategy == UPDATE_STRATEGY_EXCHANGE) {
                exchangeRecvList.insert({{adaptionInfo.rank,adaptionInfo.current}}) ;
            } else if (partitioningUpdateStrategy == UPDATE_STRATEGY_EVALUATE) {
                updateList.insert(updateList.end(), adaptionInfo.current.begin(), adaptionInfo.current.end()) ;
            }
            break;
#endif

        default:
            updateList.insert(updateList.end(), adaptionInfo.current.begin(), adaptionInfo.current.end()) ;
            break;

        }

        pruneList.insert(pruneList.end(), adaptionInfo.previous.begin(), adaptionInfo.previous.end()) ;
    }

#if BITPIT_ENABLE_MPI
    // Initialize data exchange
    bool exchangeData = (!exchangeSendList.empty() || !exchangeRecvList.empty());
    if (m_kernel->getMesh()->isPartitioned()) {
        MPI_Allreduce(MPI_IN_PLACE, &exchangeData, 1, MPI_C_BOOL, MPI_LOR, m_kernel->getCommunicator()) ;
    }

    std::unique_ptr<DataCommunicator> dataCommunicator;
    if (exchangeData) {
        dataCommunicator = m_kernel->createDataCommunicator() ;
    }

    // Start data exchange
    if (exchangeData) {
        startExchange( exchangeSendList, dataCommunicator.get() ) ;
    }
#endif

    // Prune narrow band data structures
    if (!pruneList.empty()) {
        pruneNarrowBand( pruneList ) ;
    }

#if BITPIT_ENABLE_MPI
    // Complete data exchange
    if (exchangeData) {
        completeExchange( exchangeRecvList, dataCommunicator.get() ) ;
    }
#endif

    // Update narrow band
    if (!updateList.empty()) {
        updateNarrowBand( updateList, signedDistance ) ;
    }

#if BITPIT_ENABLE_MPI
    // Update data on ghost cells
    exchangeGhosts() ;
#endif
}

#if BITPIT_ENABLE_MPI
/*!
 * Get the strategy that should be used to update the object after a partitioning.
 * @result The strategy that should be used to update the object after a partitioning.
 */
LevelSetObject::UpdateStrategy LevelSetObject::getPartitioningUpdateStrategy() const {
    return UPDATE_STRATEGY_EXCHANGE;
}
#endif

/*!
 * Calculates the value and gradient of the levelset function within the narrow band
 * @param[in] signd if signed distances should be calculted
 */
void LevelSetObject::computeNarrowBand(bool signd){
    BITPIT_UNUSED(signd);
}

/*!
 * Updates the narrow band levelset function of the specified cells.
 * @param[in] cellIds are the ids of the cells that will be updated
 * @param[in] signd if signed distances should be calculted
 */
void LevelSetObject::updateNarrowBand(const std::vector<long> &cellIds, bool signd){
    BITPIT_UNUSED(cellIds);
    BITPIT_UNUSED(signd);
}

/*! 
 * Clear narrow band information associated with the specified cells.
 * @param[in] cellIds are the ids of the cells for which narrow band information will be deleted
 */
void LevelSetObject::pruneNarrowBand(const std::vector<long> &cellIds){
    BITPIT_UNUSED(cellIds);
}

/*! 
 * Clears all levelset information
 */
void LevelSetObject::clear( ){
    _clear() ;
}

/*! 
 * Clears all levelset information stored in derived class
 */
void LevelSetObject::_clear( ){
}

/*!
 * Writes LevelSetObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::dump( std::ostream &stream ){
    // Identifier
    utils::binary::write(stream, m_id) ;

    // Narroband size
    utils::binary::write(stream, m_narrowBandSize);

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
    BITPIT_UNUSED(stream);
}

/*!
 * Reads LevelSetObject from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::restore( std::istream &stream ){
    // Identifier
    utils::binary::read(stream, m_id) ;

    // Narroband size
    utils::binary::read(stream, m_narrowBandSize);

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

        case LevelSetField::GRADIENT:
            vtkWriter.addData<double>( name, VTKFieldType::VECTOR, VTKLocation::CELL, this);
            break;

        default:
            throw std::runtime_error ("Unsupported value of field in LevelSetObject::addDataToVTK() ");
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

        case LevelSetField::GRADIENT:
            return "Gradient";

        default:
            throw std::runtime_error ("Unsupported value of field in LevelSetObject::addDataToVTK() ");
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
            flushVTKOutputData(field, stream, format);
        }
    }

}

/*!
 * Write the VTK data associated with the specified field to the given stream.
 *
 * @param[in] field is the field
 * @param[in] stream is the output stream
 * @param[in] format is the format which must be used. Supported options
 * are "ascii" or "appended". For "appended" type an unformatted binary
 * stream must be used
 */
void LevelSetObject::flushVTKOutputData(LevelSetField field, std::fstream &stream, VTKFormat format) const {

    switch(field) {

    case LevelSetField::VALUE:
    {
        for( const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange() ){
            long cellId = cell.getId();
            double value = getValue(cellId);
            flushValue(stream, format, value);
        }

        break;
    }

    case LevelSetField::GRADIENT:
    {
        for( const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange() ){
            long cellId = cell.getId();
            const std::array<double,3> &value = getGradient(cellId);
            flushValue(stream, format, value);
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
    BITPIT_UNUSED(stream);
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
    BITPIT_UNUSED(sendList) ;
    BITPIT_UNUSED(dataBuffer) ;
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
    BITPIT_UNUSED(recvList) ;
    BITPIT_UNUSED(dataBuffer) ;
}

#endif 

}
