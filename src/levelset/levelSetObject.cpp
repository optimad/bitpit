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
      m_enabledVTKOutputs(other.m_enabledVTKOutputs),
      m_narrowBandSize(other.m_narrowBandSize)
{
    setKernel(other.m_kernel);
}

/*!
 * Move constructor
 * @param[in] other is another object whose content is copied in this object
 */
LevelSetObject::LevelSetObject(LevelSetObject &&other)
    : m_id(other.m_id),
      m_nReferences(other.m_nReferences),
      m_enabledVTKOutputs(other.m_enabledVTKOutputs),
      m_narrowBandSize(other.m_narrowBandSize)
{
    for ( LevelSetWriteField field : other.m_enabledVTKOutputs ) {
        other.enableVTKOutput(field, false);
    }

    setKernel(other.m_kernel);
}

/*!
 * Destructor.
 */
LevelSetObject::~LevelSetObject() {
    // Disable all output for the object
    if (m_kernel) {
        enableVTKOutput(LevelSetWriteField::ALL, false);
    }
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

    for ( LevelSetWriteField field : m_enabledVTKOutputs ) {
        enableVTKOutput( field, true ) ;
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
 * Get the part id of projection point
 * @param[in] id cell id
 * @return part id 
 */
int LevelSetObject::getPart(long id) const {
    BITPIT_UNUSED(id) ;
    return levelSetDefaults::PART ;
}

/*!
 * Get the surface normal at the projection point. 
 * The base implementation will return the levelset gradient.
 * @param[in] id cell id
 * @return surface normal
 */
std::array<double,3> LevelSetObject::getNormal(long id) const {
    return getGradient(id);
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
 * value to tangent and bounding radius of a cell. If the value is smaler than the
 * tangent radius LevelSetIntersectionStatus::TRUE is returned, if it is larger than the
 * bounding radius LevelSetIntersectionStatus::FALSE is returned. If it is inbetwee
 * LevelSetIntersectionStatus::CLOSE is returned.
 *
 * If mode==LevelSetIntersectionMode::FAST_GUARANTEE_TRUE and the levelset value is 
 * smaller than the rangent radius LevelSetIntersectionStatus::TRUE is retuned,
 * otherwise LevelSetIntersectionStatus::FALSE.
 *
 * If mode==LevelSetIntersectionMode::FAST_GURANTEE_FALSE and the levelset value is 
 * larger than the bounding radius LevelSetIntersectionStatus::FALSE is retuned,
 * otherwise LevelSetIntersectionStatus::TRUE.
 *
 * If mode==LevelSetIntersectionMode::ACCURATE, if LevelSetIntersectionMode::FUZZY 
 * returns LevelSetIntersectionStatus::CLOSE, the intersection between the tangent
 * plane at the projection point and the cell is performed additionally. 
 * LevelSetIntersectionStatus::TRUE/::FALSE is returned accordingly. 
 * Errors of the method are related to the ratio of surface curvature over cell size.
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
 * Returns the characterstic size at the support
 * @param[in] id cell id
 * @return feature size
 */
double LevelSetObject::getSurfaceFeatureSize(long id) const{
    BITPIT_UNUSED(id);
    return (- levelSetDefaults::SIZE);
}

/*!
 * Returns the minimum surface feature size
 * @return feature size
 */
double LevelSetObject::getMinSurfaceFeatureSize() const{
    return (- levelSetDefaults::SIZE);
}

/*!
 * Returns the maximum surface feature size
 * @return feature size
 */
double LevelSetObject::getMaxSurfaceFeatureSize() const{
    return (- levelSetDefaults::SIZE);
}

/*!
 * Calculates the value and gradient of the levelset function within the narrow band
 * @param[in] signd if signed distances should be calculted
 */
void LevelSetObject::computeNarrowBand(bool signd){
    BITPIT_UNUSED(signd);
}

/*!
 * Updates the value and gradient of the levelset function within the narrow band
 * @param[in] adaptionData are the information about the adaption
 * @param[in] signd if signed distances should be calculted
 */
void LevelSetObject::updateNarrowBand(const std::vector<adaption::Info> &adaptionData, bool signd){
    BITPIT_UNUSED(adaptionData);
    BITPIT_UNUSED(signd);
}

/*! 
 * Deletes non-existing items and items outside the narrow band after grid adaption.
 * @param[in] adaptionData are the information about the adaption
 */
void LevelSetObject::clearAfterMeshAdaption( const std::vector<adaption::Info> &adaptionData ){
    _clearAfterMeshAdaption( adaptionData ) ;
}

/*!
 * Clears data structure after mesh modification
 * @param[in] adaptionData are the information about the adaption
 */
void LevelSetObject::_clearAfterMeshAdaption( const std::vector<adaption::Info> &adaptionData ){
    BITPIT_UNUSED(adaptionData) ;
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

    // Enabled VTK outputs
    std::size_t nEnabledVTKOutputs = m_enabledVTKOutputs.size() ;
    utils::binary::write(stream, nEnabledVTKOutputs) ;
    for (LevelSetWriteField field : m_enabledVTKOutputs) {
        utils::binary::write(stream, field) ;
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

    // Enabled VTK outputs
    std::size_t nEnabledVTKOutputs ;
    utils::binary::read(stream, nEnabledVTKOutputs) ;
    for (std::size_t i = 0; i < nEnabledVTKOutputs; ++i) {
        LevelSetWriteField field ;
        utils::binary::read(stream, field) ;
        m_enabledVTKOutputs.insert(field) ;
    }

    // Additional information
    _restore(stream) ;
}

/*!
 * Enables or disables the VTK output
 * @param[in] fieldset describes the field(s) that should be written
 * @param[in] enable true for enabling, false for diabling
 */
void LevelSetObject::enableVTKOutput( LevelSetWriteField fieldset, bool enable) {

    std::stringstream objectNameStream;
    objectNameStream << getId();

    enableVTKOutput(fieldset, objectNameStream.str(), enable);

}

/*!
 * Enables or disables the VTK output
 * @param[in] fieldset describes the field(s) that should be written
 * @param[in] objectName is the name that will be associated with the object
 * @param[in] enable true for enabling, false for diabling
 */
void LevelSetObject::enableVTKOutput( LevelSetWriteField fieldset, const std::string &objectName, bool enable) {

    std::vector<LevelSetWriteField> fields;

    if( fieldset==LevelSetWriteField::ALL){
        fields.push_back(LevelSetWriteField::VALUE);
        fields.push_back(LevelSetWriteField::GRADIENT);
        fields.push_back(LevelSetWriteField::NORMAL);
        fields.push_back(LevelSetWriteField::PART);

    } else if ( fieldset==LevelSetWriteField::DEFAULT){
        fields.push_back(LevelSetWriteField::VALUE);
        fields.push_back(LevelSetWriteField::GRADIENT);

    } else {
        fields.push_back(fieldset);
    }

    VTK &vtkWriter = m_kernel->getMesh()->getVTK() ;
    for( const LevelSetWriteField &field : fields){
        // Update the list of enabled VTK outputs
        if (enable) {
            m_enabledVTKOutputs.insert(field) ;
        } else {
            m_enabledVTKOutputs.erase(field) ;
        }


        // Get name of the field
        std::stringstream nameStream;
        nameStream << "levelset";
        switch(field){

            case LevelSetWriteField::VALUE:
                nameStream << "Value_" << objectName;
                break;

            case LevelSetWriteField::GRADIENT:
                nameStream << "Gradient_" << objectName;
                break;

            case LevelSetWriteField::NORMAL:
                nameStream << "Normal_" << objectName;
                break;

            case LevelSetWriteField::PART:
                nameStream << "PartId_" << objectName;
                break;

            default:
                std::cout << " field " << (int) field << std::endl;
                throw std::runtime_error ("Unsupported value of field in LevelSetObject::addDataToVTK() ");
                break;
        }

        std::string name = nameStream.str();

        // Check if the state of the filed is already the requested one
        if (enable == vtkWriter.hasData(name)) {
            continue;
        }

        // Process the field
        if (!enable) {
            vtkWriter.removeData( name);
        } else {
            switch(field){

                case LevelSetWriteField::VALUE:
                    vtkWriter.addData<double>( name, VTKFieldType::SCALAR, VTKLocation::CELL, this);
                    break;

                case LevelSetWriteField::GRADIENT:
                    vtkWriter.addData<double>( name, VTKFieldType::VECTOR, VTKLocation::CELL, this);
                    break;

                case LevelSetWriteField::NORMAL:
                    vtkWriter.addData<double>( name, VTKFieldType::VECTOR, VTKLocation::CELL, this);
                    break;

                case LevelSetWriteField::PART:
                    vtkWriter.addData<int>( name, VTKFieldType::SCALAR, VTKLocation::CELL, this);
                    break;

                default:
                    throw std::runtime_error ("Unsupported value of field in LevelSetObject::addDataToVTK() ");
                    break;

            }
        }
    }

}

/*!
 * Reads LevelSetObject from stream in binary format
 * @param[in] stream output stream
 * @param[in] name is the name of the data to be written. Either user
 * data or patch data
 * @param[in] format is the format which must be used. Supported options
 * are "ascii" or "appended". For "appended" type an unformatted binary
 * stream must be used
 */
void LevelSetObject::flushData( std::fstream &stream, const std::string &name, VTKFormat format){


    if(utils::string::keywordInString(name,"levelsetValue")){

        void (*writeFunctionPtr)(std::fstream &, const double &) = nullptr;

        if(format==VTKFormat::APPENDED){
            writeFunctionPtr = genericIO::flushBINARY<double>;
        } else if(format==VTKFormat::ASCII){
            writeFunctionPtr = genericIO::flushASCII<double>;
        } else {
            BITPIT_UNREACHABLE("Non-existent VTK format.");
        }

        for( const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange() ){
            long cellId = cell.getId();
            double value = getValue(cellId);
            (*writeFunctionPtr)(stream,value);
        }

    } else if( utils::string::keywordInString(name,"levelsetGradient")){

        void (*writeFunctionPtr)(std::fstream &, const std::array<double,3> &) = nullptr;

        if(format==VTKFormat::APPENDED){
            writeFunctionPtr = genericIO::flushBINARY<std::array<double,3>>;
        } else if(format==VTKFormat::ASCII){
            writeFunctionPtr = genericIO::flushASCII<std::array<double,3>>;
        } else {
            BITPIT_UNREACHABLE("Non-existent VTK format.");
        }

        for( const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange() ){
            long cellId = cell.getId();
            const std::array<double,3> &value = getGradient(cellId);
            (*writeFunctionPtr)(stream,value);
        }

    } else if( utils::string::keywordInString(name,"levelsetNormal")){

        void (*writeFunctionPtr)(std::fstream &, const std::array<double,3> &) = nullptr;

        if(format==VTKFormat::APPENDED){
            writeFunctionPtr = genericIO::flushBINARY<std::array<double,3>>;
        } else if(format==VTKFormat::ASCII){
            writeFunctionPtr = genericIO::flushASCII<std::array<double,3>>;
        } else {
            BITPIT_UNREACHABLE("Non-existent VTK format.");
        }

        for( const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange() ){
            long cellId = cell.getId();
            const std::array<double,3> &value = getNormal(cellId);
            (*writeFunctionPtr)(stream,value);
        }

    } else if( utils::string::keywordInString(name,"levelsetPart")){

        void (*writeFunctionPtr)(std::fstream &, const int &) = nullptr;

        if(format==VTKFormat::APPENDED){
            writeFunctionPtr = genericIO::flushBINARY<int>;
        } else if(format==VTKFormat::ASCII){
            writeFunctionPtr = genericIO::flushASCII<int>;
        } else {
            BITPIT_UNREACHABLE("Non-existent VTK format.");
        }

        for( const Cell &cell : m_kernel->getMesh()->getVTKCellWriteRange() ){
            long cellId = cell.getId();
            int value = getPart(cellId);
            (*writeFunctionPtr)(stream,value);
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
