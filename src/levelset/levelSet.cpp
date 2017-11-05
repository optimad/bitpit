/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

# include "bitpit_operators.hpp"
# include "bitpit_patchkernel.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_volcartesian.hpp"
# include "bitpit_voloctree.hpp"

# include "levelSetCommon.hpp"
# include "levelSetKernel.hpp"
# include "levelSetCartesian.hpp"
# include "levelSetOctree.hpp"
# include "levelSetObject.hpp"
# include "levelSetMetaObject.hpp"
# include "levelSetCachedObject.hpp"
# include "levelSetBoolean.hpp"
# include "levelSetSegmentation.hpp"
# include "levelSetMask.hpp"

# include "levelSet.hpp"

namespace bitpit {

/*!
 * @class LevelSet
 * @ingroup levelset
 *
 * @brief  Level Set driver class
 *
 * LevelSet is the main user interface class for computing signed- or unsigned- distance functions on Cartesian or Octree meshes
 * with respect to geometrical objects. The user needs to define the computional kernel by calling setMesh() and the objects which define the zero 
 * levelset via addObject().
 *
 * LevelSet will calculate the exact distance with respect the objects within a narrow band.
 * Outside this narrow band an approximate value will be calculated.
 *
 * The user may set the size of the narrow band explicitly.
 * Alternatively LevelSet will guarantee a least on cell center with exact levelset values across the zero-levelset iso-surface.
 *
 * LevelSet will test if the underlying mesh can provide a MPI communicator.
 * In case LevelSet is parallelized according the underlying mesh partitioning.
*/

/*!
 * Default constructor
 */
LevelSet::LevelSet() {

    m_objects.clear() ;

    m_useNarrowBand = false ;

    m_signedDF    = true ;
    m_propagateS  = false;

}

/*!
 * Destructor of LevelSet
*/
LevelSet::~LevelSet(){
    clear() ;
}

/*!
 * Sets the grid on which the levelset function should be computed.
 * Only cartesian and octree patches are supported at this moment.
 * @param[in] mesh computational grid
 */
void LevelSet::setMesh( VolumeKernel* mesh ) {

    if( VolCartesian* cartesian = dynamic_cast<VolCartesian*> (mesh) ){
        setMesh(cartesian) ;

    } else if( VolOctree* octree = dynamic_cast<VolOctree*> (mesh) ){
        setMesh(octree) ;
    
    } else{
        throw std::runtime_error ("Mesh non supported in LevelSet::setMesh()");
    } 

    for( auto &obj : m_objects){
        obj.second->setKernel(m_kernel.get());
    }

}

/*!
 * Sets the grid on which the levelset function should be computed.
 * @param[in] cartesian cartesian patch
 */
void LevelSet::setMesh( VolCartesian* cartesian ) {
    LevelSetKernel *kernel = new LevelSetCartesian( *cartesian) ;
    m_kernel = unique_ptr<LevelSetKernel>(kernel);

# if BITPIT_ENABLE_MPI
    // Initialize the communicator
    if (m_kernel->getMesh()->isPartitioned()) {
        m_kernel->initializeCommunicator();
    }
# endif
}

/*!
 * Sets the grid on which the levelset function should be computed.
 * @param[in] octree octree patch
 */
void LevelSet::setMesh( VolOctree* octree ) {
    LevelSetKernel *kernel = new LevelSetOctree( *octree) ;
    m_kernel = unique_ptr<LevelSetKernel>(kernel);

# if BITPIT_ENABLE_MPI
    // Initialize the communicator
    if (m_kernel->getMesh()->isPartitioned()) {
        m_kernel->initializeCommunicator();
    }
# endif
}

/*!
 * Adds a LevelSetSegmentation object
 * @param[in] segmentation surface segmentation
 * @param[in] angle feature angle
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( std::unique_ptr<SurfUnstructured> &&segmentation, double angle, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_objects.size();
    }

    LevelSetSegmentation* lsSeg = new LevelSetSegmentation(id, std::move(segmentation), angle ) ;
    LevelSetObject *object = static_cast<LevelSetObject *>(lsSeg);

    return registerObject(std::unique_ptr<LevelSetObject>(object));
}

/*!
 * Adds a LevelSetSegmentation object
 * @param[in] segmentation surface segmentation
 * @param[in] angle feature angle
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( SurfUnstructured *segmentation, double angle, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_objects.size();
    }

    LevelSetSegmentation* lsSeg = new LevelSetSegmentation(id, segmentation, angle) ;
    LevelSetObject *object = static_cast<LevelSetObject *>(lsSeg);

    return registerObject(std::unique_ptr<LevelSetObject>(object));
}

/*!
 * Adds a LevelSetSegmentation object
 * @param[in] segmentation surface segmentation
 * @param[in] angle feature angle
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( std::unique_ptr<SurfaceKernel> &&segmentation, double angle, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_objects.size();
    }

    SurfUnstructured *surfUnstructured = dynamic_cast<SurfUnstructured *>(segmentation.get());
    if (!surfUnstructured) {
        throw std::runtime_error ("Segmentation type not supported");
    }

    segmentation.release();
    std::unique_ptr<SurfUnstructured> surfUnstructuredUPtr = std::unique_ptr<SurfUnstructured>(surfUnstructured) ;
    LevelSetSegmentation* lsSeg = new LevelSetSegmentation(id, std::move(surfUnstructuredUPtr), angle) ;

    LevelSetObject *object = static_cast<LevelSetObject *>(lsSeg);
    return registerObject(std::unique_ptr<LevelSetObject>(object));
}

/*!
 * Adds a LevelSetSegmentation object
 * @param[in] segmentation surface segmentation
 * @param[in] angle feature angle
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( SurfaceKernel *segmentation, double angle, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_objects.size();
    }

    SurfUnstructured *surfUnstructured = dynamic_cast<SurfUnstructured *>(segmentation);
    if (!surfUnstructured) {
        throw std::runtime_error ("Segmentation type not supported");
    }

    LevelSetSegmentation* lsSeg = new LevelSetSegmentation(id,surfUnstructured, angle) ;

    LevelSetObject *object = static_cast<LevelSetObject *>(lsSeg);
    return registerObject(std::unique_ptr<LevelSetObject>(object));
}

/*!
 * Adds a boolean operation between two objects
 * @param[in] operation boolean operation
 * @param[in] id1 id of first operand
 * @param[in] id2 id of second operand
 * @param[in] id id to be assigned to object. In case default value is passed the insertion order will be used as identifier
 * @return identifier of new object
 */
int LevelSet::addObject( const LevelSetBooleanOperation &operation, const int &id1, const int &id2, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_objects.size();
    }

    LevelSetObject *ptr1 = m_objects.at(id1).get() ;
    LevelSetObject *ptr2 = m_objects.at(id2).get() ;

    return registerObject( std::unique_ptr<LevelSetObject>( new LevelSetBoolean(id, operation, ptr1, ptr2 ) ));
}

/*!
 * Adds a boolean operation between that will be applied recursivly to a series of objects
 * @param[in] operation boolean operation
 * @param[in] ids vector with indices of operand objects
 * @param[in] id id to be assigned to object. In case default value is passed the insertion order will be used as identifier
 * @return identifier of new object
 */
int LevelSet::addObject( const LevelSetBooleanOperation &operation, const std::vector<int> &ids, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_objects.size();
    }

    std::vector<LevelSetObject*> ptr;
    for( const int &id : ids){
        ptr.push_back( m_objects.at(id).get() );
    }

    return registerObject( std::unique_ptr<LevelSetObject>( new LevelSetBoolean(id, operation, ptr) ));
}
/*!
 * Adds a LevelSetMask object composed of the external envelope of a list of mesh cells.
 * The function setMesh() must have been called prior.
 * @param[in] list list of indices of cells
 * @param[in] id id to be assigned to object. In case default value is passed the insertion order will be used as identifier
 * @return identifier of new object
 */
int LevelSet::addObject( const std::unordered_set<long> &list, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_objects.size();
    }

    assert(m_kernel && " levelset: setMesh must be called befor adding a LevelSetMask object ");

    return registerObject( std::unique_ptr<LevelSetObject>( new LevelSetMask(id, list, *m_kernel->getMesh()) ) );
}

/*!
 * Adds a LevelSetMask object composed of a list of interfaces
 * The function setMesh() must have been called prior.
 * @param[in] list list of indices of interfaces
 * @param[in] id id to be assigned to object. In case default value is passed the insertion order will be used as identifier
 * @return identifier of new object
 */
int LevelSet::addObject( const std::vector<long> &list, const long &refInterface, const bool &invert, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_objects.size();
    }

    assert(m_kernel && " levelset: setMesh must be called befor adding a LevelSetMask object ");

    return registerObject( std::move(std::unique_ptr<LevelSetObject>( new LevelSetMask(id, list, refInterface, invert, *m_kernel->getMesh())) )  );
};

/*!
 * Adds a generic LevelSetObject
 * @param[in] object generic object
 * @return the index associated to the object
 */
int LevelSet::addObject( std::unique_ptr<LevelSetObject> &&object ) {

    return registerObject(std::move(object));
};

/*!
 * Adds a generic LevelSetObject
 * @param[in] object generic object
 * @return the index associated to the object
 */
int LevelSet::registerObject( std::unique_ptr<LevelSetObject> &&object ) {

    if( m_kernel){
        object->setKernel(m_kernel.get());
    }

    int objectId = object->getId();

    m_objects[objectId] = std::move(object) ;

    addProcessingOrder(objectId);

    return objectId;
}

/*!
 * Remove all levelset objects
 */
void LevelSet::removeObjects() {
    m_objects.clear();
    m_order.clear();
}

/*!
 * Remove a levelset object
 * @param[in] id id of object to be removed
 * @return true if object has been found and removed
 */
bool LevelSet::removeObject(int id) {
    if( m_objects.count(id) != 0){
        m_objects.erase(id);
        bool found = removeProcessingOrder(id);
        BITPIT_UNUSED(found);
        assert(found);
        return true;
    } 

    return false;
}

/*!
 * Adds the object to the processing order.
 *
 * The insertion order determines the processing order 
 * but priority is given to primary objects. 
 *
 * This function must be called whan a new object is inserted into m_objectss.
 *
 * @param[in] objectId the id of the newly added object
 */
void LevelSet::addProcessingOrder( int objectId ) {

    bool primary = m_objects.at(objectId)->isPrimary() ;

    if(primary){
        std::vector<int>::iterator orderItr = m_order.begin() ;
        bool iterate( orderItr != m_order.end()) ;

        while(iterate){
            int id = *orderItr ;
            if( m_objects.at(id)->isPrimary() ){
                ++orderItr ;
                iterate = orderItr != m_order.end() ;
            } else {
                iterate = false;
            }
        }

        m_order.insert(orderItr,objectId) ;

    } else {
        m_order.push_back(objectId) ;

    }

}

/*!
 * Removes the object from the processing order.
 * This function must be called whan a object is removed fromm_objectss.
 * @param[in] objectId the id of the newly added object
 * @return true if object has been found and removed
 */
bool LevelSet::removeProcessingOrder(int objectId){

    std::vector<int>::iterator orderItr;

    for(orderItr=m_order.begin(); orderItr!=m_order.end(); ++orderItr){
        if(*orderItr==objectId){
            m_order.erase(orderItr);
            return true;
        }
    }
    
    return false ;
}

/*!
 * Get a constant reference to the specified object.
 * If the specified id does not exist an exception is thrown.
 * @param id is the object id
 * @return reference to levelset object
 */
const LevelSetObject & LevelSet::getObject( int id) const{
    return getObject<LevelSetObject>(id);
}

/*!
 * Get a constant pointer to the specified object.
 * @param id is the object id
 * @return pointer to levelset object
 */
const LevelSetObject * LevelSet::getObjectPtr( int id) const{
    return getObjectPtr<LevelSetObject>(id);
}

/*!
 * Get vector of pointers to all object.
 * @return vector of pointers to levelset objects
 */
std::vector<LevelSetObject const *>  LevelSet::getObjectPtrs( ) const{
    return getObjectPtrs<LevelSetObject>();
}

/*!
 * Get the number of levelset objects
 * @return number of objects
 */
int LevelSet::getObjectCount( ) const{
    return m_objects.size() ;
}

/*!
 * Get the ids of the bodies.
 * @return a list of the body ids
*/
std::vector<int> LevelSet::getObjectIds( ) const{
    std::vector<int> ids ;
    ids.reserve(m_objects.size()) ;
    for(const auto &entry : m_objects) {
        ids.push_back(entry.first) ;
    }

    return ids ;
}

/*!
 * Clear LevelSet entirely, deleteing kernel and objects
 */
void LevelSet::clear(){
    m_kernel.reset();
    removeObjects();
}

/*!
 * Set if the signed or unsigned levelset function should be computed.
 * @param[in] flag true/false for signed /unsigned Level-Set function .
 */
void LevelSet::setSign(bool flag){
    m_signedDF = flag;

}

/*!
 * Set if the levelset sign has to be propagated from the narrow band to the whole domain.
 * @param[in] flag True/false to active/disable the propagation .
 */
void LevelSet::setPropagateSign(bool flag){
    m_propagateS = flag;
}

/*!
 * Manually set the physical size of the narrow band.
 * @param[in] r Size of the narrow band.
 */
void LevelSet::setSizeNarrowBand(double r){
    m_useNarrowBand = true ;
    for( auto &object :m_objects){
        object.second->setSizeNarrowBand(r) ;
    }
}

/*!
 * Computes levelset on given mesh with respect to the objects.
 * This routines needs to be called at least once.
 * Each object should compute the levelset and associated
 * information on both internal and ghost cells.
 */
void LevelSet::compute(){

    assert(m_kernel && "LevelSet::setMesh() must be called prior to LevelSet::compute()");

    for( int objectId : m_order){
        auto &visitor = *(m_objects.at(objectId)) ;
        visitor.computeLSInNarrowBand(m_signedDF) ;
#if BITPIT_ENABLE_MPI
        visitor.exchangeGhosts();
#endif
        if(m_propagateS) visitor.propagateSign() ;
    }
}

/*!
 * Updates the levelset after mesh adaptation.
 * Each object should compute the levelset and associated
 * information on both internal and ghost cells.
 *
 * @param[in] mapper mapper conatining mesh modifications
 */
void LevelSet::update( const std::vector<adaption::Info> &mapper ){

    assert(m_kernel && "LevelSet::setMesh() must be called prior to LevelSet::update()");

    // Udate the cache of the kernel
    m_kernel->updateGeometryCache( mapper ) ;

    // Update ls in narrow band
    for( int objectId : m_order){
        auto &visitor = *(m_objects.at(objectId)) ;
        visitor.clearAfterMeshAdaption(mapper) ;
        visitor.updateLSInNarrowBand( mapper, m_signedDF ) ;
#if BITPIT_ENABLE_MPI
        visitor.exchangeGhosts();
#endif
        if(m_propagateS) visitor.propagateSign() ;
    }
}

#if BITPIT_ENABLE_MPI
/*!
 * Updates the levelset after mesh partitioning.
 * @param[in] mapper mapper conatining mesh modifications
 */
void LevelSet::partition( const std::vector<adaption::Info> &mapper ){

    assert(m_kernel && "LevelSet::setMesh() must be called prior to LevelSet::partition()");

    // Set the communicator
    m_kernel->initializeCommunicator();

    // Udate the cache of the kernel
    m_kernel->updateGeometryCache( mapper ) ;

    // Compile send and receive lists
    std::unordered_map<int,std::vector<long>> sendList, recvList ;
    for( const auto &event : mapper){
        if( event.entity != adaption::Entity::ENTITY_CELL){
            continue;
        }

        if( event.type == adaption::Type::TYPE_PARTITION_SEND){
            sendList.insert({{event.rank,event.previous}}) ;
        } else if( event.type == adaption::Type::TYPE_PARTITION_RECV){
            recvList.insert({{event.rank,event.current}}) ;
        }
    }

    // Communicate according to new partitioning
    for( int objectId : m_order){
        auto &visitor = *(m_objects.at(objectId)) ;

        visitor.communicate( sendList, recvList, &mapper ) ;
    }

    for( int objectId : m_order){
        auto &visitor = *(m_objects.at(objectId)) ;

        visitor.exchangeGhosts();
    }

}
#endif

/*! 
 * Writes LevelSetKernel to stream in binary format
 * @param[in] stream output stream
 */
void LevelSet::dump( std::ostream &stream ){

    utils::binary::write(stream, m_order);
    utils::binary::write(stream, m_useNarrowBand);
    utils::binary::write(stream, m_signedDF);
    utils::binary::write(stream, m_propagateS);

    for( const auto &object : m_objects ){
        object.second->dump( stream ) ;
    }
}

/*! 
 * Reads LevelSetKernel from stream in binary format
 * @param[in] stream output stream
 */
void LevelSet::restore( std::istream &stream ){

    utils::binary::read(stream, m_order);
    utils::binary::read(stream, m_useNarrowBand);
    utils::binary::read(stream, m_signedDF);
    utils::binary::read(stream, m_propagateS);

    for( const auto &object : m_objects ){
        object.second->restore( stream ) ;
    }
}

}

