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

# include "bitpit_operators.hpp"
# include "bitpit_patchkernel.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_volcartesian.hpp"
# include "bitpit_voloctree.hpp"
# include "bitpit_volunstructured.hpp"

# include "levelSetCommon.hpp"
# include "levelSetKernel.hpp"
# include "levelSetCartesianKernel.hpp"
# include "levelSetOctreeKernel.hpp"
# include "levelSetObject.hpp"
# include "levelSetProxyObject.hpp"
# include "levelSetCachedObject.hpp"
# include "levelSetImmutableObject.hpp"
# include "levelSetBooleanObject.hpp"
# include "levelSetSegmentationObject.hpp"
# include "levelSetSignedObject.hpp"
# include "levelSetSignPropagator.hpp"
# include "levelSetMaskObject.hpp"
# include "levelSetObjectFactory.hpp"
# include "levelSetUnstructuredKernel.hpp"

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
 * LevelSet can use two types of storages: sparse or dense. When sparse storage
 * is used, objects will only allocate space for storing information of cells
 * that are inside the narrow band. Instead, when dense storage is used, objects
 * will still evaluate information only for cells inside the narrow band, but
 * they will allocate space for all cells of the mesh. If the narrow band covers
 * a large portion of the domain, it will be more efficient to allocate space
 * for all the cells of the mesh, rather than keep track of the cells inside
 * the narrow band and store space only for those cells.
 *
 * LevelSet will test if the underlying mesh can provide a MPI communicator.
 * In case LevelSet is parallelized according the underlying mesh partitioning.
*/

/*!
 * Default constructor
 *
 * \param storageType is the storage type that will be used for storing
 * levelset information.
 */
LevelSet::LevelSet(LevelSetStorageType storageType) : m_storageType(storageType) {

    m_objects.clear() ;

    m_narrowBandSize = levelSetDefaults::NARROWBAND_SIZE;

    m_signedDistance = true ;
    m_propagateSign  = false;

}

/*!
 * Get the storage type that will be used for storing levelset information.
 * @return storage type that will be used for storing levelset information
*/
LevelSetStorageType LevelSet::getStorageType() const{
    return m_storageType ;
}

/*!
 * Sets the mesh on which the levelset function should be computed.
 *
 * Only cartesian, octree and unstructured patches are supported. If the specified
 * mesh type is not among the supported types, an exception is thrown.
 *
 * @param[in] mesh computational mesh
 */
void LevelSet::setMesh( VolumeKernel* mesh ) {

    // Mesh can be set only once
    if (m_kernel) {
        throw std::runtime_error ("Mesh can be set just once.");
    }

    // Create the kernel
    if (VolCartesian *cartesian = dynamic_cast<VolCartesian *>(mesh)) {
        m_kernel = std::unique_ptr<LevelSetKernel>(new LevelSetCartesianKernel(*cartesian));
    } else if (VolOctree *octree = dynamic_cast<VolOctree *>(mesh)) {
        m_kernel = std::unique_ptr<LevelSetKernel>(new LevelSetOctreeKernel(*octree));
    } else if (VolUnstructured *unstructured = dynamic_cast<VolUnstructured*>(mesh)) {
        m_kernel = std::unique_ptr<LevelSetKernel>(new LevelSetUnstructuredKernel(*unstructured));
    } else {
        throw std::runtime_error ("Unable to create the levelset kernel. Mesh type non supported.");
    }

# if BITPIT_ENABLE_MPI
    // Initialize the communicator
    if (m_kernel->getMesh()->isPartitioned()) {
        m_kernel->initializeCommunicator();
    }
# endif

    for( auto &obj : m_objects){
        obj.second->setKernel(m_kernel.get());
    }

}

/*!
 * Adds a segmentation object
 * Objects can be added to the levelset only after setting the mesh.
 * @param[in] segmentation surface segmentation
 * @param[in] angle feature angle
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( std::unique_ptr<SurfUnstructured> &&segmentation, double angle, int id ) {

    std::unique_ptr<LevelSetObject> object = LevelSetObjectFactory::createSegmentationObject<LevelSetSegmentationNarrowBandCache, int &, std::unique_ptr<SurfUnstructured> &&, double &>(m_kernel.get(), m_storageType, id, std::move(segmentation), angle ) ;

    return registerObject(std::move(object));
}

/*!
 * Adds a segmentation object
 * Objects can be added to the levelset only after setting the mesh.
 * @param[in] segmentation surface segmentation
 * @param[in] angle feature angle
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( SurfUnstructured *segmentation, double angle, int id ) {

    std::unique_ptr<LevelSetObject> object = LevelSetObjectFactory::createSegmentationObject<LevelSetSegmentationNarrowBandCache, int &, SurfUnstructured * &, double &>( m_kernel.get(), m_storageType, id, segmentation, angle ) ;

    return registerObject(std::move(object));
}

/*!
 * Adds a segmentation object
 * Objects can be added to the levelset only after setting the mesh.
 * @param[in] segmentation surface segmentation
 * @param[in] angle feature angle
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( std::unique_ptr<SurfaceKernel> &&segmentation, double angle, int id ) {

    SurfUnstructured *surfUnstructured = dynamic_cast<SurfUnstructured *>(segmentation.get());
    if (!surfUnstructured) {
        throw std::runtime_error ("Segmentation type not supported");
    }

    segmentation.release();
    std::unique_ptr<SurfUnstructured> surfUnstructuredUPtr = std::unique_ptr<SurfUnstructured>(surfUnstructured) ;

    std::unique_ptr<LevelSetObject> object = LevelSetObjectFactory::createSegmentationObject<LevelSetSegmentationNarrowBandCache, int &, std::unique_ptr<SurfUnstructured> &&, double &>( m_kernel.get(), m_storageType, id, std::move(surfUnstructuredUPtr), angle ) ;

    return registerObject(std::move(object));
}

/*!
 * Adds a segmentation object
 * Objects can be added to the levelset only after setting the mesh.
 * @param[in] segmentation surface segmentation
 * @param[in] angle feature angle
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( SurfaceKernel *segmentation, double angle, int id ) {

    SurfUnstructured *surfUnstructured = dynamic_cast<SurfUnstructured *>(segmentation);
    if (!surfUnstructured) {
        throw std::runtime_error ("Segmentation type not supported");
    }

    std::unique_ptr<LevelSetObject> object = LevelSetObjectFactory::createSegmentationObject<LevelSetSegmentationNarrowBandCache, int &, SurfUnstructured * &, double &>( m_kernel.get(), m_storageType, id, surfUnstructured, angle ) ;

    return registerObject(std::move(object));
}

/*!
 * Adds a boolean operation between two objects
 * Objects can be added to the levelset only after setting the mesh.
 * @param[in] operation boolean operation
 * @param[in] id1 id of first operand
 * @param[in] id2 id of second operand
 * @param[in] id id to be assigned to object. In case default value is passed the insertion order will be used as identifier
 * @return identifier of new object
 */
int LevelSet::addObject( LevelSetBooleanOperation operation, int id1, int id2, int id ) {

    LevelSetObject *ptr1 = m_objects.at(id1).get() ;
    LevelSetObject *ptr2 = m_objects.at(id2).get() ;

    std::unique_ptr<LevelSetObject> object = LevelSetObjectFactory::createBooleanObject( m_kernel.get(), m_storageType, id, operation, ptr1, ptr2 ) ;

    return registerObject(std::move(object));
}

/*!
 * Adds a boolean operation between that will be applied recursivly to a series of objects
 * Objects can be added to the levelset only after setting the mesh.
 * @param[in] operation boolean operation
 * @param[in] ids vector with indices of operand objects
 * @param[in] id id to be assigned to object. In case default value is passed the insertion order will be used as identifier
 * @return identifier of new object
 */
int LevelSet::addObject( LevelSetBooleanOperation operation, const std::vector<int> &ids, int id ) {

    std::vector<const LevelSetObject*> ptr;
    for( int id : ids){
        ptr.push_back( m_objects.at(id).get() );
    }

    std::unique_ptr<LevelSetObject> object = LevelSetObjectFactory::createBooleanObject( m_kernel.get(), m_storageType, id, operation, ptr ) ;

    return registerObject(std::move(object));
}
/*!
 * Adds a LevelSetMask object composed of the external envelope of a list of mesh cells.
 * Objects can be added to the levelset only after setting the mesh.
 * @param[in] list list of indices of cells
 * @param[in] id id to be assigned to object. In case default value is passed the insertion order will be used as identifier
 * @return identifier of new object
 */
int LevelSet::addObject( const std::unordered_set<long> &list, int id ) {

    assert(m_kernel && " levelset: setMesh must be called befor adding a mask object ");

    std::unique_ptr<LevelSetObject> object = LevelSetObjectFactory::createMaskObject<LevelSetMaskNarrowBandCache, int &, const std::unordered_set<long> &, const VolumeKernel &>( m_kernel.get(), m_storageType, id, list, *m_kernel->getMesh() ) ;

    return registerObject(std::move(object));
}

/*!
 * Adds a LevelSetMask object composed of a list of interfaces
 * Objects can be added to the levelset only after setting the mesh.
 * @param[in] list list of indices of interfaces
 * @param[in] refInterface id of reference interface
 * @param[in] invert if orientation should be inverted with respect to the reference interface
 * @param[in] id id to be assigned to object. In case default value is passed the insertion order will be used as identifier
 * @return identifier of new object
 */
int LevelSet::addObject( const std::vector<long> &list, long refInterface, bool invert, int id ) {

    assert(m_kernel && " levelset: setMesh must be called befor adding a mask object ");

    std::unique_ptr<LevelSetObject> object = LevelSetObjectFactory::createMaskObject<LevelSetMaskNarrowBandCache, int &, const std::vector<long> &, long &, bool &, const VolumeKernel &>( m_kernel.get(), m_storageType, id, list, refInterface, invert, *m_kernel->getMesh() ) ;

    return registerObject(std::move(object));
};

/*!
 * Adds a generic LevelSetObject
 * Objects can be added to the levelset only after setting the mesh.
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

    int id = object->getId();
    if (id == levelSetDefaults::OBJECT) {
        id = m_objectIdentifierGenerator.generate();
        object->setId(id);
    } else {
        m_objectIdentifierGenerator.setAssigned(id);
    }

    if( m_kernel){
        object->setKernel(m_kernel.get());
    }

    if (m_narrowBandSize > 0.) {
        object->setSizeNarrowBand(m_narrowBandSize);
    }

    m_objects[id] = std::move(object) ;

    setObjectProcessingOrder(id);

    incrementObjectsReferenceCount(id);

    return id;
}

/*!
 * Remove all levelset objects
 */
void LevelSet::removeObjects() {
    m_objectsProcessingOrder.clear();
    m_objectIdentifierGenerator.reset();
    m_objects.clear();
}

/*!
 * Remove a levelset object
 * Objects that are sources for other objects cannot be deleted.
 * @param[in] id id of object to be removed
 * @return true if object has been found and removed
 */
bool LevelSet::removeObject(int id) {
    return removeObject(id, false);

}

/*!
 * Remove a levelset object
 * Objects that are sources for other objects cannot be deleted unless the
 * removal is forced.
 * @param[in] id id of object to be removed
 * @param[in] force if set to true the object will be deleted also if
 * it's a non-remobable object
 * @return true if object has been found and removed
 */
bool LevelSet::removeObject(int id, bool force) {
    if( !force && !isObjectRemovable(id) ) {
        return false;
    }

    decrementObjectsReferenceCount(id);
    unsetObjectProcessingOrder(id);
    m_objectIdentifierGenerator.trash(id);
    m_objects.erase(id);

    return true;

}

/*!
 * Check if an object can be deleted.
 * Non-existent objects or objects that are sources for other objects cannot
 * be deleted.
 * @param[in] id id of object
 * @return true if object can be deleted
 */
bool LevelSet::isObjectRemovable(int id) {
    auto objectItr = m_objects.find(id) ;
    if( objectItr == m_objects.end() ){
        return false;
    }

    const LevelSetObject *object = objectItr->second.get() ;
    if( object->getReferenceCount() > 0 ) {
        return false;
    }

    return true;

}

/*!
 * Convert the specified objct into an immutable object.
 *
 * The original object will be deleted and replaced with an immutable
 * object.
 *
 * @param[in] id the id of the object
 */
void LevelSet::makeObjectImmutable( int id ) {

    LevelSetObject *object = getObjectPtr( id ) ;

    // Get the list of referring objects
    //
    // Referring objects are objects for which this object is a source.
    std::vector<int> referringProxyObjectIds;
    if( object->getReferenceCount() > 0 ) {
        for( const auto &entry : m_objects ) {
            const LevelSetProxyObject *proxyObject = dynamic_cast<const LevelSetProxyObject*>(entry.second.get());
            if (!proxyObject) {
                continue;
            }

            std::vector<const LevelSetObject *> sourceObjects = proxyObject->getSourceObjects() ;
            std::size_t nSourceObjects = sourceObjects.size() ;
            for ( std::size_t i = 0; i < nSourceObjects; ++i ){
                const LevelSetObject *sourceObject = sourceObjects[i];
                if (sourceObject == object) {
                    referringProxyObjectIds.push_back(proxyObject->getId());
                }
            }
        }
    }

    // Create the immutable object
    std::unique_ptr<LevelSetObject> immutableObject = LevelSetObjectFactory::createImmutableObject<LevelSetImmutableNarrowBandCache>( m_kernel.get(), m_storageType, object ) ;
    assert( immutableObject );

    // Update referring objects
    for (int referringProxyObjectId : referringProxyObjectIds) {
        LevelSetProxyObject *referringProxyObject = static_cast<LevelSetProxyObject *>(getObjectPtr( referringProxyObjectId ));
        referringProxyObject->replaceSourceObject(object, immutableObject.get());
    }

    // Replace the object with the newly created immutable object
    removeObject(id, true);
    addObject(std::move(immutableObject));
}

/*!
 * Set the processing order of the specified object.
 *
 * The insertion order determines the processing order, however priority is
 * given to primary objects.
 *
 * Immutable objects doesn't need to be processed.
 *
 * This function must be called when a new object is added.
 *
 * @param[in] id the id of the object
 */
void LevelSet::setObjectProcessingOrder( int id ) {

    // Immutable objects are not processed.
    if (dynamic_cast<const LevelSetImmutableObjectBase*>(getObjectPtr(id))) {
        return;
    }

    // Define the processing order for the object
    //
    // The insertion order determines the processing order, however priority
    // is given to primary objects.
    std::vector<int>::iterator processingOrderItr;
    if(getObjectPtr(id)->isPrimary()){
        std::vector<int>::iterator processingOrderBegin = m_objectsProcessingOrder.begin();
        std::vector<int>::iterator processingOrderEnd   = m_objectsProcessingOrder.end();

        for (processingOrderItr = processingOrderBegin; processingOrderItr != processingOrderEnd; ++processingOrderItr) {
            int candidateId = *processingOrderItr ;
            const LevelSetObject *candidateObject = getObjectPtr(candidateId) ;
            if( !candidateObject->isPrimary() ){
                break;
            }
        }
    } else {
        processingOrderItr = m_objectsProcessingOrder.end();
    }

    m_objectsProcessingOrder.insert(processingOrderItr,id) ;

}

/*!
 * Unset the processing order of the specified object.
 * Immutable objects are not processed.
 * This function must be called whan a object is removed.
 * @param[in] id the id of the object
 */
void LevelSet::unsetObjectProcessingOrder(int id){

    // Immutable objects are not processed.
    if (dynamic_cast<const LevelSetImmutableObjectBase*>(getObjectPtr(id))) {
        return;
    }

    // Remove the object from the list of processed objects
    std::vector<int>::iterator processingOrderBegin = m_objectsProcessingOrder.begin();
    std::vector<int>::iterator processingOrderEnd   = m_objectsProcessingOrder.end();

    std::vector<int>::iterator processingOrderItr = std::find(processingOrderBegin, processingOrderEnd, id);
    assert(processingOrderItr != processingOrderEnd);
    m_objectsProcessingOrder.erase(processingOrderItr);

}

/*!
 * Increment reference counting for the sources of the specified object.
 * @param parentId is the id of the parent object
 */
void LevelSet::incrementObjectsReferenceCount(int parentId) {

    const LevelSetObject *parentObject = getObjectPtr(parentId);
    if( parentObject->isPrimary() ){
        return;
    }

    if( const LevelSetProxyObject *parentProxyObject = dynamic_cast<const LevelSetProxyObject*>(parentObject) ){
        for ( const LevelSetObject *sourceObject : parentProxyObject->getSourceObjects() ){
            int sourceObjectId = sourceObject->getId();
            getObjectPtr(sourceObjectId)->incrementReferenceCount();
        }
    }

}

/*!
 * Decrement reference counting for the sources of the specified object.
 * @param parentId is the id of the parent object
 */
void LevelSet::decrementObjectsReferenceCount(int parentId) {

    const LevelSetObject *parentObject = getObjectPtr(parentId);
    if( parentObject->isPrimary() ){
        return;
    }

    if( const LevelSetProxyObject *parentProxyObject = dynamic_cast<const LevelSetProxyObject*>(parentObject) ){
        for ( const LevelSetObject *sourceObject : parentProxyObject->getSourceObjects() ){
            int sourceObjectId = sourceObject->getId();
            getObjectPtr(sourceObjectId)->decrementReferenceCount();
        }
    }

}

/*!
 * Get a constant reference to the specified object.
 * If the specified id does not exist an exception is thrown.
 * @param id is the object id
 * @return reference to levelset object
 */
LevelSetObject & LevelSet::getObject( int id) const{
    return getObject<LevelSetObject>(id);
}

/*!
 * Get a constant pointer to the specified object.
 * @param id is the object id
 * @return pointer to levelset object
 */
LevelSetObject * LevelSet::getObjectPtr( int id) const{
    return getObjectPtr<LevelSetObject>(id);
}

/*!
 * Get vector of pointers to all object.
 * @return vector of pointers to levelset objects
 */
std::vector<LevelSetObject *>  LevelSet::getObjectPtrs( ) const{
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
    removeObjects();
    m_kernel.reset();
}

/*!
 * Set if the signed or unsigned levelset distance should be computed.
 * @param[in] flag true/false for signed /unsigned levelset distance.
 */
void LevelSet::setSign(bool flag){
    m_signedDistance = flag;

}

/*!
 * Set if the levelset sign has to be propagated from the narrow band to the whole domain.
 * @param[in] flag True/false to active/disable the propagation .
 */
void LevelSet::setPropagateSign(bool flag){
    m_propagateSign = flag;
}

/*!
 * Get the physical size of the narrow band.
 * A size equal or less than zero means that the levelset will be evaluated
 * only on cells that intersect the surface.
 * @return the physical size of the narrow band
 */
double LevelSet::getSizeNarrowBand() const{
    return m_narrowBandSize;
}

/*!
 * Manually set the physical size of the narrow band.
 * Setting a size equal or less than zero, levelset will be evaluated only on
 * the cells that intersect the surface and on all their first neighbours.
 * After setting the size of the narrowband, the levelset is not automatically
 * updated. It's up to the caller to make sure the levelset will be properly
 * updated if the size of the narrowband changes.
 * @param[in] r Size of the narrow band.
 */
void LevelSet::setSizeNarrowBand(double r){
    m_narrowBandSize = r;

    for( auto &object :m_objects){
        object.second->setSizeNarrowBand(m_narrowBandSize);
    }
}

/*!
 * Computes levelset on given mesh with respect to the objects.
 * Levelset and associated information will be computed on both internal and
 * ghost cells.
 */
void LevelSet::compute(){

    std::unordered_set<LevelSetObject *> objectProcessList = getObjectProcessList();

    compute(objectProcessList);

}

/*!
 * Computes levelset on given mesh with respect to the specified object.
 * All source objects needed for evaluating the levelset on the specified
 * object will be computed as well.
 * Levelset and associated information will be computed on both internal and
 * ghost cells.
 * @param[in] id identifier of object.
 */
void LevelSet::compute( int id ){

    std::unordered_set<LevelSetObject *> objectProcessList = getObjectProcessList(1, &id);

    compute(objectProcessList);

}

/*!
 * Computes levelset on given mesh with respect to the specified objects.
 * All source objects needed for evaluating the levelset on the specified
 * objects will be computed as well.
 * Levelset and associated information will be computed on both internal and
 * ghost cells.
 * @param[in] ids identifiers of objects.
 */
void LevelSet::compute( const std::vector<int> &ids ){

    std::unordered_set<LevelSetObject *> objectProcessList = getObjectProcessList(ids.size(), ids.data());

    compute(objectProcessList);

}

/*!
 * Computes levelset on given mesh with respect to the specified objects.
 * Only the specified objects will be computed, it's up to the caller to
 * provide an object process list that contains all the needed source
 * objects.
 * Levelset and associated information will be computed on both internal and
 * ghost cells.
 * @param[in] nObjects identifiers of objects.
 * @param[in] objectIds identifiers of objects.
 */
void LevelSet::compute( const std::unordered_set<LevelSetObject *> &objectProcessList ){

    assert(m_kernel && "LevelSet::setMesh() must be called prior to LevelSet::compute()");

    std::unique_ptr<LevelSetSignPropagator> signPropagator ;
    if (m_propagateSign) {
        signPropagator = m_kernel->createSignPropagator() ;
    }

    for( int id : m_objectsProcessingOrder){
        LevelSetObject *object = m_objects.at(id).get() ;
        if (objectProcessList.count(object) == 0) {
            continue;
        }

        LevelSetSignedObjectInterface *signPropagationObject;
        if (m_propagateSign) {
            signPropagationObject = dynamic_cast<LevelSetSignedObjectInterface *>(object);
        } else {
            signPropagationObject = nullptr;
        }

        // Clear the object
        object->clear();

        // Compute levelset inside the narrowband
        object->computeNarrowBand(m_signedDistance) ;
#if BITPIT_ENABLE_MPI
        object->exchangeGhosts();
#endif

        // Propagate sign
        if(signPropagationObject) {
            signPropagator->execute(signPropagationObject);
        }
    }
}

/*!
 * Updates the levelset after a mesh update.
 * Before being able to update the levelset, it has to be computed.
 * Levelset and associated information will be updated on both internal and
 * ghost cells.
 * @param[in] adaptionData are the information about the adaption
 */
void LevelSet::update( const std::vector<adaption::Info> &adaptionData ){

    std::unordered_set<LevelSetObject *> objectProcessList = getObjectProcessList();

    update(adaptionData, objectProcessList);

}

/*!
 * Computes levelset on given mesh with respect to the specified object.
 * Before being able to update the levelset, it has to be computed.
 * All source objects needed for evaluating the levelset on the specified
 * object will be updated as well.
 * Levelset and associated information will be updated on both internal and
 * ghost cells.
 * @param[in] adaptionData are the information about the adaption
 * @param[in] id identifier of object.
 */
void LevelSet::update( const std::vector<adaption::Info> &adaptionData, int id){

    std::unordered_set<LevelSetObject *> objectProcessList = getObjectProcessList(1, &id);

    update(adaptionData, objectProcessList);

}

/*!
 * Computes levelset on given mesh with respect to the specified objects.
 * Before being able to update the levelset, it has to be computed.
 * All source objects needed for evaluating the levelset on the specified
 * objects will be updated as well.
 * Levelset and associated information will be updated on both internal and
 * ghost cells.
 * @param[in] adaptionData are the information about the adaption
 * @param[in] ids identifiers of objects.
 */
void LevelSet::update( const std::vector<adaption::Info> &adaptionData, const std::vector<int> &ids ){

    std::unordered_set<LevelSetObject *> objectProcessList = getObjectProcessList(ids.size(), ids.data());

    update(adaptionData, objectProcessList);

}

/*!
 * Updates the levelset after a mesh update.
 * Before being able to update the levelset, it has to be computed.
 * Only the specified objects will be updated, it's up to the caller to
 * provide an object process list that contains all the needed source
 * objects.
 * Levelset and associated information will be updated on both internal and
 * ghost cells.
 * @param[in] adaptionData are the information about the adaption
 */
void LevelSet::update( const std::vector<adaption::Info> &adaptionData, const std::unordered_set<LevelSetObject *> &objectProcessList ){

    assert(m_kernel && "LevelSet::setMesh() must be called prior to LevelSet::partition()");

    VolumeKernel *mesh = m_kernel->getMesh();
#if BITPIT_ENABLE_MPI
    bool isMeshPartitioned = mesh->isPartitioned();
#endif

#if BITPIT_ENABLE_MPI
    // Set the communicator
    //
    // The mesh may have been partitioned after being set.
    if (isMeshPartitioned && !m_kernel->isCommunicatorSet()) {
        m_kernel->initializeCommunicator();
    }
#endif

    // Inspect adaption data to detect what needs to be updated
    std::unordered_map<int,std::vector<long>> partitioningSendList ;
    std::unordered_map<int,std::vector<long>> partitioningRecvList ;

    bool updateNarrowBand   = false;
    bool updatePartitioning = false;
    for( const adaption::Info &adaptionInfo : adaptionData){
        if( adaptionInfo.entity != adaption::Entity::ENTITY_CELL){
            continue;
        }

        if( adaptionInfo.type == adaption::Type::TYPE_PARTITION_SEND){
            partitioningSendList.insert({{adaptionInfo.rank,adaptionInfo.previous}}) ;
            updatePartitioning = true;
        } else if( adaptionInfo.type == adaption::Type::TYPE_PARTITION_RECV){
            partitioningRecvList.insert({{adaptionInfo.rank,adaptionInfo.current}}) ;
            updatePartitioning = true;
        } else if (adaptionInfo.type == adaption::Type::TYPE_DELETION) {
            updateNarrowBand = true;
        } else {
            if (!updateNarrowBand) {
                for (long cellId : adaptionInfo.current) {
                    const Cell &cell = mesh->getCell(cellId);
                    if (cell.isInterior()) {
                        updateNarrowBand = true;
                        break;
                    }
                }
            }
        }
    }

#if BITPIT_ENABLE_MPI
    if (isMeshPartitioned) {
        MPI_Allreduce(MPI_IN_PLACE, &updateNarrowBand, 1, MPI_C_BOOL, MPI_LOR, m_kernel->getCommunicator());
        MPI_Allreduce(MPI_IN_PLACE, &updatePartitioning, 1, MPI_C_BOOL, MPI_LOR, m_kernel->getCommunicator());
    }
#endif

    // Early return if no update is needed
    if (!updateNarrowBand && !updatePartitioning) {
        return;
    }

#if BITPIT_ENABLE_MPI
    // Get data communicator for updating partitioning
    std::unique_ptr<DataCommunicator> dataCommunicator;
    if (updatePartitioning) {
        dataCommunicator = m_kernel->createDataCommunicator();
    }
#endif

    // Update kernel
    m_kernel->update( adaptionData ) ;

    // Create sign propagator
    std::unique_ptr<LevelSetSignPropagator> signPropagator ;
    if (m_propagateSign) {
        signPropagator = m_kernel->createSignPropagator() ;
    }

    // Update objects
    for( int id : m_objectsProcessingOrder){
        LevelSetObject *object = m_objects.at(id).get() ;
        if (objectProcessList.count(object) == 0) {
            continue;
        }

        LevelSetSignedObjectInterface *signPropagationObject;
        if (m_propagateSign) {
            signPropagationObject = dynamic_cast<LevelSetSignedObjectInterface *>(object);
        } else {
            signPropagationObject = nullptr;
        }

#if BITPIT_ENABLE_MPI
        // Start partitioning update
        if (updatePartitioning) {
            object->startExchange( partitioningSendList, dataCommunicator.get() );
        }
#endif

        // Propagated sign is now dirty
        //
        // The propagated sign will be set as non-dirty by the sign propagator.
        if (signPropagationObject) {
            signPropagationObject->setSignStorageDirty(true);
        }

        // Clear data structures after mesh update
        object->clearAfterMeshAdaption( adaptionData ) ;

#if BITPIT_ENABLE_MPI
        // Complete partitioning update
        if (updatePartitioning) {
            object->completeExchange( partitioningRecvList, dataCommunicator.get() );
        }
#endif

        // Update levelset inside narrow band
        if (updateNarrowBand) {
            object->updateNarrowBand( adaptionData, m_signedDistance ) ;
        }

#if BITPIT_ENABLE_MPI
        // Update data on ghost cells
        object->exchangeGhosts();
#endif

        // Propagate sign
        //
        // It's not possible to communicate sign information, therefore sign
        // needs to be propagated also when the mesh is only partitioned.
        if (signPropagationObject) {
            signPropagator->execute(adaptionData, signPropagationObject);
        }
    }

}

#if BITPIT_ENABLE_MPI
/*!
 * Updates the levelset after mesh partitioning.
 * @param[in] adaptionData are the information about the adaption
 */
void LevelSet::partition( const std::vector<adaption::Info> &adaptionData ){

    update ( adaptionData ) ;

}
#endif

/*!
 * Return the list of all objects that can be processed.
 * @return the list of all objects that can be processed.
 */
std::unordered_set<LevelSetObject *> LevelSet::getObjectProcessList() const{

    std::unordered_set<LevelSetObject *> objectProcessList;
    for (const auto &entry : m_objects) {
        objectProcessList.insert(entry.second.get()) ;
    }

    return objectProcessList;

}

/*!
 * Return the list of objects that need to be processed in order to evaluate
 * the levelset of the specified objects.
 * The levelset of an objects may depend on the levelset of other objects,
 * therefore to evaluate the levelset of the specified object it may be
 * needed to process also other objects.
 * @param[in] nObjects is the number of objects.
 * @param[in] objectIds are the object identifiers.
 * @return the list of objects that need to be processed in order to
 * evaluate the lsevelset of the specified objects.
 */
std::unordered_set<LevelSetObject *> LevelSet::getObjectProcessList(std::size_t nObjects, const int *objectIds) const{

    std::unordered_set<LevelSetObject *> objectProcessList;
    for (std::size_t i = 0; i < nObjects; ++i) {
        int objectId = objectIds[i] ;
        LevelSetObject *object = getObjectPtr(objectId) ;
        if (objectProcessList.count(object) != 0) {
            continue;
        }

        objectProcessList.insert(object);
        if( const LevelSetProxyObject *proxyObject = dynamic_cast<const LevelSetProxyObject*>(object) ){
            std::vector<int> sourceObjectIds = proxyObject->getSourceObjectIds();
            std::unordered_set<LevelSetObject *> sourceObjectProcessList = getObjectProcessList(sourceObjectIds.size(), sourceObjectIds.data()) ;
            objectProcessList.insert(sourceObjectProcessList.begin(), sourceObjectProcessList.end());
        }
    }

    return objectProcessList;

}

/*! 
 * Writes LevelSetKernel to stream in binary format
 * @param[in] stream output stream
 */
void LevelSet::dump( std::ostream &stream ){

    m_objectIdentifierGenerator.dump(stream);

    utils::binary::write(stream, m_storageType);

    utils::binary::write(stream, m_objectsProcessingOrder);
    utils::binary::write(stream, m_narrowBandSize);
    utils::binary::write(stream, m_signedDistance);
    utils::binary::write(stream, m_propagateSign);

    for( const auto &object : m_objects ){
        object.second->dump( stream ) ;
    }
}

/*! 
 * Reads LevelSetKernel from stream in binary format
 * @param[in] stream output stream
 */
void LevelSet::restore( std::istream &stream ){

    m_objectIdentifierGenerator.restore(stream);

    utils::binary::read(stream, m_storageType);

    utils::binary::read(stream, m_objectsProcessingOrder);
    utils::binary::read(stream, m_narrowBandSize);
    utils::binary::read(stream, m_signedDistance);
    utils::binary::read(stream, m_propagateSign);

    for( const auto &object : m_objects ){
        object.second->restore( stream ) ;
    }
}

}

