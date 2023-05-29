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
# include "levelSetBooleanObject.hpp"
# include "levelSetSegmentationObject.hpp"
# include "levelSetMaskObject.hpp"
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
 */
LevelSet::LevelSet(LevelSetFillIn expectedFillIn) {

    m_expectedFillIn = expectedFillIn ;

    m_objects.clear() ;

    m_signedDistance = true ;
    m_propagateSign  = false;

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
        m_kernel = std::unique_ptr<LevelSetKernel>(new LevelSetCartesianKernel(*cartesian, m_expectedFillIn));
    } else if (VolOctree *octree = dynamic_cast<VolOctree *>(mesh)) {
        m_kernel = std::unique_ptr<LevelSetKernel>(new LevelSetOctreeKernel(*octree, m_expectedFillIn));
    } else if (VolUnstructured *unstructured = dynamic_cast<VolUnstructured*>(mesh)) {
        m_kernel = std::unique_ptr<LevelSetKernel>(new LevelSetUnstructuredKernel(*unstructured, m_expectedFillIn));
    } else {
        throw std::runtime_error ("Unable to create the levelset kernel. Mesh type non supported.");
    }

    // Assign the kernel to the existing objects
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

    auto object = std::unique_ptr<LevelSetObject>(new LevelSetSegmentationObject(id, std::move(segmentation), angle));

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

    auto object = std::unique_ptr<LevelSetObject>(new LevelSetSegmentationObject(id, segmentation, angle));

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

    auto object = std::unique_ptr<LevelSetObject>(new LevelSetSegmentationObject(id, std::move(surfUnstructuredUPtr), angle));

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

    auto object = std::unique_ptr<LevelSetObject>(new LevelSetSegmentationObject(id, surfUnstructured, angle));

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

    auto object = std::unique_ptr<LevelSetObject>(new LevelSetBooleanObject(id, operation, ptr1, ptr2));

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

    auto object = std::unique_ptr<LevelSetObject>(new LevelSetBooleanObject(id, operation, ptr));

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

    auto object = std::unique_ptr<LevelSetObject>(new LevelSetMaskObject(id, list, *m_kernel->getMesh()));

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

    auto object = std::unique_ptr<LevelSetObject>(new LevelSetMaskObject(id, list, refInterface, invert, *m_kernel->getMesh()));

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

    object->setDefaultLevelSetSigned(m_signedDistance);

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
    return unregisterObject(id, false);

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
bool LevelSet::unregisterObject(int id, bool force) {
    if( !force && !isObjectRemovable(id) ) {
        return false;
    }

    m_cacheFilledObjects.erase(id);

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
 * Set the processing order of the specified object.
 *
 * The insertion order determines the processing order, however priority is
 * given to primary objects.
 *
 * This function must be called when a new object is added.
 *
 * @param[in] id the id of the object
 */
void LevelSet::setObjectProcessingOrder( int id ) {

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
 * This function must be called whan a object is removed.
 * @param[in] id the id of the object
 */
void LevelSet::unsetObjectProcessingOrder(int id){

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
 * Fill cache.
 *
 * Cache information will be computed on both internal and ghost cells.
 */
void LevelSet::fillCache( ) {

    fillCache(m_objectsProcessingOrder);

}


/*!
 * Fill the cache of the specified object.
 *
 * Cache information will be computed on both internal and ghost cells.
 *
 * \param id is the id of the object
 */
void LevelSet::fillCache( int id ) {

    std::vector<int> ids(1, id);
    fillCache(ids);

}

/*!
 * Fill the cache of the specified objects.
 *
 * Cache information will be computed on both internal and ghost cells.
 *
 * \param ids are the ids of the objects
 */
void LevelSet::fillCache( const std::vector<int> &ids ) {

    for( int id : ids){
        LevelSetObject *object = m_objects.at(id).get() ;

        // Set cache needed for sign propagation
        if (m_propagateSign) {
            object->setFieldCache(LevelSetField::VALUE, LevelSetCacheMode::NARROW_BAND);
            object->setFieldCache(LevelSetField::SIGN, LevelSetCacheMode::FULL);
        }

        // Fill object cache
        object->fillCache();
        m_cacheFilledObjects.insert(id);
    }

}

/*!
 * Updates cached information after a mesh update.
 * Before being able to update the cached, it has to be filled.
 * Levelset and associated information will be updated on both internal and
 * ghost cells.
 * @param[in] adaptionData are the information about the adaption
 */
void LevelSet::updateCache( const std::vector<adaption::Info> &adaptionData ){

    assert(m_kernel && "LevelSet::setMesh() must be called prior to LevelSet::partition()");

    VolumeKernel *mesh = m_kernel->getMesh();
#if BITPIT_ENABLE_MPI
    bool isMeshPartitioned = mesh->isPartitioned();
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

    // Update levelset
    if (updateNarrowBand || updatePartitioning) {
#if BITPIT_ENABLE_MPI
        // Get data communicator for updating partitioning
        std::unique_ptr<DataCommunicator> dataCommunicator;
        if (updatePartitioning) {
            dataCommunicator = m_kernel->createDataCommunicator();
        }
#endif

        // Update kernel
        m_kernel->update( adaptionData ) ;

        // Update objects whose cache is already filled
        for( int id : m_cacheFilledObjects){
            LevelSetObject *object = m_objects.at(id).get() ;
#if BITPIT_ENABLE_MPI
            // Start partitioning update
            if (updatePartitioning) {
                object->startExchange( partitioningSendList, dataCommunicator.get() );
            }
#endif

            // Clear data structures after mesh update
            object->updateCache( adaptionData ) ;

#if BITPIT_ENABLE_MPI
            // Complete partitioning update
            if (updatePartitioning) {
                object->completeExchange( partitioningRecvList, dataCommunicator.get() );
            }
#endif

#if BITPIT_ENABLE_MPI
            // Update data on ghost cells
            object->exchangeGhosts();
#endif
        }
    }

    // Fill cache of objects whose cache is empty
    std::vector<int> emptyCacheObjects;
    for (int id : m_objectsProcessingOrder) {
        if (m_cacheFilledObjects.count(id) == 0) {
            emptyCacheObjects.emplace_back(id);
        }
    }

    if (!emptyCacheObjects.empty()) {
        fillCache(emptyCacheObjects);
    }

}

/*!
 * Clear cached information.
 *
 * Both information stored in objects and in the kernel will be cleared.
 *
 * \param release if set to true the memory hold by the cache will be released, otherwise
 * the cache will be cleared but its memory may not be released
 */
void LevelSet::clearCache( bool release ) {

    // Clear cache of all the objects
    for( int id : m_objectsProcessingOrder){
        LevelSetObject *object = m_objects.at(id).get() ;
        object->clearCache(release);
    }

    // Clear mesh cache
    LevelSetCachedKernel *cachedKernel = dynamic_cast<LevelSetCachedKernel *>(m_kernel.get());
    if (!cachedKernel) {
        return;
    }

    cachedKernel->clearCache(release);

}

/*!
 * Set if the levelset sign has to be propagated from the narrow band to the whole domain.
 *
 * This function is provided only for compatibility with older versions of bitpit. It set
 * the sign cache to "full" mode and the value cache to "narrow band" for all the objects.
 * The recommended way to setup sign propagation is to set caches of the objects directly.
 *
 * @param[in] flag True/false to active/disable the propagation .
 */

void LevelSet::setPropagateSign(bool flag){
    m_propagateSign = flag;
}

/*!
 * Get the physical size of the narrow band.
 * The function will return the minimum narrow band size among all the objects.
 * @return the physical size of the narrow band
 */
double LevelSet::getSizeNarrowBand() const{
    double size = std::numeric_limits<int>::max();
    for( int id : m_objectsProcessingOrder){
        LevelSetObject *object = m_objects.at(id).get() ;
        size = std::min(object->getNarrowBandSize(), size);
    }

    return size;
}

/*!
 * Manually set the physical size of the narrow band.
 * @param[in] r Size of the narrow band.
 */
void LevelSet::setSizeNarrowBand(double r){
    for( int id : m_objectsProcessingOrder){
        LevelSetObject *object = m_objects.at(id).get() ;
        object->setNarrowBandSize(r);
    }
}

/*!
 * Computes levelset on given mesh with respect to the objects.
 * Levelset and associated information will be computed on both internal and
 * ghost cells.
 */
void LevelSet::compute(){

    fillCache();

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

    BITPIT_UNUSED(id);

    log::warning() << " It is not possible to compute the levelset for a single object." << std::endl;
    log::warning() << " Levelset will be computed for all the objects." << std::endl;

    fillCache();


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

    BITPIT_UNUSED(ids);

    log::warning() << " It is not possible to compute the levelset for a single object." << std::endl;
    log::warning() << " Levelset will be computed for all the objects." << std::endl;

    fillCache();

}

/*!
 * Updates the levelset after a mesh update.
 * Before being able to update the levelset, it has to be computed.
 * Levelset and associated information will be updated on both internal and
 * ghost cells.
 * @param[in] adaptionData are the information about the adaption
 */
void LevelSet::update( const std::vector<adaption::Info> &adaptionData ){

    log::warning() << " It is not possible to compute the levelset for a single object." << std::endl;
    log::warning() << " Levelset will be computed for all the objects." << std::endl;

    updateCache(adaptionData);

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

    BITPIT_UNUSED(id);

    log::warning() << " It is not possible to update the levelset for a single object." << std::endl;
    log::warning() << " Levelset will be computed for all the objects." << std::endl;

    updateCache(adaptionData);

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

    BITPIT_UNUSED(ids);

    log::warning() << " It is not possible to update the levelset for a single object." << std::endl;
    log::warning() << " Levelset will be computed for all the objects." << std::endl;

    updateCache(adaptionData);

}

#if BITPIT_ENABLE_MPI
/*!
 * Updates the levelset after mesh partitioning.
 * @param[in] adaptionData are the information about the adaption
 */
void LevelSet::partition( const std::vector<adaption::Info> &adaptionData ){

    updateCache ( adaptionData ) ;

}
#endif

/*! 
 * Writes LevelSetKernel to stream in binary format
 * @param[in] stream output stream
 */
void LevelSet::dump( std::ostream &stream ){

    utils::binary::write(stream, m_expectedFillIn);
    utils::binary::write(stream, m_signedDistance);
    utils::binary::write(stream, m_propagateSign);

    m_objectIdentifierGenerator.dump(stream);
    for( const auto &object : m_objects ){
        object.second->dump( stream ) ;
    }
    utils::binary::write(stream, m_objectsProcessingOrder);

    utils::binary::write(stream, m_cacheFilledObjects.size());
    for( long id : m_cacheFilledObjects ){
        utils::binary::write(stream, id);
    }

}

/*! 
 * Reads LevelSetKernel from stream in binary format
 * @param[in] stream output stream
 */
void LevelSet::restore( std::istream &stream ){

    utils::binary::read(stream, m_expectedFillIn);
    utils::binary::read(stream, m_signedDistance);
    utils::binary::read(stream, m_propagateSign);

    m_objectIdentifierGenerator.restore(stream);
    for( const auto &object : m_objects ){
        object.second->restore( stream ) ;
    }
    utils::binary::read(stream, m_objectsProcessingOrder);

    std::size_t nCacheFilledObjects;
    utils::binary::read(stream, nCacheFilledObjects);
    for (std::size_t i = 0; i < nCacheFilledObjects; ++i) {
        long id;
        utils::binary::read(stream, id);
        m_cacheFilledObjects.insert(id);
    }

}

}

