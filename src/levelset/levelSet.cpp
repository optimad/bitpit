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
# include "levelSetComplementObject.hpp"
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
 * Evaluation of the fields inside the narrow band is always performed using an exact algorithm,
 * on the other hand evaluation of the fields in the bulk can be performed choosing different
 * algorithms.
 *
 * With respect to the levelset value, the domain can be divided in two regions: the narrow band
 * and the bulk. The narrow band defines the portion of the domain close to the zero-levelset iso
 * surface whereas the bulk is everything else. Regardless of the specified narrow bans size, the
 * narrow band will always contain the intersected cells and their neighbours.
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

    m_forceSignPropagation   = false;
    m_signPropagationEnabled = false;

    m_narrowBandSize = 0;

}

/*!
 * Clear the levelset entirely, deleting kernel and objects
 */
void LevelSet::clear(){
    removeObjects();
    m_kernel.reset();
}

/*!
 * Updates the levelset a mesh update.
 * Before being able to update the cached, it has to be filled.
 * Levelset and associated information will be updated on both internal and
 * ghost cells.
 * @param[in] adaptionData are the information about the adaption
 */
void LevelSet::update( const std::vector<adaption::Info> &adaptionData ){

    assert(m_kernel && "LevelSet::setMesh() must be called prior to LevelSet::update()");

    // Update kernel
    bool updated = m_kernel->update( adaptionData ) ;
    if (!updated) {
        return;
    }

    // Update objects
    for( int id : m_orderedObjectsIds){
        LevelSetObject *object = m_objects.at(id).get() ;
        object->update( adaptionData ) ;
    }
}

/*!
 * Get the size of the narrow band.
 *
 * With respect to the levelset value, the domain can be divided in two regions: the narrow band
 * and the bulk. The narrow band defines the portion of the domain close to the zero-levelset iso
 * surface whereas the bulk is everything else.
 *
 * The size of the narrow band is the absolute distance from the zero-levelset iso surface below
 * which a point is considered belonging to the narrow band. Setting the size of the narrow band
 * to LEVELSET_NARROW_BAND_UNLIMITED means that the whole domain belongs to the narrow band.
 * Regardless of the specified size, the narrow band will always contain the intersected cells
 * and their neighbours.
 *
 * @return The size of the narrow band.
 */
double LevelSet::getNarrowBandSize() const {
    return m_narrowBandSize;
}

/*!
 * Set the size of the narrow band.
 *
 * With respect to the levelset value, the domain can be divided in two regions: the narrow band
 * and the bulk. The narrow band defines the portion of the domain close to the zero-levelset iso
 * surface whereas the bulk is everything else.
 *
 * The size of the narrow band is the absolute distance from the zero-levelset iso surface below
 * which a point is considered belonging to the narrow band. Setting the size of the narrow band
 * to LEVELSET_NARROW_BAND_UNLIMITED means that the whole domain belongs to the narrow band.
 * Regardless of the specified size, the narrow band will always contain the intersected cells
 * and their neighbours.
 *
 * @param[in] size is the size of the narrow band.
 */
void LevelSet::setNarrowBandSize(double size) {
    // Set narrow band size
    m_narrowBandSize = std::max(size, 0.);

    // Set object narrow band size
    for (auto &objectEntry : m_objects) {
        LevelSetObject *object = objectEntry.second.get();
        object->setNarrowBandSize(size);
    }
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
    for( int id : m_orderedObjectsIds){
        LevelSetObject *object = m_objects.at(id).get() ;
        object->setKernel(m_kernel.get());
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

    auto surfUnstructured = std::unique_ptr<SurfUnstructured>(dynamic_cast<SurfUnstructured *>(segmentation.release())) ;
    if (!surfUnstructured) {
        throw std::runtime_error ("Segmentation type not supported");
    }

    auto object = std::unique_ptr<LevelSetObject>(new LevelSetSegmentationObject(id, std::move(surfUnstructured), angle));

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
 * Adds a LevelSetMask object composed of the external envelope of a list of mesh cells.
 * Objects can be added to the levelset only after setting the mesh.
 * @param[in] list list of indices of cells
 * @param[in] id id to be assigned to object. In case default value is passed the insertion order will be used as identifier
 * @return identifier of new object
 */
int LevelSet::addObject( const std::unordered_set<long> &list, int id ) {

    assert(m_kernel && " levelset: setMesh must be called before adding a mask object ");

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

    assert(m_kernel && " levelset: setMesh must be called before adding a mask object ");

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

    // Set object id
    int id = object->getId();
    if (id == levelSetDefaults::OBJECT) {
        id = m_objectIdentifierGenerator.generate();
        object->setId(id);
    } else {
        m_objectIdentifierGenerator.setAssigned(id);
    }

    // Set object properties
    object->setDefaultLevelSetSigndness(m_signedDistance);
    object->setNarrowBandSize(m_narrowBandSize);
    if (m_forceSignPropagation) {
        if (m_signPropagationEnabled) {
            object->setCellBulkEvaluationMode(LevelSetBulkEvaluationMode::SIGN_PROPAGATION);
        } else {
            object->setCellBulkEvaluationMode(LevelSetBulkEvaluationMode::NONE);
        }
    }

    // Set object kernel
    if (m_kernel) {
        object->setKernel(m_kernel.get());
    }

    // Add the object to the levelset
    m_objects[id] = std::move(object) ;
    registerObjectId(id);
    incrementObjectsReferenceCount(id);

    return id;
}

/*!
 * Remove all levelset objects
 */
void LevelSet::removeObjects() {
    m_orderedObjectsIds.clear();
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

    decrementObjectsReferenceCount(id);
    unregisterObjectId(id);
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
 * Register to specified object id.
 *
 * Objects should be processed in a specific order: first the primary objects (in the order they
 * were added to the levelset) and then the other objects (in the order they were added to the
 * levelset). The registration process updates the list of ids to guarantee that the objects
 * will be processed in the correct order.
 *
 * This function must be called when a new object is added.
 *
 * @param[in] id the object id that will be registered
 */
void LevelSet::registerObjectId( int id ) {

    // Add the id from the list
    //
    // Ids should be sorted according to the order in which they corresponding objects should be
    // processed: first the ids of the primary objects (in the order the objects were added to
    // the levelset) and then the ids of the other objects (in the order the objects were added
    // to the levelset).
    std::vector<int>::iterator idItr;
    if(getObjectPtr(id)->isPrimary()){
        std::vector<int>::iterator orderedIdsBegin = m_orderedObjectsIds.begin();
        std::vector<int>::iterator orderedIdsEnd   = m_orderedObjectsIds.end();

        for (idItr = orderedIdsBegin; idItr != orderedIdsEnd; ++idItr) {
            int candidateId = *idItr ;
            const LevelSetObject *candidateObject = getObjectPtr(candidateId) ;
            if( !candidateObject->isPrimary() ){
                break;
            }
        }
    } else {
        idItr = m_orderedObjectsIds.end();
    }

    m_orderedObjectsIds.insert(idItr, id) ;

}

/*!
 * Register to specified object id.
 *
 * This function must be called when a object is removed.
 *
 * @param[in] id the object id that will be unregistered
 */
void LevelSet::unregisterObjectId(int id){

    // Remove the id from the list
    std::vector<int>::iterator orderedIdsBegin = m_orderedObjectsIds.begin();
    std::vector<int>::iterator orderedIdsEnd   = m_orderedObjectsIds.end();

    std::vector<int>::iterator idItr = std::find(orderedIdsBegin, orderedIdsEnd, id);
    assert(idItr != orderedIdsEnd);
    m_orderedObjectsIds.erase(idItr);

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

    if( const LevelSetProxyBaseObject *parentProxyObject = dynamic_cast<const LevelSetProxyBaseObject *>(parentObject) ){
        for ( int sourceObjectId : parentProxyObject->getSourceObjectIds() ){
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

    if( const LevelSetProxyBaseObject *parentProxyObject = dynamic_cast<const LevelSetProxyBaseObject *>(parentObject) ){
        for ( int sourceObjectId : parentProxyObject->getSourceObjectIds() ){
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
 * Writes LevelSetKernel to stream in binary format
 * @param[in] stream output stream
 */
void LevelSet::dump( std::ostream &stream ) const{

    utils::binary::write(stream, m_expectedFillIn);
    utils::binary::write(stream, m_signedDistance);
    utils::binary::write(stream, m_forceSignPropagation);
    utils::binary::write(stream, m_signPropagationEnabled);
    utils::binary::write(stream, m_narrowBandSize);

    m_objectIdentifierGenerator.dump(stream);
    for( const auto &object : m_objects ){
        object.second->dump( stream ) ;
    }
    utils::binary::write(stream, m_orderedObjectsIds);

}

/*!
 * Reads LevelSetKernel from stream in binary format
 * @param[in] stream output stream
 */
void LevelSet::restore( std::istream &stream ){

    utils::binary::read(stream, m_expectedFillIn);
    utils::binary::read(stream, m_signedDistance);
    utils::binary::read(stream, m_forceSignPropagation);
    utils::binary::read(stream, m_signPropagationEnabled);
    utils::binary::read(stream, m_narrowBandSize);

    m_objectIdentifierGenerator.restore(stream);
    for( const auto &object : m_objects ){
        object.second->restore( stream ) ;
    }
    utils::binary::read(stream, m_orderedObjectsIds);

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
 *
 * This function is provided only for compatibility with older versions of bitpit. It sets
 * the bulk evaluation mode that matches the behaviour of the older levelset versions:
 *  - if sign propagation is enabled, the bulk evaluation mode is set to "sign propagation";
 *  - if sign propagation is disabled, the bulk evaluation mode is set to "none".
 * The recommended way to setup sign propagation is to manually set the bulk evaluation mode
 * of the relevant objects to "sign propagation".
 *
 * @param[in] flag True/false to active/disable the propagation .
 */
void LevelSet::setPropagateSign(bool flag){
    // Set propagation
    m_forceSignPropagation   = true;
    m_signPropagationEnabled = flag;

    // Set object bulk evaluation mode
    for (auto &objectEntry : m_objects) {
        LevelSetObject *object = objectEntry.second.get();
        if (m_forceSignPropagation) {
            if (m_signPropagationEnabled) {
                object->setCellBulkEvaluationMode(LevelSetBulkEvaluationMode::SIGN_PROPAGATION);
            } else {
                object->setCellBulkEvaluationMode(LevelSetBulkEvaluationMode::NONE);
            }
        }
    }
}

/*!
 * Get the size of the narrow band.
 *
 * With respect to the levelset value, the domain can be divided in two regions: the narrow band
 * and the bulk. The narrow band defines the portion of the domain close to the zero-levelset iso
 * surface whereas the bulk is everything else.
 *
 * The size of the narrow band is the absolute distance from the zero-levelset iso surface below
 * which a point is considered belonging to the narrow band. Setting the size of the narrow band
 * to LEVELSET_NARROW_BAND_UNLIMITED means that the whole domain belongs to the narrow band.
 * Regardless of the specified size, the narrow band will always contain the intersected cells
 * and their neighbours.
 *
 * @return The size of the narrow band.
 */
double LevelSet::getSizeNarrowBand() const{
    return getNarrowBandSize();
}

/*!
 * Set the size of the narrow band.
 *
 * With respect to the levelset value, the domain can be divided in two regions: the narrow band
 * and the bulk. The narrow band defines the portion of the domain close to the zero-levelset iso
 * surface whereas the bulk is everything else.
 *
 * The size of the narrow band is the absolute distance from the zero-levelset iso surface below
 * which a point is considered belonging to the narrow band. Setting the size of the narrow band
 * to LEVELSET_NARROW_BAND_UNLIMITED means that the whole domain belongs to the narrow band.
 * Regardless of the specified size, the narrow band will always contain the intersected cells
 * and their neighbours.
 *
 * @param[in] size is the size of the narrow band.
 */
void LevelSet::setSizeNarrowBand(double size){
    return setNarrowBandSize(size);
}

/*!
 * Computes levelset on given mesh with respect to the objects.
 *
 * This function is deprecated, there is no need to explicitly compute levelset information.
 */
void LevelSet::compute(){

    // Nothing to do

}

/*!
 * Computes levelset on given mesh with respect to the specified object.
 *
 * This function is deprecated, there is no need to explicitly compute levelset information.
 *
 * @param[in] id identifier of object.
 */
void LevelSet::compute( int id ){

    BITPIT_UNUSED(id);

    // Nothing to do

}

/*!
 * Computes levelset on given mesh with respect to the specified objects.
 *
 * This function is deprecated, there is no need to explicitly compute levelset information.
 *
 * @param[in] ids identifiers of objects.
 */
void LevelSet::compute( const std::vector<int> &ids ){

    BITPIT_UNUSED(ids);

    // Nothing to do

}

/*!
 * Computes levelset on given mesh with respect to the specified object.
 *
 * It is not possible to update the levelset for a specific objects, levelset will be
 * computed for all the objects. This function is provided only for compatibility
 * with older versions of bitpit.
 *
 * @param[in] adaptionData are the information about the adaption
 * @param[in] id identifier of object.
 */
void LevelSet::update( const std::vector<adaption::Info> &adaptionData, int id){

    BITPIT_UNUSED(id);

    log::warning() << " It is not possible to update the levelset for a specific objects." << std::endl;
    log::warning() << " Levelset will be computed for all the objects." << std::endl;

    update(adaptionData);

}

/*!
 * Computes levelset on given mesh with respect to the specified objects.
 *
 * It is not possible to update the levelset for a specific objects, levelset will be
 * computed for all the objects. This function is provided only for compatibility
 * with older versions of bitpit.
 *
 * @param[in] adaptionData are the information about the adaption
 * @param[in] ids identifiers of objects.
 */
void LevelSet::update( const std::vector<adaption::Info> &adaptionData, const std::vector<int> &ids ){

    BITPIT_UNUSED(ids);

    log::warning() << " It is not possible to update the levelset for a specific objects." << std::endl;
    log::warning() << " Levelset will be computed for all the objects." << std::endl;

    update(adaptionData);

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

}

