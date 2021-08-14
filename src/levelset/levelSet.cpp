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

# include "levelSetCommon.hpp"
# include "levelSetKernel.hpp"
# include "levelSetCartesian.hpp"
# include "levelSetOctree.hpp"
# include "levelSetObject.hpp"
# include "levelSetMetaObject.hpp"
# include "levelSetCachedObject.hpp"
# include "levelSetBoolean.hpp"
# include "levelSetSegmentation.hpp"
# include "levelSetSignPropagator.hpp"
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

    m_narrowBandSize = levelSetDefaults::NARROWBAND_SIZE;

    m_signedDistance = true ;
    m_propagateSign  = false;

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
    if (m_kernel) {
        throw std::runtime_error ("Mesh can be set just once.");
    }

    LevelSetKernel *kernel = new LevelSetCartesian( *cartesian) ;
    m_kernel = std::unique_ptr<LevelSetKernel>(kernel);

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
    if (m_kernel) {
        throw std::runtime_error ("Mesh can be set just once.");
    }

    LevelSetKernel *kernel = new LevelSetOctree( *octree) ;
    m_kernel = std::unique_ptr<LevelSetKernel>(kernel);

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
int LevelSet::addObject( LevelSetBooleanOperation operation, int id1, int id2, int id ) {

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
int LevelSet::addObject( LevelSetBooleanOperation operation, const std::vector<int> &ids, int id ) {

    std::vector<const LevelSetObject*> ptr;
    for( int id : ids){
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

    assert(m_kernel && " levelset: setMesh must be called befor adding a LevelSetMask object ");

    return registerObject( std::unique_ptr<LevelSetObject>( new LevelSetMask(id, list, *m_kernel->getMesh()) ) );
}

/*!
 * Adds a LevelSetMask object composed of a list of interfaces
 * The function setMesh() must have been called prior.
 * @param[in] list list of indices of interfaces
 * @param[in] refInterface id of reference interface
 * @param[in] invert if orientation should be inverted with respect to the reference interface
 * @param[in] id id to be assigned to object. In case default value is passed the insertion order will be used as identifier
 * @return identifier of new object
 */
int LevelSet::addObject( const std::vector<long> &list, long refInterface, bool invert, int id ) {

    assert(m_kernel && " levelset: setMesh must be called befor adding a LevelSetMask object ");

    return registerObject( std::unique_ptr<LevelSetObject>( new LevelSetMask(id, list, refInterface, invert, *m_kernel->getMesh())) );
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
 * If an object is referenced by other objects it cannot be deleted.
 *
 * @param[in] id id of object to be removed
 * @return true if object has been found and removed
 */
bool LevelSet::removeObject(int id) {
    if( m_objects.count(id) != 0){
        m_objectIdentifierGenerator.trash(id);
        m_objects.erase(id);
        unsetObjectProcessingOrder(id);
        return true;
    } 

    return false;
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

    std::vector<int>::iterator processingOrderBegin = m_objectsProcessingOrder.begin();
    std::vector<int>::iterator processingOrderEnd   = m_objectsProcessingOrder.end();

    std::vector<int>::iterator processingOrderItr = std::find(processingOrderBegin, processingOrderEnd, id);
    assert(processingOrderItr != processingOrderEnd);
    m_objectsProcessingOrder.erase(processingOrderItr);

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
    m_kernel.reset();
    removeObjects();
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
 * This routines needs to be called at least once.
 * Each object should compute the levelset and associated
 * information on both internal and ghost cells.
 */
void LevelSet::compute(){

    assert(m_kernel && "LevelSet::setMesh() must be called prior to LevelSet::compute()");

    std::unique_ptr<LevelSetSignPropagator> signPropagator ;
    if (m_propagateSign) {
        signPropagator = std::unique_ptr<LevelSetSignPropagator>(new LevelSetSignPropagator(m_kernel->getMesh())) ;
    }

    for( int id : m_objectsProcessingOrder){
        LevelSetObject *object = m_objects.at(id).get() ;
        object->computeLSInNarrowBand(m_signedDistance) ;
#if BITPIT_ENABLE_MPI
        object->exchangeGhosts();
#endif
        if(m_propagateSign) {
            LevelSetSignStorage *signStorage = dynamic_cast<LevelSetSignStorage *>(object);
            if (signStorage) {
                signPropagator->execute(object, signStorage);
            }
        }
    }
}

/*!
 * Updates the levelset after a mesh update.
 * Each object should compute the levelset and associated
 * information on both internal and ghost cells.
 *
 * @param[in] mapper mapper conatining mesh modifications
 */
void LevelSet::update( const std::vector<adaption::Info> &mapper ){

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

    // Inspect mapper to detect what needs to be updated
    std::unordered_map<int,std::vector<long>> partitioningSendList ;
    std::unordered_map<int,std::vector<long>> partitioningRecvList ;

    bool updateNarrowBand   = false;
    bool updatePartitioning = false;
    for( const auto &event : mapper){
        if( event.entity != adaption::Entity::ENTITY_CELL){
            continue;
        }

        if( event.type == adaption::Type::TYPE_PARTITION_SEND){
            partitioningSendList.insert({{event.rank,event.previous}}) ;
            updatePartitioning = true;
        } else if( event.type == adaption::Type::TYPE_PARTITION_RECV){
            partitioningRecvList.insert({{event.rank,event.current}}) ;
            updatePartitioning = true;
        } else {
            if (!updateNarrowBand) {
                for (long cellId : event.current) {
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

#if BITPIT_ENABLE_MPI
    // Get data communicator for updating partitioning
    std::unique_ptr<DataCommunicator> dataCommunicator;
    if (updatePartitioning) {
        dataCommunicator = m_kernel->createDataCommunicator();
    }
#endif

    // Update kernel
    m_kernel->updateGeometryCache( mapper ) ;

    // Create sign propagator
    std::unique_ptr<LevelSetSignPropagator> signPropagator ;
    if (m_propagateSign) {
        signPropagator = std::unique_ptr<LevelSetSignPropagator>(new LevelSetSignPropagator(m_kernel->getMesh())) ;
    }

    // Update objects
    for( int id : m_objectsProcessingOrder){
        LevelSetObject *object = m_objects.at(id).get() ;

#if BITPIT_ENABLE_MPI
        // Start partitioning update
        if (updatePartitioning) {
            object->startExchange( partitioningSendList, dataCommunicator.get() );
        }
#endif

        // Clear data structures after mesh update
        object->clearAfterMeshAdaption( mapper ) ;

#if BITPIT_ENABLE_MPI
        // Complete partitioning update
        if (updatePartitioning) {
            object->completeExchange( partitioningRecvList, dataCommunicator.get() );
        }
#endif

        // Update levelset inside narrow band
        if (updateNarrowBand) {
            object->updateLSInNarrowBand( mapper, m_signedDistance ) ;
        }

#if BITPIT_ENABLE_MPI
        // Update data on ghost cells
        object->exchangeGhosts();
#endif

        // Propagate sign
        //
        // It's not possible to communicate sign information, therefore sign
        // needs to be propagated also when the mesh is only partitioned.
        if (m_propagateSign) {
            LevelSetSignStorage *signStorage = dynamic_cast<LevelSetSignStorage *>(object);
            if (signStorage) {
                signPropagator->execute(mapper, object, signStorage);
            }
        }
    }

}

#if BITPIT_ENABLE_MPI
/*!
 * Updates the levelset after mesh partitioning.
 * @param[in] mapper mapper conatining mesh modifications
 */
void LevelSet::partition( const std::vector<adaption::Info> &mapper ){

    update ( mapper ) ;

}
#endif

/*! 
 * Writes LevelSetKernel to stream in binary format
 * @param[in] stream output stream
 */
void LevelSet::dump( std::ostream &stream ){

    m_objectIdentifierGenerator.dump(stream);

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

    utils::binary::read(stream, m_objectsProcessingOrder);
    utils::binary::read(stream, m_narrowBandSize);
    utils::binary::read(stream, m_signedDistance);
    utils::binary::read(stream, m_propagateSign);

    for( const auto &object : m_objects ){
        object.second->restore( stream ) ;
    }
}

}

