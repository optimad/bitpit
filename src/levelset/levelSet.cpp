/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

# if BITPIT_ENABLE_MPI
# include <mpi.h>
# endif

# include <unordered_set>

# include "bitpit_SA.hpp"
# include "bitpit_operators.hpp"

# include "levelSet.hpp"

namespace bitpit {

/*!
 * @ingroup levelset
 * @class  LevelSetInfo
 *
 * @brief  A public container which includes all information provided by LevelSet
 *
 * LevelSetInfo conatins the following information
 * - distance to closest object
 * - gradient of level set function
 * - the id of the closest object
 * - the patch of closest object which contains the projection point
 * - the id of the surface element closest to the projection point 
 *
 * If a grid point lies within the narrow band of an object, all of these information are available, provided that the clsest object may provide them.
 * On the contrary, if the grid point lies were the levelset has been propagated not all of these information may be available; 
 * in particular the object, patch and segment may be set to the default values.
 *
*/

/*!
 * Default constructor
 */
LevelSetInfo::LevelSetInfo() :value(levelSetDefaults::VALUE), gradient(levelSetDefaults::GRADIENT) {
}

/*!
 * Complete constructor
 * @param[in] v value of levelset function
 * @param[in] g gradient of levelset function
 */
LevelSetInfo::LevelSetInfo(const double &v, const std::array<double,3> &g) :value(v), gradient(g) {
}

/*!
 * @ingroup levelset
 * @class  LevelSet
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

    m_object.clear() ;

    m_userRSearch = false ;

    m_signedDF    = true ;
    m_propagateS  = false;

};

/*!
 * Destructor of LevelSet
*/
LevelSet::~LevelSet(){

    m_object.clear() ;

};

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
        log::cout() << "Mesh non supported in LevelSet::setMesh()" << std::endl ;
    }; 

    return;
};

/*!
 * Sets the grid on which the levelset function should be computed.
 * @param[in] cartesian cartesian patch
 */
void LevelSet::setMesh( VolCartesian* cartesian ) {

    LevelSetKernel *kernel = new LevelSetCartesian( *cartesian) ;

    m_kernel = unique_ptr<LevelSetKernel>(kernel);

    return;
};

/*!
 * Sets the grid on which the levelset function should be computed.
 * @param[in] octree octree patch
 */
void LevelSet::setMesh( VolOctree* octree ) {

    LevelSetKernel *kernel = new LevelSetOctree( *octree) ;

    m_kernel = unique_ptr<LevelSetKernel>(kernel);

    return;
};

/*!
 * Adds a surface segmentation
 * @param[in] segmentation surface segmentation
 * @param[in] angle feature angle
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( std::unique_ptr<SurfUnstructured> &&segmentation, double angle, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_object.size();
    }

    LevelSetSegmentation* lsSeg = new LevelSetSegmentation(id, std::move(segmentation), angle ) ;

    LevelSetObject *object = static_cast<LevelSetObject *>(lsSeg);

    return addObject(std::unique_ptr<LevelSetObject>(object));
};

/*!
 * Adds a surface segmentation
 * @param[in] segmentation surface segmentation
 * @param[in] angle feature angle
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( SurfUnstructured *segmentation, double angle, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_object.size();
    }

    LevelSetSegmentation* lsSeg = new LevelSetSegmentation(id, segmentation, angle) ;

    LevelSetObject *object = static_cast<LevelSetObject *>(lsSeg);

    return addObject(std::unique_ptr<LevelSetObject>(object));
};

/*!
 * Adds a surface segmentation
 * @param[in] segmentation surface segmentation
 * @param[in] angle feature angle
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( std::unique_ptr<SurfaceKernel> &&segmentation, double angle, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_object.size();
    }

    LevelSetSegmentation* lsSeg = new LevelSetSegmentation(id,angle) ;
    if( SurfUnstructured* surfUnstructured = dynamic_cast<SurfUnstructured*>(segmentation.get()) ){
        segmentation.release();
        std::unique_ptr<SurfUnstructured> surfUnstructuredUPtr = std::unique_ptr<SurfUnstructured>(surfUnstructured) ;
        lsSeg->setSegmentation( std::move(surfUnstructuredUPtr) );
    } else {
        throw std::runtime_error ("Segmentation type not supported");
    }

    LevelSetObject *object = static_cast<LevelSetObject *>(lsSeg);

    return addObject(std::unique_ptr<LevelSetObject>(object));
};

/*!
 * Adds a surface segmentation
 * @param[in] segmentation surface segmentation
 * @param[in] angle feature angle
 * @param[in] id identifier of object; in case no id is provided the insertion
 * order will be used as identifier
 */
int LevelSet::addObject( SurfaceKernel *segmentation, double angle, int id ) {

    if (id == levelSetDefaults::OBJECT) {
        id = m_object.size();
    }

    LevelSetSegmentation* lsSeg = new LevelSetSegmentation(id,angle) ;
    if( SurfUnstructured* surfUnstructured = dynamic_cast<SurfUnstructured*>(segmentation)){
        lsSeg->setSegmentation( surfUnstructured );
    } else {
        throw std::runtime_error ("Segmentation type not supported");
    }

    LevelSetObject *object = static_cast<LevelSetObject *>(lsSeg);

    return addObject(std::unique_ptr<LevelSetObject>(object));
};

/*!
 * Adds a generic LevelSetObject
 * @param[in] object generic object
 * @return the index associated to the object
 */
int LevelSet::addObject( std::unique_ptr<LevelSetObject> &&object ) {

    int objectId = object->getId();
    m_object[objectId] = std::move(object) ;

    addProcessingOrder(objectId) ;

    return objectId;
};

/*!
 * Adds a generic LevelSetObject
 * @param[in] object generic object
 * @return the index associated to the object
 */
int LevelSet::addObject( const std::unique_ptr<LevelSetObject> &object ) {

    int objectId = object->getId();
    m_object[objectId] = std::unique_ptr<LevelSetObject>(object->clone())  ;

    addProcessingOrder(objectId) ;

    return objectId;
};


/*!
 * Adds the object to the processing order.
 *
 * The insertion order determines the processing order 
 * but priority is given to primary objects. 
 *
 * This function must be called whan a new object is inserted into m_objects.
 *
 * @param[in] objectId the id of the newly added object
 */
void LevelSet::addProcessingOrder( int objectId ) {

    bool primary = m_object.at(objectId)->isPrimary() ;

    if(primary){
        std::vector<int>::iterator orderItr = m_order.begin() ;
        bool iterate( orderItr != m_order.end()) ;

        while(iterate){
            int id = *orderItr ;
            if( m_object.at(id)->isPrimary() ){
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

    return ;
};


/*!
 * Get a constant reference to the specified object.
 * If hte specified id does not exist an exception is thrown.
 * @param id is the object id
 * @return pointer to levelset object
 */
const LevelSetObject & LevelSet::getObject( int id) const{
    return *(m_object.at(id)) ;
};

/*!
 * Get the number of levelset objects
 * @return number of objects
 */
int LevelSet::getObjectCount( ) const{
    return m_object.size() ;
};

/*!
 * Get the ids of the bodies.
 * @return a list of the body ids
*/
std::vector<int> LevelSet::getObjectIds( ) const{
    std::vector<int> ids ;
    ids.reserve(m_object.size()) ;
    for(const auto &entry : m_object) {
        ids.push_back(entry.first) ;
    }

    return ids ;
};

/*!
 * Clear LevelSet entirely
 */
void LevelSet::clear(){

    m_kernel.reset() ;

    m_object.clear() ;
    return ;
};

/*!
 * Clear all LevelSet Objects but leave kernel unaltered
 */
void LevelSet::clearObject(){
    m_object.clear() ;
    return ;
};

/*!
 * Clear only one LevelSet Object
 * @param[in] id Id of object to be deleted
 */
void LevelSet::clearObject(int id){

    m_object.erase(id) ;
    return ;
};

/*!
 * Get levelset information.
 * If levelSetDefaults::OBJECT is passed as objectId, the information of the closest object are returned, otherwise those of the indicated object.
 * @param[in] objectId index of object
 * @param[in] cellId index of cell
 * @return levelset information
 */
LevelSetInfo LevelSet::getLevelSetInfo( const long &cellId, int objectId)const {
    if( objectId == levelSetDefaults::OBJECT){
        objectId = getClosestObject(cellId) ;
    }

    return( m_object.at(objectId)->getLevelSetInfo(cellId) ) ;
};

/*!
 * Get the levelset value of the i-th local cell.
 * If levelSetDefaults::OBJECT is passed as objectId, the information of the closest object are returned, otherwise those of the indicated object.
 * @param[in] cellId index of cell
 * @param[in] objectId index of object
 * @return levelset value 
 */
double LevelSet::getLS( const long &cellId, int objectId)const {

    if( objectId == levelSetDefaults::OBJECT){
        objectId = getClosestObject(cellId) ;
    }

    if( objectId != levelSetDefaults::OBJECT ) {
        return m_object.at(objectId)->getLS(cellId) ;
    } else {
        return levelSetDefaults::VALUE ;
    }
};

/*!
 * Get the levelset gradient of the i-th local cell.
 * If levelSetDefaults::OBJECT is passed as objectId, the information of the closest object are returned, otherwise those of the indicated object.
 * @param[in] cellId index of cell
 * @param[in] objectId index of object
 * @return Array with components of the gradient 
 */
std::array<double,3> LevelSet::getGradient(const long &cellId, int objectId) const {
    if( objectId == levelSetDefaults::OBJECT){
        objectId = getClosestObject(cellId) ;
    }
    
    if( objectId != levelSetDefaults::OBJECT ) {
        return m_object.at(objectId)->getGradient(cellId) ;
    } else {
        return levelSetDefaults::GRADIENT ;
    }
};

/*!
 * Get the id of closest object
 * @param[in] cellId index of cell
 * @return id of closest object
 */
int LevelSet::getClosestObject(const long &cellId) const {
    double value(levelSetDefaults::VALUE) ;
    int objectId(levelSetDefaults::OBJECT) ;

    for( auto &object : m_object){
        if( object.second->getLS(cellId) < value){
            value = object.second->getLS(cellId) ;
            objectId = object.first ;
        }
    }

    return objectId ;
};

/*!
 * Get the object and part id of the projection point
 * If levelSetDefaults::OBJECT is passed as objectId, the information of the closest object are returned, otherwise those of the indicated object.
 * @param[in] cellId index of cell
 * @param[in] objectId index of object
 * @return pair containing object and part id
 */
std::pair<int,int> LevelSet::getClosestPart(const long &cellId, int objectId) const {

    if( objectId == levelSetDefaults::OBJECT){
        objectId = getClosestObject(cellId) ;
    }
    long partId = levelSetDefaults::PART ;

    if(objectId!=levelSetDefaults::OBJECT){
        partId = m_object.at(objectId)->getPart(cellId) ;
    }

    return ( std::make_pair(objectId, partId) );

};

/*!
 * Get the object and support id of the projection point
 * If levelSetDefaults::OBJECT is passed as objectId, the information of the closest object are returned, otherwise those of the indicated object.
 * @param[in] cellId index of cell
 * @param[in] objectId index of object
 * @return pair containing object and support id
 */
std::pair<int,long> LevelSet::getClosestSupport(const long &cellId, int objectId) const {

    if( objectId == levelSetDefaults::OBJECT){
        objectId = getClosestObject(cellId) ;
    }
    long supportId = levelSetDefaults::SUPPORT ;

    if(objectId!=levelSetDefaults::OBJECT){
        supportId = m_object.at(objectId)->getSupport(cellId) ;
    }

    return ( std::make_pair(objectId, supportId) );
};

/*!
 * Get the number of items which fall in the support radius of the cell
 * If levelSetDefaults::OBJECT is passed as objectId, the information of the closest object are returned, otherwise those of the indicated object.
 * @param[in] cellId index of cell
 * @param[in] objectId index of object
 * @return total number of items in narrow band
 */
int LevelSet::getSupportCount(const long &cellId, int objectId) const {


    if( objectId == levelSetDefaults::OBJECT){
        objectId = getClosestObject(cellId) ;
    }

    if( objectId != levelSetDefaults::OBJECT){
        return m_object.at(objectId)->getSupportCount(cellId) ;

    } else {
        return 0 ;
    }

};

/*!
 * Get the sign of the levelset function
 * If levelSetDefaults::OBJECT is passed as objectId, the sign of the composed levelset is returned, otherwise those of the indicated object.
 * @param[in] cellId index of cell
 * @param[in] objectId index of object
 * @return sign
 */
short LevelSet::getSign(const long &cellId, int objectId)const{

    if( objectId == levelSetDefaults::OBJECT){
        return( sign( getLS(cellId) ) );
    } else {
        return( sign(m_object.at(objectId)->getSign(cellId)) )  ;

    }

};

/*!
 * Get the sign of the levelset function
 * @param[in] cellId index of cell
 * @param[in] objectId index of object
 * @return sign
 */
bool LevelSet::isInNarrowBand(const long &cellId, int objectId)const{

    return( m_object.at(objectId)->isInNarrowBand(cellId) )  ;

};
/*!
 * Get the current size of the narrow band.
 * @param[in] objectId index of object
 * @return Physical size of the current narrow band.
 */
double LevelSet::getSizeNarrowBand( const int &objectId)const{
    return( m_object.at(objectId)->getSizeNarrowBand() ) ;
};

/*!
 * Set if the signed or unsigned levelset function should be computed.
 * @param[in] flag true/false for signed /unsigned Level-Set function .
 */
void LevelSet::setSign(bool flag){
    m_signedDF = flag;

};

/*!
 * Set if the levelset sign has to be propagated from the narrow band to the whole domain.
 * @param[in] flag True/false to active/disable the propagation .
 */
void LevelSet::setPropagateSign(bool flag){
    m_propagateS = flag;
};

/*!
 * Manually set the physical size of the narrow band.
 * @param[in] r Size of the narrow band.
 */
void LevelSet::setSizeNarrowBand(double r){
    m_userRSearch = true ;
    for( auto &object :m_object){
        object.second->setSizeNarrowBand(r) ;
    }
};

/*!
 * Computes levelset on given mesh with respect to the objects.
 * This routines needs to be called at least once.
 */
void LevelSet::compute(){

    double RSearch ;

    for( int objectId : m_order){
        auto &visitor = *(m_object.at(objectId)) ;
        if( !m_userRSearch){
            RSearch = m_kernel->computeSizeNarrowBand(&visitor)  ;
            visitor.setSizeNarrowBand(RSearch) ;
        }

        visitor.computeLSInNarrowBand( m_kernel.get(), RSearch, m_signedDF) ;
        if( m_propagateS ) visitor.propagateSign( m_kernel.get() ) ;
    }


    return ;
}

/*!
 * Updates the levelset after mesh adaptation.
 * @param[in] mapper mapper conatining mesh modifications
 */
void LevelSet::update( const std::vector<adaption::Info> &mapper ){

    // Udate the cache of the kernel
    m_kernel->updateGeometryCache( mapper ) ;

    // Check the mapper to detect the operations to perform
    bool updateNarrowBand = false;
    std::unordered_map<int,std::vector<long>> sendList, recvList ;

    for( const auto &event : mapper){
        if( event.entity != adaption::Entity::ENTITY_CELL){
            continue;
        }

        if( event.type == adaption::Type::TYPE_PARTITION_SEND){
            sendList.insert({{event.rank,event.previous}}) ;
        } else if( event.type == adaption::Type::TYPE_PARTITION_RECV){
            recvList.insert({{event.rank,event.current}}) ;
        } else if( event.type != adaption::Type::TYPE_PARTITION_NOTICE){
            updateNarrowBand = true ;
        }
    }

    // Evaluate new narrow band size
    double newRSearch ;
    for( int objectId : m_order){
        auto &visitor = *(m_object.at(objectId)) ;

        if (updateNarrowBand) {

            if(m_userRSearch){
                newRSearch = visitor.getSizeNarrowBand() ;
            } else {
                newRSearch = m_kernel->updateSizeNarrowBand( mapper, &visitor )  ;
            };

            visitor.updateLSInNarrowBand( m_kernel.get(), mapper, newRSearch, m_signedDF ) ;

        } else {
            visitor.clearAfterMeshAdaption(mapper) ;
        }


#if BITPIT_ENABLE_MPI
        // Parallel communications
        if (sendList.size() > 0 || recvList.size() > 0) {
            visitor.communicate( m_kernel.get(), sendList, recvList, &mapper ) ;
        }

        visitor.exchangeGhosts( m_kernel.get() ) ;
#endif

    // Finish narrow band update
        if (updateNarrowBand) {
            newRSearch = m_kernel->computeSizeNarrowBandFromLS( &visitor, m_signedDF );
            visitor.filterOutsideNarrowBand(newRSearch) ;
            visitor.setSizeNarrowBand(newRSearch) ;

            if( m_propagateS ) visitor.propagateSign( m_kernel.get() ) ;
        }

    }

    return;

};

/*! 
 * Writes LevelSetKernel to stream in binary format
 * @param[in] stream output stream
 */
void LevelSet::dump( std::ostream &stream ){

    IO::binary::write(stream, m_order);
    IO::binary::write(stream, m_userRSearch);
    IO::binary::write(stream, m_signedDF);
    IO::binary::write(stream, m_propagateS);

    for( const auto &object : m_object ){
        object.second->dump( stream ) ;
    }


    return ;
};

/*! 
 * Reads LevelSetKernel from stream in binary format
 * @param[in] stream output stream
 */
void LevelSet::restore( std::istream &stream ){

    IO::binary::read(stream, m_order);
    IO::binary::read(stream, m_userRSearch);
    IO::binary::read(stream, m_signedDF);
    IO::binary::read(stream, m_propagateS);

    for( const auto &object : m_object ){
        object.second->restore( stream ) ;
    }


    return ;
}

#if BITPIT_ENABLE_MPI
/*!
 * Checks if MPI communicator is available in underlying mesh.
 * If available, MPI communicator is retreived from mesh (and duplicated if necessary) and parallel processing can be done.
 * If not serial processing is necessary
 * @return true if parallel
 */
bool LevelSet::assureMPI( ){

    return(m_kernel->assureMPI() ) ;
}

#endif 

}

