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

# ifndef __BITPIT_LEVELSET_PROXY_OBJECT_TPP__
# define __BITPIT_LEVELSET_PROXY_OBJECT_TPP__

namespace bitpit {

/*!
	@interface LevelSetProxyObject
	@ingroup levelset
	@brief Interface class for all objects, which depend on other LevelSetObjects
*/

/*!
 * Constructor
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::LevelSetProxyObject(int id) : BaseLevelSetObject(id){
}

/*!
 * Checks if the object is a primary object (e.g. of a surface triangulation)
 * or not (e.g. derived by boolean operations between two levelsets).
 *
 * Proxy-objects are non-primary objects by definition.
 *
 * @return Returns always true, because proxy-objects are non-primary objects
 * by definition.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
bool LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::isPrimary() const{
    return false;
} 

/*!
 * Fill the cache that contains the zone associated to the cells.
 *
 * A cell can be either in the narrow band or in the bulk. It will be considered inside the narrow
 * band if one of the following conditions holds:
 *  - its distance from the surface is less than the narrow band size;
 *  - it intersects the zero-levelset iso-surface (intersections are checked using the
 *    FAST_GUARANTEE_FALSE criterium);
 *  - one of its neighbors intersects the zero-levelset iso-surface.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
void LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::fillCellLocationCache()
{
    // Mesh information
    const VolumeKernel &mesh = *(this->getKernel()->getMesh()) ;
    VolumeKernel::CellConstIterator cellBegin = mesh.cellConstBegin();
    VolumeKernel::CellConstIterator cellEnd   = mesh.cellConstEnd();

    // Get cell zone cache
    typedef typename BaseLevelSetObject::CellCacheCollection::template ValueCache<char> BaseZoneCache;
    BaseZoneCache *locationCache = this->template getCellCache<char>(this->m_cellLocationCacheId);

    // Get zone information from reference object
    for (VolumeKernel::CellConstIterator cellItr = cellBegin; cellItr != cellEnd; ++cellItr) {
        long cellId = cellItr.getId();
        LevelSetCellLocation cellLocation = getCellReferenceObject(cellId)->getCellLocation(cellId);
        locationCache->insertEntry(cellId, static_cast<char>(cellLocation));
    }
}

/*!
 * Fill the cache that contains the zone associated to the cells.
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
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
void LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::fillCellLocationCache(const std::vector<adaption::Info> &adaptionData)
{
    // Get cell location cache
    typedef typename BaseLevelSetObject::CellCacheCollection::template ValueCache<char> BaseZoneCache;
    BaseZoneCache *locationCache = this->template getCellCache<char>(this->m_cellLocationCacheId);

    // Get location information from reference object
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
            LevelSetCellLocation cellLocation = getCellReferenceObject(cellId)->getCellLocation(cellId);
            locationCache->insertEntry(cellId, static_cast<char>(cellLocation));
        }
    }
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
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
bool LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::isCellInNarrowBand(long id)const
{
    return getCellReferenceObject(id)->isCellInNarrowBand(id);
}

/*!
 * Check if the specified point lies inside the narrow band.
 *
 * The value of the levelset is evaluated and compared with the specified narrow band size.
 *
 * \param[in] id is the cell id
 * \result Return true if the cell is in the narrow band, false otherwise.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
bool LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::isInNarrowBand(const std::array<double,3> &point)const
{
    return getReferenceObject(point)->isInNarrowBand(point);
}

/*!
 * Get the the primary object that defines the levelset information for the
 * specified cell.
 * @param[in] id cell index
 * @return The the primary object that defines the levelset information for
 * the specified cell.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
const SourceLevelSetObject * LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::getCellReferencePrimaryObject(long id) const{

    const SourceLevelSetObject *referenceObject = getCellReferenceObject(id);
    if (!referenceObject) {
        return nullptr;
    }

    if (referenceObject->isPrimary()) {
        return referenceObject;
    }

    if( const LevelSetProxyObject *referenceProxyObject = dynamic_cast<const LevelSetProxyObject*>(referenceObject) ){
        return referenceProxyObject->getCellReferencePrimaryObject(id);
    }

    return nullptr;

}

/*!
 * Get the the primary object that defines the levelset information for the
 * specified point.
 * @param[in] point are the coordinates of the point
 * @return The the primary object that defines the levelset information for
 * the specified cell.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
const SourceLevelSetObject * LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::getReferencePrimaryObject(const std::array<double, 3> &point) const{

    const SourceLevelSetObject *referenceObject = getReferenceObject(point);
    if (!referenceObject) {
        return nullptr;
    }

    if (referenceObject->isPrimary()) {
        return referenceObject;
    }

    if( const LevelSetProxyObject *referenceProxyObject = dynamic_cast<const LevelSetProxyObject*>(referenceObject) ){
        return referenceProxyObject->getReferencePrimaryObject(point);
    }

    return nullptr;

}

/*!
 * Get the id of the object that defines the levelset information for the
 * specified cell.
 * @param[in] id cell index
 * @return The id of the object that defines the levelset information for the
 * specified cell.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
int LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::getCellReferenceObjectId(long id) const{

    const SourceLevelSetObject *referenceObject = getCellReferenceObject(id);
    if (!referenceObject) {
        return levelSetDefaults::OBJECT;
    }

    return referenceObject->getId();

}

/*!
 * Get the the primary object that defines the levelset information for the
 * specified cell.
 * @param[in] id cell index
 * @return The the primary object that defines the levelset information for
 * the specified cell.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
int LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::getCellReferencePrimaryObjectId(long id) const{

    const SourceLevelSetObject *referenceObject = getCellReferenceObject(id);
    if (!referenceObject) {
        return levelSetDefaults::OBJECT;
    }

    if (referenceObject->isPrimary()) {
        return referenceObject->getId();
    }

    const LevelSetProxyBaseObject *referenceProxyObject = dynamic_cast<const LevelSetProxyBaseObject *>(referenceObject);
    if (!referenceProxyObject) {
        return levelSetDefaults::OBJECT;
    }

    return referenceProxyObject->getCellReferencePrimaryObjectId(id);
}

/*!
 * Get the id of the object that defines the levelset information for the
 * specified point.
 * @param[in] point are the coordinates of the point
 * @return The id of the object that defines the levelset information for the
 * specified point.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
int LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::getReferenceObjectId(const std::array<double, 3> &point) const{

    const SourceLevelSetObject *referenceObject = getReferenceObject(point);
    if (!referenceObject) {
        return levelSetDefaults::OBJECT;
    }

    return referenceObject->getId();

}

/*!
 * Get the the primary object that defines the levelset information for the
 * specified point.
 * @param[in] point are the coordinates of the point
 * @return The the primary object that defines the levelset information for
 * the specified point.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
int LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::getReferencePrimaryObjectId(const std::array<double, 3> &point) const{

    const SourceLevelSetObject *referenceObject = getReferenceObject(point);
    if (!referenceObject) {
        return levelSetDefaults::OBJECT;
    }

    if (referenceObject->isPrimary()) {
        return referenceObject->getId();
    }


    const LevelSetProxyBaseObject *referenceProxyObject = dynamic_cast<const LevelSetProxyBaseObject *>(referenceObject);
    if (!referenceProxyObject) {
        return levelSetDefaults::OBJECT;
    }

    return referenceProxyObject->getReferencePrimaryObjectId(point);
}

/*!
 * Get all primary objects that compose the proxy object.
 * \return pointers to all primary objects involved in the definition of the
 * proxy object levelset information.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
std::vector<const SourceLevelSetObject *> LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::getPrimarySourceObjects() const{

    std::vector<const SourceLevelSetObject *> objects;
    for( const SourceLevelSetObject *sourceObject : getSourceObjects()){
        if( const LevelSetProxyObject *proxySourceObject = dynamic_cast<const LevelSetProxyObject*>(sourceObject) ){
            std::vector<const SourceLevelSetObject *> sourcePrimarySourceObjects = proxySourceObject->getPrimarySourceObjects();
            objects.insert(objects.end(), sourcePrimarySourceObjects.begin(), sourcePrimarySourceObjects.end());
        } else {
            objects.push_back(sourceObject);
        }
    }

    return objects;

}

/*!
 * Get all primary objects that compose the meta object.
 * \return identifiers of all primary objects involved in the definition of the
 * meta object levelset information.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
std::vector<int> LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::getSourceObjectIds() const{

    std::vector<const SourceLevelSetObject*> sourceObjects = getSourceObjects();
    std::size_t nSourceObjects = sourceObjects.size();

    std::vector<int> sourceIds(nSourceObjects);
    for (std::size_t i = 0; i < nSourceObjects; ++i) {
        const SourceLevelSetObject *sourceObject = sourceObjects[i];
        sourceIds[i] = sourceObject->getId();
    }

    return sourceIds;

}

/*!
 * Get all primary objects that compose the meta object.
 * \return identifiers of all primary objects involved in the definition of the
 * meta object levelset information.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
std::vector<int> LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::getPrimarySourceObjectIds() const{

    std::vector<const SourceLevelSetObject*> primaryObjects = getPrimarySourceObjects();
    std::size_t nPrimaryObjects = primaryObjects.size();

    std::vector<int> primaryIds(nPrimaryObjects);
    for (std::size_t i = 0; i < nPrimaryObjects; ++i) {
        const SourceLevelSetObject *primaryObject = primaryObjects[i];
        primaryIds[i] = primaryObject->getId();
    }

    return primaryIds;

}

/*!
 * Get the id of the object that defines the levelset information for the
 * specified cell.
 * @param[in] id cell index
 * @return The id of the object that defines the levelset information for the
 * specified cell.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
int LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::getPrimaryObjectId(long id) const{

    return getReferencePrimaryObjectId(id);

}

}

#endif
