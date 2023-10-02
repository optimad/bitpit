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
 * If cell centroid lies within the narrow band and hence levelset is computet exactly
 * @param[in] id cell id
 * @return true/false if the centroid is in narrow band
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
bool LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::isInNarrowBand(long id)const{

    const SourceLevelSetObject *referenceObject = getReferenceObject(id) ;
    if ( referenceObject ) {
        return referenceObject->isInNarrowBand(id);
    }

    return false;

}

/*!
 * Get the the primary object that defines the levelset information for the
 * specified cell.
 * @param[in] id cell index
 * @return The the primary object that defines the levelset information for
 * the specified cell.
 */
template<typename SourceLevelSetObject, typename BaseLevelSetObject>
const SourceLevelSetObject * LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::getReferencePrimaryObject(long id) const{

    const SourceLevelSetObject *referenceObject = getReferenceObject(id);
    if (!referenceObject) {
        return nullptr;
    }

    if (referenceObject->isPrimary()) {
        return referenceObject;
    }

    if( const LevelSetProxyObject *referenceProxyObject = dynamic_cast<const LevelSetProxyObject*>(referenceObject) ){
        return referenceProxyObject->getReferencePrimaryObject(id);
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
int LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::getReferenceObjectId(long id) const{

    const SourceLevelSetObject *referenceObject = getReferenceObject(id);
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
int LevelSetProxyObject<SourceLevelSetObject, BaseLevelSetObject>::getReferencePrimaryObjectId(long id) const{

    const SourceLevelSetObject *referenceObject = getReferenceObject(id);
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

    return referenceProxyObject->getReferencePrimaryObjectId(id);
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

}

#endif
