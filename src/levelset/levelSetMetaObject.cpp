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

# include "bitpit_common.hpp"

# include "levelSetObject.hpp"
# include "levelSetMetaObject.hpp"

namespace bitpit {

/*!
	@interface LevelSetMetaObject
	@ingroup levelset
	@brief Interface class for all objects, which depend on other LevelSetObjects
*/

/*!
 * Constructor
 */
LevelSetMetaObject::LevelSetMetaObject(int id) : LevelSetObject(id){
}

/*!
 * Checks if the object is a primary object (e.g. of a surface triangulation)
 * or not (e.g. derived by boolean operations between two levelsets).
 *
 * Meta-objects are non-primary objects by definition.
 *
 * @return Returns always true, because meta-objects are non-primary objects
 * by definition.
 */
bool LevelSetMetaObject::isPrimary() const{
    return false;
} 

/*!
 * Get the the primary object that defines the levelset information for the
 * specified cell.
 * @param[in] id cell index
 * @return The the primary object that defines the levelset information for
 * the specified cell.
 */
const LevelSetObject * LevelSetMetaObject::getReferencePrimaryObject(long id) const{

    const LevelSetObject *referenceObject = getReferenceObject(id);
    if (!referenceObject) {
        return nullptr;
    }

    if (referenceObject->isPrimary()) {
        return referenceObject;
    }

    if( const LevelSetMetaObject *referenceMetaObject = dynamic_cast<const LevelSetMetaObject*>(referenceObject) ){
        return referenceMetaObject->getReferencePrimaryObject(id);
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
int LevelSetMetaObject::getReferenceObjectId(long id) const{

    const LevelSetObject *referenceObject = getReferenceObject(id);
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
int LevelSetMetaObject::getReferencePrimaryObjectId(long id) const{

    const LevelSetObject *referencePrimaryObject = getReferenceObject(id);
    if (!referencePrimaryObject) {
        return levelSetDefaults::OBJECT;
    }

    return referencePrimaryObject->getId();

}

/*!
 * Get the id of the object that defines the levelset information for the
 * specified cell.
 * @param[in] id cell index
 * @return The id of the object that defines the levelset information for the
 * specified cell.
 */
int LevelSetMetaObject::getPrimaryObjectId(long id) const{

    return getReferencePrimaryObjectId(id);

}

/*!
 * Get all primary objects that compose the meta object.
 * \return pointers to all primary objects involved in the definition of the
 * meta object levelset information.
 */
std::vector<const LevelSetObject*> LevelSetMetaObject::getPrimarySourceObjects() const{

    std::vector<const LevelSetObject*> objects;
    for( const LevelSetObject *sourceObject : getSourceObjects()){
        if( const LevelSetMetaObject *metaSourceObject = dynamic_cast<const LevelSetMetaObject*>(sourceObject) ){
            std::vector<const LevelSetObject*> sourcePrimarySourceObjects = metaSourceObject->getPrimarySourceObjects();
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
std::vector<int> LevelSetMetaObject::getSourceObjectIds() const{

    std::vector<const LevelSetObject*> sourceObjects = getSourceObjects();
    std::size_t nSourceObjects = sourceObjects.size();

    std::vector<int> sourceIds(nSourceObjects);
    for (std::size_t i = 0; i < nSourceObjects; ++i) {
        const LevelSetObject *sourceObject = sourceObjects[i];
        sourceIds[i] = sourceObject->getId();
    }

    return sourceIds;

}

/*!
 * Get all primary objects that compose the meta object.
 * \return identifiers of all primary objects involved in the definition of the
 * meta object levelset information.
 */
std::vector<int> LevelSetMetaObject::getPrimarySourceObjectIds() const{

    std::vector<const LevelSetObject*> primaryObjects = getPrimarySourceObjects();
    std::size_t nPrimaryObjects = primaryObjects.size();

    std::vector<int> primaryIds(nPrimaryObjects);
    for (std::size_t i = 0; i < nPrimaryObjects; ++i) {
        const LevelSetObject *primaryObject = primaryObjects[i];
        primaryIds[i] = primaryObject->getId();
    }

    return primaryIds;

}

}
