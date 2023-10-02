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

# ifndef __BITPIT_LEVELSET_TPP__
# define __BITPIT_LEVELSET_TPP__

namespace bitpit{

/*!
 * Adds the complement of the specified object.
 * Objects can be added to the levelset only after setting the mesh.
 * @param[in] sourceId id of source object
 * @param[in] id id to be assigned to object. In case default value is passed the insertion order will be used as identifier
 * @return identifier of new object
 */
template<typename LevelSetSourceObject>
int LevelSet::addObjectComplement( int sourceId, int id ) {

    const LevelSetSourceObject *sourceObject = getObjectPtr<LevelSetSourceObject>(sourceId) ;
    if (!sourceObject) {
        throw std::runtime_error("The type of the object does not match the type of the complement object!");
    }

    auto object = std::unique_ptr<LevelSetObject>(new LevelSetComplementObject<LevelSetSourceObject>(id, sourceObject));

    return registerObject(std::move(object));
};

/*!
 * Adds a boolean operation between two objects
 * Objects can be added to the levelset only after setting the mesh.
 * @param[in] operation boolean operation
 * @param[in] id1 id of first operand
 * @param[in] id2 id of second operand
 * @param[in] id id to be assigned to object. In case default value is passed the insertion order will be used as identifier
 * @return identifier of new object
 */
template<typename LevelSetSourceObject>
int LevelSet::addObject( LevelSetBooleanOperation operation, int id1, int id2, int id ) {

    const LevelSetObject *object1 = getObjectPtr<LevelSetSourceObject>(id1) ;
    if (!object1) {
        throw std::runtime_error("The type of the object does not match the type of the boolean object!");
    }

    const LevelSetObject *object2 = getObjectPtr<LevelSetSourceObject>(id2) ;
    if (!object2) {
        throw std::runtime_error("The type of the object does not match the type of the boolean object!");
    }

    auto object = std::unique_ptr<LevelSetObject>(new LevelSetBooleanObject<LevelSetSourceObject>(id, operation, object1, object2));

    return registerObject(std::move(object));
}

/*!
 * Adds a boolean operation between that will be applied recursively to a series of objects
 * Objects can be added to the levelset only after setting the mesh.
 * @param[in] operation boolean operation
 * @param[in] ids vector with indices of operand objects
 * @param[in] id id to be assigned to object. In case default value is passed the insertion order will be used as identifier
 * @return identifier of new object
 */
template<typename LevelSetSourceObject>
int LevelSet::addObject( LevelSetBooleanOperation operation, const std::vector<int> &ids, int id ) {

    std::vector<const LevelSetSourceObject*> objects;
    for( int id : ids){
        const LevelSetSourceObject *object = getObjectPtr<LevelSetSourceObject>(id) ;
        if (!object) {
            throw std::runtime_error("The type of the object does not match the type of the boolean object!");
        }

        objects.push_back( object );
    }

    auto object = std::unique_ptr<LevelSetObject>(new LevelSetBooleanObject<LevelSetSourceObject>(id, operation, objects));

    return registerObject(std::move(object));
}

/*!
 * Get a constant reference to the specified object.
 * If the specified id does not exist an exception is thrown.
 * @param id is the object id
 * @return reference to levelset object
 */
template<typename T>
T & LevelSet::getObject( int id) const{
    return dynamic_cast<T &>(*m_objects.at(id)) ;
}

/*!
 * Get a constant pointer to the specified object.
 * If the specified id does not exist an exception is thrown.
 * @param id is the object id
 * @return pointer to levelset object
 */
template<typename T>
T * LevelSet::getObjectPtr( int id) const{
    return dynamic_cast<T *>(m_objects.at(id).get()) ;
}

/*!
 * Get vector of pointers to all object.
 * @return vector of pointers to levelset objects
 */
template<typename T>
std::vector<T *>  LevelSet::getObjectPtrs( ) const{
    std::vector<T *> objects;
    objects.reserve(m_objects.size());
    for(const auto &entry : m_objects){
        T *object = getObjectPtr(entry.first) ;
        if(object){
            objects.push_back( object ) ;
        }
    }

    return objects;
}

}

#endif
