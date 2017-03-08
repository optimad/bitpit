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

# ifndef __BITPIT_LEVELSET_TPP__
# define __BITPIT_LEVELSET_TPP__

namespace bitpit{

/*!
 * Get a constant reference to the specified object.
 * If the specified id does not exist an exception is thrown.
 * @param id is the object id
 * @return reference to levelset object
 */
template<typename T>
const T & LevelSet::getObject( int id) const{
    return dynamic_cast<const T &>(*m_objects.at(id)) ;
}

/*!
 * Get a constant pointer to the specified object.
 * If the specified id does not exist an exception is thrown.
 * @param id is the object id
 * @return pointer to levelset object
 */
template<typename T>
const T * LevelSet::getObjectPtr( int id) const{
    return dynamic_cast<const T *>(m_objects.at(id).get()) ;
}

/*!
 * Get vector of pointers to all object.
 * @return vector of pointers to levelset objects
 */
template<typename T>
std::vector<T const *>  LevelSet::getObjectPtrs( ) const{
    std::vector<T const *> objects;
    objects.reserve(m_objects.size());
    for( auto const &entry : m_objects){
        T const *object = getObjectPtr(entry.first) ;
        if(object){
            objects.push_back( object ) ;
        }
    }

    return objects;
}

}

#endif
