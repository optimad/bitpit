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

#include "levelSetBooleanObject.hpp"
#include "levelSetBooleanObject.tpp"
#include "levelSetObject.hpp"

namespace bitpit {

/*!
	@class      LevelSetBooleanObject
	@ingroup    levelset
	@brief      Class which deals with boolean operation between two LevelSetObjects
*/

/*!
 * Constructor taking two objects.
 * @param[in] id identifier of object
 * @param[in] op type of boolean operation
 * @param[in] source1 pointer to first source object
 * @param[in] source2 pointer to second source object
 */
LevelSetBooleanObject<LevelSetObject>::LevelSetBooleanObject( int id, LevelSetBooleanOperation op, const LevelSetObject *source1, const LevelSetObject *source2  )
    : LevelSetBooleanBaseObject<LevelSetObject>(id, op, source1, source2) {
}

/*!
 * Constructor taking a vector of objects.
 * The boolean operation will be applied recursively on each entry.
 * @param[in] id identifier of object
 * @param[in] op type of boolean operation
 * @param[in] sourceObjects pointers to source objects
 */
LevelSetBooleanObject<LevelSetObject>::LevelSetBooleanObject( int id, LevelSetBooleanOperation op, const std::vector<const LevelSetObject *> &sourceObjects )
    : LevelSetBooleanBaseObject<LevelSetObject>(id, op, sourceObjects) {
}

/*!
 * Clones the object
 * @return pointer to cloned object
 */
LevelSetBooleanObject<LevelSetObject> * LevelSetBooleanObject<LevelSetObject>::clone() const {
    return new LevelSetBooleanObject<LevelSetObject>( *this );
}

}
