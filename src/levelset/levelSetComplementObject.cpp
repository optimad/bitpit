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

#include "levelSetComplementObject.hpp"
#include "levelSetComplementObject.tpp"

namespace bitpit{

/*!
 * \class LevelSetComplementObject
 * \ingroup levelset
 *
 * \brief Class that allows to evaluate the complement of a LevelSetObjects.
 */

/*!
 * Constructor.
 *
 * \param[in] id identifier of object
 * \param[in] source pointer to source object
 */
LevelSetComplementObject<LevelSetObject>::LevelSetComplementObject(int id, const LevelSetObject *source)
    : LevelSetComplementBaseObject<LevelSetObject>(id, source)
{
}

/*!
 * Clones the object
 * @return pointer to cloned object
 */
LevelSetComplementObject<LevelSetObject> * LevelSetComplementObject<LevelSetObject>::clone() const
{
    return new LevelSetComplementObject<LevelSetObject>(*this);
}

}
