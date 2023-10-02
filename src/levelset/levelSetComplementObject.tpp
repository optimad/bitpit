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

# include <cassert>

# include "bitpit_IO.hpp"
# include "bitpit_common.hpp"

# include "levelSetObject.hpp"
# include "levelSetProxyObject.hpp"
# include "levelSetComplementObject.hpp"

namespace bitpit {

/*!
 * \class LevelSetComplementBaseObject
 * \ingroup levelset
 *
 * \brief Base class that allows to evaluate the complement of a LevelSetObjects.
 */

/*!
 * Constructor.
 *
 * \param[in] id identifier of object
 * \param[in] source pointer to source object
 */
template<typename SourceLevelSetObject>
LevelSetComplementBaseObject<SourceLevelSetObject>::LevelSetComplementBaseObject(int id, const SourceLevelSetObject *source)
    : LevelSetProxyObject<SourceLevelSetObject, SourceLevelSetObject>(id),
      m_sourceObject(source)
{
}

/*!
 * Replace a source object.
 *
 * \param[in] current current source object
 * \param[in] updated updated source object
 */
template<typename SourceLevelSetObject>
void LevelSetComplementBaseObject<SourceLevelSetObject>::replaceSourceObject(const SourceLevelSetObject *current, const SourceLevelSetObject *updated)
{
    if (current != m_sourceObject) {
        throw std::runtime_error("Unable to find the source that should be replaced.");
    }

    m_sourceObject = updated;
}

/*!
 * Get the levelset value.
 *
 * \param[in] id cell id
 * \return levelset value in cell
 */
template<typename SourceLevelSetObject>
double LevelSetComplementBaseObject<SourceLevelSetObject>::getValue(long id) const
{
    return (- m_sourceObject->getValue(id));
}

/*!
 * Get the levelset gradient.
 *
 * \param[in] id cell id
 * \return levelset gradient in cell
 */
template<typename SourceLevelSetObject>
std::array<double,3> LevelSetComplementBaseObject<SourceLevelSetObject>::getGradient(long id) const
{
    return (- 1. * m_sourceObject->getGradient(id));
}

/*!
 * Computes the LevelSetInfo in a point.
 *
 * \param[in] coords point coordinates
 * \return LevelSetInfo
*/
template<typename SourceLevelSetObject>
LevelSetInfo LevelSetComplementBaseObject<SourceLevelSetObject>::computeLevelSetInfo(const std::array<double,3> &coords) const
{
    LevelSetInfo levelSetInfo = m_sourceObject->computeLevelSetInfo(coords);
    levelSetInfo.value    *= -1.;
    levelSetInfo.gradient *= -1.;

    return levelSetInfo;
}

/*!
 * Get the object that defines the levelset information for the specified cell.
 *
 * \param[in] id cell index
 * \return The object that defines the levelset information for the specified
 * cell.
 */
template<typename SourceLevelSetObject>
const SourceLevelSetObject * LevelSetComplementBaseObject<SourceLevelSetObject>::getReferenceObject(long id) const
{
    BITPIT_UNUSED(id);

    return m_sourceObject;
}

/*!
 * Get the objects that define
 *
 * \return pointers to all primary objects involved in the definition of the
 * boolean object levelset information.
 */
template<typename SourceLevelSetObject>
const SourceLevelSetObject * LevelSetComplementBaseObject<SourceLevelSetObject>::getSourceObject() const
{
    return m_sourceObject;
}

/*!
 * Get all objects that compose the boolean object.
 *
 * \return pointers to all primary objects involved in the definition of the
 * boolean object levelset information.
 */
template<typename SourceLevelSetObject>
std::vector<const SourceLevelSetObject *> LevelSetComplementBaseObject<SourceLevelSetObject>::getSourceObjects() const
{
    return std::vector<const SourceLevelSetObject *>(1, m_sourceObject);
}

}