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
LevelSetComplementObject::LevelSetComplementObject(int id, const LevelSetObject *source)
    : LevelSetProxyObject(id),
      m_sourceObject(source)
{
}

/*!
 * Copy constructor.
 *
 * Assigns same id to new object;
 * \param[in] other object to be copied
 */
LevelSetComplementObject::LevelSetComplementObject(const LevelSetComplementObject &other)
    : LevelSetProxyObject(other),
      m_sourceObject(other.m_sourceObject)
{
}

/*!
 * Get the levelset value.
 *
 * \param[in] id cell id
 * \return levelset value in cell
 */
double LevelSetComplementObject::getValue(long id) const
{
    return (- m_sourceObject->getValue(id));
}

/*!
 * Get the levelset gradient.
 *
 * \param[in] id cell id
 * \return levelset gradient in cell
 */
std::array<double,3> LevelSetComplementObject::getGradient(long id) const
{
    return (- 1. * m_sourceObject->getGradient(id));
}

/*!
 * Computes the LevelSetInfo in a point.
 *
 * \param[in] coords point coordinates
 * \return LevelSetInfo
*/
LevelSetInfo LevelSetComplementObject::computeLevelSetInfo(const std::array<double,3> &coords) const
{
    LevelSetInfo levelSetInfo = m_sourceObject->computeLevelSetInfo(coords);
    levelSetInfo.value    *= -1.;
    levelSetInfo.gradient *= -1.;

    return levelSetInfo;
}

/*!
 * Replace a source object.
 *
 * \param[in] current current source object
 * \param[in] updated updated source object
 */
void LevelSetComplementObject::replaceSourceObject(const LevelSetObject *current, const LevelSetObject *updated)
{
    if (current != m_sourceObject) {
        throw std::runtime_error("Unable to find the source that should be replaced.");
    }

    m_sourceObject = updated;
}

/*!
 * Clones the object.
 *
 * \return pointer to cloned object
 */
LevelSetComplementObject* LevelSetComplementObject::clone() const
{
    return new LevelSetComplementObject(*this);
}

/*!
 * Gets the surface normal at the projection point.
 *
 * \param[in] id cell index
 * \return closest part
 */
std::array<double,3> LevelSetComplementObject::getNormal(long id) const
{
    return (-1. * m_sourceObject->getNormal(id));
}

/*!
 * Gets the closest part index.
 *
 * \param[in] id cell index
 * \return closest part
 */
int LevelSetComplementObject::getPart(long id) const
{
    return m_sourceObject->getPart(id);
}

/*!
 * Get surface feature size.
 *
 * \param[in] id cell index
 * \return characteristic size
 */
double LevelSetComplementObject::getSurfaceFeatureSize(long id) const
{
    return m_sourceObject->getSurfaceFeatureSize(id);
}

/*!
 * Get the smallest surface feature size.
 *
 * \return characteristic size
 */
double LevelSetComplementObject::getMinSurfaceFeatureSize() const
{
    return m_sourceObject->getMinSurfaceFeatureSize();
}

/*!
 * Get the largest surface feature size.
 *
 * \return characteristic size
 */
double LevelSetComplementObject::getMaxSurfaceFeatureSize() const
{
    return m_sourceObject->getMaxSurfaceFeatureSize();
}

/*!
 * Get the object that defines the levelset information for the specified cell.
 *
 * \param[in] id cell index
 * \return The object that defines the levelset information for the specified
 * cell.
 */
const LevelSetObject * LevelSetComplementObject::getReferenceObject(long id) const
{
    BITPIT_UNUSED(id);

    return m_sourceObject;
}

/*!
 * Get all objects that compose the boolean object.
 *
 * \return pointers to all primary objects involved in the definition of the
 * boolean object levelset information.
 */
std::vector<const LevelSetObject *> LevelSetComplementObject::getSourceObjects() const
{
    return std::vector<const LevelSetObject *>(1, m_sourceObject);
}

}
