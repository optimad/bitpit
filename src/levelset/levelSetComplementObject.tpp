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
 * Checks if the object is empty.
 *
 * \result Returns true is the object is empty, false otherwise.
 */
template<typename SourceLevelSetObject>
bool LevelSetComplementBaseObject<SourceLevelSetObject>::empty() const
{
    return m_sourceObject->empty();
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
 * Fill the cache that contains the propagated cell sign.
 */
template<typename SourceLevelSetObject>
void LevelSetComplementBaseObject<SourceLevelSetObject>::fillCellPropagatedSignCache()
{
    // Early return if propagated sign cannot be copied from the source object
    LevelSetBulkEvaluationMode sourceBulkEvaluationMode = m_sourceObject->getCellBulkEvaluationMode();

    bool useSourceSign = false;
    if (sourceBulkEvaluationMode == LevelSetBulkEvaluationMode::SIGN_PROPAGATION) {
        useSourceSign = true;
    } else if (sourceBulkEvaluationMode == LevelSetBulkEvaluationMode::EXACT) {
        LevelSetCacheMode sourceSignCacheMode = m_sourceObject->getFieldCellCacheMode(LevelSetField::SIGN);
        if (sourceSignCacheMode == LevelSetCacheMode::FULL) {
            useSourceSign = true;
        } else {
            LevelSetCacheMode sourceValueCacheMode = m_sourceObject->getFieldCellCacheMode(LevelSetField::VALUE);
            if (sourceValueCacheMode == LevelSetCacheMode::FULL) {
                useSourceSign = true;
            }
        }
    }

    if (!useSourceSign) {
        SourceLevelSetObject::fillCellPropagatedSignCache();
        return;
    }

    // Mesh information
    const VolumeKernel &mesh = *(this->getKernel()->getMesh()) ;
    VolumeKernel::CellConstIterator cellBegin = mesh.cellConstBegin();
    VolumeKernel::CellConstIterator cellEnd   = mesh.cellConstEnd();

    // Get cache for sign propagation
    typedef typename SourceLevelSetObject::CellCacheCollection::template ValueCache<char> ZoneCache;
    ZoneCache *propagatedSignCache = this->template getCellCache<char>(this->m_cellPropagatedSignCacheId);

    // Get propagated sign from source object
    for (VolumeKernel::CellConstIterator cellItr = cellBegin; cellItr != cellEnd; ++cellItr) {
        long cellId = cellItr.getId();
        short cellSign = - m_sourceObject->evalCellSign(cellId);

        propagatedSignCache->insertEntry(cellId, static_cast<char>(cellSign));
    }
}

/*!
 * Evaluate levelset sign at the specified cell.
 *
 * \param id is the id of the cell
 * \result The sign of the levelset at the specified cell.
 */
template<typename SourceLevelSetObject>
short LevelSetComplementBaseObject<SourceLevelSetObject>::_evalCellSign(long id) const
{
    return (-1 * m_sourceObject->evalCellSign(id));
}

/*!
 * Evaluate levelset value at the specified cell.
 *
 * \param id is the id of the cell
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The value of the levelset at the specified cell.
 */
template<typename SourceLevelSetObject>
double LevelSetComplementBaseObject<SourceLevelSetObject>::_evalCellValue(long id, bool signedLevelSet) const
{
    double value = m_sourceObject->evalCellValue(id, signedLevelSet);
    if (signedLevelSet) {
        value *= -1.;
    }

    return value;
}

/*!
 * Evaluate levelset gradient at the specified cell.
 *
 * \param id is the id of the cell
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The gradient of the levelset at the specified cell.
 */
template<typename SourceLevelSetObject>
std::array<double,3> LevelSetComplementBaseObject<SourceLevelSetObject>::_evalCellGradient(long id, bool signedLevelSet) const
{
    std::array<double,3> gradient = m_sourceObject->evalCellGradient(id, signedLevelSet);
    if (signedLevelSet) {
        gradient *= -1.;
    }

    return gradient;
}

/*!
 * Evaluate levelset sign at the specified point.
 *
 * \param point are the coordinates of the point
 * \result The sign of the levelset at the specified point.
 */
template<typename SourceLevelSetObject>
short LevelSetComplementBaseObject<SourceLevelSetObject>::_evalSign(const std::array<double,3> &point) const
{
    return (-1 * m_sourceObject->evalSign(point));
}

/*!
 * Evaluate levelset value at the specified point.
 *
 * \param point are the coordinates of the point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The value of the levelset at the specified point.
 */
template<typename SourceLevelSetObject>
double LevelSetComplementBaseObject<SourceLevelSetObject>::_evalValue(const std::array<double,3> &point, bool signedLevelSet) const
{
    double value = m_sourceObject->evalValue(point, signedLevelSet);
    if (signedLevelSet) {
        value *= -1.;
    }

    return value;
}

/*!
 * Evaluate levelset gradient at the specified point.
 *
 * \param point are the coordinates of the point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The gradient of the levelset at the specified point.
 */
template<typename SourceLevelSetObject>
std::array<double,3> LevelSetComplementBaseObject<SourceLevelSetObject>::_evalGradient(const std::array<double,3> &point, bool signedLevelSet) const
{
    std::array<double,3> gradient = m_sourceObject->evalGradient(point, signedLevelSet);
    if (signedLevelSet) {
        gradient *= -1.;
    }

    return gradient;
}

/*!
 * Get the object that defines the levelset information for the specified cell.
 *
 * \param[in] id cell index
 * \return The object that defines the levelset information for the specified
 * cell.
 */
template<typename SourceLevelSetObject>
const SourceLevelSetObject * LevelSetComplementBaseObject<SourceLevelSetObject>::getCellReferenceObject(long id) const
{
    BITPIT_UNUSED(id);

    return m_sourceObject;
}

/*!
 * Get the object that defines the levelset information for the specified point.
 *
 * @param[in] point are the coordinates of the point
 * \return The object that defines the levelset information for the specified
 * point.
 */
template<typename SourceLevelSetObject>
const SourceLevelSetObject * LevelSetComplementBaseObject<SourceLevelSetObject>::getReferenceObject(const std::array<double, 3> &point) const
{
    BITPIT_UNUSED(point);

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
