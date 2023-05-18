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

#include "levelSetCache.hpp"

namespace bitpit {

/*!
 * Dummy vector to be used when the cache entry stores a real boolean reference
 */
const std::vector<bool> LevelSetValueCacheEntry<bool>::m_dummyVector = std::vector<bool>(1, false);

/*!
 * Constructor for invalid entries.
 */
LevelSetValueCacheEntry<bool>::LevelSetValueCacheEntry()
    : LevelSetValueCacheBaseEntry<bool>(false),
      m_value(LevelSetValueCacheBaseEntry<bool>::m_dummyValue),
      m_useVectorValue(false), m_vectorValue(m_dummyVector[0])
{
}

/*!
 * Constructor for valid entries.
 *
 * \param value is the cached value
 */
LevelSetValueCacheEntry<bool>::LevelSetValueCacheEntry(const bool &value)
    : LevelSetValueCacheBaseEntry<bool>(true),
      m_value(value),
      m_useVectorValue(false), m_vectorValue(m_dummyVector[0])
{
}

/*!
 * Constructor for valid entries.
 *
 * \param value is the cached value
 */
LevelSetValueCacheEntry<bool>::LevelSetValueCacheEntry(const std::vector<bool>::reference &value)
    : LevelSetValueCacheBaseEntry<bool>(true),
      m_value(LevelSetValueCacheBaseEntry<bool>::m_dummyValue),
      m_useVectorValue(true), m_vectorValue(value)
{
}

/*!
 * Get the stored value.
 *
 * Trying to deference an invalid entry results in undefined behaviour.
 */
bool LevelSetValueCacheEntry<bool>::operator*() const
{
    if (m_useVectorValue) {
        return m_vectorValue;
    } else {
        return m_value;
    }
}

}
