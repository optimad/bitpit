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

#include <cassert>
#include <limits>

#include "bitpit_IO.hpp"

#include "index_generator.hpp"

namespace bitpit {

/*!
    \class IndexGenerator
    \ingroup GenericIO

    \brief The IndexGenerator class allows to generate unique ids.
*/

const long IndexGenerator::NULL_ID = std::numeric_limits<long>::min();

/*!
    Creates a new generator.
*/
IndexGenerator::IndexGenerator()
    : m_latest(NULL_ID), m_highest(NULL_ID)
{

}

/*!
    Generates a unique index.

    If the trash is empty a new index is generated, otherwise an index taken
    from the trash is recycled.

    \return A new unique index.
*/
long IndexGenerator::generate()
{
    // If the trash is empty generate a new id otherwise recycle the first id
    // in the trash.
    if (m_trash.empty()) {
        if (m_highest == NULL_ID) {
            m_highest = 0;
        } else {
            assert(m_highest < std::numeric_limits<long>::max());
            ++m_highest;
        }
        m_latest = m_highest;
    } else {
        m_latest = m_trash.front();
        m_trash.pop_front();
    }

    return m_latest;
}

/*!
    Gets the latest assigned id.

    \return The latest assigned index.
*/
long IndexGenerator::getLatest()
{
    return m_latest;
}

/*!
    Gets the highest assigned id.

    \return The highest assigned index.
*/
long IndexGenerator::getHighest()
{
    return m_highest;
}

/*!
    Checks if an id is currently assigned.

    \return True if the id has already been assigned, false otherwise.
*/
bool IndexGenerator::isAssigned(long id)
{
    // If the latest id is valid, it is assigned by definition
    if (m_latest >= 0 && id == m_latest) {
        return true;
    }

    // Ids past the highest one are not assigned
    if (id > m_highest) {
        return false;
    }

    // The id is assigned only if is not in the trash
    for (long trashedId : m_trash) {
        if (trashedId == id) {
            return false;
        }
    }

    return true;
}

/*!
    Marks the specified id as currently assigned.

    \param id is the id to be marked as assigned
*/
void IndexGenerator::setAssigned(long id)
{
    // If the id is past the highest assigned id we need to trash all the ids
    // from the highest assigned to this one (the generator only handles
    // contiguous ids), otherwise look for the id in the trash and, if found,
    // recycle it (if the id is not in the trash, this means it was already
    // an assigned id).
    if (id > m_highest) {
        for (long wasteId = m_highest + 1; wasteId < id; ++wasteId) {
            trash(wasteId);
        }

        m_highest = id;
    } else {
        for (auto trashItr = m_trash.begin(); trashItr != m_trash.end(); ++trashItr) {
            long trashedId = *trashItr;
            if (trashedId == id) {
                m_trash.erase(trashItr);
                break;
            }
        }
    }

    // This is the latest assigned id
    m_latest = id;
}

/*!
    Trashes an index.

    A trashed index is an index no more used that can be recycled.

    \param id is the index that will be trashed
*/
void IndexGenerator::trash(const long &id)
{
    m_trash.push_back(id);

    // We only keep track of the latest assigned id, if we trash that id we
    // have no information of the previous assigned id.
    if (id == m_latest) {
        m_latest = NULL_ID;
    }
}

/*!
    Reset the generator.
*/
void IndexGenerator::reset()
{
    m_latest  = NULL_ID;
    m_highest = NULL_ID;
    m_trash.clear();
}

/*!
    Gets the version of the binary archive.

    \result The version of the binary archive.
*/
int IndexGenerator::getBinaryArchiveVersion()
{
    return 1;
}

/*!
 *  Write the index generator to the specified stream.
 *
 *  \param stream is the stream to write to
 */
void IndexGenerator::dump(std::ostream &stream)
{
    utils::binary::write(stream, getBinaryArchiveVersion());
    utils::binary::write(stream, m_latest);
    utils::binary::write(stream, m_highest);

    utils::binary::write(stream, m_trash.size());
    for (long id : m_trash) {
        utils::binary::write(stream, id);
    }
}

/*!
 *  Restore the index generator from the specified stream.
 *
 *  \param stream is the stream to read from
 */
void IndexGenerator::restore(std::istream &stream)
{
    // Version
    int version;
    utils::binary::read(stream, version);
    if (version != getBinaryArchiveVersion()) {
        throw std::runtime_error ("The version of the file does not match the required version");
    }

    // Generator data
    utils::binary::read(stream, m_latest);
    utils::binary::read(stream, m_highest);

    size_t nTrashedIds;
    utils::binary::read(stream, nTrashedIds);

    m_trash.resize(nTrashedIds);
    for (size_t i = 0; i < nTrashedIds; ++i) {
        utils::binary::read(stream, m_trash[i]);
    }
}

/*!
    @}
*/

}
