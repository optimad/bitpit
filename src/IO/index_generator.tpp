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

#ifndef __BITPIT_INDEX_GENERATOR_TPP__
#define __BITPIT_INDEX_GENERATOR_TPP__

namespace bitpit {

template<typename id_t>
const typename IndexGenerator<id_t>::id_type IndexGenerator<id_t>::NULL_ID = std::numeric_limits<typename IndexGenerator<id_t>::id_type>::min();

/*!
    Creates a new generator.
*/
template<typename id_t>
IndexGenerator<id_t>::IndexGenerator()
    : m_latest(NULL_ID), m_lowest(NULL_ID), m_highest(NULL_ID)
{

}

/*!
    Generates a unique index.

    If the trash is empty a new index is generated, otherwise an index taken
    from the trash is recycled.

    \return A new unique index.
*/
template<typename id_t>
typename IndexGenerator<id_t>::id_type IndexGenerator<id_t>::generate()
{
    // If the trash is empty generate a new id otherwise recycle the first id
    // in the trash.
    if (m_trash.empty()) {
        if (m_lowest == NULL_ID) {
            m_lowest  = 0;
            m_highest = 0;
            m_latest = m_lowest;
        } else if (m_lowest > 0) {
            --m_lowest;
            m_latest = m_lowest;
        } else {
            assert(m_highest < std::numeric_limits<id_type>::max());
            ++m_highest;
            m_latest = m_highest;
        }
    } else {
        m_latest = *(m_trash.begin());
        eraseFromTrash(m_latest);
    }

    return m_latest;
}

/*!
    Gets the latest assigned id.

    \return The latest assigned index.
*/
template<typename id_t>
typename IndexGenerator<id_t>::id_type IndexGenerator<id_t>::getLatest()
{
    return m_latest;
}

/*!
    Gets the highest assigned id.

    \return The highest assigned index.
*/
template<typename id_t>
typename IndexGenerator<id_t>::id_type IndexGenerator<id_t>::getHighest()
{
    return m_highest;
}

/*!
    Checks if an id is currently assigned.

    \return True if the id has already been assigned, false otherwise.
*/
template<typename id_t>
bool IndexGenerator<id_t>::isAssigned(id_type id)
{
    // If the latest id is valid, it is assigned by definition
    if (m_latest >= 0 && id == m_latest) {
        return true;
    }

    // Ids past the highest one are not assigned
    if (id > m_highest) {
        return false;
    }

    // Ids before the lowest one are not assigned
    if (id < m_lowest) {
        return false;
    }

    // The id is assigned only if is not in the trash
    if (m_trash.count(id) > 0) {
        return false;
    }

    return true;
}

/*!
    Marks the specified id as currently assigned.

    \param id is the id to be marked as assigned
*/
template<typename id_t>
void IndexGenerator<id_t>::setAssigned(typename IndexGenerator<id_t>::id_type id)
{
    // Check if the id has already been assigned
    if (isAssigned(id)) {
        throw std::runtime_error("Requested id has already been assigned.");
    }

    // If the id is past the highest assigned id we need to trash all the ids
    // from the highest assigned to this one (the generator only handles
    // contiguous ids), otherwise look for the id in the trash and, if found,
    // recycle it (if the id is not in the trash, this means it was already
    // an assigned id).
    if (m_highest == NULL_ID && m_lowest == NULL_ID) {
        m_lowest  = id;
        m_highest = id;
    } else if (id > m_highest) {
        for (id_type wasteId = m_highest + 1; wasteId < id; ++wasteId) {
            trash(wasteId);
        }

        m_highest = id;
    } else if (id < m_lowest) {
        for (id_type wasteId = id + 1; wasteId < m_lowest; ++wasteId) {
            trash(wasteId);
        }

        m_lowest = id;
    } else {
        assert(m_trash.count(id) > 0);
        eraseFromTrash(id);
    }

    // This is the latest assigned id
    m_latest = id;
}

/*!
    Trashes an index.

    A trashed index is an index no more used that can be recycled.

    \param id is the index that will be trashed
*/
template<typename id_t>
void IndexGenerator<id_t>::trash(typename IndexGenerator<id_t>::id_type id)
{
    assert(m_trash.count(id) == 0);

    // If we are trashing the highest or the lowest id we can update the
    // limits of the generator, otherwise we add the id to the trash.
    if (id == m_highest && id == m_lowest) {
        // If highest and lowest values are equal it means that we have
        // generated only one id. Since we are now trashing the only id
        // that has been generated, the trash must be emtpy.
        assert(m_trash.empty());

        m_lowest  = NULL_ID;
        m_highest = NULL_ID;
    } else if (id == m_highest) {
        --m_highest;
        while (m_trash.count(m_highest) > 0) {
            eraseFromTrash(m_highest);
            --m_highest;
        }
    } else if (id == m_lowest) {
        ++m_lowest;
        while (m_trash.count(m_lowest) > 0) {
            eraseFromTrash(m_lowest);
            ++m_lowest;
        }
    } else {
        pushToTrash(id);
    }

    // We only keep track of the latest assigned id, if we trash that id we
    // have no information of the previous assigned id.
    if (id == m_latest) {
        m_latest = NULL_ID;
    }
}

/*!
    Reset the generator.
*/
template<typename id_t>
void IndexGenerator<id_t>::reset()
{
    m_latest  = NULL_ID;
    m_lowest  = NULL_ID;
    m_highest = NULL_ID;

    m_trash.clear();
}

/*!
    Add the specified index to the trash.

    \param id is the id to be added
*/
template<typename id_t>
void IndexGenerator<id_t>::pushToTrash(id_type id)
{
    m_trash.insert(id);
}

/*!
    Squeeze the trash to keep the unused size in the trash limited.

    \param id is the id to be erased
*/
template<typename id_t>
void IndexGenerator<id_t>::eraseFromTrash(id_type id)
{
    // Erase the index
    m_trash.erase(id);
}

/*!
    Gets the version of the binary archive.

    \result The version of the binary archive.
*/
template<typename id_t>
int IndexGenerator<id_t>::getBinaryArchiveVersion() const
{
    return 1;
}

/*!
 *  Write the index generator to the specified stream.
 *
 *  \param stream is the stream to write to
 */
template<typename id_t>
void IndexGenerator<id_t>::dump(std::ostream &stream) const
{
    utils::binary::write(stream, getBinaryArchiveVersion());
    utils::binary::write(stream, m_latest);
    utils::binary::write(stream, m_lowest);
    utils::binary::write(stream, m_highest);

    utils::binary::write(stream, m_trash.size());
    for (id_type id : m_trash) {
        utils::binary::write(stream, id);
    }
}

/*!
 *  Restore the index generator from the specified stream.
 *
 *  \param stream is the stream to read from
 */
template<typename id_t>
void IndexGenerator<id_t>::restore(std::istream &stream)
{
    // Version
    int version;
    utils::binary::read(stream, version);
    if (version != getBinaryArchiveVersion()) {
        throw std::runtime_error ("The version of the file does not match the required version");
    }

    // Generator data
    utils::binary::read(stream, m_latest);
    utils::binary::read(stream, m_lowest);
    utils::binary::read(stream, m_highest);

    size_t nTrashedIds;
    utils::binary::read(stream, nTrashedIds);

    std::set<id_type>().swap(m_trash);

    for (size_t i = 0; i < nTrashedIds; ++i) {
        id_type id;
        utils::binary::read(stream, id);
        pushToTrash(id);
    }
}

}

#endif
