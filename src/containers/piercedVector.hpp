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

#ifndef __BITPIT_PIERCED_VECTOR_HPP__
#define __BITPIT_PIERCED_VECTOR_HPP__

#define __PV_REFERENCE__       typename PiercedVector<value_t, id_t>::reference
#define __PV_CONST_REFERENCE__ typename PiercedVector<value_t, id_t>::const_reference
#define __PV_POINTER__         typename PiercedVector<value_t, id_t>::pointer
#define __PV_CONST_POINTER__   typename PiercedVector<value_t, id_t>::const_pointer

#include <cassert>
#include <limits>
#include <unordered_map>
#include <type_traits>
#include <vector>

#include <bitpit_common.hpp>

#include <piercedIterator.hpp>

namespace bitpit {

/*!
    \ingroup containers

    \brief Base class for the pierced vectors.
*/
class BasePiercedVector {
};

/*!
    \ingroup containers

    \brief Metafunction for generating a pierced vector.

    \details
    Usage: use <tt>PiercedVector<value_t, id_t></tt> to declare a pierced
    vector.

    Internally all the holes are stored in a single vector. The first part of
    this vector contains the "regular" holes, whereas the last part contains
    the "pending" holes. The space reserved to the pending holes is fixed.
    To track begin and end of pending/regular holes section, two couple of
    begin/end iterators are used.

             /------ Regular holes begin
             |
             |                     /------ Regular holes end
             |                     |
             v                     v
        |-+-+R+R+R+R+R+R+R+R+R+R+R+-+-+-+-+-+P+P+P+P+P+P+P+-+-+-+-+-+-|
         <   REGULAR HOLES SECTION   >   <     MAX PENDING HOLES     >
                                             ^             ^
                                             |             |
                   Pending holes begin ------/             |
                                                           |
                                   Pending holes end ------/

    At the beginning, the vector of the holes will have a size equal to the
    maximum number of pending holes. When an element is deleted, its position
    is marked as a pending hole and it is inserted at the end of the existing
    pending holes (or after the regular holes if there are no pending holes).
    When the maximum number of pending holes is reached, the hole's vector
    is flushed. First, pending holes are then converted to regular holes.
    The difference between pending and regular holes is that the position of
    a pending hole is just marked as empty, whereas the position of a regular
    hole contains the distance from the next non-epty element (to speed-up
    vector traversal). Once the positions associated to the pending holes are
    updated, pending holes are moved into the regular holes and all the holes
    are compacted at the beginning of the holes' vector. Finally, the vector
    is resized: the new size of the vector is the number of regular holes plus
    the maximum number of allowed pending holes.

    When a new element needs to be inserted in the element, first a suitable
    position is searched in the pending holes, if no suitable positions is
    found, the search is extended to regular holes. If, among the holes, there
    is no suitable position, a new element is added in the container.

    \tparam value_t The type of the elements stored in the vector
    \tparam id_t The type of the ids to associate to the elements
*/
template<typename value_t, typename id_t = long>
class PiercedVector : public BasePiercedVector {
    static_assert(std::is_integral<id_t>::value, "Signed integer required for id.");
    static_assert(std::numeric_limits<id_t>::is_signed, "Signed integer required for id.");

private:
    /*!
        Maximum number of pending deletes before the changes are flushed.
    */
    static const std::size_t MAX_PENDING_HOLES;

public:
    // Friendships
    template<typename PI_value_t, typename PI_id_t, typename PI_value_no_cv_t>
    friend class PiercedIterator;

    /*!
        Type of data stored in the container
    */
    typedef value_t value_type;

    /*!
        Type of ids stored in the container
    */
    typedef id_t id_type;

    /*!
        Iterator for the pierced array.
    */
    typedef PiercedIterator<value_t, id_t> iterator;

    /*!
        Constant iterator for the pierced array.
    */
    typedef PiercedIterator<const value_t, id_t> const_iterator;

    /*!
        Iterator for the pierced array raw container.
    */
    typedef typename std::vector<value_t>::iterator raw_iterator;

    /*!
        Constant iterator for the pierced array raw container.
    */
    typedef typename std::vector<value_t>::const_iterator raw_const_iterator;

    /*!
        Reference
    */
    typedef typename std::vector<value_t>::reference reference;

    /*!
        Constant reference
    */
    typedef typename std::vector<value_t>::const_reference const_reference;

    /*!
        Pointer
    */
    typedef typename std::vector<value_t>::pointer pointer;

    /*!
        Constant pointer
    */
    typedef typename std::vector<value_t>::const_pointer const_pointer;

    /*!
        Functional for compare the position of two elements
    */
    struct positionLess
    {
        positionLess(PiercedVector<value_t> &vector)
        {
            m_vector = &vector;
        }

        bool operator()(const id_t &id_1, const id_t &id_2) const
        {
            return m_vector->getPosFromId(id_1) < m_vector->getPosFromId(id_2);
        }

        PiercedVector<value_t> *m_vector;
    };

    /*!
        Functional for compare the position of two elements
    */
    struct positionGreater
    {
        positionGreater(PiercedVector<value_t> &vector)
        {
            m_vector = &vector;
        }

        bool operator()(const id_t &id_1, const id_t &id_2) const
        {
            return m_vector->getPosFromId(id_1) > m_vector->getPosFromId(id_2);
        }

        PiercedVector<value_t> *m_vector;
    };

    // Contructors
    PiercedVector();
    PiercedVector(std::size_t n);

    // Methods that modify the contents of the container
    iterator pushBack(const id_t &id, value_t &&value);

    iterator reclaim(const id_t &id);
    iterator reclaimAfter(const id_t &referenceId, const id_t &id);
    iterator reclaimBack(const id_t &id);
    iterator reclaimBefore(const id_t &referenceId, const id_t &id);

    iterator moveAfter(const id_t &referenceId, const id_t &id, bool delayed = false);
    iterator moveBefore(const id_t &referenceId, const id_t &id, bool delayed = false);

    iterator insert(const id_t &id, const value_t &value);
    iterator insertAfter(const id_t &referenceId, const id_t &id, const value_t &value);
    iterator insertBefore(const id_t &referenceId, const id_t &id, const value_t &value);

    iterator replace(id_t id, value_t &&value);

    void updateId(const id_t &currentId, const id_t &updatedId);

    template<typename... Args>
    typename PiercedVector<value_t, id_t>::iterator emplace(const id_t &id, Args&&... args);
    template<typename... Args>
    typename PiercedVector<value_t, id_t>::iterator emplaceAfter(const id_t &referenceId, const id_t &id, Args&&... args);
    template<typename... Args>
    void emplaceBack(const id_t &id, Args&&... args);
    template<typename... Args>
    typename PiercedVector<value_t, id_t>::iterator emplaceBefore(const id_t &referenceId, const id_t &id, Args&&... args);

    template<typename... Args>
    typename PiercedVector<value_t, id_t>::iterator emreplace(id_t id, Args&&... args);

    iterator erase(id_t id, bool delayed = false);

    void popBack();

    void swap(const id_t &id_first, const id_t &id_second);

    // Methods that modify the container as a whole
    void clear(bool release = true);
    void flush();
    void reserve(std::size_t n);
    void resize(std::size_t n);
    void sort();
    void squeeze();
    void shrinkToFit();
    void swap(PiercedVector& x) noexcept;

    // Methods that extract information on the container
    std::size_t capacity() const;
    bool contiguous() const;
    void dump();
    bool empty() const;
    bool isIteratorSlow();
    std::size_t maxSize() const;
    std::size_t size() const;

    // Methods that extract information on the contents of the container
    bool exists(id_t id) const;
    const_iterator find(id_t id) const;
    iterator find(id_t id);
    std::size_t evalFlatIndex(id_t id);
    std::size_t getRawIndex(id_t id) const;

    std::vector<id_t> getIds(bool ordered = true) const;
    id_t getSizeMarker(const size_t &targetSize, const id_t &fallback = -1);

    // Methods that extract the contents of the container
    __PV_POINTER__ data() noexcept;

    __PV_REFERENCE__ back();
    __PV_CONST_REFERENCE__ back() const;

    __PV_REFERENCE__ front();
    __PV_CONST_REFERENCE__ front() const;

    __PV_REFERENCE__ at(const id_t &id);
    __PV_CONST_REFERENCE__ at(const id_t &id) const;

    __PV_REFERENCE__ rawAt(const std::size_t &pos);
    __PV_CONST_REFERENCE__ rawAt(const std::size_t &pos) const;
    std::size_t rawIndex(id_t id) const;

    __PV_CONST_REFERENCE__ operator[](const id_t &id) const;
    __PV_REFERENCE__ operator[](const id_t &id);

    // Iterators
    iterator getIterator(const id_t &id) noexcept;
    const_iterator getConstIterator(const id_t &id) const noexcept;

    iterator getIteratorFromRawIndex(const std::size_t &pos) noexcept;
    const_iterator getConstIteratorFromRawIndex(const std::size_t &pos) const noexcept;

    iterator begin() noexcept;
    iterator end() noexcept;
    const_iterator begin() const noexcept;
    const_iterator end() const noexcept;
    const_iterator cbegin() const noexcept;
    const_iterator cend() const noexcept;

    raw_iterator rawBegin() noexcept;
    raw_iterator rawEnd() noexcept;
    raw_const_iterator rawBegin() const noexcept;
    raw_const_iterator rawEnd() const noexcept;
    raw_const_iterator rawCbegin() const noexcept;
    raw_const_iterator rawCend() const noexcept;

private:
    /*!
        Hasher for the id map.

        Since the id are uniques, the hasher can be a function that
        takes the id and cast it to a size_t.

        The hasher is defined as a struct, because a struct can be
        passed as an object into metafunctions (meaning that the type
        deduction for the template paramenters can take place, and
        also meaning that inlining is easier for the compiler). A bare
        function would have to be passed as a function pointer.
        To transform a function template into a function pointer,
        the template would have to be manually instantiated (with a
        perhaps unknown type argument).

    */
    struct PiercedHasher {
        /*!
            Function call operator that casts the specified
            value to a size_t.

            \tparam U type of the value
            \param value is the value to be casted
            \result Returns the value casted to a size_t.
        */
        template<typename U>
        constexpr std::size_t operator()(U&& value) const noexcept
        {
            return static_cast<std::size_t>(std::forward<U>(value));
        }
    };

    /*!
        Storage position
    */
    struct StoragePosition {
        StoragePosition(size_t _pos, bool _ready)
            : pos(_pos), ready(_ready)
        {
        }

        size_t pos;
        bool ready;
    };

    /*!
        Container used for storing holes
    */
    typedef std::vector<std::size_t> hole_container;

    /*!
        Hole iterator
    */
    typedef hole_container::iterator hole_iterator;

    /*!
        Vector that will hold the elements.
    */
    std::vector<value_t>m_v;

    /*!
        Vector that will hold the ids.
    */
    std::vector<id_t> m_ids;

    /*!
        Container that will hold a list of the holes present in
        the piecrecd vector.
    */
    hole_container m_holes;

    /*!
        Iterator pointing to the first regular hole
    */
    hole_iterator m_holes_regular_begin;

    /*!
        Iterator pointing to the last regular hole
    */
    hole_iterator m_holes_regular_end;

    /*!
        Iterator pointing to the first pending hole
    */
    hole_iterator m_holes_pending_begin;

    /*!
        Iterator pointing to the last pending hole
    */
    hole_iterator m_holes_pending_end;

    /*!
        Tracks if the regular holes are sorted
    */
    bool m_holes_regular_sorted;

    /*!
        Tracks if the pending holes are sorted
    */
    bool m_holes_pending_sorted;

    /*!
        Map that links the id of the elements and their position
        inside the internal vector.
    */
    std::unordered_map<id_t, std::size_t, PiercedHasher> m_pos;

    /*!
        Position of the first element in the internal vector.
    */
    std::size_t m_begin_pos;

    /*!
        Position of the last element in the internal vector.
    */
    std::size_t m_end_pos;

    /*!
        Position of the first dirty element.

        After the first dirty position the id of the holes can not be
        properly defined, meaning that the iterator can take longer to
        iterate through the elements.
    */
    std::size_t m_dirty_begin_pos;

    /*!
        Compares the id of the elements in the specified position.

        \param pos_x is the position to the first element to compare
        \param y is the position to the second element to compare
        \result Returns true if the element x has an id lower than the element
        y, false otherwise. Negative ids are special ids and are considered
        higher than positive ids.
    */
    struct idLess
    {
        const std::vector<id_t> &m_ids;

        idLess(const std::vector<id_t> &ids)
            : m_ids(ids)
        {
        }

        inline bool operator() (const std::size_t &pos_x, const std::size_t &pos_y)
        {
            id_t id_x = m_ids[pos_x];
            id_t id_y = m_ids[pos_y];

            if (id_x >= 0 && id_y < 0) {
                return true;
            } else if (id_x < 0 && id_y >= 0) {
                return false;
            } else if (id_x >= 0) {
                return (id_x < id_y);
            } else {
                return (id_x > id_y);
            }
        }
    };

    iterator getIteratorFromPos(const std::size_t &pos) noexcept;
    const_iterator getConstIteratorFromPos(const std::size_t &pos) const noexcept;

    StoragePosition fillHead(const id_t &id);
    StoragePosition fillTail(const id_t &id);
    StoragePosition fillBefore(const std::size_t &referencePos, const id_t &id);
    StoragePosition fillAfter(const std::size_t &referencePos, const id_t &id);
    StoragePosition fillHole(const hole_iterator &holeItr, const id_t &id);
    StoragePosition fillAppend(const id_t &id);
    StoragePosition fillInsert(const std::size_t &pos, const id_t &id);

    void pierce(const std::size_t &pos, bool flush = true);

    void holesClear(bool release = false);
    void holesClearRegular(bool release = false);
    void holesClearPending(bool release = false);
    void holesResize(size_t offset, size_t nRegulars, size_t nPendings, bool release = true);
    std::size_t holesCount() const;
    std::size_t holesCountPending() const;
    std::size_t holesCountRegular() const;
    void holesFlush();
    void holesSortPending();
    void holesSortRegular();

    bool validateId(const id_t &id);

    std::size_t findPrevUsedPos(std::size_t pos) const;
    std::size_t findNextUsedPos(std::size_t pos) const;
    bool isPosEmpty(std::size_t pos) const;
    std::size_t getPosFromId(id_t id) const;
    void setBeginPos(const std::size_t &pos);
    void setEndPos(const std::size_t &pos);

    void setPosId(const std::size_t &pos, const id_t &id);
    void setPosEmptyId(const std::size_t &pos, const std::size_t &nextUsedPos);
    void swapPosIds(const std::size_t &pos_1, const id_t &id_1, const std::size_t &pos_2, const id_t &id_2);

    void storageShrink(size_t n, bool force = false);
    std::size_t storageSize() const;

};

}

// Include the implementation
#include "piercedVector.tpp"

#endif
