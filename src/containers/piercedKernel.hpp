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

#ifndef __BITPIT_PIERCED_KERNEL_HPP__
#define __BITPIT_PIERCED_KERNEL_HPP__

#include <algorithm>
#include <cassert>
#include <limits>
#include <unordered_map>
#include <type_traits>
#include <vector>

#include <bitpit_common.hpp>

#include "piercedSync.hpp"

namespace bitpit {

class BasePiercedKernel : protected virtual PiercedSyncMaster {

protected:
    BasePiercedKernel();

};

template<typename PS_value_t, typename PS_id_t>
class PiercedStorage;


/**
* \ingroup containers
*
* \brief Metafunction for generating a pierced kernel.
*
* \details
* Usage: use <tt>PiercedKernel<id_t></tt> to declare a pierced kernel.
*
* Internally all the holes are stored in a single vector. The first part of
* this vector contains the "regular" holes, whereas the last part contains
* the "pending" holes. The space reserved to the pending holes is fixed.
* To track begin and end of pending/regular holes section, two couple of
* begin/end iterators are used.
*
*          /------ Regular holes begin
*          |
*          |                     /------ Regular holes end
*          |                     |
*          v                     v
*     |-+-+R+R+R+R+R+R+R+R+R+R+R+-+-+-+-+-+P+P+P+P+P+P+P+-+-+-+-+-+-|
*      <   REGULAR HOLES SECTION   >   <     MAX PENDING HOLES     >
*                                          ^             ^
*                                          |             |
*                Pending holes begin ------/             |
*                                                        |
*                                Pending holes end ------/
*
* At the beginning, the vector of the holes will have a size equal to the
* maximum number of pending holes. When an element is deleted, its position
* is marked as a pending hole and it is inserted at the end of the existing
* pending holes (or after the regular holes if there are no pending holes).
* When the maximum number of pending holes is reached, the hole's vector
* is flushed. First, pending holes are then converted to regular holes.
* The difference between pending and regular holes is that the position of
* a pending hole is just marked as empty, whereas the position of a regular
* hole contains the distance from the next non-epty element (to speed-up
* vector traversal). Once the positions associated to the pending holes are
* updated, pending holes are moved into the regular holes and all the holes
* are compacted at the beginning of the holes' vector. Finally, the vector
* is resized: the new size of the vector is the number of regular holes plus
* the maximum number of allowed pending holes.
*
* When a new element needs to be inserted in the element, first a suitable
* position is searched in the pending holes, if no suitable positions is
* found, the search is extended to regular holes. If, among the holes, there
* is no suitable position, a new element is added in the container.
*
* \tparam id_t The type of the ids to associate to the elements
*/
template<typename id_t = long>
class PiercedKernel : public BasePiercedKernel {

static_assert(std::is_integral<id_t>::value, "Signed integer required for id.");
static_assert(std::numeric_limits<id_t>::is_signed, "Signed integer required for id.");

// Friendships
template<typename PI_value_t, typename PI_id_t, typename PI_value_no_cv_t>
friend class PiercedIterator;

template<typename PS_value_t, typename PS_id_t>
friend class PiercedStorage;

public:
    /**
    * Type of ids stored in the kernel
    */
    typedef id_t id_type;

    /**
    * Functional for compare the position of two elements
    */
    struct positionLess
    {
        positionLess(const PiercedKernel<id_t> &kernel)
        {
            m_kernel = &kernel;
        }

        bool operator()(id_t id_1, id_t id_2) const
        {
            return m_kernel->getPos(id_1) < m_kernel->getPos(id_2);
        }

        const PiercedKernel<id_t> *m_kernel;
    };

    /**
    * Functional for compare the position of two elements
    */
    struct positionGreater
    {
        positionGreater(const PiercedKernel<id_t> &kernel)
        {
            m_kernel = &kernel;
        }

        bool operator()(id_t id_1, id_t id_2) const
        {
            return m_kernel->getPos(id_1) > m_kernel->getPos(id_2);
        }

        const PiercedKernel<id_t> *m_kernel;
    };

    // Sync mode
    using PiercedSyncMaster::SyncMode;

    // Contructors
    PiercedKernel();
    PiercedKernel(std::size_t n);

    // Methods that modify the kernel as a whole
    void clear(bool release = true);
    void flush();
    void reserve(std::size_t n);
    void resize(std::size_t n);
    std::vector<std::size_t> sort();
    std::vector<std::size_t> squeeze();
    void shrinkToFit();
    void swap(PiercedKernel &x) noexcept;

    // Methods that extract information about the kernel
    bool contiguous() const;
    void dump() const;
    bool empty() const;
    bool isIteratorSlow();
    std::size_t maxSize() const;
    std::size_t size() const;
    std::size_t capacity() const;

    // Methods that extract information about the elements of the kernel
    bool contains(id_t id) const;
    std::size_t getRawIndex(id_t id) const;
    std::size_t evalFlatIndex(id_t id);

    std::vector<id_t> getIds(bool ordered = true) const;
    id_t getSizeMarker(std::size_t targetSize, const id_t &fallback = -1);

    // Methods to handle the storage
    void registerStorage(PiercedSyncSlave *storage, PiercedSyncMaster::SyncMode syncMode);
    void unregisterStorage(const PiercedSyncSlave *storage);
    bool isStorageRegistered(const PiercedSyncSlave *storage) const;
    PiercedSyncMaster::SyncMode getStorageSyncMode(const PiercedSyncSlave *storage) const;

    using PiercedSyncMaster::sync;

protected:
    /**
    * Fill action
    */
    class FillAction : public PiercedSyncAction {

    public:
        enum FillActionType {
            TYPE_UNDEFINED = PiercedSyncAction::TYPE_UNDEFINED,
            TYPE_OVERWRITE = PiercedSyncAction::TYPE_OVERWRITE,
            TYPE_INSERT    = PiercedSyncAction::TYPE_INSERT,
            TYPE_APPEND    = PiercedSyncAction::TYPE_APPEND
        };

        FillAction(FillActionType type)
            : PiercedSyncAction(static_cast<PiercedSyncAction::ActionType>(type))
        {
        };
    };

    /**
    * Move action
    */
    class MoveAction : public PiercedSyncAction {

    public:
        enum MoveActionType {
            TYPE_UNDEFINED = PiercedSyncAction::TYPE_UNDEFINED,
            TYPE_OVERWRITE = PiercedSyncAction::TYPE_MOVE_OVERWRITE,
            TYPE_INSERT    = PiercedSyncAction::TYPE_MOVE_INSERT,
            TYPE_APPEND    = PiercedSyncAction::TYPE_MOVE_APPEND
        };

        MoveAction(MoveActionType type)
            : PiercedSyncAction(static_cast<PiercedSyncAction::ActionType>(type))
        {
        };
    };

    /**
    * Swap action
    */
    class SwapAction : public PiercedSyncAction {

    public:
        enum SwapActionType {
            TYPE_UNDEFINED = PiercedSyncAction::TYPE_UNDEFINED,
            TYPE_SWAP      = PiercedSyncAction::TYPE_SWAP
        };

        SwapAction(SwapActionType type)
            : PiercedSyncAction(static_cast<PiercedSyncAction::ActionType>(type))
        {
        };
    };

    /**
    * Erase action
    */
    class EraseAction : public PiercedSyncAction {

    public:
        enum EraseActionType {
            TYPE_UNDEFINED = PiercedSyncAction::TYPE_UNDEFINED,
            TYPE_PIERCE    = PiercedSyncAction::TYPE_PIERCE,
            TYPE_SHRINK    = PiercedSyncAction::TYPE_RESIZE
        };

        EraseAction(EraseActionType type)
            : PiercedSyncAction(static_cast<PiercedSyncAction::ActionType>(type))
        {
        };
    };

    /**
    * Container used for storing holes
    */
    typedef std::vector<std::size_t> holes_container;

    /**
    * Holes iterator
    */
    typedef holes_container::iterator holes_iterator;

    /**
    * Holes const iterator
    */
    typedef holes_container::const_iterator holes_const_iterator;

    // Methods that extract information about the kernel
    void checkIntegrity() const;
    std::size_t rawSize() const;

    // Methods that modify the elements stored in the kernel
    void updateId(const id_t &currentId, const id_t &updatedId);

    FillAction fillHead(id_t id);
    FillAction fillTail(id_t id);
    FillAction fillAfter(id_t referenceId, id_t id);
    FillAction fillBefore(id_t referenceId, id_t id);
    FillAction fillAppend(id_t id);
    FillAction fillHole(const holes_iterator &holeItr, id_t id);
    FillAction fillInsert(std::size_t pos, id_t id);

    MoveAction moveAfter(id_t referenceId, id_t id, bool flush = false);
    MoveAction moveBefore(id_t referenceId, id_t id, bool flush = false);

    SwapAction swap(id_t id_first, id_t id_second);

    EraseAction erase(id_t id, bool flush = false);
    EraseAction popBack();

    // Methods that extract information about the elements of the kernel
    std::size_t back() const;
    std::size_t front() const;

    std::size_t findPrevUsedPos(std::size_t pos) const;
    std::size_t findNextUsedPos(std::size_t pos) const;
    bool isPosEmpty(std::size_t pos) const;
    std::size_t getPos(id_t id) const;
    std::size_t getFirstUsedPos() const;
    std::size_t getLastUsedPos() const;

private:
    /**
    * Hasher for the id map.
    *
    * Since the id are uniques, the hasher can be a function that
    * takes the id and cast it to a std::size_t.
    *
    * The hasher is defined as a struct, because a struct can be
    * passed as an object into metafunctions (meaning that the type
    * deduction for the template paramenters can take place, and
    * also meaning that inlining is easier for the compiler). A bare
    * function would have to be passed as a function pointer.
    * To transform a function template into a function pointer,
    * the template would have to be manually instantiated (with a
    * perhaps unknown type argument).
    */
    struct PiercedHasher {
        /**
        * Function call operator that casts the specified
        * value to a std::size_t.
        *
        * \tparam U type of the value
        * \param value is the value to be casted
        * \result Returns the value casted to a std::size_t.
        */
        template<typename U>
        constexpr std::size_t operator()(U&& value) const noexcept
        {
            return static_cast<std::size_t>(std::forward<U>(value));
        }
    };

    /**
    * Compares the id of the elements in the specified position.
    *
    * \param pos_x is the position to the first element to compare
    * \param y is the position to the second element to compare
    * \result Returns true if the element x has an id lower than the element
    * y, false otherwise. Negative ids are special ids and are considered
    * higher than positive ids.
    */
    struct idLess
    {
        const std::vector<id_t> &m_ids;

        idLess(const std::vector<id_t> &ids)
            : m_ids(ids)
        {
        }

        inline bool operator() (std::size_t pos_x, std::size_t pos_y)
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

    /**
    * Maximum number of pending deletes before the changes are flushed.
    */
    static const std::size_t MAX_PENDING_HOLES;

    /**
    * Vector that will hold the ids.
    */
    std::vector<id_t> m_ids;

    /**
    * Map that links the id of the elements and their position inside the
    * internal vector.
    */
    std::unordered_map<id_t, std::size_t, PiercedHasher> m_pos;

    /**
    * Position of the first element in the internal vector.
    */
    std::size_t m_begin_pos;

    /**
    * Position of the last element in the internal vector.
    */
    std::size_t m_end_pos;

    /**
    * Position of the first dirty element.
    *
    * After the first dirty position the id of the holes can not be properly
    * defined, meaning that the iterator can take longer to iterate through
    * the elements.
    */
    std::size_t m_dirty_begin_pos;

    /**
    * Container that will hold a list of the holes present in
    * the piecrecd vector.
    */
    holes_container m_holes;

    /**
    * Iterator pointing to the first regular hole
    */
    holes_iterator m_holes_regular_begin;

    /**
    * Iterator pointing to the last regular hole
    */
    holes_iterator m_holes_regular_end;

    /**
    * Iterator pointing to the first pending hole
    */
    holes_iterator m_holes_pending_begin;

    /**
    * Iterator pointing to the last pending hole
    */
    holes_iterator m_holes_pending_end;

    /**
    * Tracks if the regular holes are sorted
    */
    bool m_holes_regular_sorted;

    /**
    * Tracks if the pending holes are sorted
    */
    bool m_holes_pending_sorted;

    void pierce(std::size_t pos, bool flush = true);

    void holesClear(bool release = false);
    void holesClearRegular(bool release = false);
    void holesClearPending(bool release = false);
    void holesResize(std::size_t offset, std::size_t nRegulars, std::size_t nPendings, bool release = true);
    std::size_t holesCount() const;
    std::size_t holesCountPending() const;
    std::size_t holesCountRegular() const;
    void holesFlush();
    void holesSortPending();
    void holesSortRegular();

    bool validateId(id_t id);

    void setBeginPos(std::size_t pos);
    void setEndPos(std::size_t pos);

    void setPosId(std::size_t pos, id_t id);
    void setPosEmptyId(std::size_t pos, std::size_t nextUsedPos);
    void swapPosIds(std::size_t pos_1, id_t id_1, std::size_t pos_2, id_t id_2);

    void shrink(std::size_t n, bool force = false);

};

}

// Include the implementation
#include "piercedKernel.tpp"


#endif
