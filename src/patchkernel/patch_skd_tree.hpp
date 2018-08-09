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

# ifndef __BITPIT_PATCH_SKD_TREE_HPP__
# define __BITPIT_PATCH_SKD_TREE_HPP__

#include "patch_kernel.hpp"

namespace bitpit {

class PatchSkdTree;

struct SkdPatchInfo {

friend class PatchSkdTree;

public:
    void buildCache();
    void destroyCache();

    const PatchKernel & getPatch() const;
    const std::vector<std::size_t> & getCellRawIds() const;
    std::size_t getCellRawId(std::size_t n) const;

    const std::array<double, 3> & getCachedCentroid(std::size_t rawId) const;

    const std::array<double, 3> & getCachedBoxMin(std::size_t rawId) const;
    const std::array<double, 3> & getCachedBoxMax(std::size_t rawId) const;
    std::array<double, 3> evalCachedBoxMean(std::size_t rawId) const;

protected:
    typedef PiercedStorage<std::array<double, 3>, long> CentroidCache;
    typedef PiercedStorage<std::array<double, 3>, long> BoxCache;

    SkdPatchInfo(const PatchKernel *patch, const std::vector<std::size_t> *cellRawIds);

    const PatchKernel *m_patch;
    const std::vector<std::size_t> *m_cellRawIds;

    std::unique_ptr<CentroidCache> m_cellCentroids;
    std::unique_ptr<BoxCache> m_cellBoxes;

};

class SkdNode {

friend class PatchSkdTree;
friend class SkdNodeAllocator;

public:
    static const std::size_t NULL_ID;

    constexpr static const int MAX_CHILDREN = 2;

    enum ChildLocation {
        CHILD_LEFT  = 0,
        CHILD_RIGHT = 1,
        CHILD_BEGIN = 0,
        CHILD_END   = 2
    };

    std::size_t getCellCount() const;
    std::vector<long> getCells() const;
    long getCell(std::size_t n) const;

    const std::array<double, 3> & getBoxMin() const;
    const std::array<double, 3> & getBoxMax() const ;
    std::array<double,3> evalBoxWeightedMean() const;

    bool isLeaf() const;
    bool hasChild(ChildLocation child) const;
    std::size_t getChildId(ChildLocation child) const;

    bool boxContainsPoint(const std::array<double,3> &point, double offset) const;
    bool boxIntersectsSphere(const std::array<double,3> &center, double radius) const;

    double evalPointMinDistance(const std::array<double, 3> &point) const;
    double evalPointMaxDistance(const std::array<double, 3> &point) const;
    double evalPointDistance(const std::array<double, 3> &point) const;

    void findPointClosestCell(const std::array<double, 3> &point, long *id, double *distance) const;

protected:
    struct Allocator : std::allocator<SkdNode>
    {
        template<typename U, typename... Args>
        void construct(U* p, Args&&... args)
        {
            ::new((void *)p) U(std::forward<Args>(args)...);
        }

        template<typename U> struct rebind { typedef Allocator other; };

    };

    SkdNode();
    SkdNode(const SkdPatchInfo *patchInfo, std::size_t cellRangeBegin, std::size_t cellRangeEnd);

private:
    const SkdPatchInfo *m_patchInfo;

    std::size_t m_cellRangeBegin;
    std::size_t m_cellRangeEnd;

    std::array<double,3> m_boxMin;
    std::array<double,3> m_boxMax;

    std::array<std::size_t, MAX_CHILDREN> m_children;

    void initializeBoundingBox();

};

class PatchSkdTree {

public:
    virtual ~PatchSkdTree() = default;

    void build(int leafCapacity = 1);
    void clear(bool release = false);

    const PatchKernel & getPatch() const;

    int getLeafCapacity() const;

    std::size_t getNodeCount() const;
    std::size_t getLeafCount() const;

    const SkdNode & getNode(std::size_t nodeId) const;

    std::size_t evalMaxDepth(std::size_t rootId = 0) const;

protected:
    SkdPatchInfo m_patchInfo;
    std::vector<std::size_t> m_cellRawIds;

    std::size_t m_leafCapacity;

    std::size_t m_nLeafs;
    std::vector<SkdNode, SkdNode::Allocator> m_nodes;

    PatchSkdTree(const PatchKernel *patch);

    SkdNode & _getNode(std::size_t nodeId);

private:
    void setLeafCapacity(int capacity);

    void createChildren(std::size_t parentId);

};

}

#endif
