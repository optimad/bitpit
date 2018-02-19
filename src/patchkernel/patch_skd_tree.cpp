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

#include <cmath>
#include <queue>

#include "bitpit_CG.hpp"

#include "patch_skd_tree.hpp"

namespace bitpit {

/*!
* \class SkdPatchInfo
*
* \brief The SkdPatchInfo class allows to store patch information
* needed for the construction and the utilization of a skd-tree.
*/

/*!
* Constructor
*
* \param patch is the patch
* \param cellRawIds is the list of cell raw ids
*/
SkdPatchInfo::SkdPatchInfo(const PatchKernel *patch, const std::vector<std::size_t> *cellRawIds)
    : m_patch(patch), m_cellRawIds(cellRawIds)
{
    if (!m_patch) {
        throw std::runtime_error("Unable to initialize the patch info. Provided patch is not valid.");
    }
}

/*!
* Build the cache.
*/
void SkdPatchInfo::buildCache()
{
    const PiercedVector<Cell> &cells = m_patch->getCells();

    m_cellBoxes     = std::unique_ptr<BoxCache>(new BoxCache(2, &(m_patch->getCells())));
    m_cellCentroids = std::unique_ptr<BoxCache>(new CentroidCache(1, &(m_patch->getCells())));
    for (auto itr = m_patch->internalConstBegin(); itr != m_patch->internalConstEnd(); ++itr) {
        std::size_t rawCellId = itr.getRawIndex();

        // Cell info
        const Cell &cell = cells.rawAt(rawCellId);
        int nCellVertices = cell.getVertexCount();
        ConstProxyVector<long> cellConnect = cell.getVertexIds();

        // Bounding box
        std::array<double, 3> &cellBoxMin   = m_cellBoxes->rawAt(rawCellId, 0);
        std::array<double, 3> &cellBoxMax   = m_cellBoxes->rawAt(rawCellId, 1);
        std::array<double, 3> &cellCentroid = m_cellCentroids->rawAt(rawCellId, 0);

        cellBoxMin   = m_patch->getVertexCoords(cellConnect[0]);
        cellBoxMax   = cellBoxMin;
        cellCentroid = cellBoxMin;
        for (int i = 1; i < nCellVertices; ++i) {
            const std::array<double, 3> &coords = m_patch->getVertexCoords(cellConnect[i]);
            for (int d = 0; d < 3; ++d) {
                cellBoxMin[d]    = std::min(coords[d], cellBoxMin[d]);
                cellBoxMax[d]    = std::max(coords[d], cellBoxMax[d]);
                cellCentroid[d] += coords[d];
            }
        }
        cellCentroid /= double(nCellVertices);
    }
}

/*!
* Destroy the cache.
*/
void SkdPatchInfo::destroyCache()
{
    m_cellBoxes.reset();
    m_cellCentroids.reset();
}

/*!
* Get a constant reference to the patch.
*
* \result A constant reference to the patch.
*/
const PatchKernel & SkdPatchInfo::getPatch() const
{
    return *m_patch;
}

/*!
* Get a constant reference to the list of cell raw ids.
*
* \result A constant reference to the list of cell raw ids.
*/
const std::vector<std::size_t> & SkdPatchInfo::getCellRawIds() const
{
    return *m_cellRawIds;
}

/*!
* Get the n-th cell raw id.
*
* \param n is the requested raw id
* \result The n-th cell raw id.
*/
std::size_t SkdPatchInfo::getCellRawId(std::size_t n) const
{
    return (*m_cellRawIds)[n];
}

/*!
* Get the cached centroid of the specified cell.
*
* \param rawId is the raw id of the cell
* \result The cached centroid of the specified cell.
*/
const std::array<double, 3> & SkdPatchInfo::getCachedCentroid(std::size_t rawId) const
{
    return m_cellCentroids->rawAt(rawId, 0);
}

/*!
* Get the minimum point of the cached bounding box of the specified cell.
*
* \param rawId is the raw id of the cell
* \result The minimum point of the cached bounding box of the specified cell.
*/
const std::array<double, 3> & SkdPatchInfo::getCachedBoxMin(std::size_t rawId) const
{
    return m_cellBoxes->rawAt(rawId, 0);
}

/*!
* Get the maximum point of the cached bounding box of the specified cell.
*
* \param rawId is the raw id of the cell
* \result The maximum point of the cached bounding box of the specified cell.
*/
const std::array<double, 3> & SkdPatchInfo::getCachedBoxMax(std::size_t rawId) const
{
    return m_cellBoxes->rawAt(rawId, 1);
}

/*!
* Get the average point of the cached bounding box of the specified cell.
*
* \param rawId is the raw id of the cell
* \result The average point of the cached bounding box of the specified cell.
*/
std::array<double, 3> SkdPatchInfo::evalCachedBoxMean(std::size_t rawId) const
{
    return 0.5 * (m_cellBoxes->rawAt(rawId, 0) + m_cellBoxes->rawAt(rawId, 1));
}

/*!
* \class SkdNode
*
* \brief The SkdPatchInfo class defines a node of the skd-tree.
*/

/*!
* Null node id
*/
const std::size_t SkdNode::NULL_ID = std::numeric_limits<std::size_t>::max();

/*!
* Default constructor.
*/
SkdNode::SkdNode()
    : m_patchInfo(nullptr),
      m_cellRangeBegin(0), m_cellRangeEnd(0),
      m_boxMin({{0., 0., 0.}}), m_boxMax({{0., 0., 0.}}),
      m_children({{NULL_ID, NULL_ID}})
{
}

/*!
* Constructor
*
* \param patchInfo are the patch information
* \param cellRangeBegin is the index of the first cell enclosed in the
* bounding box associated to the node
* \param cellRangeEnd is the index of the past-the-end cell enclosed
* in the bounding box associated to the node
*/
SkdNode::SkdNode(const SkdPatchInfo *patchInfo, std::size_t cellRangeBegin, std::size_t cellRangeEnd)
    : m_patchInfo(patchInfo),
      m_cellRangeBegin(cellRangeBegin), m_cellRangeEnd(cellRangeEnd),
      m_children({{NULL_ID, NULL_ID}})
{
    initializeBoundingBox();
}

/*!
* Initialize the bounding box associated to the node.
*/
void SkdNode::initializeBoundingBox()
{
    const std::vector<std::size_t> &cellRawIds = m_patchInfo->getCellRawIds();

    // Evaluate the bounding box
    m_boxMin = m_patchInfo->getCachedBoxMin(cellRawIds[m_cellRangeBegin]);
    m_boxMax = m_patchInfo->getCachedBoxMax(cellRawIds[m_cellRangeBegin]);
    for (std::size_t n = m_cellRangeBegin + 1; n < m_cellRangeEnd; n++) {
        const std::size_t rawCellId = cellRawIds[n];
        const std::array<double, 3> &cellBoxMin = m_patchInfo->getCachedBoxMin(rawCellId);
        const std::array<double, 3> &cellBoxMax = m_patchInfo->getCachedBoxMax(rawCellId);
        for (int d = 0; d < 3; ++d) {
            m_boxMin[d] = std::min(cellBoxMin[d], m_boxMin[d]);
            m_boxMax[d] = std::max(cellBoxMax[d], m_boxMax[d]);
        }
    }

    // Inlfate the bounding box by a small epsilon
    //
    // The small inflation allows us to test intersection of a query
    // against a safety box and be assured that non-intersection of
    // with the box implies non-intersection with all the contents
    // in the box, even when computations are susceptible to round-off
    // error.
    double tolerance = m_patchInfo->getPatch().getTol();

    m_boxMin -= tolerance;
    m_boxMax += tolerance;
}

/*!
* Count the cells enclosed in the bounding box associated to the node.
*
* \result The number of cells enclosed in the bounding box associated to the
* node.
*/
std::size_t SkdNode::getCellCount() const
{
    return (m_cellRangeEnd - m_cellRangeBegin);
}

/*!
* Get the n-th cell enclosed in the bounding box associated to the node.
*
* \param n is the index of the requested cell
* \result The n-th cell enclosed in the bounding box associated to the node.
*/
long SkdNode::getCell(std::size_t n) const
{
    const PatchKernel &patch = m_patchInfo->getPatch();
    const PiercedKernel<long> &cellKernel = patch.getCells().getKernel();

    const std::size_t rawCellId = m_patchInfo->getCellRawId(m_cellRangeBegin + n);
    long cellId = cellKernel.rawFind(rawCellId).getId();

    return cellId;
}

/*!
* Get the list of cells enclosed in the bounding box associated to the node.
*
* \result The list of cells enclosed in the bounding box associated to the
* node.
*/
std::vector<long> SkdNode::getCells() const
{
    const PatchKernel &patch = m_patchInfo->getPatch();
    const std::vector<std::size_t> &cellRawIds = m_patchInfo->getCellRawIds();
    const PiercedKernel<long> &cellKernel = patch.getCells().getKernel();

    std::size_t nCells = getCellCount();
    std::vector<long> cellList(nCells);
    for (std::size_t n = m_cellRangeBegin; n < m_cellRangeEnd; ++n) {
        const std::size_t rawCellId = cellRawIds[n];
        long cellId = cellKernel.rawFind(rawCellId).getId();

        cellList[n - m_cellRangeBegin] = cellId;
    }

    return cellList;
}

/*!
* Get the minimum coordinate of the bounding box associated to the node.
*
* \result The minimum coordinate of the bounding box associated to the
* node.
*/
const std::array<double, 3> & SkdNode::getBoxMin() const
{
    return m_boxMin;
}

/*!
* Get the maximum coordinate of the bounding box associated to the node.
*
* \result The maximum coordinate of the bounding box associated to the
* node.
*/
const std::array<double, 3> & SkdNode::getBoxMax() const
{
    return m_boxMax;
}

/*!
* Evaluates the weighted centroid of the bounding box associated to
* the node.
*
* \result The the weighted centroid of the bounding box associated to
* the node.
*/
std::array<double, 3> SkdNode::evalBoxWeightedMean() const
{
    const std::vector<std::size_t> &cellRawIds = m_patchInfo->getCellRawIds();

    std::array<double, 3> boxWeightedMean = m_patchInfo->getCachedCentroid(cellRawIds[m_cellRangeBegin]);
    for (std::size_t n = m_cellRangeBegin + 1; n < m_cellRangeEnd; ++n) {
        const std::size_t rawCellId = cellRawIds[n];
        boxWeightedMean += m_patchInfo->getCachedCentroid(rawCellId);
    }
    boxWeightedMean /= (double) getCellCount();

    return boxWeightedMean;
}

/*!
* Checks if a node is a leaf.
*
* \result Returns true is the node is a leaf, false otherwise.
*/
bool SkdNode::isLeaf() const
{
    for (int i = CHILD_BEGIN; i != CHILD_END; ++i) {
        SkdNode::ChildLocation childLocation = static_cast<ChildLocation>(i);
        if (hasChild(childLocation)) {
            return false;
        }
    }

    return true;
}

/*!
* Checks if a node has the specified child.
*
* \param child is the child that has to be checked
* \result Returns true if the specified child exists, false otherwise.
*/
bool SkdNode::hasChild(ChildLocation child) const
{
    return (m_children[static_cast<std::size_t>(child)] != NULL_ID);
}

/*!
* Get the id of the specified child.
*
* \param child is the child that has to reteived
* \result The id of the specified child.
*/
std::size_t SkdNode::getChildId(ChildLocation child) const
{
    return m_children[static_cast<std::size_t>(child)];
}

/*!
* Evaluates the minimum distance among the specified point and the
* cells contained the node.
*
* \param point is the point
* \result The minimum distance among the specified point and the cells
* contained the node.
*/
double SkdNode::evalPointMinDistance(const std::array<double, 3> &point) const
{
    double distance = 0.;
    for (int d = 0; d < 3; ++d) {
        distance += std::pow(std::max({0., m_boxMin[d] - point[d], point[d] - m_boxMax[d]}), 2);
    }
    distance = std::sqrt(distance);

    return distance;
}

/*!
* Evaluates the maximum distance among the specified point and the
* cells contained in the bounding box associated to the node.
*
* \param point is the point
* \result The maximum distance among the specified point and the cells
* contained in the bounding box associated to the node.
*/
double SkdNode::evalPointMaxDistance(const std::array<double, 3> &point) const
{
    double distance = 0.;
    for (int d = 0; d < 3; ++d) {
        distance += std::pow(std::max(point[d] - m_boxMin[d], m_boxMax[d] - point[d]), 2);
    }
    distance = std::sqrt(distance);

    return distance;
}

/*!
* Computes the distance between the specified point and the closest
* cell contained in the bounding box associated to the node.
*
* \param point is the point
* \result The distance between the specified point and the closest
* cell contained in the bounding box associated to the node.
*/
double SkdNode::evalPointDistance(const std::array<double, 3> &point) const
{
    long id;
    double distance;

    findPointClosestCell(point, &id, &distance);

    return distance;
}

/*!
* Given the specified point find the closest cell contained in the
* bounding box associated to the node and evaluates the distance
* between that cell and the given point.
*
* \param point is the point
* \param[out] id on output it will contain the id of the closest cell
* \param[out] distance on output it will contain the distance between
* the point and the closest cell
*/
void SkdNode::findPointClosestCell(const std::array<double, 3> &point,
                                   long *id, double *distance) const
{
    *id       = Cell::NULL_ID;
    *distance = std::numeric_limits<double>::max();
    if (getCellCount() == 0) {
        return;
    }

    const PatchKernel &patch = m_patchInfo->getPatch();
    const PiercedVector<Cell> &cells = patch.getCells();
    const std::vector<std::size_t> &cellRawIds = m_patchInfo->getCellRawIds();
    std::vector<std::array<double, 3>> cellVertexCoordinates(ReferenceElementInfo::MAX_ELEM_VERTICES);

    for (std::size_t n = m_cellRangeBegin; n < m_cellRangeEnd; n++) {
        std::size_t cellRawId = cellRawIds[n];
        const Cell &cell = cells.rawAt(cellRawId);

        // Get the vertices ids
        ConstProxyVector<long> elementVertexIds = cell.getVertexIds();
        const int nElementVertices = elementVertexIds.size();

        // Get vertex coordinates
        cellVertexCoordinates.resize(nElementVertices);
        for (int i = 0; i < nElementVertices; ++i) {
            cellVertexCoordinates[i] = patch.getVertex(elementVertexIds[i]).getCoords();
        }

        // Evaluate the distance from the cell
        double cellDistance = cell.evalPointDistance(point, cellVertexCoordinates.data());

        // Update the distance
        if (cellDistance < *distance) {
            *distance = cellDistance;
            *id       = cell.getId();
        }
    }
}

/*!
* Checks if the specified point is inside the bounding box associated to the
* node. The bounding box size will be expanded by the specified offset value.
*
* \param point is the point
* \param offset is the offset that will be used to expand the bounding
* box
* \result Returns true if the point is inside the inflated bounding box,
* false otherwise.
*/
bool SkdNode::boxContainsPoint(const std::array<double, 3> &point, double offset) const
{
    for (int d = 0; d < 3; d++) {
        if (point[d] < (m_boxMin[d] - offset)) {
            return false;
        }

        if (point[d] > (m_boxMax[d] + offset)) {
            return false;
        }
    }

    return true;
}

/*!
* Checks if the bounding box associated to the node intersects the
* sphere with given center and radius.
*
* \param[in] center is the center of the sphere
* \param[in] radius is the radius of the sphere
* \result Returns true if the bounding box associated to the node
* intersects the sphere with given center and radius, false otherwise.
*/
bool SkdNode::boxIntersectsSphere(const std::array<double, 3> &center, double radius) const
{
    // If the box contains the center it will also intersect the sphere
    if (boxContainsPoint(center, 0.)) {
        return true;
    }

    // Check if the distance between the center and its projection on
    // the bounding box is smaller than the radius of the sphere.
    std::array<double, 3> delta;
    for (int d = 0; d < 3; d++) {
        delta[d] = std::min(std::max(center[d], m_boxMin[d]), m_boxMax[d]) - center[d];
    }
    double distance = dotProduct(delta, delta);

    return (distance < radius * radius);
}

/*!
* \class PatchSkdTree
*
* \brief PatchSkdTree is the class that implements a spatial kd-tree
* (skd-tree) a bitpit patch.
*
* Spatial kd-tree is a spatial access method presented bt Ooi et al.
* where successive levels are split along different dimensions.
* Tke skd-tree was proposed to handle spatial obects with extension,
* which could not be handled by the original kd-tree. Skd-trees belong
* to the Bounding Volume Hierarchy (BVH) class of geometric search
* structures. Bounding volume hierarchies come in many variations:
* swept spheres, OBB-trees, sphere trees, and skd-trees are all
* examples.
*
* Each node of a skd-tree structure contains a bounding volume for
* some subset of the boundary geometry. The nodes are generally
* arranged in an oriented tree, where child nodes bound non-empty
* subsets of their parents’ geometry. The bounding volumes are
* selected to minimize the cost of query operations (e.g. proximity,
* intersection, or containment) while providing a close fit to the
* underlying geometry.
*
* Analyzing the asymptotic behavior of skd-trees is difficult. In
* the best case, skd-tree queries can be answered in constant time.
* In the worst case, each node of the tree may have to be visited,
* leading to O(n) work for the example of a binary tree with n leaf
* nodes. On average, the effort for intersection testing queries
* appears to be O(log2 n), while nearest features queries appear to
* be considerably cheaper than O(n).
*
* See:
*
* Deterministic point inclusion methods for computational applications
* with complex geometry, Ahmed Khamayseh and Andrew Kuprat.
* Computational Science & Discovery, Volume 1, Number 1, 2008.
* [This paper contains the algorithms for building the skd-tree and
* for performing closest cell queries.]
*
* Efficient Query Processing in Geographic Information Systems (Lecture
* Notes in Computer Science vol 471), Ooi B. C., 1990.
*
* A comparative study of spatial indexing techniques for
* multidimensional scientific datasets, Beomseok Nam and Alan Sussman
* Nam, International Conference on Scientific and Statistical Database
* Management, 2004.
*/

/*!
* Constructor.
*
* \param patch is the patch that will be use to build the tree
*/
PatchSkdTree::PatchSkdTree(const PatchKernel *patch)
    : m_patchInfo(patch, &m_cellRawIds), m_cellRawIds(patch->getInternalCount()),
      m_leafCapacity(0), m_nLeafs(0)
{
}

/*!
* Gets the maximum number of elements that can be contained in the leafs of
* the tree.
*
* \result The maximum number of elements that can be contained in the leafs
* of the tree.
*/
int PatchSkdTree::getLeafCapacity() const
{
    return m_leafCapacity;
}

/*!
* Sets the maximum number of elements that can be contained in the leafs of
* the tree.
*
* \param capacity is the maximum number of elements that can be contained in
* the leafs of the tree
*/
void PatchSkdTree::setLeafCapacity(int capacity)
{
    if (capacity <= 0) {
        throw std::runtime_error("Leaf capacity should be greater than zero.");
    }

    m_leafCapacity = capacity;
}

/*!
* Build the tree.
*/
void PatchSkdTree::build(int leafCapacity)
{
    const PatchKernel &patch = m_patchInfo.getPatch();

    // Clear existing tree
    clear();

    // Set leaf size
    setLeafCapacity(leafCapacity);

    // Initialize list of cell raw ids
    std::size_t nCells = patch.getInternalCount();
    if (m_cellRawIds.size() != nCells) {
        m_cellRawIds.resize(nCells);
    }

    std::size_t k = 0;
    for (auto itr = patch.internalConstBegin(); itr != patch.internalConstEnd(); ++itr) {
        m_cellRawIds[k++] = itr.getRawIndex();
    }

    // Build patch cache
    m_patchInfo.buildCache();

    // Initialize node list
    m_nodes.reserve(1.5 * nCells / getLeafCapacity());

    // Create the root
    m_nodes.emplace_back(&m_patchInfo, 0, nCells);

    // Create the tree
    std::queue<std::size_t> nodeStack;
    nodeStack.push(0);
    while (!nodeStack.empty()) {
        std::size_t nodeId = nodeStack.front();
        nodeStack.pop();

        // Create the children of the node
        createChildren(nodeId);

        // Add the newly created childrend to the stack
        const SkdNode &node = getNode(nodeId);
        for (int i = SkdNode::CHILD_BEGIN; i != SkdNode::CHILD_END; ++i) {
            std::size_t childId = node.getChildId(static_cast<SkdNode::ChildLocation>(i));
            if (childId != SkdNode::NULL_ID) {
                nodeStack.push(childId);
            }
        }
    }

    // Patch cache is no longer needed
    m_patchInfo.destroyCache();
}

/*!
* Clear the tree.
*
* \param release if it's true the memory hold by the tree will be released,
* otherwise the treewill be cleared but its memory will not be relased
*/
void PatchSkdTree::clear(bool release)
{
    m_nLeafs       = 0;
    m_leafCapacity = 0;

    if (release) {
        std::vector<SkdNode, SkdNode::Allocator>().swap(m_nodes);
        std::vector<std::size_t>().swap(m_cellRawIds);
    } else {
        m_nodes.clear();
        m_cellRawIds.clear();
    }
}

/*!
* Get a constant reference to the patch associated to the tree.
*
* \result A a constant reference to the patch associated to the tree.
*/
const PatchKernel & PatchSkdTree::getPatch() const
{
      return m_patchInfo.getPatch();
}

/*!
* Get the number of nodes contained in the tree.
*
* \result The number of nodes contained in the tree.
*/
std::size_t PatchSkdTree::getNodeCount() const
{
      return m_nodes.size();
}

/*!
* Get the number of leafs contained in the tree.
*
* \result The number of nodes contained in the tree.
*/
std::size_t PatchSkdTree::getLeafCount() const
{
      return m_nLeafs;
}

/*!
* Get a constant reference to the specified node.
*
* \param nodeId is the id of the node
* \result A constant reference to the specified node.
*/
const SkdNode & PatchSkdTree::getNode(std::size_t nodeId) const
{
      return m_nodes[nodeId];
}

/*!
* Get a reference to the specified node.
*
* \param nodeId is the id of the node
* \result A reference to the specified node.
*/
SkdNode & PatchSkdTree::_getNode(std::size_t nodeId)
{
      return m_nodes[nodeId];
}

/*!
* Evaluate the maximum depth of the tree.
*
* \param rootId is the id of the root node
* \result The maximum depth of the tree.
*/
std::size_t PatchSkdTree::evalMaxDepth(std::size_t rootId) const
{
    if (rootId == SkdNode::NULL_ID) {
        return 0;
    }

    const SkdNode &node = getNode(rootId);

    std::size_t depth = 0;
    for (int i = SkdNode::CHILD_BEGIN; i != SkdNode::CHILD_END; ++i) {
        SkdNode::ChildLocation childLocation = static_cast<SkdNode::ChildLocation>(i);
        depth = std::max(evalMaxDepth(node.getChildId(childLocation)), depth);
    }
    ++depth;

    return depth;
}

/*!
* Create the children of the specified node.
*
* \param parentId is the index of the parent node
*/
void PatchSkdTree::createChildren(std::size_t parentId)
{
    SkdNode &parent = _getNode(parentId);

    // Check if the parent is a leaf.
    if (parent.getCellCount() <= m_leafCapacity) {
        ++m_nLeafs;
        return;
    }

    // Evaluate the preferred direction along which elements will be split.
    //
    // The elements will be split along a plane normal to the direction
    // for which the bounding box has the maximum length.
    const std::array<double, 3> &parentBoxMin = parent.getBoxMin();
    const std::array<double, 3> &parentBoxMax = parent.getBoxMax();

    int largerDirection = 0;
    double boxMaximumLength = parentBoxMax[largerDirection] - parentBoxMin[largerDirection];
    for (int d = 1; d < 3; ++d) {
        double length = parentBoxMax[d] - parentBoxMin[d];
        if (length > boxMaximumLength) {
            largerDirection  = d;
            boxMaximumLength = length;
        }
    }

    // Split the elements.
    //
    // Preferred direction is tried first, if the split along this direction
    // will produce and empty child (e.g., all the cell's box centroids have
    // the same coordinate along the split direction and therefore all the
    // cells would be assigned to the left chiled) the split will be performed
    // along one of the other directions.
    std::array<double, 3> parentWeightedCentroid = parent.evalBoxWeightedMean();
    for (int d = 0; d < 3; ++d) {
        // Update the split direction
        int splitDirection = (largerDirection + d) % 3;

        // Get the threshold for the split
        double splitThreshold = parentWeightedCentroid[splitDirection];

        // Order the elements
        //
        // All the elements with a centroid coordinate less or equal than the
        // threshold will be assigned to the left child, the others will be
        // assigned to the right child.
        std::size_t leftBegin  = parent.m_cellRangeBegin;
        std::size_t leftEnd    = parent.m_cellRangeEnd;
        std::size_t rightBegin = parent.m_cellRangeBegin;
        std::size_t rightEnd   = parent.m_cellRangeEnd;
        while (true) {
            // Update the right begin
            while (rightBegin != leftEnd && m_patchInfo.getCachedCentroid(m_cellRawIds[rightBegin])[splitDirection] <= splitThreshold) {
                rightBegin++;
            }

            // Update the left end
            while (rightBegin != leftEnd && m_patchInfo.getCachedCentroid(m_cellRawIds[leftEnd - 1])[splitDirection] > splitThreshold) {
                leftEnd--;
            }

            // If all the elements are in the right position we can exit
            if (rightBegin == leftEnd) {
                break;
            }

            // If left end and and right begin are not equal, that the two ids
            // point to misplaced elements. Swap the elements, advance the ids
            // and continue iterating.
            std::iter_swap(m_cellRawIds.begin() + (leftEnd - 1), m_cellRawIds.begin() + rightBegin);
        }

        if (leftEnd <= leftBegin || rightEnd <= rightBegin) {
            continue;
        }

        // Create the left child
        long leftId = m_nodes.size();
        m_nodes.emplace_back(&m_patchInfo, leftBegin, leftEnd);
        m_nodes[parentId].m_children[static_cast<std::size_t>(SkdNode::CHILD_LEFT)] = leftId;

        // Create the right child
        long rightId = m_nodes.size();
        m_nodes.emplace_back(&m_patchInfo, rightBegin, rightEnd);
        m_nodes[parentId].m_children[static_cast<std::size_t>(SkdNode::CHILD_RIGHT)] = rightId;

        // Done
        return;
    }

    throw std::runtime_error("Unable to create the children.");
}

}
