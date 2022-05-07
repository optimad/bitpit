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

#include "bitpit_common.hpp"

#include "surface_skd_tree.hpp"

namespace bitpit {

/*!
* \class SurfaceSkdTree
*
* \brief The SurfaceSkdTree implements a Bounding Volume Hierarchy tree for
* surface patches.
*/

/*!
* Constructor.
*
* \param patch is the surface patch that will be use to build the tree
* \param interiorCellsOnly if set to true, only interior cells will be considered
*/
SurfaceSkdTree::SurfaceSkdTree(const SurfaceKernel *patch, bool interiorCellsOnly)
    : PatchSkdTree(patch, interiorCellsOnly)
{
}

/*!
* Clear the tree.
*
* \param release if it's true the memory hold by the tree will be released,
* otherwise the treewill be cleared but its memory will not be relased
*/
void SurfaceSkdTree::clear(bool release)
{
    if (release) {
        std::vector<std::size_t>().swap(m_nodeStack);
        std::vector<std::size_t>().swap(m_candidateIds);
        std::vector<double>().swap(m_candidateMinDistances);
    }

    PatchSkdTree::clear(release);
}

/*!
* Computes the distance between the specified point and the closest
* cell contained in the tree. Only cells with a distance less than
* the specified maximum distance will be considered.
*
* \param[in] point is the point
* \result The distance between the specified point and the closest
* cell contained in the tree.
*/
double SurfaceSkdTree::evalPointDistance(const std::array<double, 3> &point) const
{
    return evalPointDistance(point, std::numeric_limits<double>::max(), false);
}

/*!
* Computes the distance between the specified point and the closest
* cell contained in the tree. Only cells with a distance less than
* the specified maximum distance will be considered. If all cells
* contained in the tree are farther than the maximum distance, the
* function will return the maximum representable distance.
*
* \param[in] point is the point
* \param[in] maxDistance all cells whose distance is greater than
* this parameters will not be considered for the evaluation of the
* distance
* \result The distance between the specified point and the closest
* cell contained in the tree. If all cells contained in the tree are
* farther than the maximum distance, the function will return the
* maximum representable distance.
*/
double SurfaceSkdTree::evalPointDistance(const std::array<double, 3> &point, double maxDistance) const
{
    return evalPointDistance(point, maxDistance, false);
}

/*!
* Computes the distance between the specified point and the closest
* cell contained in the tree. Only cells with a distance less than
* the specified maximum distance will be considered. If all cells
* contained in the tree are farther than the maximum distance, the
* function will return the maximum representable distance.
*
* \param[in] point is the point
* \param[in] maxDistance all cells whose distance is greater than
* this parameters will not be considered for the evaluation of the
* distance
* \param[in] interiorCellsOnly if set to true, only interior cells will be considered,
* it will be possible to consider non-interior cells only if the tree has been
* instantiated with non-interior cells support enabled
* \result The distance between the specified point and the closest
* cell contained in the tree. If all cells contained in the tree are
* farther than the maximum distance, the function will return the
* maximum representable distance.
*/
double SurfaceSkdTree::evalPointDistance(const std::array<double, 3> &point, double maxDistance, bool interiorCellsOnly) const
{
    long id;
    double distance = maxDistance;

    findPointClosestCell(point, interiorCellsOnly, &id, &distance);

    return distance;
}

/*!
* Computes the distance between each of the the specified points and the
* closest cell contained in the tree. Only cells with a distance less than
* the specified maximum distance will be considered.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[out] distances on output it will contain the distances
* between the points and closest cells
*/
void SurfaceSkdTree::evalPointDistance(int nPoints, const std::array<double, 3> *points, double *distances) const
{
    std::vector<double> maxDistances(nPoints, std::numeric_limits<double>::max());

    evalPointDistance(nPoints, points, maxDistances.data(), false, distances);
}

/*!
* Computes the distance between each of the the specified points and the
* closest cell contained in the tree. Only cells with a distance less than
* the specified maximum distance will be considered. If all cells contained
* in the tree are farther than the maximum distance, the function will return
* the maximum representable distance.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[in] maxDistance all cells whose distance is greater than
* this parameters will not be considered for the evaluation of the
* distance
* \param[out] distances on output it will contain the distances
* between the points and closest cells. If all cells contained in the tree are
* farther than the maximum distance, the related argument will be set to the
* maximum representable distance.
*/
void SurfaceSkdTree::evalPointDistance(int nPoints, const std::array<double, 3> *points, double maxDistance, double *distances) const
{
    std::vector<double> maxDistances(nPoints, maxDistance);

    evalPointDistance(nPoints, points, maxDistances.data(), false, distances);
}

/*!
* Computes the distance between the specified points and the closest
* cell contained in the tree. Only cells with a distance less than
* the specified maximum distance will be considered. If all cells
* contained in the tree are farther than the maximum distance, the
* function will return the maximum representable distance.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[in] maxDistances are the maximum allowed distances, all cells whose
* distance is greater than this parameter will not be considered for the
* evaluation of the distance with respect to the related point
* \param[out] distances on output it will contain the distances
* between the points and closest cells. If all cells contained in the tree are
* farther than the maximum distance, the related argument will be set to the
* maximum representable distance.
*/
void SurfaceSkdTree::evalPointDistance(int nPoints, const std::array<double, 3> *points, const double *maxDistances, double *distances) const
{
    evalPointDistance(nPoints, points, maxDistances, false, distances);
}


/*!
* Computes the distance between the specified points and the closest
* cell contained in the tree. Only cells with a distance less than
* the specified maximum distance will be considered. If all cells
* contained in the tree are farther than the maximum distance, the
* function will return the maximum representable distance.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[in] maxDistances are the maximum allowed distances, all cells whose
* distance is greater than this parameter will not be considered for the
* evaluation of the distance with respect to the related point
* \param[in] interiorCellsOnly if set to true, only interior cells will be
* considered
* \param[out] distances on output it will contain the distances
* between the points and closest cells. If all cells contained in the tree are
* farther than the maximum distance, the related argument will be set to the
* maximum representable distance.
*/
void SurfaceSkdTree::evalPointDistance(int nPoints, const std::array<double, 3> *points, const double *maxDistances, bool interiorCellsOnly, double *distances) const
{
    std::vector<long> ids(nPoints, Cell::NULL_ID);

    findPointClosestCell(nPoints, points, maxDistances, interiorCellsOnly, ids.data(), distances);
}

#if BITPIT_ENABLE_MPI
/*!
* Computes the distance between the specified point, considered distributed
* on the processes, and the closest cells contained in the tree.
* Only cells with a distance less than the specified maximum distance will
* be considered.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[out] distances on output it will contain the distances
* between the points and closest cells
*/
void SurfaceSkdTree::evalPointGlobalDistance(int nPoints, const std::array<double, 3> *points, double *distances) const
{
    std::vector<double> maxDistances(nPoints, std::numeric_limits<double>::max());

    evalPointGlobalDistance(nPoints, points, maxDistances.data(), distances);
}

/*!
* Computes the distance between the specified point, considered distributed
* on the processes, and the closest cells contained in the tree.
* Only cells with a distance less than the specified maximum distance will
* be considered. If all cells contained in the tree are farther
* than the maximum distance, the function will return the maximum representable
* distance.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[in] maxDistance all cells whose distance is greater than
* this parameters will not be considered for the evaluation of the
* distance
* \param[out] distances on output it will contain the distances
* between the points and closest cells. If all cells contained in the tree are
* farther than the maximum distance, the related argument will be set to the
* maximum representable distance.
*/
void SurfaceSkdTree::evalPointGlobalDistance(int nPoints, const std::array<double, 3> *points, double maxDistance, double *distances) const
{
    std::vector<double> maxDistances(nPoints, maxDistance);

    evalPointGlobalDistance(nPoints, points, maxDistances.data(), distances);
}

/*!
* Computes the distance between the specified point, considered distributed
* on the processes, and the closest cells contained in the tree.
* Only cells with a distance less than the specified maximum distance will
* be considered. If all cells contained in the tree are farther
* than the maximum distance, the function will return the maximum representable
* distance.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[in] maxDistances are the maximum allowed distances, all cells whose
* distance is greater than this parameter will not be considered for the
* evaluation of the distance with respect to the related point
* \param[out] distances on output it will contain the distances
* between the points and closest cells. If all cells contained in the tree are
* farther than the maximum distance, the related argument will be set to the
* maximum representable distance.
*/
void SurfaceSkdTree::evalPointGlobalDistance(int nPoints, const std::array<double, 3> *points, const double *maxDistances, double *distances) const
{
    std::vector<long> ids(nPoints, Cell::NULL_ID);
    std::vector<int> ranks(nPoints, -1);

    findPointClosestGlobalCell(nPoints, points, maxDistances, ids.data(), ranks.data(), distances);
}
#endif

/*!
* Given the specified point find the closest cell contained in the
* three and evaluates the distance between that cell and the given
* point.
*
* \param[in] point is the point
* \param[out] id on output it will contain the id of the cell closest
* to the point, if the distance between the point and the surface is
* greater than the specified maximum distance, the id parameter will
* be set to the null id
* \param[in,out] distance on output it will contain the distance
* between the point and closest cell. If all cells contained in the
* tree are farther than the maximum distance, the function will return
* the maximum representable distance.
*/
long SurfaceSkdTree::findPointClosestCell(const std::array<double, 3> &point, long *id, double *distance) const
{
    return findPointClosestCell(point, std::numeric_limits<double>::max(), false, id, distance);
}

/*!
* Given the specified point find the closest cell contained in the
* three and evaluates the distance between that cell and the given
* point.
*
* \param[in] point is the point
* \param[in] maxDistance all cells whose distance is greater than
* this parameters will not be considered for the evaluation of the
* distance
* \param[out] id on output it will contain the id of the closest cell.
* If all cells contained in the tree are farther than the maximum
* distance, the argument will be set to the null id
* \param[out] distance on output it will contain the distance between
* the point and closest cell. If all cells contained in the tree are
* farther than the maximum distance, the argument will be set to the
* maximum representable distance
*/
long SurfaceSkdTree::findPointClosestCell(const std::array<double, 3> &point, double maxDistance,
                                          long *id, double *distance) const
{
    return findPointClosestCell(point, maxDistance, false, id, distance);
}

/*!
* Given the specified point find the closest cell contained in the tree and
* evaluate the distance between that cell and the given point.
*
* \param[in] point is the point
* \param[in] maxDistance all cells whose distance is greater than
* this parameters will not be considered for the evaluation of the
* distance
* \param[in] interiorCellsOnly if set to true, only interior cells will be considered,
* it will be possible to consider non-interior cells only if the tree has been
* instantiated with non-interior cells support enabled
* \param[out] id on output it will contain the id of the closest cell.
* If all cells contained in the tree are farther than the maximum
* distance, the argument will be set to the null id
* \param[out] distance on output it will contain the distance between
* the point and closest cell. If all cells contained in the tree are
* farther than the maximum distance, the argument will be set to the
* maximum representable distance
*/
long SurfaceSkdTree::findPointClosestCell(const std::array<double, 3> &point, double maxDistance,
                                          bool interiorCellsOnly, long *id, double *distance) const
{
    // Tolerance for distance evaluations
    const PatchKernel &patch = getPatch();
    double tolerance = patch.getTol();

    // Initialize the cell id
    *id = Cell::NULL_ID;

    // Get the root of the tree
    std::size_t rootId = 0;
    const SkdNode &root = m_nodes[rootId];
    if (root.isEmpty()) {
        *distance = std::numeric_limits<double>::max();

        return 0;
    }

    // Initialize a distance estimate
    //
    // The real distance will be lesser than or equal to the estimate.
    //
    // Care must be taken to avoid overflow when performing the multiplication.
    double squaredMaxDistance;
    if (maxDistance <= 1. || maxDistance < std::numeric_limits<double>::max() / maxDistance) {
        squaredMaxDistance = maxDistance * maxDistance;
    } else {
        squaredMaxDistance = std::numeric_limits<double>::max();
    }

    double squareDistanceEstimate = std::min(root.evalPointMaxSquareDistance(point), squaredMaxDistance);

    // Get a list of candidates nodes
    //
    // If threads safe lookups are not needed, some temporary data structures
    // are declared as member of the class to avoid their reallocation every
    // time the function is called.
    //
    // First, we gather all the candidates and then we evaluate the distance
    // of each candidate. Since distance estimate is constantly updated when
    // new nodes are processed, the final estimate may be smaller than the
    // minimum distance of some candidates. Processing the candidates after
    // scanning all the tree, allows to discard some of them without the need
    // of evaluating the exact distance.
    std::unique_ptr<std::vector<std::size_t>> privateNodeStack;
    std::unique_ptr<std::vector<std::size_t>> privateCandidateIds;
    std::unique_ptr<std::vector<double>> privateCandidateMinDistances;

    std::vector<std::size_t> *nodeStack;
    std::vector<std::size_t> *candidateIds;
    std::vector<double> *candidateMinDistances;
    if (areLookupsThreadSafe()) {
        privateNodeStack = std::unique_ptr<std::vector<std::size_t>>(new std::vector<std::size_t>());
        privateCandidateIds = std::unique_ptr<std::vector<std::size_t>>(new std::vector<std::size_t>());
        privateCandidateMinDistances = std::unique_ptr<std::vector<double>>(new std::vector<double>());

        nodeStack = privateNodeStack.get();
        candidateIds = privateCandidateIds.get();
        candidateMinDistances = privateCandidateMinDistances.get();
    } else {
        m_nodeStack.clear();
        m_candidateIds.clear();
        m_candidateMinDistances.clear();

        nodeStack = &m_nodeStack;
        candidateIds = &m_candidateIds;
        candidateMinDistances = &m_candidateMinDistances;
    }

    nodeStack->push_back(rootId);
    while (!nodeStack->empty()) {
        std::size_t nodeId = nodeStack->back();
        const SkdNode &node = m_nodes[nodeId];
        nodeStack->pop_back();

        // Do not consider nodes with a minimum distance greater than
        // the distance estimate
        double nodeMinSquareDistance = node.evalPointMinSquareDistance(point);
        if (utils::DoubleFloatingGreater()(nodeMinSquareDistance, squareDistanceEstimate, tolerance, tolerance)) {
            continue;
        }

        // Update the distance estimate
        //
        // The real distance will be less than or equal to the estimate.
        double nodeMaxSquareDistance = node.evalPointMaxSquareDistance(point);
        squareDistanceEstimate = std::min(nodeMaxSquareDistance, squareDistanceEstimate);

        // If the node is a leaf add it to the candidates, otherwise add its
        // children to the stack.
        bool isLeaf = true;
        for (int i = SkdNode::CHILD_BEGIN; i != SkdNode::CHILD_END; ++i) {
            SkdNode::ChildLocation childLocation = static_cast<SkdNode::ChildLocation>(i);
            std::size_t childId = node.getChildId(childLocation);
            if (childId != SkdNode::NULL_ID) {
                isLeaf = false;
                nodeStack->push_back(childId);
            }
        }

        if (isLeaf) {
            candidateIds->push_back(nodeId);
            candidateMinDistances->push_back(std::sqrt(nodeMinSquareDistance));
        }
    }

    // Process the candidates and find the closest cell
    if (!candidateIds->empty()) {
        *distance = std::sqrt(squareDistanceEstimate);
    } else {
        *distance = std::numeric_limits<double>::max();
    }

    long nDistanceEvaluations = 0;
    for (std::size_t k = 0; k < candidateIds->size(); ++k) {
        // Do not consider nodes with a minimum distance greater than the
        // distance estimate
        if (utils::DoubleFloatingGreater()(candidateMinDistances->at(k), *distance, tolerance, tolerance)) {
            continue;
        }

        // Evaluate the distance
        std::size_t nodeId = candidateIds->at(k);
        const SkdNode &node = m_nodes[nodeId];

        node.updatePointClosestCell(point, interiorCellsOnly, id, distance);
        ++nDistanceEvaluations;
    }

    return nDistanceEvaluations;
}

/*!
* For each of the specified points find the closest cell contained in the
* three and evaluates the distance between that cell and the point.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[out] ids on output it will contain the ids of the cells closest
* to the local points
* \param[out] distances on output it will contain the distances
* between the points and closest cells
*/
long SurfaceSkdTree::findPointClosestCell(int nPoints, const std::array<double, 3> *points, long *ids, double *distances) const
{
    long nDistanceEvaluations = 0;
    for (int i = 0; i < nPoints; ++i) {
        nDistanceEvaluations += findPointClosestCell(points[i], std::numeric_limits<double>::max(), false, ids + i, distances + i);
    }

    return nDistanceEvaluations;
}

/*!
* For each of the specified points find the closest cell contained in the
* three and evaluates the distance between that cell and the point.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[in] maxDistance all cells whose distance is greater than this
* parameters will not be considered for the evaluation of the distance
* \param[out] ids on output it will contain the ids of the cells closest
* to the points. If all cells contained in the tree are farther from a point
* than the maximum distance, the related id will be set to the null id
* \param[out] distances on output it will contain the distances
* between the points and closest cells. If all cells contained in the tree are
* farther than the maximum distance, the related argument will be set to the
* maximum representable distance.
*/
long SurfaceSkdTree::findPointClosestCell(int nPoints, const std::array<double, 3> *points, double maxDistance, long *ids, double *distances) const
{
    long nDistanceEvaluations = 0;
    for (int i = 0; i < nPoints; ++i) {
        nDistanceEvaluations += findPointClosestCell(points[i], maxDistance, false, ids + i, distances + i);
    }

    return nDistanceEvaluations;
}

/*!
* For each of the specified points find the closest cell contained in the
* three and evaluates the distance between that cell and the point.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[in] maxDistances are the maximum allowed distances, all cells whose
* distance is greater than this parameter will not be considered for the
* evaluation of the distance with respect to the related point
* \param[out] ids on output it will contain the ids of the cells closest
* to the points. If all cells contained in the tree are farther from a point
* than the maximum distance, the related id will be set to the null id
* \param[out] distances on output it will contain the distances
* between the points and closest cells. If all cells contained in the tree are
* farther than the maximum distance, the related argument will be set to the
* maximum representable distance.
*/
long SurfaceSkdTree::findPointClosestCell(int nPoints, const std::array<double, 3> *points, const double *maxDistances, long *ids, double *distances) const
{
    return findPointClosestCell(nPoints, points, maxDistances, false, ids, distances);
}

/*!
* For each of the specified points find the closest cell contained in the
* three and evaluates the distance between that cell and the point.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[in] maxDistances are the maximum allowed distances, all cells whose
* distance is greater than this parameter will not be considered for the
* evaluation of the distance with respect to the related point
* \param[in] interiorCellsOnly if set to true, only interior cells will be
* considered
* \param[out] ids on output it will contain the ids of the cells closest
* to the points. If all cells contained in the tree are farther from a point
* than the maximum distance, the related id will be set to the null id
* \param[out] distances on output it will contain the distances
* between the points and closest cells. If all cells contained in the tree are
* farther than the maximum distance, the related argument will be set to the
* maximum representable distance.
*/
long SurfaceSkdTree::findPointClosestCell(int nPoints, const std::array<double, 3> *points, const double *maxDistances, bool interiorCellsOnly, long *ids, double *distances) const
{
    long nDistanceEvaluations = 0;
    for (int i = 0; i < nPoints; ++i) {
        nDistanceEvaluations += findPointClosestCell(points[i], maxDistances[i], interiorCellsOnly, ids + i, distances + i);
    }

    return nDistanceEvaluations;
}

#if BITPIT_ENABLE_MPI
/*!
* Given the specified points, considered distributed on the processes, find
* the closest cells contained in the tree and evaluates the distance values
* between those cells and the given points.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[out] ids on output it will contain the ids of the cells closest
* to the local points
* \param[out] ranks on output it will contain the rank indices of the processes
* owner of the cells closest to the points
* \param[out] distances on output it will contain the distances
* between the points and closest cells
*/
long SurfaceSkdTree::findPointClosestGlobalCell(int nPoints, const std::array<double, 3> *points,
                                                long *ids, int *ranks, double *distances) const
{
    return findPointClosestGlobalCell(nPoints, points, std::numeric_limits<double>::max(), ids, ranks, distances);
}

/*!
* Given the specified points, considered distributed on the processes, find the
* closest cells contained in the tree and evaluates the distance values
* between those cells and the given points.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[in] maxDistance all cells whose distance is greater than this
* parameters will not be considered for the evaluation of the distance
* \param[out] ids on output it will contain the ids of the cells closest
* to the points. If all cells contained in the tree are farther from a point
* than the maximum distance, the related id will be set to the null id
* \param[out] ranks on output it will contain the rank indices of the processes
* owner of the cells closest to the points
* \param[out] distances on output it will contain the distances
* between the points and closest cells. If all cells contained in the tree are
* farther than the maximum distance, the related argument will be set to the
* maximum representable distance.
*/
long SurfaceSkdTree::findPointClosestGlobalCell(int nPoints, const std::array<double, 3> *points, double maxDistance,
                                                long *ids, int *ranks, double *distances) const
{
    std::vector<double> maxDistances(nPoints, maxDistance);

    return findPointClosestGlobalCell(nPoints, points, maxDistances.data(), ids, ranks, distances);
}

/*!
* Given the specified points, considered distributed on the processes, find the
* closest cells contained in the tree and evaluates the distance values
* between those cells and the given points.
*
* \param[in] nPoints is the number of the points
* \param[in] points are the points coordinates
* \param[in] maxDistances are the maximum allowed distances, all cells whose
* distance is greater than this parameter will not be considered for the
* evaluation of the distance with respect to the related point
* \param[out] ids on output it will contain the ids of the cells closest
* to the points. If all cells contained in the tree are farther from a point
* than the maximum distance, the related id will be set to the null id
* \param[out] ranks on output it will contain the rank indices of the processes
* owner of the cells closest to the points
* \param[out] distances on output it will contain the distances
* between the points and closest cells. If all cells contained in the tree are
* farther than the maximum distance, the related argument will be set to the
* maximum representable distance.
*/
long SurfaceSkdTree::findPointClosestGlobalCell(int nPoints, const std::array<double, 3> *points, const double *maxDistances,
                                                long *ids, int *ranks, double *distances) const
{
    long nDistanceEvaluations = 0;

    // Early return is the patch is not partitioned
    const PatchKernel &patch = getPatch();
    if (!patch.isPartitioned()) {
        for (int i = 0; i < nPoints; ++i) {
            // Evaluate distance
            nDistanceEvaluations += findPointClosestCell(points[i], maxDistances[i], ids + i, distances + i);

            // The patch is not partitioned, all cells are local
            ranks[i] = patch.getRank();
        }

        return nDistanceEvaluations;
    }

    // Get MPI communicator
    if (!isCommunicatorSet()) {
        throw std::runtime_error("Skd-tree communicator has not been set.");
    }

    MPI_Comm communicator = getCommunicator();

    // Gather the number of points associated to each process
    std::vector<int> pointsCount(m_nProcessors);
    MPI_Allgather(&nPoints, 1, MPI_INT, pointsCount.data(), 1, MPI_INT, communicator);

    // Evaluate information for data communications
    std::vector<int> globalPointsDispls(m_nProcessors, 0);
    std::vector<int> globalPointsOffsets(m_nProcessors, 0);
    std::vector<int> globalPointsDataCount(m_nProcessors, 0);

    globalPointsDataCount[0] = 3 * pointsCount[0];
    for (int i = 1; i < m_nProcessors; ++i) {
        globalPointsDispls[i]     = globalPointsDispls[i - 1] + 3 * pointsCount[i - 1];
        globalPointsOffsets[i]    = globalPointsOffsets[i - 1] + pointsCount[i - 1];
        globalPointsDataCount[i]  = 3 * pointsCount[i];
    }

    int nGlobalPoints = globalPointsDispls.back() + pointsCount.back();

    // Gather point coordinates
    std::vector<std::array<double,3>> globalPoints(nGlobalPoints);
    int pointsDataCount = 3 * nPoints;
    MPI_Allgatherv(points, pointsDataCount, MPI_DOUBLE, globalPoints.data(),
                globalPointsDataCount.data(), globalPointsDispls.data(), MPI_DOUBLE, communicator);

    // Gather vector with all maximum distances
    std::vector<double> globalMaxDistances(nGlobalPoints);
    MPI_Allgatherv(maxDistances, nPoints, MPI_DOUBLE, globalMaxDistances.data(),
                    pointsCount.data(), globalPointsOffsets.data(), MPI_DOUBLE, communicator);

    // Initialize distance information
    std::vector<SkdGlobalCellDistance> globalCellDistances(nGlobalPoints);

    // Call local find point closest cell for each global point collected
    for (int i = 0; i < nGlobalPoints; ++i) {
        // Get point information
        const std::array<double, 3> &point = globalPoints[i];

        // Use a maximum distance for each point given by an estimation
        // based on partition bounding boxes. The distance will be lesser
        // than or equal to the point maximum distance.
        double pointMaxDistance = globalMaxDistances[i];
        for (int rank = 0; rank < m_nProcessors; ++rank) {
            pointMaxDistance = std::min(getPartitionBox(rank).evalPointMaxDistance(point, std::numeric_limits<double>::max()), pointMaxDistance);
        }

        // Get cell distance information
        SkdGlobalCellDistance &globalCellDistance = globalCellDistances[i];
        int &cellRank = globalCellDistance.getRank();
        long &cellId = globalCellDistance.getId();
        double &cellDistance = globalCellDistance.getDistance();

        // Evaluate local distance from the point
        bool interiorCellsOnly = true;
        nDistanceEvaluations += findPointClosestCell(point, pointMaxDistance, interiorCellsOnly, &cellId, &cellDistance);

        // Set cell rank
        if (cellId != Cell::NULL_ID) {
            cellRank = patch.getCellOwner(cellId);
        }
    }

    // Exchange distance information
    MPI_Datatype globalCellDistanceDatatype = SkdGlobalCellDistance::getMPIDatatype();
    MPI_Op globalCellDistanceMinOp = SkdGlobalCellDistance::getMPIMinOperation();
    for (int rank = 0; rank < m_nProcessors; ++rank) {
        SkdGlobalCellDistance *globalCellDistance = globalCellDistances.data() + globalPointsOffsets[rank];
        if (m_rank == rank) {
            MPI_Reduce(MPI_IN_PLACE, globalCellDistance, pointsCount[rank], globalCellDistanceDatatype, globalCellDistanceMinOp, rank, communicator);
        } else {
            MPI_Reduce(globalCellDistance, globalCellDistance, pointsCount[rank], globalCellDistanceDatatype, globalCellDistanceMinOp, rank, communicator);
        }
    }

    // Update output arguments
    for (int i = 0; i < nPoints; ++i) {
        int globalIndex = i + globalPointsOffsets[m_rank];
        SkdGlobalCellDistance &globalCellDistance = globalCellDistances[globalIndex];
        globalCellDistance.exportData(ranks + i, ids + i, distances + i);
    }

    return nDistanceEvaluations;
}
#endif

}
