/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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
* \param includeGhosts if set to true (default value) the ghost cells are included in the tree
*/
SurfaceSkdTree::SurfaceSkdTree(const SurfaceKernel *patch, bool includeGhosts)
    : PatchSkdTree(patch,includeGhosts)
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
    return evalPointDistance(point, std::numeric_limits<double>::max());
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
    long id;
    double distance = maxDistance;

    findPointClosestCell(point, &id, &distance);

    return distance;
}

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
    return findPointClosestCell(point, std::numeric_limits<double>::max(), id, distance);
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
    // Initialize the cell id
    *id = Cell::NULL_ID;

    // Initialize the distance with an estimate
    //
    // The real distance will be lesser than or equal to the estimate.
    std::size_t rootId = 0;
    const SkdNode &root = m_nodes[rootId];
    *distance = std::min(root.evalPointMaxDistance(point), maxDistance);

    // Get a list of candidates nodes
    //
    // The list of candidates is a memeber of the class to avoid its
    // reallocation every time the function is called.
    m_candidateIds.clear();
    m_candidateMinDistances.clear();

    std::vector<std::size_t> nodeStack;
    nodeStack.push_back(rootId);
    while (!nodeStack.empty()) {
        std::size_t nodeId = nodeStack.back();
        const SkdNode &node = m_nodes[nodeId];
        nodeStack.pop_back();

        // Do not consider nodes with a minimum distance greater than
        // the distance estimate
        double nodeMinDistance = node.evalPointMinDistance(point);
        if (nodeMinDistance > *distance) {
            continue;
        }

        // Update the distance estimate
        //
        // The real distance will be lesser than or equal to the
        // estimate.
        double nodeMaxDistance = node.evalPointMaxDistance(point);
        *distance = std::min(nodeMaxDistance, *distance);

        // If the node is a leaf add it to the candidates, otherwise
        // add its children to the stack.
        bool isLeaf = true;
        for (int i = SkdNode::CHILD_BEGIN; i != SkdNode::CHILD_END; ++i) {
            SkdNode::ChildLocation childLocation = static_cast<SkdNode::ChildLocation>(i);
            std::size_t childId = node.getChildId(childLocation);
            if (childId != SkdNode::NULL_ID) {
                isLeaf = false;
                nodeStack.push_back(childId);
            }
        }

        if (isLeaf) {
            m_candidateIds.push_back(nodeId);
            m_candidateMinDistances.push_back(nodeMinDistance);
        }
    }

    // Process the candidates and find the closest cell
    long nDistanceEvaluations = 0;
    for (std::size_t k = 0; k < m_candidateIds.size(); ++k) {
        // Do not consider nodes with a minimum distance greater than
        // the distance estimate
        if (m_candidateMinDistances[k] > *distance) {
            continue;
        }

        // Evaluate the distance
        std::size_t nodeId = m_candidateIds[k];
        const SkdNode &node = m_nodes[nodeId];

        node.updatePointClosestCell(point, id, distance);
        ++nDistanceEvaluations;
    }

    // If no closest cell was found set the distance to the maximum
    // representable distance.
    if (*id == Cell::NULL_ID) {
        *distance = std::numeric_limits<double>::max();
    }

    return nDistanceEvaluations;
}

#if BITPIT_ENABLE_MPI
/*!
* Computes the distance between the specified point, considered distributed
* on the processes, and the closest cells contained in the tree.
* Only cells with a distance less than the specified maximum distance will
* be considered.
*
* \param[in] nPoints number of the points
* \param[in] points points coordinates
* \param[out] distances on output it will contain the distances
* between the points and closest cells
*/
void SurfaceSkdTree::evalPointGlobalDistance(const std::size_t nPoints, const std::array<double, 3> *points, double *distances) const
{
    evalPointGlobalDistance(nPoints, points, std::numeric_limits<double>::max(), distances);
}

/*!
* Computes the distance between the specified point, considered distributed
* on the processes, and the closest cells contained in the tree.
* Only cells with a distance less than the specified maximum distance will
* be considered. If all cells contained in the tree are farther
* than the maximum distance, the function will return the maximum representable
* distance.
*
* \param[in] nPoints number of the points
* \param[in] points points coordinates
* \param[in] maxDistance all cells whose distance is greater than
* this parameters will not be considered for the evaluation of the
* distance
* \param[out] distances on output it will contain the distances
* between the points and closest cells
* \param[out] distances on output it will contain the distances
* between the points and closest cells. If all cells contained in the tree are
* farther than the maximum distance, the related argument will be set to the
* maximum representable distance.
*/
void SurfaceSkdTree::evalPointGlobalDistance(const std::size_t nPoints, const std::array<double, 3> *points, double maxDistance, double *distances) const
{
    std::vector<long> ids(nPoints, Cell::NULL_ID);
    std::vector<int> ranks(nPoints, -1);

    findPointClosestGlobalCell(nPoints, points, maxDistance, ids.data(), ranks.data(), distances);
}

/*!
 * Minimum operation to reduce the found distances and related cell info (rank and id)
 */
void minSkdCellInfo(SurfaceSkdTree::SkdCellInfo * in, SurfaceSkdTree::SkdCellInfo * inout, int * len, MPI_Datatype * datatype)
{
    for (std::size_t i = 0; i < *len; i++){
        if (std::abs(in[i].distance) < std::abs(inout[i].distance)){
            inout[i] = in[i];
        }
    }
}

/*!
* Given the specified points, considered distributed on the processes, find the
* closest cells contained in the tree and evaluates the distance values
* between those cells and the given points.
*
* \param[in] nPoints number of the points
* \param[in] points points coordinates
* \param[out] ids on output it will contain the ids of the cells closest
* to the local points
* \param[out] ranks on output it will contain the rank indices of the processes
* owner of the cells closest to the points
* \param[out] distances on output it will contain the distances
* between the points and closest cells
*/
void SurfaceSkdTree::findPointClosestGlobalCell(const std::size_t nPoints, const std::array<double, 3> *points,
                                          long *ids, int *ranks, double *distances) const
{
    findPointClosestGlobalCell(nPoints, points, std::numeric_limits<double>::max(), ids, ranks, distances);
}

/*!
* Given the specified points, considered distributed on the processes, find the
* closest cells contained in the tree and evaluates the distance values
* between those cells and the given points.
*
* \param[in] nPoints number of the points
* \param[in] points points coordinates
* \param[in] maxDistance all cells whose distance is greater than
* this parameters will not be considered for the evaluation of the
* distance
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
void SurfaceSkdTree::findPointClosestGlobalCell(const std::size_t nPoints, const std::array<double, 3> *points,
                                          double maxDistance, long *ids, int *ranks, double *distances) const
{
    // Patch is partitioned call the parallel method
    if (getPatch().isPartitioned()){

        // Initialize the cell ids and ranks
        for (std::size_t i = 0; i < nPoints; i++){
            ids[i] = Cell::NULL_ID;
            ranks[i] = -1;
        }

        // Communicate all the points to all the processes
        // Gather number of points multiplied by number of coordinates,
        // i.e. number of data to communicate
        std::vector<int> pointsCount(m_nProcessors);
        MPI_Allgather(&nPoints, 1, MPI_INT, pointsCount.data(), 1, MPI_INT, m_communicator);

        // Recover points displacements, offset and number of points data per process
        std::vector<int> pointsDispls(m_nProcessors, 0);
        std::vector<int> pointsOffsets(m_nProcessors, 0);
        std::vector<int> pointsCountData(m_nProcessors, 0);
        pointsCountData[0] = pointsCount[0] * 3;
        for (int i = 1; i < m_nProcessors; ++i) {
            pointsDispls[i] = pointsDispls[i - 1] + pointsCount[i - 1] * 3;
            pointsOffsets[i] = pointsOffsets[i - 1] + pointsCount[i - 1];
            pointsCountData[i] = pointsCount[i] * 3;
        }

        // Sum to global number of points
        int nGlobalPoints = 0;
        for (int np : pointsCount){
            nGlobalPoints += np;
        }

        // Gather vector with all global points
        std::vector<std::array<double,3>> globalPoints(nGlobalPoints);
        std::size_t nPointsCountData = nPoints * 3;
        MPI_Allgatherv(points, nPointsCountData, MPI_DOUBLE, globalPoints.data(),
                pointsCountData.data(), pointsDispls.data(), MPI_DOUBLE, m_communicator);

        // Call local find point closest cell for each global point collected

        // Instantiate global container for distances, ids and ranks (SkdCellInfo)
        std::vector<SkdCellInfo> dri_data(nGlobalPoints, SkdCellInfo(std::numeric_limits<double>::max(), m_rank, Cell::NULL_ID));

        // Call local find point closest cell for each global point collected
        for (std::size_t ip = 0; ip < nGlobalPoints; ip++){

            const std::array<double,3> & point = globalPoints[ip];
            double & distance = dri_data[ip].distance;
            int & rank = dri_data[ip].rank;
            long & id = dri_data[ip].id;

            // Use a maximum distance for each point given by an estimation based on partition
            // bounding boxes. The distance will be lesser than or equal to the point maximum distance
            double pointMaxDistance = maxDistance;
            for (int irank = 0; irank < m_nProcessors; irank++){
                pointMaxDistance = std::min(getPartitionBox(irank).evalPointMaxDistance(point), pointMaxDistance);
            }

            // Call local find point closest cell with estimated distance used as maximum distance
            long nDistanceEvaluations = findPointClosestCell(point, pointMaxDistance, &id, &distance);

        }

        // Force distance to numeric limits maximum value for point projected on ghost cells
        // The desired result id is the local cell id on the ghost owner rank
        for (std::size_t ip = 0; ip < nGlobalPoints; ip++){
            long cellId = dri_data[ip].id;
            if(cellId != Cell::NULL_ID && !getPatch().getCell(cellId).isInterior()) {
                dri_data[ip].distance = std::numeric_limits<double>::max();
            }
        }

        // Communicate the computed distances of the distributed input points to all processes
        // and retain the indices of the rank owner and the id of the closest cell

        // Prepare MPI custom data type and Operation
        // The data are of MPI custom data type  with distance, ids and rank as member
        int blocklengths[3] = {1,1,1};
        MPI_Aint displacements[3] = {offsetof(SkdCellInfo, distance),
                offsetof(SkdCellInfo, rank),
                offsetof(SkdCellInfo, id)};
        MPI_Datatype types[3] = {MPI_DOUBLE, MPI_INT, MPI_LONG};
        MPI_Datatype MPI_DRI;
        MPI_Type_create_struct(3, blocklengths, displacements, types, &MPI_DRI);
        MPI_Type_commit(&MPI_DRI);

        MPI_Op MPI_MIN_DRI;
        MPI_Op_create((MPI_User_function *) minSkdCellInfo, false, &MPI_MIN_DRI);

        // Communicate the closest cells distances, ranks and ids by reduce with custom minimum distance operation
        // Reduce only the right portion of data to the right process
        for (int irank = 0; irank < m_nProcessors; irank++){
            SkdCellInfo *buffer = (dri_data.data() + pointsOffsets[irank]);
            if (m_rank == irank){
                MPI_Reduce(MPI_IN_PLACE, buffer, pointsCount[irank], MPI_DRI, MPI_MIN_DRI, irank, getCommunicator());
            } else {
                MPI_Reduce(buffer, buffer, pointsCount[irank], MPI_DRI, MPI_MIN_DRI, irank, getCommunicator());
            }
        }

        // Update distances, rank indices and cell ids
        for (std::size_t ip = 0; ip < nPoints; ip++){
            std::size_t globalIndex = ip + pointsOffsets[m_rank];
            distances[ip] = dri_data[globalIndex].distance;
            ranks[ip] = dri_data[globalIndex].rank;
            ids[ip] = dri_data[globalIndex].id;
        }

    } else {

        // Call the serial method and fill the rank output with -1
        for (std::size_t ip = 0; ip < nPoints; ip++){

            const std::array<double,3> & point = points[ip];
            double & distance = distances[ip];
            long & id = ids[ip];
            int & rank = ranks[ip];

            // Call serial find point closest cell with maximum allowed distance
            long nDistanceEvaluations = findPointClosestCell(point, maxDistance, &id, &distance);

            // Fix rank with the default value given by getRank method pf the patch
            rank = getPatch().getRank();

        }
    }
}
#endif

}
