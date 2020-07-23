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
* \param interiorOnly if set to true, only interior cells will be considered
*/
SurfaceSkdTree::SurfaceSkdTree(const SurfaceKernel *patch, bool interiorOnly)
    : PatchSkdTree(patch, interiorOnly)
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
* \param[in] interiorOnly if set to true, only interior cells will be considered,
* it will be possible to consider non-interior cells only if the tree has been
* instantiated with non-interior cells support enabled
* \result The distance between the specified point and the closest
* cell contained in the tree. If all cells contained in the tree are
* farther than the maximum distance, the function will return the
* maximum representable distance.
*/
double SurfaceSkdTree::evalPointDistance(const std::array<double, 3> &point, double maxDistance, bool interiorOnly) const
{
    long id;
    double distance = maxDistance;

    findPointClosestCell(point, interiorOnly, &id, &distance);

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
* \param[in] interiorOnly if set to true, only interior cells will be
* considered
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
* \param[in] interiorOnly if set to true, only interior cells will be
* considered
* \param[out] distances on output it will contain the distances
* between the points and closest cells. If all cells contained in the tree are
* farther than the maximum distance, the related argument will be set to the
* maximum representable distance.
*/
void SurfaceSkdTree::evalPointDistance(int nPoints, const std::array<double, 3> *points, const double *maxDistances, bool interiorOnly, double *distances) const
{
    std::vector<long> ids(nPoints, Cell::NULL_ID);

    findPointClosestCell(nPoints, points, maxDistances, interiorOnly, ids.data(), distances);
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
* Given the specified point find the closest cell contained in the
* three and evaluates the distance between that cell and the given
* point.
*
* \param[in] point is the point
* \param[in] maxDistance all cells whose distance is greater than
* this parameters will not be considered for the evaluation of the
* distance
* \param[in] interiorOnly if set to true, only interior cells will be considered,
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
                                          bool interiorOnly, long *id, double *distance) const
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

        node.updatePointClosestCell(point, interiorOnly, id, distance);
        ++nDistanceEvaluations;
    }

    // If no closest cell was found set the distance to the maximum
    // representable distance.
    if (*id == Cell::NULL_ID) {
        *distance = std::numeric_limits<double>::max();
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
* \param[in] interiorOnly if set to true, only interior cells will be
* considered
* \param[out] ids on output it will contain the ids of the cells closest
* to the points. If all cells contained in the tree are farther from a point
* than the maximum distance, the related id will be set to the null id
* \param[out] distances on output it will contain the distances
* between the points and closest cells. If all cells contained in the tree are
* farther than the maximum distance, the related argument will be set to the
* maximum representable distance.
*/
long SurfaceSkdTree::findPointClosestCell(int nPoints, const std::array<double, 3> *points, const double *maxDistances, bool interiorOnly, long *ids, double *distances) const
{
    long nDistanceEvaluations = 0;
    for (int i = 0; i < nPoints; ++i) {
        nDistanceEvaluations += findPointClosestCell(points[i], maxDistances[i], interiorOnly, ids + i, distances + i);
    }

    return nDistanceEvaluations;
}

}
