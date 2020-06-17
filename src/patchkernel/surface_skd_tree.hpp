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

# ifndef __BITPIT_SURFACE_SKD_TREE_HPP__
# define __BITPIT_SURFACE_SKD_TREE_HPP__

#include "patch_skd_tree.hpp"
#include "surface_kernel.hpp"

namespace bitpit {

class SurfaceSkdTree : public PatchSkdTree {

public:
#if BITPIT_ENABLE_MPI
    /*!
     * Info structure with distance, rank owner and cell id attributes.
     * The structure is used to communicate info about a cell between processes.
     */
    struct SkdCellInfo {
        double distance = 1.0e+18;
        int rank = -1;
        long id = Cell::NULL_ID;

        SkdCellInfo(){};
        SkdCellInfo(double distance_, int rank_, long id_):distance(distance_),rank(rank_),id(id_){};
    };
#endif

    SurfaceSkdTree(const SurfaceKernel *patch, bool includeGhosts = true);

    void clear(bool release);

    double evalPointDistance(const std::array<double,3> &point) const;
    double evalPointDistance(const std::array<double,3> &point, double maxDistance) const;

    long findPointClosestCell(const std::array<double,3> &point, long *id, double *distance) const;
    long findPointClosestCell(const std::array<double, 3> &point, double maxDistance, long *id, double *distance) const;

#if BITPIT_ENABLE_MPI
    void evalPointGlobalDistance(const std::size_t nPoints, const std::array<double, 3> *points, double *distances) const;
    void evalPointGlobalDistance(const std::size_t nPoints, const std::array<double, 3> *points, double maxDistance, double *distances) const;
    void findPointClosestGlobalCell(const std::size_t nPoints, const std::array<double, 3> *points, long *ids, int *ranks, double *distances) const;
    void findPointClosestGlobalCell(const std::size_t nPoints, const std::array<double, 3> *points, double maxDistance, long *ids, int *ranks, double *distances) const;
#endif

private:
    mutable std::vector<std::size_t> m_candidateIds;
    mutable std::vector<double> m_candidateMinDistances;

};

#if BITPIT_ENABLE_MPI
void minSkdCellInfo(SurfaceSkdTree::SkdCellInfo * in, SurfaceSkdTree::SkdCellInfo * inout, int * len, MPI_Datatype * datatype);
#endif

}

#endif
