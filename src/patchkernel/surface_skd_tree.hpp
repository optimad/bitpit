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
    SurfaceSkdTree(const SurfaceKernel *patch, bool interorOnly = false);

    void clear(bool release);

    double evalPointDistance(const std::array<double,3> &point) const;
    double evalPointDistance(const std::array<double,3> &point, double maxDistance) const;
    double evalPointDistance(const std::array<double,3> &point, double maxDistance, bool interorOnly) const;
    void evalPointDistance(int nPoints, const std::array<double, 3> *points, double *distances) const;
    void evalPointDistance(int nPoints, const std::array<double, 3> *points, double maxDistance, double *distances) const;
    void evalPointDistance(int nPoints, const std::array<double, 3> *points, const double *maxDistances, double *distances) const;
    void evalPointDistance(int nPoints, const std::array<double, 3> *points, const double *maxDistances, bool interorOnly, double *distances) const;

    long findPointClosestCell(const std::array<double,3> &point, long *id, double *distance) const;
    long findPointClosestCell(const std::array<double, 3> &point, double maxDistance, long *id, double *distance) const;
    long findPointClosestCell(const std::array<double, 3> &point, double maxDistance, bool interorOnly, long *id, double *distance) const;
    long findPointClosestCell(int nPoints, const std::array<double, 3> *points, long *ids, double *distances) const;
    long findPointClosestCell(int nPoints, const std::array<double, 3> *points, double maxDistance, long *ids, double *distances) const;
    long findPointClosestCell(int nPoints, const std::array<double, 3> *points, const double *maxDistances, long *ids, double *distances) const;
    long findPointClosestCell(int nPoints, const std::array<double, 3> *points, const double *maxDistances, bool interorOnly, long *ids, double *distances) const;

private:
    mutable std::vector<std::size_t> m_candidateIds;
    mutable std::vector<double> m_candidateMinDistances;

};

}

#endif
