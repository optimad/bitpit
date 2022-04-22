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

# ifndef __BITPIT_LEVELSET_OCTREE_KERNEL_HPP__
# define __BITPIT_LEVELSET_OCTREE_KERNEL_HPP__

#include "levelSetKernel.hpp"

#include "bitpit_voloctree.hpp"

namespace bitpit{

class LevelSetExternalPiercedStorageManager;
class LevelSetInternalPiercedStorageManager;

class LevelSetOctreeKernel : public LevelSetKernel{

    private:
    std::vector<double>                         m_levelToCellIncircle ;        /**< Incircles associated with cell levels*/
    std::vector<double>                         m_levelToCellCircumcircle ;    /**< Circumcircles associated with cell levels*/

    void                                        clearCellCirclesCache();
    void                                        updateCellCirclesCache();

    public:
    typedef LevelSetExternalPiercedStorageManager DenseStorageManager;
    typedef LevelSetInternalPiercedStorageManager SparseStorageManager;

    LevelSetOctreeKernel( VolOctree & );

    VolOctree *                                 getMesh() const override;

    double                                      computeCellIncircle(long) const override;
    double                                      computeCellCircumcircle(long) const override;

    void                                        clearGeometryCache() override;
    void                                        updateGeometryCache(const std::vector<adaption::Info> &) override;

    bool                                        intersectCellPlane(long, const std::array<double,3> &, const std::array<double,3> &, double) override;

};

// Typdefs for compatibility with older versions
typedef LevelSetOctreeKernel LevelSetOctree;

}

#endif
