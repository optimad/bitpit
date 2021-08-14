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

# ifndef __BITPIT_LEVELSET_OCTREE_HPP__
# define __BITPIT_LEVELSET_OCTREE_HPP__

namespace bitpit{

class VolOctree;
class LevelSetKernel ;

class LevelSetOctree : public LevelSetKernel{

    private:
    VolOctree*                                  m_octree ;       /**< Pointer to underlying octree mesh*/

    std::vector<double>                         m_levelToCellIncircle ;        /**< Incircles associated with cell levels*/
    std::vector<double>                         m_levelToCellCircumcircle ;    /**< Circumcircles associated with cell levels*/

    void                                        clearCellCirclesCache();
    void                                        updateCellCirclesCache();

    public:
    LevelSetOctree( VolOctree & );

    VolOctree *                                 getOctreeMesh() const;

    double                                      computeCellIncircle(long) const override;
    double                                      computeCellCircumcircle(long) const override;

    void                                        clearGeometryCache() override;
    void                                        updateGeometryCache(const std::vector<adaption::Info> &) override;

    bool                                        intersectCellPlane(long, const std::array<double,3> &, const std::array<double,3> &, double) override;

};

}

#endif
