/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2022 OPTIMAD engineering Srl
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

# ifndef __BITPIT_LEVELSET_UNSTRUCTURED_KERNEL_HPP__
# define __BITPIT_LEVELSET_UNSTRUCTURED_KERNEL_HPP__

#include "levelSetKernel.hpp"

# include "bitpit_volunstructured.hpp"

namespace bitpit{

class LevelSetExternalPiercedStorageManager;
class LevelSetInternalPiercedStorageManager;

struct LevelSetUnstructuredCellCacheEntry {

    std::array<double, 3> centroid;

    double tangentRadius;
    double boundingRadius;

    LevelSetUnstructuredCellCacheEntry(const VolumeKernel &patch, long cellId);

};

class LevelSetUnstructuredKernel : public LevelSetCachedKernel<LevelSetUnstructuredCellCacheEntry> {

    public:
    typedef LevelSetExternalPiercedStorageManager DenseStorageManager;
    typedef LevelSetInternalPiercedStorageManager SparseStorageManager;

    LevelSetUnstructuredKernel( VolUnstructured & );

    VolUnstructured *                           getMesh() const;

    std::array<double, 3>                       computeCellCentroid(long) const override;
    double                                      computeCellTangentRadius(long) const override;
    double                                      computeCellBoundingRadius(long) const override;

};

// Typdefs for compatibility with older versions
typedef LevelSetUnstructuredKernel LevelSetUnstructured;

}

#endif
