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

# ifndef __BITPIT_LEVELSET_MASK_OBJECT_HPP__
# define __BITPIT_LEVELSET_MASK_OBJECT_HPP__

// Standard Template Library
# include <vector>
# include <set>

# include "bitpit_operators.hpp"
# include "bitpit_surfunstructured.hpp"

# include "levelSetSegmentationObject.hpp"

namespace bitpit{

class PatchKernel;
class SurfUnstructured;

template<typename storage_manager_t>
using LevelSetMaskNarrowBandCache = LevelSetSegmentationNarrowBandCache<storage_manager_t>;

template<typename storage_manager_t>
class LevelSetMaskObject : public LevelSetSegmentationObject<storage_manager_t> {

    private:
    std::unique_ptr<SurfUnstructured> extractCellEnvelope(const std::unordered_set<long> &, const VolumeKernel &, std::unordered_map<long,long> &);
    std::unique_ptr<SurfUnstructured> extractFaceEnvelope(const std::vector<long> &, const VolumeKernel &, std::unordered_map<long,long> &);
    bool sameInterfaceEnvelopeOrientation(const VolumeKernel &, long, SurfUnstructured &, long );


    public:
    LevelSetMaskObject(int, const std::unordered_set<long> &, const VolumeKernel &);
    LevelSetMaskObject(int, const std::vector<long> &, long, bool, const VolumeKernel &);
};

// Typdefs for compatibility with older versions
typedef LevelSetMaskObject<LevelSetInternalPiercedStorageManager> LevelSetMask;

}

// Include template implementations
#include "levelSetMaskObject.tpp"

// Explicit instantization
#ifndef __BITPIT_LEVELSET_MASK_OBJECT_SRC__
namespace bitpit {

extern template class LevelSetMaskObject<LevelSetMaskNarrowBandCache<LevelSetExternalPiercedStorageManager>>;
extern template class LevelSetMaskObject<LevelSetMaskNarrowBandCache<LevelSetInternalPiercedStorageManager>>;

}
#endif

#endif
