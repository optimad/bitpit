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

#ifndef __BITPIT_LEVELSET_IMMUTABLE_OBJECT_HPP__
#define __BITPIT_LEVELSET_IMMUTABLE_OBJECT_HPP__

#include "levelSetObject.hpp"
#include "levelSetBooleanObject.hpp"
#include "levelSetCachedObject.hpp"
#include "levelSetSignedObject.hpp"

namespace bitpit {

class LevelSetImmutableObjectBase {

protected:
    LevelSetImmutableObjectBase() = default;

};

template<typename storage_manager_t>
using LevelSetImmutableNarrowBandCache = LevelSetNarrowBandCache<storage_manager_t>;

template<typename narrow_band_cache_t>
class LevelSetImmutableObject : public LevelSetObject, public LevelSetImmutableObjectBase, public LevelSetCachedObjectInterface<narrow_band_cache_t>, public LevelSetSignedObjectInterface
{

public:
    LevelSetImmutableObject(LevelSetCachedObject<narrow_band_cache_t> *object);
    LevelSetImmutableObject(LevelSetBooleanObject *object);

    LevelSetImmutableObject * clone() const override;

    double getValue(long id) const override;
    short getSign(long id) const override;
    std::array<double,3> getGradient(long id) const override;

    LevelSetInfo computeLevelSetInfo(const std::array<double,3> &coords) const override;

protected:
    std::shared_ptr<LevelSetSignStorage> createSignStorage() override;

    LevelSetImmutableObject(int);

    void _clear() override;

    void _dump(std::ostream &stream) override;
    void _restore(std::istream &stream) override;

};

}

// Include template implementations
#include "levelSetImmutableObject.tpp"

// Explicit instantization
#ifndef __BITPIT_LEVELSET_IMMUTABLE_OBJECT_SRC__
namespace bitpit {

extern template class LevelSetImmutableObject<LevelSetImmutableNarrowBandCache<LevelSetExternalPiercedStorageManager>>;
extern template class LevelSetImmutableObject<LevelSetImmutableNarrowBandCache<LevelSetInternalPiercedStorageManager>>;

}
#endif

#endif
