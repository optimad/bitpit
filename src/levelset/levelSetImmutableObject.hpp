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

#ifndef __BITPIT_LEVELSET_IMMUTABLE_HPP__
#define __BITPIT_LEVELSET_IMMUTABLE_HPP__

#include "levelSetObject.hpp"
#include "levelSetBoolean.hpp"
#include "levelSetCachedObject.hpp"
#include "levelSetSignedObject.hpp"

namespace bitpit {

class LevelSetImmutableObject : public LevelSetObject, public LevelSetCachedObjectInterface, public LevelSetSignedObjectInterface
{

public:
    LevelSetImmutableObject(LevelSetCachedObject *object);
    LevelSetImmutableObject(LevelSetBoolean *object);

    LevelSetImmutableObject * clone() const override;

    LevelSetInfo getLevelSetInfo(long id) const override;
    double getLS(long id) const override;
    double getValue(long id) const override;
    short getSign(long id) const override;
    std::array<double,3> getGradient(long id) const override;

    LevelSetInfo computeLevelSetInfo(const std::array<double,3> &coords) const override;

protected:
    std::shared_ptr<LevelSetNarrowBandCache> createNarrowBandCache() override;

    std::shared_ptr<LevelSetSignStorage> createSignStorage() override;

    LevelSetImmutableObject(int);

    void _clear() override;

    void _dump(std::ostream &stream) override;
    void _restore(std::istream &stream) override;

};

}

#endif
