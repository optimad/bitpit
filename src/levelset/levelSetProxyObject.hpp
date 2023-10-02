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

# ifndef __BITPIT_LEVELSET_PROXY_OBJECT_HPP__
# define __BITPIT_LEVELSET_PROXY_OBJECT_HPP__

#include <levelSetObject.hpp>

#include <bitpit_common.hpp>

namespace bitpit{

class LevelSetProxyBaseObject {
    public:
    virtual int getReferenceObjectId( long ) const = 0;
    virtual int getReferencePrimaryObjectId( long ) const = 0;

    virtual std::vector<int> getSourceObjectIds() const = 0;
    virtual std::vector<int> getPrimarySourceObjectIds() const = 0;

};

template<typename SourceLevelSetObject, typename BaseLevelSetObject = LevelSetObject>
class LevelSetProxyObject : public BaseLevelSetObject, public LevelSetProxyBaseObject {
    protected:
    LevelSetProxyObject(int);

    virtual void replaceSourceObject(const SourceLevelSetObject *current, const SourceLevelSetObject *updated) = 0;

    public:
    bool            isPrimary() const override;

    bool            isInNarrowBand(long id) const override;

    virtual const SourceLevelSetObject *    getReferenceObject( long ) const =0;
    virtual const SourceLevelSetObject *    getReferencePrimaryObject( long ) const;

    int getReferenceObjectId( long ) const override;
    int getReferencePrimaryObjectId( long ) const override;
    BITPIT_DEPRECATED(int getPrimaryObjectId( long ) const);

    virtual std::vector<const SourceLevelSetObject *> getSourceObjects() const =0;
    virtual std::vector<const SourceLevelSetObject *> getPrimarySourceObjects() const;

    std::vector<int> getSourceObjectIds() const override;
    std::vector<int> getPrimarySourceObjectIds() const override;

};

// Compatibility with older versions
typedef LevelSetProxyObject<LevelSetObject> LevelSetMetaObject;

}

// Include template implementations
#include "levelSetProxyObject.tpp"

#endif 
