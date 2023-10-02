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

#include <array>

namespace bitpit{

class LevelSetProxyBaseObject {
    public:
    virtual int getCellReferenceObjectId( long id ) const = 0;
    virtual int getCellReferencePrimaryObjectId( long id ) const = 0;

    virtual int getReferenceObjectId( const std::array<double, 3> &point ) const = 0;
    virtual int getReferencePrimaryObjectId( const std::array<double, 3> &point ) const = 0;

    virtual std::vector<int> getSourceObjectIds() const = 0;
    virtual std::vector<int> getPrimarySourceObjectIds() const = 0;

};

template<typename SourceLevelSetObject, typename BaseLevelSetObject = LevelSetObject>
class LevelSetProxyObject : public BaseLevelSetObject, public LevelSetProxyBaseObject {
    protected:
    LevelSetProxyObject(int);

    void fillCellLocationCache() override;
    void fillCellLocationCache(const std::vector<adaption::Info> &adaptionData) override;

    virtual void replaceSourceObject(const SourceLevelSetObject *current, const SourceLevelSetObject *updated) = 0;

    public:
    bool isPrimary() const override;

    bool isCellInNarrowBand(long id) const override;
    bool isInNarrowBand(const std::array<double,3> &point) const override;

    virtual const SourceLevelSetObject * getCellReferenceObject( long id ) const =0;
    virtual const SourceLevelSetObject * getCellReferencePrimaryObject( long id ) const;

    virtual const SourceLevelSetObject * getReferenceObject( const std::array<double, 3> &point ) const =0;
    virtual const SourceLevelSetObject * getReferencePrimaryObject( const std::array<double, 3> &point ) const;

    int getCellReferenceObjectId( long id ) const override;
    int getCellReferencePrimaryObjectId( long id ) const override;

    int getReferenceObjectId( const std::array<double, 3> &point ) const override;
    int getReferencePrimaryObjectId( const std::array<double, 3> &point ) const override;

    virtual std::vector<const SourceLevelSetObject *> getSourceObjects() const =0;
    virtual std::vector<const SourceLevelSetObject *> getPrimarySourceObjects() const;

    std::vector<int> getSourceObjectIds() const override;
    std::vector<int> getPrimarySourceObjectIds() const override;

    BITPIT_DEPRECATED(int getPrimaryObjectId( long ) const);

};

// Compatibility with older versions
typedef LevelSetProxyObject<LevelSetObject> LevelSetMetaObject;

}

// Include template implementations
#include "levelSetProxyObject.tpp"

#endif 
