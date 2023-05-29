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

#include <bitpit_common.hpp>

namespace bitpit{

class LevelSet;
class LevelSetObject;
class LevelSetObject;

class LevelSetProxyObject : public LevelSetObject {
    friend class LevelSet;

    protected:
    virtual void replaceSourceObject(const LevelSetObject *current, const LevelSetObject *updated) = 0;

    public:
    LevelSetProxyObject(int);

    bool            isPrimary() const override;

    virtual const LevelSetObject *    getReferenceObject( long ) const =0;
    virtual const LevelSetObject *    getReferencePrimaryObject( long ) const;

    int getReferenceObjectId( long ) const;
    int getReferencePrimaryObjectId( long ) const;
    BITPIT_DEPRECATED(int getPrimaryObjectId( long ) const);

    virtual std::vector<const LevelSetObject *> getSourceObjects() const =0;
    virtual std::vector<const LevelSetObject *> getPrimarySourceObjects() const;

    std::vector<int> getSourceObjectIds() const;
    std::vector<int> getPrimarySourceObjectIds() const;

};

// Compatibility with older versions
typedef LevelSetProxyObject LevelSetMetaObject;

}

#endif 
