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

# ifndef __BITPIT_LEVELSET_COMPLEMENT_OBJECT_HPP__
# define __BITPIT_LEVELSET_COMPLEMENT_OBJECT_HPP__

# include <array>

# include "levelSetCommon.hpp"
# include "levelSetProxyObject.hpp"

namespace bitpit{

class LevelSetObject ;

template<typename SourceLevelSetObject>
class LevelSetComplementBaseObject : public LevelSetProxyObject<SourceLevelSetObject, SourceLevelSetObject> {

    private:
    const SourceLevelSetObject *                         m_sourceObject;        /**< Pointers to source object */

    protected:
    LevelSetComplementBaseObject(int id, const SourceLevelSetObject *source);

    void                                                 replaceSourceObject(const SourceLevelSetObject *current, const SourceLevelSetObject *updated) override;

    void                                                 fillCellPropagatedSignCache() override;

    short                                                _evalCellSign(long id) const override;
    double                                               _evalCellValue(long id, bool signedLevelSet) const override;
    std::array<double,3>                                 _evalCellGradient(long id, bool signedLevelSet) const override;

    short                                                _evalSign(const std::array<double,3> &point) const override;
    double                                               _evalValue(const std::array<double,3> &point, bool signedLevelSet) const override;
    std::array<double,3>                                 _evalGradient(const std::array<double,3> &point, bool signedLevelSet) const override;

    public:
    bool                                                 empty() const override;

    const SourceLevelSetObject *                         getCellReferenceObject(long id) const override;

    const SourceLevelSetObject *                         getReferenceObject(const std::array<double, 3> &point) const override;

    virtual const SourceLevelSetObject *                 getSourceObject() const;
    std::vector<const SourceLevelSetObject *>            getSourceObjects() const override;

};

template<typename SourceLevelSetObject>
class LevelSetComplementObject : public LevelSetComplementBaseObject<SourceLevelSetObject> {

};

template<>
class LevelSetComplementObject<LevelSetObject> : public LevelSetComplementBaseObject<LevelSetObject> {

public:
    LevelSetComplementObject(int id, const LevelSetObject *source);

    LevelSetComplementObject<LevelSetObject> * clone() const override;

};

}

#endif
