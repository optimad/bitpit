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
class LevelSetProxyObject ;

class LevelSetComplementObject: public LevelSetProxyObject {

    private:
    const LevelSetObject*                       m_sourceObject;        /**< Pointers to source object */

    protected:
    void                                        replaceSourceObject(const LevelSetObject *current, const LevelSetObject *updated) override ;

    public:
    LevelSetComplementObject(int id, const LevelSetObject*);
    LevelSetComplementObject(const LevelSetComplementObject &);

    LevelSetComplementObject*                   clone() const override;

    double                                      getValue(long ) const override;
    std::array<double,3>                        getGradient(long ) const override;

    std::array<double,3>                        getNormal(long ) const override;
    int                                         getPart(long ) const override;
    double                                      getSurfaceFeatureSize(long ) const override;

    LevelSetInfo                                computeLevelSetInfo(const std::array<double,3> &) const override;

    double                                      getMinSurfaceFeatureSize() const override;
    double                                      getMaxSurfaceFeatureSize() const override;

    const LevelSetObject *                      getReferenceObject( long ) const override;

    std::vector<const LevelSetObject *>         getSourceObjects() const override;

};

}

#endif
