/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

# ifndef __BITPIT_LEVELSET_BOOLEAN_HPP__
# define __BITPIT_LEVELSET_BOOLEAN_HPP__

// Standard Template Library
# include <array>
# include <vector>

# include "levelSetCommon.hpp"

namespace bitpit{

namespace adaption{
    struct Info;
}

class LevelSetObject ;
class LevelSetMetaObject ;

class LevelSetBoolean: public LevelSetMetaObject {

    private:
    LevelSetBooleanOperation                    m_operation;            /**< identifier of operation */
    std::vector<LevelSetObject*>                m_objPtr;               /**< pointers to objects */

    LevelSetInfo                                booleanOperation(const long &) const ;
    LevelSetBooleanOperation                    getBooleanOperation() const;
    LevelSetObject*                             getCompetentObject(const long &, double *factorPtr=nullptr) const ;

    protected:
    void                                        _dump( std::ostream &);
    void                                        _restore( std::istream &);

    public:
    ~LevelSetBoolean();
    LevelSetBoolean(int, LevelSetBooleanOperation, LevelSetObject*, LevelSetObject*);
    LevelSetBoolean(int, LevelSetBooleanOperation, std::vector<LevelSetObject*>);
    LevelSetBoolean(const LevelSetBoolean &);

    LevelSetBoolean*                            clone() const ;

    LevelSetInfo                                getLevelSetInfo(const long &) const;
    double                                      getLS(const long &) const; 
    std::array<double,3>                        getGradient(const long &) const; 

    std::array<double,3>                        getNormal(const long &) const ;
    int                                         getPart(const long &) const ;
    double                                      getSurfaceFeatureSize(const long &) const;
    double                                      getMinSurfaceFeatureSize() const;
    double                                      getMaxSurfaceFeatureSize() const;

    void                                        computeLSInNarrowBand(bool) ;
    void                                        updateLSInNarrowBand(const std::vector<adaption::Info> &, bool);

    int                                         getPrimaryObjectId(const long &) const;

};

}

#endif /*__BITPIT_LEVELSET_BOOLEAN_HPP__ */
