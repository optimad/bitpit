/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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
    class Info ;
}
class LevelSetObject ;
class LevelSetKernel ;

class LevelSetBoolean: public LevelSetObject{

    private:
    LevelSetBooleanOperation                    m_operation;            /**< identifier of operation */
    LevelSetObject*                             m_objPtr1;              /**< pointer to first object */
    LevelSetObject*                             m_objPtr2;              /**< pointer to second object */
    int                                         m_objId1;               /**< identifier of first object */
    int                                         m_objId2;               /**< identifier of second object */

    protected:
    void                                        _dump( std::ostream &);
    void                                        _restore( std::istream &);

    public:
    ~LevelSetBoolean();
    LevelSetBoolean(int, LevelSetBooleanOperation, LevelSetObject*, LevelSetObject*);
    LevelSetBoolean(const LevelSetBoolean &);

    LevelSetBoolean*                            clone() const ;

    LevelSetInfo                                getLevelSetInfo(const long &) const;
    double                                      getLS(const long &) const; 
    std::array<double,3>                        getGradient(const long &) const; 

    int                                         getPart(const long &) const ;
    long                                        getSupport(const long &) const;
    int                                         getSupportCount(const long &) const ;

    void                                        setSizeNarrowBand(double) ;

    double                                      computeSizeNarrowBand(LevelSetKernel*);
    double                                      updateSizeNarrowBand(LevelSetKernel*,const std::vector<adaption::Info> &);
    void                                        computeLSInNarrowBand( LevelSetKernel *, const double &, const bool &) ;
    void                                        updateLSInNarrowBand( LevelSetKernel *, const std::vector<adaption::Info> &, const double &, const bool &);

    LevelSetBooleanOperation                    getBooleanOperation() const;
    LevelSetObject*                             getClosestObject(const long &) const ;
    LevelSetInfo                                booleanOperation(const long &) const ;
};

}

#endif /*__BITPIT_LEVELSET_BOOLEAN_HPP__ */
