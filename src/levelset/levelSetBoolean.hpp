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

class LevelSetBooleanResult {

    private:
    LevelSetBooleanOperation                    m_operation;            /**< Boolean operation */

    LevelSetObject                             *m_object;               /**< Object that defines the result */
    int                                         m_objectSign;           /**< Sign associated with the object */

    double                                      m_value;                /**< Value of the results */

    public:
    LevelSetBooleanResult(LevelSetBooleanOperation operation) ;
    LevelSetBooleanResult(LevelSetBooleanOperation operation, LevelSetObject *object, double value) ;

    void update(LevelSetObject *object, double value) ;

    LevelSetObject * getObject() const ;
    int getObjectSign () const ;

    double getValue() const ;

};

class LevelSetBoolean: public LevelSetMetaObject {

    private:
    LevelSetBooleanOperation                    m_operation;            /**< identifier of operation */
    std::vector<LevelSetObject*>                m_objects;              /**< pointers to objects */

    LevelSetBooleanOperation                    getBooleanOperation() const;

    LevelSetBooleanResult                       computeBooleanResult( long ) const ;
    LevelSetBooleanResult                       computeBooleanResult( const std::array<double,3> &coords ) const ;

    protected:
    void                                        _dump( std::ostream &) override;
    void                                        _restore( std::istream &) override;

    public:
    LevelSetBoolean(int, LevelSetBooleanOperation, LevelSetObject*, LevelSetObject*);
    LevelSetBoolean(int, LevelSetBooleanOperation, const std::vector<LevelSetObject*> &);
    LevelSetBoolean(const LevelSetBoolean &);

    LevelSetBoolean*                            clone() const override;

    LevelSetInfo                                getLevelSetInfo(long ) const override;
    double                                      getLS(long ) const override;
    double                                      getValue(long ) const override;
    std::array<double,3>                        getGradient(long ) const override;

    std::array<double,3>                        getNormal(long ) const override;
    int                                         getPart(long ) const override;
    double                                      getSurfaceFeatureSize(long ) const override;

    LevelSetInfo                                computeLevelSetInfo(const std::array<double,3> &) const override;

    double                                      getMinSurfaceFeatureSize() const override;
    double                                      getMaxSurfaceFeatureSize() const override;

    void                                        computeLSInNarrowBand(bool) override;
    void                                        updateLSInNarrowBand(const std::vector<adaption::Info> &, bool) override;

    int                                         getPrimaryObjectId(long ) const override;
    std::vector<const LevelSetObject*>          getPrimaryObjects() const override;

    bool                                        isInNarrowBand(long ) const override;

};

}

#endif /*__BITPIT_LEVELSET_BOOLEAN_HPP__ */
