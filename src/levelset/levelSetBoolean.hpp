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
# include "levelSetProxyObject.hpp"

namespace bitpit{

namespace adaption{
    struct Info;
}

class LevelSetObject ;
class LevelSetProxyObject ;

class LevelSetBooleanResult {

    private:
    LevelSetBooleanOperation                    m_operation;            /**< Boolean operation */

    const LevelSetObject                       *m_object;               /**< Object that defines the result */
    int                                         m_objectSign;           /**< Sign associated with the object */

    double                                      m_value;                /**< Value of the results */

    public:
    LevelSetBooleanResult(LevelSetBooleanOperation operation) ;
    LevelSetBooleanResult(LevelSetBooleanOperation operation, const LevelSetObject *object, double value) ;

    void update(const LevelSetObject *object, double value) ;

    const LevelSetObject * getObject() const ;
    int getObjectSign () const ;

    double getValue() const ;

};

class LevelSetBoolean: public LevelSetProxyObject {

    private:
    LevelSetBooleanOperation                    m_operation;            /**< identifier of operation */
    std::vector<const LevelSetObject*>          m_sourceObjects;        /**< Pointers to source objects */

    LevelSetBooleanOperation                    getBooleanOperation() const;

    LevelSetBooleanResult                       computeBooleanResult( long ) const ;
    LevelSetBooleanResult                       computeBooleanResult( const std::array<double,3> &coords ) const ;

    protected:
    void                                        _dump( std::ostream &) override;
    void                                        _restore( std::istream &) override;

    public:
    LevelSetBoolean(int, LevelSetBooleanOperation, const LevelSetObject*, const LevelSetObject*);
    LevelSetBoolean(int, LevelSetBooleanOperation, const std::vector<const LevelSetObject*> &);
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

    const LevelSetObject *                      getReferenceObject( long ) const override;

    std::vector<const LevelSetObject *>         getSourceObjects() const override;

    bool                                        isInNarrowBand(long ) const override;

};

}

#endif /*__BITPIT_LEVELSET_BOOLEAN_HPP__ */
