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

# ifndef __BITPIT_LEVELSET_BOOLEAN_OBJECT_HPP__
# define __BITPIT_LEVELSET_BOOLEAN_OBJECT_HPP__

// Standard Template Library
# include <array>
# include <vector>

# include "levelSetCommon.hpp"
# include "levelSetProxyObject.hpp"

namespace bitpit{

class LevelSetObject ;

namespace adaption{
    struct Info;
}

template<typename SourceLevelSetObject>
class LevelSetBooleanResult {

    private:
    LevelSetBooleanOperation                    m_operation;            /**< Boolean operation */

    const SourceLevelSetObject                 *m_object;               /**< Object that defines the result */
    int                                         m_objectSign;           /**< Sign associated with the object */

    double                                      m_value;                /**< Value of the results */

    public:
    LevelSetBooleanResult(LevelSetBooleanOperation operation) ;
    LevelSetBooleanResult(LevelSetBooleanOperation operation, const SourceLevelSetObject *object, double value) ;

    void update(const SourceLevelSetObject *object, double value) ;

    const SourceLevelSetObject * getObject() const ;
    int getObjectSign () const ;

    double getValue() const ;

};

template<typename SourceLevelSetObject>
class LevelSetBooleanBaseObject : public LevelSetProxyObject<SourceLevelSetObject, SourceLevelSetObject> {

    private:
    LevelSetBooleanOperation                    m_operation;            /**< identifier of operation */
    std::vector<const SourceLevelSetObject *>   m_sourceObjects;        /**< Pointers to source objects */

    LevelSetBooleanOperation                    getBooleanOperation() const;

    LevelSetBooleanResult<SourceLevelSetObject> computeBooleanResult( long, bool signedLevelSet ) const ;
    LevelSetBooleanResult<SourceLevelSetObject> computeBooleanResult( const std::array<double,3> &coords, bool signedLevelSet ) const ;

    protected:
    LevelSetBooleanBaseObject(int, LevelSetBooleanOperation, const SourceLevelSetObject *, const SourceLevelSetObject *);
    LevelSetBooleanBaseObject(int, LevelSetBooleanOperation, const std::vector<const SourceLevelSetObject *> &);

    void                                        replaceSourceObject(const SourceLevelSetObject *current, const SourceLevelSetObject *updated) override ;

    void                                        fillCellPropagatedSignCache() override;

    short                                       _evalCellSign(long id) const override;
    double                                      _evalCellValue(long id, bool signedLevelSet) const override;
    std::array<double,3>                        _evalCellGradient(long id, bool signedLevelSet) const override;

    double                                      _evalValue(const std::array<double,3> &point, bool signedLevelSet) const override;
    std::array<double,3>                        _evalGradient(const std::array<double,3> &point, bool signedLevelSet) const override;

    public:
    bool                                        empty() const override;

    const SourceLevelSetObject *                getCellReferenceObject(long id) const override;

    const SourceLevelSetObject *                getReferenceObject(const std::array<double,3> &point) const override;

    std::vector<const SourceLevelSetObject *>   getSourceObjects() const override;

protected:
    template<typename data_t, typename function_t>
    data_t _evalCellFunction(long id, bool signedLevelSet, const function_t &function) const;

    template<typename data_t, typename function_t>
    data_t _evalFunction(const std::array<double,3> &point, bool signedLevelSet, const function_t &function) const;

};

template<typename SourceLevelSetObject>
class LevelSetBooleanObject : public LevelSetBooleanBaseObject<SourceLevelSetObject> {

};

template<>
class LevelSetBooleanObject<LevelSetObject> : public LevelSetBooleanBaseObject<LevelSetObject> {

public:
    LevelSetBooleanObject(int, LevelSetBooleanOperation, const LevelSetObject *, const LevelSetObject *);
    LevelSetBooleanObject(int, LevelSetBooleanOperation, const std::vector<const LevelSetObject *> &);

    LevelSetBooleanObject<LevelSetObject> * clone() const override;

};

// Typdefs for compatibility with older versions
typedef LevelSetBooleanObject<LevelSetObject> LevelSetBoolean;

}

#endif
