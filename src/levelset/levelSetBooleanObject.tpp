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

# ifndef __BITPIT_LEVELSET_BOOLEAN_OBJECT_TPP__
# define __BITPIT_LEVELSET_BOOLEAN_OBJECT_TPP__

namespace bitpit {

/*!
	@class      LevelSetBooleanResult
	@ingroup    levelset
	@brief      Allow to evaluate the result of a boolean operation between two LevelSetObjects.
*/

/*!
 * Constructor
 *
 * @param[in] operation type of boolean operation
 */
template<typename SourceLevelSetObject>
LevelSetBooleanResult<SourceLevelSetObject>::LevelSetBooleanResult(LevelSetBooleanOperation operation)
    : m_operation(operation), m_object(nullptr), m_objectSign(0), m_value(levelSetDefaults::VALUE)
{
}

/*!
 * Constructor
 *
 * @param[in] operation type of boolean operation
 * @param[in] object is the object that will be used to initialize the result
 * @param[in] value is the value that will be used to initialize the result
 */
template<typename SourceLevelSetObject>
LevelSetBooleanResult<SourceLevelSetObject>::LevelSetBooleanResult(LevelSetBooleanOperation operation, const SourceLevelSetObject *object, double value)
    : m_operation(operation), m_object(object), m_objectSign(1), m_value(value)
{
}

/*!
 * Update the result.
 *
 * @param[in] object is the object that will be used to update the result
 * @param[in] value is the value that will be used to update the result
 */
template<typename SourceLevelSetObject>
void LevelSetBooleanResult<SourceLevelSetObject>::update(const SourceLevelSetObject *object, double value)
{
    // Early return if the result was not initialized
    if (!m_object) {
            m_object     = object;
            m_objectSign = 1;

            m_value = value;

            return;
    }

    // Update the result
    if( m_operation == LevelSetBooleanOperation::UNION){
        if(m_value > value) {
            m_object     = object;
            m_objectSign = 1;

            m_value = m_objectSign * value;
        }

    } else if ( m_operation == LevelSetBooleanOperation::INTERSECTION){
        if(m_value < value) {
            m_object     = object;
            m_objectSign = 1;

            m_value = m_objectSign * value;
        }

    } else if ( m_operation == LevelSetBooleanOperation::SUBTRACTION){
        if(m_value < - value) {
            m_object     = object;
            m_objectSign = -1;

            m_value = m_objectSign * value;
        }
    }
}

/*!
 * Get the object associated with the results.
 */
template<typename SourceLevelSetObject>
const SourceLevelSetObject * LevelSetBooleanResult<SourceLevelSetObject>::getObject() const {
    return m_object;
}

/*!
 * Get the object sign associated with the results.
 */
template<typename SourceLevelSetObject>
int LevelSetBooleanResult<SourceLevelSetObject>::getObjectSign() const {
    return m_objectSign;
}

/*!
 * Get the value associated with the results.
 */
template<typename SourceLevelSetObject>
double LevelSetBooleanResult<SourceLevelSetObject>::getValue() const {
    return m_value;
}

/*!
	@class      LevelSetBooleanBaseObject
	@ingroup    levelset
	@brief      Base class which deals with boolean operation between two LevelSetObjects
*/

/*!
 * Constructor taking two objects.
 * @param[in] id identifier of object
 * @param[in] op type of boolean operation
 * @param[in] source1 pointer to first source object
 * @param[in] source2 pointer to second source object
 */
template<typename SourceLevelSetObject>
LevelSetBooleanBaseObject<SourceLevelSetObject>::LevelSetBooleanBaseObject( int id, LevelSetBooleanOperation op, const SourceLevelSetObject *source1, const SourceLevelSetObject *source2  )
    : LevelSetProxyObject<SourceLevelSetObject, SourceLevelSetObject>(id) {

    m_operation = op;
    m_sourceObjects.push_back(source1);
    m_sourceObjects.push_back(source2);
}

/*!
 * Constructor taking a vector of objects.
 * The boolean operation will be applied recursivly on each entry.
 * @param[in] id identifier of object
 * @param[in] op type of boolean operation
 * @param[in] sourceObjects pointers to source objects
 */
template<typename SourceLevelSetObject>
LevelSetBooleanBaseObject<SourceLevelSetObject>::LevelSetBooleanBaseObject( int id, LevelSetBooleanOperation op, const std::vector<const SourceLevelSetObject *> &sourceObjects )
    : LevelSetProxyObject<SourceLevelSetObject, SourceLevelSetObject>(id),
      m_operation(op), m_sourceObjects(sourceObjects) {

}

/*!
 * Checks if the object is empty.
 *
 * \result Returns true is the object is empty, false otherwise.
 */
template<typename SourceLevelSetObject>
bool LevelSetBooleanBaseObject<SourceLevelSetObject>::empty() const{

    for (const SourceLevelSetObject *sourceObject : m_sourceObjects) {
        if (!sourceObject->empty()) {
            return false;
        }
    }

    return true;
}

/*
 * Returns the boolean operation
 * @return boolean operation
 */
template<typename SourceLevelSetObject>
LevelSetBooleanOperation LevelSetBooleanBaseObject<SourceLevelSetObject>::getBooleanOperation() const{
    return m_operation;
}

/*!
 * Compute the result of the boolean operation.
 * @param[in] id cell index
 * @param signedLevelSet controls if signed levelset function will be used
 * @return result of the boolean operation.
 */
template<typename SourceLevelSetObject>
LevelSetBooleanResult<SourceLevelSetObject> LevelSetBooleanBaseObject<SourceLevelSetObject>::computeBooleanResult( long id, bool signedLevelSet ) const{

    if (m_sourceObjects.empty()) {
        return LevelSetBooleanResult<SourceLevelSetObject>(getBooleanOperation());
    }

    LevelSetBooleanResult<SourceLevelSetObject> result( getBooleanOperation(), m_sourceObjects[0], m_sourceObjects[0]->evalCellValue(id, signedLevelSet) );
    for( size_t n=1; n<m_sourceObjects.size(); ++n){
        result.update(m_sourceObjects[n], m_sourceObjects[n]->evalCellValue(id, signedLevelSet));
    }

    return result;
}

/*!
 * Compute the result of the boolean operation.
 * @param[in] coords point coordinates
 * @param signedLevelSet controls if signed levelset function will be used
 * @return result of the boolean operation.
 */
template<typename SourceLevelSetObject>
LevelSetBooleanResult<SourceLevelSetObject> LevelSetBooleanBaseObject<SourceLevelSetObject>::computeBooleanResult( const std::array<double,3> &coords, bool signedLevelSet ) const{

    if (m_sourceObjects.empty()) {
        return LevelSetBooleanResult<SourceLevelSetObject>(getBooleanOperation());
    }

    LevelSetBooleanResult<SourceLevelSetObject> result( getBooleanOperation(), m_sourceObjects[0], m_sourceObjects[0]->evalValue(coords, signedLevelSet) );
    for( size_t n=1; n<m_sourceObjects.size(); ++n){
        result.update(m_sourceObjects[n], m_sourceObjects[n]->evalValue(coords, signedLevelSet));
    }

    return result;
}

/*!
 * Replace a source object.
 * @param[in] current current source object
 * @param[in] updated updated source object
 */
template<typename SourceLevelSetObject>
void LevelSetBooleanBaseObject<SourceLevelSetObject>::replaceSourceObject(const SourceLevelSetObject *current, const SourceLevelSetObject *updated){

    std::size_t nSources = m_sourceObjects.size();
    for (std::size_t i = 0; i < nSources; ++i) {
        const SourceLevelSetObject *source = m_sourceObjects[i];
        if (source == current) {
            m_sourceObjects[i] = updated;
            return;
        }
    }

    throw std::runtime_error("Unable to find the source that should be replaced.");
}

/*!
 * Fill the cache that contains the propagated cell sign.
 */
template<typename SourceLevelSetObject>
void LevelSetBooleanBaseObject<SourceLevelSetObject>::fillCellPropagatedSignCache()
{
    // Early return if propagated sign cannot be copied from the source object
    for (const SourceLevelSetObject *sourceObject : m_sourceObjects) {
        LevelSetBulkEvaluationMode sourceBulkEvaluationMode = sourceObject->getCellBulkEvaluationMode();

        bool useSourceSign = false;
        if (sourceBulkEvaluationMode == LevelSetBulkEvaluationMode::SIGN_PROPAGATION) {
            useSourceSign = true;
        } else if (sourceBulkEvaluationMode == LevelSetBulkEvaluationMode::EXACT) {
            LevelSetCacheMode sourceSignCacheMode = sourceObject->getFieldCellCacheMode(LevelSetField::SIGN);
            if (sourceSignCacheMode == LevelSetCacheMode::FULL) {
                useSourceSign = true;
            } else {
                LevelSetCacheMode sourceValueCacheMode = sourceObject->getFieldCellCacheMode(LevelSetField::VALUE);
                if (sourceValueCacheMode == LevelSetCacheMode::FULL) {
                    useSourceSign = true;
                }
            }
        }

        if (!useSourceSign) {
            SourceLevelSetObject::fillCellPropagatedSignCache();
            return;
        }
    }

    // Mesh information
    const VolumeKernel &mesh = *(this->getKernel()->getMesh()) ;
    VolumeKernel::CellConstIterator cellBegin = mesh.cellConstBegin();
    VolumeKernel::CellConstIterator cellEnd   = mesh.cellConstEnd();

    // Get cache for sign propagation
    typedef typename SourceLevelSetObject::CellCacheCollection::template ValueCache<char> ZoneCache;
    ZoneCache *propagatedSignCache = this->template getCellCache<char>(this->m_cellPropagatedSignCacheId);

    // Get propagated sign from source objects
    for (VolumeKernel::CellConstIterator cellItr = cellBegin; cellItr != cellEnd; ++cellItr) {
        long cellId = cellItr.getId();
        const LevelSetBooleanResult<SourceLevelSetObject> cellResult = computeBooleanResult(cellId, true);
        short cellSign = this->evalValueSign(cellResult.getValue());

        propagatedSignCache->insertEntry(cellId, static_cast<char>(cellSign));
    }
}

/*!
 * Evaluate levelset sign at the specified cell.
 *
 * \param id is the id of the cell
 * \result The sign of the levelset at the specified cell.
 */
template<typename SourceLevelSetObject>
short LevelSetBooleanBaseObject<SourceLevelSetObject>::_evalCellSign(long id) const {
    return this->evalValueSign(_evalCellValue(id, true));
}

/*!
 * Evaluate levelset value at the specified cell.
 *
 * \param id is the id of the cell
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The value of the levelset at the specified cell.
 */
template<typename SourceLevelSetObject>
double LevelSetBooleanBaseObject<SourceLevelSetObject>::_evalCellValue(long id, bool signedLevelSet) const {
    return _evalCellFunction<double>(id, signedLevelSet, [] (const LevelSetBooleanResult<SourceLevelSetObject> &result)
        {
            const LevelSetObject *resultObject = result.getObject();
            if ( !resultObject ) {
                return levelSetDefaults::VALUE ;
            }

            return result.getValue();
        });
}

/*!
 * Evaluate levelset gradient at the specified cell.
 *
 * \param id is the id of the cell
 * \result The gradient of the levelset at the specified cell.
 */
template<typename SourceLevelSetObject>
std::array<double,3> LevelSetBooleanBaseObject<SourceLevelSetObject>::_evalCellGradient(long id, bool signedLevelSet) const {
    return _evalCellFunction<std::array<double,3>>(id, signedLevelSet, [id, signedLevelSet] (const LevelSetBooleanResult<SourceLevelSetObject> &result)
        {
            const LevelSetObject *resultObject = result.getObject();
            if ( !resultObject ) {
                return levelSetDefaults::GRADIENT ;
            }

            std::array<double,3> gradient = resultObject->evalCellGradient(id, signedLevelSet);
            if (signedLevelSet) {
                return gradient;
            }

            return static_cast<double>(result.getObjectSign()) * gradient;
        });
}

/*!
 * Evaluate levelset value at the specified cell.
 *
 * \param point are the coordinates of the point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The value of the levelset at the specified cell.
 */
template<typename SourceLevelSetObject>
double LevelSetBooleanBaseObject<SourceLevelSetObject>::_evalValue(const std::array<double,3> &point, bool signedLevelSet) const {
    return _evalFunction<double>(point, signedLevelSet, [] (const LevelSetBooleanResult<SourceLevelSetObject> &result)
        {
            const SourceLevelSetObject *resultObject = result.getObject();
            if ( !resultObject ) {
                return levelSetDefaults::VALUE ;
            }

            return result.getValue();
        });
}

/*!
 * Evaluate levelset gradient at the specified cell.
 *
 * \param point are the coordinates of the point
 * \param signedLevelSet controls if signed levelset function will be used
 * \result The gradient of the levelset at the specified cell.
 */
template<typename SourceLevelSetObject>
std::array<double,3> LevelSetBooleanBaseObject<SourceLevelSetObject>::_evalGradient(const std::array<double,3> &point, bool signedLevelSet) const {
    return _evalFunction<std::array<double,3>>(point, signedLevelSet, [&point, signedLevelSet] (const LevelSetBooleanResult<SourceLevelSetObject> &result)
        {
            const LevelSetObject *resultObject = result.getObject();
            if ( !resultObject ) {
                return levelSetDefaults::GRADIENT ;
            }

            std::array<double,3> gradient = resultObject->evalGradient(point, signedLevelSet);
            if (signedLevelSet) {
                return gradient;
            }

            return static_cast<double>(result.getObjectSign()) * gradient;
        });
}

/*!
 * Evaluate the specified function at the specified cell.
 * @param[in] id cell index
 * @param[in] function is the function that will be evaluated
 * @return The value of the function at specified cell.
 */
template<typename SourceLevelSetObject>
template<typename data_t, typename function_t>
data_t LevelSetBooleanBaseObject<SourceLevelSetObject>::_evalCellFunction(long id, bool signedLevelSet,
                                                                          const function_t &function) const
{
    const LevelSetBooleanResult<SourceLevelSetObject> result = computeBooleanResult( id, signedLevelSet ) ;

    return function(result) ;
}

/*!
 * Evaluate the specified function at the specified point.
 * @param point are the coordinates of the point
 * @param[in] function is the function that will be evaluated
 * @return The value of the function at specified point.
 */
template<typename SourceLevelSetObject>
template<typename data_t, typename function_t>
data_t LevelSetBooleanBaseObject<SourceLevelSetObject>::_evalFunction(const std::array<double,3> &point, bool signedLevelSet,
                                                                      const function_t &function) const
{
    const LevelSetBooleanResult<SourceLevelSetObject> result = computeBooleanResult( point, signedLevelSet ) ;

    return function(result) ;
}

/*!
 * Get the object that defines the levelset information for the specified cell.
 * @param[in] id cell index
 * @return The object that defines the levelset information for the specified
 * cell.
 */
template<typename SourceLevelSetObject>
const SourceLevelSetObject * LevelSetBooleanBaseObject<SourceLevelSetObject>::getCellReferenceObject(long id) const{

    // Early return if the object has no sources
    if (m_sourceObjects.empty()) {
        return nullptr;
    }

    // Early return if the object has only one source
    if (m_sourceObjects.size() == 1) {
        return m_sourceObjects.front();
    }

    // Evaluate reference object from boolean result
    const LevelSetBooleanResult<SourceLevelSetObject> result = computeBooleanResult(id, true) ;

    return result.getObject();

}

/*!
 * Get the object that defines the levelset information for the specified point.
 * @param[in] point are the coordinates of the point
 * @return The object that defines the levelset information for the specified
 * point.
 */
template<typename SourceLevelSetObject>
const SourceLevelSetObject * LevelSetBooleanBaseObject<SourceLevelSetObject>::getReferenceObject(const std::array<double, 3> &point) const{

    // Early return if the object has no sources
    if (m_sourceObjects.empty()) {
        return nullptr;
    }

    // Early return if the object has only one source
    if (m_sourceObjects.size() == 1) {
        return m_sourceObjects.front();
    }

    // Evaluate reference object from boolean result
    const LevelSetBooleanResult<SourceLevelSetObject> result = computeBooleanResult(point, true) ;

    return result.getObject();

}

/*!
 * Get all objects that compose the boolean object
 * \return pointers to all primary objects involved in the definition of the
 * boolean object levelset information.
 */
template<typename SourceLevelSetObject>
std::vector<const SourceLevelSetObject *> LevelSetBooleanBaseObject<SourceLevelSetObject>::getSourceObjects() const{

    return m_sourceObjects;

}

}

#endif
