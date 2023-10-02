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
    : LevelSetProxyObject<SourceLevelSetObject>(id) {

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
    : LevelSetProxyObject<SourceLevelSetObject>(id),
      m_operation(op), m_sourceObjects(sourceObjects) {

}

/*!
 * Get the levelset value
 * @param[in] id cell id
 * @return levelset value in cell
 */
template<typename SourceLevelSetObject>
double LevelSetBooleanBaseObject<SourceLevelSetObject>::getValue( long id)const {
    const LevelSetBooleanResult<SourceLevelSetObject> result = computeBooleanResult( id ) ;
    const SourceLevelSetObject *resultObject = result.getObject();
    if ( resultObject ) {
        double value = result.getValue();

        return value ;
    }

    return levelSetDefaults::VALUE ;
}

/*!
 * Get the levelset gradient
 * @param[in] id cell id
 * @return levelset gradient in cell
 */
template<typename SourceLevelSetObject>
std::array<double,3> LevelSetBooleanBaseObject<SourceLevelSetObject>::getGradient(long id) const {
    const LevelSetBooleanResult<SourceLevelSetObject> result = computeBooleanResult( id ) ;
    const SourceLevelSetObject *resultObject = result.getObject();
    if ( resultObject ) {
        std::array<double, 3> gradient = static_cast<double>(result.getObjectSign()) * resultObject->getGradient( id ) ;

        return gradient ;
    }

    return levelSetDefaults::GRADIENT ;
}

/*!
 * Computes the LevelSetInfo in a point
 * @param[in] coords point coordinates
 * @return LevelSetInfo
*/
template<typename SourceLevelSetObject>
LevelSetInfo LevelSetBooleanBaseObject<SourceLevelSetObject>::computeLevelSetInfo( const std::array<double,3> &coords) const{
    const LevelSetBooleanResult<SourceLevelSetObject> result = computeBooleanResult( coords ) ;
    const SourceLevelSetObject *componentObject = result.getObject();
    if ( componentObject ) {
        LevelSetInfo levelSetInfo = componentObject->computeLevelSetInfo( coords ) ;
        levelSetInfo.value    *= result.getObjectSign();
        levelSetInfo.gradient *= static_cast<double>(result.getObjectSign());

        return levelSetInfo;
    }

    return LevelSetInfo() ;
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

/*
 * Returns the boolean operation
 * @return boolean operation
 */
template<typename SourceLevelSetObject>
LevelSetBooleanOperation LevelSetBooleanBaseObject<SourceLevelSetObject>::getBooleanOperation() const{
    return m_operation;
}

/*!
 * Get the object that defines the levelset information for the specified cell.
 * @param[in] id cell index
 * @return The object that defines the levelset information for the specified
 * cell.
 */
template<typename SourceLevelSetObject>
const SourceLevelSetObject * LevelSetBooleanBaseObject<SourceLevelSetObject>::getReferenceObject(long id) const{

    const LevelSetBooleanResult<SourceLevelSetObject> result = computeBooleanResult(id) ;

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

/*!
 * Compute the result of the boolean operation.
 * @param[in] id cell index
 * @return result of the boolean operation.
 */
template<typename SourceLevelSetObject>
LevelSetBooleanResult<SourceLevelSetObject> LevelSetBooleanBaseObject<SourceLevelSetObject>::computeBooleanResult( long id ) const{

    // Early return if the are no objects
    if (m_sourceObjects.empty()) {
        return LevelSetBooleanResult<SourceLevelSetObject>( getBooleanOperation() );
    }

    // Identify information about the source
    LevelSetBooleanResult<SourceLevelSetObject> result( getBooleanOperation(), m_sourceObjects[0], m_sourceObjects[0]->getValue(id) ) ;
    for( size_t n=1; n<m_sourceObjects.size(); ++n){
        result.update(m_sourceObjects[n], m_sourceObjects[n]->getValue(id));
    }

    return result;
}

/*!
 * Compute the result of the boolean operation.
 * @param[in] coords point coordinates
 * @return result of the boolean operation.
 */
template<typename SourceLevelSetObject>
LevelSetBooleanResult<SourceLevelSetObject> LevelSetBooleanBaseObject<SourceLevelSetObject>::computeBooleanResult( const std::array<double,3> &coords ) const{

    // Early return if the are no objects
    if (m_sourceObjects.empty()) {
        return LevelSetBooleanResult<SourceLevelSetObject>( getBooleanOperation() );
    }

    // Identify information about the source
    LevelSetBooleanResult<SourceLevelSetObject> result( getBooleanOperation(), m_sourceObjects[0], m_sourceObjects[0]->computeLevelSetInfo(coords).value);
    for( size_t n=1; n<m_sourceObjects.size(); ++n){
        result.update(m_sourceObjects[n], m_sourceObjects[n]->computeLevelSetInfo( coords ).value );
    }

    return result;
}

}

#endif
