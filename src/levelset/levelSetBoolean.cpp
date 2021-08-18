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

# include <cassert>

# include "bitpit_common.hpp"
# include "bitpit_operators.hpp"

# include "logger.hpp"
# include "adaption.hpp"

# include "levelSetObject.hpp"
# include "levelSetProxyObject.hpp"
# include "levelSetBoolean.hpp"

namespace bitpit {

/*!
	@class      LevelSetBooleanResult
	@ingroup    levelset
	@brief      Class for evaluating boolean results.
*/

/*!
 * Constructor
 *
 * @param[in] operation type of boolean operation
 */
LevelSetBooleanResult::LevelSetBooleanResult(LevelSetBooleanOperation operation)
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
LevelSetBooleanResult::LevelSetBooleanResult(LevelSetBooleanOperation operation, const LevelSetObject *object, double value)
    : m_operation(operation), m_object(object), m_objectSign(1), m_value(value)
{
}

/*!
 * Update the result.
 *
 * @param[in] object is the object that will be used to update the result
 * @param[in] value is the value that will be used to update the result
 */
void LevelSetBooleanResult::update(const LevelSetObject *object, double value)
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
const LevelSetObject * LevelSetBooleanResult::getObject() const {
    return m_object;
}

/*!
 * Get the object sign associated with the results.
 */
int LevelSetBooleanResult::getObjectSign() const {
    return m_objectSign;
}

/*!
 * Get the value associated with the results.
 */
double LevelSetBooleanResult::getValue() const {
    return m_value;
}

/*!
	@class      LevelSetBoolean
	@ingroup    levelset
	@brief      Class which deals with boolean operation between two LevelSetObjects
*/

/*!
 * Constructor taking two objects.
 * @param[in] id identifier of object
 * @param[in] op type of boolean operation
 * @param[in] source1 pointer to first source object
 * @param[in] source2 pointer to second source object
 */
LevelSetBoolean::LevelSetBoolean( int id, LevelSetBooleanOperation op, const LevelSetObject *source1, const LevelSetObject *source2  ) :LevelSetProxyObject(id) {
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
LevelSetBoolean::LevelSetBoolean( int id, LevelSetBooleanOperation op, const std::vector<const LevelSetObject*> &sourceObjects ) :LevelSetProxyObject(id) {
    m_operation     = op;
    m_sourceObjects = sourceObjects;
}

/*!
 * Copy constructor.
 * Assigns same id to new object;
 * @param[in] other object to be coppied
 */
LevelSetBoolean::LevelSetBoolean( const LevelSetBoolean &other) :LevelSetProxyObject(other) {
    m_operation     = other.m_operation;
    m_sourceObjects = other.m_sourceObjects;
}

/*!
 * Returns LevelSetInfo 
 * @param[in] id cell id
 * @return LevelSetInfo
*/
LevelSetInfo LevelSetBoolean::getLevelSetInfo( long id)const{
    LevelSetBooleanResult result = computeBooleanResult( id ) ;
    const LevelSetObject *resultObject = result.getObject();
    if ( resultObject ) {
        double value = result.getValue();
        std::array<double, 3> gradient = static_cast<double>(result.getObjectSign()) * resultObject->getGradient( id ) ;

        return LevelSetInfo(value, gradient) ;
    }

    return LevelSetInfo() ;
} 

/*!
 * Get the levelset value
 * @param[in] id cell id
 * @return levelset value in cell
 */
double LevelSetBoolean::getLS( long id)const {
    return getValue(id) ;
}

/*!
 * Get the levelset value
 * @param[in] id cell id
 * @return levelset value in cell
 */
double LevelSetBoolean::getValue( long id)const {
    const LevelSetBooleanResult result = computeBooleanResult( id ) ;
    const LevelSetObject *resultObject = result.getObject();
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
std::array<double,3> LevelSetBoolean::getGradient(long id) const {
    const LevelSetBooleanResult result = computeBooleanResult( id ) ;
    const LevelSetObject *resultObject = result.getObject();
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
LevelSetInfo LevelSetBoolean::computeLevelSetInfo( const std::array<double,3> &coords) const{
    const LevelSetBooleanResult result = computeBooleanResult( coords ) ;
    const LevelSetObject *componentObject = result.getObject();
    if ( componentObject ) {
        LevelSetInfo levelSetInfo = componentObject->computeLevelSetInfo( coords ) ;
        levelSetInfo.value    *= result.getObjectSign();
        levelSetInfo.gradient *= static_cast<double>(result.getObjectSign());

        return levelSetInfo;
    }

    return LevelSetInfo() ;
}

/*!
 * Writes LevelSetBoolean to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetBoolean::_dump( std::ostream &stream ){
    BITPIT_UNUSED(stream);
}

/*!
 * Reads LevelSetBoolean from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetBoolean::_restore( std::istream &stream ){
    BITPIT_UNUSED(stream);
}

/*!
 * Clones the object
 * @return pointer to cloned object
 */
LevelSetBoolean* LevelSetBoolean::clone() const {
    return new LevelSetBoolean( *this ); 
}

/*!
 * Gets the surface normal at the projection point
 * @param[in] id cell index
 * @return closest part
 */
std::array<double,3> LevelSetBoolean::getNormal( long id ) const{

    const LevelSetBooleanResult result = computeBooleanResult(id) ;
    const LevelSetObject *resultObject = result.getObject();
    if ( resultObject ) {
        return (static_cast<double>(result.getObjectSign()) * resultObject->getNormal(id)) ;
    }

    return levelSetDefaults::GRADIENT;
}

/*!
 * Gets the closest part index
 * @param[in] id cell index
 * @return closest part
 */
int LevelSetBoolean::getPart( long id ) const{

    const LevelSetBooleanResult result = computeBooleanResult(id) ;
    const LevelSetObject *resultObject = result.getObject();
    if ( resultObject ) {
        return resultObject->getPart(id) ;
    }

    return levelSetDefaults::PART;
}

/*!
 * Computes the levelset function within the narrow band
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 */
void LevelSetBoolean::computeNarrowBand(bool signd){
    BITPIT_UNUSED(signd) ;
    log::cout() << "Computing levelset within the narrow band... " << std::endl;
}

/*!
 * Updates the levelset function within the narrow band after mesh adaptation.
 * @param[in] adaptionData are the information about the adaption
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 */
void LevelSetBoolean::updateNarrowBand( const std::vector<adaption::Info> &adaptionData, bool signd){
    BITPIT_UNUSED(adaptionData);
    BITPIT_UNUSED(signd);
    log::cout() << "Updating levelset within the narrow band... " << std::endl;
}

/*
 * Returns the boolean operation
 * @return boolean operation
 */
LevelSetBooleanOperation LevelSetBoolean::getBooleanOperation() const{
    return m_operation;
}

/*!
 * Get surface feature size
 * @param[in] id cell index
 * @return charcteristic size
 */
double LevelSetBoolean::getSurfaceFeatureSize( long id ) const {

    const LevelSetBooleanResult result = computeBooleanResult(id) ;
    const LevelSetObject *resultObject = result.getObject();
    if ( resultObject ) {
        return resultObject->getSurfaceFeatureSize(id);
    }

    return (- levelSetDefaults::SIZE);

}

/*!
 * Get the smallest surface feature size
 * @return charcteristic size
 */
double LevelSetBoolean::getMinSurfaceFeatureSize() const {

    bool   minimumValid = false;
    double minimumSize  = levelSetDefaults::SIZE;
    for( const LevelSetObject *object : m_sourceObjects ){
        double objectMinimumSize = object->getMinSurfaceFeatureSize();
        if (objectMinimumSize < 0) {
            continue;
        }

        minimumValid = true;
        minimumSize  = std::min(objectMinimumSize, minimumSize);
    }

    if (!minimumValid) {
        minimumSize = - levelSetDefaults::SIZE;
    }

    return minimumSize;
}

/*!
 * Get the largest surface feature size
 * @return charcteristic size
 */
double LevelSetBoolean::getMaxSurfaceFeatureSize() const {

    double maximumSize = - levelSetDefaults::SIZE;
    for( const LevelSetObject *object : m_sourceObjects ){
        double objectMaximumSize = object->getMaxSurfaceFeatureSize();
        maximumSize = std::max(objectMaximumSize, maximumSize);
    }

    return maximumSize;
}

/*!
 * Get the object that defines the levelset information for the specified cell.
 * @param[in] id cell index
 * @return The object that defines the levelset information for the specified
 * cell.
 */
const LevelSetObject * LevelSetBoolean::getReferenceObject(long id) const{

    const LevelSetBooleanResult result = computeBooleanResult(id) ;

    return result.getObject();

}

/*!
 * Get all objects that compose the boolean object
 * \return pointers to all primary objects involved in the definition of the
 * boolean object levelset information.
 */
std::vector<const LevelSetObject*> LevelSetBoolean::getSourceObjects() const{

    return m_sourceObjects;

}

/*!
 * If cell centroid lies within the narrow band and hence levelset is computet exactly
 * @param[in] id cell id
 * @return true/false if the centroid is in narrow band
 */
bool LevelSetBoolean::isInNarrowBand(long id)const{

    const LevelSetBooleanResult result = computeBooleanResult(id) ;
    const LevelSetObject *resultObject = result.getObject();
    if ( resultObject ) {
        return resultObject->isInNarrowBand(id);
    }

    return false;

}

/*!
 * Compute the result of the boolean operation.
 * @param[in] id cell index
 * @return result of the boolean operation.
 */
LevelSetBooleanResult LevelSetBoolean::computeBooleanResult( long id ) const{

    // Early return if the are no objects
    if (m_sourceObjects.empty()) {
        return LevelSetBooleanResult( getBooleanOperation() );
    }

    // Identify informaiton about the source
    LevelSetBooleanResult result( getBooleanOperation(), m_sourceObjects[0], m_sourceObjects[0]->getValue(id) ) ;
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
LevelSetBooleanResult LevelSetBoolean::computeBooleanResult( const std::array<double,3> &coords ) const{

    // Early return if the are no objects
    if (m_sourceObjects.empty()) {
        return LevelSetBooleanResult( getBooleanOperation() );
    }

    // Identify informaiton about the source
    LevelSetBooleanResult result( getBooleanOperation(), m_sourceObjects[0], m_sourceObjects[0]->computeLevelSetInfo(coords).value);
    for( size_t n=1; n<m_sourceObjects.size(); ++n){
        result.update(m_sourceObjects[n], m_sourceObjects[n]->computeLevelSetInfo( coords ).value );
    }

    return result;
}

}
