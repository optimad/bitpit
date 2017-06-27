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

# include <cassert>

# include "bitpit_common.hpp"
# include "bitpit_operators.hpp"

# include "logger.hpp"
# include "adaption.hpp"

# include "levelSetObject.hpp"
# include "levelSetMetaObject.hpp"
# include "levelSetBoolean.hpp"

namespace bitpit {

/*!
	@class      LevelSetBoolean
	@ingroup    levelset
	@brief      Class which deals with boolean operation between two LevelSetObjects
*/

/*!
 * Destructor
 */
LevelSetBoolean::~LevelSetBoolean(){
    m_objPtr.clear();
}

/*!
 * Constructor taking two objects.
 * @param[in] id identifier of object
 * @param[in] op type of boolean operation
 * @param[in] ptr1 pointer to first object
 * @param[in] ptr2 pointer to second object
 */
LevelSetBoolean::LevelSetBoolean( int id, LevelSetBooleanOperation op, LevelSetObject *ptr1, LevelSetObject *ptr2  ) :LevelSetMetaObject(id) {
    m_operation = op;
    m_objPtr.push_back(ptr1);
    m_objPtr.push_back(ptr2);
}

/*!
 * Constructor taking a vector of objects.
 * The boolean operation will be applied recursivly on each entry.
 * @param[in] id identifier of object
 * @param[in] op type of boolean operation
 * @param[in] objPtr vector of pointers to objects
 */
LevelSetBoolean::LevelSetBoolean( int id, LevelSetBooleanOperation op, std::vector<LevelSetObject*> objPtr ) :LevelSetMetaObject(id) {
    m_operation = op;
    m_objPtr = objPtr;
}

/*!
 * Copy constructor.
 * Assigns same id to new object;
 * @param[in] other object to be coppied
 */
LevelSetBoolean::LevelSetBoolean( const LevelSetBoolean &other) :LevelSetMetaObject(other.getId()) {
    m_operation = other.m_operation;
    m_objPtr = other.m_objPtr;
}

/*!
 * Returns LevelSetInfo 
 * @param[in] i cell index
 * @return LevelSetInfo
*/
LevelSetInfo LevelSetBoolean::getLevelSetInfo( const long &i)const{
    return booleanOperation(i) ;
} 

/*!
 * Get the levelset value
 * @param[in] i cell index
 * @return levelset value in cell
 */
double LevelSetBoolean::getLS( const long &i)const {
    return booleanOperation(i).value ;
}

/*!
 * Get the levelset gradient
 * @param[in] i cell index
 * @return levelset gradient in cell 
 */
std::array<double,3> LevelSetBoolean::getGradient(const long &i) const {
    return booleanOperation(i).gradient ;
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
 * Gets the closest part index
 * @param[in] id cell index
 * @return closest part
 */
int LevelSetBoolean::getPart( const long &id ) const{
    LevelSetObject *objPtr = getCompetentObject(id) ;

    if( objPtr == nullptr){
        return levelSetDefaults::PART;
    }

    return objPtr->getPart(id) ;
}

/*!
 * Manually set the size of the narrow band.
 * @param[in] r size of the narrow band.
 */
void LevelSetBoolean::setSizeNarrowBand(double r){
    m_RSearch = r;
    for( auto objPtr : m_objPtr){
        objPtr->setSizeNarrowBand(r);
    }
}

/*!
 * Computes the size of the narrow band
 * @return size of the narrow band.
 */
double LevelSetBoolean::computeSizeNarrowBand(){
    double r = -std::numeric_limits<double>::max() ;
    for( auto objPtr : m_objPtr){
        r = std::max(r, objPtr->getSizeNarrowBand());
    }

    return r;
}

/*!
 * Updates the size of the narrow band after mesh has been modified
 * @param[in] mapper descriptor of mesh modifications
 * @return size of the narrow band.
 */
double LevelSetBoolean::updateSizeNarrowBand(const std::vector<adaption::Info> &mapper){
    BITPIT_UNUSED(mapper);
    return computeSizeNarrowBand();
}

/*!
 * Computes the levelset function within the narrow band
 * @param[in] RSearch size of narrow band
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 */
void LevelSetBoolean::computeLSInNarrowBand( const double &RSearch, const bool &signd ){

    BITPIT_UNUSED(RSearch) ;
    BITPIT_UNUSED(signd) ;

    log::cout() << "Computing levelset within the narrow band... " << std::endl;
}

/*!
 * Updates the levelset function within the narrow band after mesh adaptation.
 * @param[in] mapper information concerning mesh adaption 
 * @param[in] RSearch size of narrow band
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 */
void LevelSetBoolean::updateLSInNarrowBand( const std::vector<adaption::Info> &mapper, const double &RSearch, const bool &signd ){

    BITPIT_UNUSED(mapper);
    BITPIT_UNUSED(RSearch);
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
double LevelSetBoolean::getSurfaceFeatureSize( const long &id ) const {

    LevelSetObject *objectPtr = getCompetentObject(id);
    if (!objectPtr) {
        return (- levelSetDefaults::SIZE);
    }

    return objectPtr->getSurfaceFeatureSize(id);
}

/*!
 * Get the smallest surface feature size
 * @return charcteristic size
 */
double LevelSetBoolean::getMinSurfaceFeatureSize() const {

    bool   minimumValid = false;
    double minimumSize  = levelSetDefaults::SIZE;
    for( const auto & object : m_objPtr ){
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
    for( const auto & object : m_objPtr ){
        double objectMaximumSize = object->getMaxSurfaceFeatureSize();
        maximumSize = std::max(objectMaximumSize, maximumSize);
    }

    return maximumSize;
}

/*
 * Determines the relevant object which determines the levelset value in the cell
 * Taken from http://www.iue.tuwien.ac.at/phd/ertl/node57.html
 * @param[in] id cell index
 * @return pointer to competent LevelSetObject
 */
LevelSetObject* LevelSetBoolean::getCompetentObject( const long &id) const{

    if(m_objPtr.size()==0){ 
        return nullptr;
    }

    double result, second;
    LevelSetObject *resPtr, *secPtr;

    result = m_objPtr[0]->getLS(id);
    resPtr = m_objPtr[0];

    for( size_t n=1; n<m_objPtr.size(); ++n){
        second = m_objPtr[n]->getLS(id) ;
        secPtr = m_objPtr[n];

        if( getBooleanOperation() == LevelSetBooleanOperation::UNION){
            if(result>second) {
                result = second;
                resPtr = secPtr;
            }

        } else if ( getBooleanOperation() == LevelSetBooleanOperation::INTERSECTION){
            if(result<second) {
                result = second;
                resPtr = secPtr;
            }

        } else if ( getBooleanOperation() == LevelSetBooleanOperation::SUBTRACTION){
            if(result<-1.*second) {
                result = -1.*second;
                resPtr = secPtr;
            }
        }
    }

    return resPtr;
}

/*!
 * Performs the bolean operation
 * Taken from http://www.iue.tuwien.ac.at/phd/ertl/node57.html
 * @param[in] id cell index
 * @return resulting levelset value and gradient in LevelSetInfo
 */
LevelSetInfo LevelSetBoolean::booleanOperation(const long &id) const{

    if(m_objPtr.size()==0){ 
        return LevelSetInfo();
    }

    LevelSetInfo result, second;
    result = m_objPtr[0]->getLevelSetInfo(id);

    for( size_t n=1; n<m_objPtr.size(); ++n){
        second = m_objPtr[n]->getLevelSetInfo(id) ;

        if( getBooleanOperation() == LevelSetBooleanOperation::UNION){
            if(result.value>second.value) {
                result = second;
            }

        } else if ( getBooleanOperation() == LevelSetBooleanOperation::INTERSECTION){
            if(result.value<second.value) {
                result = second;
            }

        } else if ( getBooleanOperation() == LevelSetBooleanOperation::SUBTRACTION){
            if(result.value<-1.*second.value) {
                second.value *= -1.;
                second.gradient *= -1.;
                result = second;
            }
        }
    }

    return result;
}

/*!
 * Get the index of the primary object
 * @param[in] cellId cell index
 * @return primary object
 */
int LevelSetBoolean::getPrimaryObjectId(const long &cellId) const{

    LevelSetObject *obj = getCompetentObject(cellId);

    if( LevelSetMetaObject *meta = dynamic_cast<LevelSetMetaObject*>(obj) ){
        return meta->getPrimaryObjectId(cellId);
    }

    return obj->getId();
    
}

}
