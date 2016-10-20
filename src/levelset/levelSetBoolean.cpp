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

# include <cassert>

# include "levelSet.hpp"

# include "bitpit_common.hpp"
# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"

namespace bitpit {

/*!
	@ingroup    levelset
	@class      LevelSetBoolean
	@brief      Class which deals with boolean operation between two LevelSetObjects
*/

/*!
 * Destructor
 */
LevelSetBoolean::~LevelSetBoolean(){
    m_objPtr1 = nullptr;
    m_objPtr2 = nullptr;
};

/*!
 * Constructor
 * @param[in] id identifier of object
 * @param[in] angle feature angle; if the angle between two segments is bigger than this angle, the enclosed edge is considered as a sharp edge.
 */
LevelSetBoolean::LevelSetBoolean( int id, LevelSetBooleanOperation op, LevelSetObject *ptr1, LevelSetObject *ptr2  ) :LevelSetObject(id,false) {
    m_operation = op;
    m_objPtr1 = ptr1;
    m_objPtr2 = ptr2;
    m_objId1 = ptr1->getId() ;
    m_objId2 = ptr2->getId() ;
};

/*!
 * Copy constructor.
 * Assigns same id to new object;
 * @param[in] other object to be coppied
 */
LevelSetBoolean::LevelSetBoolean( const LevelSetBoolean &other) :LevelSetObject(other.getId(),false) {
    m_operation = other.m_operation;
    m_objPtr1 = other.m_objPtr1;
    m_objPtr2 = other.m_objPtr2;
    m_objId1 = other.m_objId1 ;
    m_objId2 = other.m_objId2 ;
};

/*!
 * Returns reference to LevelSetInfo
*/
LevelSetInfo LevelSetBoolean::getLevelSetInfo( const long &i)const{
    return booleanOperation(i) ;
} 

/*!
 * Get the Sdf value of the i-th local element of the octree mesh.
 * @param[in] i cell index
 * @return levelset value in cell
 */
double LevelSetBoolean::getLS( const long &i)const {
    return booleanOperation(i).value ;
};

/*!
 * Get the Sdf gradient vector of the i-th local element of the octree mesh.
 * @param[in] i cell index
 * @return levelset gradient in cell 
 */
std::array<double,3> LevelSetBoolean::getGradient(const long &i) const {
    return booleanOperation(i).gradient ;
};

/*!
 * Writes LevelSetBoolean to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetBoolean::_dump( std::ostream &stream ){

    IO::binary::write( stream, m_objId1 ) ;
    IO::binary::write( stream, m_objId2 ) ;
};

/*!
 * Reads LevelSetBoolean from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetBoolean::_restore( std::istream &stream ){

    IO::binary::read( stream, m_objId1 ) ;
    IO::binary::read( stream, m_objId2 ) ;
};

/*!
 * Clones the object
 * @return pointer to cloned object
 */
LevelSetBoolean* LevelSetBoolean::clone() const {
    return new LevelSetBoolean( *this ); 
}

/*!
 * Gets the closest support within the narrow band of cell
 * @param[in] id index of cell
 * @return closest segment in narrow band
 */
int LevelSetBoolean::getPart( const long &id ) const{
    LevelSetObject *objPtr = getClosestObject(id) ;
    return objPtr->getPart(id) ;
};

/*!
 * Gets the closest support within the narrow band of cell
 * @param[in] id index of cell
 * @return closest segment in narrow band
 */
long LevelSetBoolean::getSupport( const long &id ) const{
    LevelSetObject *objPtr = getClosestObject(id) ;
    return objPtr->getSupport(id) ;
};

/*!
 * Gets the number of support items within the narrow band of cell
 * @param[in] id index of cell
 * @return number of segments in narrow band 
 */
int LevelSetBoolean::getSupportCount( const long &id ) const{
    return m_objPtr1->getSupportCount(id)+m_objPtr2->getSupportCount(id);
};

/*!
 * Manually set the size of the narrow band.
 * @param[in] r size of the narrow band.
 */
void LevelSetBoolean::setSizeNarrowBand(double r){
    m_RSearch = r;
    m_objPtr1->setSizeNarrowBand(r);
    m_objPtr2->setSizeNarrowBand(r);
};

/*!
 * Computes the size of the narrow band
 * @return size of the narrow band.
 */
double LevelSetBoolean::computeSizeNarrowBand(LevelSetKernel* visitee){
    BITPIT_UNUSED(visitee);
    return std::max( m_objPtr1->getSizeNarrowBand(), m_objPtr2->getSizeNarrowBand() );
};

/*!
 * Computes the size of the narrow band
 * @return size of the narrow band.
 */
double LevelSetBoolean::updateSizeNarrowBand(LevelSetKernel* visitee, const std::vector<adaption::Info> &mapper){
    BITPIT_UNUSED(visitee);
    BITPIT_UNUSED(mapper);
    return std::max( m_objPtr1->getSizeNarrowBand(), m_objPtr2->getSizeNarrowBand() );
};

/*!
 * Computes the levelset function within the narrow band
 * @param[in] visitee pointer to mesh
 * @param[in] RSearch size of narrow band
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 */
void LevelSetBoolean::computeLSInNarrowBand( LevelSetKernel *visitee, const double &RSearch, const bool &signd ){

    BITPIT_UNUSED(visitee) ;
    BITPIT_UNUSED(RSearch) ;
    BITPIT_UNUSED(signd) ;

    log::cout() << "Computing levelset within the narrow band... " << std::endl;
};

/*!
 * Updates the levelset function within the narrow band after mesh adaptation.
 * @param[in] visitee pointer to mesh
 * @param[in] mapper information concerning mesh adaption 
 * @param[in] RSearch size of narrow band
 * @param[in] signd if signed- or unsigned- distance function should be calculated
 */
void LevelSetBoolean::updateLSInNarrowBand( LevelSetKernel *visitee, const std::vector<adaption::Info> &mapper, const double &RSearch, const bool &signd ){

    BITPIT_UNUSED(visitee);
    BITPIT_UNUSED(mapper);
    BITPIT_UNUSED(RSearch);
    BITPIT_UNUSED(signd);

    log::cout() << "Updating levelset within the narrow band... " << std::endl;
};

/*
 * Returns the boolean operation
 * @return boolean operation
 */
LevelSetBooleanOperation LevelSetBoolean::getBooleanOperation() const{
    return m_operation;
}

/*
 * Determines the closest object
 * @param[in] id id of cell
 * @return pointer to closest LevelSetObject
 */
LevelSetObject* LevelSetBoolean::getClosestObject( const long &id) const{

    double value1 = m_objPtr1->getLS(id) ;
    double value2 = m_objPtr2->getLS(id) ;
    
    if( getBooleanOperation() == LevelSetBooleanOperation::UNION){
        return (value1<=value2) ? m_objPtr1 : m_objPtr2 ;

    } else if ( getBooleanOperation() == LevelSetBooleanOperation::INTERSECTION){
        return (value1>=value2) ? m_objPtr1 : m_objPtr2 ;

    } else if ( getBooleanOperation() == LevelSetBooleanOperation::SUBTRACTION){
        return (value1>=-1.*value2) ? m_objPtr1 : m_objPtr2 ;

    }

    return nullptr;
}

/*
 * Performs the bolean operation
 * @param[in] id id of cell
 * @return resulting levelset value and gradient in LevelSetInfo
 */
LevelSetInfo LevelSetBoolean::booleanOperation(const long &id) const{

    LevelSetInfo info1 = m_objPtr1->getLevelSetInfo(id); 
    LevelSetInfo info2 = m_objPtr2->getLevelSetInfo(id); 

    if( getBooleanOperation() == LevelSetBooleanOperation::UNION){
        return (info1.value<=info2.value) ? info1 : info2 ;

    } else if ( getBooleanOperation() == LevelSetBooleanOperation::INTERSECTION){
        return (info1.value>=info2.value) ? info1 : info2 ;
    
    } else if ( getBooleanOperation() == LevelSetBooleanOperation::SUBTRACTION){
        return (info1.value>=-1.*info2.value) ? info1 : LevelSetInfo(-1.*info2.value,-1.*-1.*info2.gradient) ;

    }

    assert(false);
    return LevelSetInfo() ;
}

}
