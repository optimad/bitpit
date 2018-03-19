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

#include "VTK.hpp"

namespace bitpit{

/*!
 * @class VTKField
 * @ingroup VisualizationToolKit
 * @brief VTKField handles geometry and data field information for the VTK format
 *
 */

/*!
 * Destructor
 */
VTKField::~VTKField(){
}


/*!
 * Default constructor
 */
VTKField::VTKField(){

    m_name            = "undefined" ;  
    m_dataType        = VTKDataType::UNDEFINED ;
    m_location        = VTKLocation::UNDEFINED ;
    m_codification    = VTKFormat::UNDEFINED ;
    m_offset          = 0 ;
    m_fieldType       = VTKFieldType::UNDEFINED ;
    m_position        = 0 ;
    m_streamer        = nullptr ;
    m_enabled         = true ;

}

/*!
 * Copy constructor
 * @param[in] other object to be copied
 */
VTKField::VTKField( const VTKField &other){
    *this = other ;
}

/*!
 * Constructor
 * @param[in] name name of data field
 */
VTKField::VTKField( std::string name ): VTKField() {
    setName(name) ;
}

/*!
 * Assignment operator
 */
VTKField& VTKField::operator=( const VTKField & other){

    m_name  = other.m_name;
    m_fieldType = other.m_fieldType ;
    m_dataType = other.m_dataType ;
    m_codification = other.m_codification;
    m_location = other.m_location ;
    m_offset = other.m_offset ;
    m_position = other.m_position ;
    m_streamer = other.m_streamer ;
    m_enabled = other.m_enabled ;

    return *this;
}

/*!
 * set name of data field
 * @param[in] name name of data field
 */
void VTKField::setName( std::string  name){ 
    m_name= name; 
}

/*!
 * set type of data field
 * @param[in] type type of data [ VTKDataType::[[U]Int[8/16/32/64] / Float[32/64] ] ]
 */
void VTKField::setDataType( VTKDataType type){
    m_dataType= type; 
}

/*!
 * set location of data field
 * @param[in] loc location of data field [VTKLocation::CELL/VTKLocation::POINT]
 */
void VTKField::setLocation( VTKLocation loc ){ 
    m_location= loc; 
}

/*!
 * set codification of data field
 * @param[in] cod codification [VTKFormat::APPENDED/VTKFormat::ASCII]
 */
void VTKField::setCodification( VTKFormat cod ){ 
    m_codification= cod; 
}

/*!
 * set type of data field
 * @param[in] type type of data field [VTKFieldType::SCALAR/VECTOR/KNOWN_BY_CLASS]
 */
void VTKField::setFieldType( VTKFieldType type ){ 
    m_fieldType= type; 
}

/*!
 * set position of data field
 * @param[in] pos position of data field [VTKLocation::CELL/VTKLocation::POINT]
 */
void VTKField::setPosition( std::fstream::pos_type pos ){ 
    m_position =pos ; 
}

/*!
 * set offset of data field for appended output
 * @param[in] off offset from "_" character of appended section
 */
void VTKField::setOffset( uint64_t off){ 
    m_offset= off; 
}

/*!
 * set streamer for reading and writing 
 * @param[in] streamer streamer 
 */
void VTKField::setStreamer(VTKBaseStreamer& streamer ){ 
    m_streamer= &streamer; 
}

/*!
 * Enables the field for writing/reading
 */
void VTKField::enable(){
    m_enabled=true;
}

/*!
 * Disables the field for writing/reading
 */
void VTKField::disable(){
    m_enabled=false;
}

/*!
 * get name of data field
 * @return  name of data field
 */
std::string VTKField::getName() const{ 
    return m_name; 
}

/*!
 * get type of field
 * @return type of data field [ VTKFieldType ]
 */
VTKFieldType VTKField::getFieldType() const{ 
    return m_fieldType; 
}

/*!
 * get type of data field
 * @return type of data field [ VTKDataType ]
 */
VTKDataType VTKField::getDataType() const{ 
    return m_dataType; 
}

/*!
 * get location of data field
 * @return location [VTKLocation::CELL/VTKLocation::POINT]
 */
VTKLocation VTKField::getLocation() const{ 
    return m_location; 
}

/*!
 * get codification of data field
 * @return codification [VTKFormat::APPENDED/VTKFormat::ASCII]
 */
VTKFormat VTKField::getCodification() const{ 
    return m_codification; 
}

/*!
 * get offset in appended section
 * @return offset from "_" character in appended section
 */
uint64_t VTKField::getOffset() const{ 
    return m_offset; 
}

/*!
 * get position of data field in VTK file. 
 * This information is available after VTK::ReadMetaData() has been called.
 * @return position in VTK file.
 */
std::fstream::pos_type   VTKField::getPosition() const{ 
    return m_position; 
}

/*!
 * get streamer of data field
 * @return streamer of data field
 */
const VTKBaseStreamer & VTKField::getStreamer() const {
    return *m_streamer;
}

/*!
 * Returns if field is enabled for readind/writing
 * @return true if enabled
 */
bool VTKField::isEnabled() const{ 
    return m_enabled; 
}

/*!
 * Check if all information of field has been set
 * @return true if all information regarding field is available
 */
bool VTKField::hasAllMetaData() const{ 

    bool allData(true);

    allData &= m_name != "undefined" ;
    allData &= m_dataType != VTKDataType::UNDEFINED ;
    allData &= m_location != VTKLocation::UNDEFINED ;
    allData &= m_codification != VTKFormat::UNDEFINED ;
    allData &= m_fieldType != VTKFieldType::UNDEFINED ;
    allData &= m_streamer != nullptr ;

    return allData;

}

/*!
 * Writes the field through its streamer to file
 * @param[in] str file stream
 */
void  VTKField::write( std::fstream &str) const{ 
    if(m_streamer){
        m_streamer->flushData( str, m_name, m_codification) ;
    }
}

/*!
 * Reads the field through its streamer from file
 * @param[in] str file stream
 * @param[in] entries total number of entries to be read
 * @param[in] components size of subgroup
 */
void VTKField::read( std::fstream &str, uint64_t entries, uint8_t components ) const{ 
    if(m_streamer){
        m_streamer->absorbData( str, m_name, m_codification, entries, components, m_dataType) ;
    }
}

}
