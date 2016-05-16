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

#include "VTK.hpp"

namespace bitpit{
/*!
 * @ingroup VisualizationToolKit
 * @{
 * @class        VTKField
 * @brief        VTKField handles geometry and data field information for the VTK format
 *
 */

/*!
 * Destructor
 */
VTKField::~VTKField(){
};


/*!
 * Default constructor
 */
VTKField::VTKField(){

    name            = "undefined" ;  
    dataType        = VTKDataType::UNDEFINED ;
    location        = VTKLocation::UNDEFINED ;
    codification    = VTKFormat::UNDEFINED ;
    fieldType       = VTKFieldType::UNDEFINED ;
    components      = 0 ;
    nr_elements     = 0 ;
    position        = 0 ;

    derived         = false;
    implicitKnown   = false;

};

/*!
 * Copy constructor
 * @param[in]   other   object to be copied
 */
VTKField::VTKField( const VTKField &other){

    *this = other ;
};

/*!
 * Constructor
 * @param[in]   name_   name of data field
 */
VTKField::VTKField( std::string name_ ): VTKField() {

    name = name_ ;
};

/*!
 * Assignment operator
 */
VTKField& VTKField::operator=( const VTKField & other){

    name  = other.name;
    fieldType = other.fieldType ;
    components = other.components ;
    dataType = other.dataType ;
    codification = other.codification;
    location = other.location ;
    nr_elements = other.nr_elements ;
    offset = other.offset ;
    position = other.position ;
    derived = other.derived ;
    implicitKnown = other.implicitKnown ;

    return *this;
};

/*!
 * set name of data field
 * @param[in]   name_   name of data field
 */
void      VTKField::setName( std::string  name_){ 
    name= name_; 
    return; 
};

/*!
 * set type of data field
 * @param[in]  type_    type of data [ VTKDataType::[[U]Int[8/16/32/64] / Float[32/64] ] ]
 */
void      VTKField::setDataType( VTKDataType  type_){
    dataType= type_; 
    return; 
};

/*!
 * set location of data field
 * @param[in]   loc_   location of data field [VTKLocation::CELL/VTKLocation::POINT]
 */
void      VTKField::setLocation( VTKLocation  loc_ ){ 
    location= loc_; 
    return; 
};

/*!
 * set codification of data field
 * @param[in]   code_  codification [VTKFormat::APPENDED/VTKFormat::ASCII]
 */
void      VTKField::setCodification( VTKFormat  code_ ){ 
    codification= code_; 
    return; 
};

/*!
 * set type of data field
 * @param[in]   type_   type of data field [VTKFieldType::SCALAR/VECTOR/CONSTANT/VARIABLE]
 */
void      VTKField::setFieldType( VTKFieldType type_){ 
    fieldType= type_; 
    return; 
};

/*!
 * set number of components of data field to be used if VTKFieldType::CONSTANT
 * @param[in]   comp_   number of coponents of data field 
 */
void      VTKField::setComponents( uint8_t comp_){ 
    components= comp_; 
    return; 
};

/*!
 * set numer of elements of data field
 * @param[in]   ele_   number of elements
 */
void      VTKField::setElements( uint64_t ele_){ 
    nr_elements= ele_; 
    return; 
};

/*!
 * set position of data field
 * @param[in]   pos_    position of data field [VTKLocation::CELL/VTKLocation::POINT]
 */
void      VTKField::setPosition( std::fstream::pos_type pos_ ){ 
    position =pos_ ; 
    return; 
} ;

/*!
 * set offset of data field for appended output
 * @param[in]   offs_   offset from "_" character of appended section
 */
void      VTKField::setOffset( uint64_t offs_){ 
    offset= offs_; 
    return; 
};

/*!
 * set if class storing the field is aware about data to be written
 * @param[in]   known   true if aware
 */
void      VTKField::setImplicit( bool known ){ 
    implicitKnown= known; 
    return; 
};

/*!
 * get name of data field
 * @return  name of data field
 */
std::string    VTKField::getName() const{ 
    return name; 
};

/*!
 * get type of field
 * @return   type of data field [ VTKFieldType ]
 */
VTKFieldType    VTKField::getFieldType() const{ 
    return fieldType; 
};

/*!
 * get type of data field
 * @return   type of data field [ VTKDataType ]
 */
VTKDataType    VTKField::getDataType() const{ 
    return dataType; 
};

/*!
 * get location of data field
 * @return  location [VTKLocation::CELL/VTKLocation::POINT]
 */
VTKLocation    VTKField::getLocation() const{ 
    return location; 
};

/*!
 * get codification of data field
 * @return  codification [VTKFormat::APPENDED/VTKFormat::ASCII]
 */
VTKFormat    VTKField::getCodification() const{ 
    return codification; 
};

/*!
 * get number of components of data field
 * @return  number of components 
 */
uint8_t   VTKField::getComponents() const{ 
    return components; 
};

/*!
 * get number of elements of data field. 
 * Generally number of points or number of cells.
 * @return  number of elements
 */
uint64_t  VTKField::getElements() const{ 
    return nr_elements; 
};

/*!
 * get total number of entries data field
 * @return  total size = number of elements *nuber of components
 */
uint64_t  VTKField::getSize() const{ 

    if( fieldType == VTKFieldType::SCALAR || fieldType == VTKFieldType::VARIABLE ){
        return nr_elements ; 

    } else if( fieldType == VTKFieldType::VECTOR){
        return static_cast<int>(fieldType) *nr_elements ; 

    } else if( fieldType == VTKFieldType::CONSTANT ){
        return components *nr_elements ; 
    };

    assert(false);

};

/*!
 * get offset in appended section
 * @return  offset from "_" character in appended section
 */
uint64_t  VTKField::getOffset() const{ 
    return offset; 
};

/*!
 * get bytes of  data field
 * @return  memory size of data field
 */
uint64_t  VTKField::getNbytes() const{ 
    return getSize() *VTKTypes::sizeOfType( dataType ) ;
};

/*!
 * get position of data field in VTK file. 
 * This information is available after VTK::ReadMetaData() has been called.
 * @return      position in VTK file.
 */
std::fstream::pos_type   VTKField::getPosition() const{ 
    return position; 
};

/*!
 * get position of data field in VTK file. 
 * This information is available after VTK::ReadMetaData() has been called.
 * @return      position in VTK file.
 */
bool   VTKField::hasAllMetaData() const{ 

    bool    allData(true);

    allData = allData && name != "undefined" ;  
    allData = allData && dataType != VTKDataType::UNDEFINED ;
    allData = allData && location != VTKLocation::UNDEFINED ;
    allData = allData && codification != VTKFormat::UNDEFINED ;
    allData = allData && fieldType != VTKFieldType::UNDEFINED ;   ;
    allData = allData && nr_elements != 0 ;

    return allData;
};

/*!
 * Checks which mechanism is used to read/write data
 * @return  true if flushData is used, false if internal pointer to data is used
 */
bool VTKField::usesInterface( )const{ 

    return( !autoWrite() && !implicitKnown ) ;
};

/*!
 * Checks which mechanism is used to read/write data
 * @return  true if a derived class with explicit data storage is used
 */
bool VTKField::autoWrite( )const{ 

    return( derived ) ;
};

/*!
 * Checks which mechanism is used to read/write data
 * @return  true if flushData is used, false if internal pointer to data is used
 */
void VTKField::flushData( std::fstream &str )const{ 
    BITPIT_UNUSED( str);

    return  ;
};

/*!
 * Checks which mechanism is used to read/write data
 * @return  true if flushData is used, false if internal pointer to data is used
 */
void VTKField::absorbData( std::fstream &str )const{ 
    BITPIT_UNUSED( str);

    return  ;
};
/*!
 * @}
 */

}
