#include "VTK.hpp"

/*!
 * @class        VTKField
 * @brief        VTKField handles geometry and data field information for the VTK format
 *
 */

/*!
 * Default constructor
 */
VTKField::VTKField(){

    name            = "undefined" ;  
    type            = VTKDataType::UNDEFINED ;
    location        = VTKLocation::UNDEFINED ;
    codification    = VTKFormat::UNDEFINED ;
    components      = VTKFieldType::UNDEFINED ;   ;
    nr_elements     = 0 ;
    position        = 0 ;

};

/*!
 * Constructor
 * @param[in]   name_   name of data field
 * @param[in]   comp_   number of coponents of data field [VTKFieldType::SCALAR/VTKFieldType::VECTOR]
 * @param[in]   loc_    location of data field [VTKLocation::CELL/VTKLocation::POINT]
 */
VTKField::VTKField( std::string name_, VTKFieldType comp_, VTKLocation loc_ ): name(name_), components(comp_), location(loc_) {

    type            = VTKDataType::UNDEFINED ;
    nr_elements     = 0 ;
    codification    = VTKFormat::UNDEFINED ;
    position        = 0 ;

};

/*!
 * Constructor
 * @param[in]   name_   name of data field
 * @param[in]   comp_   number of coponents of data field [VTKFieldType::SCALAR/VTKFieldType::VECTOR]
 * @param[in]   loc_    location of data field [VTKLocation::CELL/VTKLocation::POINT]
 * @param[in]   type_   type of data field [ VTKDataType::[[U]Int[8/16/32/64] / Float[32/64] ] ]
 */
VTKField::VTKField( std::string name_, VTKFieldType comp_, VTKLocation loc_, VTKDataType type_ ): name(name_), components(comp_), type(type_), location(loc_) {

    nr_elements     = 0 ;
    codification    = VTKFormat::UNDEFINED ;
    position        = 0 ;

};

/*!
 * Constructor
 * @param[in]   name_   name of data field
 * @param[in]   comp_   number of coponents of data field [VTKFieldType::SCALAR/VTKFieldType::VECTOR]
 * @param[in]   loc_    location of data field [VTKLocation::CELL/VTKLocation::POINT]
 * @param[in]   type_   type of data field [ VTKDataType::[[U]Int[8/16/32/64] / Float[32/64] ] ]
 * @param[in]   cod_    codex [VTKFormat::APPENDED/VTKFormat::ASCII]
 * @param[in]   nr_elements_    number of elements of data field
 */
VTKField::VTKField( std::string name_, VTKFieldType comp_, VTKLocation loc_, VTKDataType type_, VTKFormat cod_, uint64_t nr_elements_):
    name(name_), components(comp_), type(type_), location(loc_), codification(cod_), nr_elements(nr_elements_){
};

/*!
 * Destructor
 */
VTKField::~VTKField(){
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
void      VTKField::setType( VTKDataType  type_){
    type= type_; 
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
 * set number of components of data field
 * @param[in]   comp_   number of coponents of data field [VTKFieldType::SCALAR/VTKFieldType::VECTOR]
 */
void      VTKField::setComponents( VTKFieldType comp_){ 
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
 * get name of data field
 * @return  name of data field
 */
std::string    VTKField::getName() const{ 
    return name; 
};

/*!
 * get type of data field
 * @return   type of data field [ "[U]Int8", "[U]Int16", "[U]Int32", "[U]Int64", "Float32", "Float64" ]
 */
VTKDataType    VTKField::getType() const{ 
    return type; 
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
 * @return  number of components [1/3]
 */
VTKFieldType   VTKField::getComponents() const{ 
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

    return static_cast<int>(components) *nr_elements ; 
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
    return static_cast<int>(components) *nr_elements *VTKTypes::sizeOfType( type ) ;
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
    allData = allData && type != VTKDataType::UNDEFINED ;
    allData = allData && location != VTKLocation::UNDEFINED ;
    allData = allData && codification != VTKFormat::UNDEFINED ;
    allData = allData && components != VTKFieldType::UNDEFINED ;   ;
    allData = allData && nr_elements != 0 ;

    return allData;
};


/*!
 * Update field information using VTKFieldMetaData class
 */
void VTKField::importMetaData( const VTKFieldMetaData &data){ 

    setType( VTKTypes::whichType( data.getType() ) );

    if( getComponents() == VTKFieldType::UNDEFINED)
        setComponents( VTKFieldType::SCALAR);
    setElements( data.getSize() / static_cast<int>(components)  );

    return ;
};


/*!
 * @}
 */
