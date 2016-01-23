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
    type            = "undefined" ;
    location        = "undefined" ;
    codification    = "undefined" ;
    components      = 1 ;   ;
    nr_elements     = 0 ;
    position        = 0 ;

};

/*!
 * Constructor
 * @param[in]   name_   name of data field
 * @param[in]   comp_   number of coponents of data field [1/3]
 * @param[in]   loc_    location of data field ["Cell"/"Point"]
 */
VTKField::VTKField( std::string name_, uint8_t comp_, std::string loc_ ): name(name_), components(comp_), location(loc_) {

    type            = "undefined" ;
    nr_elements     = 0 ;
    codification    = "undefined" ;
    position        = 0 ;

};

/*!
 * Constructor
 * @param[in]   name_   name of data field
 * @param[in]   comp_   number of coponents of data field [1/3]
 * @param[in]   type_   type of data field [ "[U]Int8", "[U]Int16", "[U]Int32", "[U]Int64", "Float32", "Float64" ]
 * @param[in]   loc_    location of data field ["Cell"/"Point"]
 */
VTKField::VTKField( std::string name_, uint8_t comp_, std::string type_, std::string loc_ ): name(name_), components(comp_), type(type_), location(loc_) {

    nr_elements     = 0 ;
    codification    = "undefined" ;
    position        = 0 ;

};

/*!
 * Constructor
 * @param[in]   name_   name of data field
 * @param[in]   comp_   number of coponents of data field [1/3]
 * @param[in]   type_   type of data field [ "[U]Int8", "[U]Int16", "[U]Int32", "[U]Int64", "Float32", "Float64" ]
 * @param[in]   loc_    location of data field ["Cell"/"Point"]
 * @param[in]   cod_    codex ["ascii"/ "appended"]
 * @param[in]   nr_elements_    number of elements of data field
 */
VTKField::VTKField( std::string name_, uint8_t comp_, std::string type_, std::string loc_, std::string cod_, uint64_t nr_elements_):
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
 * @param[in]   type_   type of data field [ "[U]Int8", "[U]Int16", "[U]Int32", "[U]Int64", "Float32", "Float64" ]
 */
void      VTKField::setType( std::string  type_){

    if( type_ == "Int8"    || 
            type_ == "UInt8"   ||
            type_ == "Int16"   || 
            type_ == "UInt16"  ||
            type_ == "Int32"   || 
            type_ == "UInt32"  ||
            type_ == "Int64"   || 
            type_ == "UInt64"  ||
            type_ == "Float32" ||
            type_ == "Float64" ){

        type= type_; 

    }

    else{

        std::cout << type_ << " not admitted in VTKField::setType." << std::endl ;
    };


    return; 
};

/*!
 * set location of data field
 * @param[in]   loc_   location of data field [ "Cell", "Point" ] 
 */
void      VTKField::setLocation( std::string  loc_ ){ 
    location= loc_; 
    return; 
};

/*!
 * set codification of data field
 * @param[in]   code_  codification ["appended"/ "ascii"]
 */
void      VTKField::setCodification( std::string  code_ ){ 
    codification= code_; 
    return; 
};

/*!
 * set number of components of data field
 * @param[in]   comp_   number of components [1/3]
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
 * @param[in]   pos_    position of data field ["Cell"/"Point"]
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
std::string    VTKField::getName(){ 
    return name; 
};

/*!
 * get type of data field
 * @return   type of data field [ "[U]Int8", "[U]Int16", "[U]Int32", "[U]Int64", "Float32", "Float64" ]
 */
std::string    VTKField::getType(){ 
    return type; 
};

/*!
 * get location of data field
 * @return  location ["Cell"/"Point"]
 */
std::string    VTKField::getLocation(){ 
    return location; 
};

/*!
 * get codification of data field
 * @return  codification ["ascii"/"appended"]
 */
std::string    VTKField::getCodification(){ 
    return codification; 
};

/*!
 * get number of components of data field
 * @return  number of components [1/3]
 */
uint8_t   VTKField::getComponents(){ 
    return components; 
};

/*!
 * get number of elements of data field. 
 * Generally number of points or number of cells.
 * @return  number of elements
 */
uint64_t  VTKField::getElements(){ 
    return nr_elements; 
};

/*!
 * get total number of entries data field
 * @return  total size = number of elements *nuber of components
 */
uint64_t  VTKField::getSize(){ 
    return components *nr_elements ; 
};

/*!
 * get offset in appended section
 * @return  offset for "_" character in appended section
 */
uint64_t  VTKField::getOffset(){ 
    return offset; 
};

/*!
 * get bytes of  data field
 * @return  memory size of data field
 */
uint64_t  VTKField::getNbytes(){ 
    return components *nr_elements *VTKUtils::sizeOfType( type ) ;
};

/*!
 * get position of data field in VTK file. 
 * This information is available after VTK::ReadMetaData() has been called.
 * @return      position in VTK file.
 */
std::fstream::pos_type   VTKField::getPosition(){ 
    return position; 
};



/*!
 * @}
 */
