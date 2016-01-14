#include"Class_VTK.hpp"

/*!
 * @class        VTK::Field_C
 * @brief        Field_C handles geometry and data field information for the VTK format
 *
 */

/*! ------------------------------------------------------------------
 * Default constructor
 */
VTK::Field_C::Field_C(){

    name            = "undefined" ;  
    type            = "undefined" ;
    location        = "undefined" ;
    codification    = "undefined" ;
    components      = 1 ;   ;
    nr_elements     = 0 ;
    position        = 0 ;

};

/*! ------------------------------------------------------------------
 * Constructor
 * @param[in]   name_   name of data field
 * @param[in]   comp_   number of coponents of data field [1/3]
 * @param[in]   loc_    location of data field ["Cell"/"Point"]
 */
VTK::Field_C::Field_C( std::string name_, uint8_t comp_, std::string loc_ ): name(name_), components(comp_), location(loc_) {

    type            = "undefined" ;
    nr_elements     = 0 ;
    codification    = "undefined" ;
    position        = 0 ;

};

/*! ------------------------------------------------------------------
 * Constructor
 * @param[in]   name_   name of data field
 * @param[in]   comp_   number of coponents of data field [1/3]
 * @param[in]   type_   type of data field [ "[U]Int8", "[U]Int16", "[U]Int32", "[U]Int64", "Float32", "Float64" ]
 * @param[in]   loc_    location of data field ["Cell"/"Point"]
 */
VTK::Field_C::Field_C( std::string name_, uint8_t comp_, std::string type_, std::string loc_ ): name(name_), components(comp_), type(type_), location(loc_) {

    nr_elements     = 0 ;
    codification    = "undefined" ;
    position        = 0 ;

};

/*! ------------------------------------------------------------------
 * Constructor
 * @param[in]   name_   name of data field
 * @param[in]   comp_   number of coponents of data field [1/3]
 * @param[in]   type_   type of data field [ "[U]Int8", "[U]Int16", "[U]Int32", "[U]Int64", "Float32", "Float64" ]
 * @param[in]   loc_    location of data field ["Cell"/"Point"]
 * @param[in]   cod_    codex ["ascii"/ "appended"]
 * @param[in]   nr_elements_    number of elements of data field
 */
VTK::Field_C::Field_C( std::string name_, uint8_t comp_, std::string type_, std::string loc_, std::string cod_, uint64_t nr_elements_):
    name(name_), components(comp_), type(type_), location(loc_), codification(cod_), nr_elements(nr_elements_){

};

/*! ------------------------------------------------------------------
 * Destructor
 * */
VTK::Field_C::~Field_C(){
};


/*! ------------------------------------------------------------------
 * Set name of data field
 * @param[in]   name_   name of data field
 */
void      VTK::Field_C::SetName( std::string  name_){ 
    name= name_; 
    return; 
};

/*! ------------------------------------------------------------------
 * Set type of data field
 * @param[in]   type_   type of data field [ "[U]Int8", "[U]Int16", "[U]Int32", "[U]Int64", "Float32", "Float64" ]
 */
void      VTK::Field_C::SetType( std::string  type_){

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

        std::cout << type_ << " not admitted in VTK::Field_C::SetType." << std::endl ;
    };


    return; 
};

/*! ------------------------------------------------------------------
 * Set location of data field
 * @param[in]   loc_   location of data field [ "Cell", "Point" ] 
 */
void      VTK::Field_C::SetLocation( std::string  loc_ ){ 
    location= loc_; 
    return; 
};

/*! ------------------------------------------------------------------
 * Set codification of data field
 * @param[in]   code_  codification ["appended"/ "ascii"]
 */
void      VTK::Field_C::SetCodification( std::string  code_ ){ 
    codification= code_; 
    return; 
};

/*! ------------------------------------------------------------------
 * Set number of components of data field
 * @param[in]   comp_   number of components [1/3]
 */
void      VTK::Field_C::SetComponents( uint8_t comp_){ 
    components= comp_; 
    return; 
};

/*! ------------------------------------------------------------------
 * Set numer of elements of data field
 * @param[in]   ele_   number of elements
 */
void      VTK::Field_C::SetElements( uint64_t ele_){ 
    nr_elements= ele_; 
    return; 
};

/*! ------------------------------------------------------------------
 * Set position of data field
 * @param[in]   pos_    position of data field ["Cell"/"Point"]
 */
void      VTK::Field_C::SetPosition( std::fstream::pos_type pos_ ){ 
    position =pos_ ; 
    return; 
} ;

/*! ------------------------------------------------------------------
 * Set offset of data field for appended output
 * @param[in]   offs_   offset from "_" character of appended section
 */
void      VTK::Field_C::SetOffset( uint64_t offs_){ 
    offset= offs_; 
    return; 
};

/*! ------------------------------------------------------------------
 * Get name of data field
 * @return  name of data field
 */
std::string    VTK::Field_C::GetName(){ 
    return name; 
};

/*! ------------------------------------------------------------------
 * Get type of data field
 * @return   type of data field [ "[U]Int8", "[U]Int16", "[U]Int32", "[U]Int64", "Float32", "Float64" ]
 */
std::string    VTK::Field_C::GetType(){ return type; };

/*! ------------------------------------------------------------------
 * Get location of data field
 * @return  location ["Cell"/"Point"]
 */
std::string    VTK::Field_C::GetLocation(){ 
    return location; 
};

/*! ------------------------------------------------------------------
 * Get codification of data field
 * @return  codification ["ascii"/"appended"]
 */
std::string    VTK::Field_C::GetCodification(){ 
    return codification; 
};

/*! ------------------------------------------------------------------
 * Get number of components of data field
 * @return  number of components [1/3]
 */
uint8_t   VTK::Field_C::GetComponents(){ 
    return components; 
};

/*! ------------------------------------------------------------------
 * Get number of elements of data field. 
 * Generally number of points or number of cells.
 * @return  number of elements
 */
uint64_t  VTK::Field_C::GetElements(){ return nr_elements; };

/*! ------------------------------------------------------------------
 * Get total number of entries data field
 * @return  total size = number of elements *nuber of components
 */
uint64_t  VTK::Field_C::GetSize(){ 
    return components *nr_elements ; 
};

/*! ------------------------------------------------------------------
 * Get offset in appended section
 * @return  offset for "_" character in appended section
 */
uint64_t  VTK::Field_C::GetOffset(){ 
    return offset; 
};

/*! ------------------------------------------------------------------
 * Get bytes of  data field
 * @return  memory size of data field
 */
uint64_t  VTK::Field_C::GetNbytes(){ 
    return components *nr_elements *SizeOfType( type ) ;
};

/*! ------------------------------------------------------------------
 * Get position of data field in VTK file. 
 * This information is available after VTK::ReadMetaData() has been called.
 * @return      position in VTK file.
 */
std::fstream::pos_type   VTK::Field_C::GetPosition(){ 
    return position; 
};



/*!
 * @}
 */
