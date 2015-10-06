#include"Class_VTK.hpp"

//------------------------------------------------------------------
VTK::Field_C::Field_C(){

    name            = "undefined" ;  
    type            = "undefined" ;
    location        = "undefined" ;
    codification    = "undefined" ;
    components      = 1 ;   ;
    nr_elements     = 0 ;
    position        = 0 ;

};

//------------------------------------------------------------------
VTK::Field_C::Field_C( string name_, uint8_t comp_, string loc_ ):
       name(name_), components(comp_), location(loc_) {

    type            = "undefined" ;
    nr_elements     = 0 ;
    codification    = "undefined" ;
    position        = 0 ;

};

//------------------------------------------------------------------
VTK::Field_C::Field_C( string name_, uint8_t comp_, string type_, string loc_ ):
       name(name_), components(comp_), type(type_), location(loc_) {

    nr_elements     = 0 ;
    codification    = "undefined" ;
    position        = 0 ;

};

//------------------------------------------------------------------
VTK::Field_C::Field_C( string name_, uint8_t comp_, string type_, string loc_, string cod_, uint64_t nr_elements_):
       name(name_), components(comp_), type(type_), location(loc_), codification(cod_), nr_elements(nr_elements_){

};

//------------------------------------------------------------------
VTK::Field_C::~Field_C(){};


//------------------------------------------------------------------
void      VTK::Field_C::SetName( string  name_){ name= name_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::SetType( string  type_){

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

        cout << type_ << " not admitted in VTK::Field_C::SetType." << endl ;
    };


    return; 
};

//------------------------------------------------------------------
void      VTK::Field_C::SetLocation( string  loc_ ){ location= loc_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::SetCodification( string  code_ ){ codification= code_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::SetComponents( uint8_t comp_){ components= comp_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::SetElements( uint64_t ele_){ nr_elements= ele_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::SetPosition( fstream::pos_type pos_ ){ position =pos_ ; return; } ;

//------------------------------------------------------------------
void      VTK::Field_C::SetOffset( uint64_t offs_){ offset= offs_; return; };

//------------------------------------------------------------------
string    VTK::Field_C::GetName(){ return name; };

//------------------------------------------------------------------
string    VTK::Field_C::GetType(){ return type; };

//------------------------------------------------------------------
string    VTK::Field_C::GetLocation(){ return location; };

//------------------------------------------------------------------
string    VTK::Field_C::GetCodification(){ return codification; };

//------------------------------------------------------------------
uint8_t   VTK::Field_C::GetComponents(){ return components; };

//------------------------------------------------------------------
uint64_t  VTK::Field_C::GetElements(){ return nr_elements; };

//------------------------------------------------------------------
uint64_t  VTK::Field_C::GetSize(){ 
    return components *nr_elements ; 
};

//------------------------------------------------------------------
uint64_t  VTK::Field_C::GetOffset(){ return offset; };

//------------------------------------------------------------------
uint64_t  VTK::Field_C::GetNbytes(){ 
    return components *nr_elements *SizeOfType( type ) ;
};

//------------------------------------------------------------------
fstream::pos_type   VTK::Field_C::GetPosition(){ return position; };



