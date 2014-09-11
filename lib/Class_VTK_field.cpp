#include"Class_VTK.hpp"

//------------------------------------------------------------------
VTK::Field_C::Field_C(){};

//------------------------------------------------------------------
VTK::Field_C::Field_C( string name_, int comp_, string type_, string loc_ ):
       name(name_), components(comp_), type(type_), location(loc_) {

       nr_elements  = -9999 ;
       size         = -9999 ;
       nbytes       = -9999 ;
       codification = "undefined" ;


};


//------------------------------------------------------------------
VTK::Field_C::Field_C( string name_, int comp_, string type_, string loc_, string cod_, int nr_elements_):
       name(name_), components(comp_), type(type_), location(loc_), codification(cod_), nr_elements(nr_elements_){

       size   = components *nr_elements ;

       nbytes = size *size_of_type( type ) ;

};

//------------------------------------------------------------------
VTK::Field_C::~Field_C(){};

//------------------------------------------------------------------
int VTK::Field_C::size_of_type( string type ){

  int nbytes ;

  if( type == "Int8" || type == "UInt8"){
    nbytes = sizeof(int8_t) ;
  }

  if( type == "Int16" || type == "UInt16"){
    nbytes = sizeof(int16_t) ;
  }

  else if( type == "Int32" || type == "UInt32"){
    nbytes = sizeof(int32_t) ;
  }

  else if( type == "Int64" || type == "UInt64"){
    nbytes = sizeof(int64_t) ;
  }

  else if( type == "Float32" ){
    nbytes = sizeof(float) ;
  }

  else if( type == "Float64" ){
    nbytes = sizeof(double) ;
  };

  return nbytes ;
};

//------------------------------------------------------------------
void      VTK::Field_C::Set_Name( string  name_){ name= name_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::Set_Type( string  type_){ type= type_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::Set_Location( string  loc_ ){ location= loc_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::Set_Codification( string  code_ ){ codification= code_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::Set_Components( int comp_){ components= comp_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::Set_Elements( int ele_){ 
  nr_elements= ele_; 
  size   = components *nr_elements ;

  nbytes = size *size_of_type( type ) ;
 
  return; 

};

//------------------------------------------------------------------
void      VTK::Field_C::Set_Position( fstream::pos_type pos_ ){ position =pos_ ; return; } ;

//------------------------------------------------------------------
void      VTK::Field_C::Set_Offset( int offs_){ offset= offs_; return; };

//------------------------------------------------------------------
string    VTK::Field_C::Get_Name(){ return name; };

//------------------------------------------------------------------
string    VTK::Field_C::Get_Type(){ return type; };

//------------------------------------------------------------------
string    VTK::Field_C::Get_Location(){ return location; };

//------------------------------------------------------------------
string    VTK::Field_C::Get_Codification(){ return codification; };

//------------------------------------------------------------------
int       VTK::Field_C::Get_Components(){ return components; };

//------------------------------------------------------------------
int       VTK::Field_C::Get_Elements(){ return nr_elements; };

//------------------------------------------------------------------
int       VTK::Field_C::Get_Size(){ return size; };

//------------------------------------------------------------------
int       VTK::Field_C::Get_Offset(){ return offset; };

//------------------------------------------------------------------
int       VTK::Field_C::Get_Nbytes(){ return nbytes; };

//------------------------------------------------------------------
fstream::pos_type   VTK::Field_C::Get_Position(){ return position; };



