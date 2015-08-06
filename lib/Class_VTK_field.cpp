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

       nbytes = size *SizeOfType( type ) ;

};

//------------------------------------------------------------------
VTK::Field_C::~Field_C(){};

//------------------------------------------------------------------
int VTK::Field_C::SizeOfType( string type ){

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
void      VTK::Field_C::SetName( string  name_){ name= name_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::SetType( string  type_){ type= type_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::SetLocation( string  loc_ ){ location= loc_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::SetCodification( string  code_ ){ codification= code_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::SetComponents( int comp_){ components= comp_; return; };

//------------------------------------------------------------------
void      VTK::Field_C::SetElements( int ele_){ 
  nr_elements= ele_; 
  size   = components *nr_elements ;

  nbytes = size *SizeOfType( type ) ;
 
  return; 

};

//------------------------------------------------------------------
void      VTK::Field_C::SetPosition( fstream::pos_type pos_ ){ position =pos_ ; return; } ;

//------------------------------------------------------------------
void      VTK::Field_C::SetOffset( int offs_){ offset= offs_; return; };

//------------------------------------------------------------------
string    VTK::Field_C::GetName(){ return name; };

//------------------------------------------------------------------
string    VTK::Field_C::GetType(){ return type; };

//------------------------------------------------------------------
string    VTK::Field_C::GetLocation(){ return location; };

//------------------------------------------------------------------
string    VTK::Field_C::GetCodification(){ return codification; };

//------------------------------------------------------------------
int       VTK::Field_C::GetComponents(){ return components; };

//------------------------------------------------------------------
int       VTK::Field_C::GetElements(){ return nr_elements; };

//------------------------------------------------------------------
int       VTK::Field_C::GetSize(){ return size; };

//------------------------------------------------------------------
int       VTK::Field_C::GetOffset(){ return offset; };

//------------------------------------------------------------------
int       VTK::Field_C::GetNbytes(){ return nbytes; };

//------------------------------------------------------------------
fstream::pos_type   VTK::Field_C::GetPosition(){ return position; };



