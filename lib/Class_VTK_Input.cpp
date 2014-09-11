

#include "Class_VTK.hpp"

using namespace std;


// ==================================================================================== //
void VTK::Read_Data_Header( fstream &str ){


  string                   location ;
  string                   line, loc_;
  stringstream             ss;

  bool                     read ;

  Field_C                  temp ;

  data.clear() ;
  nr_data = 0  ;

  for( int i=0; i<2; i++){

    ss.str("") ;
    if( i== 0) location = "Point" ;
    if( i== 1) location = "Cell" ;

    temp.Set_Location( location ) ;

    ss << "</" << location << "Data>" ;
    loc_ = ss.str();
 
    read= true ; 
    if( ! getline( str, line) ) read = false ;
    if( Keyword_In_String( line, loc_) ) read=false ;
    

//    cout << line << endl ;
//    cout << loc_ << endl ;
//    cout << read << endl ;

    while( read ){
      if( Read_DataArray( line, temp  ) ) {
        data.push_back( temp ) ;
        nr_data++ ;
      };

      if( ! getline( str, line) ) read = false ;
      if( Keyword_In_String( line, loc_) ) read=false ;
    }; 

  };

  return ;
};

// =================================================================================== //
void VTK::Read_FieldValues( fstream &str ){

  fstream::pos_type        position_appended;
  string                   line;
  char                     c_ ;

  //Read appended data
  //Go to the initial position of the appended section
  while( getline(str, line) && (! Keyword_In_String( line, "<AppendedData")) ){} ;

  str >> c_;
  while( c_ != '_') str >> c_;

  position_appended = str.tellg();


  str.close();
  str.clear();

  //Open in binary for read
  str.open( fh.Get_Name( ), ios::in | ios::binary);

  //Read fields
  for( int i=0; i< nr_data; i++){
    if( data[i].Get_Codification() == "appended"){
      str.seekg( position_appended) ;
      str.seekg( data[i].Get_Offset(), ios::cur) ;
      Absorb( str, "binary", data[i].Get_Name() ) ;
    };
  };

  //Read geometry
  for(int i=0; i<geometry.size(); i++){
    if( geometry[i].Get_Codification() == "appended"){
      str.seekg( position_appended) ;
      str.seekg( geometry[i].Get_Offset(), ios::cur) ;
      Absorb( str, "binary", geometry[i].Get_Name() ) ;
    };
  };

  //Read ascii data
  for( int i=0; i< nr_data; i++){
    if(  data[i].Get_Codification() == "ascii"){
      str.seekg( data[i].Get_Position() ) ;
      Absorb( str, "ascii", data[i].Get_Name() ) ;
    };
  };

  //Read geometry
  for(int i=0; i<geometry.size(); i++){
    if(  geometry[i].Get_Codification() == "ascii"){
      str.seekg( geometry[i].Get_Position() ) ;
      Absorb( str, "ascii", geometry[i].Get_Name() ) ;
    };
  };



  return ; 
};

// =================================================================================== //
bool VTK::Read_DataArray( string &line_, VTK::Field_C &data_  ){

  string type_, name_, code_ ;
  int    comp_, offs_ ;


  if( Keyword_In_String( line_, "<DataArray ") ){  
    type_ = Get_After_Keyword( line_, "type=", '\"') ;
    name_ = Get_After_Keyword( line_, "Name=", '\"') ;
    code_ = Get_After_Keyword( line_, "format=", '\"') ;
    convert_string( Get_After_Keyword( line_, "NumberOfComponents=", '\"'), comp_ ) ;
  
    data_.Set_Type(type_) ;
    data_.Set_Name(name_) ;
    data_.Set_Components(comp_) ;
    data_.Set_Codification(code_) ;
  
    if(code_=="appended") {
      convert_string( Get_After_Keyword( line_, "offset=", '\"'), offs_ ) ;
      data_.Set_Offset(offs_) ;
    };

    return true ;
  }

  else{
    return false ;
  };

 
};


// =================================================================================== //
bool  VTK::Seek_and_Read( fstream &str, const string &name_, VTK::Field_C &field_  ){

  bool                read_, found_ ;
  string              line_ ;
  fstream::pos_type   pos_ ;

  found_ = false ;
  read_  = true ;

  if( !getline(str, line_) ) read_ =false ;

  while( read_  ){


    if( Keyword_In_String( line_, name_) ){
      Read_DataArray( line_, field_  ) ;
      if( field_.Get_Codification() == "ascii") {
        pos_ = str.tellg() ;
        field_.Set_Position( pos_ ) ;
      };
      found_ = true ;
      read_  = false ;
    };

    if(read_){
      if( !getline(str, line_) ) read_ =false ;
    };

  };
 
  return found_ ; 

};
