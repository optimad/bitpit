

#include "Class_VTK.hpp"

using namespace std;

// ==================================================================================== //
void VTK::Read( ){

    ReadMetaData( );
    ReadData( ) ;

    return ;
};

// ==================================================================================== //
void VTK::ReadDataHeader( fstream &str ){


  fstream::pos_type   pos_ ;

  string                   location ;
  string                   line, loc_;
  stringstream             ss;

  int                      size ;
  bool                     read ;

  Field_C                  temp ;
  Field_C*                 ptemp ;


  for( int i=0; i<2; i++){

    ss.str("") ;
    if( i== 0) location = "Point" ;
    if( i== 1) location = "Cell" ;

    temp.SetLocation( location ) ;

    ss << "</" << location << "Data>" ;
    loc_ = ss.str();
 
    read= true ; 
    if( ! getline( str, line) ) read = false ;
    if( Keyword_In_String( line, loc_) ) read=false ;
    

    while( read ){
      if( StringToDataArray( line, temp  ) ) {

        if( temp.GetCodification() == "ascii") {
            pos_ = str.tellg() ;
        }

        else{
            pos_ =  0 ; 
        };

        temp.SetPosition( pos_ ) ;

        if( ! GetFieldByName( temp.GetName(), ptemp )) {
            data.push_back( temp ) ;
            nr_data++ ;
        }
        
        else{
            *ptemp = temp ;
        };

      };

      if( ! getline( str, line) ) read = false ;
      if( Keyword_In_String( line, loc_) ) read=false ;
    }; 

  };

  return ;
};

// =================================================================================== //
void VTK::ReadData( ){

  fstream                  str  ;
  fstream::pos_type        position_appended;
  string                   line;
  char                     c_ ;
  uint32_t                 nbytes32 ;
  uint64_t                 nbytes64 ;

  str.open( fh.GetName( ), ios::in ) ;

  //Read appended data
  //Go to the initial position of the appended section
  while( getline(str, line) && (! Keyword_In_String( line, "<AppendedData")) ){} ;

  str >> c_;
  while( c_ != '_') str >> c_;

  position_appended = str.tellg();


  str.close();
  str.clear();

  //Open in binary for read
  str.open( fh.GetName( ), ios::in | ios::binary);

  //Read fields
  for( int i=0; i< nr_data; i++){
    if( data[i].GetCodification() == "appended"){
      str.seekg( position_appended) ;
      str.seekg( data[i].GetOffset(), ios::cur) ;
      if( HeaderType== "UInt32") absorb_binary( str, nbytes32 ) ;
      if( HeaderType== "UInt64") absorb_binary( str, nbytes64 ) ;
      Absorb( str, "binary", data[i].GetName() ) ;
    };
  };

  //Read geometry
  for(int i=0; i<geometry.size(); i++){
    if( geometry[i].GetCodification() == "appended"){
      str.seekg( position_appended) ;
      str.seekg( geometry[i].GetOffset(), ios::cur) ;
      if( HeaderType== "UInt32") absorb_binary( str, nbytes32 ) ;
      if( HeaderType== "UInt64") absorb_binary( str, nbytes64 ) ;
      Absorb( str, "binary", geometry[i].GetName() ) ;
    };
  };

  //Read ascii data
  for( int i=0; i< nr_data; i++){
    if(  data[i].GetCodification() == "ascii"){
      str.seekg( data[i].GetPosition() ) ;
      Absorb( str, "ascii", data[i].GetName() ) ;
    };
  };

  //Read geometry
  for(int i=0; i<geometry.size(); i++){
    if(  geometry[i].GetCodification() == "ascii"){
      str.seekg( geometry[i].GetPosition() ) ;
      Absorb( str, "ascii", geometry[i].GetName() ) ;
    };
  };

  str.close();

  return ; 
};

// =================================================================================== //
bool  VTK::ReadDataArray( fstream &str, VTK::Field_C &field_  ){

  string              line_ ;

  while( getline(str, line_)  ){

    if( Keyword_In_String( line_, field_.GetName() ) ){
      if( StringToDataArray( line_, field_  ) ){

        if( field_.GetCodification() == "ascii") {
          field_.SetPosition( str.tellg() ) ;
        };

        return true ;
      };

    };

  };
 
  return false ; 

};
