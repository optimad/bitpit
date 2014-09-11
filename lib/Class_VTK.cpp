
#include "Class_VTK.hpp"

using namespace std;

//------------------------------------------------------------------
VTK::VTK(){

  nr_data = 0;
  nr_procs = 0;
  my_proc  = 0;

  fh.Set_Series( false ) ;
  fh.Set_Parallel( false ) ;

};

//------------------------------------------------------------------
VTK::VTK(string dir_, string name_ ):
     VTK(){

//  nr_data = 0;
//  nr_procs = 0;
//  my_proc  = 0;
//
//  fh.Set_Directory(dir_);
//  fh.Set_Name(name_);
//
//  fh.Set_Series( false ) ;
//  fh.Set_Parallel( false ) ;

  Set_Names(dir_, name_) ;

  return;
};

//------------------------------------------------------------------
VTK::~VTK(){} ;


//-------------------------------------------------------------------
void  VTK::Set_Names( string dir_, string name_ ){

  fh.Set_Directory(dir_);
  fh.Set_Name(name_);

  return ;
};

//------------------------------------------------------------------
void  VTK::Set_Counter( int c_){ 

  fh.Set_Series(true) ;
  fh.Set_Counter(c_) ;

  return; 
} ;


//------------------------------------------------------------------
void  VTK::Set_Parallel( int nr, int my){ 

  if( nr <  1 ) cout << " Numer of processes must be greater than 0" << endl ;
  if( my >= nr) cout << " my_process is not in valid range " << endl ;

  nr_procs = nr; 
  my_proc  = my; 

  if(nr_procs == 0) {
   fh.Set_Parallel(false) ;
  }

  else {
   fh.Set_Parallel(true) ;
   fh.Set_Block( my ) ;
  };

  return; 
} ;


// =================================================================================== //
void VTK::Add_Data( string name_, int comp_, string type_, string loc_, string cod_ ){

  int          size_ ;
  bool         allocate(true) ;

  // Check location -------------------------------------------------------------
  if( loc_ == "Point"){
     size_ = nr_points ;
  }

  else if( loc_ == "Cell" ){
     size_ = nr_cells ;
  }

  else{
     cout << "Location: " << loc_ << " not supported in VTK::Add_Data"  << endl ;
     allocate = false ;
  };

  // Check type -----------------------------------------------------------------
  if(      type_ == "Int8"    || type_ == "UInt8"    ||
           type_ == "Int16"   || type_ == "UInt16"   ||
           type_ == "Int32"   || type_ == "UInt32"   || 
           type_ == "Int64"   || type_ == "UInt64"   || 
           type_ == "Float32" || type_ == "Float64"  ){
  }

  else{
    cout << "Type: " << type_ << " not supported in VTK::Add_Data"  << endl ;
    allocate = false ;
  };

  // Check codification ---------------------------------------------------------
  if( cod_ == "ascii"  || cod_ == "appended" ) {
  } 
  
  else{
    cout << "Codification: " << cod_ << " not supported in VTK::Add_Data"  << endl ;
    allocate = false ;
  };

  // Do allocation if everything ok ---------------------------------------------
  if( allocate) {
    data.push_back( Field_C( name_, comp_, type_, loc_, cod_, size_ ) ) ;
    nr_data++ ;
  };

  return ;

};

// =================================================================================== //
void VTK::Remove_Data( string name_ ){

  vector<VTK::Field_C>::iterator  i_, j_ ;
  bool                            found ;

  found     = false ;

  for ( i_ = data.begin(); i_ != data.end(); i_++){
    if( (i_)->Get_Name() == name_){
       j_ = i_ ;
       found     = true ;
    };
  };

  if( found ){
    data.erase(j_);
  }

  else{
    cout << "did not find field for removing: " << name_ << endl;
  };


  return ;

};

// =================================================================================== //
void  VTK::Flush( fstream &str, string codex_, string name  ){ ;}

// =================================================================================== //
void  VTK::Absorb( fstream &str, string codex_, string name  ){ ;}


// =================================================================================== //
VTK::Field_C* VTK::get_field_by_name( vector<VTK::Field_C> &fields_, const string &name_ ){

  VTK::Field_C*   the_field ;
  bool            found ;

  the_field = NULL ;
  found     = false ;

  for (vector<VTK::Field_C>::iterator i = fields_.begin(); i != fields_.end(); i++){
    if( (i)->Get_Name() == name_){

       the_field = &(*i) ;
       found     = true ;
    };
  };

  if( ! found) cout << " did not find field called: " << name_ << endl;

  return the_field;
};

// =================================================================================== //
void VTK::Calc_Appended_Offsets(){

  int offset(0) ;

  for( int i=0; i< nr_data; i++){
    if( data[i].Get_Codification() == "appended" && data[i].Get_Location() == "Point") {
      data[i].Set_Offset( offset) ;
      offset += sizeof(int) + data[i].Get_Nbytes()  ;
    };
  };

  for( int i=0; i< nr_data; i++){
    if( data[i].Get_Codification() == "appended" && data[i].Get_Location() == "Cell") {
      data[i].Set_Offset( offset) ;
      offset += sizeof(int) + data[i].Get_Nbytes()  ;
    };
  };

  for( int i=0; i< geometry.size(); i++){
    if( geometry[i].Get_Codification() == "appended" ) {
      geometry[i].Set_Offset( offset) ;
      offset += sizeof(int) + geometry[i].Get_Nbytes()  ;
    }; 
  };


  return ;
};


