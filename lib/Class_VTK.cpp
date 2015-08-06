
#include "Class_VTK.hpp"

using namespace std;

//------------------------------------------------------------------
VTK::VTK(){

  nr_data = 0;
  nr_procs = 0;
  my_proc  = 0;

  fh.SetSeries( false ) ;
  fh.SetParallel( false ) ;

};

//------------------------------------------------------------------
VTK::VTK(string dir_, string name_ ):
     VTK(){

  SetNames(dir_, name_) ;

  return;
};

//------------------------------------------------------------------
VTK::~VTK(){} ;


//-------------------------------------------------------------------
void  VTK::SetNames( string dir_, string name_ ){

  fh.SetDirectory(dir_);
  fh.SetName(name_);

  return ;
};

//------------------------------------------------------------------
void  VTK::SetCounter( int c_){ 

  fh.SetSeries(true) ;
  fh.SetCounter(c_) ;

  return; 
} ;

//------------------------------------------------------------------
void  VTK::SetParallel( int nr, int my){ 

  if( nr <  1 ) cout << " Numer of processes must be greater than 0" << endl ;
  if( my >= nr) cout << " my_process is not in valid range " << endl ;

  nr_procs = nr; 
  my_proc  = my; 

  if(nr_procs == 0) {
   fh.SetParallel(false) ;
  }

  else {
   fh.SetParallel(true) ;
   fh.SetBlock( my ) ;
  };

  return; 
} ;

//------------------------------------------------------------------
void  VTK::SetGeomCodex( string cod_ ) {

    int i, nf( geometry.size() ) ;

    for( i=0; i<nf ; i++) geometry[i].SetCodification( cod_ ) ;
};

//------------------------------------------------------------------
void  VTK::SetDataCodex( string cod_ ) {

    int i, nf( data.size() ) ;

    for( i=0; i<nf ; i++) data[i].SetCodification( cod_ ) ;
};

// =================================================================================== //
void VTK::AddData( string name_, int comp_, string type_, string loc_, string cod_ ){

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
void VTK::RemoveData( string name_ ){

  vector<VTK::Field_C>::iterator  i_, j_ ;
  bool                            found ;

  found     = false ;

  for ( i_ = data.begin(); i_ != data.end(); i_++){
    if( (i_)->GetName() == name_){
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
bool VTK::GetFieldByName( vector<VTK::Field_C> &fields_, const string &name_, VTK::Field_C *&the_field ){


  vector<VTK::Field_C>::iterator    it_ ;

  for( it_=fields_.begin(); it_!=fields_.end(); ++it_){
      if( (*it_).GetName() == name_ ){
          the_field = &(*it_) ;
          return true ;
      };
  };

  return false ;
};

// =================================================================================== //
void VTK::CalcAppendedOffsets(){

  int offset(0) ;

  for( int i=0; i< nr_data; i++){
    if( data[i].GetCodification() == "appended" && data[i].GetLocation() == "Point") {
      data[i].SetOffset( offset) ;
      offset += sizeof(int) + data[i].GetNbytes()  ;
    };
  };

  for( int i=0; i< nr_data; i++){
    if( data[i].GetCodification() == "appended" && data[i].GetLocation() == "Cell") {
      data[i].SetOffset( offset) ;
      offset += sizeof(int) + data[i].GetNbytes()  ;
    };
  };

  for( int i=0; i< geometry.size(); i++){
    if( geometry[i].GetCodification() == "appended" ) {
      geometry[i].SetOffset( offset) ;
      offset += sizeof(int) + geometry[i].GetNbytes()  ;
    }; 
  };


  return ;
};

// =================================================================================== //
bool VTK::StringToDataArray( string &line_, VTK::Field_C &data_  ){

  string type_, name_, code_ ;
  int    comp_, offs_ ;


  if( Keyword_In_String( line_, "<DataArray ") ){  
    type_ = Get_After_Keyword( line_, "type=", '\"') ;
    name_ = Get_After_Keyword( line_, "Name=", '\"') ;
    code_ = Get_After_Keyword( line_, "format=", '\"') ;
    convert_string( Get_After_Keyword( line_, "NumberOfComponents=", '\"'), comp_ ) ;
  
    data_.SetType(type_) ;
    data_.SetName(name_) ;
    data_.SetComponents(comp_) ;
    data_.SetCodification(code_) ;
  
    if(code_=="appended") {
      convert_string( Get_After_Keyword( line_, "offset=", '\"'), offs_ ) ;
      data_.SetOffset(offs_) ;
    }

    return true ;
  }

  else{
    return false ;
  };

 
};

// =================================================================================== //
void  VTK::DataArrayToString( string &str, VTK::Field_C &field_ ){

  stringstream  os("") ;

  os << "        <DataArray "
       << "type=\"" << field_.GetType() << "\" "
       << "Name=\"" << field_.GetName() << "\" "
       << "NumberOfComponents=\""<< field_.GetComponents() << "\" "
       << "format=\"" << field_.GetCodification() << "\" ";

  if( field_.GetCodification() == "appended"){
    os << "offset=\"" << field_.GetOffset() << "\" " ;
  };
  
  os << ">" ;

  str = os.str() ;       
       
  return ;
  
};

// =================================================================================== //
void  VTK::PDataArrayToString( string &str, VTK::Field_C &field_ ){

  stringstream  os("") ;

  os << "        <PDataArray "
      << "type=\"" << field_.GetType() << "\" "
      << "Name=\"" << field_.GetName() << "\" "
      << "NumberOfComponents=\""<< field_.GetComponents() << "\" " 
      << "/>" ;

  str = os.str() ;

  return ;
  
};


