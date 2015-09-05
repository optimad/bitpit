
#include "Class_VTK.hpp"

using namespace std;

//------------------------------------------------------------------
VTK::VTK(){

    nr_data = 0;
    nr_procs = 0;
    my_proc  = 0;

    data.reserve(50) ;

    HeaderType = "UInt32" ;
    
    fh.SetSeries( false ) ;
    fh.SetParallel( false ) ;

    GeomCodex = "undefined" ;
    DataCodex = "undefined" ;

};

//------------------------------------------------------------------
VTK::VTK(string dir_, string name_ ):
     VTK(){

    SetNames(dir_, name_) ;

};

//------------------------------------------------------------------
VTK::~VTK(){} ;

//------------------------------------------------------------------
void  VTK::SetHeaderType( string st_){ 

    if( st_ == "UInt32" || st_ == "UInt64"){
        HeaderType = st_ ;
    }

    else{
        cout << "Unsupported HeaderType " << st_ << endl ;
    };

    return; 
} ;

//------------------------------------------------------------------
string  VTK::GetHeaderType( ){ 
    return HeaderType ;
};

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
void  VTK::SetCodex( string cod_ ) {

    SetGeomCodex( cod_) ;
    SetDataCodex( cod_) ;

    return ;
};

//------------------------------------------------------------------
void  VTK::SetGeomCodex( string cod_ ) {

    int i, nf( geometry.size() ) ;

    if( cod_ == "ascii"  || cod_ == "appended" ) {
        GeomCodex = cod_ ;
        for( i=0; i<nf ; i++) geometry[i].SetCodification( cod_ ) ;
    } 
    
    else{
      cout << "Codification: " << cod_ << " not supported in VTK::SetGeomCodex"  << endl ;
    };

    return;
};

//------------------------------------------------------------------
void  VTK::SetDataCodex( string cod_ ) {

    int i, nf( data.size() ) ;


    if( cod_ == "ascii"  || cod_ == "appended" ) {
        DataCodex = cod_ ;
        for( i=0; i<nf ; i++) data[i].SetCodification( cod_ ) ;
    } 
    
    else{
      cout << "Codification: " << cod_ << " not supported in VTK::SetDataCodex"  << endl ;
    };

 
    return ;
};

// =================================================================================== //
VTK::Field_C* VTK::AddData( string name_, int comp_, string type_, string loc_ ){

    int          size_ ;
    bool         allocate(true) ;
    VTK::Field_C* ptr(NULL) ;
    
    if(  GetFieldByName( name_, ptr ) ) {
        ptr->SetComponents(comp_) ;
        ptr->SetType(type_) ;
        ptr->SetLocation(loc_) ;
    
    }
    
    else{
    
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
        
        // Do allocation if everything ok ---------------------------------------------
        if( allocate) {
          data.push_back( Field_C( name_, comp_, type_, loc_, DataCodex, size_ ) ) ;
          ptr = &data[nr_data] ;
          nr_data++ ;
        };

    };
    
    return ptr ;

};

// =================================================================================== //
VTK::Field_C* VTK::AddData( string name_, int comp_, string type_, string loc_, string cod_ ){

    int           size_ ;
    bool          allocate(true) ;
    VTK::Field_C* ptr(NULL) ;
    
    if( ! GetFieldByName( name_, ptr ) ){
        ptr->SetComponents(comp_) ;
        ptr->SetType(type_) ;
        ptr->SetLocation(loc_) ;
        ptr->SetCodification(cod_) ;
    }

    else{
    
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
          ptr = &data[nr_data] ;
          nr_data++ ;
        };

    };
    
    return ptr ;

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
        nr_data--;
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
bool VTK::GetFieldByName( const string &name_, VTK::Field_C *&the_field ){


    vector<VTK::Field_C>::iterator    it_ ;

    for( it_=data.begin(); it_!=data.end(); ++it_){
        if( (*it_).GetName() == name_ ){
            the_field = &(*it_) ;
            return true ;
        };
    };

    for( it_=geometry.begin(); it_!=geometry.end(); ++it_){
        if( (*it_).GetName() == name_ ){
            the_field = &(*it_) ;
            return true ;
        };
    };

    return false ;
};

// =================================================================================== //
void VTK::CalcAppendedOffsets(){

    uint64_t offset(0) ;
    uint64_t  HeaderByte ;
    
    if( GetHeaderType() == "UInt32"){
        HeaderByte = sizeof(uint32_t) ;
    }

    else if( GetHeaderType() == "UInt64") {
        HeaderByte = sizeof(uint64_t) ;
    };
    
    
    for( int i=0; i< nr_data; i++){
        if( data[i].GetCodification() == "appended" && data[i].GetLocation() == "Point") {
            data[i].SetOffset( offset) ;
            offset += HeaderByte + data[i].GetNbytes()  ;
        };
    };
    
    for( int i=0; i< nr_data; i++){
        if( data[i].GetCodification() == "appended" && data[i].GetLocation() == "Cell") {
            data[i].SetOffset( offset) ;
            offset += HeaderByte + data[i].GetNbytes()  ;
        };
    };
    
    for( int i=0; i< geometry.size(); i++){
        if( geometry[i].GetCodification() == "appended" ) {
            geometry[i].SetOffset( offset) ;
            offset += HeaderByte + geometry[i].GetNbytes()  ;
        }; 
    };
    
    
    return ;
};

// =================================================================================== //
bool VTK::StringToDataArray( string &line_, VTK::Field_C &data_  ){

    string type_, name_, code_, comp_, offs_ ;
    int    components_(1), offset_ ;
    
    bool  success(true) ;
    
    
    if( Keyword_In_String( line_, "<DataArray ") ){  
        success = success && Get_After_Keyword( line_, "type=", '\"', type_) ;
        success = success && Get_After_Keyword( line_, "Name=", '\"', name_) ;
        success = success && Get_After_Keyword( line_, "format=", '\"', code_) ;
    
        if( Get_After_Keyword( line_, "NumberOfComponents=", '\"', comp_)  ){
            convert_string( comp_, components_ ) ;
        };
    
        data_.SetType(type_) ;
        data_.SetName(name_) ;
        data_.SetComponents(components_) ;
        data_.SetCodification(code_) ;
    
        if(code_=="appended") {
            if( Get_After_Keyword( line_, "offset=", '\"', offs_) ){
                convert_string( offs_, offset_ ) ;
                data_.SetOffset(offset_) ;
            }
            else{
                success = false ;
            };
        }
    
        return success ;
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
       << "NumberOfComponents=\""<< unsigned(field_.GetComponents()) << "\" "
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

//------------------------------------------------------------------
uint8_t SizeOfType( string type ){

  uint8_t nbytes ;

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
