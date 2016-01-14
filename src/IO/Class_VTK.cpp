#include "Class_VTK.hpp"

/*! ========================================================================================
 *
 * \ingroup    VisualizationToolKit
 * \{
 *
 * \class      VTK
 * \brief A base class for VTK input output. 
 *
 * VTK provides all basic methods for reading and writing VTK files.
 * ASCII and APPENDED mode are supported.
 *
 */

/*! ------------------------------------------------------------------
 * Default constructor referes to a serial VTK file with appended binary data.
 */
VTK::VTK(){

    nr_data = 0;
    nr_procs = 1;
    my_proc  = 0;

    HeaderType = "UInt32" ;
    
    fh.SetSeries( false ) ;
    fh.SetParallel( false ) ;

    GeomCodex = "appended" ;
    DataCodex = "appended" ;

};

/*! ------------------------------------------------------------------
 * Constructor referes to a serial VTK file with appended binary data.
 * \param[in]  dir_    directory of file with final "/"
 * \param[in]  name_   file name without suffix
 */
VTK::VTK( std::string dir_,  std::string name_ ):
     VTK(){

    SetNames(dir_, name_) ;

};

/*! ------------------------------------------------------------------
 * Destructor
 */
VTK::~VTK(){} ;

/*! ------------------------------------------------------------------
 * Set header type for appended binary output
 * \param[in]  st_    header type ["UInt32"/"UInt64"]
 */
void  VTK::SetHeaderType( std::string st_){ 

    if( st_ == "UInt32" || st_ == "UInt64"){
        HeaderType = st_ ;
    }

    else{
        std::cout << "Unsupported HeaderType " << st_ << std::endl ;
    };

    return; 
} ;

/*! ------------------------------------------------------------------
 * Get header type for appended binary output
 * \return    header type ["UInt32"/"UInt64"]
 */
std::string  VTK::GetHeaderType( ){ 
    return HeaderType ;
};

/*! ------------------------------------------------------------------
 * Set name for VTK file
 * \param[in]  dir_    directory of file with final "/"
 * \param[in]  name_   file name without suffix
 */
void  VTK::SetNames( std::string dir_, std::string name_ ){

  fh.SetDirectory(dir_);
  fh.SetName(name_);

  return ;
};

/*! ------------------------------------------------------------------
 * Activates output for time series. Sets series to true first output index to input
 * \param[in]  c_    first output index
 */
void  VTK::SetCounter( int c_){ 

  fh.SetSeries(true) ;
  fh.SetCounter(c_) ;

  return; 
} ;

/*! ------------------------------------------------------------------
 * Activates parallel output
 * \param[in]  nr    number of processes
 * \param[in]  my    my rank
 */
void  VTK::SetParallel( int nr, int my){ 

  if( nr <  1 ) std::cout << " Numer of processes must be greater than 0" << std::endl ;
  if( my >= nr) std::cout << " my_process is not in valid range " << std::endl ;

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

/*! ------------------------------------------------------------------
 * Sets codex for geometry and data
 * \param[in]  cod_    codex for VTK output ["ascii"/"appended"]
 */
void  VTK::SetCodex( std::string cod_ ) {

    SetGeomCodex( cod_) ;
    SetDataCodex( cod_) ;

    return ;
};

/*! ------------------------------------------------------------------
 * Sets codex for geometry only
 * \param[in]  cod_    codex for VTK output ["ascii"/"appended"]
 */
void  VTK::SetGeomCodex( std::string cod_ ) {

    int i, nf( geometry.size() ) ;

    if( cod_ == "ascii"  || cod_ == "appended" ) {
        GeomCodex = cod_ ;
        for( i=0; i<nf ; i++) geometry[i].SetCodification( cod_ ) ;
    } 
    
    else{
      std::cout << "Codification: " << cod_ << " not supported in VTK::SetGeomCodex"  << std::endl ;
    };

    return;
};

/*! ------------------------------------------------------------------
 * Sets codex for data only
 * \param[in]  cod_    codex for VTK output ["ascii"/"appended"]
 */
void  VTK::SetDataCodex( std::string cod_ ) {

    int i, nf( data.size() ) ;


    if( cod_ == "ascii"  || cod_ == "appended" ) {
        DataCodex = cod_ ;
        for( i=0; i<nf ; i++) data[i].SetCodification( cod_ ) ;
    } 
    
    else{
      std::cout << "Codification: " << cod_ << " not supported in VTK::SetDataCodex"  << std::endl ;
    };

 
    return ;
};

/*! ------------------------------------------------------------------
 * Add user data for input or output. 
 * Codification will be set according to default value [appended] or to value set by VTK::SetDataCodex( std::string ) or VTK::SetCodex( std::string )
 * \param[in]  name_    name of field
 * \param[in]  comp_    nr of components [1/3]
 * \param[in]  type_    type of data [ "[U]Int8", "[U]Int16", "[U]Int32", "[U]Int64", "Float32", "Float6"4 ]
 * \param[in]  loc_     location of data ["Cell"/"Point"]
 */
VTK::Field_C* VTK::AddData( std::string name_, int comp_, std::string type_, std::string loc_ ){

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
           std::cout << "Location: " << loc_ << " not supported in VTK::Add_Data"  << std::endl ;
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
          std::cout << "Type: " << type_ << " not supported in VTK::Add_Data"  << std::endl ;
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

/*! ------------------------------------------------------------------
 * Add user data for input or output. 
 * \param[in]  name_    name of field
 * \param[in]  comp_    nr of components [1/3]
 * \param[in]  type_    type of data [ [U]Int8, [U]Int16, [U]Int32, [U]Int64, Float32, Float64 ]
 * \param[in]  loc_     location of data ["Cell"/"Point"]
 * \param[in]  cod_     codification of data ["ascii"/"appended"]
 */
VTK::Field_C* VTK::AddData( std::string name_, int comp_, std::string type_, std::string loc_, std::string cod_ ){

    int           size_ ;
    bool          allocate(true) ;
    VTK::Field_C* ptr(NULL) ;
    
    if( GetFieldByName( name_, ptr ) ){
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
           std::cout << "Location: " << loc_ << " not supported in VTK::Add_Data"  << std::endl ;
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
          std::cout << "Type: " << type_ << " not supported in VTK::Add_Data"  << std::endl ;
          allocate = false ;
        };
        
        // Check codification ---------------------------------------------------------
        if( cod_ == "ascii"  || cod_ == "appended" ) {
        } 
        
        else{
          std::cout << "Codification: " << cod_ << " not supported in VTK::Add_Data"  << std::endl ;
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

/*! ------------------------------------------------------------------
 * Removes user data from input or output 
 * \param[in]  name_    name of field to be removed
 */
void VTK::RemoveData( std::string name_ ){

    std::vector<VTK::Field_C>::iterator  i_, j_ ;
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
        std::cout << "did not find field for removing: " << name_ << std::endl;
    };
    
    
    return ;

};

/*! ------------------------------------------------------------------
 * Finds data field through name. 
 * Search is performed in user data and geometry
 * \param[in]  name_    name of field to be found
 * \param[out]  the_field    pointer to the filed if found; if not found value is not altered
 * \return      true if field has been found, else false
 */
bool VTK::GetFieldByName( const std::string &name_, VTK::Field_C *&the_field ){


    std::vector<VTK::Field_C>::iterator    it_ ;

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

/*! ------------------------------------------------------------------
 * Calculates the offsets of all geometry and data fields for appended output.
 * Offsets are stored in each field.
 */
void VTK::CalcAppendedOffsets(){

    uint64_t offset(0) ;
    uint64_t  HeaderByte ;
    
    if( GetHeaderType() == "UInt32"){
        HeaderByte = sizeof(uint32_t) ;
    }

    else if( GetHeaderType() == "UInt64") {
        HeaderByte = sizeof(uint64_t) ;
    };
    
    
    for( unsigned i=0; i< nr_data; i++){
        if( data[i].GetCodification() == "appended" && data[i].GetLocation() == "Point") {
            data[i].SetOffset( offset) ;
            offset += HeaderByte + data[i].GetNbytes()  ;
        };
    };
    
    for( unsigned i=0; i< nr_data; i++){
        if( data[i].GetCodification() == "appended" && data[i].GetLocation() == "Cell") {
            data[i].SetOffset( offset) ;
            offset += HeaderByte + data[i].GetNbytes()  ;
        };
    };
    
    for( unsigned i=0; i< geometry.size(); i++){
        if( geometry[i].GetCodification() == "appended" ) {
            geometry[i].SetOffset( offset) ;
            offset += HeaderByte + geometry[i].GetNbytes()  ;
        }; 
    };
    
    
    return ;
};

/*! ------------------------------------------------------------------
 * Converts a string into a Field information.
 * \param[in]   line_       string to be converted
 * \param[out]  data_       Field information
 * \return      true if successful
 */
bool VTK::StringToDataArray( std::string &line_, VTK::Field_C &data_  ){

    std::string type_, name_, code_, comp_, offs_ ;
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

/*! ------------------------------------------------------------------
 * Converts a Field information to string as requested by VTK format.
 * \param[out]   str        string in VTK format
 * \param[in]  field_       Field information
 * \return      true if successful
 */
void  VTK::DataArrayToString( std::string &str, VTK::Field_C &field_ ){

  std::stringstream  os("") ;

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

/*! ------------------------------------------------------------------
 * Converts a parallel field information to string as requested by VTK format.
 * \param[out]   str        string in VTK format
 * \param[in]  field_       Field information
 * \return      true if successful
 */
void  VTK::PDataArrayToString( std::string &str, VTK::Field_C &field_ ){

  std::stringstream  os("") ;

  os << "        <PDataArray "
      << "type=\"" << field_.GetType() << "\" "
      << "Name=\"" << field_.GetName() << "\" "
      << "NumberOfComponents=\""<< unsigned(field_.GetComponents()) << "\" " 
      << ">" ;

  str = os.str() ;

  return ;
  
};

/*! ------------------------------------------------------------------
 * Converts a parallel field information to string as requested by VTK format.
 * \param[in]   type       basic types supported by VTK [ [U]Int8, [U]Int16, [U]Int32, [U]Int64, Float32, Float64 ]
 * \return      size of basic type
 */
uint8_t SizeOfType( std::string type ){

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

/*! ------------------------------------------------------------------
 * Writes entire VTK file (headers and data).
 */
void VTK::Write( ){

    CalcAppendedOffsets() ;

    WriteMetaData() ;
    WriteData() ;

    if( nr_procs != 0  && my_proc == 0)  WriteCollection() ;

    fh.IncrementCounter() ;

    return ;
};

/*! ------------------------------------------------------------------
 * Writes data only in VTK file
 */
void VTK::WriteData( ){

    std::fstream             str ;
    std::fstream::pos_type   position_insert, position_eof ;

    int                 length;
    char*               buffer ;

    str.open( fh.GetName( ), std::ios::in | std::ios::out ) ;

    { // Write Ascii

        position_insert = str.tellg();
        VTK::Field_C    temp ;

        //Writing first point data then cell data
        for( unsigned i=0; i< nr_data; i++){
            if( data[i].GetCodification() == "ascii" && data[i].GetLocation() == "Point") {
            str.seekg( position_insert);
            ReadDataArray( str, data[i] ) ;

            str.seekg( data[i].GetPosition() ) ;
            CopyUntilEOFInString( str, buffer, length );

            Flush( str, "ascii", data[i].GetName() ) ;

            position_insert = str.tellg();
            str << std::endl ;
            flush_binary( str, buffer, length) ;

            delete [] buffer ;

            };
        }; 

        for( unsigned i=0; i< nr_data; i++){
            if( data[i].GetCodification() == "ascii" && data[i].GetLocation() == "Cell") {
            str.seekg( position_insert);
            ReadDataArray( str, data[i] ) ;

            str.seekg( data[i].GetPosition() ) ;
            CopyUntilEOFInString( str, buffer, length );

            Flush( str, "ascii", data[i].GetName() ) ;

            position_insert = str.tellg();
            str << std::endl ;
            flush_binary( str, buffer, length) ;

            delete [] buffer ;
            };
        }; 

        for( unsigned i=0; i< geometry.size(); i++){
            if( geometry[i].GetCodification() == "ascii" ) {
            str.seekg( position_insert);
            ReadDataArray( str, geometry[i] ) ;

            str.seekg( geometry[i].GetPosition() ) ;
            CopyUntilEOFInString( str, buffer, length );

            Flush( str, "ascii", geometry[i].GetName() ) ;

            position_insert = str.tellg();
            str << std::endl ;
            flush_binary( str, buffer, length) ;

            delete [] buffer ;
            };
        }; 

        str.seekg( temp.GetPosition() ) ;


    }
    
    { // Write Appended

        char                c_;
        std::string              line ;
        std::fstream::pos_type   position_appended ;

        //Go to the initial position of the appended section
        while( getline(str, line) && (! Keyword_In_String( line, "<AppendedData")) ){} ;
        
        str >> c_;
        while( c_ != '_') str >> c_;
        
        position_insert = str.tellg();
        CopyUntilEOFInString( str, buffer, length );

        str.close();
        str.clear();
        
        
        //Reopening in binary mode
        str.open( fh.GetName( ), std::ios::out | std::ios::in | std::ios::binary);
        str.seekg( position_insert) ;

        //str.open( "data.dat", std::ios::out | std::ios::binary);

        //Writing first point data then cell data
        for( unsigned i=0; i< nr_data; i++){
            if( data[i].GetCodification() == "appended" && data[i].GetLocation() == "Point") {
                if( GetHeaderType() == "UInt32"){
                    uint32_t    nbytes = data[i].GetNbytes() ;
                    flush_binary(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = data[i].GetNbytes() ;
                    flush_binary(str, nbytes) ;
                };
                Flush( str, "binary", data[i].GetName() ) ;
            };
        } 
        
        for( unsigned i=0; i< nr_data; i++){
            if( data[i].GetCodification() == "appended" && data[i].GetLocation() == "Cell") {
                if( GetHeaderType() == "UInt32"){
                    uint32_t    nbytes = data[i].GetNbytes() ;
                    flush_binary(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = data[i].GetNbytes() ;
                    flush_binary(str, nbytes) ;
                };
                Flush( str, "binary", data[i].GetName() ) ;
            };
        } 
        
        //Writing Geometry Data
        for(unsigned i=0; i<geometry.size(); i++){
            if( geometry[i].GetCodification() == "appended" ) {
                if( GetHeaderType() == "UInt32"){
                    uint32_t    nbytes = geometry[i].GetNbytes() ;
                    flush_binary(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = geometry[i].GetNbytes() ;
                    flush_binary(str, nbytes) ;
                };
                Flush( str, "binary", geometry[i].GetName() ) ;           
            };
        };
   
        // { 
        // std::fstream             str2 ;
        // str2.open( "test2.dat", std::ios::out | std::ios::binary);
        // flush_binary( str2, buffer, length) ;
        // str2.close();
        // }

        flush_binary( str, buffer, length) ;

        delete [] buffer ;
    };
    
    // Closing Appended Secyion
    str.close();
    
    return ;

};

/*! ------------------------------------------------------------------
 * Writes data headers in strean
 * \param[in]   str     output stream
 * \param[in]   parallel     flag for parallel data headers for collection files [true/false]
 */
void VTK::WriteDataHeader( std::fstream &str, bool parallel ){

  std::string        location, line ;
  std::stringstream  scalars, vectors ;

  for( int j=0; j<2; j++){

    if( j==0 ) location = "Point" ;
    if( j==1 ) location = "Cell" ;

//    scalars.clear();
//    vectors.clear();

    scalars.str("");
    vectors.str("");

    //Creating Scalar and Vector Lists
    scalars << "\"" ;
    vectors << "\"" ;

    for( unsigned i=0; i< nr_data; i++ ){

      if( data[i].GetLocation() == location){
        if(      data[i].GetComponents() == 1 ) scalars <<  data[i].GetName() << " " ;
        else if( data[i].GetComponents() == 3 ) vectors <<  data[i].GetName() << " " ;
      };

    };

    scalars << "\"" ;
    vectors << "\"" ;

    if(      location == "Point") {
      str << "      <" ;
      if( parallel )  str << "P" ;
      str << "PointData " ;
    }

    else if( location == "Cell")  {
      str << "      <" ;
      if( parallel )  str << "P" ;
      str << "CellData " ;
    };

    str << " Scalars=" << scalars.str()
        << " Vectors=" << vectors.str()
        << ">" << std::endl;

    //Writing DataArray
    for( unsigned i=0; i< nr_data; i++){
      if( data[i].GetLocation() == location && !parallel) WriteDataArray( str, data[i] ) ;
      if( data[i].GetLocation() == location &&  parallel) WritePDataArray( str, data[i] ); 
    };

    str << "      </" ;
    if( parallel )  str << "P" ;

    if( location == "Point") str << "PointData> " << std::endl;
    if( location == "Cell")  str << "CellData> "  << std::endl;

  };

  return ;

};

/*! ------------------------------------------------------------------
 * Writes data array related to a field in strean
 * \param[in]   str     output stream
 * \param[in]   field_   field to be written
 */
void VTK::WriteDataArray( std::fstream &str, VTK::Field_C &field_ ){

    std::string          line;

    DataArrayToString( line, field_ ) ;

    str << line << std::endl ;
    
    str << "        </DataArray>" << std::endl ;

    return ;

};

/*! ------------------------------------------------------------------
 * Writes parallel data array related to a field in strean
 * \param[in]   str     output stream
 * \param[in]   field_   field to be written
 */
void VTK::WritePDataArray( std::fstream &str, VTK::Field_C &field_ ){

    std::string          line;

    PDataArrayToString( line, field_ ) ;

    str << line << std::endl ;
    str << "        </PDataArray>" << std::endl ;

    return ;

};

/*! ------------------------------------------------------------------
 * Reads entire VTK file (headers and data).
 */
void VTK::Read( ){

    ReadMetaData( );
    ReadData( ) ;

    return ;
};

/*! ------------------------------------------------------------------
 * Reads data only from VTK file
 */
void VTK::ReadData( ){

  std::fstream                  str  ;
  std::fstream::pos_type        position_appended;
  std::string                   line;
  char                     c_ ;
  uint32_t                 nbytes32 ;
  uint64_t                 nbytes64 ;

  str.open( fh.GetName( ), std::ios::in ) ;

  //Read appended data
  //Go to the initial position of the appended section
  while( getline(str, line) && (! Keyword_In_String( line, "<AppendedData")) ){} ;

  str >> c_;
  while( c_ != '_') str >> c_;

  position_appended = str.tellg();


  str.close();
  str.clear();

  //Open in binary for read
  str.open( fh.GetName( ), std::ios::in | std::ios::binary);

  //Read fields
  for( unsigned i=0; i< nr_data; i++){
    if( data[i].GetCodification() == "appended"){
      str.seekg( position_appended) ;
      str.seekg( data[i].GetOffset(), std::ios::cur) ;
      if( HeaderType== "UInt32") absorb_binary( str, nbytes32 ) ;
      if( HeaderType== "UInt64") absorb_binary( str, nbytes64 ) ;
      Absorb( str, "binary", data[i].GetName() ) ;
    };
  };

  //Read geometry
  for(unsigned i=0; i<geometry.size(); i++){
    if( geometry[i].GetCodification() == "appended"){
      str.seekg( position_appended) ;
      str.seekg( geometry[i].GetOffset(), std::ios::cur) ;
      if( HeaderType== "UInt32") absorb_binary( str, nbytes32 ) ;
      if( HeaderType== "UInt64") absorb_binary( str, nbytes64 ) ;
      Absorb( str, "binary", geometry[i].GetName() ) ;
    };
  };

  //Read ascii data
  for( unsigned i=0; i< nr_data; i++){
    if(  data[i].GetCodification() == "ascii"){
      str.seekg( data[i].GetPosition() ) ;
      Absorb( str, "ascii", data[i].GetName() ) ;
    };
  };

  //Read geometry
  for(unsigned i=0; i<geometry.size(); i++){
    if(  geometry[i].GetCodification() == "ascii"){
      str.seekg( geometry[i].GetPosition() ) ;
      Absorb( str, "ascii", geometry[i].GetName() ) ;
    };
  };

  str.close();

  return ; 
};

/*! ------------------------------------------------------------------
 * Reads data headers from strean.
 * All field information available in file are stored.
 * \param[in]   str     output stream
 */
void VTK::ReadDataHeader( std::fstream &str ){


  std::fstream::pos_type   pos_ ;

  std::string                   location ;
  std::string                   line, loc_;
  std::stringstream             ss;

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

/*! ------------------------------------------------------------------
 * Reads data array from stream and stores in field information
 * \param[in]   str     output stream
 * \param[out]   field_   field information
 */
bool  VTK::ReadDataArray( std::fstream &str, VTK::Field_C &field_  ){

  std::string              line_ ;

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

/*! 
 * \} 
 */
