#include "VTK.hpp"

/*!
 * @ingroup    VisualizationToolKit
 * @{
 *
 * @class      VTK
 * @brief A base class for VTK input output. 
 *
 * VTK provides all basic methods for reading and writing VTK files.
 * ASCII and APPENDED mode are supported.
 *
 */

/*! 
 * Default constructor referes to a serial VTK file with appended binary data.
 */
VTK::VTK(){

    nr_data = 0;
    nr_procs = 1;
    my_proc  = 0;

    HeaderType = "UInt32" ;
    
    fh.setSeries( false ) ;
    fh.setParallel( false ) ;

    GeomCodex = "appended" ;
    DataCodex = "appended" ;

};

/*! 
 * Constructor referes to a serial VTK file with appended binary data.
 * @param[in]  dir_    directory of file with final "/"
 * @param[in]  name_   file name without suffix
 */
VTK::VTK( std::string dir_,  std::string name_ ):
     VTK(){

    setNames(dir_, name_) ;

};

/*!
 * Destructor
 */
VTK::~VTK(){} ;

/*! 
 * set header type for appended binary output
 * @param[in]  st_    header type ["UInt32"/"UInt64"]
 */
void  VTK::setHeaderType( std::string st_){ 

    if( st_ == "UInt32" || st_ == "UInt64"){
        HeaderType = st_ ;
    }

    else{
        std::cout << "Unsupported HeaderType " << st_ << std::endl ;
    };

    return; 
} ;

/*! 
 * Get header type for appended binary output
 * @return    header type ["UInt32"/"UInt64"]
 */
std::string  VTK::getHeaderType( ){ 
    return HeaderType ;
};

/*! 
 * set name for VTK file
 * @param[in]  dir_    directory of file with final "/"
 * @param[in]  name_   file name without suffix
 */
void  VTK::setNames( std::string dir_, std::string name_ ){

  fh.setDirectory(dir_);
  fh.setName(name_);

  return ;
};

/*!
 * Activates output for time series. sets series to true first output index to input
 * @param[in]  c_    first output index
 */
void  VTK::setCounter( int c_){ 

  fh.setSeries(true) ;
  fh.setCounter(c_) ;

  return; 
} ;

/*!
 * Activates parallel output
 * @param[in]  nr    number of processes
 * @param[in]  my    my rank
 */
void  VTK::setParallel( int nr, int my){ 

  if( nr <  1 ) std::cout << " Numer of processes must be greater than 0" << std::endl ;
  if( my >= nr) std::cout << " my_process is not in valid range " << std::endl ;

  nr_procs = nr; 
  my_proc  = my; 

  if(nr_procs == 0) {
   fh.setParallel(false) ;
  }

  else {
   fh.setParallel(true) ;
   fh.setBlock( my ) ;
  };

  return; 
} ;

/*!
 * sets codex for geometry and data
 * @param[in]  cod_    codex for VTK output ["ascii"/"appended"]
 */
void  VTK::setCodex( std::string cod_ ) {

    setGeomCodex( cod_) ;
    setDataCodex( cod_) ;

    return ;
};

/*!
 * sets codex for geometry only
 * @param[in]  cod_    codex for VTK output ["ascii"/"appended"]
 */
void  VTK::setGeomCodex( std::string cod_ ) {

    int i, nf( geometry.size() ) ;

    if( cod_ == "ascii"  || cod_ == "appended" ) {
        GeomCodex = cod_ ;
        for( i=0; i<nf ; i++) geometry[i].setCodification( cod_ ) ;
    } 
    
    else{
      std::cout << "Codification: " << cod_ << " not supported in VTK::setGeomCodex"  << std::endl ;
    };

    return;
};

/*!
 * sets codex for data only
 * @param[in]  cod_    codex for VTK output ["ascii"/"appended"]
 */
void  VTK::setDataCodex( std::string cod_ ) {

    int i, nf( data.size() ) ;


    if( cod_ == "ascii"  || cod_ == "appended" ) {
        DataCodex = cod_ ;
        for( i=0; i<nf ; i++) data[i].setCodification( cod_ ) ;
    } 
    
    else{
      std::cout << "Codification: " << cod_ << " not supported in VTK::setDataCodex"  << std::endl ;
    };

 
    return ;
};

/*!
 * Add user data for input or output. 
 * Codification will be set according to default value [appended] or to value set by VTK::setDataCodex( std::string ) or VTK::setCodex( std::string )
 * @param[in]  name_    name of field
 * @param[in]  comp_    nr of components [1/3]
 * @param[in]  type_    type of data [ "[U]Int8", "[U]Int16", "[U]Int32", "[U]Int64", "Float32", "Float6"4 ]
 * @param[in]  loc_     location of data ["Cell"/"Point"]
 */
VTKField* VTK::addData( std::string name_, int comp_, std::string type_, std::string loc_ ){

    int         size_ ;
    bool        allocate(true) ;
    VTKField*   ptr(NULL) ;
    
    if(  getFieldByName( name_, ptr ) ) {
        ptr->setComponents(comp_) ;
        ptr->setType(type_) ;
        ptr->setLocation(loc_) ;
    
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
           std::cout << "Location: " << loc_ << " not supported in VTK::addData"  << std::endl ;
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
          std::cout << "Type: " << type_ << " not supported in VTK::add_Data"  << std::endl ;
          allocate = false ;
        };
        
        // Do allocation if everything ok ---------------------------------------------
        if( allocate) {
          data.push_back( VTKField( name_, comp_, type_, loc_, DataCodex, size_ ) ) ;
          ptr = &data[nr_data] ;
          nr_data++ ;
        };

    };
    
    return ptr ;

};

/*!
 * add user data for input or output. 
 * @param[in]  name_    name of field
 * @param[in]  comp_    nr of components [1/3]
 * @param[in]  type_    type of data [ [U]Int8, [U]Int16, [U]Int32, [U]Int64, Float32, Float64 ]
 * @param[in]  loc_     location of data ["Cell"/"Point"]
 * @param[in]  cod_     codification of data ["ascii"/"appended"]
 */
VTKField* VTK::addData( std::string name_, int comp_, std::string type_, std::string loc_, std::string cod_ ){

    int         size_ ;
    bool        allocate(true) ;
    VTKField*   ptr(NULL) ;
    
    if( getFieldByName( name_, ptr ) ){
        ptr->setComponents(comp_) ;
        ptr->setType(type_) ;
        ptr->setLocation(loc_) ;
        ptr->setCodification(cod_) ;
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
           std::cout << "Location: " << loc_ << " not supported in VTK::add_Data"  << std::endl ;
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
          std::cout << "Type: " << type_ << " not supported in VTK::add_Data"  << std::endl ;
          allocate = false ;
        };
        
        // Check codification ---------------------------------------------------------
        if( cod_ == "ascii"  || cod_ == "appended" ) {
        } 
        
        else{
          std::cout << "Codification: " << cod_ << " not supported in VTK::add_Data"  << std::endl ;
          allocate = false ;
        };
        
        // Do allocation if everything ok ---------------------------------------------
        if( allocate) {
          data.push_back( VTKField( name_, comp_, type_, loc_, cod_, size_ ) ) ;
          ptr = &data[nr_data] ;
          nr_data++ ;
        };

    };
    
    return ptr ;

};

/*!
 * Removes user data from input or output 
 * @param[in]  name_    name of field to be removed
 */
void VTK::removeData( std::string name_ ){

    std::vector<VTKField>::iterator  i_, j_ ;
    bool                            found ;
    
    found     = false ;
    
    for ( i_ = data.begin(); i_ != data.end(); i_++){
        if( (i_)->getName() == name_){
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

/*!
 * Finds data field through name. 
 * Search is performed in user data and geometry
 * @param[in]  name_    name of field to be found
 * @param[out]  the_field    pointer to the filed if found; if not found value is not altered
 * \return      true if field has been found, else false
 */
bool VTK::getFieldByName( const std::string &name_, VTKField *&the_field ){


    std::vector<VTKField>::iterator    it_ ;

    for( it_=data.begin(); it_!=data.end(); ++it_){
        if( (*it_).getName() == name_ ){
            the_field = &(*it_) ;
            return true ;
        };
    };

    for( it_=geometry.begin(); it_!=geometry.end(); ++it_){
        if( (*it_).getName() == name_ ){
            the_field = &(*it_) ;
            return true ;
        };
    };

    return false ;
};

/*!
 * Calculates the offsets of all geometry and data fields for appended output.
 * Offsets are stored in each field.
 */
void VTK::calcAppendedOffsets(){

    uint64_t offset(0) ;
    uint64_t  HeaderByte ;
    
    if( getHeaderType() == "UInt32"){
        HeaderByte = sizeof(uint32_t) ;
    }

    else if( getHeaderType() == "UInt64") {
        HeaderByte = sizeof(uint64_t) ;
    };
    
    
    for( unsigned i=0; i< nr_data; i++){
        if( data[i].getCodification() == "appended" && data[i].getLocation() == "Point") {
            data[i].setOffset( offset) ;
            offset += HeaderByte + data[i].getNbytes()  ;
        };
    };
    
    for( unsigned i=0; i< nr_data; i++){
        if( data[i].getCodification() == "appended" && data[i].getLocation() == "Cell") {
            data[i].setOffset( offset) ;
            offset += HeaderByte + data[i].getNbytes()  ;
        };
    };
    
    for( unsigned i=0; i< geometry.size(); i++){
        if( geometry[i].getCodification() == "appended" ) {
            geometry[i].setOffset( offset) ;
            offset += HeaderByte + geometry[i].getNbytes()  ;
        }; 
    };
    
    
    return ;
};

/*!
 * Writes entire VTK file (headers and data).
 */
void VTK::write( ){

    calcAppendedOffsets() ;

    writeMetaData() ;
    writeData() ;

    if( nr_procs > 1  && my_proc == 0)  writeCollection() ;

    fh.incrementCounter() ;

    return ;
};

/*!
 * Writes data only in VTK file
 */
void VTK::writeData( ){

    std::fstream             str ;
    std::fstream::pos_type   position_insert, position_eof ;

    int                 length;
    char*               buffer ;

    str.open( fh.getName( ), std::ios::in | std::ios::out ) ;

    { // Write Ascii

        position_insert = str.tellg();
        VTKField    temp ;

        //Writing first point data then cell data
        for( unsigned i=0; i< nr_data; i++){
            if( data[i].getCodification() == "ascii" && data[i].getLocation() == "Point") {
            str.seekg( position_insert);
            readDataArray( str, data[i] ) ;

            str.seekg( data[i].getPosition() ) ;
            CopyUntilEOFInString( str, buffer, length );

            flush( str, "ascii", data[i].getName() ) ;

            position_insert = str.tellg();
            str << std::endl ;
            flush_binary( str, buffer, length) ;

            delete [] buffer ;

            };
        }; 

        for( unsigned i=0; i< nr_data; i++){
            if( data[i].getCodification() == "ascii" && data[i].getLocation() == "Cell") {
            str.seekg( position_insert);
            readDataArray( str, data[i] ) ;

            str.seekg( data[i].getPosition() ) ;
            CopyUntilEOFInString( str, buffer, length );

            flush( str, "ascii", data[i].getName() ) ;

            position_insert = str.tellg();
            str << std::endl ;
            flush_binary( str, buffer, length) ;

            delete [] buffer ;
            };
        }; 

        for( unsigned i=0; i< geometry.size(); i++){
            if( geometry[i].getCodification() == "ascii" ) {
            str.seekg( position_insert);
            readDataArray( str, geometry[i] ) ;

            str.seekg( geometry[i].getPosition() ) ;
            CopyUntilEOFInString( str, buffer, length );

            flush( str, "ascii", geometry[i].getName() ) ;

            position_insert = str.tellg();
            str << std::endl ;
            flush_binary( str, buffer, length) ;

            delete [] buffer ;
            };
        }; 

        str.seekg( temp.getPosition() ) ;


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
        str.open( fh.getName( ), std::ios::out | std::ios::in | std::ios::binary);
        str.seekg( position_insert) ;

        //str.open( "data.dat", std::ios::out | std::ios::binary);

        //Writing first point data then cell data
        for( unsigned i=0; i< nr_data; i++){
            if( data[i].getCodification() == "appended" && data[i].getLocation() == "Point") {
                if( getHeaderType() == "UInt32"){
                    uint32_t    nbytes = data[i].getNbytes() ;
                    flush_binary(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = data[i].getNbytes() ;
                    flush_binary(str, nbytes) ;
                };
                flush( str, "binary", data[i].getName() ) ;
            };
        } 
        
        for( unsigned i=0; i< nr_data; i++){
            if( data[i].getCodification() == "appended" && data[i].getLocation() == "Cell") {
                if( getHeaderType() == "UInt32"){
                    uint32_t    nbytes = data[i].getNbytes() ;
                    flush_binary(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = data[i].getNbytes() ;
                    flush_binary(str, nbytes) ;
                };
                flush( str, "binary", data[i].getName() ) ;
            };
        } 
        
        //Writing Geometry Data
        for(unsigned i=0; i<geometry.size(); i++){
            if( geometry[i].getCodification() == "appended" ) {
                if( getHeaderType() == "UInt32"){
                    uint32_t    nbytes = geometry[i].getNbytes() ;
                    flush_binary(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = geometry[i].getNbytes() ;
                    flush_binary(str, nbytes) ;
                };
                flush( str, "binary", geometry[i].getName() ) ;           
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

/*!
 * Writes data headers in strean
 * @param[in]   str     output stream
 * @param[in]   parallel     flag for parallel data headers for collection files [true/false]
 */
void VTK::writeDataHeader( std::fstream &str, bool parallel ){

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

      if( data[i].getLocation() == location){
        if(      data[i].getComponents() == 1 ) scalars <<  data[i].getName() << " " ;
        else if( data[i].getComponents() == 3 ) vectors <<  data[i].getName() << " " ;
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
      if( data[i].getLocation() == location && !parallel) writeDataArray( str, data[i] ) ;
      if( data[i].getLocation() == location &&  parallel) writePDataArray( str, data[i] ); 
    };

    str << "      </" ;
    if( parallel )  str << "P" ;

    if( location == "Point") str << "PointData> " << std::endl;
    if( location == "Cell")  str << "CellData> "  << std::endl;

  };

  return ;

};

/*!
 * Writes data array related to a field in strean
 * @param[in]   str     output stream
 * @param[in]   field_   field to be written
 */
void VTK::writeDataArray( std::fstream &str, VTKField &field_ ){

    std::string          line;

    VTKUtils::convertDataArrayToString( line, field_ ) ;

    str << line << std::endl ;
    
    str << "        </DataArray>" << std::endl ;

    return ;

};

/*!
 * Writes parallel data array related to a field in strean
 * @param[in]   str     output stream
 * @param[in]   field_   field to be written
 */
void VTK::writePDataArray( std::fstream &str, VTKField &field_ ){

    std::string          line;

    VTKUtils::convertPDataArrayToString( line, field_ ) ;

    str << line << std::endl ;
    str << "        </PDataArray>" << std::endl ;

    return ;

};

/*!
 * Reads entire VTK file (headers and data).
 */
void VTK::read( ){

    readMetaData( );
    readData( ) ;

    return ;
};

/*!
 * Reads data only from VTK file
 */
void VTK::readData( ){

  std::fstream                  str  ;
  std::fstream::pos_type        position_appended;
  std::string                   line;
  char                     c_ ;
  uint32_t                 nbytes32 ;
  uint64_t                 nbytes64 ;

  str.open( fh.getName( ), std::ios::in ) ;

  //Read appended data
  //Go to the initial position of the appended section
  while( getline(str, line) && (! Keyword_In_String( line, "<AppendedData")) ){} ;

  str >> c_;
  while( c_ != '_') str >> c_;

  position_appended = str.tellg();


  str.close();
  str.clear();

  //Open in binary for read
  str.open( fh.getName( ), std::ios::in | std::ios::binary);

  //Read fields
  for( unsigned i=0; i< nr_data; i++){
    if( data[i].getCodification() == "appended"){
      str.seekg( position_appended) ;
      str.seekg( data[i].getOffset(), std::ios::cur) ;
      if( HeaderType== "UInt32") absorb_binary( str, nbytes32 ) ;
      if( HeaderType== "UInt64") absorb_binary( str, nbytes64 ) ;
      absorb( str, "binary", data[i].getName() ) ;
    };
  };

  //Read geometry
  for(unsigned i=0; i<geometry.size(); i++){
    if( geometry[i].getCodification() == "appended"){
      str.seekg( position_appended) ;
      str.seekg( geometry[i].getOffset(), std::ios::cur) ;
      if( HeaderType== "UInt32") absorb_binary( str, nbytes32 ) ;
      if( HeaderType== "UInt64") absorb_binary( str, nbytes64 ) ;
      absorb( str, "binary", geometry[i].getName() ) ;
    };
  };

  //Read ascii data
  for( unsigned i=0; i< nr_data; i++){
    if(  data[i].getCodification() == "ascii"){
      str.seekg( data[i].getPosition() ) ;
      absorb( str, "ascii", data[i].getName() ) ;
    };
  };

  //Read geometry
  for(unsigned i=0; i<geometry.size(); i++){
    if(  geometry[i].getCodification() == "ascii"){
      str.seekg( geometry[i].getPosition() ) ;
      absorb( str, "ascii", geometry[i].getName() ) ;
    };
  };

  str.close();

  return ; 
};

/*!
 * Reads data headers from strean.
 * All field information available in file are stored.
 * @param[in]   str     output stream
 */
void VTK::readDataHeader( std::fstream &str ){


  std::fstream::pos_type   pos_ ;

  std::string                   location ;
  std::string                   line, loc_;
  std::stringstream             ss;

  bool                     read ;

  VTKField                  temp ;
  VTKField*                 ptemp ;


  for( int i=0; i<2; i++){

    ss.str("") ;
    if( i== 0) location = "Point" ;
    if( i== 1) location = "Cell" ;

    temp.setLocation( location ) ;

    ss << "</" << location << "Data>" ;
    loc_ = ss.str();
 
    read= true ; 
    if( ! getline( str, line) ) read = false ;
    if( Keyword_In_String( line, loc_) ) read=false ;
    

    while( read ){
      if( VTKUtils::convertStringToDataArray( line, temp  ) ) {

        if( temp.getCodification() == "ascii") {
            pos_ = str.tellg() ;
        }

        else{
            pos_ =  0 ; 
        };

        temp.setPosition( pos_ ) ;

        if( ! getFieldByName( temp.getName(), ptemp )) {
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

/*!
 * Reads data array from stream and stores in field information
 * @param[in]   str     output stream
 * @param[out]   field_   field information
 */
bool  VTK::readDataArray( std::fstream &str, VTKField &field_  ){

  std::string              line_ ;

  while( getline(str, line_)  ){

    if( Keyword_In_String( line_, field_.getName() ) ){
      if( VTKUtils::convertStringToDataArray( line_, field_  ) ){

        if( field_.getCodification() == "ascii") {
          field_.setPosition( str.tellg() ) ;
        };

        return true ;
      };

    };

  };
 
  return false ; 

};

/*! 
 * @} 
 */
