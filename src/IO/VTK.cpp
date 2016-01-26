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

    nr_procs = 1;
    my_proc  = 0;

    HeaderType = "UInt32" ;
    
    fh.setSeries( false ) ;
    fh.setParallel( false ) ;

    GeomCodex = VTKFormat::APPENDED ;
    DataCodex = VTKFormat::APPENDED ;

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
 * @param[in]  cod_    codex for VTK output [VTKFormat::APPENDED/VTKFormat::ASCII]
 */
void  VTK::setCodex( VTKFormat cod_) {

    setGeomCodex( cod_) ;
    setDataCodex( cod_) ;

    return ;
};

/*!
 * sets codex for geometry only
 * @param[in]  cod_    codex for VTK output [VTKFormat::APPENDED/VTKFormat::ASCII]
 */
void  VTK::setGeomCodex( VTKFormat cod_ ) {

    GeomCodex = cod_ ;
    for( auto &field : geometry)
        field.setCodification( cod_ ) ;
    return;
};

/*!
 * sets codex for data only
 * @param[in]  cod_    codex for VTK output [VTKFormat::APPENDED/VTKFormat::ASCII]
 */
void  VTK::setDataCodex( VTKFormat cod_ ) {

    DataCodex = cod_ ;
    for( auto &field : data)
        field.setCodification( cod_ ) ;

    return ;
};

/*!
 * Add user data for input or output. 
 * Codification will be set according to default value [appended] or to value set by VTK::setDataCodex( VTKFormat ) or VTK::setCodex( VTKFormat )
 * @param[in]  name_    name of field
 * @param[in]  comp_    type of data field [ VTKFieldType::SCALAR/ VTKFieldType::VECTOR ] 
 * @param[in]  loc_     location of data [VTKLocation::CELL/VTKLocation::POINT]
 */
VTKField* VTK::addData( std::string name_, VTKFieldType comp_,  VTKLocation loc_ ){

    int         size_ ;
    bool        allocate(true) ;
    VTKField*   ptr(NULL) ;

    if(  getFieldByName( name_, ptr ) ) {
        ptr->setComponents(comp_) ;
        ptr->setLocation(loc_) ;
    }

    else{

        if( loc_ == VTKLocation::UNDEFINED ){
            std::cout << " VTK::addData needs defined VTKLocation "  << std::endl ;
            allocate = false ;

        } else if( loc_ == VTKLocation::CELL ){
            size_ = nr_cells;

        } else if( loc_ == VTKLocation::POINT ){
            size_ = nr_points;
        }

        // Do allocation if everything ok ---------------------------------------------
        if( allocate) {
            data.push_back( VTKField( name_, comp_, loc_, VTKDataType::UNDEFINED, DataCodex, size_ ) ) ;
            ptr = &data.back() ;
        };

    };

    return ptr ;

};

/*!
 * Add user data for input or output. 
 * Codification will be set according to default value [appended] or to value set by VTK::setDataCodex( VTKFormat ) or VTK::setCodex( VTKFormat )
 * @param[in]  name_    name of field
 * @param[in]  comp_    type of data field [ VTKFieldType::SCALAR/ VTKFieldType::VECTOR ] 
 * @param[in]  loc_     location of data [VTKLocation::CELL/VTKLocation::POINT]
 * @param[in]  type_    type of data [ VTKDataType::[[U]Int[8/16/32/64] / Float[32/64] ] ]
 */
VTKField* VTK::addData( std::string name_, VTKFieldType comp_,  VTKLocation loc_, VTKDataType type_ ){

    int         size_ ;
    bool        allocate(true) ;
    VTKField*   ptr(NULL) ;

    if(  getFieldByName( name_, ptr ) ) {
        ptr->setComponents(comp_) ;
        ptr->setType(type_) ;
        ptr->setLocation(loc_) ;

    }

    else{

        if( loc_ == VTKLocation::UNDEFINED ){
            std::cout << " VTK::addData needs defined VTKLocation "  << std::endl ;
            allocate = false ;

        } else if( loc_ == VTKLocation::CELL ){
            size_ = nr_cells;

        } else if( loc_ == VTKLocation::POINT ){
            size_ = nr_points;
        }

        if( type_ == VTKDataType::UNDEFINED){
            std::cout << " VTK::addData needs defined VTKDataType "  << std::endl ;
            allocate = false ;
        }


        // Do allocation if everything ok ---------------------------------------------
        if( allocate) {
            data.push_back( VTKField( name_, comp_, loc_, type_, DataCodex, size_ ) ) ;
            ptr = &data.back() ;
        };

    };

    return ptr ;

};

/*!
 * add user data for input or output. 
 * @param[in]  name_    name of field
 * @param[in]  comp_    type of data field [ VTKFieldType::SCALAR/ VTKFieldType::VECTOR ] 
 * @param[in]  loc_     location of data [VTKLocation::CELL/VTKLocation::POINT]
 * @param[in]  type_    type of data [ VTKDataType::[[U]Int[8/16/32/64] / Float[32/64] ] ]
 * @param[in]  cod_     codification of data [VTKFormat::APPENDED/VTKFormat::ASCII]
 */
VTKField* VTK::addData( std::string name_, VTKFieldType comp_, VTKLocation loc_, VTKDataType type_, VTKFormat cod_ ){

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

        if( loc_ == VTKLocation::UNDEFINED ){
            std::cout << " VTK::addData needs defined VTKLocation "  << std::endl ;
            allocate = false ;

        } else if( loc_ == VTKLocation::CELL ){
            size_ = nr_cells;

        } else if( loc_ == VTKLocation::POINT ){
            size_ = nr_points;
        }

        if( type_ == VTKDataType::UNDEFINED){
            std::cout << " VTK::addData needs defined VTKDataType "  << std::endl ;
            allocate = false ;
        }

        // Do allocation if everything ok ---------------------------------------------
        if( allocate) {
            data.push_back( VTKField( name_, comp_, loc_, type_, cod_, size_ ) ) ;
            ptr = &data.back() ;
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

    for( auto &field : data ){
        if( field.getName() == name_ ){
            the_field = &(field) ;
            return true ;
        };
    };

    for( auto &field : geometry ){
        if( field.getName() == name_ ){
            the_field = &(field) ;
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

    uint64_t    offset(0) ;
    uint64_t    HeaderByte ;

    if( getHeaderType() == "UInt32"){
        HeaderByte = sizeof(uint32_t) ;
    }

    else if( getHeaderType() == "UInt64") {
        HeaderByte = sizeof(uint64_t) ;
    };


    for( auto & field : data ){
        if( field.getCodification() == VTKFormat::APPENDED && field.getLocation() == VTKLocation::POINT ) {
            field.setOffset( offset) ;
            offset += HeaderByte + field.getNbytes()  ;
        };
    };

    for( auto & field : data ){
        if( field.getCodification() == VTKFormat::APPENDED && field.getLocation() == VTKLocation::CELL) {
            field.setOffset( offset) ;
            offset += HeaderByte + field.getNbytes()  ;
        };
    };

    for( auto & field : geometry ){
        if( field.getCodification() == VTKFormat::APPENDED ) {
            field.setOffset( offset) ;
            offset += HeaderByte + field.getNbytes()  ;
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
        for( auto &field : data ){
            if( field.getCodification() == VTKFormat::ASCII && field.getLocation() == VTKLocation::POINT ) {
                str.seekg( position_insert);
                readDataArray( str, field ) ;

                str.seekg( field.getPosition() ) ;
                CopyUntilEOFInString( str, buffer, length );

                flushData( str, VTKFormat::ASCII, field.getName() ) ;

                position_insert = str.tellg();
                str << std::endl ;
                flush_binary( str, buffer, length) ;

                delete [] buffer ;

            };
        }; 

        for( auto &field : data ){
            if( field.getCodification() == VTKFormat::ASCII && field.getLocation() == VTKLocation::CELL ) {
                str.seekg( position_insert);
                readDataArray( str, field ) ;

                str.seekg( field.getPosition() ) ;
                CopyUntilEOFInString( str, buffer, length );

                flushData( str, VTKFormat::ASCII, field.getName() ) ;

                position_insert = str.tellg();
                str << std::endl ;
                flush_binary( str, buffer, length) ;

                delete [] buffer ;
            };
        }; 

        for( auto &field : geometry ){
            if( field.getCodification() == VTKFormat::ASCII ) {
                str.seekg( position_insert);
                readDataArray( str, field ) ;

                str.seekg( field.getPosition() ) ;
                CopyUntilEOFInString( str, buffer, length );

                flushData( str, VTKFormat::ASCII, field.getName() ) ;

                position_insert = str.tellg();
                str << std::endl ;
                flush_binary( str, buffer, length) ;

                delete [] buffer ;
            };
        }; 

        str.seekg( temp.getPosition() ) ;


    }

    { // Write Appended

        char                    c_;
        std::string             line ;
        std::fstream::pos_type  position_appended ;

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
        for( auto &field : data ){
            if( field.getCodification() == VTKFormat::APPENDED && field.getLocation() == VTKLocation::POINT ) {
                if( getHeaderType() == "UInt32"){
                    uint32_t    nbytes = field.getNbytes() ;
                    flush_binary(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = field.getNbytes() ;
                    flush_binary(str, nbytes) ;
                };
                flushData( str, VTKFormat::APPENDED, field.getName() ) ;
            };
        } 

        for( auto &field : data ){
            if( field.getCodification() == VTKFormat::APPENDED && field.getLocation() == VTKLocation::CELL ) {

                if( getHeaderType() == "UInt32"){
                    uint32_t    nbytes = field.getNbytes() ;
                    flush_binary(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = field.getNbytes() ;
                    flush_binary(str, nbytes) ;
                };
                flushData( str, VTKFormat::APPENDED, field.getName() ) ;

            };
        } 

        //Writing Geometry Data
        for( auto &field : geometry ){
            if( field.getCodification() == VTKFormat::APPENDED ) {
                if( getHeaderType() == "UInt32"){
                    uint32_t    nbytes = field.getNbytes() ;
                    flush_binary(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = field.getNbytes() ;
                    flush_binary(str, nbytes) ;
                };
                flushData( str, VTKFormat::APPENDED, field.getName() ) ;           
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

    std::string           line ;
    VTKLocation           location ;
    std::stringstream     scalars, vectors ;

    for( int j=0; j<2; j++){

        if( j==0 ) location = VTKLocation::POINT ;
        if( j==1 ) location = VTKLocation::CELL ;

        //    scalars.clear();
        //    vectors.clear();

        scalars.str("");
        vectors.str("");

        //Creating Scalar and Vector Lists
        scalars << "\"" ;
        vectors << "\"" ;

        for( auto &field : data ){
            if( field.getLocation() == location){
                if(      field.getComponents() == VTKFieldType::SCALAR ) scalars <<  field.getName() << " " ;
                else if( field.getComponents() == VTKFieldType::VECTOR ) vectors <<  field.getName() << " " ;
            };

        };

        scalars << "\"" ;
        vectors << "\"" ;

        if(      location == VTKLocation::POINT) {
            str << "      <" ;
            if( parallel )  str << "P" ;
            str << "PointData " ;
        }

        else if( location == VTKLocation::CELL )  {
            str << "      <" ;
            if( parallel )  str << "P" ;
            str << "CellData " ;
        };

        str << " Scalars=" << scalars.str()
            << " Vectors=" << vectors.str()
            << ">" << std::endl;

        //Writing DataArray
        for( auto &field : data ){
            if( field.getLocation() == location && !parallel) writeDataArray( str, field ) ;
            if( field.getLocation() == location &&  parallel) writePDataArray( str, field ); 
        };

        str << "      </" ;
        if( parallel )  str << "P" ;

        if( location == VTKLocation::POINT) str << "PointData> " << std::endl;
        if( location == VTKLocation::CELL)  str << "CellData> "  << std::endl;

    };

    return ;

};

/*!
 * Writes data array related to a field in strean
 * @param[in]   str     output stream
 * @param[in]   field_   field to be written
 */
void VTK::writeDataArray( std::fstream &str, VTKField &field_ ){

    str << VTKUtils::convertDataArrayToString( field_ )  << std::endl ;
    str << "        </DataArray>" << std::endl ;

    return ;

};

/*!
 * Writes parallel data array related to a field in strean
 * @param[in]   str     output stream
 * @param[in]   field_   field to be written
 */
void VTK::writePDataArray( std::fstream &str, VTKField &field_ ){

    str << VTKUtils::convertPDataArrayToString( field_ ) << std::endl ;
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

    std::fstream              str  ;
    std::fstream::pos_type    position_appended;
    std::string               line;
    char                      c_ ;
    uint32_t                  nbytes32 ;
    uint64_t                  nbytes64 ;

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
    for( auto & field : data){
        if( field.getCodification() == VTKFormat::APPENDED){
            str.seekg( position_appended) ;
            str.seekg( field.getOffset(), std::ios::cur) ;
            if( HeaderType== "UInt32") absorb_binary( str, nbytes32 ) ;
            if( HeaderType== "UInt64") absorb_binary( str, nbytes64 ) ;
            absorbData( str, VTKFormat::APPENDED, field.getName() ) ;
        };
    };

    //Read geometry
    for( auto & field : geometry ){
        if( field.getCodification() == VTKFormat::APPENDED){
            str.seekg( position_appended) ;
            str.seekg( field.getOffset(), std::ios::cur) ;
            if( HeaderType== "UInt32") absorb_binary( str, nbytes32 ) ;
            if( HeaderType== "UInt64") absorb_binary( str, nbytes64 ) ;
            absorbData( str, VTKFormat::APPENDED, field.getName() ) ;
        };
    };

    //Read ascii data
    for( auto & field : data ){
        if(  field.getCodification() == VTKFormat::ASCII){
            str.seekg( field.getPosition() ) ;
            absorbData( str, VTKFormat::ASCII, field.getName() ) ;
        };
    };

    //Read geometry
    for( auto & field : geometry ){
        if( field.getCodification() == VTKFormat::ASCII){
            str.seekg( field.getPosition() ) ;
            absorbData( str, VTKFormat::ASCII, field.getName() ) ;
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

    VTKLocation             location;
    std::string             locationString ;
    std::string             line, loc_;
    std::stringstream       ss;

    bool                    read ;

    VTKField                temp ;
    VTKField*               ptemp ;


    for( int i=0; i<2; i++){

        ss.str("") ;
        if( i== 0) {
            location = VTKLocation::POINT;
            locationString = "Point" ;
        } else if( i== 1) {
            location = VTKLocation::CELL;
            locationString = "Cell" ;
        }

        temp.setLocation( location ) ;

        ss << "</" << locationString << "Data>" ;
        loc_ = ss.str();

        read= true ; 
        if( ! getline( str, line) ) read = false ;
        if( Keyword_In_String( line, loc_) ) read=false ;


        while( read ){
            if( VTKUtils::convertStringToDataArray( line, temp  ) ) {

                if( temp.getCodification() == VTKFormat::ASCII) {
                    pos_ = str.tellg() ;
                }

                else{
                    pos_ =  0 ; 
                };

                temp.setPosition( pos_ ) ;

                if( ! getFieldByName( temp.getName(), ptemp )) {
                    data.push_back( temp ) ;
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

                if( field_.getCodification() == VTKFormat::ASCII) {
                    field_.setPosition( str.tellg() ) ;
                };

                return true ;
            };

        };

    };

    return false ; 

};

/*!
 * Interface method for writing field data
 * @param[in] str output stream
 * @param[in] format codex to be used [VTKFormat::ASCII/VKFormat::APPENDED]
 * @param[in] name name of the field to be written
 */
void VTK::flushData( std::fstream &str, VTKFormat format, std::string name ){
    BITPIT_UNUSED( str ) ;
    BITPIT_UNUSED( format ) ;
    BITPIT_UNUSED( name ) ;
};

/*!
 * Interface method for reading field data
 * @param[in] str output stream
 * @param[in] format codex to be used [VTKFormat::ASCII/VKFormat::APPENDED]
 * @param[in] name name of the field to be read
 */
void VTK::absorbData( std::fstream &str, VTKFormat format, std::string name ){
    BITPIT_UNUSED( str ) ;
    BITPIT_UNUSED( format ) ;
    BITPIT_UNUSED( name ) ;
};
/*! 
 * @} 
 */
