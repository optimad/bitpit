/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include "VTK.hpp"

namespace bitpit{
/*!
 * @ingroup VisualizationToolKit
 * @interface VTK
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

    m_procs = 1;
    m_rank  = 0;

    m_headerType = "UInt32" ;
    
    m_fh.setDirectory( "." ) ;
    m_fh.setSeries( false ) ;
    m_fh.setParallel( false ) ;

    m_geomCodex = VTKFormat::APPENDED ;
    m_dataCodex = VTKFormat::APPENDED ;

    m_cells = 0;
    m_points= 0;

}

/*! 
 * Constructor referes to a serial VTK file with appended binary data.
 * @param[in]  dir    directory of file with final "/"
 * @param[in]  name   file name without suffix
 */
VTK::VTK( std::string dir,  std::string name ):
     VTK(){

    setNames(dir, name) ;

}

/*!
 * Destructor
 */
VTK::~VTK(){

    m_data.clear();
    m_geometry.clear();

}

/*! 
 * set header type for appended binary output
 * @param[in] st header type ["UInt32"/"UInt64"]
 */
void  VTK::setHeaderType( std::string st){ 

    if( st == "UInt32" || st == "UInt64"){
        m_headerType = st ;
    }

    else{
        log::cout() << "Unsupported HeaderType " << st << std::endl ;
    }

}

/*! 
 * Get header type for appended binary output
 * @return header type ["UInt32"/"UInt64"]
 */
std::string  VTK::getHeaderType( ) const{ 
    return m_headerType ;
}

/*! 
 * set directory and name for VTK file
 * @param[in] dir directory of file with final "/"
 * @param[in] name file name without suffix
 */
void  VTK::setNames( std::string dir, std::string name ){

  setDirectory(dir);
  setName(name);

}

/*!
 * set name for VTK file
 * @param[in] name file name without suffix
 */
void  VTK::setName( std::string name ){

  m_fh.setName(name);

}

/*!
 * set directory for VTK file
 * @param[in] dir directory of file with final "/"
 */
void  VTK::setDirectory( std::string dir ){

  m_fh.setDirectory(dir);

}

/*!
 * Get the name of the VTK file
 * @return The name of the VTK file.
 */
std::string  VTK::getName() const {

  return m_fh.getName();
}

/*!
 * Get the directory where the VTK file will be saved
 * @return The directory where the VTK file will be saved.
 */
std::string  VTK::getDirectory() const {

  return m_fh.getDirectory();
}

/*!
 * Activates output for time series. sets series to true first output index to input
 * @param[in] counter first output index
 */
void  VTK::setCounter( int counter){ 

  m_fh.setSeries(true) ;
  m_fh.setCounter(counter) ;

}

/*!
 * De-activates output for time series. 
 * @return last value of counter
 */
int  VTK::unsetCounter( ){ 

  int counter = m_fh.getCounter() ;
  m_fh.setSeries(false) ;

  return counter; 
}

/*!
 * Returns the time index of the following file
 * @return counter 
 */
int  VTK::getCounter( ) const{ 

  return m_fh.getCounter( ) ;

}

/*!
 * Activates parallel output
 * @param[in] procs number of processes
 * @param[in] rank my rank
 */
void  VTK::setParallel( uint16_t procs, uint16_t rank){ 

  if( procs <  1 ) log::cout() << " Numer of processes must be greater than 0" << std::endl ;
  if( rank  >= procs) log::cout() << " m_rankess is not in valid range " << std::endl ;

  m_procs = procs; 
  m_rank  = rank; 

  if(m_procs == 0) {
   m_fh.setParallel(false) ;
  }

  else {
   m_fh.setParallel(true) ;
   m_fh.setBlock( rank ) ;
  }

}

/*!
 * sets codex for geometry and data
 * @param[in] cod codex for VTK output [VTKFormat::APPENDED/VTKFormat::ASCII]
 */
void  VTK::setCodex( VTKFormat cod) {

    setGeomCodex( cod) ;
    setDataCodex( cod) ;

}

/*!
 * sets codex for geometry only
 * @param[in] cod codex for VTK output [VTKFormat::APPENDED/VTKFormat::ASCII]
 */
void  VTK::setGeomCodex( VTKFormat cod ) {

    m_geomCodex = cod ;
    for( auto &field : m_geometry)
        field.setCodification( cod ) ;

}

/*!
 * sets codex for data only
 * @param[in] cod codex for VTK output [VTKFormat::APPENDED/VTKFormat::ASCII]
 */
void  VTK::setDataCodex( VTKFormat cod ) {

    m_dataCodex = cod ;
    for( auto &field : m_data)
        field.setCodification( cod ) ;

}

/*!
 * Add user data for input or output. 
 * Codification will be set according to default value [appended] or to value set by VTK::setDataCodex( VTKFormat ) or VTK::setCodex( VTKFormat )
 * @param[in] name name of field
 * @param[in] streamer data streamer
 */
VTKField& VTK::addData( std::string name, VTKBaseStreamer* streamer ){

    VTKField *ptr(NULL), *test(NULL) ;

    if( ! getGeomByName(name, test) ){

        if(  ! getDataByName( name, ptr ) ) {
            m_data.push_back( VTKField( name ) ) ;
            ptr = &(m_data.back()) ;
        }

        ptr->setCodification(m_dataCodex) ;
        ptr->setStreamer(*streamer) ;

    } else {
        log::cout() << "Not admissible to add user data with same name as geometry field " << name << std::endl ;
    }

    return *ptr ;

}

/*!
 * Add user data for input or output. 
 * Codification will be set according to default value [appended] or to value set by VTK::setDataCodex( VTKFormat ) or VTK::setCodex( VTKFormat )
 * @param[in] name name of field
 * @param[in] comp type of data field [ VTKFieldType::SCALAR/ VTKFieldType::VECTOR ] 
 * @param[in] loc location of data [VTKLocation::CELL/VTKLocation::POINT]
 * @param[in] type type of data [ VTKDataType::[[U]Int[8/16/32/64] / Float[32/64] ] ]
 * @param[in] streamer data streamer
 */
VTKField& VTK::addData( std::string name, VTKFieldType comp,  VTKLocation loc, VTKDataType type, VTKBaseStreamer *streamer ){

    VTKField&   field = addData( name, streamer ) ;

    field.setFieldType(comp) ;
    field.setLocation(loc) ;
    field.setDataType(type) ;

    return field ;

}

/*!
 * Removes user data from input or output 
 * @param[in] name name of field to be removed
 */
void VTK::removeData( std::string name ){

    std::vector<VTKField>::iterator  fieldItr = m_data.begin();

    while( fieldItr != m_data.end()){

        if( fieldItr->getName() == name){
            fieldItr = m_data.erase(fieldItr);
            return;
        } else {
            fieldItr++;
        }
    } 

    log::cout() << "did not find field for removing in VTK: " << name << std::endl;

}

/*!
 * Enables field for reading and writing
 * @param[in] name name of field to be enabled
 */
void VTK::enableData( std::string name ){

    VTKField* field ;

    if( getDataByName(name,field) ) {
        field->enable() ;

    } else{
        log::cout() << "did not find field for enabling: " << name << std::endl;
    }


}

/*!
 * Disables field for reading and writing
 * @param[in] name name of field to be disabled
 */
void VTK::disableData( std::string name ){

    VTKField* field ;

    if( getDataByName(name,field) ) {
        field->disable() ;

    } else{
        log::cout() << "did not find field for disabling: " << name << std::endl;
    }

}

/*!
 * Finds user data field through name. 
 * @param[in] name name of field to be found
 * @param[out] the_field pointer to the filed if found; if not found value is not altered
 * \return true if field has been found, else false
 */
bool VTK::getDataByName( const std::string &name, VTKField*& the_field ){

    the_field = NULL ;

    for( auto &field : m_data ){
        if( field.getName() == name ){
            the_field = &field ;
            return true ;
        }
    }

    return false ;
}

/*!
 * Finds geometry data field through name. 
 * @param[in] name name of field to be found
 * @param[out] the_field pointer to the filed if found; if not found value is not altered
 * \return true if field has been found, else false
 */
bool VTK::getGeomByName( const std::string &name, VTKField*& the_field ){

    the_field = NULL ;

    for( auto &field : m_geometry ){
        if( field.getName() == name ){
            the_field = &field ;
            return true ;
        }
    }

    return false ;
}

/*!
 * Calculates the offsets of all geometry and data fields for appended output.
 * Offsets are stored in each field.
 */
void VTK::calcAppendedOffsets(){

    uint64_t    offset(0) ;
    uint64_t    HeaderByte(0) ;

    if( getHeaderType() == "UInt32"){
        HeaderByte = sizeof(uint32_t) ;
    }

    else if( getHeaderType() == "UInt64") {
        HeaderByte = sizeof(uint64_t) ;
    }

    for( auto & field : m_data ){
        if( field.isEnabled() && field.getCodification() == VTKFormat::APPENDED && field.getLocation() == VTKLocation::POINT ) {
            field.setOffset( offset) ;
            offset += HeaderByte + calcFieldSize(field) ;
        }
    }

    for( auto & field : m_data ){
        if( field.isEnabled() && field.getCodification() == VTKFormat::APPENDED && field.getLocation() == VTKLocation::CELL) {
            field.setOffset( offset) ;
            offset += HeaderByte + calcFieldSize(field)  ;
        }
    }

    for( auto & field : m_geometry ){
        if( field.isEnabled() && field.getCodification() == VTKFormat::APPENDED ) {
            field.setOffset( offset) ;
            offset += HeaderByte + calcFieldSize(field)  ;
        }
    }

}

/*!
 * Check all geometry fields and user data for available metadata.
 * If metadata is missing field is disabled for writing or reading
 */
void VTK::checkAllFields(){

    for( auto & field : m_data ){
        if( field.isEnabled() && field.hasAllMetaData() ) {
            field.enable() ;
        } else {
            field.disable() ;
            log::cout() << "Data field " << field.getName() << " has not all metadata and has been disabled from reading/writing in VTK" << std::endl ;
        }
    }

    for( auto & field : m_geometry ){
        if( field.hasAllMetaData()){
            field.enable() ;
        } else {
            field.disable() ;
            log::cout() << "Geometry field " << field.getName() << " has not all metadata and has been disabled from reading/writing in VTK" << std::endl ;
        }
    }

}

/*!
 * Writes entire VTK file (headers and data).
 * @param[in] writeMode if writeMode == VTKWriteMode::DEFAULT the default write setting will be used according to setCounter();
 * if writeMode == VTKWriteMode::NO_SERIES no time stamp will be added and the counter will not be increased;
 * if writeMode == VTKWriteMode::NO_INCREMENT the output file will have the same time stamp like the previous one ;
 */
void VTK::write( VTKWriteMode writeMode ){

    int     counter(0);

    if( writeMode == VTKWriteMode::NO_SERIES ){
        counter = unsetCounter() ;
    } 

    if( writeMode == VTKWriteMode::NO_INCREMENT ){
        counter = getCounter() -1 ;
        setCounter( counter) ;
    } 

    checkAllFields() ;
    calcAppendedOffsets() ;

    writeMetaInformation() ;
    writeData() ;

    if( m_procs > 1  && m_rank == 0)  writeCollection() ;

    if( writeMode == VTKWriteMode::DEFAULT || writeMode == VTKWriteMode::NO_INCREMENT ){
        m_fh.incrementCounter() ;
    }

    if( writeMode == VTKWriteMode::NO_SERIES ){
        setCounter(counter) ;
    }

}

/*!
 * Writes entire VTK file (headers and data).
 * @param[in] name filename to be set for this output only
 * @param[in] writeMode if writeMode == VTKWriteMode::DEFAULT the default write setting will be used according to setCounter();
 * if writeMode == VTKWriteMode::NO_SERIES no time stamp will be added and the counter will not be increased;
 * if writeMode == VTKWriteMode::NO_INCREMENT the output file will have the same time stamp like the previous one ;
 */
void VTK::write( std::string name, VTKWriteMode writeMode ){

    std::string oldName = getName() ;

    setName(name) ;
    write(writeMode) ;
    setName(oldName) ;

}

/*!
 * Writes data only in VTK file
 */
void VTK::writeData( ){

    std::fstream             str ;
    std::fstream::pos_type   position_insert, position_eof ;

    int                 length;
    char*               buffer ;

    str.open( m_fh.getPath( ), std::ios::in | std::ios::out ) ;

    { // Write Ascii

        position_insert = str.tellg();
        VTKField    temp ;

        //Writing first point data then cell data
        for( auto &field : m_data ){
            if( field.isEnabled() && field.getCodification() == VTKFormat::ASCII && field.getLocation() == VTKLocation::POINT ) {
                str.seekg( position_insert);
                readDataArray( str, field ) ;

                str.seekg( field.getPosition() ) ;
                genericIO::copyUntilEOFInString( str, buffer, length );

                field.write( str ) ;
                position_insert = str.tellg();

                if(position_insert==field.getPosition()){
                    log::cout() << "Error VTK: No data has been written for field " << field.getName() << std::endl;
                    assert(false);
                }

                str << std::endl ;
                genericIO::flushBINARY( str, buffer, length) ;

                delete [] buffer ;

            }
        }

        for( auto &field : m_data ){
            if( field.isEnabled() && field.getCodification() == VTKFormat::ASCII && field.getLocation() == VTKLocation::CELL ) {
                str.seekg( position_insert);
                readDataArray( str, field ) ;

                str.seekg( field.getPosition() ) ;
                genericIO::copyUntilEOFInString( str, buffer, length );

                field.write( str ) ;
                position_insert = str.tellg();

                if(position_insert==field.getPosition()){
                    log::cout() << "Error VTK: No data has been written for field " << field.getName() << std::endl;
                    assert(false);
                }

                str << std::endl ;
                genericIO::flushBINARY( str, buffer, length) ;

                delete [] buffer ;
            }
        }

        for( auto &field : m_geometry ){
            if( field.isEnabled() && field.getCodification() == VTKFormat::ASCII ) {
                str.seekg( position_insert);
                readDataArray( str, field ) ;

                str.seekg( field.getPosition() ) ;
                genericIO::copyUntilEOFInString( str, buffer, length );

                field.write( str ) ;
                position_insert = str.tellg();

                if(position_insert==field.getPosition()){
                    log::cout() << "Error VTK: No data has been written for field " << field.getName() << std::endl;
                    assert(false);
                }

                str << std::endl ;
                genericIO::flushBINARY( str, buffer, length) ;

                delete [] buffer ;
            }
        }

        str.seekg( temp.getPosition() ) ;


    }

    { // Write Appended

        char                    c_;
        std::string             line ;
        std::fstream::pos_type  position_appended ;
        std::fstream::pos_type  position_before_write ;

        //Go to the initial position of the appended section
        while( getline(str, line) && (! bitpit::utils::string::keywordInString( line, "<AppendedData")) ){}

        str >> c_;
        while( c_ != '_') str >> c_;

        position_insert = str.tellg();
        genericIO::copyUntilEOFInString( str, buffer, length );

        str.close();
        str.clear();


        //Reopening in binary mode
        str.open( m_fh.getPath( ), std::ios::out | std::ios::in | std::ios::binary);
        str.seekg( position_insert) ;

        //str.open( "data.dat", std::ios::out | std::ios::binary);

        //Writing first point data then cell data
        for( auto &field : m_data ){
            if( field.isEnabled() && field.getCodification() == VTKFormat::APPENDED && field.getLocation() == VTKLocation::POINT ) {
                if( getHeaderType() == "UInt32"){
                    uint32_t    nbytes = calcFieldSize(field) ;
                    genericIO::flushBINARY(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = calcFieldSize(field) ;
                    genericIO::flushBINARY(str, nbytes) ;
                }

                position_before_write = str.tellg();
                field.write(str) ;
                if( (uint64_t) str.tellg()-position_before_write != calcFieldSize(field) ){
                    log::cout() << "Error VTK: Data written do not corrispond to size of field " << field.getName() << std::endl;
                    assert(false);
                }
            }
        } 

        for( auto &field : m_data ){
            if( field.isEnabled() && field.getCodification() == VTKFormat::APPENDED && field.getLocation() == VTKLocation::CELL ) {

                if( getHeaderType() == "UInt32"){
                    uint32_t    nbytes = calcFieldSize(field) ;
                    genericIO::flushBINARY(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = calcFieldSize(field) ;
                    genericIO::flushBINARY(str, nbytes) ;
                }

                position_before_write = str.tellg();
                field.write(str) ;
                if( (uint64_t) str.tellg()-position_before_write != calcFieldSize(field) ){
                    log::cout() << "Error VTK: Data written do not corrispond to size of field " << field.getName() << std::endl;
                    assert(false);
                }

            }
        } 

        //Writing Geometry Data
        for( auto &field : m_geometry ){
            if( field.isEnabled() && field.getCodification() == VTKFormat::APPENDED ) {
                if( getHeaderType() == "UInt32"){
                    uint32_t    nbytes = calcFieldSize(field) ;
                    genericIO::flushBINARY(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = calcFieldSize(field) ;
                    genericIO::flushBINARY(str, nbytes) ;
                }

                position_before_write = str.tellg();
                field.write(str) ;
                if( (uint64_t) str.tellg()-position_before_write != calcFieldSize(field) ){
                    log::cout() << "Error VTK: Data written do not corrispond to size of field " << field.getName() << std::endl;
                    assert(false);
                }
            }
        }

        genericIO::flushBINARY( str, buffer, length) ;

        delete [] buffer ;
    }

    // Closing Appended Secyion
    str.close();

}

/*!
 * Writes data headers in strean
 * @param[in] str output stream
 * @param[in] parallel flag for parallel data headers for collection files [true/false]
 */
void VTK::writeDataHeader( std::fstream &str, bool parallel ){

    std::string           line ;
    VTKLocation           location ;
    std::stringstream     scalars, vectors ;

    for( int j=0; j<2; j++){

        if( j==0 ) location = VTKLocation::POINT ;
        if( j==1 ) location = VTKLocation::CELL ;

        scalars.str("");
        vectors.str("");

        //Creating Scalar and Vector Lists
        scalars << "\"" ;
        vectors << "\"" ;

        for( auto &field : m_data ){
            if( field.isEnabled() && field.getLocation() == location){
                if(      field.getFieldType() == VTKFieldType::SCALAR ) scalars <<  field.getName() << " " ;
                else if( field.getFieldType() == VTKFieldType::VECTOR ) vectors <<  field.getName() << " " ;
            }

        }

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
        }

        str << " Scalars=" << scalars.str()
            << " Vectors=" << vectors.str()
            << ">" << std::endl;

        //Writing DataArray
        for( auto &field : m_data ){
            if( field.isEnabled() && field.getLocation() == location && !parallel) writeDataArray( str, field ) ;
            if( field.isEnabled() && field.getLocation() == location &&  parallel) writePDataArray( str, field ); 
        }

        str << "      </" ;
        if( parallel )  str << "P" ;

        if( location == VTKLocation::POINT) str << "PointData> " << std::endl;
        if( location == VTKLocation::CELL)  str << "CellData> "  << std::endl;

    }

}

/*!
 * Writes data array related to a field in strean
 * @param[in] str output stream
 * @param[in] field field to be written
 */
void VTK::writeDataArray( std::fstream &str, VTKField &field ){

    str << vtk::convertDataArrayToString( field )  << std::endl ;
    str << "        </DataArray>" << std::endl ;

}

/*!
 * Writes parallel data array related to a field in strean
 * @param[in] str output stream
 * @param[in] field field to be written
 */
void VTK::writePDataArray( std::fstream &str, VTKField &field ){

    str << vtk::convertPDataArrayToString( field ) << std::endl ;
    str << "        </PDataArray>" << std::endl ;

}

/*!
 * Reads entire VTK file (headers and data).
 */
void VTK::read( ){

    readMetaInformation( );
    checkAllFields() ;
    readData( ) ;

}

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

    str.open( m_fh.getPath( ), std::ios::in ) ;

    //Read appended data
    //Go to the initial position of the appended section
    bool foundAppendedSection = false;
    while( !foundAppendedSection && getline(str, line) ){
        foundAppendedSection = bitpit::utils::string::keywordInString( line, "<AppendedData");
    }

    if( foundAppendedSection){
        str >> c_;
        while( c_ != '_') str >> c_;

        position_appended = str.tellg();

        str.close();
        str.clear();

        //Open in binary for read
        str.open( m_fh.getPath( ), std::ios::in | std::ios::binary);

        //Read appended data
        for( auto & field : m_data){
            if( field.isEnabled() && field.getCodification() == VTKFormat::APPENDED){
                str.seekg( position_appended) ;
                str.seekg( field.getOffset(), std::ios::cur) ;
                if( m_headerType== "UInt32") genericIO::absorbBINARY( str, nbytes32 ) ;
                if( m_headerType== "UInt64") genericIO::absorbBINARY( str, nbytes64 ) ;

                std::fstream::pos_type position_before = str.tellg();

                field.read( str, calcFieldEntries(field), calcFieldComponents(field) ) ;

                if( (uint64_t) str.tellg()-position_before != calcFieldSize(field) ){
                    log::cout() << "Warning VTK: Size of data read does not corrispond to size of field " << field.getName() << std::endl;
                }

            }
        }

        //Read appended m_geometry
        for( auto & field : m_geometry ){
            if( field.isEnabled() && field.getCodification() == VTKFormat::APPENDED){
                str.seekg( position_appended) ;
                str.seekg( field.getOffset(), std::ios::cur) ;
                if( m_headerType== "UInt32") genericIO::absorbBINARY( str, nbytes32 ) ;
                if( m_headerType== "UInt64") genericIO::absorbBINARY( str, nbytes64 ) ;

                std::fstream::pos_type position_before = str.tellg();

                field.read( str, calcFieldEntries(field), calcFieldComponents(field) ) ;

                if( (uint64_t) str.tellg()-position_before != calcFieldSize(field) ){
                    log::cout() << "Warning VTK: Size of data read does not corrispond to size of field " << field.getName() << std::endl;
                }
            }
        }

        str.close();
        str.clear();

        //ReOpen in for ascii read
        str.open( m_fh.getPath( ), std::ios::in ) ;
    }

    //Read ascii data
    for( auto & field : m_data ){
        if( field.isEnabled() &&  field.getCodification() == VTKFormat::ASCII){
            str.seekg( field.getPosition() ) ;
            field.read( str, calcFieldEntries(field), calcFieldComponents(field) ) ;

            if(str.tellg()==field.getPosition()){
                log::cout() << "Warning VTK: No data have been read for field " << field.getName() << std::endl;
            }
        }
    }

    //Read ascii geometry
    for( auto & field : m_geometry ){
        if( field.isEnabled() && field.getCodification() == VTKFormat::ASCII){
            str.seekg( field.getPosition() ) ;
            field.read( str, calcFieldEntries(field), calcFieldComponents(field) ) ;

            if(str.tellg()==field.getPosition()){
                log::cout() << "Warning VTK: No data have been read for field " << field.getName() << std::endl;
            }
        }
    }

    str.close();

}

/*!
 * Reads data headers from strean.
 * All field information available in file are stored.
 * @param[in] str output stream
 */
void VTK::readDataHeader( std::fstream &str ){


    std::fstream::pos_type   pos_ ;

    VTKLocation             location;
    std::string             locationString ;
    std::string             line, loc;
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
        loc = ss.str();

        read= true ; 
        if( ! getline( str, line) ) read = false ;
        if( bitpit::utils::string::keywordInString( line, loc) ) read=false ;


        while( read ){
            if( vtk::convertStringToDataArray( line, temp  ) ) {

                if( temp.getCodification() == VTKFormat::ASCII) {
                    pos_ = str.tellg() ;
                }

                else{
                    pos_ =  0 ; 
                }

                temp.setPosition( pos_ ) ;

                if( ! getDataByName( temp.getName(), ptemp )) {
                    m_data.push_back( VTKField( temp ) ) ;
                }

                else{
                    ptemp->setOffset( temp.getOffset() ) ;
                    ptemp->setLocation( temp.getLocation() ) ;
                    ptemp->setDataType( temp.getDataType() ) ;
                    ptemp->setFieldType( temp.getFieldType() ) ;
                    ptemp->setCodification( temp.getCodification() ) ;

                }

            }

            if( ! getline( str, line) ) read = false ;
            if( bitpit::utils::string::keywordInString( line, loc) ) read=false ;
        }

    }

}

/*!
 * Reads data array from stream and stores in field information
 * @param[in]  str output stream
 * @param[out] field field information
 */
bool VTK::readDataArray( std::fstream &str, VTKField &field  ){

    std::string line ;

    while( getline(str, line)  ){

        if( bitpit::utils::string::keywordInString( line, field.getName() ) ){
            if( vtk::convertStringToDataArray( line, field  ) ){

                if( field.getCodification() == VTKFormat::ASCII) {
                    field.setPosition( str.tellg() ) ;
                }

                return true ;
            }
        }
    }

    return false ; 

}

}
