/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

#include <cassert>
#include <sstream>
#include <string>
#include <unordered_set>

#include "VTK.hpp"

#include "logger.hpp"

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
VTK::VTK(const std::string &dir, const std::string &name ):
     VTK(){

    setNames(dir, name) ;

}

/*!
 * set header type for appended binary output
 * @param[in] st header type ["UInt32"/"UInt64"]
 */
void  VTK::setHeaderType(const std::string &st){

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
const std::string &  VTK::getHeaderType( ) const{
    return m_headerType ;
}

/*!
 * set directory and name for VTK file
 * @param[in] dir directory of file with final "/"
 * @param[in] name file name without suffix
 */
void  VTK::setNames( const std::string &dir, const std::string &name ){

  setDirectory(dir);
  setName(name);

}

/*!
 * set name for VTK file
 * @param[in] name file name without suffix
 */
void  VTK::setName( const std::string &name ){

  m_fh.setName(name);

}

/*!
 * set directory for VTK file
 * @param[in] dir directory of file with final "/"
 */
void  VTK::setDirectory( const std::string &dir ){

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
 * Increment the time series' counter. If the time series is not active, it will
 * be activated and the counter will be reset.
 */
void  VTK::incrementCounter( ){

  if (!m_fh.isSeries()) {
    setCounter(0);
  } else {
    m_fh.incrementCounter() ;
  }

}

/*!
 * De-activates output for time series.
 * @return last value of counter
 */
int  VTK::unsetCounter( ){

  int counter = getCounter() ;
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
 * Set the specified geometry field.
 * @param[in] field is the field that will be set
 */
void VTK::setGeomData( VTKField &&field ){

    const std::string &name = field.getName();
    int id = _findFieldIndex( name, m_geometry ) ;
    if( id < 0 ){
        throw std::runtime_error("Not admissible to add geometry field with same name as user data " + name );
    }

    m_geometry[id] = std::move(field);

}

/*!
 * Add user data for input or output.
 * Codification will be set according to default value [appended] or to value
 * set by VTK::setDataCodex( VTKFormat ) or VTK::setCodex( VTKFormat )
 * @param[in] field is the field that will be added
 */
VTKField& VTK::addData( VTKField &&field ){

    const std::string &name = field.getName();
    if( findGeomData(name) ){
        throw std::runtime_error("Not admissible to add user data with same name as geometry field " + name );
    }

    int id = _findFieldIndex( name, m_data ) ;
    if (id >= 0) {
        m_data[id] = std::move(field);
    } else {
        id = m_data.size();
        m_data.push_back(std::move(field));
    }

    return m_data[id] ;

}

/*!
 * Add user data for input or output.
 * Codification will be set according to default value [appended] or to value
 * set by VTK::setDataCodex( VTKFormat ) or VTK::setCodex( VTKFormat )
 * @param[in] name name of field
 * @param[in] streamer data streamer
 */
VTKField& VTK::addData( const std::string &name, VTKBaseStreamer* streamer ){

    if( findGeomData(name) ){
        throw std::runtime_error("Not admissible to add user data with same name as geometry field " + name);
    }

    VTKField *ptr = _findData( name ) ;
    if ( ! ptr ) {
        m_data.push_back( VTKField( name ) ) ;
        ptr = &(m_data.back()) ;
    }

    ptr->setCodification(m_dataCodex) ;
    ptr->setStreamer(*streamer) ;

    return *ptr ;

}

/*!
 * Removes user data from input or output
 * @param[in] name name of field to be removed
 */
void VTK::removeData( const std::string &name ){

    std::vector<VTKField>::iterator  fieldItr = m_data.begin();

    while( fieldItr != m_data.end()){

        if( fieldItr->getName() == name){
            fieldItr = m_data.erase(fieldItr);
            return;
        } else {
            ++fieldItr;
        }
    }

    log::cout() << "did not find field for removing in VTK: " << name << std::endl;

}

/*!
 * Check if the VTK contains the specified user data
 * @param[in] name name of field
 */
bool VTK::hasData( const std::string &name ) const {

    for ( const VTKField &field : m_data ) {
        if( field.getName() == name){
            return true;
        }
    }

    return false;

}

/*!
 * Enables field for reading and writing
 * @param[in] name name of field to be enabled
 */
void VTK::enableData( const std::string &name ){

    VTKField *field = _findData(name);

    if( field ) {
        field->enable() ;

    } else{
        log::cout() << "did not find field for enabling: " << name << std::endl;
    }


}

/*!
 * Disables field for reading and writing
 * @param[in] name name of field to be disabled
 */
void VTK::disableData( const std::string &name ){

    VTKField *field = _findData(name);

    if( field ) {
        field->disable() ;

    } else{
        log::cout() << "did not find field for disabling: " << name << std::endl;
    }

}

/*!
 * Get the list of the names of the fields.
 * \return the list of the names of the fields.
 */
std::vector<std::string> VTK::getFieldNames( const std::vector<VTKField> &fields ) const {

    std::size_t nFields = fields.size();
    std::vector<std::string> names(nFields);
    for (std::size_t i = 0; i < nFields; ++i) {
        names[i] = fields[i].getName();
    }

    return names;
}

/*!
 * Get the a constant iterator pointing to the first data field.
 * \return a constant iterator pointing to the first data field.
 */
std::vector<VTKField>::const_iterator VTK::getDataBegin( ) const {

    return m_data.cbegin();
}

/*!
 * Get the a constant iterator pointing to the past-the-element data field.
 * \return a constant iterator pointing to the past-the-element data field.
 */
std::vector<VTKField>::const_iterator VTK::getDataEnd( ) const {

    return m_data.cend();
}

/*!
 * Get the a constant iterator pointing to the first geometry field.
 * \return a constant iterator pointing to the first geometry field.
 */
std::vector<VTKField>::const_iterator VTK::getGeomDataBegin( ) const {

    return m_geometry.cbegin();
}

/*!
 * Get the a constant iterator pointing to the past-the-element geometry field.
 * \return a constant iterator pointing to the past-the-element geometry field.
 */
std::vector<VTKField>::const_iterator VTK::getGeomDataEnd( ) const {

    return m_geometry.cend();
}

/*!
 * Get the number of data fields.
 * \return the number of data
 */
std::size_t VTK::getDataCount( ) const {

    return m_data.size();
}

/*!
 * Get the number of geometry fields.
 * \return the number of geometry fields.
 */
std::size_t VTK::getGeomDataCount( ) const {

    return m_geometry.size();
}

/*!
 * Finds the data field with the specified name.
 * @param[in] name name of field to be found
 * \return A constant pointer to the field with the specified name, if such a
 * field does not exists a null pointer is returned.
 */
const VTKField * VTK::findData( const std::string &name ) const {

    int id = _findFieldIndex(name, m_data) ;
    if (id >= 0) {
        return getData(id);
    } else {
        return nullptr;
    }
}

/*!
 * Finds the data field with the specified name.
 * @param[in] name name of field to be found
 * \return A pointer to the field with the specified name, if such a field does
 * not exists a null pointer is returned.
 */
VTKField * VTK::_findData( const std::string &name ) {

    int id = _findFieldIndex(name, m_data) ;
    if (id >= 0) {
        return getData(id);
    } else {
        return nullptr;
    }
}

/*!
 * Finds the geometry field with the specified name.
 * @param[in] name name of field to be found
 * \return A constant pointer to the field with the specified name, if such a
 * field does not exists a null pointer is returned.
 */
const VTKField * VTK::findGeomData( const std::string &name ) const {

    int id = _findFieldIndex(name, m_geometry) ;
    if (id >= 0) {
        return getGeomData(id);
    } else {
        return nullptr;
    }
}

/*!
 * Finds the geometry field with the specified name.
 * @param[in] name name of field to be found
 * \return A pointer to the field with the specified name, if such a field does
 * not exists a null pointer is returned.
 */
VTKField * VTK::_findGeomData( const std::string &name ) {

    int id = _findFieldIndex(name, m_geometry) ;
    if (id >= 0) {
        return getGeomData(id);
    } else {
        return nullptr;
    }
}

/*!
 * Gets a constant pointer to the data field with the specified index.
 * @param[in] id index of the field
 * \return A constant pointer to the field with the specified index.
 */
const VTKField * VTK::getData( std::size_t id ) const {

    assert(id < m_data.size());

    return m_data.data() + id;
}

/*!
 * Gets a pointer to the data field with the specified index.
 * @param[in] id index of the field
 * \return A pointer to the field with the specified index.
 */
VTKField * VTK::getData( std::size_t id ) {

    assert(id < m_data.size());

    return m_data.data() + id;
}

/*!
 * Gets a constant pointer to the geometry field with the specified index.
 * @param[in] id index of the field
 * \return A constant pointer to the field with the specified index.
 */
const VTKField * VTK::getGeomData( std::size_t id ) const {

    assert(id < m_geometry.size());

    return m_geometry.data() + id;
}

/*!
 * Gets a pointer to the geometry field with the specified index.
 * @param[in] id index of the field
 * \return A pointer to the field with the specified index.
 */
VTKField * VTK::getGeomData( std::size_t id ) {

    assert(id < m_geometry.size());

    return m_geometry.data() + id;
}

/*!
 * Finds the field with the specified name.
 * @param[in] name name of field to be found
 * @param[in] fields is the list of fields among which the requested field will
 * be searched
 * \return The index of the specified field, if such a field does not exists,
 * -1 is returned.
 */
int VTK::_findFieldIndex( const std::string &name, const std::vector<VTKField> &fields ) const {

    for( std::size_t i = 0; i < fields.size(); ++i ) {
        if( fields[i].getName() == name ){
            return i ;
        }
    }

    return -1 ;
}

/*!
 * Check if the VTK file contains data in appended format.
 * \return True if the VTK file contains data in appended format, false otherwise.
 */
bool VTK::isAppendedActive() const {

    for( auto & field : m_data ){
        if( field.isEnabled() && field.getCodification() == VTKFormat::APPENDED) {
            return true;
        }
    }

    return false;

}

/*!
 * Check if the VTK file contains data in ASCII format.
 * \return True if the VTK file contains data in ASCII format, false otherwise.
 */
bool VTK::isASCIIActive() const {

    for( auto & field : m_data ){
        if( field.isEnabled() && field.getCodification() == VTKFormat::ASCII) {
            return true;
        }
    }

    return false;

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
        if( field.isEnabled() && !field.hasAllMetaData() ) {
            field.disable() ;
#if BITPIT_ENABLE_DEBUG
            log::cout() << "Data field " << field.getName() << " has not all metadata and has been disabled from reading/writing in VTK" << std::endl ;
#endif
        }
    }

    for( auto & field : m_geometry ){
        if( field.isEnabled() && !field.hasAllMetaData() ) {
            field.disable() ;
#if BITPIT_ENABLE_DEBUG
            log::cout() << "Geometry field " << field.getName() << " has not all metadata and has been disabled from reading/writing in VTK" << std::endl ;
#endif
        }
    }

}

/*!
 * Writes entire VTK file (headers and data).
 * @param[in] writeMode if writeMode == VTKWriteMode::DEFAULT the default write setting will be used according to setCounter();
 * if writeMode == VTKWriteMode::NO_SERIES no time stamp will be added and the counter will not be increased;
 * if writeMode == VTKWriteMode::NO_INCREMENT the output file will have the same time stamp like the previous one ;
 * @param[in] time is the time associated with the file
 */
void VTK::write( VTKWriteMode writeMode, double time ){

    // Write VTK file
    write(writeMode);

    // Write time series file
    if( writeMode == VTKWriteMode::DEFAULT || writeMode == VTKWriteMode::NO_INCREMENT ){
        writeTimeSeries(time);
    }

}

/*!
 * Writes entire VTK file (headers and data).
 * @param[in] writeMode if writeMode == VTKWriteMode::DEFAULT the default write setting will be used according to setCounter();
 * if writeMode == VTKWriteMode::NO_SERIES no time stamp will be added and the counter will not be increased;
 * if writeMode == VTKWriteMode::NO_INCREMENT the output file will have the same time stamp like the previous one ;
 */
void VTK::write( VTKWriteMode writeMode ){

    // Initialize counter
    int counter = getCounter();
    if( writeMode == VTKWriteMode::NO_SERIES ){
        unsetCounter() ;
    } else if( writeMode == VTKWriteMode::NO_INCREMENT ){
        setCounter(counter - 1) ;
    }

    // Write data
    checkAllFields() ;
    calcAppendedOffsets() ;

    writeMetaInformation() ;
    writeData() ;

    // Write collection
    writeCollection() ;

    // Update counter
    if( writeMode == VTKWriteMode::DEFAULT || writeMode == VTKWriteMode::NO_INCREMENT ){
        incrementCounter() ;
    } else if( writeMode == VTKWriteMode::NO_SERIES ){
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
void VTK::write( const std::string &name, VTKWriteMode writeMode ){

    std::string oldName = getName() ;

    setName(name) ;
    write(writeMode) ;
    setName(oldName) ;
}

/*!
 * Writes entire VTK file (headers and data).
 * @param[in] name filename to be set for this output only
 * @param[in] writeMode if writeMode == VTKWriteMode::DEFAULT the default write setting will be used according to setCounter();
 * if writeMode == VTKWriteMode::NO_SERIES no time stamp will be added and the counter will not be increased;
 * if writeMode == VTKWriteMode::NO_INCREMENT the output file will have the same time stamp like the previous one ;
 * @param[in] time is the time associated with the file
 */
void VTK::write( const std::string &name, VTKWriteMode writeMode, double time ){

    std::string oldName = getName() ;

    setName(name) ;
    write(writeMode, time) ;
    setName(oldName) ;

}

/*!
 *  Writes collection file for parallel output.
 *  Is called by rank 0 in VTK::Write()
 */
void VTK::writeCollection( ) const {
    writeCollection(getName(), getName()) ;
}

/*!
 *  Writes collection file for parallel output.
 *  Is called by rank 0 in VTK::Write()
 *
 *  \param outputName filename to be set for this output only
 */
void VTK::writeCollection( const std::string &outputName ) const {
    writeCollection(outputName, outputName) ;
}

/*!
 *  Writes time series file.
 *
 *  If the file doesn't exists or if this is the first data set of the series
 *  (i.e., the counter is equal to zero), the time series will be written from
 *  scratch and will contain only the current data set. Otherwise, the existing
 *  time series will be inspected and the current data set will be either added
 *  at the end of the series or will it will replace the entry with the same
 *  file name it such an entry exists.
 *
 *  \param time is the time associated with the data set
 */
void VTK::writeTimeSeries( double time ) const {
    writeTimeSeries(getName(), getName(), time) ;
}

/*!
 *  Writes time series file.
 *
 *  If the file doesn't exists or if this is the first data set of the series
 *  (i.e., the counter is equal to zero), the time series will be written from
 *  scratch and will contain only the current data set. Otherwise, the existing
 *  time series will be inspected and the current data set will be either added
 *  at the end of the series or will it will replace the entry with the same
 *  file name it such an entry exists.
 *
 *  \param outputName filename to be set for this output only
 *  \param time is the time associated with the data set
 */
void VTK::writeTimeSeries( const std::string &outputName, double time ) const {
    writeTimeSeries(outputName, outputName, time) ;
}

/*!
 *  Writes time series file.
 *
 *  If the file doesn't exists or if this is the first data set of the series
 *  (i.e., the counter is equal to zero), the time series will be written from
 *  scratch and will contain only the current data set. Otherwise, the existing
 *  time series will be inspected and the current data set will be either added
 *  at the end of the series or will it will replace the entry with the same
 *  file name it such an entry exists.
 *
 *  \param outputName filename to be set for this output only
 *  \param seriesName series filename to be set for this output only
 *  \param time is the time associated with the data set
 */
void VTK::writeTimeSeries( const std::string &outputName, const std::string &seriesName, double time ) const {

    // Only one process should write the time series
    if (m_rank != 0) {
        return;
    }

    // Initialize series file handle
    FileHandler seriesFileHandler = m_fh ;
    seriesFileHandler.setSeries(false) ;
    seriesFileHandler.setParallel(false) ;
    seriesFileHandler.setName(seriesName) ;
    seriesFileHandler.setAppendix("pvd") ;

    // Initialize series
    if (!seriesFileHandler.exists() || getCounter() == 0) {
        std::fstream emptySeriesFileStream;
        emptySeriesFileStream.open( seriesFileHandler.getPath(), std::ios::out | std::ios::trunc ) ;
        if (!emptySeriesFileStream.is_open()) {
            throw std::runtime_error("Cannot create file \"" + seriesFileHandler.getName() + "\"" + " inside the directory \"" + seriesFileHandler.getDirectory() + "\"");
        }

        emptySeriesFileStream << "<?xml version=\"1.0\"?>" << std::endl;
        emptySeriesFileStream << "<VTKFile type=\"Collection\" version=\"0.1\">" << std::endl;
        emptySeriesFileStream << "    <Collection>" << std::endl;
        emptySeriesFileStream << "    </Collection>" << std::endl;
        emptySeriesFileStream << "</VTKFile>" << std::endl;
        emptySeriesFileStream.close();
    }

    // Create updated series
    FileHandler dataSetFileHandler;
    if (m_procs <= 1) {
        dataSetFileHandler = m_fh ;
        dataSetFileHandler.setSeries(false) ;
        dataSetFileHandler.setParallel(false) ;
        dataSetFileHandler.setName(outputName);
    } else {
        dataSetFileHandler = createCollectionHandler(outputName) ;
    }
    dataSetFileHandler.setDirectory(".");

    std::string line;
    std::stringstream series;

    std::fstream seriesFileStream;
    seriesFileStream.open( seriesFileHandler.getPath( ), std::ios::in ) ;
    if (!seriesFileStream.is_open()) {
        throw std::runtime_error("Cannot read file \"" + seriesFileHandler.getName() + "\"" + " inside the directory \"" + seriesFileHandler.getDirectory() + "\"");
    }

    while (std::getline(seriesFileStream, line)) {
        bool duplicateEntry = (line.find("file=\"" + dataSetFileHandler.getPath() + "\"") != std::string::npos);
        bool collectionEnd  = (line.find("</Collection>") != std::string::npos);

        if (collectionEnd || duplicateEntry) {
            series << "        <DataSet";
            series << " timestep=\"" << time << "\"";
            series << " file=\"" << dataSetFileHandler.getPath() << "\"";
            series << "/>" << std::endl;
        }

        if (!duplicateEntry) {
            series << line << std::endl;
        }
    }
    seriesFileStream.close();

    // Write updated series
    std::ofstream updatedSeriesFileStream;
    updatedSeriesFileStream.open( seriesFileHandler.getPath( )) ;
    if (!updatedSeriesFileStream.is_open()) {
        throw std::runtime_error("Cannot create file \"" + seriesFileHandler.getName() + "\"" + " inside the directory \"" + seriesFileHandler.getDirectory() + "\"");
    }

    updatedSeriesFileStream << series.rdbuf();
    updatedSeriesFileStream.close();
}

/*!
 *  Creates a handler for the collection.
 *
 *  \param collectionName collection filename to be set for this output only
 */
FileHandler VTK::createCollectionHandler( const std::string &collectionName ) const {
    FileHandler handler(m_fh) ;
    handler.setSeries(false) ;
    handler.setParallel(false) ;
    handler.setName(collectionName) ;
    handler.setAppendix(getCollectionExtension()) ;

    return handler;
}

/*!
 * Writes data only in VTK file
 */
void VTK::writeData( ){

    int                 length;
    char*               buffer ;


    //We need to open the file in binary mode in order for the seekg/tellg
    //methods to work properly. ASCII mode does not guarantee a correct
    //stream positioning, because of automatic jump of escape characters
    // (e.g., "\n" on Unix or "\r\n" on Windows).
    std::fstream str ;
    str.open( m_fh.getPath( ), std::ios::in | std::ios::out | std::ios::binary ) ;
    if (!str.is_open()) {
        throw std::runtime_error("Cannot create file \"" + m_fh.getName() + "\"" + " inside the directory \"" + m_fh.getDirectory() + "\"");
    }

    //Position the input stream at the beginning
    str.seekg(0, std::ios::beg);

    if (isASCIIActive()) { // Write Ascii
      //
      // It is fine to write ASCII on a binary opened file because the method
      // genericIO::flushBinary is basically a reinterpret_cast of chars into
      // chars and ensures consistent results.

        std::fstream::pos_type position_insert = str.tellg();

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

        // Position the stream at the beginning of the stream
        str.seekg(0, std::ios::beg);

    }

    if (isAppendedActive()) { // Write Appended

        char                    c_;
        std::string             line ;
        std::fstream::pos_type  position_before_write ;

        //Go to the initial position of the appended section
        while( getline(str, line) && (! bitpit::utils::string::keywordInString( line, "<AppendedData")) ){}

        str >> c_;
        while( c_ != '_') str >> c_;

        genericIO::copyUntilEOFInString( str, buffer, length );

        //Writing first point data then cell data
        for( auto &field : m_data ){
            if( field.isEnabled() && field.getCodification() == VTKFormat::APPENDED && field.getLocation() == VTKLocation::POINT ) {
                if( getHeaderType() == "UInt32"){
                    uint32_t    nbytes = static_cast<uint32_t>(calcFieldSize(field)) ;
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
                    uint32_t    nbytes = static_cast<uint32_t>(calcFieldSize(field)) ;
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
                    uint32_t    nbytes = static_cast<uint32_t>(calcFieldSize(field)) ;
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

    // Closing Appended Section
    str.close();

}

/*!
 * Writes data headers in strean
 * @param[in] str output stream
 * @param[in] parallel flag for parallel data headers for collection files [true/false]
 */
void VTK::writeDataHeader( std::fstream &str, bool parallel ) const {

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
                else if( field.getFieldType() == VTKFieldType::TENSOR ) vectors <<  field.getName() << " " ;
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
void VTK::writeDataArray( std::fstream &str, const VTKField &field ) const {

    str << vtk::convertDataArrayToString( field )  << std::endl ;
    str << "        </DataArray>" << std::endl ;

}

/*!
 * Writes parallel data array related to a field in strean
 * @param[in] str output stream
 * @param[in] field field to be written
 */
void VTK::writePDataArray( std::fstream &str, const VTKField &field ) const {

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

    str.open( m_fh.getPath( ), std::ios::in | std::ios::binary) ;

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


        //Read appended data
        for( auto & field : m_data){
            if( field.isEnabled() && field.getCodification() == VTKFormat::APPENDED){
                str.clear();
                str.seekg( position_appended) ;
                str.seekg( field.getOffset(), std::ios::cur) ;
                if( m_headerType== "UInt32") genericIO::absorbBINARY( str, nbytes32 ) ;
                if( m_headerType== "UInt64") genericIO::absorbBINARY( str, nbytes64 ) ;

#if BITPIT_ENABLE_DEBUG
                std::fstream::pos_type position_before = str.tellg();
#endif

                field.read( str, calcFieldEntries(field), calcFieldComponents(field) ) ;

#if BITPIT_ENABLE_DEBUG
                if( uint64_t(str.tellg()-position_before) != calcFieldSize(field) ){
                    log::cout() << "Warning VTK: Size of data read does not corrispond to size of field " << field.getName() << std::endl;
                    log::cout() << "Found data chunk dimension of : "<<  uint64_t(str.tellg()-position_before)<<" different from estimated field size : " <<calcFieldSize(field)<<std::endl;
                }
#endif

            }
        }

        //Read appended m_geometry
        for( auto & field : m_geometry ){
            if( field.isEnabled() && field.getCodification() == VTKFormat::APPENDED){
                str.clear();
                str.seekg( position_appended) ;
                str.seekg( field.getOffset(), std::ios::cur) ;
                if( m_headerType== "UInt32") genericIO::absorbBINARY( str, nbytes32 ) ;
                if( m_headerType== "UInt64") genericIO::absorbBINARY( str, nbytes64 ) ;

#if BITPIT_ENABLE_DEBUG
                std::fstream::pos_type position_before = str.tellg();
#endif

                field.read( str, calcFieldEntries(field), calcFieldComponents(field) ) ;

#if BITPIT_ENABLE_DEBUG
                if( uint64_t(str.tellg()-position_before) != calcFieldSize(field) ){
                    log::cout() << "Warning VTK: Size of data read does not corrispond to size of field " << field.getName() << std::endl;
                    log::cout() << "Found data chunk dimension of : "<< uint64_t(str.tellg()-position_before)<<" different from estimated field size : " <<calcFieldSize(field)<<std::endl;
                }
#endif
            }
        }

        //Leave the file opened in binary mode and reposition it to the start
        //to allow using seekg/tellg later.
        str.clear();
        str.seekg(0, std::ios::beg) ;

    }


    //Read ascii data
    for( auto & field : m_data ){
        if( field.isEnabled() &&  field.getCodification() == VTKFormat::ASCII){

            str.clear();
            str.seekg( field.getPosition() ) ;
            field.read( str, calcFieldEntries(field), calcFieldComponents(field) ) ;

#if BITPIT_ENABLE_DEBUG
            if(str.tellg()==field.getPosition()){
                log::cout() << "Warning VTK: No data have been read for field " << field.getName() << std::endl;
            }
#endif
        }
    }

    //Read ascii geometry
    for( auto & field : m_geometry ){
        if( field.isEnabled() && field.getCodification() == VTKFormat::ASCII){

            str.clear();
            str.seekg( field.getPosition() ) ;

            field.read( str, calcFieldEntries(field), calcFieldComponents(field) ) ;

#if BITPIT_ENABLE_DEBUG
            if(str.tellg()==field.getPosition()){
                log::cout() << "Warning VTK: No data have been read for field " << field.getName() << std::endl;
            }
#endif
        }
    }

    str.close();

}

/*!
 * Reads data headers from stream.
 * All field information available in file are stored.
 * All pre-existent fields not found in the header will be disabled.
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


    // Mark all fields as unused
    std::unordered_set<std::string> unusedFields;
    for (VTKField &field : m_data) {
        unusedFields.insert(field.getName());
    }

    // Read header
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

                ptemp = _findData( temp.getName() );
                if( !ptemp ) {
                    m_data.push_back( VTKField( temp ) ) ;
                }

                else{
                    ptemp->setOffset( temp.getOffset() ) ;
                    ptemp->setLocation( temp.getLocation() ) ;
                    ptemp->setDataType( temp.getDataType() ) ;
                    ptemp->setFieldType( temp.getFieldType() ) ;
                    ptemp->setCodification( temp.getCodification() ) ;
                    ptemp->setPosition( temp.getPosition() ) ;

                    // The field should not be marked as unused
                    if( unusedFields.count(ptemp->getName()) > 0 ){
                        unusedFields.erase(ptemp->getName());
                    }
                }

            }

            if( ! getline( str, line) ) read = false ;
            if( bitpit::utils::string::keywordInString( line, loc) ) read=false ;
        }

    }

    // Disable all unused fields
    for (const std::string &fieldname : unusedFields) {
        VTKField *field = _findData(fieldname);
        assert(field);
        field->disable();
    }
}

/*!
 * Reads data array information from stream and stores in field.
 * @param[in]  str output stream
 * @param[out] field field information
 */
bool VTK::readDataArray( std::fstream &str, VTKField &field  ) const {

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

/*!
 *  Gets the extension of the collection file for parallel output.
 *
 *  \result The extension of the collection file for parallel output.
 */
std::string VTK::getCollectionExtension() const {
    return "p" + getExtension();
}

}
