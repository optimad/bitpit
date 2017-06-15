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
 * @class VTKRectilinearGrid
 * @ingroup VisualizationToolKit
 * @brief VTK input output for Rectilinear Meshes
 *
 * VTKRectilinearGrid provides methods to read and write parallel and serial rectlinear meshes and data. 
 * The class is agnostic with respect to the container used for the data and provides an interface through the CRTP mechanism.
 * The numbering of nodes start with 0. Different numbering scheme is not supported.
 *
 */

/*!  
 *  Default constructor.
 *  Allocates three geometry fields called "x_Coord", "y_Coord" and "z_Coord".
 */
VTKRectilinearGrid::VTKRectilinearGrid( ) :VTK() {

    m_fh.setAppendix( "vtr" );

    m_geometry.push_back( VTKField("x_Coord") ) ;
    m_geometry.push_back( VTKField("y_Coord") ) ;
    m_geometry.push_back( VTKField("z_Coord") ) ;

    for( auto & field : m_geometry ){
        field.setLocation( VTKLocation::POINT ) ;
        field.setFieldType( VTKFieldType::KNOWN_BY_CLASS ) ;
        field.setDataType( VTKDataType::Float64 ) ;
        field.setCodification(m_geomCodex);
    }

}

/*!  
 *  Constructor for parallel 3D grid.
 *  Calls default constructor and sets provided input information.
 *  @param[in] dir directory of VTK file with final "/"
 *  @param[in] name name of VTK without suffix
 */
VTKRectilinearGrid::VTKRectilinearGrid( std::string dir, std::string name ) :VTKRectilinearGrid( ) {

    setNames( dir, name ) ;

}

/*!  
 *  Constructor for parallel 3D grid.
 *  Calls default constructor and sets provided input information.
 *  @param[in] dir directory of VTK file with final "/"
 *  @param[in] name name of VTK without suffix
 *  @param[in] codex codex of data [VTKFormat::ASCII/VTKFormat::APPENDED]
 *  @param[in] n1 min node index in first direction 
 *  @param[in] n2 max node index in first direction 
 *  @param[in] m1 min node index in second direction 
 *  @param[in] m2 max node index in second direction 
 *  @param[in] l1 min node index in third direction 
 *  @param[in] l2 max node index in third direction 
 */
VTKRectilinearGrid::VTKRectilinearGrid( std::string dir, std::string name, VTKFormat codex, int n1, int n2, int m1, int m2, int l1, int l2 ) :VTKRectilinearGrid( ) {

    setNames( dir, name ) ;

    setDimensions( n1, n2, m1, m2, l1, l2) ;
    setGeomCodex( codex ) ;

}

/*!  
 *  Constructor for serial 3D grids.
 *  Min index is set automatically to zero.
 *  Calls default constructor and sets provided input information.
 *  @param[in] dir directory of VTK file with final "/"
 *  @param[in] name name of VTK without suffix
 *  @param[in] codex codex of data [VTKFormat::ASCII/VTKFormat::APPENDED]
 *  @param[in] n number of nodes in first direction 
 *  @param[in] m number of nodes in second direction 
 *  @param[in] l number of nodes in third direction 
 */
VTKRectilinearGrid::VTKRectilinearGrid( std::string dir, std::string name, VTKFormat codex, int n, int m, int l ) : VTKRectilinearGrid( ) {

    setNames( dir, name ) ;

    setDimensions( 0, n-1, 0, m-1, 0, l-1) ;
    setGeomCodex( codex ) ;

}

/*!  
 *  Constructor for parallel 2D grid.
 *  Calls default constructor and sets provided input information.
 *  @param[in] dir directory of VTK file with final "/"
 *  @param[in] name name of VTK without suffix
 *  @param[in] codex codex of data [VTKFormat::ASCII/VTKFormat::APPENDED]
 *  @param[in] n1 min node index in first direction 
 *  @param[in] n2 max node index in first direction 
 *  @param[in] m1 min node index in second direction 
 *  @param[in] m2 max node index in second direction 
 */
VTKRectilinearGrid::VTKRectilinearGrid( std::string dir, std::string name, VTKFormat codex, int n1, int n2, int m1, int m2 ) :VTKRectilinearGrid( ) {

    setNames( dir, name ) ;

    setDimensions( n1, n2, m1, m2, 0, 0) ;
    setGeomCodex( codex ) ;

}

/*!  
 *  Constructor for serial 2D grid.
 *  Calls default constructor and sets provided input information.
 *  @param[in] dir directory of VTK file with final "/"
 *  @param[in] name name of VTK without suffix
 *  @param[in] codex codex of data [VTKFormat::ASCII/VTKFormat::APPENDED]
 *  @param[in] n number of nodes in first direction 
 *  @param[in] m number of nodes in second direction 
 */
VTKRectilinearGrid::VTKRectilinearGrid( std::string dir, std::string name, VTKFormat codex, int n, int m ) : VTKRectilinearGrid( ) {

    setNames( dir, name ) ;

    setDimensions( 0, n-1, 0, m-1, 0, 0) ;
    setGeomCodex( codex ) ;

}

/*!  
 *  Destructor 
 */
VTKRectilinearGrid::~VTKRectilinearGrid( ){
}

/*!  
 *  Reads meta data of VTR file (grid size, data fields, codex, position of data within file).
 *  Calls setDimension.
 */
void VTKRectilinearGrid::readMetaInformation( ){

    std::fstream str;
    std::string line, temp;

    std::fstream::pos_type        position;
    std::array<int,6>             extensions ;


    str.open( m_fh.getPath( ), std::ios::in ) ;

    getline( str, line);
    while( ! bitpit::utils::string::keywordInString( line, "<VTKFile")){
        getline(str, line);
    }

    if( bitpit::utils::string::getAfterKeyword( line, "header_type", '\"', temp) ){
        setHeaderType( temp) ;
    }

    while( ! bitpit::utils::string::keywordInString( line, "<Piece")){
        getline(str, line);
    }

    bitpit::utils::string::getAfterKeyword( line, "Extent", '\"', temp) ;
    bitpit::utils::string::convertString( temp, extensions );

    m_localIndex[0][0] = extensions[0] ;
    m_localIndex[0][1] = extensions[1] ;
    m_localIndex[1][0] = extensions[2] ;
    m_localIndex[1][1] = extensions[3] ;
    m_localIndex[2][0] = extensions[4] ;
    m_localIndex[2][1] = extensions[5] ;

    position = str.tellg() ;

    readDataHeader( str ) ;

    for( auto &field : m_geometry ){ //int i=0; i<geometry.size(); ++i){
        str.seekg( position) ;
        if( ! readDataArray( str, field ) ) {
            log::cout() << field.getName() << " DataArray not found" << std::endl ;
        }
    }


    setDimensions( m_localIndex[0][0], m_localIndex[0][1], m_localIndex[1][0], m_localIndex[1][1], m_localIndex[2][0], m_localIndex[2][1] ) ;
    str.close() ; 

}

/*!  
 *  Writes entire VTR but the data.
 */
void VTKRectilinearGrid::writeMetaInformation( ){

    std::fstream str;

    str.open( m_fh.getPath( ), std::ios::out ) ;

    //Writing XML header
    str << "<?xml version=\"1.0\"?>" << std::endl;

    //Writing Piece Information
    str << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\"  header_type=\"" << m_headerType << "\">" << std::endl; 
    str << "  <RectilinearGrid WholeExtent= \"" 
        << m_globalIndex[0][0] << " " << m_globalIndex[0][1]<< " "
        << m_globalIndex[1][0] << " " << m_globalIndex[1][1]<< " "
        << m_globalIndex[2][0] << " " << m_globalIndex[2][1]<< " "
        << "\" >" << std::endl;

    str << "    <Piece Extent= \" " 
        << m_localIndex[0][0] << " " << m_localIndex[0][1]<< " "
        << m_localIndex[1][0] << " " << m_localIndex[1][1]<< " "
        << m_localIndex[2][0] << " " << m_localIndex[2][1]<< " "
        << "\" >" << std::endl;



    //Header for Data
    writeDataHeader( str, false ) ;

    //Wring Geometry Information   
    str << "       <Coordinates>" << std::endl;
    writeDataArray( str, m_geometry[0] ) ;
    writeDataArray( str, m_geometry[1] ) ;
    writeDataArray( str, m_geometry[2] ) ;
    str << "       </Coordinates>" << std::endl;

    //Closing Piece
    str << "    </Piece>" << std::endl;
    str << "  </RectilinearGrid>" << std::endl;

    //Write Appended Section
    str << "  <AppendedData encoding=\"raw\">" << std::endl;
    str << "_" ;
    str << std::endl ;

    //Closing XML
    str << "</VTKFile>" << std::endl;

    str.close() ;

}

/*!  
 *  Writes collection file for parallel output. 
 *  Is called by rank 0 in VTK::Write()
 */
void VTKRectilinearGrid::writeCollection( ){

    std::fstream str ;

    FileHandler     fhp, fho ;

    fhp = m_fh ;
    fho = m_fh ;

    fhp.setParallel(false) ;
    fhp.setAppendix("pvtr") ;

    fho.setDirectory(".") ;

    str.open( fhp.getPath( ), std::ios::out ) ;

    //Writing XML header
    str << "<?xml version=\"1.0\"?>" << std::endl;

    //Writing Piece Information
    str << "<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    str << "  <PRectilinearGrid WholeExtent= \"" 
        << m_globalIndex[0][0] << " " << m_globalIndex[0][1]<< " "
        << m_globalIndex[1][0] << " " << m_globalIndex[1][1]<< " "
        << m_globalIndex[2][0] << " " << m_globalIndex[2][1]<< " "
        << "GhostLevel=\"0\">" << std::endl;



    //Header for Data
    writeDataHeader( str, true );

    //Wring Geometry Information
    str << "      <PCoordinates>" << std::endl;
    writePDataArray( str, m_geometry[0] ) ;
    writePDataArray( str, m_geometry[1] ) ;
    writePDataArray( str, m_geometry[2] ) ;
    str << "      </PCoordinates>" << std::endl;


    for( int i=0; i<m_procs; i++){
        fho.setBlock(i) ;
        extension3D_t &index = m_procIndex[i] ;

        str << "    <Piece Extent= \" " 
            << index[0][0] << " " << index[0][1] << " "
            << index[1][0] << " " << index[1][1] << " "
            << index[2][0] << " " << index[2][1] << " "
            << "\" Source= \"" << fho.getPath() << "\"/>" << std::endl;
    }

    str << "  </PRectilinearGrid>"  << std::endl;
    str << "</VTKFile>" << std::endl;

    str.close() ;

}

/*!  
 *  sets the dimension for 3D parallel grids.
 *  @param[in] n1 min node index in first direction 
 *  @param[in] n2 max node index in first direction 
 *  @param[in] m1 min node index in second direction 
 *  @param[in] m2 max node index in second direction 
 *  @param[in] l1 min node index in third direction 
 *  @param[in] l2 max node index in third direction 
 */
void VTKRectilinearGrid::setDimensions( int n1, int n2, int m1, int m2, int l1, int l2 ){

    m_dimensions = ( l1 == l2 ) ? 2 :3 ;

    m_localIndex[0][0] = n1 ;
    m_localIndex[0][1] = n2 ;
    m_localIndex[1][0] = m1 ;
    m_localIndex[1][1] = m2 ;
    m_localIndex[2][0] = l1 ;
    m_localIndex[2][1] = l2 ;

    if(m_procs==1)
        m_globalIndex = m_localIndex ;

    m_cells = 1 ;
    m_points = 1 ;

    for(int d=0; d<m_dimensions; ++d){
        m_cells = m_cells *  ( m_localIndex[d][1] -m_localIndex[d][0] ) ;
        m_points = m_points *  ( m_localIndex[d][1] -m_localIndex[d][0] +1 ) ;
    }

}

/*!  
 *  sets the dimension for 3D serial grids.
 *  @param[in] n number of nodes in first direction 
 *  @param[in] m number of nodes in second direction 
 *  @param[in] l number of nodes in third direction 
 */
void VTKRectilinearGrid::setDimensions( int n, int m, int l ){

    this->setDimensions( 0, n-1, 0, m-1, 0, l-1 );

}

/*!  
 *  sets the dimension for 2D parallel grids.
 *  @param[in] n1 min node index in first direction 
 *  @param[in] n2 max node index in first direction 
 *  @param[in] m1 min node index in second direction 
 *  @param[in] m2 max node index in second direction 
 */
void VTKRectilinearGrid::setDimensions( int n1, int n2, int m1, int m2 ){

    this->setDimensions( n1, n2, m1, m2, 0, 0 );

}

/*!  
 *  sets the dimension for 2D serial grids.
 *  @param[in] n number of nodes in first direction 
 *  @param[in] m number of nodes in second direction 
 */
void VTKRectilinearGrid::setDimensions( int n, int m ){

    this->setDimensions( 0, n-1, 0, m-1, 0, 0 );

}

/*!  
 *  sets the global 3D grid information for parallel output.
 *  Needs to be called by all processes
 *  @param[in] I numer of nodes in first direction
 *  @param[in] J numer of nodes in second direction
 *  @param[in] K numer of nodes in third direction
 */
void VTKRectilinearGrid::setGlobalDimensions( int I, int J, int K ){

    m_globalIndex[0][0] = 0;
    m_globalIndex[0][1] = I;
    m_globalIndex[1][0] = 0;
    m_globalIndex[1][1] = J;
    m_globalIndex[2][0] = 0;
    m_globalIndex[2][1] = K;

}

/*!
 * Associates streamer to a geometrical field
 * @param[in] fieldEnum which geometrical field 
 * @param[in] streamer VTKBaseStreamer
 */
void VTKRectilinearGrid::setGeomData( VTKRectilinearField fieldEnum, VTKBaseStreamer *streamer ){

    int      index = static_cast<int>(fieldEnum) ;
    VTKField& field = m_geometry[index] ;

    field.setStreamer( *streamer ) ;

}

/*!
 * Associates streamer to a geometrical field
 * @param[in] fieldEnum which geometrical field 
 * @param[in] type type of data [ VTKDataType::[[U]Int[8/16/32/64] / Float[32/64] ] ]
 * @param[in] streamer VTKBaseStreamer
 */
void VTKRectilinearGrid::setGeomData( VTKRectilinearField fieldEnum, VTKDataType type, VTKBaseStreamer *streamer ){

    int      index = static_cast<int>(fieldEnum) ;
    VTKField& field = m_geometry[index] ;

    field.setDataType( type ) ;
    field.setStreamer( *streamer ) ;

}

/*!  
 *  sets the global 2D grid information for parallel output.
 *  Needs to be called by all processes
 *  @param[in] I numer of nodes in first direction
 *  @param[in] J numer of nodes in second direction
 */
void VTKRectilinearGrid::setGlobalDimensions( int I, int J ){

    m_globalIndex[0][0] = 0;
    m_globalIndex[0][1] = I;
    m_globalIndex[1][0] = 0;
    m_globalIndex[1][1] = J;
    m_globalIndex[2][0] = 0;
    m_globalIndex[2][1] = 0;

}

/*!  
 *  sets the global 3D grid information for parallel output.
 *  Needs to be called only by rank 0.
 *  @param[in] loc min I-index, max I-index, min J-index, max J-index, min K-index, max K-index for each process
 */
void VTKRectilinearGrid::setGlobalIndex( std::vector<extension3D_t> loc ){

    if( loc.size() != m_procs ) 
        log::cout() << "Size of loc_ in VTKRectilinearGrid::setParallelIndex does not fit m_procs " << std::endl ;

    m_procIndex   = loc ;

}

/*!  
 *  sets the global 2D grid information for parallel output.
 *  Needs to be called only by rank 0.
 *  @param[in]  loc min I-index, max I-index, min J-index, max J-index for each process
 */
void VTKRectilinearGrid::setGlobalIndex( std::vector<extension2D_t> loc ){

    if( loc.size() !=m_procs ) log::cout() << "Size of loc_ in VTKRectilinearGrid::setParallelIndex does not fit m_procs " << std::endl ;

    m_procIndex.resize(m_procs) ;

    for( int i=0; i< m_procs; ++i){
        m_procIndex[i][0]   = loc[i][0] ;
        m_procIndex[i][1]   = loc[i][1] ;
        m_procIndex[i][2]   = {{0,0}} ;
    }

}

/*!
 * Calculates the size (in bytes) of a field
 * @param[in] field field 
 * @return size of the field
 */
uint64_t VTKRectilinearGrid::calcFieldSize( const VTKField &field ){

    uint64_t bytes = calcFieldEntries(field) ;
    bytes *= VTKTypes::sizeOfType( field.getDataType() ) ;

    return bytes ;

}

/*!
 * Calculates the number of entries of a field
 * @param[in] field field 
 * @return size of the field
 */
uint64_t VTKRectilinearGrid::calcFieldEntries( const VTKField &field ){

    uint64_t entries(0) ;
    std::string name( field.getName() ) ;

    if( name == "x_Coord" ){
        entries = m_localIndex[0][1] -m_localIndex[0][0] +1 ;

    } else if( name == "y_Coord" ){
        entries = m_localIndex[1][1] -m_localIndex[1][0] +1 ;

    } else if( name == "z_Coord" ){
        entries = m_localIndex[2][1] -m_localIndex[2][0] +1 ;

    } else{

        VTKLocation location( field.getLocation() ) ;
        assert( location != VTKLocation::UNDEFINED) ;

        if( location == VTKLocation::CELL ){
            entries = m_cells ;

        } else if( location == VTKLocation::POINT ){
            entries = m_points ;

        }

        VTKFieldType fieldType( field.getFieldType() ) ;
        assert( fieldType != VTKFieldType::UNDEFINED) ;

        entries *= static_cast<uint64_t>(fieldType) ;

    }

    return entries ;

}

/*!
 * Calculates the compnents of a field
 * @param[in] field field 
 * @return size of the field
 */
uint8_t VTKRectilinearGrid::calcFieldComponents( const VTKField &field ){

    uint8_t comp ;
    std::string name( field.getName() ) ;

    if( name == "x_Coord" || name == "y_Cooord" || name == "z_Coord" ){
        comp = 1 ;

    } else{
        VTKFieldType fieldType( field.getFieldType() ) ;
        assert( fieldType != VTKFieldType::UNDEFINED) ;

        comp = static_cast<uint8_t>(fieldType) ;

    }

    return comp ;

}

}
