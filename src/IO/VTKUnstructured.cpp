/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by by the Free Software Foundation.
 *
 *  BitPit is distributed in the hope that it will be useful, but WITHOUT
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
 * @ingroup     VisualizationToolKit
 * @{
 *
 * @class       VTKUnstructuredGrid
 * @brief       VTK input output for Unstructured Meshes
 *
 * VTKUnstructuredGrid provides methods to read and write parallel and serial unstructured meshes and data. 
 * The class is agnostic with respect to the container used for the data and provides an interface through the CRTP mechanism.
 *
 */

/*!  
 *  Default constructor.
 *  Allocates four geometry fields called "Points"(Float64), "offsets"(Int32), "types"(Int32) and "connectivity"(Int32).
 */
VTKUnstructuredGrid::VTKUnstructuredGrid( ) :VTK() {

  fh.setAppendix("vtu");

  geometry.push_back( VTKField( "Points",       VTKFieldType::VECTOR, VTKLocation::POINT, VTKDataType::Float64) ) ;
  geometry.push_back( VTKField( "offsets",      VTKFieldType::SCALAR, VTKLocation::CELL, VTKDataType::Int32   ) ) ;
  geometry.push_back( VTKField( "types",        VTKFieldType::SCALAR, VTKLocation::CELL, VTKDataType::Int32   ) ) ;
  geometry.push_back( VTKField( "connectivity", VTKFieldType::SCALAR, VTKLocation::CELL, VTKDataType::Int32   ) ) ;

};

/*!  
 *  Constructor.
 *  sets input parameters and calls default constructor
 *  @param[in]  dir_        Directory of vtk file with final "/"
 *  @param[in]  name_       Name of vtk file without suffix
 */
VTKUnstructuredGrid::VTKUnstructuredGrid( std::string dir_, std::string name_ ):VTKUnstructuredGrid( ){

  setNames( dir_, name_ ) ; 
  return ;

};

/*!  
 *  Destructor.
 */
VTKUnstructuredGrid::~VTKUnstructuredGrid( ) {
    
};

/*!  
 *  sets the type of the geometry variables
 *  @param[in]  Ptype       Type of "Point" geometry information [ VTKDataType::Float[32/64]]
 *  @param[in]  Otype       Type of "offset" geometry information [ VTKDataType::[U]Int[8/16/32/64] ] 
 *  @param[in]  Ttype       Type of "types" geometry information [ VTKDataType::[U]Int[8/16/32/64] ]
 *  @param[in]  Ctype       Type of "connectivity" geometry information [ VTKDataType::[U]Int[8/16/32/64] ]
 */
void VTKUnstructuredGrid::setGeomTypes( VTKDataType Ptype, VTKDataType Otype, VTKDataType Ttype, VTKDataType Ctype  ){

    geometry[0].setType(Ptype) ;
    geometry[1].setType(Otype) ;
    geometry[2].setType(Ttype) ;
    geometry[3].setType(Ctype) ;

    return ;
};

/*!  
 *  sets the size of the unstructured grid.
 *  @param[in]  ncells_     number of cells
 *  @param[in]  npoints_    number of points
 *  @param[in]  nconn_      size of the connectivity information
 */
void VTKUnstructuredGrid::setDimensions( uint64_t ncells_, uint64_t npoints_, uint64_t nconn_ ){

    nr_cells        = ncells_ ;
    nr_points       = npoints_ ;

    geometry[0].setElements(nr_points) ;
    geometry[1].setElements(nr_cells) ;
    geometry[2].setElements(nr_cells) ;
    geometry[3].setElements(nconn_) ;

    for( auto &field : data ){
        if( field.getLocation() == VTKLocation::CELL)  field.setElements(nr_cells) ;
        if( field.getLocation() == VTKLocation::POINT) field.setElements(nr_points) ;
    };

    return ;
};

/*!  
 *  Reads "type" information of existing grid and calculates the correspondng connectivity size.
 *  @return     size of the connectivity information
 */
uint64_t VTKUnstructuredGrid::calcSizeConnectivity( ){

    std::fstream                  str  ;
    std::fstream::pos_type        position_appended;
    std::string                   line;
    char                     c_ ;
    uint32_t                 nbytes32 ;
    uint64_t                 nbytes64 ;
    uint64_t                 nconn ;

    str.open( fh.getName( ), std::ios::in ) ;

    //Read appended data
    //Go to the initial position of the appended section
    while( getline(str, line) && (! keywordInString( line, "<AppendedData")) ){} ;

    str >> c_;
    while( c_ != '_') str >> c_;

    position_appended = str.tellg();


    str.close();
    str.clear();

    //Open in binary for read
    str.open( fh.getName( ), std::ios::in | std::ios::binary);

    if( geometry[3].getCodification() == VTKFormat::APPENDED ){
        str.seekg( position_appended) ;
        str.seekg( geometry[3].getOffset(), std::ios::cur) ;

        if( HeaderType== "UInt32") {
            absorbBINARY( str, nbytes32 ) ;
            nconn = nbytes32 /VTKTypes::sizeOfType( geometry[3].getType() ) ;
        }

        if( HeaderType== "UInt64") {
            absorbBINARY( str, nbytes64 ) ;
            nconn = nbytes64 /VTKTypes::sizeOfType( geometry[3].getType() ) ;
        };
    };


    //Read geometry
    if(  geometry[3].getCodification() == VTKFormat::ASCII ){
        str.seekg( geometry[3].getPosition() ) ;

        std::string              line ;
        std::vector<uint64_t>    temp;

        nconn = 0 ;

        getline( str, line) ;
        while( ! keywordInString(line,"/DataArray") ) {

            temp.clear() ;
            convertString( line, temp) ;
            nconn += temp.size() ;
        };


    };

    str.close();

    return nconn ;


};

/*!  
 *  Writes entire VTU but the data.
 */
void VTKUnstructuredGrid::writeMetaData( ){

    std::fstream str ;
    std::string line ; 
    
    str.open( fh.getName( ), std::ios::out ) ;
    
    //Writing XML header
    str << "<?xml version=\"1.0\"?>" << std::endl;
    
    //Writing Piece Information
    str << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"" << HeaderType << "\">" << std::endl;
    str << "  <UnstructuredGrid>"  << std::endl;;
    str << "    <Piece  NumberOfPoints=\"" << nr_points << "\" NumberOfCells=\"" << nr_cells << "\">" << std::endl;
    
    //Header for Data
    writeDataHeader( str, false );
    
    //Wring Geometry Information
    str << "      <Points>" << std::endl ;;
    writeDataArray( str, geometry[0] ) ;
    str << "      </Points>" << std::endl;
    
    str << "      <Cells>" << std::endl ;;
    writeDataArray( str, geometry[1] ) ;
    writeDataArray( str, geometry[2] ) ;
    writeDataArray( str, geometry[3] ) ;
    str << "      </Cells>" << std::endl;
    
    //Closing Piece
    str << "    </Piece>" << std::endl;
    str << "  </UnstructuredGrid>"  << std::endl;
    
    //Appended Section
    
    str << "  <AppendedData encoding=\"raw\">" << std::endl;
    str << "_" ;
    str << std::endl ;
    str << "</VTKFile>" << std::endl;
    
    str.close() ;
    
    return ;
};

/*!  
 *  Writes collection file for parallel output. 
 *  Is called by rank 0 in VTK::Write()
 */
void VTKUnstructuredGrid::writeCollection( ){

  std::fstream str ;

  FileHandler     fhp, fho ;

  fhp = fh ;
  fho = fh ;

  fhp.setParallel(false) ;
  fhp.setAppendix("pvtu") ;

  fho.setDirectory(".") ;

  str.open( fhp.getName( ), std::ios::out ) ;

  //Writing XML header
  str << "<?xml version=\"1.0\"?>" << std::endl;

  //Writing Piece Information
  str << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  str << "  <PUnstructuredGrid GhostLevel=\"0\">"  << std::endl;;

  //Header for Data
  writeDataHeader( str, true );

  //Wring Geometry Information
  str << "      <PPoints>" << std::endl;
  writePDataArray( str, geometry[0] ) ;
  str << std::endl ;
  str << "      </PPoints>" << std::endl;


  for( int i=0; i<nr_procs; i++){
    fho.setBlock(i) ;
    str << "    <Piece  Source=\"" << fho.getName() <<  "\"/>" << std::endl;
  };

  str << "  </PUnstructuredGrid>"  << std::endl;
  str << "</VTKFile>" << std::endl;

  str.close() ;


  return ;
};

/*!  
 *  Reads meta data of VTU file (grid size, data fields, codex, position of data within file).
 *  Calls setDimension.
 */
void VTKUnstructuredGrid::readMetaData( ){

    std::fstream str;
    std::string line, temp;
    
    std::fstream::pos_type        position;
    
    str.open( fh.getName( ), std::ios::in ) ;
    
    getline( str, line);
    while( ! keywordInString( line, "<VTKFile")){
        getline(str, line);
    };
                                              
    if( getAfterKeyword( line, "header_type", '\"', temp) ){
        setHeaderType( temp) ;
    };

    while( ! keywordInString( line, "<Piece")){
      getline(str, line);
    };
    
    getAfterKeyword( line, "NumberOfPoints", '\"', temp) ;
    convertString( temp, nr_points );
    
    getAfterKeyword( line, "NumberOfCells", '\"', temp) ;
    convertString( temp, nr_cells );
    
    
    position = str.tellg() ;
    readDataHeader( str ) ;
    
    
    for( auto &field : geometry ){ 
        str.seekg( position) ;
        if( ! readDataArray( str, field ) ) {
          std::cout << field.getName() << " DataArray not found" << std::endl ;
        };
    };
    
    
    str.close() ;
    
    setDimensions( nr_cells, nr_points, calcSizeConnectivity( )  ) ;

    return ;
};

/*!  
 *  Returns the size of the connectivity information
 *  @return     size of connectivity
 */
uint64_t VTKUnstructuredGrid::getNConnectivity( ){

  return nconnectivity ;
};

/*!  
 *  sets the dimensions of the VTKUnstructuredGrid deduced from the geometry fields
 */
void VTKUnstructuredGrid::setMissingGlobalData( ){

    uint64_t    ncells, npoints, nconn;

    ncells  = std::max( geometry[1].getElements(), geometry[2].getElements() );
    npoints = geometry[0].getElements() ;
    nconn   = geometry[3].getElements() ;

    setDimensions( ncells, npoints, nconn) ;

    return ;
};

/*!
 *   @}
 */

}
