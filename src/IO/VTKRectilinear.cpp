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
 * @ingroup     VisualizationToolKit
 * @{
 * @class       VTKRectilinearGrid
 * @brief       VTK input output for Rectilinear Meshes
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
VTKRectilinearGrid::VTKRectilinearGrid( )  :VTK() {

    fh.setAppendix( "vtr" );

    geometry.push_back( VTKField( "x_Coord", VTKFieldType::SCALAR, VTKLocation::POINT, VTKDataType::Float64 ) ) ;
    geometry.push_back( VTKField( "y_Coord", VTKFieldType::SCALAR, VTKLocation::POINT, VTKDataType::Float64 ) ) ;
    geometry.push_back( VTKField( "z_Coord", VTKFieldType::SCALAR, VTKLocation::POINT, VTKDataType::Float64 ) ) ;

} ;

/*!  
 *  Constructor for parallel 3D grid.
 *  Calls default constructor and sets provided input information.
 *  @param[in]  dir_        directory of VTK file with final "/"
 *  @param[in]  name_       name of VTK without suffix
 */
VTKRectilinearGrid::VTKRectilinearGrid( std::string dir_, std::string name_ ) :VTKRectilinearGrid( ) {

    setNames( dir_, name_ ) ;

    return ;
};

/*!  
 *  Constructor for parallel 3D grid.
 *  Calls default constructor and sets provided input information.
 *  @param[in]  dir_        directory of VTK file with final "/"
 *  @param[in]  name_       name of VTK without suffix
 *  @param[in]  codex_      codex of data [VTKFormat::ASCII/VTKFormat::APPENDED]
 *  @param[in]  n1_         min node index in first direction 
 *  @param[in]  n2_         max node index in first direction 
 *  @param[in]  m1_         min node index in second direction 
 *  @param[in]  m2_         max node index in second direction 
 *  @param[in]  l1_         min node index in third direction 
 *  @param[in]  l2_         max node index in third direction 
 */
VTKRectilinearGrid::VTKRectilinearGrid( std::string dir_, std::string name_, VTKFormat codex_, int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ ) :VTKRectilinearGrid( ) {

    setNames( dir_, name_ ) ;

    setDimensions( n1_, n2_, m1_, m2_, l1_, l2_) ;
    setGeomCodex( codex_ ) ;

    return ;
};

/*!  
 *  Constructor for serial 3D grids.
 *  Min index is set automatically to zero.
 *  Calls default constructor and sets provided input information.
 *  @param[in]  dir_        directory of VTK file with final "/"
 *  @param[in]  name_       name of VTK without suffix
 *  @param[in]  codex_      codex of data [VTKFormat::ASCII/VTKFormat::APPENDED]
 *  @param[in]  n_          number of nodes in first direction 
 *  @param[in]  m_          number of nodes in second direction 
 *  @param[in]  l_          number of nodes in third direction 
 */
VTKRectilinearGrid::VTKRectilinearGrid( std::string dir_, std::string name_, VTKFormat codex_, int n_, int m_, int l_ ) : VTKRectilinearGrid( ) {

    setNames( dir_, name_ ) ;

    setDimensions( 0, n_-1, 0, m_-1, 0, l_-1) ;
    setGeomCodex( codex_ ) ;

    return ;
};

/*!  
 *  Constructor for parallel 2D grid.
 *  Calls default constructor and sets provided input information.
 *  @param[in]  dir_        directory of VTK file with final "/"
 *  @param[in]  name_       name of VTK without suffix
 *  @param[in]  codex_      codex of data [VTKFormat::ASCII/VTKFormat::APPENDED]
 *  @param[in]  n1_         min node index in first direction 
 *  @param[in]  n2_         max node index in first direction 
 *  @param[in]  m1_         min node index in second direction 
 *  @param[in]  m2_         max node index in second direction 
 */
VTKRectilinearGrid::VTKRectilinearGrid( std::string dir_, std::string name_, VTKFormat codex_, int n1_, int n2_, int m1_, int m2_ ) :VTKRectilinearGrid( ) {

    setNames( dir_, name_ ) ;

    setDimensions( n1_, n2_, m1_, m2_, 0, 0) ;
    setGeomCodex( codex_ ) ;

    return ;
};

/*!  
 *  Constructor for serial 2D grid.
 *  Calls default constructor and sets provided input information.
 *  @param[in]  dir_        directory of VTK file with final "/"
 *  @param[in]  name_       name of VTK without suffix
 *  @param[in]  codex_      codex of data [VTKFormat::ASCII/VTKFormat::APPENDED]
 *  @param[in]  n_          number of nodes in first direction 
 *  @param[in]  m_          number of nodes in second direction 
 */
VTKRectilinearGrid::VTKRectilinearGrid( std::string dir_, std::string name_, VTKFormat codex_, int n_, int m_ ) : VTKRectilinearGrid( ) {

    setNames( dir_, name_ ) ;

    setDimensions( 0, n_-1, 0, m_-1, 0, 0) ;
    setGeomCodex( codex_ ) ;

    return ;
};

/*!  
 *  Destructor 
 */
VTKRectilinearGrid::~VTKRectilinearGrid( ){
} ;

/*!  
 *  sets the type of the geometry variables
 *  @param[in]  Ptype       Type of "x_Coord", "y_Coord" and "z_Coord" geometry information [ VTKDataType::Float[32/64] ]
 */
void VTKRectilinearGrid::setGeomTypes( VTKDataType Ptype ){

    geometry[0].setType(Ptype) ;
    geometry[1].setType(Ptype) ;
    geometry[2].setType(Ptype) ;

    return ;
};

/*!  
 *  Reads meta data of VTR file (grid size, data fields, codex, position of data within file).
 *  Calls setDimension.
 */
void VTKRectilinearGrid::readMetaData( ){

    std::fstream str;
    std::string line, temp;

    std::fstream::pos_type        position;
    std::array<int,6>             extensions ;


    str.open( fh.getName( ), std::ios::in ) ;

    getline( str, line);
    while( ! bitpit::utils::keywordInString( line, "<VTKFile")){
        getline(str, line);
    };

    if( bitpit::utils::getAfterKeyword( line, "header_type", '\"', temp) ){
        setHeaderType( temp) ;
    };

    while( ! bitpit::utils::keywordInString( line, "<Piece")){
        getline(str, line);
    };

    bitpit::utils::getAfterKeyword( line, "Extent", '\"', temp) ;
    bitpit::utils::convertString( temp, extensions );

    local_index[0][0] = extensions[0] ;
    local_index[0][1] = extensions[1] ;
    local_index[1][0] = extensions[2] ;
    local_index[1][1] = extensions[3] ;
    local_index[2][0] = extensions[4] ;
    local_index[2][1] = extensions[5] ;

    position = str.tellg() ;

    readDataHeader( str ) ;

    for( auto &field : geometry ){ //int i=0; i<geometry.size(); ++i){
        str.seekg( position) ;
        if( ! readDataArray( str, field ) ) {
            std::cout << field.getName() << " DataArray not found" << std::endl ;
        };
    };


    setDimensions( local_index[0][0], local_index[0][1], local_index[1][0], local_index[1][1], local_index[2][0], local_index[2][1] ) ;
    str.close() ; 

    return ;

};

/*!  
 *  Writes entire VTR but the data.
 */
void VTKRectilinearGrid::writeMetaData( ){

    std::fstream str;

    str.open( fh.getName( ), std::ios::out ) ;

    //Writing XML header
    str << "<?xml version=\"1.0\"?>" << std::endl;

    //Writing Piece Information
    str << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\"  header_type=\"" << HeaderType << "\">" << std::endl; 
    str << "  <RectilinearGrid WholeExtent= \"" 
        << global_index[0][0] << " " << global_index[0][1]<< " "
        << global_index[1][0] << " " << global_index[1][1]<< " "
        << global_index[2][0] << " " << global_index[2][1]<< " "
        << "\" >" << std::endl;

    str << "    <Piece Extent= \" " 
        << local_index[0][0] << " " << local_index[0][1]<< " "
        << local_index[1][0] << " " << local_index[1][1]<< " "
        << local_index[2][0] << " " << local_index[2][1]<< " "
        << "\" >" << std::endl;



    //Header for Data
    writeDataHeader( str, false ) ;

    //Wring Geometry Information   
    str << "       <Coordinates>" << std::endl;
    writeDataArray( str, geometry[0] ) ;
    writeDataArray( str, geometry[1] ) ;
    writeDataArray( str, geometry[2] ) ;
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

    return ;
};

/*!  
 *  Writes collection file for parallel output. 
 *  Is called by rank 0 in VTK::Write()
 */
void VTKRectilinearGrid::writeCollection( ){

    std::fstream str ;

    FileHandler     fhp, fho ;

    fhp = fh ;
    fho = fh ;

    fhp.setParallel(false) ;
    fhp.setAppendix("pvtr") ;

    fho.setDirectory(".") ;

    str.open( fhp.getName( ), std::ios::out ) ;

    //Writing XML header
    str << "<?xml version=\"1.0\"?>" << std::endl;

    //Writing Piece Information
    str << "<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
    str << "  <PRectilinearGrid WholeExtent= \"" 
        << global_index[0][0] << " " << global_index[0][1]<< " "
        << global_index[1][0] << " " << global_index[1][1]<< " "
        << global_index[2][0] << " " << global_index[2][1]<< " "
        << "GhostLevel=\"0\">" << std::endl;



    //Header for Data
    writeDataHeader( str, true );

    //Wring Geometry Information
    str << "      <PCoordinates>" << std::endl;
    writePDataArray( str, geometry[0] ) ;
    writePDataArray( str, geometry[1] ) ;
    writePDataArray( str, geometry[2] ) ;
    str << "      </PCoordinates>" << std::endl;


    for( int i=0; i<nr_procs; i++){
        fho.setBlock(i) ;
        extension3D_t &index = proc_index[i] ;

        str << "    <Piece Extent= \" " 
            << index[0][0] << " " << index[0][1] << " "
            << index[1][0] << " " << index[1][1] << " "
            << index[2][0] << " " << index[2][1] << " "
            << "\" Source= \"" << fho.getName() << "\"/>" << std::endl;
    };

    str << "  </PRectilinearGrid>"  << std::endl;
    str << "</VTKFile>" << std::endl;

    str.close() ;


    return ;
};

/*!  
 *  sets the dimension for 3D parallel grids.
 *  @param[in]  n1_         min node index in first direction 
 *  @param[in]  n2_         max node index in first direction 
 *  @param[in]  m1_         min node index in second direction 
 *  @param[in]  m2_         max node index in second direction 
 *  @param[in]  l1_         min node index in third direction 
 *  @param[in]  l2_         max node index in third direction 
 */
void VTKRectilinearGrid::setDimensions( int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ ){

    dimensions = ( l1_ == l2_ ) ? 2 :3 ;

    local_index[0][0] = n1_ ;
    local_index[0][1] = n2_ ;
    local_index[1][0] = m1_ ;
    local_index[1][1] = m2_ ;
    local_index[2][0] = l1_ ;
    local_index[2][1] = l2_ ;

    if(nr_procs==1)
        global_index = local_index ;

    for(int d=0; d<3; ++d){
        geometry[d].setElements( local_index[d][1] -local_index[d][0] +1 ) ;
    };


    nr_cells = 1 ;
    nr_points = 1 ;

    for(int d=0; d<dimensions; ++d){
        nr_cells = nr_cells *  ( local_index[d][1] -local_index[d][0] ) ;
        nr_points = nr_points *  ( local_index[d][1] -local_index[d][0] +1 ) ;
    };

    for( auto &field : data ){
        if( field.getLocation() == VTKLocation::CELL)  field.setElements(nr_cells) ;
        if( field.getLocation() == VTKLocation::POINT) field.setElements(nr_points) ;
    };

    return ;
};

/*!  
 *  sets the dimension for 3D serial grids.
 *  @param[in]  n_          number of nodes in first direction 
 *  @param[in]  m_          number of nodes in second direction 
 *  @param[in]  l_          number of nodes in third direction 
 */
void VTKRectilinearGrid::setDimensions( int n_, int m_, int l_ ){

    this->setDimensions( 0, n_-1, 0, m_-1, 0, l_-1 );

    return ;
};

/*!  
 *  sets the dimension for 2D parallel grids.
 *  @param[in]  n1_         min node index in first direction 
 *  @param[in]  n2_         max node index in first direction 
 *  @param[in]  m1_         min node index in second direction 
 *  @param[in]  m2_         max node index in second direction 
 */
void VTKRectilinearGrid::setDimensions( int n1_, int n2_, int m1_, int m2_ ){

    this->setDimensions( n1_, n2_, m1_, m2_, 0, 0 );

    return ;
};

/*!  
 *  sets the dimension for 2D serial grids.
 *  @param[in]  n_          number of nodes in first direction 
 *  @param[in]  m_          number of nodes in second direction 
 */
void VTKRectilinearGrid::setDimensions( int n_, int m_ ){

    this->setDimensions( 0, n_-1, 0, m_-1, 0, 0 );

    return ;
};

/*!  
 *  sets the global 3D grid information for parallel output.
 *  Needs to be called by all processes
 *  @param[in]  I           numer of nodes in first direction
 *  @param[in]  J           numer of nodes in second direction
 *  @param[in]  K           numer of nodes in third direction
 */
void VTKRectilinearGrid::setGlobalDimensions( int I, int J, int K ){


    global_index[0][0] = 0;
    global_index[0][1] = I;
    global_index[1][0] = 0;
    global_index[1][1] = J;
    global_index[2][0] = 0;
    global_index[2][1] = K;

    return;
};

/*!  
 *  sets the global 2D grid information for parallel output.
 *  Needs to be called by all processes
 *  @param[in]  I           numer of nodes in first direction
 *  @param[in]  J           numer of nodes in second direction
 */
void VTKRectilinearGrid::setGlobalDimensions( int I, int J ){


    global_index[0][0] = 0;
    global_index[0][1] = I;
    global_index[1][0] = 0;
    global_index[1][1] = J;
    global_index[2][0] = 0;
    global_index[2][1] = 0;

    return;
};

/*!  
 *  sets the global 3D grid information for parallel output.
 *  Needs to be called only by rank 0.
 *  @param[in]  loc_        min I-index, max I-index, min J-index, max J-index, min K-index, max K-index for each process
 */
void VTKRectilinearGrid::setGlobalIndex( std::vector<extension3D_t> loc_ ){

    if( loc_.size() != nr_procs ) 
        std::cout << "Size of loc_ in VTKRectilinearGrid::setParallelIndex does not fit nr_procs " << std::endl ;

    proc_index   = loc_ ;

    return;
};

/*!  
 *  sets the global 2D grid information for parallel output.
 *  Needs to be called only by rank 0.
 *  @param[in]  loc_        min I-index, max I-index, min J-index, max J-index for each process
 */
void VTKRectilinearGrid::setGlobalIndex( std::vector<extension2D_t> loc_ ){

    if( loc_.size() !=nr_procs ) std::cout << "Size of loc_ in VTKRectilinearGrid::setParallelIndex does not fit nr_procs " << std::endl ;

    proc_index.resize(nr_procs) ;

    for( int i=0; i< nr_procs; ++i){
        proc_index[i][0]   = loc_[i][0] ;
        proc_index[i][1]   = loc_[i][1] ;
        proc_index[i][2]   = {{0,0}} ;
    }

    return;
};

/*!  
 *  sets the dimensions of the VTKRectilinerGrid deduced from the geometry fields
 */
void VTKRectilinearGrid::setMissingGlobalData( ){

    int    nX, nY, nZ;

    nX = geometry[0].getElements() ;
    nY = geometry[1].getElements() ;
    nZ = geometry[2].getElements() ;

    setDimensions( nX, nY, nZ ) ;

    return ;
};
/*!
 *   @}
 */

}
