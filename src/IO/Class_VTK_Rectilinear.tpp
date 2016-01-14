/*! ========================================================================================
 * @ingroup     VisualizationToolKit
 * @{
 * @class       VTK_RectilinearGrid
 * @brief       VTK input output for Rectilinear Meshes
 * @tparam      Derived     this argument is used for the static-dispatch interface through CRTP
 *
 * VTK_RectilinearGrid provides methods to read and write parallel and serial rectlinear meshes and data. 
 * The class is agnostic with respect to the container used for the data and provides an interface through the CRTP mechanism.
 * The numbering of nodes start with 0. Different numbering scheme is not supported.
 *
 */

/*!  ===================================================================================
 *  Default constructor.
 *  Allocates three geometry fields called "x_Coord", "y_Coord" and "z_Coord".
 *  @tparam     Derived     Derived class for CRTP
 */
template <class Derived>
VTK_RectilinearGrid<Derived>::VTK_RectilinearGrid( )  :VTK() {

    fh.SetAppendix( "vtr" );

    geometry.push_back( VTK::Field_C( "x_Coord", 1, "Float64", "Point" ) ) ;
    geometry.push_back( VTK::Field_C( "y_Coord", 1, "Float64", "Point" ) ) ;
    geometry.push_back( VTK::Field_C( "z_Coord", 1, "Float64", "Point" ) ) ;

} ;

/*!  ===================================================================================
 *  Constructor for parallel 3D grid.
 *  Calls default constructor and sets provided input information.
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  dir_        directory of VTK file with final "/"
 *  @param[in]  name_       name of VTK without suffix
 */
template <class Derived>
VTK_RectilinearGrid<Derived>::VTK_RectilinearGrid( std::string dir_, std::string name_ ) :VTK_RectilinearGrid( ) {

    SetNames( dir_, name_ ) ;

    return ;
};

/*!  ===================================================================================
 *  Constructor for parallel 3D grid.
 *  Calls default constructor and sets provided input information.
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  dir_        directory of VTK file with final "/"
 *  @param[in]  name_       name of VTK without suffix
 *  @param[in]  codex_      codex of data ["ascii"/"appended"]
 *  @param[in]  n1_         min node index in first direction 
 *  @param[in]  n2_         max node index in first direction 
 *  @param[in]  m1_         min node index in second direction 
 *  @param[in]  m2_         max node index in second direction 
 *  @param[in]  l1_         min node index in third direction 
 *  @param[in]  l2_         max node index in third direction 
 */
template <class Derived>
VTK_RectilinearGrid<Derived>::VTK_RectilinearGrid( std::string dir_, std::string name_, std::string codex_, int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ ) :VTK_RectilinearGrid( ) {

    SetNames( dir_, name_ ) ;

    SetDimensions( n1_, n2_, m1_, m2_, l1_, l2_) ;
    SetGeomCodex( codex_ ) ;

    return ;
};

/*!  ===================================================================================
 *  Constructor for serial 3D grids.
 *  Min index is set automatically to zero.
 *  Calls default constructor and sets provided input information.
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  dir_        directory of VTK file with final "/"
 *  @param[in]  name_       name of VTK without suffix
 *  @param[in]  codex_      codex of data ["ascii"/"appended"]
 *  @param[in]  n_          number of nodes in first direction 
 *  @param[in]  m_          number of nodes in second direction 
 *  @param[in]  l_          number of nodes in third direction 
 */
template <class Derived>
VTK_RectilinearGrid<Derived>::VTK_RectilinearGrid( std::string dir_, std::string name_, std::string codex_, int n_, int m_, int l_ ) : VTK_RectilinearGrid( ) {

    SetNames( dir_, name_ ) ;

    SetDimensions( 0, n_-1, 0, m_-1, 0, l_-1) ;
    SetGeomCodex( codex_ ) ;

    return ;
};

/*!  ===================================================================================
 *  Constructor for parallel 2D grid.
 *  Calls default constructor and sets provided input information.
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  dir_        directory of VTK file with final "/"
 *  @param[in]  name_       name of VTK without suffix
 *  @param[in]  codex_      codex of data ["ascii"/"appended"]
 *  @param[in]  n1_         min node index in first direction 
 *  @param[in]  n2_         max node index in first direction 
 *  @param[in]  m1_         min node index in second direction 
 *  @param[in]  m2_         max node index in second direction 
 */
template <class Derived>
VTK_RectilinearGrid<Derived>::VTK_RectilinearGrid( std::string dir_, std::string name_, std::string codex_, int n1_, int n2_, int m1_, int m2_ ) :VTK_RectilinearGrid( ) {

    SetNames( dir_, name_ ) ;

    SetDimensions( n1_, n2_, m1_, m2_, 0, 0) ;
    SetGeomCodex( codex_ ) ;

    return ;
};

/*!  ===================================================================================
 *  Constructor for serial 2D grid.
 *  Calls default constructor and sets provided input information.
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  dir_        directory of VTK file with final "/"
 *  @param[in]  name_       name of VTK without suffix
 *  @param[in]  codex_      codex of data ["ascii"/"appended"]
 *  @param[in]  n_          number of nodes in first direction 
 *  @param[in]  m_          number of nodes in second direction 
 */
template <class Derived>
VTK_RectilinearGrid<Derived>::VTK_RectilinearGrid( std::string dir_, std::string name_, std::string codex_, int n_, int m_ ) : VTK_RectilinearGrid( ) {

    SetNames( dir_, name_ ) ;

    SetDimensions( 0, n_-1, 0, m_-1, 0, 0) ;
    SetGeomCodex( codex_ ) ;

    return ;
};


/*!  ===================================================================================
 *  Destructor 
 *  @tparam     Derived     Derived class for CRTP
 */
template <class Derived>
VTK_RectilinearGrid<Derived>::~VTK_RectilinearGrid( ) {} ;

/*!  ===================================================================================
 *  Sets the type of the geometry variables
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  Ptype       Type of "x_Coord", "y_Coord" and "z_Coord" geometry information [ "Float[32/64]"]
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::SetGeomTypes( std::string Ptype ){

    geometry[0].SetType(Ptype) ;
    geometry[1].SetType(Ptype) ;
    geometry[2].SetType(Ptype) ;

    return ;
};

/*!  ===================================================================================
 *  Reads meta data of VTR file (grid size, data fields, codex, position of data within file).
 *  Calls SetDimension.
 *  @tparam     Derived     Derived class for CRTP
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::ReadMetaData( ){

    std::fstream str;
    std::string line, temp;

    std::fstream::pos_type        position;

    bool                     read ;

    std::array<int,6>             extensions ;


    str.open( fh.GetName( ), ios::in ) ;

    getline( str, line);
    while( ! Keyword_In_String( line, "<VTKFile")){
        getline(str, line);
    };

    if( Get_After_Keyword( line, "header_type", '\"', temp) ){
        SetHeaderType( temp) ;
    };

    while( ! Keyword_In_String( line, "<Piece")){
        getline(str, line);
    };

    Get_After_Keyword( line, "Extent", '\"', temp) ;
    convert_string( temp, extensions );

    local_index[0][0] = extensions[0] ;
    local_index[0][1] = extensions[1] ;
    local_index[1][0] = extensions[2] ;
    local_index[1][1] = extensions[3] ;
    local_index[2][0] = extensions[4] ;
    local_index[2][1] = extensions[5] ;

    position = str.tellg() ;

    ReadDataHeader( str ) ;

    for( int i=0; i<geometry.size(); ++i){
        str.seekg( position) ;
        if( ! ReadDataArray( str, geometry[i] ) ) {
            cout << geometry[i].GetName() << " DataArray not found" << endl ;
        };
    };


    SetDimensions( local_index[0][0], local_index[0][1], local_index[1][0], local_index[1][1], local_index[2][0], local_index[2][1] ) ;
    str.close() ; 

    return ;

};

/*!  ===================================================================================
 *  Writes entire VTR but the data.
 *  @tparam     Derived     Derived class for CRTP
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::WriteMetaData( ){

    std::fstream str;

    str.open( fh.GetName( ), ios::out ) ;

    //Writing XML header
    str << "<?xml version=\"1.0\"?>" << endl;

    //Writing Piece Information
    str << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\"  header_type=\"" << HeaderType << "\">" << endl; 
    str << "  <RectilinearGrid WholeExtent= \"" 
        << global_index[0][0] << " " << global_index[0][1]<< " "
        << global_index[1][0] << " " << global_index[1][1]<< " "
        << global_index[2][0] << " " << global_index[2][1]<< " "
        << "\" >" << endl;

    str << "    <Piece Extent= \" " 
        << local_index[0][0] << " " << local_index[0][1]<< " "
        << local_index[1][0] << " " << local_index[1][1]<< " "
        << local_index[2][0] << " " << local_index[2][1]<< " "
        << "\" >" << endl;



    //Header for Data
    WriteDataHeader( str, false ) ;

    //Wring Geometry Information   
    str << "       <Coordinates>" << endl;
    WriteDataArray( str, geometry[0] ) ;
    WriteDataArray( str, geometry[1] ) ;
    WriteDataArray( str, geometry[2] ) ;
    str << "       </Coordinates>" << endl;

    //Closing Piece
    str << "    </Piece>" << endl;
    str << "  </RectilinearGrid>" << endl;

    //Write Appended Section
    str << "  <AppendedData encoding=\"raw\">" << endl;
    str << "_" ;
    str << endl ;

    //Closing XML
    str << "</VTKFile>" << endl;

    str.close() ;

    return ;
};

/*!  ===================================================================================
 *  Writes collection file for parallel output. 
 *  Is called by rank 0 in VTK::Write()
 *  @tparam     Derived     Derived class for CRTP
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::WriteCollection( ){

    std::fstream str ;

    FileHandler_C   fhp, fho ;

    fhp = fh ;
    fho = fh ;

    fhp.SetParallel(false) ;
    fhp.SetAppendix("pvtr") ;

    fho.SetDirectory(".") ;

    str.open( fhp.GetName( ), ios::out ) ;

    //Writing XML header
    str << "<?xml version=\"1.0\"?>" << endl;

    //Writing Piece Information
    str << "<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    str << "  <PRectilinearGrid WholeExtent= \"" 
        << global_index[0][0] << " " << global_index[0][1]<< " "
        << global_index[1][0] << " " << global_index[1][1]<< " "
        << global_index[2][0] << " " << global_index[2][1]<< " "
        << "GhostLevel=\"0\">" << endl;



    //Header for Data
    WriteDataHeader( str, true );

    //Wring Geometry Information
    str << "      <PCoordinates>" << endl;
    WritePDataArray( str, geometry[0] ) ;
    WritePDataArray( str, geometry[1] ) ;
    WritePDataArray( str, geometry[2] ) ;
    str << "      </PCoordinates>" << endl;


    for( int i=0; i<nr_procs; i++){
        fho.SetBlock(i) ;
        extension3D_t &index = proc_index[i] ;

        str << "    <Piece Extent= \" " 
            << index[0][0] << " " << index[0][1] << " "
            << index[1][0] << " " << index[1][1] << " "
            << index[2][0] << " " << index[2][1] << " "
            << "\" Source= \"" << fho.GetName() << "\"/>" << endl;
    };

    str << "  </PRectilinearGrid>"  << endl;
    str << "</VTKFile>" << endl;

    str.close() ;


    return ;
};

/*!  ===================================================================================
 *  Sets the dimension for 3D parallel grids.
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  n1_         min node index in first direction 
 *  @param[in]  n2_         max node index in first direction 
 *  @param[in]  m1_         min node index in second direction 
 *  @param[in]  m2_         max node index in second direction 
 *  @param[in]  l1_         min node index in third direction 
 *  @param[in]  l2_         max node index in third direction 
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::SetDimensions( int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ ){

    dimensions = ( l1_ == l2_ ) ? 2 :3 ;

    local_index[0][0] = n1_ ;
    local_index[0][1] = n2_ ;
    local_index[1][0] = m1_ ;
    local_index[1][1] = m2_ ;
    local_index[2][0] = l1_ ;
    local_index[2][1] = l2_ ;

    if(nr_procs=1)
        global_index = local_index ;

    for(int d=0; d<3; ++d){
        geometry[d].SetElements( local_index[d][1] -local_index[d][0] +1 ) ;
    };


    nr_cells = 1 ;
    nr_points = 1 ;

    for(int d=0; d<dimensions; ++d){
        nr_cells = nr_cells *  ( local_index[d][1] -local_index[d][0] ) ;
        nr_points = nr_points *  ( local_index[d][1] -local_index[d][0] +1 ) ;
    };

    for( int i=0; i< nr_data; i++){
        if( data[i].GetLocation() == "Cell")  data[i].SetElements(nr_cells) ;
        if( data[i].GetLocation() == "Point") data[i].SetElements(nr_points) ;
    };

    return ;
};

/*!  ===================================================================================
 *  Sets the dimension for 3D serial grids.
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  n_          number of nodes in first direction 
 *  @param[in]  m_          number of nodes in second direction 
 *  @param[in]  l_          number of nodes in third direction 
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::SetDimensions( int n_, int m_, int l_ ){

    this->SetDimensions( 0, n_-1, 0, m_-1, 0, l_-1 );

    return ;
};

/*!  ===================================================================================
 *  Sets the dimension for 2D parallel grids.
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  n1_         min node index in first direction 
 *  @param[in]  n2_         max node index in first direction 
 *  @param[in]  m1_         min node index in second direction 
 *  @param[in]  m2_         max node index in second direction 
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::SetDimensions( int n1_, int n2_, int m1_, int m2_ ){

    this->SetDimensions( n1_, n2_, m1_, m2_, 0, 0 );

    return ;
};

/*!  ===================================================================================
 *  Sets the dimension for 2D serial grids.
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  n_          number of nodes in first direction 
 *  @param[in]  m_          number of nodes in second direction 
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::SetDimensions( int n_, int m_ ){

    this->SetDimensions( 0, n_-1, 0, m_-1, 0, 0 );

    return ;
};

/*!  ===================================================================================
 *  Sets the global 3D grid information for parallel output.
 *  Needs to be called by all processes
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  I           numer of nodes in first direction
 *  @param[in]  J           numer of nodes in second direction
 *  @param[in]  K           numer of nodes in third direction
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::SetGlobalDimensions( int I, int J, int K ){


    global_index[0][0] = 0;
    global_index[0][1] = I;
    global_index[1][0] = 0;
    global_index[1][1] = J;
    global_index[2][0] = 0;
    global_index[2][1] = K;

    return;
};

/*!  ===================================================================================
 *  Sets the global 2D grid information for parallel output.
 *  Needs to be called by all processes
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  I           numer of nodes in first direction
 *  @param[in]  J           numer of nodes in second direction
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::SetGlobalDimensions( int I, int J ){


    global_index[0][0] = 0;
    global_index[0][1] = I;
    global_index[1][0] = 0;
    global_index[1][1] = J;
    global_index[2][0] = 0;
    global_index[2][1] = 0;

    return;
};

/*!  ===================================================================================
 *  Sets the global 3D grid information for parallel output.
 *  Needs to be called only by rank 0.
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  loc_        min I-index, max I-index, min J-index, max J-index, min K-index, max K-index for each process
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::SetGlobalIndex( std::vector<extension3D_t> loc_ ){

    if( loc_.size() !=nr_procs ) cout << "Size of loc_ in VTK_RectilinearGrid<Derived>::SetParallelIndex does not fit nr_procs " << endl ;

    proc_index   = loc_ ;

    return;
};

/*!  ===================================================================================
 *  Sets the global 2D grid information for parallel output.
 *  Needs to be called only by rank 0.
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  loc_        min I-index, max I-index, min J-index, max J-index for each process
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::SetGlobalIndex( std::vector<extension2D_t> loc_ ){

    if( loc_.size() !=nr_procs ) cout << "Size of loc_ in VTK_RectilinearGrid<Derived>::SetParallelIndex does not fit nr_procs " << endl ;

    proc_index.resize(nr_procs) ;

    for( int i=0; i< nr_procs; ++i){
        proc_index[i][0]   = loc_[i][0] ;
        proc_index[i][1]   = loc_[i][1] ;
        proc_index[i][2]   = {0,0} ;
    }

    return;
};

/*!  ===================================================================================
 *  CRPT interface for writing data to stream.
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  str         stream to write to
 *  @param[in]  codex       codex which must be used ["ascii"/"appended"]. If "appended" a unformatted binary stream must be used
 *  @param[in]  name        name of the data to be written. Either user data or grid data
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::Flush( std::fstream &str, std::string codex, std::string name ){

    static_cast<Derived *>(this)->Flush( str, codex, name );
    return ;

};

/*!  ===================================================================================
 *  CRPT interface for reading data from stream.
 *  @tparam     Derived     Derived class for CRTP
 *  @param[in]  str         stream to read from
 *  @param[in]  codex       codex which must be used ["ascii"/"appended"]. If "appended" a unformatted binary stream must be used
 *  @param[in]  name        name of the data to be read. Either user data or grid data
 */
template <class Derived>
void VTK_RectilinearGrid<Derived>::Absorb( std::fstream &str, std::string codex, std::string name ){

    static_cast<Derived *>(this)->Absorb( str, codex, name );
    return ;

};

/*!
 *   @}
 */
