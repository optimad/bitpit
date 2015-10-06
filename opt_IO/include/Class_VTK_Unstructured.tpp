//------------------------------------------------------------------
template <class Derived>
VTK_UnstructuredGrid<Derived>::VTK_UnstructuredGrid( )
                              :VTK() {

  fh.SetAppendix("vtu");

  geometry.push_back( VTK::Field_C( "Points",       3,  "Point" ) ) ;
  geometry.push_back( VTK::Field_C( "offsets",      1,  "Cell"  ) ) ;
  geometry.push_back( VTK::Field_C( "types",        1,  "Cell"  ) ) ;
  geometry.push_back( VTK::Field_C( "connectivity", 1,  "Cell"  ) ) ;

};
//------------------------------------------------------------------
template <class Derived>
VTK_UnstructuredGrid<Derived>::VTK_UnstructuredGrid( string dir_, string name_  )
                              :VTK_UnstructuredGrid( ){

  SetNames( dir_, name_ ) ; 
  return ;

};

//------------------------------------------------------------------
template <class Derived>
VTK_UnstructuredGrid<Derived>::~VTK_UnstructuredGrid( ) {
    
};

//------------------------------------------------------------------
template <class Derived>
void VTK_UnstructuredGrid<Derived>::SetGeomTypes(string Ptype, string Otype, string Ttype, string Ctype  ){

    geometry[0].SetType(Ptype) ;
    geometry[1].SetType(Otype) ;
    geometry[2].SetType(Ttype) ;
    geometry[3].SetType(Ctype) ;

    return ;
};

//------------------------------------------------------------------
template <class Derived>
void VTK_UnstructuredGrid<Derived>::SetDimensions( uint64_t ncells_, uint64_t npoints_, uint64_t nconn_ ){

    nr_cells        = ncells_ ;
    nr_points       = npoints_ ;

    geometry[0].SetElements(nr_points) ;
    geometry[1].SetElements(nr_cells) ;
    geometry[2].SetElements(nr_cells) ;
    geometry[3].SetElements(nconn_) ;

    for( int i=0; i< nr_data; i++){
        if( data[i].GetLocation() == "Cell")  data[i].SetElements(nr_cells) ;
        if( data[i].GetLocation() == "Point") data[i].SetElements(nr_points) ;
    };

    return ;
};

// =================================================================================== //
template <class Derived>
uint64_t VTK_UnstructuredGrid<Derived>::CalcSizeConnectivity( ){

    fstream                  str  ;
    fstream::pos_type        position_appended;
    string                   line;
    char                     c_ ;
    uint32_t                 nbytes32 ;
    uint64_t                 nbytes64 ;
    uint64_t                 nconn ;

    str.open( fh.GetName( ), ios::in ) ;

    //Read appended data
    //Go to the initial position of the appended section
    while( getline(str, line) && (! Keyword_In_String( line, "<AppendedData")) ){} ;

    str >> c_;
    while( c_ != '_') str >> c_;

    position_appended = str.tellg();


    str.close();
    str.clear();

    //Open in binary for read
    str.open( fh.GetName( ), ios::in | ios::binary);

    if( geometry[3].GetCodification() == "appended" ){
        str.seekg( position_appended) ;
        str.seekg( geometry[3].GetOffset(), ios::cur) ;

        if( HeaderType== "UInt32") {
            absorb_binary( str, nbytes32 ) ;
            nconn = nbytes32 /SizeOfType( geometry[3].GetType() ) ;
        }

        if( HeaderType== "UInt64") {
            absorb_binary( str, nbytes64 ) ;
            nconn = nbytes64 /SizeOfType( geometry[3].GetType() ) ;
        };
    };


    //Read geometry
    if(  geometry[3].GetCodification() == "ascii"){
        str.seekg( geometry[3].GetPosition() ) ;

        string              line ;
        vector<uint64_t>    temp;

        nconn = 0 ;

        getline( str, line) ;
        while( ! Keyword_In_String(line,"/DataArray") ) {

            temp.clear() ;
            convert_string( line, temp) ;
            nconn += temp.size() ;
        };


    };

    str.close();

    return nconn ;


};

// =================================================================================== //
template <class Derived>
uint8_t VTK_UnstructuredGrid<Derived>::NumberOfElements( uint8_t t){

  int e;

  switch (t){

    case 1: e= 1; break;
    case 3: e= 2; break;
    case 5: case 21: e =3; break;
    case 8: case 9:  case 10: e= 4; break;
    case 11: case 12: case 23: case 24: e=8; break;
    case 13: case 22: e= 6; break;
    case 14: e=5; break;
    case 25: e=20; break;
    default:
      cout << "Element type not supported: " << t << endl;

  };

  return e;

};

// =================================================================================== //
template <class Derived>
void VTK_UnstructuredGrid<Derived>::WriteMetaData( ){

    fstream str ;
    string line ; 
    
    str.open( fh.GetName( ), ios::out ) ;
    
    //Writing XML header
    str << "<?xml version=\"1.0\"?>" << endl;
    
    //Writing Piece Information
    str << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"" << HeaderType << "\">" << endl;
    str << "  <UnstructuredGrid>"  << endl;;
    str << "    <Piece  NumberOfPoints=\"" << nr_points << "\" NumberOfCells=\"" << nr_cells << "\">" << endl;
    
    //Header for Data
    WriteDataHeader( str, false );
    
    //Wring Geometry Information
    str << "      <Points>" << endl ;;
    WriteDataArray( str, geometry[0] ) ;
    str << "      </Points>" << endl;
    
    str << "      <Cells>" << endl ;;
    WriteDataArray( str, geometry[1] ) ;
    WriteDataArray( str, geometry[2] ) ;
    WriteDataArray( str, geometry[3] ) ;
    str << "      </Cells>" << endl;
    
    //Closing Piece
    str << "    </Piece>" << endl;
    str << "  </UnstructuredGrid>"  << endl;
    
    //Appended Section
    
    str << "  <AppendedData encoding=\"raw\">" << endl;
    str << "_" ;
    str << endl ;
    str << "</VTKFile>" << endl;
    
    str.close() ;
    
    return ;
};

// =================================================================================== //
template <class Derived>
void VTK_UnstructuredGrid<Derived>::WriteCollection( ){

  fstream str ;

  FileHandler_C   fhp, fho ;

  fhp = fh ;
  fho = fh ;

  fhp.SetParallel(false) ;
  fhp.SetAppendix("pvtu") ;

  fho.SetDirectory(".") ;

  str.open( fhp.GetName( ), ios::out ) ;

  //Writing XML header
  str << "<?xml version=\"1.0\"?>" << endl;

  //Writing Piece Information
  str << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
  str << "  <PUnstructuredGrid GhostLevel=\"0\">"  << endl;;

  //Header for Data
  WriteDataHeader( str, true );

  //Wring Geometry Information
  str << "      <PPoints>" << endl;
  WritePDataArray( str, geometry[0] ) ;
  str << endl ;
  str << "      </PPoints>" << endl;


  for( int i=0; i<nr_procs; i++){
    fho.SetBlock(i) ;
    str << "    <Piece  Source=\"" << fho.GetName() <<  "\"/>" << endl;
  };

  str << "  </PUnstructuredGrid>"  << endl;
  str << "</VTKFile>" << endl;

  str.close() ;


  return ;
};

// =================================================================================== //
template <class Derived>
void VTK_UnstructuredGrid<Derived>::ReadMetaData( ){

    fstream str;
    string line, temp;
    
    fstream::pos_type        position;
    
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
    
    Get_After_Keyword( line, "NumberOfPoints", '\"', temp) ;
    convert_string( temp, nr_points );
    
    Get_After_Keyword( line, "NumberOfCells", '\"', temp) ;
    convert_string( temp, nr_cells );
    
    
    position = str.tellg() ;
    ReadDataHeader( str ) ;
    
    
    for( int i=0; i<geometry.size(); ++i){
        str.seekg( position) ;
        if( ! ReadDataArray( str, geometry[i] ) ) {
          cout << geometry[i].GetName() << " DataArray not found" << endl ;
        };
    };
    
    
    str.close() ;
    
    SetDimensions( nr_cells, nr_points, CalcSizeConnectivity( )  ) ;

    return ;
};

// =================================================================================== //
template <class Derived>
uint64_t VTK_UnstructuredGrid<Derived>::GetNConnectivity( ){

  return nconnectivity ;
};

// =================================================================================== //
template <class Derived>
void VTK_UnstructuredGrid<Derived>::Flush( fstream &str, string codex, string name ){

  static_cast<Derived *>(this)->Flush( str, codex, name );
  return ;
};

// =================================================================================== //
template <class Derived>
void VTK_UnstructuredGrid<Derived>::Absorb( fstream &str, string codex, string name ){

  static_cast<Derived *>(this)->Absorb( str, codex, name );
  return ;
};
