//------------------------------------------------------------------
template <class Derived>
VTK_UnstructuredGrid<Derived>::VTK_UnstructuredGrid( )
                              :VTK() {

  fh.SetAppendix("vtu");

  geometry.push_back( VTK::Field_C( "xyz",          3, "Float64", "Point" ) ) ;
  geometry.push_back( VTK::Field_C( "offsets",      1, "Int32"  , "Cell"  ) ) ;
  geometry.push_back( VTK::Field_C( "types",        1, "Int32"  , "Cell"  ) ) ;
  geometry.push_back( VTK::Field_C( "connectivity", 1, "Int32"  , "Cell"  ) ) ;

};

//------------------------------------------------------------------
template <class Derived>
VTK_UnstructuredGrid<Derived>::VTK_UnstructuredGrid( string dir_, string name_, string cod_, int ncell_, int npoints_, int nconn_  )
                              :VTK_UnstructuredGrid( ){

  SetNames( dir_, name_ ) ; 
  SetDimensions( ncell_, npoints_, nconn_ ) ; 

  SetGeomCodex( cod_ ) ;

  return ;

};

//------------------------------------------------------------------
template <class Derived>
VTK_UnstructuredGrid<Derived>::~VTK_UnstructuredGrid( ) {
    
    SetDimensions( -1, -1, -1) ;    
};

//------------------------------------------------------------------
template <class Derived>
void VTK_UnstructuredGrid<Derived>::SetDimensions( int ncells_, int npoints_, int nconn_ ){

  ncells        = ncells_ ;
  npoints       = npoints_ ;
  nconnectivity = nconn_ ;

  geometry[0].SetElements(npoints) ;
  geometry[1].SetElements(ncells) ;
  geometry[2].SetElements(ncells) ;
  geometry[3].SetElements(nconnectivity) ;

  nr_points = npoints ;
  nr_cells  = ncells ;

  for( int i=0; i< nr_data; i++){
      if( data[i].GetLocation() == "Cell")  data[i].SetElements(nr_cells) ;
      if( data[i].GetLocation() == "Point") data[i].SetElements(nr_points) ;
  };

  return ;
};

// =================================================================================== //
template <class Derived>
int VTK_UnstructuredGrid<Derived>::NumberOfElements( int t){

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
  str << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
  str << "  <UnstructuredGrid>"  << endl;;
  str << "    <Piece  NumberOfPoints=\"" << npoints << "\" NumberOfCells=\"" << ncells << "\">" << endl;

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
void VTK_UnstructuredGrid<Derived>::Flush( fstream &str, string codex, string name ){

  static_cast<Derived *>(this)->Flush( str, codex, name );
  return ;
};

// =================================================================================== //
template <class Derived>
void VTK_UnstructuredGrid<Derived>::ReadMetaData( ){

  fstream str;
  string line;

  fstream::pos_type        position;

  bool                     read ;
 

  str.open( fh.GetName( ), ios::in ) ;

//  Read_Piece( str ) ;
  getline( str, line);
  while( ! Keyword_In_String( line, "<Piece")){
    getline(str, line);
  };

  convert_string( Get_After_Keyword( line, "NumberOfPoints", '\"'), npoints );
  convert_string( Get_After_Keyword( line, "NumberOfCells", '\"') , ncells );


  position = str.tellg() ;
  ReadDataHeader( str ) ;


  for( int i=0; i<geometry.size(); ++i){
      str.seekg( position) ;
      if( ! ReadDataArray( str, geometry[i] ) ) {
        cout << geometry[i].GetName() << " DataArray not found" << endl ;
      };
  };

  SetDimensions( ncells, npoints, -1 ) ;

  str.close() ;

  return ;
};

// =================================================================================== //
template <class Derived>
void VTK_UnstructuredGrid<Derived>::Absorb( fstream &str, string codex, string name ){

  static_cast<Derived *>(this)->Absorb( str, codex, name );
  return ;
};
