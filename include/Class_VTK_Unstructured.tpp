//------------------------------------------------------------------
template <class Derived>
VTK_UnstructuredGrid<Derived>::VTK_UnstructuredGrid( )
                              :VTK() {

  fh.Set_Appendix("vtu");

  geometry.push_back( VTK::Field_C( "xyz",          3, "Float64", "Point" ) ) ;
  geometry.push_back( VTK::Field_C( "offsets",      1, "Int32"  , "Cell"  ) ) ;
  geometry.push_back( VTK::Field_C( "types",        1, "Int32"  , "Cell"  ) ) ;
  geometry.push_back( VTK::Field_C( "connectivity", 1, "Int32"  , "Cell"  ) ) ;

};

//------------------------------------------------------------------
template <class Derived>
VTK_UnstructuredGrid<Derived>::VTK_UnstructuredGrid( string dir_, string name_, string cod_, int ncell_, int npoints_, int nconn_  )
                              :VTK_UnstructuredGrid( ){

  Set_Names( dir_, name_ ) ; 
  Set_Dimensions( ncell_, npoints_, nconn_ ) ; 

  for( int i=0; i<4 ; i++) geometry[i].Set_Codification( cod_ ) ;

  return ;

};


//------------------------------------------------------------------
template <class Derived>
VTK_UnstructuredGrid<Derived>::~VTK_UnstructuredGrid( ) {
    
    Set_Dimensions( -1, -1, -1) ;    
};

//------------------------------------------------------------------
template <class Derived>
void VTK_UnstructuredGrid<Derived>::Set_Dimensions( int ncells_, int npoints_, int nconn_ ){

  ncells        = ncells_ ;
  npoints       = npoints_ ;
  nconnectivity = nconn_ ;

  geometry[0].Set_Elements(npoints) ;
  geometry[1].Set_Elements(ncells) ;
  geometry[2].Set_Elements(ncells) ;
  geometry[3].Set_Elements(nconnectivity) ;

  nr_points = npoints ;
  nr_cells  = ncells ;

  for( int i=0; i< nr_data; i++){
      if( data[i].Get_Location() == "Cell")  data[i].Set_Elements(nr_cells) ;
      if( data[i].Get_Location() == "Point") data[i].Set_Elements(nr_points) ;
  };

  return ;
};

// =================================================================================== //
template <class Derived>
int VTK_UnstructuredGrid<Derived>::numberofelements( int t){

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
void VTK_UnstructuredGrid<Derived>::Write( ){

  fstream str ;

  Calc_Appended_Offsets() ;

  str.open( fh.Get_Name( ), ios::out ) ;

  //Writing XML header
  str << "<?xml version=\"1.0\"?>" << endl;

  //Writing Piece Information
  str << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
  str << "  <UnstructuredGrid>"  << endl;;
  str << "    <Piece  NumberOfPoints=\"" << npoints << "\" NumberOfCells=\"" << ncells << "\">" << endl;

  //Header for Data
  Write_Data_Header( str, false );

  //Wring Geometry Information
  str << "      <Points>" << endl ;;
  Write_DataArray( str, geometry[0] ) ;
  str << "      </Points>" << endl;

  str << "      <Cells>" << endl ;;
  Write_DataArray( str, geometry[1] ) ;
  Write_DataArray( str, geometry[2] ) ;
  Write_DataArray( str, geometry[3] ) ;
  str << "      </Cells>" << endl;

  //Closing Piece
  str << "    </Piece>" << endl;
  str << "  </UnstructuredGrid>"  << endl;

  //Write Appended Section
  Write_All_Appended( str ) ;

  str << "</VTKFile>" << endl;

  str.close() ;

  if( nr_procs != 0  && my_proc == 0)  Write_pvtu() ;

  fh.Increment_Counter() ;

  return ;
};

// =================================================================================== //
template <class Derived>
void VTK_UnstructuredGrid<Derived>::Write_pvtu( ){

  fstream str ;

  FileHandler_C   fhp, fho ;

  fhp = fh ;
  fho = fh ;

  fhp.Set_Parallel(false) ;
  fhp.Set_Appendix("pvtu") ;

  fho.Set_Directory(".") ;

  str.open( fhp.Get_Name( ), ios::out ) ;

  //Writing XML header
  str << "<?xml version=\"1.0\"?>" << endl;

  //Writing Piece Information
  str << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
  str << "  <PUnstructuredGrid GhostLevel=\"0\">"  << endl;;

  //Header for Data
  Write_Data_Header( str, true );

  //Wring Geometry Information
  str << "      <PPoints>" << endl;
  Write_PDataArray( str, geometry[0] ) ;
  str << endl ;
  str << "      </PPoints>" << endl;


  for( int i=0; i<nr_procs; i++){
    fho.Set_Block(i) ;
    str << "    <Piece  Source=\"" << fho.Get_Name() <<  "\"/>" << endl;
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
void VTK_UnstructuredGrid<Derived>::Read( ){

  fstream str;
  string line;

  fstream::pos_type        position;

  bool                     read ;
 

  str.open( fh.Get_Name( ), ios::in ) ;

//  Read_Piece( str ) ;
  getline( str, line);
  while( ! Keyword_In_String( line, "<Piece")){
    getline(str, line);
  };

  convert_string( Get_After_Keyword( line, "NumberOfPoints", '\"'), npoints );
  convert_string( Get_After_Keyword( line, "NumberOfCells", '\"') , ncells );


  Read_Data_Header( str ) ;
  position = str.tellg() ;

//  Read_Geometry_Header( str ) ;
  str.seekg( position) ;
  if( Seek_and_Read( str, "xyz", geometry[0] ) ) {
    geometry[0].Set_Elements(npoints) ;
  }

  else{
    cout << "xyz DataArray not found" << endl ;
  };

  str.seekg( position) ;
  if( Seek_and_Read( str, "offsets", geometry[1] ) ){
    geometry[1].Set_Elements(ncells) ;
  }

  else{
    cout << "offsets DataArray not found" << endl ;
  };

  str.seekg( position) ;
  if( Seek_and_Read( str, "types", geometry[2] ) ){
    geometry[2].Set_Elements(ncells) ;
  }

  else{
    cout << "types DataArray not found" << endl ;
  };

  str.seekg( position) ;
  if( Seek_and_Read( str, "connectivity", geometry[3] ) ){
    geometry[3].Set_Elements(-1) ;
  }

  else{
    cout << "connectivity DataArray not found" << endl ;
  };


  Read_FieldValues( str ) ;

  str.close() ;

  return ;
};

// =================================================================================== //
template <class Derived>
void VTK_UnstructuredGrid<Derived>::Absorb( fstream &str, string codex, string name ){

  static_cast<Derived *>(this)->Flush( str, codex, name );
  return ;
};
