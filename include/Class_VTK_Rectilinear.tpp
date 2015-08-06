
//------------------------------------------------------------------
template <class Derived>
VTK_RectilinearGrid<Derived>::VTK_RectilinearGrid( )
                             :VTK() {

  fh.SetAppendix( "vtr" );

  geometry.push_back( VTK::Field_C( "x_Coord", 1, "Float64", "Point" ) ) ;
  geometry.push_back( VTK::Field_C( "y_Coord", 1, "Float64", "Point" ) ) ;
  geometry.push_back( VTK::Field_C( "z_Coord", 1, "Float64", "Point" ) ) ;

} ;

//------------------------------------------------------------------
template <class Derived>
VTK_RectilinearGrid<Derived>::VTK_RectilinearGrid( string dir_, string name_, string codex_, int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ )
                             :VTK_RectilinearGrid( ) {

  SetNames( dir_, name_ ) ;

  SetDimensions( n1_, n2_, m1_, m2_, l1_, l2_) ;
  SetGeomCodex( codex_ ) ;

  return ;
};

//------------------------------------------------------------------
template <class Derived>
VTK_RectilinearGrid<Derived>::~VTK_RectilinearGrid( ) {} ;

// =================================================================================== //
template <class Derived>
void VTK_RectilinearGrid<Derived>::ReadMetaData( ){

  fstream str;
  string line;

  fstream::pos_type        position;

  bool                     read ;

  array<int,6>             extensions ;


  str.open( fh.GetName( ), ios::in ) ;

  getline( str, line);
  while( ! Keyword_In_String( line, "<Piece")){
    getline(str, line);
  };

  convert_string( Get_After_Keyword( line, "Extent", '\"'), extensions );
 
  n1 = extensions[0] ;
  n2 = extensions[1] ;
  m1 = extensions[2] ;
  m2 = extensions[3] ;
  l1 = extensions[4] ;
  l2 = extensions[5] ;

  position = str.tellg() ;

  ReadDataHeader( str ) ;

    for( int i=0; i<geometry.size(); ++i){
        str.seekg( position) ;
        if( ! ReadDataArray( str, geometry[i] ) ) {
            cout << geometry[i].GetName() << " DataArray not found" << endl ;
        };
    };


  SetDimensions( n1, n2, m1, m2, l1, l2 ) ;
  str.close() ; 

  return ;
 
};

// =================================================================================== //
template <class Derived>
void VTK_RectilinearGrid<Derived>::WriteMetaData( ){

  fstream str;

  str.open( fh.GetName( ), ios::out ) ;

  //Writing XML header
  str << "<?xml version=\"1.0\"?>" << endl;

  //Writing Piece Information
  str << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
  str << "  <RectilinearGrid WholeExtent= \"" <<n1<<" "<<n2<<" "<<m1<<" "<<m2<<" "<<l1<<" "<<l2<< "\" >" << endl;
  str << "    <Piece Extent= \" " <<n1<<" "<<n2<<" "<<m1<<" "<<m2<<" "<<l1<<" "<<l2<< "\" >" << endl;


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

// =================================================================================== //
template <class Derived>
void VTK_RectilinearGrid<Derived>::WriteCollection( ){

  fstream str ;

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
      << global_index[0] << " " << global_index[1]<< " "
      << global_index[2] << " " << global_index[3]<< " "
      << global_index[4] << " " << global_index[5]<< " "
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
    str << "    <Piece Extent= \" " 
        << proc_index[i][0] << " " << proc_index[i][1] << " "
        << proc_index[i][2] << " " << proc_index[i][3] << " "
        << proc_index[i][4] << " " << proc_index[i][5] << " "
        << "\" Source= \"" << fho.GetName() << "\"/>" << endl;
  };

  str << "  </PRectilinearGrid>"  << endl;
  str << "</VTKFile>" << endl;

  str.close() ;


  return ;
};

// =================================================================================== //
template <class Derived>
void VTK_RectilinearGrid<Derived>::SetParallelIndex( array<int,6> glo_, vector<array<int,6>> loc_ ){

  if( loc_.size() !=nr_procs ) cout << "Size of loc_ in VTK_RectilinearGrid<Derived>::SetParallelIndex does not fit nr_procs " << endl ;

  global_index = glo_ ;
  proc_index   = loc_ ;
  
  return;
};

// =================================================================================== //
template <class Derived>
void VTK_RectilinearGrid<Derived>::SetDimensions( int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ ){

  n1 = n1_ ;
  n2 = n2_ ;
  m1 = m1_ ;
  m2 = m2_ ;
  l1 = l1_ ;
  l2 = l2_ ;

  geometry[0].SetElements(n2-n1+1) ;
  geometry[1].SetElements(m2-m1+1) ;
  geometry[2].SetElements(l2-l1+1) ;

  nr_cells =  (n2-n1+0)*(m2-m1+0)*(l2-l1+0) ;
  nr_points = (n2-n1+1)*(m2-m1+1)*(l2-l1+1) ;

  for( int i=0; i< nr_data; i++){
    if( data[i].GetLocation() == "Cell")  data[i].SetElements(nr_cells) ;
    if( data[i].GetLocation() == "Point") data[i].SetElements(nr_points) ;
  };

  return ;
};

// =================================================================================== //
template <class Derived>
void VTK_RectilinearGrid<Derived>::Flush( fstream &str, string codex, string name ){

  static_cast<Derived *>(this)->Flush( str, codex, name );
  return ;

};

// =================================================================================== //
template <class Derived>
void VTK_RectilinearGrid<Derived>::Absorb( fstream &str, string codex, string name ){

  static_cast<Derived *>(this)->Absorb( str, codex, name );
  return ;

};

