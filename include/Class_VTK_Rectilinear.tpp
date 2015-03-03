
//------------------------------------------------------------------
template <class Derived>
VTK_RectilinearGrid<Derived>::VTK_RectilinearGrid( )
                             :VTK() {

  fh.Set_Appendix( "vtr" );

  geometry.push_back( VTK::Field_C( "x_Coord", 1, "Float64", "Point" ) ) ;
  geometry.push_back( VTK::Field_C( "y_Coord", 1, "Float64", "Point" ) ) ;
  geometry.push_back( VTK::Field_C( "z_Coord", 1, "Float64", "Point" ) ) ;

} ;

//------------------------------------------------------------------
template <class Derived>
VTK_RectilinearGrid<Derived>::VTK_RectilinearGrid( string dir_, string name_ )
                             :VTK(dir_, name_ ) {

  fh.Set_Appendix( "vtr" );

  geometry.push_back( VTK::Field_C( "x_Coord", 1, "Float64", "Point" ) ) ;
  geometry.push_back( VTK::Field_C( "y_Coord", 1, "Float64", "Point" ) ) ;
  geometry.push_back( VTK::Field_C( "z_Coord", 1, "Float64", "Point" ) ) ;

  return ;
};


//------------------------------------------------------------------
template <class Derived>
VTK_RectilinearGrid<Derived>::VTK_RectilinearGrid( string dir_, string name_, string codex_, int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ )
                             :VTK(dir_, name_ ) {

  fh.Set_Appendix( "vtr" );

  n1 = n1_ ;
  n2 = n2_ ;
  m1 = m1_ ;
  m2 = m2_ ;
  l1 = l1_ ;
  l2 = l2_ ;

  nr_cells =  (n2-n1+0)*(m2-m1+0)*(l2-l1+0) ;
  nr_points = (n2-n1+1)*(m2-m1+1)*(l2-l1+1) ;

  geometry.push_back( VTK::Field_C( "x_Coord", 1, "Float64", "Point", codex_, n2-n1+1 ) ) ;
  geometry.push_back( VTK::Field_C( "y_Coord", 1, "Float64", "Point", codex_, m2-m1+1 ) ) ;
  geometry.push_back( VTK::Field_C( "z_Coord", 1, "Float64", "Point", codex_, l2-l1+1 ) ) ;

  return ;
};

//------------------------------------------------------------------
template <class Derived>
VTK_RectilinearGrid<Derived>::~VTK_RectilinearGrid( ) {} ;

// =================================================================================== //
template <class Derived>
void VTK_RectilinearGrid<Derived>::Read( ){

  fstream str;
  string line;

  fstream::pos_type        position;

  bool                     read ;

  array<int,6>             extensions ;


  str.open( fh.Get_Name( ), ios::in ) ;

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

  Read_Data_Header( str ) ;
  position = str.tellg() ;

//  Read_Geometry_Header( str ) ;
  str.seekg( position) ;
  if( Seek_and_Read( str, "x_Coord", geometry[0] ) ) {
    geometry[0].Set_Elements(n2-n1+1) ;
  }

  else{
    cout << "x_Coord DataArray not found" << endl ;
  };

  str.seekg( position) ;
  if( Seek_and_Read( str, "y_Coord", geometry[1] ) ){
    geometry[1].Set_Elements(m2-m1+1) ;
  }

  else{
    cout << "y_Coord DataArray not found" << endl ;
  };

  str.seekg( position) ;
  if( Seek_and_Read( str, "z_Coord", geometry[2] ) ){
    geometry[2].Set_Elements(l2-l1+1) ;
  }

  else{
    cout << "z_Coord DataArray not found" << endl ;
  };

  Read_FieldValues( str ) ;

  str.close() ; 

  return ;
 
};

// =================================================================================== //
template <class Derived>
void VTK_RectilinearGrid<Derived>::Write( ){

  fstream str;

  Calc_Appended_Offsets();

  str.open( fh.Get_Name( ), ios::out | ios::app ) ;

  //Writing XML header
  str << "<?xml version=\"1.0\"?>" << endl;

  //Writing Piece Information
  str << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
  str << "  <RectilinearGrid WholeExtent= \"" <<n1<<" "<<n2<<" "<<m1<<" "<<m2<<" "<<l1<<" "<<l2<< "\" >" << endl;
  str << "    <Piece Extent= \" " <<n1<<" "<<n2<<" "<<m1<<" "<<m2<<" "<<l1<<" "<<l2<< "\" >" << endl;


  //Header for Data
  Write_Data_Header( str, false ) ;

  //Wring Geometry Information   
  str << "       <Coordinates>" << endl;
  Write_DataArray( str, geometry[0] ) ;
  Write_DataArray( str, geometry[1] ) ;
  Write_DataArray( str, geometry[2] ) ;
  str << "       </Coordinates>" << endl;

  //Closing Piece
  str << "    </Piece>" << endl;
  str << "  </RectilinearGrid>" << endl;

  //Write Appended Section
  Write_All_Appended( str ) ;

  //Closing XML
  str << "</VTKFile>" << endl;

  str.close() ;

  fh.Increment_Counter() ;
  if( nr_procs != 0  && my_proc == 0)  Write_pvtr() ;

  return ;
};

// =================================================================================== //
template <class Derived>
void VTK_RectilinearGrid<Derived>::Write_pvtr( ){

  fstream str ;

  FileHandler_C   fhp, fho ;

  fhp = fh ;
  fho = fh ;

  fhp.Set_Parallel(false) ;
  fhp.Set_Appendix("pvtr") ;

  fho.Set_Directory(".") ;

  str.open( fhp.Get_Name( ), ios::out ) ;

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
  Write_Data_Header( str, true );

  //Wring Geometry Information
  str << "      <PCoordinates>" << endl;
  Write_PDataArray( str, geometry[0] ) ;
  Write_PDataArray( str, geometry[1] ) ;
  Write_PDataArray( str, geometry[2] ) ;
  str << "      </PCoordinates>" << endl;


  for( int i=0; i<nr_procs; i++){
    fho.Set_Block(i) ;
    str << "    <Piece Extent= \" " 
        << proc_index[i][0] << " " << proc_index[i][1] << " "
        << proc_index[i][2] << " " << proc_index[i][3] << " "
        << proc_index[i][4] << " " << proc_index[i][5] << " "
        << "\" Source= \"" << fho.Get_Name() << "\"/>" << endl;
  };

  str << "  </PRectilinearGrid>"  << endl;
  str << "</VTKFile>" << endl;

  str.close() ;


  return ;
};

// =================================================================================== //
template <class Derived>
void VTK_RectilinearGrid<Derived>::Set_Parallel_Index( array<int,6> glo_, vector<array<int,6>> loc_ ){

  if( loc_.size() !=nr_procs ) cout << "Size of loc_ in VTK_RectilinearGrid<Derived>::Set_Parallel_Index does not fit nr_procs " << endl ;

  global_index = glo_ ;
  proc_index   = loc_ ;
  
  return;
};

// =================================================================================== //
template <class Derived>
void VTK_RectilinearGrid<Derived>::Set_Dimensions( int n1_, int n2_, int m1_, int m2_, int l1_, int l2_ ){

  n1 = n1_ ;
  n2 = n2_ ;
  m1 = m1_ ;
  m2 = m2_ ;
  l1 = l1_ ;
  l2 = l2_ ;

  geometry[0].Set_Elements(n2-n1+1) ;
  geometry[1].Set_Elements(m2-m1+1) ;
  geometry[2].Set_Elements(l2-l1+1) ;

  nr_cells =  (n2-n1+0)*(m2-m1+0)*(l2-l1+0) ;
  nr_points = (n2-n1+1)*(m2-m1+1)*(l2-l1+1) ;

  for( int i=0; i< nr_data; i++){
    if( data[i].Get_Location() == "Cell")  data[i].Set_Elements(nr_cells) ;
    if( data[i].Get_Location() == "Point") data[i].Set_Elements(nr_points) ;
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

