#include"my_vtk.hpp"

//====================================================================================================
my_vtk_unstr::my_vtk_unstr()
             :VTK_UnstructuredGrid<my_vtk_unstr>(){};

//====================================================================================================
my_vtk_unstr::my_vtk_unstr( string dir_, string name_, string codex_, int ncells_, int npoints_, int nconn_ )
             :VTK_UnstructuredGrid<my_vtk_unstr>( dir_, name_, codex_, ncells_, npoints_, nconn_){ };

//====================================================================================================
my_vtk_unstr::~my_vtk_unstr(){

  points       = NULL;
  connectivity = NULL;
};

//====================================================================================================
void my_vtk_unstr::Link_Data( vector<vec1> &points_ext, vector<node_v> &connectivity_ext ){

  points       = &points_ext;
  connectivity = &connectivity_ext ;

  return;
};

//====================================================================================================
void my_vtk_unstr::Flush( fstream &str, string codex_, string name ) {

  int n;
  string indent("         ") ;


    if( codex_ == "ascii"){

      if( name == "xyz"){
        for( n=0; n<8; n++) {
          flush_ascii( str, indent ) ;
          flush_ascii( str, 3, points->at(n)  ) ;
          str << endl ;
        };
      };

      if( name == "connectivity"){
        for( n=0; n<1; n++) {
          flush_ascii( str, indent ) ;
          flush_ascii( str, 8, connectivity->at(n)  ) ;
          str << endl ;
        };
      };

      if( name == "types"){
        int type_(11) ;
        for( n=0; n<1; n++) {
          flush_ascii( str, indent ) ;
          flush_ascii( str, type_  ) ;
          str << endl ;
        };
      };

      if( name == "offsets"){
        int off_(0) ;
        for( n=0; n<1; n++) {
          off_ += numberofelements( 11 ) ; 
          flush_ascii( str, indent ) ;
          flush_ascii( str, off_  ) ;
          str << endl ;
        };
      };

    }

    else{

      if( name == "xyz"){
        for( n=0; n<8; n++) flush_binary( str, points->at(n)  ) ;
      };

      if( name == "connectivity"){
        for( n=0; n<1; n++) flush_binary( str, connectivity->at(n)  ) ;
      };

      if( name == "types"){
        int type_(11) ;
        for( n=0; n<1; n++) flush_binary( str, type_  ) ;
      };

      if( name == "offsets"){
        int off_(0) ;
        for( n=0; n<1; n++) {
          off_ += numberofelements( 11 ) ; 
          flush_binary( str, off_  ) ;
        };
      };

    };

  return ;

};

//====================================================================================================
void my_vtk_unstr::Absorb( fstream &str, string codex_, string name ) {

  int n;


    if( codex_ == "ascii"){


      if( name == "xyz"){
        for( n=0; n<8; n++) {
          absorb_ascii( str, points->at(n)  ) ;
        };
      };

      if( name == "connectivity"){
        for( n=0; n<1; n++) {
          absorb_ascii( str, connectivity->at(n)  ) ;
        };
      };

      if( name == "types"){
        int type_ ;
        for( n=0; n<1; n++) {
          absorb_ascii( str, type_  ) ;
        };
      };

      if( name == "offsets"){
        int off_ ;
        for( n=0; n<1; n++) {
          absorb_ascii( str, off_  ) ;
        };
      };

    }

    else{

      if( name == "xyz"){
        for( n=0; n<8; n++) absorb_binary( str, points->at(n)  ) ;
      };

      if( name == "connectivity"){
        for( n=0; n<1; n++) absorb_binary( str, connectivity->at(n)  ) ;
      };

      if( name == "types"){
        int type_ ;
        for( n=0; n<1; n++) absorb_binary( str, type_  ) ;
      };

      if( name == "offsets"){
        int off_ ;
        for( n=0; n<1; n++) {
          absorb_binary( str, off_  ) ;
        };
      };

    };

  return ;

};
