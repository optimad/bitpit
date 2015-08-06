#include"Class_VTK_Wrappers.hpp"

//====================================================================================================
VtkUnstrVec::VtkUnstrVec()
             :VTK_UnstructuredGrid<VtkUnstrVec>(){};

//====================================================================================================
VtkUnstrVec::VtkUnstrVec( string dir_, string name_, string codex_, int type_, dvecarr3E &points_ext, ivector2D &connectivity_ext )
             :VTK_UnstructuredGrid<VtkUnstrVec>( ){ 


    int ncells_, npoints_, nconn_;

    SetNames( dir_, name_ );
    SetGeomCodex( codex_ ) ;

    type = type_ ;
    codex= codex_ ;

    points = &points_ext ;
    connectivity = &connectivity_ext ;

    ncells_ = connectivity_ext.size() ;
    npoints_ = points_ext.size() ;

    nconn_ = ncells_ * NumberOfElements( type ) ;

    SetDimensions( ncells_, npoints_, nconn_ ) ;
};

//====================================================================================================
VtkUnstrVec::~VtkUnstrVec(){

  points       = NULL;
  connectivity = NULL;
};

//====================================================================================================
void VtkUnstrVec::AddData( dvector1D &data_, string name_, string loc_ ){

    int n = scalar_data.size() ;
    VTK::Field_C*     f_ ;

    scalar_data.push_back( sfield() ) ;

    scalar_data[n].name = name_ ;
    scalar_data[n].data = &data_ ;
  
    if( ! GetFieldByName( data, name_, f_ )) {
        VTK_UnstructuredGrid<VtkUnstrVec>::AddData( name_, 1, "Float64", loc_, codex ) ;
    };

    return;
};

//====================================================================================================
void VtkUnstrVec::AddData( dvecarr3E &data_, string name_, string loc_ ){


    int n = vector_data.size() ;
    VTK::Field_C*     f_ ;

    vector_data.push_back( vfield() ) ;

    vector_data[n].name = name_ ;
    vector_data[n].data = &data_ ;

    if( !GetFieldByName( data, name_, f_ ) ) {
        VTK_UnstructuredGrid<VtkUnstrVec>::AddData( name_, 3, "Float64", loc_, codex ) ;
    };

    return;
};

//====================================================================================================
void VtkUnstrVec::Flush( fstream &str, string codex_, string name ) {

  int n;
  string indent("         ") ;


    if( codex_ == "ascii"){

      if( name == "xyz"){
          flush_ascii( str, 1, (*points)  ) ;
      }

      else if( name == "connectivity"){
         flush_ascii( str, 1, (*connectivity)  ) ;
      }

      else if( name == "types"){
        for( n=0; n<ncells-1; n++) {
          flush_ascii( str, type  ) ;
          str << endl ;
        };
        flush_ascii( str, type  ) ;
      }

      else if( name == "offsets"){
        int off_(0) ;
        for( n=0; n<ncells-1; n++) {
          off_ += NumberOfElements( type ) ; 
          flush_ascii( str, off_  ) ;
          str << endl ;
        };
        off_ += NumberOfElements( type ) ;
        flush_ascii( str, off_  ) ;
      }

      else{

          for( int i=0; i<scalar_data.size(); ++i ){
              if( scalar_data[i].name == name) flush_ascii( str, 1, *(scalar_data[i].data) ) ;
          };

          for( int i=0; i<vector_data.size(); ++i ){
              if( vector_data[i].name == name) flush_ascii( str, 1, *(vector_data[i].data) ) ;
          };
      };

    }

    else{

      if( name == "xyz"){
        flush_binary( str, *points  ) ;
      }

      else if( name == "connectivity"){
        flush_binary( str, *connectivity  ) ;
      }

      else if( name == "types"){
        for( n=0; n<ncells; n++) flush_binary( str, type  ) ;
      }

      else if( name == "offsets"){
        int off_(0) ;
        for( n=0; n<ncells; n++) {
          off_ += NumberOfElements( type ) ; 
          flush_binary( str, off_  ) ;
        };
      }

      else{

          for( int i=0; i<scalar_data.size(); ++i ){
              if( scalar_data[i].name == name) flush_binary( str, *(scalar_data[i].data) ) ;
          };

          for( int i=0; i<vector_data.size(); ++i ){
              if( vector_data[i].name == name) flush_binary( str, *(vector_data[i].data) ) ;
          };
      };

    };

  return ;

};

//====================================================================================================
void VtkUnstrVec::Absorb( fstream &str, string codex_, string name ) {

  int n;



    if( codex_ == "ascii"){

      if( name == "xyz"){
        (*points).resize( npoints ) ;
         absorb_ascii( str, (*points) ) ;
      }

      else if( name == "connectivity"){
        (*connectivity).resize( ncells ) ;
        for( n=0; n<ncells; n++) {
          (*connectivity)[n].resize( NumberOfElements(type) ) ;
        };
          absorb_ascii( str, (*connectivity)  ) ;
      }

      else if( name == "types"){
        for( n=0; n<ncells; n++) {
          absorb_ascii( str, type  ) ;
        };
      }

      else if( name == "offsets"){
        int off_ ;
        for( n=0; n<ncells; n++) {
          absorb_ascii( str, off_  ) ;
        };
      }

      else{
          VTK::Field_C*     f_ ;
          string            loc; 
          for( int i=0; i<scalar_data.size(); ++i ){
              if( scalar_data[i].name == name) {
                  GetFieldByName( data, name, f_ );
                  loc = f_->GetLocation() ;
                  
                  if( loc == "Point") (scalar_data[i].data)->resize(npoints) ;
                  if( loc == "Cell")  (scalar_data[i].data)->resize(ncells) ;

                  absorb_ascii( str, *(scalar_data[i].data) ) ;
              };
          };

          for( int i=0; i<vector_data.size(); ++i ){
              if( vector_data[i].name == name) {
                  GetFieldByName( data, name, f_ );
                  loc = f_->GetLocation() ;
                  
                  if( loc == "Point") (vector_data[i].data)->resize(npoints) ;
                  if( loc == "Cell")  (vector_data[i].data)->resize(ncells) ;

                  absorb_ascii( str, *(vector_data[i].data) ) ;
              };

          };


      };

    }

    else{

      if( name == "xyz"){
        (*points).resize( npoints ) ;
        absorb_binary( str, (*points) ) ;
      }

      else if( name == "connectivity"){
        (*connectivity).resize( ncells ) ;
        for( n=0; n<ncells; n++) {
          (*connectivity)[n].resize( NumberOfElements(type) ) ;
        };
        absorb_binary( str, (*connectivity)  ) ;
      }

      else if( name == "types"){
      }

      else if( name == "offsets"){
      }

      else{
          VTK::Field_C*     f_ ;
          string            loc; 
          for( int i=0; i<scalar_data.size(); ++i ){
              if( scalar_data[i].name == name) {
                  GetFieldByName( data, name, f_ );
                  loc = f_->GetLocation() ;
                  
                  if( loc == "Point") (scalar_data[i].data)->resize(npoints) ;
                  if( loc == "Cell")  (scalar_data[i].data)->resize(ncells) ;

                  absorb_binary( str, *(scalar_data[i].data) ) ;
              };
          };

          for( int i=0; i<vector_data.size(); ++i ){
              if( vector_data[i].name == name) {
                  GetFieldByName( data, name, f_ );
                  loc = f_->GetLocation() ;
                  
                  if( loc == "Point") (vector_data[i].data)->resize(npoints) ;
                  if( loc == "Cell")  (vector_data[i].data)->resize(ncells) ;

                  absorb_binary( str, *(vector_data[i].data) ) ;
              };

          };


      };

    };

  return ;

};
