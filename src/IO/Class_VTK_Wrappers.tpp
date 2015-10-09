#include"Class_VTK_Wrappers.hpp"

//====================================================================================================
template< class T0, class T1>
VtkUnstrVec::VtkUnstrVec( string dir_, string name_, string codex_, uint8_t type_, vector<T0> &points_ext, vector<T1> &connectivity_ext )
             :VTK_UnstructuredGrid<VtkUnstrVec>( ){ 


    int ncells_, npoints_, nconn_;

    T0  dum0 ;
    T1  dum1 ;

    SetNames( dir_, name_ );
    SetCodex( codex_ ) ;

    type = type_ ;

    ncells_ = connectivity_ext.size() ;
    npoints_ = points_ext.size() ;

    nconn_ = ncells_ * NumberOfElements( type ) ;

    SetGeomTypes( WhichType(dum0), "UInt64", "UInt8", WhichType(dum1)  ) ;
    SetDimensions( ncells_, npoints_, nconn_ ) ;

    adata.resize(2) ;

    adata[0].DPtr = &points_ext ;
    adata[0].name = "Points" ;

    adata[1].DPtr = &connectivity_ext ;
    adata[1].name = "connectivity" ;

};

//====================================================================================================
template<class T>
void VtkUnstrVec::AddData( vector<T> &data_, string name_, string loc_ ){

    T               dum_ ;

    adata.push_back( ufield(name_,data_) ) ;

    VTK_UnstructuredGrid<VtkUnstrVec>::AddData( name_, 1, WhichType(dum_), loc_ ) ;

    return;
};

//====================================================================================================
template<class T>
void VtkUnstrVec::AddData( vector< array<T,3> > &data_, string name_, string loc_ ){

    T               dum_ ;

    adata.push_back( ufield(name_,data_) ) ;

    VTK_UnstructuredGrid<VtkUnstrVec>::AddData( name_, 3, WhichType(dum_), loc_ ) ;

    return;
};

//====================================================================================================
template<class T>
void VtkUnstrVec::AddData( vector< vector<T> > &data_, string name_, string loc_ ){

    T               dum_ ;

    adata.push_back( ufield(name_,data_) ) ;

    VTK_UnstructuredGrid<VtkUnstrVec>::AddData( name_, 3, WhichType(dum_), loc_ ) ;


    return;
};

//====================================================================================================
template<class T>
void VtkUnstrVec::LinkData( vector<T> &data_, string name_ ){

    adata.push_back( ufield(name_,data_) ) ;

    return;
};
