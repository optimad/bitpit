#include"Class_VTK_Wrappers.hpp"

//====================================================================================================
VtkUnstrVec::VtkUnstrVec()
             :VTK_UnstructuredGrid<VtkUnstrVec>(){};

// //====================================================================================================
// VtkUnstrVec::VtkUnstrVec( string dir_, string name_, string codex_, uint8_t type_, dvecarr3E &points_ext, ivector2D &connectivity_ext )
//              :VTK_UnstructuredGrid<VtkUnstrVec>( ){ 
// 
// 
//     int ncells_, npoints_, nconn_;
// 
//     data.reserve(50) ;
// 
//     SetNames( dir_, name_ );
//     SetCodex( codex_ ) ;
// 
//     type = type_ ;
// 
//     points = &points_ext ;
//     connectivity = &connectivity_ext ;
// 
//     ncells_ = connectivity_ext.size() ;
//     npoints_ = points_ext.size() ;
// 
//     nconn_ = ncells_ * NumberOfElements( type ) ;
// 
//     geometry[0].SetType("Float64") ;
//     geometry[1].SetType("UInt32") ;
//     geometry[2].SetType("UInt8") ;
//     geometry[3].SetType("UInt32") ;
// 
//     SetDimensions( connectivity_ext.size(), points_ext.size(), nconn_ ) ;
// };

//====================================================================================================
template< class T0, class T1>
VtkUnstrVec::VtkUnstrVec( string dir_, string name_, string codex_, uint8_t type_, vector<T0> &points_ext, vector<T1> &connectivity_ext )
             :VTK_UnstructuredGrid<VtkUnstrVec>( ){ 


    int ncells_, npoints_, nconn_;

    T0  dum0 ;
    T1  dum1 ;


    data.reserve(50) ;

    SetNames( dir_, name_ );
    SetCodex( codex_ ) ;

    type = type_ ;

    ncells_ = connectivity_ext.size() ;
    npoints_ = points_ext.size() ;

    nconn_ = ncells_ * NumberOfElements( type ) ;

    geometry[0].SetType( WhichType(dum0) ) ;
    geometry[1].SetType("UInt64") ;
    geometry[2].SetType("UInt8") ;
    geometry[3].SetType( WhichType(dum1) ) ;

    SetDimensions( connectivity_ext.size(), points_ext.size(), nconn_ ) ;

    adata.resize(4) ;

    adata[0].DPtr = &points_ext ;
    adata[1].DPtr = static_cast<vector<int32_t>*>(NULL) ;
    adata[2].DPtr = static_cast<vector<int32_t>*>(NULL) ;
    adata[3].DPtr = &connectivity_ext ;

    adata[0].FPtr = &geometry[0] ;
    adata[1].FPtr = &geometry[1] ;
    adata[2].FPtr = &geometry[2] ;
    adata[3].FPtr = &geometry[3] ;

};

//====================================================================================================
VtkUnstrVec::VtkUnstrVec( string dir_, string name_, string codex_, uint8_t type_ )
             :VTK_UnstructuredGrid<VtkUnstrVec>( ){ 

    data.reserve(50) ;

    SetNames( dir_, name_ );
    SetCodex( codex_ ) ;

    type = type_ ;

};


//====================================================================================================
VtkUnstrVec::~VtkUnstrVec(){

};

//====================================================================================================
template<class T>
void VtkUnstrVec::AddData( vector<T> &data_, string name_, string loc_ ){

    int n = adata.size() ;
    VTK::Field_C*     f_ ;
    T               dum_ ;

    adata.push_back( ufield() ) ;

    adata[n].DPtr = &data_ ;

    if( name_ == "Points" ){
        adata[n].FPtr = &geometry[0] ;
    }

    else if( name_ == "types"){
        adata[n].FPtr = &geometry[1] ;
    }

    else if( name_ == "offsets"){
        adata[n].FPtr = &geometry[2] ;
    }

    else if( name_ == "connectivity"){
        adata[n].FPtr = &geometry[3] ;
    }

    else{
        adata[n].FPtr = VTK_UnstructuredGrid<VtkUnstrVec>::AddData( name_, 1, WhichType(dum_), loc_ ) ;
    };

    return;
};

//====================================================================================================
template<class T>
void VtkUnstrVec::AddData( vector< array<T,3> > &data_, string name_, string loc_ ){

    int n = adata.size() ;
    VTK::Field_C*     f_ ;
    T               dum_ ;

    adata.push_back( ufield() ) ;

    adata[n].DPtr = &data_ ;

    if( name_ == "Points" ){
        adata[n].FPtr = &geometry[0] ;
    }

    else if( name_ == "types"){
        adata[n].FPtr = &geometry[1] ;
    }

    else if( name_ == "offsets"){
        adata[n].FPtr = &geometry[2] ;
    }

    else if( name_ == "connectivity"){
        adata[n].FPtr = &geometry[3] ;
    }

    else{

        adata[n].FPtr = VTK_UnstructuredGrid<VtkUnstrVec>::AddData( name_, 3, WhichType(dum_), loc_ ) ;
    };

    return;
};

//====================================================================================================
template<class T>
void VtkUnstrVec::AddData( vector< vector<T> > &data_, string name_, string loc_ ){

    int n = adata.size() ;
    VTK::Field_C*     f_ ;
    T               dum_ ;

    adata.push_back( ufield() ) ;

    adata[n].DPtr = &data_ ;

    if( name_ == "Points" ){
        adata[n].FPtr = &geometry[0] ;
    }

    else if( name_ == "types"){
        adata[n].FPtr = &geometry[1] ;
    }

    else if( name_ == "offsets"){
        adata[n].FPtr = &geometry[2] ;
    }

    else if( name_ == "connectivity"){
        adata[n].FPtr = &geometry[3] ;
    }

    else{
        adata[n].FPtr = VTK_UnstructuredGrid<VtkUnstrVec>::AddData( name_, 3, WhichType(dum_), loc_ ) ;
    };


    return;
};

// =================================================================================== //
bool VtkUnstrVec::GetFieldByName( const string &name_, VtkUnstrVec::ufield *&the_field ){


    vector<VtkUnstrVec::ufield>::iterator    it_ ;

    for( it_=adata.begin(); it_!=adata.end(); ++it_){
        if( (*it_).FPtr->GetName() == name_ ){
            the_field = &(*it_) ;
            return true ;
        };
    };

    return false ;
};

//====================================================================================================
void VtkUnstrVec::Flush( fstream &str, string codex, string name ) {


    if( codex == "ascii" && name == "types"){
        for( uint64_t n=0; n<nr_cells-1; n++) {
            flush_ascii( str, type  ) ;
            str << endl ;
        };
        flush_ascii( str, type  ) ;
    }

    else if( codex == "binary" && name == "types"){
        for( uint64_t n=0; n<nr_cells; n++) flush_binary( str, type  ) ;
    }

    else if( codex == "ascii" && name == "offsets"){
        uint64_t off_(0) ;
        for(uint64_t  n=0; n<nr_cells-1; n++) {
            off_ += NumberOfElements( type ) ; 
            flush_ascii( str, off_  ) ;
            str << endl ;
        };
        off_ += NumberOfElements( type ) ;
        flush_ascii( str, off_  ) ;
    }

    else if( codex == "binary" && name == "offsets"){
        uint64_t off_(0), nT( NumberOfElements( type ) ) ;
        for( uint64_t n=0; n<nr_cells; n++) {
            off_ += nT ; 
            flush_binary( str, off_  ) ;
        };
    }

    else{

        ufield  *f_ ;
        stream_visitor  visitor;
        visitor.SetStream( str ) ;
        visitor.SetCodex( codex ) ;
        visitor.SetName( name ) ;
        visitor.SetTask( "write" ) ;

        GetFieldByName( name, f_) ;
        boost::apply_visitor(visitor, f_->DPtr ); 

    };

    return ;

};

//====================================================================================================
void VtkUnstrVec::Absorb( fstream &str, string codex, string name ) {

    ufield  *f_ ;

    if( GetFieldByName( name, f_)){

        stream_visitor  visitor;
        visitor.SetStream( str ) ;
        visitor.SetCodex( codex ) ;
        visitor.SetName( name ) ;
        visitor.SetTask( "read" ) ;

        visitor.SetSize( f_->FPtr->GetElements() ) ;
        visitor.SetComponents( f_->FPtr->GetComponents() ) ;

        if( name == "connectivity") visitor.SetComponents( NumberOfElements(type) ) ;

        if( ! (boost::get<vector<int32_t>*>(&f_->DPtr) && boost::get<vector<int32_t>*>(f_->DPtr) == NULL ) ){
            boost::apply_visitor(visitor, f_->DPtr ); 
        };
    }

    return ;

};
