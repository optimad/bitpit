
#include"Class_VTK_Wrappers.hpp"

using namespace std;

VtkUnstrVec::VtkUnstrVec()
             :VTK_UnstructuredGrid<VtkUnstrVec>(){};

//====================================================================================================
VtkUnstrVec::VtkUnstrVec( string dir_, string name_, string codex_, uint8_t type_ )
             :VTK_UnstructuredGrid<VtkUnstrVec>( ){ 

    SetNames( dir_, name_ );
    SetCodex( codex_ ) ;

    type = type_ ;

};

//====================================================================================================
VtkUnstrVec::~VtkUnstrVec(){

};

// =================================================================================== //
bool VtkUnstrVec::GetFieldByName( const string &name_, VtkUnstrVec::ufield *&the_field ){

    vector<VtkUnstrVec::ufield>::iterator    it_ ;

    for( it_=adata.begin(); it_!=adata.end(); ++it_){
        if( (*it_).name == name_ ){
            the_field = &(*it_) ;
            return true ;
        };
    };

    return false ;
};

//====================================================================================================
void VtkUnstrVec::Write(  ) {

    bool            CanWrite(true) ;
    ufield*         FPtr ; 
    string          name ;
    vector<string>  TBD;

    CanWrite = CanWrite && GetFieldByName( "Points", FPtr ) ; 
    CanWrite = CanWrite && GetFieldByName( "connectivity", FPtr ) ; 

    if( CanWrite){

        for( unsigned i=0; i<nr_data; ++i){
            name = data[i].GetName() ;
            if( !GetFieldByName(name,FPtr) ){
                TBD.push_back(name) ;
            };
        };

        for( unsigned i=0; i<TBD.size(); ++i){
            VTK::RemoveData( TBD[i] ) ;
        };

        VTK::Write() ;
    }

    else{
        cout << " VtkUnstrVec::Write cannot write file since \"Points\" or \"connectivity\" vectors are missing" << endl ;

    };

    return;
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

        if( GetFieldByName( name, f_)){
            stream_visitor  visitor;
            visitor.SetStream( str ) ;
            visitor.SetCodex( codex ) ;
            visitor.SetName( name ) ;
            visitor.SetTask( "write" ) ;

            boost::apply_visitor(visitor, f_->DPtr ); 
        };

    };

    return ;

};

//====================================================================================================
void VtkUnstrVec::Absorb( fstream &str, string codex, string name ) {

    ufield  *f_ ;

    if( GetFieldByName( name, f_)){

        VTK::Field_C*   FPtr ;
        if( VTK::GetFieldByName(name, FPtr) ){

            stream_visitor  visitor;
            visitor.SetStream( str ) ;
            visitor.SetCodex( codex ) ;
            visitor.SetName( name ) ;
            visitor.SetTask( "read" ) ;

            visitor.SetSize( FPtr->GetElements() ) ;
            visitor.SetComponents( FPtr->GetComponents() ) ;

            if( name == "connectivity") visitor.SetComponents( NumberOfElements(type) ) ;

            boost::apply_visitor(visitor, f_->DPtr ); 
        };
    }

    return ;

};
