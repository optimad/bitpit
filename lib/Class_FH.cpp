#include"Class_FH.hpp"

//---------------------------------------------------------------
FileHandler_C::FileHandler_C(){
   
    directory  = "./" ; 
    name       = "file" ;
    appendix   = "dat" ;
    series     = false ;
    parallel   = false ;
    counter    = 999 ;
    block      = 999 ;
    
} ;

//---------------------------------------------------------------
FileHandler_C::FileHandler_C( string dir_, string name_, string app_):
    directory(dir_), name(name_), appendix(app_){
    series     = false ;
    parallel   = false ;
    counter    = 999 ;
    block      = 999 ;
};

//---------------------------------------------------------------
FileHandler_C::FileHandler_C(const FileHandler_C& other){

    (*this) = other ;

};

//---------------------------------------------------------------
FileHandler_C::~FileHandler_C(){

    if (str.is_open()) {
        str.close() ;
        
    };
} ;

//---------------------------------------------------------------
FileHandler_C&    FileHandler_C::operator=(const FileHandler_C& other){

    directory   = other.directory ;
    name        = other.name ;
    appendix    = other.appendix ;
    series      = other.series ;
    parallel    = other.parallel ;
    counter     = other.counter ;
    block       = other.block ;

    str.clear() ;

};

//---------------------------------------------------------------
void    FileHandler_C::Set_Directory( string d_){ directory=d_; return;} ;

//---------------------------------------------------------------
void    FileHandler_C::Set_Name( string n_){ name=n_; return;} ;

//---------------------------------------------------------------
void    FileHandler_C::Set_Appendix( string a_){ appendix=a_; return;} ;

//---------------------------------------------------------------
void    FileHandler_C::Set_Series( bool s_){ series=s_; return;} ;

//---------------------------------------------------------------
void    FileHandler_C::Set_Parallel( bool p_){ parallel=p_; return;} ;

//---------------------------------------------------------------
void    FileHandler_C::Set_Counter(int c_){ counter=c_; return; };

//---------------------------------------------------------------
void    FileHandler_C::Set_Block( int b_){ block=b_; return;} ;

//---------------------------------------------------------------
void    FileHandler_C::Increment_Counter(){ if(series) counter++; return; };

//---------------------------------------------------------------
string  FileHandler_C::Get_Name(){
  stringstream filename ;

  filename << directory << "/"<<name ;
  if(parallel) filename <<".b"<< ZeroPadNumber(4, block)    ;
  if(series)   filename <<"." << ZeroPadNumber(4, counter)  ;
  filename <<"."<< appendix  ;

  return filename.str() ;

};
//---------------------------------------------------------------
bool   FileHandler_C::Exists() {
    ifstream f( Get_Name() );
    return f.good() ;
};
