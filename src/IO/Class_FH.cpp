#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include"Class_FH.hpp"
#include"Operators.hpp"

/*! ---------------------------------------------------------------
 * Default constructor. 
 * Refers to a serial file in the current directory called "file.dat"
 * Default values of counters and blocks are set to 0.
 */
FileHandler_C::FileHandler_C(){
   
    directory  = "./" ; 
    name       = "file" ;
    appendix   = "dat" ;
    series     = false ;
    parallel   = false ;
    counter    = 0 ;
    block      = 0 ;
    
} ;

/*! ---------------------------------------------------------------
 * Constructor. 
 * Refers to a serial file named according to input parameters.
 * Default values of counters and blocks are set to 0.
 * \param[in]   dir_    directory of file
 * \param[in]   name_   name of the file
 * \param[in]   app_    appendix of the file
 */
FileHandler_C::FileHandler_C( std::string dir_, std::string name_, std::string app_):
    directory(dir_), name(name_), appendix(app_){
    series     = false ;
    parallel   = false ;
    counter    = 0 ;
    block      = 0 ;
};

/*! ---------------------------------------------------------------
 * Copy Constructor. 
 * \param[in]   other    FileHandler_C to be copied.
 */
FileHandler_C::FileHandler_C(const FileHandler_C& other){

    (*this) = other ;

};

/*! ---------------------------------------------------------------
 * Destructor.
 */
FileHandler_C::~FileHandler_C(){

} ;

/*! ---------------------------------------------------------------
 * Assignment operator.
 */
FileHandler_C&    FileHandler_C::operator=(const FileHandler_C& other){

    directory   = other.directory ;
    name        = other.name ;
    appendix    = other.appendix ;
    series      = other.series ;
    parallel    = other.parallel ;
    counter     = other.counter ;
    block       = other.block ;


    return *this ;
};

/*! ---------------------------------------------------------------
 * Sets the directory
 * \param[in]   d_      name of directory
 */
void    FileHandler_C::SetDirectory( std::string d_){ 
    directory=d_; 
    return;
} ;

/*! ---------------------------------------------------------------
 * Sets the name
 * \param[in]   n_      name of file
 */
void    FileHandler_C::SetName( std::string n_){ 
    name=n_; 
    return;
} ;

/*! ---------------------------------------------------------------
 * Sets the appendix
 * \param[in]   a_      appendix of file
 */
void    FileHandler_C::SetAppendix( std::string a_){ 
    appendix=a_; 
    return;
} ;

/*! ---------------------------------------------------------------
 * Sets if file belongs to a time series
 * \param[in]   s_      [true/false] if series
 */
void    FileHandler_C::SetSeries( bool s_){ 
    series=s_; 
    return;
} ;

/*! ---------------------------------------------------------------
 * Sets the counter of the series
 * \param[in]   c_      index of series
 */
void    FileHandler_C::SetCounter(int c_){ 
    counter=c_; 
    return; 
};

/*! ---------------------------------------------------------------
 * Sets if file belongs to a parallel output
 * \param[in]   s_      [true/false] if parallel
 */
void    FileHandler_C::SetParallel( bool p_){ 
    parallel=p_; 
    return;
} ;

/*! ---------------------------------------------------------------
 * Sets the counter of the series
 * \param[in]   v_      index of block of parallel set
 */
void    FileHandler_C::SetBlock( int b_){ 
    block=b_; 
    return;
} ;

/*! ---------------------------------------------------------------
 * Increments the counter used for series.
 */
void    FileHandler_C::IncrementCounter(){ 
    counter++; 
    return; 
};

/*! ---------------------------------------------------------------
 * Composes the filename.
 * \return  complete filename
 */
std::string  FileHandler_C::GetName(){
  std::stringstream filename ;

  filename << directory << "/"<<name ;
  if(parallel) filename <<".b"<< ZeroPadNumber(4, block)    ;
  if(series)   filename <<"." << ZeroPadNumber(4, counter)  ;
  filename <<"."<< appendix  ;

  return filename.str() ;

};

/*! ---------------------------------------------------------------
 * Checks if file exists.
 * \return  [true/false] if file exists
 */
bool   FileHandler_C::Exists() {
    std::ifstream f( GetName() );
    return f.good() ;
};
