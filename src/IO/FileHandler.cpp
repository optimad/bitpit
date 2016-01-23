#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include"FileHandler.hpp"
#include"Operators.hpp"

/*!
 * Default constructor. 
 * Refers to a serial file in the current directory called "file.dat"
 * Default values of counters and blocks are set to 0.
 */
FileHandler::FileHandler(){
   
    directory  = "./" ; 
    name       = "file" ;
    appendix   = "dat" ;
    series     = false ;
    parallel   = false ;
    counter    = 0 ;
    block      = 0 ;
    
} ;

/*!
 * Constructor. 
 * Refers to a serial file named according to input parameters.
 * Default values of counters and blocks are set to 0.
 * @param[in]   dir_    directory of file
 * @param[in]   name_   name of the file
 * @param[in]   app_    appendix of the file
 */
FileHandler::FileHandler( std::string dir_, std::string name_, std::string app_):
    directory(dir_), name(name_), appendix(app_){
    series     = false ;
    parallel   = false ;
    counter    = 0 ;
    block      = 0 ;
};

/*!
 * Copy Constructor. 
 * @param[in]   other    FileHandler to be copied.
 */
FileHandler::FileHandler(const FileHandler& other){

    (*this) = other ;

};

/*!
 * Destructor.
 */
FileHandler::~FileHandler(){

} ;

/*!
 * Assignment operator.
 */
FileHandler&    FileHandler::operator=(const FileHandler& other){

    directory   = other.directory ;
    name        = other.name ;
    appendix    = other.appendix ;
    series      = other.series ;
    parallel    = other.parallel ;
    counter     = other.counter ;
    block       = other.block ;


    return *this ;
};

/*!
 * sets the directory
 * @param[in]   d_      name of directory
 */
void    FileHandler::setDirectory( std::string d_){ 
    directory=d_; 
    return;
} ;

/*!
 * sets the name
 * @param[in]   n_      name of file
 */
void    FileHandler::setName( std::string n_){ 
    name=n_; 
    return;
} ;

/*!
 * sets the appendix
 * @param[in]   a_      appendix of file
 */
void    FileHandler::setAppendix( std::string a_){ 
    appendix=a_; 
    return;
} ;

/*!
 * sets if file belongs to a time series
 * @param[in]   s_      [true/false] if series
 */
void    FileHandler::setSeries( bool s_){ 
    series=s_; 
    return;
} ;

/*!
 * sets the counter of the series
 * @param[in]   c_      index of series
 */
void    FileHandler::setCounter(int c_){ 
    counter=c_; 
    return; 
};

/*!
 * sets if file belongs to a parallel output
 * @param[in]   s_      [true/false] if parallel
 */
void    FileHandler::setParallel( bool p_){ 
    parallel=p_; 
    return;
} ;

/*!
 * sets the counter of the series
 * @param[in]   v_      index of block of parallel set
 */
void    FileHandler::setBlock( int b_){ 
    block=b_; 
    return;
} ;

/*!
 * Increments the counter used for series.
 */
void    FileHandler::incrementCounter(){ 
    counter++; 
    return; 
};

/*!
 * Composes the filename.
 * @return  complete filename
 */
std::string  FileHandler::getName(){
  std::stringstream filename ;

  filename << directory << "/"<<name ;
  if(parallel) filename <<".b"<< ZeroPadNumber(4, block)    ;
  if(series)   filename <<"." << ZeroPadNumber(4, counter)  ;
  filename <<"."<< appendix  ;

  return filename.str() ;

};

/*!
 * Checks if file exists.
 * @return  [true/false] if file exists
 */
bool   FileHandler::exists() {
    std::ifstream f( getName() );
    return f.good() ;
};
