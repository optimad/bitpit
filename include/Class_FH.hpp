#ifndef __CLASS_FH__HH__
#define __CLASS_FH__HH__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "Operators.hpp"

using namespace std;

class FileHandler_C{

    public:
    fstream               str;

    protected:
    string                directory;       // name od directory where file resides
    string                name;            // name of file;
    string                appendix ;       // appendix;
    bool                  series ;         // is time series?
    bool                  parallel ;       // is part of distributed data?
    int                   counter ;        // counter for time series ;
    int                   block  ;         // index of block if parallel data


    public:
    FileHandler_C() ;
    FileHandler_C( string dir_, string name_, string app_) ;
    FileHandler_C(  const FileHandler_C& other ) ;

    ~FileHandler_C() ;
    
    FileHandler_C& operator=( const FileHandler_C& other) ;

    void                SetDirectory( string d_) ;
    void                SetName( string n_) ;
    void                SetAppendix( string a_) ;
    void                SetSeries( bool s_) ;
    void                SetParallel( bool p_) ;
    void                SetCounter(int c_);
    void                SetBlock( int b_) ;

    void                IncrementCounter();
    string              GetName() ;
    bool                Exists() ;

};


#endif
