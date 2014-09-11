#ifndef __CLASS_FH__HH__
#define __CLASS_FH__HH__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "Operators.hpp"

using namespace std;

class FileHandler_C{

    protected:
      string              directory;       // name od directory where file resides
      string              name;            // name of file;
      string              appendix ;       // appendix;
      bool                series ;         // is time series?
      bool                parallel ;       // is part of distributed data?
      int                 counter ;        // counter for time series ;
      int                 block  ;         // index of block if parallel data

    public:
       FileHandler_C() ;
       FileHandler_C( string dir_, string name_, string app_) ;
      ~FileHandler_C() ;

      void                Set_Directory( string d_) ;
      void                Set_Name( string n_) ;
      void                Set_Appendix( string a_) ;
      void                Set_Series( bool s_) ;
      void                Set_Parallel( bool p_) ;
      void                Set_Counter(int c_);
      void                Set_Block( int b_) ;

      void                Increment_Counter();
      string              Get_Name() ;
      bool                Exists() ;

};


#endif
