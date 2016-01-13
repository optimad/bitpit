#ifndef __CLASS_FH__HH__
#define __CLASS_FH__HH__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>


/*!
 * \class   FileHandler_C
 * \brief   Creates file names and checks status
 *
 * Permits to create easily filenames for parallel output and for writing time series.
 * Typical filenames are *directory/name.bxxxx.yyyy.appendix* and can be retirieved via GetName().
 * - *directory* set via SetDirectory( std::string)
 * - *name* set via SetName(std::string)
 * - *bxxxx* is added if SetParallel(true) and SetBlock(xxxx) have been called
 * - *yyyy* is added if SetSeries(true) and SetCounter(yyyy) have been called. yyyy is incremented by one if IncrementCounter() is called.
 */
class FileHandler_C{

    private:
        std::string           directory;       /**< name od directory where file resides */
        std::string           name;            /**< name of file; */
        std::string           appendix ;       /**< appendix; */
        bool                  series ;         /**< is time series? */
        bool                  parallel ;       /**< is part of distributed data? */
        int                   counter ;        /**< counter for time series ; */
        int                   block  ;         /**< index of block if parallel data */


    public:
        FileHandler_C() ;
        FileHandler_C( std::string dir_, std::string name_, std::string app_) ;
        FileHandler_C(  const FileHandler_C& other ) ;

        ~FileHandler_C() ;

        FileHandler_C& operator=( const FileHandler_C& other) ;

        void                SetDirectory( std::string d_) ;
        void                SetName( std::string n_) ;
        void                SetAppendix( std::string a_) ;
        void                SetSeries( bool s_) ;
        void                SetParallel( bool p_) ;
        void                SetCounter(int c_);
        void                SetBlock( int b_) ;

        void                IncrementCounter();
        std::string         GetName() ;
        bool                Exists() ;

};


#endif
