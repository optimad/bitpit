#ifndef __FileHandler__HH__
#define __FileHandler__HH__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>


/*!
 * \class   FileHandler
 * \brief   Creates file names and checks status
 *
 * Permits to create easily filenames for parallel output and for writing time series.
 * Typical filenames are *directory/name.bxxxx.yyyy.appendix* and can be retirieved via GetName().
 * - *directory* set via SetDirectory( std::string)
 * - *name* set via SetName(std::string)
 * - *bxxxx* is added if SetParallel(true) and SetBlock(xxxx) have been called
 * - *yyyy* is added if SetSeries(true) and SetCounter(yyyy) have been called. yyyy is incremented by one if IncrementCounter() is called.
 */
class FileHandler{

    private:
        std::string           directory;       /**< name od directory where file resides */
        std::string           name;            /**< name of file; */
        std::string           appendix ;       /**< appendix; */
        bool                  series ;         /**< is time series? */
        bool                  parallel ;       /**< is part of distributed data? */
        int                   counter ;        /**< counter for time series ; */
        int                   block  ;         /**< index of block if parallel data */


    public:
        FileHandler() ;
        FileHandler( std::string dir_, std::string name_, std::string app_) ;
        FileHandler(  const FileHandler& other ) ;

        ~FileHandler() ;

        FileHandler& operator=( const FileHandler& other) ;

        void                setDirectory( std::string d_) ;
        void                setName( std::string n_) ;
        void                setAppendix( std::string a_) ;
        void                setSeries( bool s_) ;
        void                setParallel( bool p_) ;
        void                setCounter(int c_);
        void                setBlock( int b_) ;

        std::string         getName() ;

        void                incrementCounter();
        bool                exists() ;

};


#endif
