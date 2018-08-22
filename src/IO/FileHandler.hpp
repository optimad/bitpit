/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#ifndef __BITPIT_FILEHANDLER_HPP__
#define __BITPIT_FILEHANDLER_HPP__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

namespace bitpit{

/*!
 * \class   FileHandler
 * \brief   Creates file names and checks status
 *
 * Permits to create easily filenames for parallel output and for writing time series.
 * Typical filenames are *directory/name.bxxxx.yyyy.appendix* and can be retirieved via GetName().
 * - *directory* set via SetDirectory( const std::string &)
 * - *name* set via SetName(const std::string &)
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
        FileHandler( const std::string & dir_, const std::string & name_, const std::string & app_) ;
        FileHandler(  const FileHandler& other ) ;

        ~FileHandler() ;

        FileHandler& operator=( const FileHandler& other) ;

        void                setDirectory( const std::string & d_) ;
        void                setName( const std::string & n_) ;
        void                setAppendix( const std::string & a_) ;
        void                setSeries( bool s_) ;
        void                setParallel( bool p_) ;
        void                setCounter(int c_);
        void                setBlock( int b_) ;

        std::string         getPath() ;
        std::string         getName() const ;
        std::string         getDirectory() const ;
        std::string         getAppendix() const ;
        int                 getCounter() const ;

        void                incrementCounter();
        bool                exists() ;

};

}

#endif
