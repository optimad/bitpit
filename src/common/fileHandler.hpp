/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

#include "bitpit_common.hpp"

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
        int                   counter ;        /**< counter for time series ; */

        bool                  parallel ;       /**< is part of distributed data? */
        int                   block  ;         /**< index of the parallel block for distributed data */

    public:
        BITPIT_PUBLIC_API static const int INVALID_COUNTER = -1; /**< dummy counter associate with files that doesn't belong to a time series*/

        FileHandler() ;
        FileHandler( const std::string & dir_, const std::string & name_, const std::string & app_) ;
        FileHandler(  const FileHandler& other ) ;

        FileHandler& operator=( const FileHandler& other) ;

        bool                exists() const ;

        std::string         getPath() const ;
        const std::string & getDirectory() const ;
        const std::string & getName() const ;
        const std::string & getAppendix() const ;
        void                setDirectory( const std::string & d_) ;
        void                setName( const std::string & n_) ;
        void                setAppendix( const std::string & a_) ;

        bool                isSeries( ) const ;
        void                setSeries( bool s_) ;

        int                 getCounter() const ;
        void                setCounter(int c_);
        void                incrementCounter();

        bool                isParallel( ) const ;
        void                setParallel( bool p_) ;

        int                 getBlock( ) const ;
        void                setBlock( int b_) ;

};

}

#endif
