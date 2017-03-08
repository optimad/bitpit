/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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


#include <iostream>

#include "GenericIO.hpp"


int main()
{

    int exitStatus(0) ;

    { //write ASCII
        bitpit::PiercedVector<int>  dataOut ;    

        dataOut.insert(0,4) ;
        dataOut.insert(4,7) ;
        dataOut.insert(2,10) ;

        std::fstream    file;
        file.open( "writeASCII.dat", std::ios::out ) ;
        bitpit::genericIO::flushASCII(file, 3, dataOut) ;
        file.close() ;

        file.open( "writeASCII_INDEX.dat", std::ios::out ) ;
        bitpit::genericIO::flushASCII(file, 3, dataOut, true) ;
        file.close() ;

    }

    { //read ASCII and write BINARY
        bitpit::PiercedVector<int>  dataNoId, dataId ;    

        dataNoId.reclaim(0) ;
        dataNoId.reclaim(4) ;
        dataNoId.reclaim(2) ;

        std::fstream    file;
        file.open( "writeASCII.dat", std::ios::in ) ;
        bitpit::genericIO::absorbASCII(file, dataNoId) ;
        file.close() ;

        file.open( "writeASCII_INDEX.dat", std::ios::in ) ;
        bitpit::genericIO::absorbASCII(file, dataId, 3) ;
        file.close() ;

        file.open( "writeBINARY.dat", std::ios::out | std::ios::binary ) ;
        bitpit::genericIO::flushBINARY(file, dataNoId) ;
        file.close() ;

        file.open( "writeBINARY_INDEX.dat", std::ios::out | std::ios::binary ) ;
        bitpit::genericIO::flushBINARY(file, dataId, true) ;
        file.close() ;
    }

    { //read BINARY and compare with original data
        bitpit::PiercedVector<int>  dataNoId, dataId ;    

        dataNoId.reclaim(0) ;
        dataNoId.reclaim(4) ;
        dataNoId.reclaim(2) ;

        std::fstream    file;
        file.open( "writeBINARY.dat", std::ios::in | std::ios::binary ) ;
        bitpit::genericIO::absorbBINARY(file, dataNoId ) ;
        file.close() ;

        file.open( "writeBINARY_INDEX.dat", std::ios::in | std::ios::binary ) ;
        bitpit::genericIO::absorbBINARY(file, dataId, 3) ;
        file.close() ;

        if( dataNoId[0] != 4 || dataNoId[4] != 7 || dataNoId[2] != 10 ){
            std::cout << "read /write of bitpit::PiercedVector<long> without Ids not successful" << std::endl;
            exitStatus++;

        } else if( dataId[0] != 4 || dataId[4] != 7 || dataId[2] != 10 ){
            std::cout << "read /write of bitpit::PiercedVector<long> with Ids not successful" << std::endl ;
            exitStatus++;

        }


    }


    return exitStatus;

}
