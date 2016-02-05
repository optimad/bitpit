/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

#include"Operators.hpp"

namespace bitpit{

namespace genericIO{

/*!
 * @ingroup GenericIO
 * @{
 */

/*!
 * Writes a POD data type to file stream in ascii format using scientific format
 * @tparam  data_T  type of POD data
 * @param[in]   str     file stream
 * @param[in]   data    data to be written
 */
template< class data_T >
void flushASCII( std::fstream &str, const data_T &data ){

    str << std::setprecision(8) << std::scientific ;
    str << data << " ";

    return ;
};

/*!
 * Writes a vector of POD data type to file stream in ascii format using scientific format on a single line
 * @tparam  data_T  type of POD data
 * @param[in]   str     file stream
 * @param[in]   data    data to be written
 */
template< class data_T >
void flushASCII( std::fstream &str, const std::vector<data_T> &data ){

    flushASCII( str, data.size(), data) ;

    return ;
};

/*!
 * Writes a vector of POD data type to file stream in ascii format using scientific format putting a fixed number of elements per line
 * @tparam  data_T  type of POD data
 * @param[in]   str     file stream
 * @param[in]   elements_per_line     number of entries per line
 * @param[in]   data    data to be written
 */
template< class data_T >
void flushASCII( std::fstream &str, int elements_per_line, const std::vector<data_T> &data ){

    int i(0), j(0), k(0) ;
    int nr ;

    int lines, this_line ;

    bool next(true) ;

    nr = data.size() ;
    lines = (nr-1) /elements_per_line + 1;

    str << std::setprecision(8) << std::scientific ;


    while( next ) {

        this_line = std::min( elements_per_line, nr - k ) ;

        for( j=0; j<this_line; j++){
            flushASCII( str, data[k] ) ;
            k++ ;
        };

        i++ ;
        if( i<lines){
            str << std::endl;
        }
        else{
            next = false;
        };

    };

    return ;
};

/*!
 * Writes a array of POD data type to file stream in ascii format using scientific format on one single line
 * @tparam  data_T  type of POD data
 * @tparam  d       size of array
 * @param[in]       str     file stream
 * @param[in]       data    data to be written
 */
template< class data_T, size_t d >
void flushASCII( std::fstream &str, const std::array<data_T,d> &data ){

    flushASCII( str, d, data ) ;
    return ;
};

/*!
 * Writes a std::array of POD data type to file stream in ascii format using scientific format putting a fixed number of elements per line
 * @tparam          data_T  type of POD data
 * @tparam          d       size of array
 * @param[in]       str     file stream
 * @param[in]       elements_per_line     number of entries per line
 * @param[in]       data    data to be written
 */
template< class data_T, size_t d >
void flushASCII( std::fstream &str, int elements_per_line, const std::array<data_T,d> &data ){

    int i(0), j(0), k(0) ;
    int nr ;

    int lines, this_line ;

    bool next(true) ;

    nr = d ;
    lines = nr /elements_per_line ;

    str << std::setprecision(8) << std::scientific ;


    while( next ) {

        this_line = std::min( elements_per_line, nr - k ) ;
        for( j=0; j<this_line; j++){
            flushASCII( str, data[k] ) ;
            k++ ;
        };

        i++ ;
        if( i<lines){
            str << std::endl;
        }
        else{
            next = false;
        };

    };

    return ;
};

/*!
 * Writes a C array of POD data type to file stream in ascii format using scientific format putting a fixed number of elements per line
 * @tparam          data_T  type of POD data
 * @param[in]       str     file stream
 * @param[in]       elements_per_line     number of entries per line
 * @param[in]       data    data to be written
 * @param[in]       nr      size of the C array
 */
template< class data_T >
void flushASCII( std::fstream &str, int elements_per_line, const data_T *data, int nr ){

    int i(0), j(0), k(0) ;

    int lines, this_line ;

    bool next(true) ;

    lines = nr /elements_per_line ;

    str << std::setprecision(8) << std::scientific ;


    while( next ) {

        this_line = std::min( elements_per_line, nr - k ) ;
        for( j=0; j<this_line; j++){
            flushASCII( str, data[k] ) ;
            k++ ;
        };

        i++ ;
        if( i<lines){
            str << std::endl;
        }
        else{
            next = false;
        };

    };

    return ;
};

/*!
 * Writes a POD data type to file stream in binary format.
 * The function makes uses of memory contigiuty.
 * @tparam          data_T  type of POD data
 * @param[in]       str     file stream
 * @param[in]       data    data to be written
 */
template< class data_T >
void flushBINARY( std::fstream &str, const data_T &data ){

    int nbytes;
    nbytes = sizeof(data_T) ;

    str.write( reinterpret_cast<const char*>(&data), nbytes ) ;

    return ;
};

/*!
 * Writes a std::vector of POD data type to file stream in binary format.
 * The function makes uses of memory contigiuty.
 * @tparam          data_T  type of POD data
 * @param[in]       str     file stream
 * @param[in]       data    data to be written
 */
template< class data_T >
void flushBINARY( std::fstream &str, const std::vector<data_T> &data ){

    int nbytes, nr;
    nr = data.size() ;
    nbytes = sizeof(data_T) *nr ;

    str.write( reinterpret_cast<const char*>(&data[0]), nbytes ) ;

    return ;
};

/*!
 * Writes a std::vector<std::vector> of POD data type to file stream in binary format.
 * @tparam          data_T  type of POD data
 * @param[in]       str     file stream
 * @param[in]       data    data to be written
 */
template< class data_T >
void flushBINARY( std::fstream &str, const std::vector< std::vector<data_T> > &data ){


    for( const auto &item : data){
        flushBINARY( str, item ) ;
    };

    return ;
};

/*!
 * Writes a std::vector<std::array> of POD data type to file stream in binary format.
 * The function makes uses of memory contigiuty.
 * @tparam          data_T  type of POD data
 * @tparam          d       size of the array
 * @param[in]       str     file stream
 * @param[in]       data    data to be written
 */
template< class data_T, size_t d >
void flushBINARY( std::fstream &str, const std::vector< std::array<data_T,d> > &data ){

    int nbytes, nr;
    nr = data.size() ;
    nbytes = sizeof(data_T) *nr *d ;

    str.write( reinterpret_cast<const char*>(&data[0]), nbytes ) ;

    return ;
};

/*!
 * Writes a std::array of POD data type to file stream in binary format 
 * @tparam          data_T  type of POD data
 * @tparam          d       size of the array
 * @param[in]       str     file stream
 * @param[in]       data    data to be written
 */
template< class data_T, size_t d >
void flushBINARY( std::fstream &str, const std::array<data_T,d> &data ){

    int nbytes;
    nbytes = sizeof(data_T)*d ;

    str.write( reinterpret_cast<const char*>(&data[0]), nbytes ) ;

    return ;
};

/*!
 * Writes a C array of POD data type to file stream in binary format 
 * @tparam          data_T  type of POD data
 * @param[in]       str     file stream
 * @param[in]       data    data to be written
 * @param[in]       nr      size of the C array
 */
template< class data_T >
void flushBINARY( std::fstream &str, const data_T *data, int nr ){

    int nbytes;
    nbytes = sizeof(data_T) *nr ;

    str.write( reinterpret_cast<const char*>(data), nbytes ) ;

    return ;
};

/*!
 * Reads one line into templated data type.
 * Relies on operator ">>" of the templated data type.
 * In case the information on the line is not sufficient or exceeds the data type and error message is displayed std::cout
 * @tparam          data_T  type of POD data
 * @param[in]       str     file stream
 * @param[in]       data    data to be written
 */
template< class data_T >
void  lineStream( std::fstream &str, data_T &data){

    std::vector<data_T> temp;
    data_T         x_;
    std::string         line;
    int            expected, read(0) ;

    expected = 1;

    getline( str, line );
    bitpit::utils::trim( line );

    std::stringstream ss( line );

    while( ss.good() ){
        ss >> x_ ;
        temp.push_back(x_);
        read++ ;
    };

    if( read != expected){
        std::cout << " Not expected nr of element in line" << std::endl;
        std::cout << " Expected number: "<< expected << std::endl ; 
        std::cout << " Actual number: "<< read << std::endl ; 
    }

    else{
        data=temp[0];
    };

    return;

};

/*!
 * Reads one line into std::vector of templated data type.
 * Relies on operator ">>" of the templated data type.
 * In case data.size() == 0, data will be resized to match all information available on the line.
 * Otherwise the information on the line should fit exactly the size of data.
 * If not an error message is displayed std::cout
 * @tparam          data_T  type of POD data
 * @param[in]       str     file stream
 * @param[in]       data    data to be written
 */
template< class data_T >
void  lineStream( std::fstream &str, std::vector<data_T> &data){

    std::vector<data_T>   temp;
    data_T                x_;
    std::string           line;
    int                   expected, read(0) ;

    expected = data.size() ;

    getline( str, line );
    bitpit::utils::trim( line );

    std::stringstream ss( line );

    while( ss.good() ){
        ss >> x_ ;
        temp.push_back(x_);
        read++ ;
    };

    if( expected == 0) {
        data = temp ;
    }

    else{
        if( expected == read){
            data = temp ;
        }
        else{
            std::cout << " Not expected nr of element in line" << std::endl;
            std::cout << " Expected number: "<< expected << std::endl ; 
            std::cout << " Actual number: "<< read << std::endl ; 
        };
    };

    return;

};

/*!
 * Reads one line into std::array of templated data type and templated size
 * Relies on operator ">>" of the templated data type.
 * The information on the line should fit exactly the size of the array.
 * If not an error message is displayed std::cout
 * @tparam          data_T  type of POD data
 * @tparam          d       size of std::array
 * @param[in]       str     file stream
 * @param[in]       data    data to be written
 */
template< class data_T, size_t d >
void  lineStream( std::fstream &str, std::array<data_T,d> &data){

    std::vector<data_T> temp;
    data_T              x_;
    std::string         line;
    int                 expected, read(0), i ;

    expected = d ;

    getline( str, line );
    bitpit::utils::trim( line );

    std::stringstream ss( line );

    while( ss.good() ){
        ss >> x_ ;
        temp.push_back(x_);
        read++ ;
    };

    if( expected == read){
        for(i=0; i<read; i++) data[i] = temp[i] ;
    }
    else{
        std::cout << " Not expected nr of element in line" << std::endl;
        std::cout << " Expected number: "<< expected << std::endl ; 
        std::cout << " Actual number: "<< read << std::endl ; 
    };

    return;

};

/*!
 * Reads one line into C array of templated data type and given size.
 * Relies on operator ">>" of the templated data type.
 * The information on the line should fit exactly the size .
 * If not an error message is displayed std::cout
 * @tparam          data_T  type of POD data
 * @param[in]       str     file stream
 * @param[in]       data    data to be written
 * @param[in]       nr      number of elements to be read
 */
template< class data_T >
void  lineStream( std::fstream &str, data_T *data, int nr ){

    std::vector<data_T> temp;
    data_T              x_;
    std::string         line;
    int                 expected, read(0), i ;

    expected = nr ;

    getline( str, line );
    bitpit::utils::trim( line );

    std::stringstream ss( line );

    while( ss.good() ){
        ss >> x_ ;
        temp.push_back(x_);
        read++ ;
    };

    if( expected == read){
        for(i=0; i<read; i++) data[i] = temp[i] ;
    }
    else{
        std::cout << " Not expected nr of element in line" << std::endl;
        std::cout << " Expected number: "<< expected << std::endl ; 
        std::cout << " Actual number: "<< read << std::endl ; 
    };

    return;

};

/*!
 * Reads a templated data type from file stream in ascii 
 * Relies on the ">>" operator.
 * @tparam  data_T  type of POD data
 * @param[in]   str     file stream
 * @param[in]   data    data to be read
 */
template< class data_T >
void absorbASCII( std::fstream &str, data_T &data ){

    str >> data ;

    return ;
};

/*!
 * Reads a std::vector of data type from file stream in ascii.
 * The size of the vector defines the number of elements to be read.
 * If not enough elements are present in the file an error message is displayed on std::cout
 * Relies on the function lineStream.
 * @tparam  data_T  type of POD data
 * @param[in]   str     file stream
 * @param[in]   data    data to be read
 */
template< class data_T >
void absorbASCII( std::fstream &str, std::vector<data_T> &data ){

    std::vector<data_T>             temp;

    typename std::vector<data_T>::iterator   itrData, begData, endData;
    typename std::vector<data_T>::iterator   itrTemp, begTemp, endTemp;

    begData = data.begin() ;
    endData = data.end() ;

    itrData = begData ;

    while( str.good() && itrData!=endData ) {

        temp.clear() ;
        lineStream( str, temp) ;

        begTemp = temp.begin() ;
        endTemp = temp.end() ;

        for( itrTemp=begTemp; (itrTemp!=endTemp && itrData!=endData); ){
            *itrData = *itrTemp ;

            ++itrTemp;
            ++itrData;

        };


    };

    if( itrData != endData ) {
        std::cout << "Not enough elements found to fill vector" << std::endl ;
    };


    return ;
};

/*!
 * Reads a std::array of data type from file stream in ascii.
 * The size of the array defines the number of elements to be read.
 * If not enough elements are present in the file an error message is displayed on std::cout
 * Relies on the function lineStream.
 * @tparam      data_T  class stored in std::array
 * @tparam      d       size of std::array
 * @param[in]   str     file stream
 * @param[in]   data    data to be read
 */
template< class data_T, size_t d >
void absorbASCII( std::fstream &str, std::array<data_T,d> &data ){

    std::vector<data_T>             temp;

    typename std::array<data_T,d>::iterator  itrData, begData, endData;
    typename std::vector<data_T>::iterator   itrTemp, begTemp, endTemp;

    begData = data.begin() ;
    endData = data.end() ;

    itrData = begData ;

    while( str.good() && itrData!=endData ) {

        temp.clear() ;
        lineStream( str, temp) ;

        begTemp = temp.begin() ;
        endTemp = temp.end() ;

        for( itrTemp=begTemp; (itrTemp!=endTemp && itrData!=endData); ){
            *itrData = *itrTemp ;

            ++itrTemp;
            ++itrData;

        };


    };

    if( itrData != endData ) {
        std::cout << "Not enough elements found to fill array" << std::endl ;
    };



    return ;
};

/*!
 * Reads a C array of data type from file stream in ascii.
 * If not enough elements are present in the file an error message is displayed on std::cout
 * Relies on the function lineStream.
 * @tparam      data_T  class stored in C array
 * @tparam      d       size of std::array
 * @param[in]   str     file stream
 * @param[in]   data    pointer to C array
 * @param[in]   nr      number of elements
 */
template< class data_T >
void absorbASCII( std::fstream &str, data_T *data, int nr ){


    std::vector<data_T>                      temp;

    data_T*                                  itrData, begData, endData;
    typename std::vector<data_T>::iterator   itrTemp, begTemp, endTemp;

    begData = &data[0] ;
    endData = &data[nr-1] ;

    itrData = begData ;

    while( str.good() && itrData!=endData ) {

        temp.clear() ;
        lineStream( str, temp) ;

        begTemp = temp.begin() ;
        endTemp = temp.end() ;

        for( itrTemp=begTemp; (itrTemp!=endTemp && itrData!=endData); ){
            *itrData = *itrTemp ;

            ++itrTemp;
            ++itrData;

        };


    };

    if( itrData != endData ) {
        std::cout << "Not enough elements found to fill array" << std::endl ;
    };


    return ;
};

/*!
 * Reads a templated data type from file stream in binary
 * @tparam      data_T  type of data
 * @param[in]   str     file stream
 * @param[in]   data    data to be read
 */
template< class data_T >
void absorbBINARY( std::fstream &str, data_T &data ){

    int nbytes ;
    nbytes = sizeof(data_T) ;

    str.read( reinterpret_cast<char*>(&data), nbytes ) ;

    return ;
};

/*!
 * Reads a std::vector of templated POD data type from file stream in binary format
 * @tparam      data_T  type of data
 * @param[in]   str     file stream
 * @param[in]   data    data to be read
 */
template< class data_T >
void absorbBINARY( std::fstream &str, std::vector<data_T> &data ){

    int nbytes, nr;
    nr = data.size() ;
    nbytes = sizeof(data_T) *nr ;

    str.read( reinterpret_cast<char*>(&data[0]), nbytes ) ;

    return ;
};

/*!
 * Reads a std::vector of std::vector of templated POD data type from file stream in binary format
 * @tparam      data_T  type of data
 * @param[in]   str     file stream
 * @param[in]   data    data to be read
 */
template< class data_T >
void absorbBINARY( std::fstream &str, std::vector< std::vector<data_T> > &data ){

    for( auto &item: data ){
        absorbBINARY( str, item ) ;
    };

    return ;
};

/*!
 * Reads a std::vector of std::array of templated POD data type of templated size from file stream in binary
 * @tparam      data_T  type of data
 * @tparam      d       size of array
 * @param[in]   str     file stream
 * @param[in]   data    data to be read
 */
template< class data_T, size_t d >
void absorbBINARY( std::fstream &str, std::vector< std::array<data_T,d> > &data ){

    int  nbytes, nr;
    nr = data.size() ;
    nbytes = sizeof(data_T) *nr *d ;

    str.read( reinterpret_cast<char*>(&data[0]), nbytes ) ;

    return ;
};

/*!
 * Reads a std::array of templated POD data type of templated size from file stream in binary format
 * @tparam      data_T  type of data
 * @tparam      d       size of array
 * @param[in]   str     file stream
 * @param[in]   data    data to be read
 */
template< class data_T, size_t d >
void absorbBINARY( std::fstream &str, std::array<data_T,d> &data ){

    int nbytes;
    nbytes = sizeof(data_T) *d ;

    str.read( reinterpret_cast<char*>(&data[0]), nbytes ) ;

    return ;
};

/*!
 * Reads a C array of templated data type from file stream in binary format
 * @tparam      data_T  type of data
 * @param[in]   str     file stream
 * @param[in]   data    data to be read
 * @param[in]   nr      number of elements to be read
 */
template< class data_T >
void absorbBINARY( std::fstream &str, data_T *data, int nr ){

    int nbytes;
    nbytes = sizeof(data_T) *nr ;

    str.read( reinterpret_cast<char*>(data), nbytes ) ;

    return ;
};

/*!
 * @}
 */
}

}
