/*!
 * @ingroup Generic
 * @{
 */

#include "GenericIO.hpp"

/*!
 * Writes uint_8 data as formatted integer 
 * @param[in]   str     file stream to be copied; file needs already to be opened
 * @param[in]   data    data to be written
 */
template< >
void flushASCII( std::fstream &str, const uint8_t &data ){

    str << std::setprecision(8) << std::scientific ;
    str << unsigned(data) << " ";

    return ;
};


/*!
 * Copies entire file into a char array
 * @param[in]   str     file stream to be copied; file nneds already to be opened
 * @param[out]  buffer  char array containing entire array
 * @param[out]  length  number of elements which hvae benn copied
 */
void copyUntilEOFInString( std::fstream &str, char*& buffer, int &length){

    std::fstream::pos_type   position_insert, position_eof ;

    position_insert = str.tellg() ;

    str.seekg(0, str.end);
    position_eof = str.tellg();

    length = position_eof - position_insert ;
    buffer = new char[length] ;

    str.seekg(position_insert);
    str.read(buffer,length ) ;

    str.seekg(position_insert);

    return;
};

/*!
 * @}
 */
