#include "Generic_IO.hpp"

void CopyUntilEOFInString( fstream &str, char*& buffer, int &length){

    fstream::pos_type   position_insert, position_eof ;

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
