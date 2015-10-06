#include "Generic_IO.hpp"

// =================================================================================== //
template< >
void flush_ascii( fstream &str, const uint8_t data ){

    str << setprecision(8) << scientific ;
    str << unsigned(data) << " ";

    return ;
};


// =================================================================================== //
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
