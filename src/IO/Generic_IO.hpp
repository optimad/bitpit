
#ifndef __Generic_IO__HH__
#define __Generic_IO__HH__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>

/*

line_stream        reads one line from stream str and check if input is sufficient to fill data. 
                   if successful data will contain input values;
                   exeption:  if data contains a zero-size vector, data will be resized to match the entire line;

flush_(format)     writes the content of data into stream str.
                   if format == ascii, elements_per_line need to be specified


absorb_(format)    reads the filestream str and fills the content of data.
                   checks if data has been filled completely 

*/

template< class data_T >
void  line_stream( std::fstream &str, data_T &data) ;

template< class data_T >
void  line_stream( std::fstream &str, std::vector<data_T> &data) ;

template< class data_T, size_t d >
void  line_stream( std::fstream &str, std::array<data_T,d> &data) ;

template< class data_T, size_t d >
void  line_stream( std::fstream &str, int nr, data_T *data) ;

template< class data_T >
void flush_ascii( std::fstream &str, const data_T &data  ) ;

template<>
void flush_ascii( std::fstream &str, const uint8_t &data  ) ;

template< class data_T >
void flush_ascii( std::fstream &str, int elements_per_line, const std::vector<data_T> &data  ) ;

template< class data_T, size_t d >
void flush_ascii( std::fstream &str, int elements_per_line, const std::array<data_T,d> &data  ) ;

template< class data_T >
void flush_ascii( std::fstream &str, int elements_per_line, const data_T *data, int nr  ) ;

template< class data_T >
void flush_binary( std::fstream &str, const data_T &data  ) ;

template< class data_T >
void flush_binary( std::fstream &str, const std::vector<data_T> &data  ) ;

template< class data_T >
void flush_binary( std::fstream &str, const std::vector< std::vector<data_T> > &data  ) ;

template< class data_T, size_t d >
void flush_binary( std::fstream &str, const std::vector< std::array<data_T,d> > &data  ) ;

template< class data_T, size_t d >
void flush_binary( std::fstream &str, const std::array<data_T,d> &data  ) ;

template< class data_T >
void flush_binary( std::fstream &str, const data_T *data, int nr  ) ;

template< class data_T >
void absorb_ascii( std::fstream &str, data_T &data  ) ;

template< class data_T >
void absorb_ascii( std::fstream &str, std::vector<data_T> &data  ) ;

template< class data_T, size_t d >
void absorb_ascii( std::fstream &str, std::array<data_T,d> &data  ) ;

template< class data_T >
void absorb_ascii( std::fstream &str, data_T *data, int nr  ) ;

template< class data_T >
void absorb_binary( std::fstream &str, data_T &data  ) ;

template< class data_T >
void absorb_binary( std::fstream &str, std::vector<data_T> &data  ) ;

template< class data_T >
void absorb_binary( std::fstream &str, std::vector< std::vector<data_T> > &data  ) ;

template< class data_T, size_t d >
void absorb_binary( std::fstream &str, std::vector< std::array<data_T,d> > &data  ) ;

template< class data_T, size_t d >
void absorb_binary( std::fstream &str, std::array<data_T,d> &data  ) ;

template< class data_T >
void absorb_binary( std::fstream &str, data_T *data, int nr  ) ;

void CopyUntilEOFInString( std::fstream &str, char*& buffer, int& length);

#include "Generic_IO.tpp"

#endif
