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
 *  as published by by the Free Software Foundation.
 *
 *  BitPit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#ifndef __Generic_IO__HH__
#define __Generic_IO__HH__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>

/*

lineStream        reads one line from stream str and check if input is sufficient to fill data. 
                   if successful data will contain input values;
                   exeption:  if data contains a zero-size vector, data will be resized to match the entire line;

flush_(format)     writes the content of data into stream str.
                   if format == ascii, elements_per_line need to be specified


absorb_(format)    reads the filestream str and fills the content of data.
                   checks if data has been filled completely 

*/

template< class data_T >
void  lineStream( std::fstream &str, data_T &data) ;

template< class data_T >
void  lineStream( std::fstream &str, std::vector<data_T> &data) ;

template< class data_T, size_t d >
void  lineStream( std::fstream &str, std::array<data_T,d> &data) ;

template< class data_T, size_t d >
void  lineStream( std::fstream &str, int nr, data_T *data) ;

template< class data_T >
void flushASCII( std::fstream &str, const data_T &data  ) ;

template<>
void flushASCII( std::fstream &str, const uint8_t &data  ) ;

template< class data_T >
void flushASCII( std::fstream &str, int elements_per_line, const std::vector<data_T> &data  ) ;

template< class data_T, size_t d >
void flushASCII( std::fstream &str, int elements_per_line, const std::array<data_T,d> &data  ) ;

template< class data_T >
void flushASCII( std::fstream &str, int elements_per_line, const data_T *data, int nr  ) ;

template< class data_T >
void flushBINARY( std::fstream &str, const data_T &data  ) ;

template< class data_T >
void flushBINARY( std::fstream &str, const std::vector<data_T> &data  ) ;

template< class data_T >
void flushBINARY( std::fstream &str, const std::vector< std::vector<data_T> > &data  ) ;

template< class data_T, size_t d >
void flushBINARY( std::fstream &str, const std::vector< std::array<data_T,d> > &data  ) ;

template< class data_T, size_t d >
void flushBINARY( std::fstream &str, const std::array<data_T,d> &data  ) ;

template< class data_T >
void flushBINARY( std::fstream &str, const data_T *data, int nr  ) ;

template< class data_T >
void absorbASCII( std::fstream &str, data_T &data  ) ;

template< class data_T >
void absorbASCII( std::fstream &str, std::vector<data_T> &data  ) ;

template< class data_T, size_t d >
void absorbASCII( std::fstream &str, std::array<data_T,d> &data  ) ;

template< class data_T >
void absorbASCII( std::fstream &str, data_T *data, int nr  ) ;

template< class data_T >
void absorbBINARY( std::fstream &str, data_T &data  ) ;

template< class data_T >
void absorbBINARY( std::fstream &str, std::vector<data_T> &data  ) ;

template< class data_T >
void absorbBINARY( std::fstream &str, std::vector< std::vector<data_T> > &data  ) ;

template< class data_T, size_t d >
void absorbBINARY( std::fstream &str, std::vector< std::array<data_T,d> > &data  ) ;

template< class data_T, size_t d >
void absorbBINARY( std::fstream &str, std::array<data_T,d> &data  ) ;

template< class data_T >
void absorbBINARY( std::fstream &str, data_T *data, int nr  ) ;

void copyUntilEOFInString( std::fstream &str, char*& buffer, int& length);

#include "GenericIO.tpp"

#endif
