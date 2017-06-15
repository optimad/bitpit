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

#ifndef __BITPIT_STRING_UTILS_HPP__
#define __BITPIT_STRING_UTILS_HPP__

/*! \file */

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace bitpit {

namespace utils {

/*!
    \ingroup common::utils::string

    \brief Namespace for string utility functions
*/
namespace string {

// Trimming operators --------------------------------------------------------------- //
inline std::string &ltrim(                                                     // STRING LEFT TRIMMING
        std::string                             &                                     // (input) std::string to be trimmed
        );
inline std::string &rtrim(                                                     // STRING RIGHT TRIMMING
        std::string                             &                                     // (input) std::string to be trimmed
        );
inline std::string &trim(                                                      // STRING TRIMMING
        std::string                             &                                     // (input) std::string to be trimmed
        );

// Padding operators ---------------------------------------------------------------- //
inline std::string lfill(                                                     // Left filler for input string
        const int                               &,                            // (input) Final string length
        std::string                             &,                            // (input) input string
        char                                                                  // (input) char used as filler
);
inline std::string rfill(                                                     // Right filler for input string
        const int                               &,                            // (input) Final string length
        std::string                             &,                            // (input) input string
        char                                                                  // (input) char used as filler
);
inline std::string zeroPadNumber(                                              // PERFORMS CONVERSION OF INTEGER INTO STRING
        int                                      ,                                    // (input) number of char in std::string
        int                                                                           // (input) integer to be padded
        );

// Input stream operator ------------------------------------------------------------ //
bool getAfterKeyword(                                                               // EXTRACT FIELD AFTER SPECIFIC KEYWORD
        std::string                              ,                                    // (input) std::string
        std::string                              ,                                    // (input) keyword
        char                                     ,                                    // (input) field delimiter
        std::string                             &                                     // (input/output) field found
        );

// returns true if key_ is present in line ------------------------------------------ //
inline bool keywordInString(                                                 // SEARCH KEYWORD IN STRING
        std::string                              ,                                    // (input) input string            
        std::string                                                                   // (input) keyword
        ) ;

// converts a string to fundamental data types and vectors or arrays of them -------- //
template <class T>
void convertString(                                                                  // EXTRACT SCALAR FROM STRING
        std::string                              ,                                    // (input) input string
        T                                       &                                     // (input/output) scalar
        );

template <class T>
void  convertString(                                                                 // EXTRACT DATA FROM STRING AND STORE THEM INTO VECTOR
        std::string                              ,                                    // (input) string
        std::vector<T>                          &                                     // (input/output) vector used to store string
        );

template <class T, size_t n>
void  convertString(                                                                 // EXTRACT DATA FROM STRING AND STORE THEM INTO ARRAY
        std::string                              ,                                    // (input) string
        std::array<T,n>                         &                                     // (input/output) array used to store data
        );

}

}

}

// Template implementation
#include "stringUtils.tpp"

#endif
