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

#include "stringUtils.hpp"

using namespace std;

namespace bitpit {

namespace utils {

namespace string {

/*!
* Given an input string containing several fields separated by a delimiter,
* returns the field after the specified search key.
* For instance, if the input string is str = "field_1 ; field_2 ; field_3"
* getAfterKeyword(str, "field_1, ';', output) returns
* output = "field_2"
*
* \param[in] line input string
* \param[in] key search key
* \param[in] del delimiter char
* \param[in,out] result on output stores the field after the search key
*
* \result boolean flag, (true) if search key has been found, (false) otherwise
*/
bool getAfterKeyword(std::string line, std::string key, char del, std::string &result)
{
    result.clear();

    std::size_t  pos = line.find(key);
    if (pos == std::string::npos) {
        return false;
    }

    std::string::iterator it = line.begin();
    advance(it, pos);
    advance(it, key.size());

    while ((*it) != del) {
        ++it;
    }
    std::size_t c1= it- line.begin() + 1;

    ++it;

    while ((*it) != del) {
        ++it;
    }
    std::size_t c2 = it - line.begin() - 1;

    pos= c2 - c1 + 1;

    result = line.substr(c1, pos);
    trim(result);

    return true;
}

}

}

}
