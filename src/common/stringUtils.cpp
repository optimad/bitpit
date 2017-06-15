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

// ---------------------------------------------------------------------------------- //
/*!
    Given an input string containing several fields separated by a delimiter,
    returns the field after the specified search key.
    For instance, if the input string is str = "field_1 ; field_2 ; field_3"
    getAfterKeyword(str, "field_1, ';', output) returns
    output = "field_2"

    \param[in] line_ input string
    \param[in] key_ search key
    \param[in] del_ delimiter char
    \param[in,out] result_ on output stores the field after the search key

    \result boolean flag, (true) if search key has been found, (false) otherwise
*/
bool getAfterKeyword( string line_, string key_, char del_, string& result_){

  std::size_t c1, c2, pos ;

  string::iterator it;

  result_.clear() ;

  pos = line_.find( key_);

  if( pos == string::npos ){
      return false;
  }

  else{

    it= line_.begin();
    advance(it, pos);
    advance(it, key_.size() ) ;

    while( (*it) != del_){
      it++;
    };
    c1= it- line_.begin() +1;

    it++;

    while( (*it) != del_){
      it++;
    };
    c2= it- line_.begin()-1;

    pos= c2 -c1 +1;

    result_= line_.substr( c1, pos) ;
    trim( result_ ) ;

    return true;
  };

};

}

}
