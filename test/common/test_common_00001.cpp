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

// ================================================================================== //
//                         OPERATORS - EXAMPLES OF USAGE -                            //
//                                                                                    //
// Examples of Operators usage.                                                       //
// ================================================================================== //
// INFO                                                                               //
// ================================================================================== //
// Author     : Alessandro Alaia                                                      //
// Version    : v3.0                                                                  //
//                                                                                    //
// All rights reserved.                                                               //
// ================================================================================== //

// ================================================================================== //
// INCLUDES                                                                           //
// ================================================================================== //
# include <cmath>
# include <array>
# include <vector>
# include <iostream>

// CC_Lib
# include "bitpit_common.hpp"

// ================================================================================== //
// NAMESPACES                                                                         //
// ================================================================================== //
using namespace std;


// ================================================================================== //
// IMPLEMENTATIONS                                                                    //
// ================================================================================== //

// ---------------------------------------------------------------------------------- //
int main(
    void
) {

// ================================================================================== //
// OPERATOR "ltrim"                                                                   //
// ================================================================================== //
{

    // Output message --------------------------------------------------------------- //
    cout << "** function ltrim()" << endl;

    // Scope variables -------------------------------------------------------------- //
    string                       s = "  test string #1  ";

    // Initialize scope variables --------------------------------------------------- //
    cout << "  std::string s =\"" << s << "\"" << endl;

    // left trimming ---------------------------------------------------------------- //
    cout << "  ltrim(s) = \"" << bitpit::utils::string::ltrim(s) << "\"" << endl;
    cout << endl;
}

// ================================================================================== //
// OPERATOR "rtrim"                                                                   //
// ================================================================================== //
{
    // Output message --------------------------------------------------------------- //
    cout << "** function rtrim()" << endl;

    // Scope variables -------------------------------------------------------------- //
    string                       s = "  test string #2  ";

    // Initialize scope variables --------------------------------------------------- //
    cout << "  std::string s =\"" << s << "\"" << endl;

    // left trimming ---------------------------------------------------------------- //
    cout << "  rtrim(s) = \"" << bitpit::utils::string::rtrim(s) << "\"" << endl;
    cout << endl;
}

// ================================================================================== //
// OPERATOR "trim"                                                                    //
// ================================================================================== //
{
    // Output message --------------------------------------------------------------- //
    cout << "** function trim()" << endl;

    // Scope variables -------------------------------------------------------------- //
    string                       s = "  test string #3  ";

    // Initialize scope variables --------------------------------------------------- //
    cout << "  std::string s =\"" << s << "\"" << endl;

    // left trimming ---------------------------------------------------------------- //
    cout << "  trim(s) = \"" << bitpit::utils::string::trim(s) << "\"" << endl;
    cout << endl;
}

// ================================================================================== //
// OPERATOR "lfill"                                                                   //
// ================================================================================== //
{

    // Output message --------------------------------------------------------------- //
    cout << "** function lfill()" << endl;

    // Scope variables -------------------------------------------------------------- //
    string                       s = "  test string #4  ";

    // Initialize scope variables --------------------------------------------------- //
    cout << "  std::string s =\"" << s << "\"" << endl;

    // left trimming ---------------------------------------------------------------- //
    cout << "  lfill(20, s, '*') = \"" << bitpit::utils::string::lfill(20, s, '*') << "\"" << endl;
    cout << endl;
}

// ================================================================================== //
// OPERATOR "rfill"                                                                   //
// ================================================================================== //
{

    // Output message --------------------------------------------------------------- //
    cout << "** function rfill()" << endl;

    // Scope variables -------------------------------------------------------------- //
    string                       s = "  test string #5  ";

    // Initialize scope variables --------------------------------------------------- //
    cout << "  std::string s =\"" << s << "\"" << endl;

    // left trimming ---------------------------------------------------------------- //
    cout << "  lfill(20, s, '*') = \"" << bitpit::utils::string::rfill(20, s, '*') << "\"" << endl;
    cout << endl;
}

// ================================================================================== //
// OPERATOR "zeroPadNumber"                                                           //
// ================================================================================== //
{
    // Output message --------------------------------------------------------------- //
    cout << "** function zeroPadNumber()" << endl;

    // Scope variables -------------------------------------------------------------- //
    int                          x = 14, y = 1283;
    string                       s;

    // Padding strings -------------------------------------------------------------- //
    cout << "  zeroPadNumber(6, 14) =\"" << bitpit::utils::string::zeroPadNumber(6, x) << "\"" << endl;
    cout << "  zeroPadNumber(3, 1283) =\"" << bitpit::utils::string::zeroPadNumber(3, y) << "\"" << endl;
    cout << endl;
}

// ================================================================================== //
// OPERATOR "ConvertString"                                                           //
// ================================================================================== //
{
    // Output message --------------------------------------------------------------- //
    cout << "** function convertString()" << endl;

    // Scope variables -------------------------------------------------------------- //
    double                       x;
    std::vector<double>          v;
    array<double, 3>             a;
    string                       s0 = " 0.12 ";
    string                       s1 = "  0.12 0.13 0.14 ";
    string                       s2 = "  0.22 0.23 0.24 0.25";

    // Padding strings -------------------------------------------------------------- //
    cout << "  std::string s0 = \"" << s0 << "\"" << endl;
    cout << "  std::string s1 = \"" << s1 << "\"" << endl;
    cout << "  std::string s2 = \"" << s2 << "\"" << endl;
    bitpit::utils::string::convertString(s0, x);
    cout << "  convertString(s0, x), x = " << x << endl;
    bitpit::utils::string::convertString(s1, v);
    cout << "  convertString(s1, v), v = ";
	for (size_t i = 0; i < v.size(); ++i) {
		cout << v[i] << "   ";
	}
    cout << endl;
    bitpit::utils::string::convertString(s2, a);
    cout << "  convertString(s2, a), a = ";
	for (size_t i = 0; i < a.size(); ++i) {
		cout << a[i] << "   ";
	}
    cout << endl;
}

// ================================================================================== //
// OPERATOR "getAfterKeyword"                                                       //
// ================================================================================== //
{
    // Output message --------------------------------------------------------------- //
    cout << "** function getAfterKeyword()" << endl;

    // Scope variables -------------------------------------------------------------- //
    string                       s = " field_1 ; field_2 ; field_3";
    string                       field;

    // Padding strings -------------------------------------------------------------- //
    cout << "  std::string s = \"" << s << "\"" << endl;
    bitpit::utils::string::getAfterKeyword(s, "field_1", ';', field);
    cout << "  getAfterKeyword(s \"field_1\", ';', field), field = \"" << field << "\"" << endl;
    cout << endl;
}

// ================================================================================== //
// OUTPUT MESSAGE                                                                     //
// ================================================================================== //
{
    // Scope variables -------------------------------------------------------------- //
    // none

    // Output message --------------------------------------------------------------- //
    cout << endl;
}

return 0;

};

