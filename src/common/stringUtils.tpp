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

#ifndef __BITPIT_STRING_UTILS_TPP__
#define __BITPIT_STRING_UTILS_TPP__

namespace bitpit {

namespace utils {

namespace string {

/*!
* Left-trim operator for std::string.
*
* Remove left trailing spaces from string. For instance, if the input string is
* "  test_string  ", on output this function returns "test_string  "
*
* \param[in] s is the input string
* \result A reference to the input string.
*/
inline std::string & ltrim(std::string &s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not1(std::ptr_fun<int, int>(std::isspace))));

    return s;
}

/*!
* Right-trim operator for std::string.
* Remove right blank spaces from string. For instance, if the input string is
* "  test_string  ", on output this function returns "  test_string"
*
* \param[in] s is the input string
* \result A reference to the input string.
*/
inline std::string &rtrim(std::string &s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), not1(std::ptr_fun<int, int>(isspace))).base(), s.end());

    return s;
}

/*!
* Trim operator for std::string.
* Remove left/right blank spaces from string. For instance, if the input string is
* "  test_string  ", on output this function returns "test_string"
*
* \param[in] s is the input string
* \result A reference to the input string.
*/
inline std::string &trim(std::string &s)
{
    return ltrim(rtrim(s));
}

/*!
* String left-filler. Create a string composed of the input string left-filled
* with a specified character. E.g.
* given the input string s = "test", lfill(10, s, '_') will return
* "______test".
*
* \param[in] nchars is the final length of the string
* \param[in] s is the input string
* \param[in] c is the char used as filler
*/
inline std::string lfill(const int &nchars, std::string &s, char c)
{
    std::stringstream ss;
    ss << std::string(nchars - s.length(), c) << s;

    return ss.str();
}

/*!
* String right-filler. Create a string composed of the input string right-filled
* with a specified character. E.g.
* given the input string s = "test", rfill(10, s, '_') will return
* "test______".
*
* \param[in] nchars is the final length of the string
* \param[in] s is the input string
* \param[in] c is the char used as filler
*/
inline std::string rfill(const int &nchars, std::string &s, char c)
{
    std::stringstream ss;
    ss << s << std::string(nchars - s.length(), c);

    return ss.str();
}

/*!
* Given an integer, returns a string of length nchars, composed of the input number
* and nchars - ndigits '0' characters (where ndigits is the number of digits of the input integer)
* in the following format "000xxx".
* If ndigits > nchars, the output string will contaiens ndigits characters storing the
* digits of the input number.
* For instance, if nchars = 4 and num = 3, this function returns the string "0003".
* If nchars = 4, and num = 12345, this function returns "12345".
*
* \param[in] nchars is the final length of the string
* \param[in] num is the input integer
* \result A string storing the input number in the format 000xxx.
*/
inline std::string zeroPadNumber(int nchars, int num)
{
    std::ostringstream ss;
    ss << std::setw(nchars) << std::setfill('0') << num;

    return ss.str();
}

/*!
* Check whether a string contains the kwyword or not.
*
* \param[in] line is the input string
* \param[in] key is the search key
* \result Return true if the keyword has been found, false otherwise.
*/
inline bool keywordInString(std::string line, std::string key)
{
    return (line.find(key) != std::string::npos);
}

/*!
* Convertes a string into fundamental data type.
*
* If no data of type T can be extracted from the input string a 0 value,
* will be stored in output.
* If multiple values can be extracted from the input string, only the first
* value will be saved in output.
*
* \param[in] input is the input string
* \param[out] output on output contains the value extracted from the input
* string
*/
template <class T>
void convertString(std::string input, T &output)
{
    trim(input);
    std::stringstream ss(input);

    T x;
    std::vector<T> tmp;
    while (ss.good()) {
        ss >> x;
        tmp.push_back(x);
    }

    if (tmp.size() == 0) {
        std::cout << " no useful information in string " << input   << std::endl;
        std::cout << " casting zero                   " <<  std::endl;

        x = static_cast<T> (0);
    } else if (tmp.size() == 1) {
        x = tmp[0];
    } else if(tmp.size() > 1) {
        std::cout << " more than one element in string " << input   << std::endl;
        std::cout << " assigning first element             "  << std::endl;

        x = tmp[0];
    }

    output = x;
}

/*!
* Convertes a string into a vector of fundamental data type.
*
* If no data of type T can be extracted from the input string a void vector is returned.
* Values extracted from string are pushed at the end of the vector.
*
* \param[in] input is the input string
* \param[out] output on output contains the values extracted from the input
* string
*/
template <class T>
void convertString(std::string input, std::vector<T> &output)
{
    output.clear();

    trim(input);
    std::stringstream ss(input);

    T x;
    std::vector<T> tmp;
    while (ss.good()) {
        ss >> x;
        tmp.push_back(x);
    }

    if (tmp.size() == 0) {
        std::cout << " no useful information in string " << input   << std::endl;
        std::cout << " returning void vector          " <<  std::endl;
    };

    output = tmp;
}

/*!
* Convertes a string into a arrayof fundamental data type.
*
* If no data of type T can be extracted from the input string a void array with null
* elements is returned.
* If the number of elements which can be extracted from the input string is larger
* than the array size, only the first n elements are saved in the array.
*
* \param[in] input is the input string
* \param[out] output on output contains the values extracted from the input
* string
*/
template <class T, size_t n>
void convertString(std::string input, std::array<T,n> &output)
{
    T x;
    std::vector<T> tmp;

    tmp.clear();

    trim(input);
    std::stringstream ss(input);

    while (ss.good()) {
        ss >> x;
        tmp.push_back(x);
    }

    if (tmp.size() < n) {
        std::cout << " not enough useful information in string " << input   << std::endl;
        std::cout << " casting zero into missing elements      " <<  std::endl;

        x = static_cast<T>(0);
        output.fill(x);

        for(size_t i=0; i<tmp.size(); i++) {
            output[i] = tmp[i];
        }
    } else if (tmp.size() == n) {
        for(size_t i = 0; i < n; i++) {
            output[i] = tmp[i];
        }
    } else if (tmp.size() > n) {
        std::cout << " more than " << n << " elements in string " << input   << std::endl;
        std::cout << " assigning first element " << n << " elements "   << std::endl;

        for(size_t i = 0; i < n; i++) {
            output[i] = tmp[i];
        }
    }
}

}

}

}

#endif
