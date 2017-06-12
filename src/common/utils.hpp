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

#ifndef __BITPIT_UTILS_HPP__
#define __BITPIT_UTILS_HPP__

/*! \file */

#include <array>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

// Stringification macro
#define BITPIT_STR2(X) #X
#define BITPIT_STR(X) BITPIT_STR2(X)

// Macros to allow using oveload in preprocessing macro
#define BITPIT_CAT(A, B) A ## B

#define BITPIT_COUNT_ARGS_MAX6(_1, _2, _3, _4, _5, _6 /* ad nauseam */, COUNT, ...) COUNT
#define BITPIT_EXPAND_ARGS_FOR_COUNT(ARGS) BITPIT_COUNT_ARGS_MAX6 ARGS
#define BITPIT_ARGS_SIZE(...) BITPIT_EXPAND_ARGS_FOR_COUNT((__VA_ARGS__, 6, 5, 4, 3, 2, 1, 0))

#define BITPIT_SELECT_OVERLOAD(NAME, NUM) BITPIT_CAT(NAME ## _, NUM)
#define BITPIT_OVERLOAD_CALL(NAME, ...) BITPIT_SELECT_OVERLOAD(NAME, BITPIT_ARGS_SIZE(__VA_ARGS__))(__VA_ARGS__)

namespace bitpit {

/*!
    \ingroup common::utils

    \brief Namespace for generic utility functions
*/
namespace utils {

template <typename T, typename Comparator = std::less<T> >
bool addToOrderedVector(const T &value, std::vector<T> &list, Comparator comparator = Comparator());

template <typename T, typename Comparator = std::less<T> >
typename std::vector<T>::const_iterator findInOrderedVector(const T &value, const std::vector<T> &list, Comparator comparator = Comparator());

template<typename T>
void reorderVector(std::vector<size_t>& order, std::vector<T>& v, const size_t &size);

template <class T>
void eraseValue(std::vector<T> &, const T&);

template <class T>
std::vector<T> intersectionVector(const std::vector<T>&, const std::vector<T>&);

#ifndef __BITPIT_UTILS_SRC__
extern template bool addToOrderedVector<>(const long&, std::vector<long>&, std::less<long>);
extern template bool addToOrderedVector<>(const unsigned long&, std::vector<unsigned long>&, std::less<unsigned long>);

extern template std::vector<long>::const_iterator findInOrderedVector<>(const long&, const std::vector<long>&, std::less<long>);
extern template std::vector<unsigned long>::const_iterator findInOrderedVector<>(const unsigned long&, const std::vector<unsigned long>&, std::less<unsigned long>);
#endif

void extractWithoutReplacement(                                               // Extract integers without replacement
    int                         ,                                             // (input) number of integers to be extracted
    int                         ,                                             // (input) upper bound of extraction interval
    std::vector<int>           &                                              // (input/output) list of extracted value
);

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

/*!
	Functor to compare two double precision floating point numbers.

	See: https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
*/
struct DoubleFloatingEqual
{
	/*!
		Compares the specified double precision floating point numbers
		and returns true if the numbers match.

		\param x if the first value to compare
		\param y if the second value to compare
		\result Returns true if the numbers match, false otherwise.
	*/
	bool operator()(const double &x, const double &y) const
	{
		const double ABS_MAX_DIFF = 1e-14;
		const double REL_MAX_DIFF = DBL_EPSILON;

		// Check if the numbers are really close (needed when comparing
		// numbers near zero).
		double diff = std::abs(x - y);
		if (diff <= ABS_MAX_DIFF) {
			return true;
		}

		// Check if the numbers have the same sign
		if ((x < 0 && y > 0) || (x > 0 && y < 0)) {
			return false;
		}

		// Compare using a relative difference
		double abs_x   = std::abs(x);
		double abs_y   = std::abs(y);
		double largest = (abs_y > abs_x) ? abs_y : abs_x;

		return (diff <= largest * REL_MAX_DIFF);
	}
};

}

/*!
	\ingroup common::utils

	\brief Functions for generating the hash of data types.
*/
namespace hashing {

namespace
{
	// Code from boost
	// Reciprocal of the golden ratio helps spread entropy
	//     and handles duplicates.
	// See Mike Seymour in magic-numbers-in-boosthash-combine:
	//     http://stackoverflow.com/questions/4948780

	template <class T>
	inline void hash_combine(std::size_t& seed, T const& v)
	{
		seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
	}

	// Recursive template code derived from Matthieu M.
	template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
	struct HashValueImpl
	{
		static void apply(size_t & seed, Tuple const & tuple)
		{
			HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
			hash_combine(seed, std::get<Index>(tuple));
		}
	};

	template <class Tuple>
	struct HashValueImpl<Tuple,0>
	{
		static void apply(size_t & seed, Tuple const & tuple)
		{
			hash_combine(seed, std::get<0>(tuple));
		}
	};
}

template <typename TT>
struct hash
{
    size_t
    operator()(TT const& tt) const
    {
        return std::hash<TT>()(tt);
    }
};

template <typename ... TT>
struct hash<std::tuple<TT...>>
{
	size_t
	operator()(std::tuple<TT...> const & tt) const
	{
		size_t seed = 0;
		HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
		return seed;
	}
};

}

}

# include "utilsString.tpp"
#endif
