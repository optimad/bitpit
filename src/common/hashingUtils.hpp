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

#ifndef __BITPIT_HASHING_UTILS_HPP__
#define __BITPIT_HASHING_UTILS_HPP__

/*! \file */

#include <functional>
#include <tuple>

namespace bitpit {

namespace utils {

/*!
	\ingroup common::utils::hashing

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

}

#endif
