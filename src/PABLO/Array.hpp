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

#ifndef __BITPIT_ARRAY_HPP__
#define __BITPIT_ARRAY_HPP__

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include <cstdint>

namespace bitpit {

// =================================================================================== //
// CLASS DEFINITION                                                                    //
// =================================================================================== //

/*! @ingroup PABLO 
 * @{ 
 *
 *	\date		03/mar/2014
 *  \authors	Marco Cisternino
 *  \authors	Edoardo Lombardi
 *	\copyright	Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par		License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *  \brief Customized array definition
 *
 *  Array contains a pointer to an array of integer values.
 *  Implemented here for fast using in PABLO.
 *
 */
class Array {

	friend class ParaTree;

	// ------------------------------------------------------------------------------- //
	// MEMBERS ----------------------------------------------------------------------- //
private:
	uint32_t m_arraySize;		/**< Size of array */
	int* m_array;				/**< Pointer to array of integers */

	// ------------------------------------------------------------------------------- //
	// CONSTRUCTORS ------------------------------------------------------------------ //
public:
	Array();
	Array(const Array& other);
	~Array();
private:
	Array(uint32_t size, int value);

	// ------------------------------------------------------------------------------- //
	// METHODS ----------------------------------------------------------------------- //
	Array& operator=(const Array& rhs);
};

/* @} */

}

#endif /* __BITPIT_ARRAY_HPP__ */
