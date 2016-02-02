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

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "Array.hpp"
#include <algorithm>

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;

// =================================================================================== //
// CLASS IMPLEMENTATION                                                                    //
// =================================================================================== //

Array::Array() {
	// TODO Auto-generated constructor stub
	m_arraySize = 0;
	m_array = NULL;
}

Array::Array(uint32_t size, int value) {

	m_arraySize = size;
	m_array = new int [size];
	for(uint32_t i = 0; i < size; ++i)
		m_array[i] = value;
}

Array::Array(
		const Array& other) {
	m_arraySize = other.m_arraySize;
	m_array = new int [m_arraySize];
	for(uint32_t i = 0; i < m_arraySize; ++i)
		m_array[i] = other.m_array[i];
}

Array::~Array() {
}

Array& Array::operator =(
		const Array& rhs) {
	if(this != &rhs)
	{
		int* new_array = new int[rhs.m_arraySize];
		copy(rhs.m_array,rhs.m_array+rhs.m_arraySize,new_array);

		delete [] m_array;

		m_array = new_array;
		m_arraySize = rhs.m_arraySize;
	}
	return *this;
}
