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
