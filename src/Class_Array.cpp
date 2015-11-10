#include "Class_Array.hpp"

Class_Array::Class_Array() {
	// TODO Auto-generated constructor stub
	arraySize = 0;
	array = NULL;
}

Class_Array::Class_Array(uint32_t size, int value) {

	arraySize = size;
	array = new int [size];
	for(uint32_t i = 0; i < size; ++i)
		array[i] = value;
}

Class_Array::Class_Array(
		const Class_Array& other) {
	arraySize = other.arraySize;
	//	if(array != NULL){
	//		delete [] array;
	//		array = NULL;
	//	}
	array = new int [arraySize];
	for(uint32_t i = 0; i < arraySize; ++i)
		array[i] = other.array[i];
}

Class_Array::~Class_Array() {
	//	delete [] array;
	//	array = NULL;
}

Class_Array& Class_Array::operator =(
		const Class_Array& rhs) {
	if(this != &rhs)
	{
		int* new_array = new int[rhs.arraySize];
		copy(rhs.array,rhs.array+rhs.arraySize,new_array);

		delete [] array;

		array = new_array;
		arraySize = rhs.arraySize;
	}
	return *this;
}
