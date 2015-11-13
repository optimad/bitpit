#if NOMPI==0
#include "Class_Comm_Buffer.hpp"

Class_Comm_Buffer::Class_Comm_Buffer(){

	commBufferSize = 0;
	commBuffer = NULL;
	pos = 0;
	comm = MPI_COMM_WORLD;
}

Class_Comm_Buffer::Class_Comm_Buffer(MPI_Comm comm_) : comm(comm_){

	commBufferSize = 0;
	commBuffer = NULL;
	pos = 0;
}

Class_Comm_Buffer::Class_Comm_Buffer(uint32_t size, char value, MPI_Comm comm_) : comm(comm_){

	commBufferSize = size;
	commBuffer = new char [size];
	for(uint32_t i = 0; i < size; ++i)
		commBuffer[i] = value;
	pos = 0;
}

Class_Comm_Buffer::Class_Comm_Buffer(const Class_Comm_Buffer& other) {

	commBufferSize = other.commBufferSize;
	//	if(commBuffer != NULL){
	//		delete [] commBuffer;
	//		commBuffer = NULL;
	//	}
	commBuffer = new char [commBufferSize];
	for(uint32_t i = 0; i < commBufferSize; ++i)
		commBuffer[i] = other.commBuffer[i];
	pos = other.pos;
	comm = other.comm;
}

Class_Comm_Buffer::~Class_Comm_Buffer() {
		delete [] commBuffer;
		commBuffer = NULL;
}

Class_Comm_Buffer& Class_Comm_Buffer::operator =(const Class_Comm_Buffer& rhs) {
	if(this != &rhs)
	{
		char* new_array = new char[rhs.commBufferSize];
		copy(rhs.commBuffer,rhs.commBuffer+rhs.commBufferSize,new_array);

		delete [] commBuffer;

		commBuffer = new_array;
		commBufferSize = rhs.commBufferSize;
		pos = rhs.pos;
		comm = rhs.comm;
	}
	return *this;

}
#endif
