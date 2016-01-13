#if ENABLE_MPI
// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "ClassCommBuffer.hpp"
#include <cstdint>

// =================================================================================== //
// NAME SPACES                                                                         //
// =================================================================================== //
using namespace std;


ClassCommBuffer::ClassCommBuffer(){

	m_commBufferSize = 0;
	m_commBuffer = NULL;
	m_pos = 0;
	m_comm = MPI_COMM_WORLD;
}

ClassCommBuffer::ClassCommBuffer(MPI_Comm comm_) : m_comm(comm_){

	m_commBufferSize = 0;
	m_commBuffer = NULL;
	m_pos = 0;
}

ClassCommBuffer::ClassCommBuffer(uint32_t size, char value, MPI_Comm comm_) : m_comm(comm_){

	m_commBufferSize = size;
	m_commBuffer = new char [size];
	for(uint32_t i = 0; i < size; ++i)
		m_commBuffer[i] = value;
	m_pos = 0;
}

ClassCommBuffer::ClassCommBuffer(const ClassCommBuffer& other) {

	m_commBufferSize = other.m_commBufferSize;
	//	if(commBuffer != NULL){
	//		delete [] commBuffer;
	//		commBuffer = NULL;
	//	}
	m_commBuffer = new char [m_commBufferSize];
	for(uint32_t i = 0; i < m_commBufferSize; ++i)
		m_commBuffer[i] = other.m_commBuffer[i];
	m_pos = other.m_pos;
	m_comm = other.m_comm;
}

ClassCommBuffer::~ClassCommBuffer() {
		delete [] m_commBuffer;
		m_commBuffer = NULL;
}

ClassCommBuffer& ClassCommBuffer::operator =(const ClassCommBuffer& rhs) {
	if(this != &rhs)
	{
		char* new_array = new char[rhs.m_commBufferSize];
		copy(rhs.m_commBuffer,rhs.m_commBuffer+rhs.m_commBufferSize,new_array);

		delete [] m_commBuffer;

		m_commBuffer = new_array;
		m_commBufferSize = rhs.m_commBufferSize;
		m_pos = rhs.m_pos;
		m_comm = rhs.m_comm;
	}
	return *this;

}
#endif
