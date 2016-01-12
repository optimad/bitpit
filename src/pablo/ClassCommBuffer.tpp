// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //

#include "mpi_datatype_conversion.hpp"

// =================================================================================== //
// CLASS IMPLEMENTATION                                                                    //
// =================================================================================== //

template<class T>
void ClassCommBuffer::write(T& val) {
	MPI_Datatype datatype = convert<T>();
	int error = MPI_Pack(&val,1,datatype,m_commBuffer,m_commBufferSize,&m_pos,m_comm);
};

template<class T>
void ClassCommBuffer::read(T& val) {
	MPI_Datatype datatype = convert<T>();
	int error = MPI_Unpack(m_commBuffer,m_commBufferSize,&m_pos,&val,1,datatype,m_comm);
};

