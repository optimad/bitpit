#if ENABLE_MPI==1
#ifndef MPI_DATATYPE_CONVERSION_HPP_
#define MPI_DATATYPE_CONVERSION_HPP_

#include <mpi.h>

template <class T>
MPI_Datatype convert();

template <>
inline MPI_Datatype convert<char>(){return MPI::CHAR;}

template <>
inline MPI_Datatype convert<int>(){return MPI::INT;}

template <>
inline MPI_Datatype convert<short>(){return MPI::SHORT;}

template <>
inline MPI_Datatype convert<long>(){return MPI::LONG;}

template <>
inline MPI_Datatype convert<unsigned int>(){return MPI::UNSIGNED;}

template <>
inline MPI_Datatype convert<unsigned short>(){return MPI::UNSIGNED_SHORT;}

template <>
inline MPI_Datatype convert<unsigned long>(){return MPI::UNSIGNED_LONG;}

template <>
inline MPI_Datatype convert<double>(){return MPI::DOUBLE;}

template <>
inline MPI_Datatype convert<float>(){return MPI::FLOAT;}

template <>
inline MPI_Datatype convert<long double>(){return MPI::LONG_DOUBLE;}

template <>
inline MPI_Datatype convert<bool>(){return MPI::BOOL;}

template <>
inline MPI_Datatype convert<uint8_t>(){return MPI_UINT8_T;}

//template <>
//inline MPI_Datatype convert<uint64_t>(){return MPI_UINT64_T;}

template <>
inline MPI_Datatype convert<unsigned long long>(){return MPI_UINT64_T;}

template <>
inline MPI_Datatype convert<int8_t>(){return MPI_INT8_T;}

#endif /* MPI_DATATYPE_CONVERSION_HPP_ */
#endif /* NOMPI */
