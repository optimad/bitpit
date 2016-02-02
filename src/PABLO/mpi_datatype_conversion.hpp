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

#if ENABLE_MPI==1
#ifndef __BITPIT_MPI_DATATYPE_CONVERSION_HPP__
#define __BITPIT_MPI_DATATYPE_CONVERSION_HPP__

#include <mpi.h>

namespace bitpit {

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

}

#endif /* __BITPIT_MPI_DATATYPE_CONVERSION_HPP__ */
#endif /* NOMPI */
