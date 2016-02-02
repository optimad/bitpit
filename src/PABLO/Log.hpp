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

/*
 * Log.hpp
 *
 *  Created on: 3 dec 2014
 *      Author: marco
 */

#ifndef __BITPIT_PABLO_LOG_HPP__
#define __BITPIT_PABLO_LOG_HPP__

#if ENABLE_MPI==1
#include <mpi.h>
#endif
#include <string>

namespace bitpit {

class Log {

	std::string m_filename;

#if ENABLE_MPI==1
	MPI_Comm m_comm;
#endif

public:
#if ENABLE_MPI==1
	Log(std::string filename_,MPI_Comm comm_ = MPI_COMM_WORLD);
#else
	Log(std::string filename_);
#endif
	~Log();

	void writeLog(std::string msg);

};

}

#endif /* __BITPIT_PABLO_LOG_HPP__ */
