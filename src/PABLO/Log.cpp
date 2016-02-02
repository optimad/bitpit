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
 * Log.cpp
 *
 *  Created on: 3 dec 2014
 *      Author: marco
 */

#include "Log.hpp"
#include <fstream>

namespace bitpit {

using namespace std;

#if ENABLE_MPI==1
Log::Log(string filename_,MPI_Comm comm_) : m_filename(filename_),m_comm(comm_) {};
#else
Log::Log(string filename_) : m_filename(filename_) {};
#endif

Log::~Log() {};

// ----------------------------------------------------------------------------------- //
void Log::writeLog(string msg) {

	// =================================================================================== //
	// void Write_Log(string msg)                                                          //
	//                                                                                     //
	// Append message into a .log file                                                     //
	// =================================================================================== //
	// INPUT                                                                               //
	// =================================================================================== //
	// - msg    : string, message to be appended into the .log file.                       //
	// =================================================================================== //
	// OUTPUT                                                                              //
	// =================================================================================== //
	// - none                                                                              //
	// =================================================================================== //

	// =================================================================================== //
	// VARIABLES DECLARATION                                                               //
	// =================================================================================== //

	// Local variables
	ofstream              file_handle;

	// Counters
	// none

	// =================================================================================== //
	// APPEND MESSAGE TO THE LOG FILE                                                      //
	// =================================================================================== //

	int rank = 0;

#if ENABLE_MPI==1
	bool flag = MPI::Is_finalized();
	if (!(flag))
		int error_flag = MPI_Comm_rank(m_comm,&rank);
#endif
	if(rank == 0){
		// Open the .log file
		file_handle.open(m_filename.c_str(), ifstream::app);
		if(!file_handle.is_open())
			exit(1);

		// Append message to the .log file
		file_handle << msg << endl;

		// Close file
		file_handle.close();
	}
//#if ENABLE_MPI==1
//	error_flag = MPI_Barrier(comm);
//#endif
	return; };

}
