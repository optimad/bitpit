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
#include "logFunct.hpp"
#if ENABLE_MPI==1
#include <mpi.h>
#endif
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>

namespace bitpit {

// =================================================================================== //
// NAMESPACES                                                                          //
// =================================================================================== //
using namespace std;

// ----------------------------------------------------------------------------------- //
void writeLog(string msg) {

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
	int error_flag = MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
	if(rank == 0){
		// Open the .log file
		stringstream ss,time;
		time_t now = chrono::system_clock::to_time_t(chrono::system_clock::now());
		time << ctime(&now);
		string timeformat;
		timeformat.append(time.str().substr(0,3));
		timeformat.append("_");
		timeformat.append(time.str().substr(4,3));
		timeformat.append("_");
		timeformat.append(time.str().substr(9,1));
		timeformat.append("_");
		timeformat.append(time.str().substr(11,2));
		timeformat.append("h");
		timeformat.append(time.str().substr(14,2));
		timeformat.append("m");
		timeformat.append(time.str().substr(17,2));
		timeformat.append("s_");
		timeformat.append(time.str().substr(20,4));

		ss << "PABLO_"<< timeformat << ".log";
		//file_handle.open("PABLO.log", ifstream::app);
		file_handle.open(ss.str().c_str(), ofstream::app);
		if(!file_handle.is_open())
			exit(1);

		// Append message to the .log file
		file_handle << msg << endl;

		// Close file
		file_handle.close();
	}
#if ENABLE_MPI==1
	error_flag = MPI_Barrier(MPI_COMM_WORLD);
#endif
	return; };

}
