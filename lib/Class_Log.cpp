/*
 * ClassLog.cpp
 *
 *  Created on: 3 déc. 2014
 *      Author: marco
 */

#include <Class_Log.hpp>

#if NOMPI==0
Class_Log::Class_Log(string filename_,MPI_Comm comm_) : filename(filename_),comm(comm_) {};
#else
Class_Log::Class_Log(string filename_) : filename(filename_) {};
#endif

Class_Log::~Class_Log() {};

// ----------------------------------------------------------------------------------- //
void Class_Log::writeLog(string msg) {

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
#if NOMPI==0
	int error_flag = MPI_Comm_rank(comm,&rank);
#endif
	if(rank == 0){
		// Open the .log file
		file_handle.open(filename.c_str(), ifstream::app);
		if(!file_handle.is_open())
			exit(1);

		// Append message to the .log file
		file_handle << msg << endl;

		// Close file
		file_handle.close();
	}
#if NOMPI==0
	error_flag = MPI_Barrier(comm);
#endif
	return; };
