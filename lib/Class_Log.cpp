/*
 * ClassLog.cpp
 *
 *  Created on: 3 déc. 2014
 *      Author: marco
 */

#include <Class_Log.hpp>

Class_Log::Class_Log(string filename_) : filename(filename_) {};

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

	int rank;
	int error_flag = MPI_Comm_rank(MPI_COMM_WORLD,&rank);
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
	error_flag = MPI_Barrier(MPI_COMM_WORLD);
	return; };
