/*
 * ClassLog.hpp
 *
 *  Created on: 3 déc. 2014
 *      Author: marco
 */

#ifndef INCLUDE_CLASS_LOG_HPP_
#define INCLUDE_CLASS_LOG_HPP_

#if NOMPI==0
#include <mpi.h>
#endif
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;

class Class_Log {

	string filename;
#if NOMPI==0
	MPI_Comm comm;
#endif

public:
#if NOMPI==0
	Class_Log(string filename_,MPI_Comm comm_ = MPI_COMM_WORLD);
#else
	Class_Log(string filename_);
#endif
	~Class_Log();

	void writeLog(string msg);

};

#endif /* INCLUDE_CLASS_LOG_HPP_ */
