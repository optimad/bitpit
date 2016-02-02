/*
 * Log.hpp
 *
 *  Created on: 3 dec 2014
 *      Author: marco
 */

#ifndef LOG_HPP_
#define LOG_HPP_

#if ENABLE_MPI==1
#include <mpi.h>
#endif
#include <string>

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

#endif /* LOG_HPP_ */
