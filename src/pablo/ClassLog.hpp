/*
 * ClassLog.hpp
 *
 *  Created on: 3 dec 2014
 *      Author: marco
 */

#ifndef INCLUDE_CLASSLOG_HPP_
#define INCLUDE_CLASSLOG_HPP_

#if ENABLE_MPI
#include <mpi.h>
#endif

class ClassLog {

	std::string m_filename;

#if ENABLE_MPI
	MPI_Comm m_comm;
#endif

public:
#if ENABLE_MPI
	ClassLog(std::string filename_,MPI_Comm comm_ = MPI_COMM_WORLD);
#else
	ClassLog(std::string filename_);
#endif
	~ClassLog();

	void writeLog(std::string msg);

};

#endif /* INCLUDE_CLASS_LOG_HPP_ */
