/*
 * ClassLog.hpp
 *
 *  Created on: 3 déc. 2014
 *      Author: marco
 */

#ifndef INCLUDE_CLASS_LOG_HPP_
#define INCLUDE_CLASS_LOG_HPP_

# include <vector>
# include <string>
# include <fstream>
# include <sstream>
#include <mpi.h>

using namespace std;

class Class_Log {

	string filename;

public:
	Class_Log(string filename_);
	~Class_Log();

	void writeLog(string msg);

};

#endif /* INCLUDE_CLASS_LOG_HPP_ */
