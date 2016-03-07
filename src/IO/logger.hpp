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
#ifndef __BITPIT_LOGGER_HPP__
#define __BITPIT_LOGGER_HPP__

#include <memory>
#include <unordered_map>
#include <vector>

#include "utils.hpp"
#include "FileHandler.hpp"

#define BITPIT_DEBUG_COUT(...) BITPIT_OVERLOAD_CALL(BITPIT_DEBUG_COUT, __VA_ARGS__)
#if ENABLE_DEBUG==1
#define BITPIT_DEBUG_COUT_0()     bitpit::log::cout()
#define BITPIT_DEBUG_COUT_1(NAME) bitpit::log::cout(NAME)
#else
#define BITPIT_DEBUG_COUT_0()     while(0) bitpit::log::cout()
#define BITPIT_DEBUG_COUT_1(NAME) while(0) bitpit::log::cout(NAME)
#endif

namespace bitpit{

namespace log {
	enum Verbosity {
	    QUIET = 0,
	    NORMAL,
	    DEBUG
	};

	enum Visibility {
	    MASTER = 0,
	    GLOBAL
	};

	typedef Verbosity Priority;
}

// Logger buffer
class LoggerBuffer : public std::streambuf
{

public:
	LoggerBuffer(std::size_t bufferSize = 256);
	~LoggerBuffer();

	void setConsoleEnabled(bool enabled);
	void setConsoleStream(std::ostream *console);
	std::ostream & getConsoleStream();
	void setConsolePrefix(const std::string &prefix);

	void setFileEnabled(bool enabled);
	void setFileStream(std::ofstream *file);
	std::ofstream & getFileStream();
	void setFilePrefix(const std::string &prefix);

	void setContext(const std::string &context);
	void setPadding(const std::string &padding);

    int flush(bool terminate);

private:
	std::vector<char> m_buffer;

	std::string m_context;
	std::string m_padding;

	bool m_consoleEnabled;
	std::ostream *m_console;
	std::string m_consolePrefix;

	bool m_fileEnabled;
	std::ofstream *m_file;
	std::string m_filePrefix;

	const std::string getTimestamp() const;

    int_type overflow(int_type ch);
    int sync();

};

// Logger
class Logger : public std::ostream
{

public:
	Logger(std::string name);
	~Logger();

	static Logger & cout(std::string name = "bitpit");

	void setParallel(int nProcessors, int rank);

	void setContext(const std::string &context);
	std::string getContext();

	void setIndentation(const int &delta);
	int getIndentation();

	void setPriority(log::Priority priority);
	log::Priority getPriority();

	void setVisibility(log::Visibility visibility);
	log::Visibility getVisibility();

	void setConsoleVerbosity(log::Verbosity verbosity);
	log::Verbosity getConsoleVerbosity();

	void setFileVerbosity(log::Verbosity verbosity);
	log::Verbosity getFileVerbosity();
	void resetLogFile();

	std::string getName();

	void println(const std::string &line);
	void println(const std::string &line, log::Priority priority);
	void println(const std::string &line, log::Visibility visibility);
	void println(const std::string &line, const log::Priority priority, log::Visibility visibility);

	void print(const std::string &line);
	void print(const std::string &line, log::Priority priority);
	void print(const std::string &line, log::Visibility visibility);
	void print(const std::string &line, log::Priority priority, log::Visibility visibility);

private:
	typedef std::unordered_map<std::string, std::unique_ptr<Logger>> Cache;
	typedef std::unordered_map<std::string, std::unique_ptr<Logger>>::iterator CacheIterator;
	static Cache m_loggers;

	bool m_isMaster;
	LoggerBuffer m_buffer;

	int m_indentation;
	std::string m_context;
	log::Priority m_priority;
	log::Visibility m_visibility;

	log::Verbosity m_consoleVerbosity;

	FileHandler m_fileHandler;
	std::ofstream m_fileStream;
	log::Verbosity m_fileVerbosity;

	Logger(Logger const&) = delete;
	Logger& operator=(Logger const&) = delete;

	void openLogFile(bool reset);

};

/*!
	\brief The namespace 'log' contains routines for interacting with the
	message logger.
 */
namespace log {

	// Generic global functions
	Logger & cout(std::string name = "bitpit");

	// Manipulators global functions

	/*!
		Struct that allows to manipulate a logger.
	*/
	template<typename T>
	struct LoggerManipulator {
		Logger& (*f) (Logger&, const T &);
		T value;

		/*!
			Creates a new logger manipulator
		*/
		LoggerManipulator(Logger& (*ff)(Logger&, const T &), const T & ss)
			: f(ff), value(ss) {

		}
	};

	/*!
		Applies a manipulator to a logger.
	*/
	template<typename T>
	Logger& operator<<(Logger& logger, LoggerManipulator<T>&& m)
	{
		return m.f(logger, m.value);
	}

	Logger& setContext(Logger& logger, const std::string &context);
	LoggerManipulator<std::string> context(const std::string &context);

	Logger& setPriority(Logger& logger, const log::Priority &priority);
	LoggerManipulator<log::Priority> priority(const log::Priority &priority);

	Logger& setVisibility(Logger& logger, const log::Visibility &visibility);
	LoggerManipulator<log::Visibility> visibility(const log::Visibility &visibility);

	Logger& setConsoleVerbosity(Logger& logger, const log::Verbosity &verbosity);
	LoggerManipulator<log::Priority> consoleVerbosity(const log::Verbosity &verbosity);

	Logger& setFileVerbosity(Logger& logger, const log::Verbosity &verbosity);
	LoggerManipulator<log::Priority> fileVerbosity(const log::Verbosity &verbosity);

	Logger& setIndentation(Logger& logger, const int &delta);
	LoggerManipulator<int> indent(const int &delta);

}

}

#endif
