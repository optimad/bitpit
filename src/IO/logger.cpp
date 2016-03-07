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

#include <cassert>
#include <cctype>
#include <cmath>
#include <chrono>
#include <ctime>
#include <ostream>
#include <functional>
#include <sys/types.h>
#include <sys/stat.h>

#include "logger.hpp"

namespace bitpit{

/*!
	@ingroup Logger
	@{
*/

/*!
    \class LoggerBuffer

    \brief Stream buffer for the message logger.

    This class implements a stream buffer used by the message logger.
*/

/*!
	Builds a new buffer.
*/
LoggerBuffer::LoggerBuffer(std::size_t bufferSize)
	: m_buffer(bufferSize + 1), m_context(""), m_padding(""),
	  m_consoleEnabled(false), m_console(&std::cout), m_consolePrefix(""),
	  m_fileEnabled(false), m_file(nullptr), m_filePrefix("")
{
	// Set the buffer
    char *bufferBegin = &m_buffer.front();
    setp(bufferBegin, bufferBegin + m_buffer.size() - 1);

	// Set console data
	setConsoleEnabled(true);

	// Set log file data
	setFileEnabled(false);
}

/*!
	Destructor
*/
LoggerBuffer::~LoggerBuffer()
{
	sync();
}

/*!
	Put character on overflow.

	It is called by public member functions such as sputc to write a character
	when there are no writing positions available at the put pointer (pptr).

	\return In case of success, the character put is returned, converted to
	a value of type int_type using member traits_type::to_int_type. Otherwise,
	it returns the end-of-file value (EOF) either if called with this value
	as argument c or to signal a failure.
*/
LoggerBuffer::int_type LoggerBuffer::overflow(int_type character)
{
    if (character != traits_type::eof()) {
		// Write the charater to the buffer and then increment pptr() by
		// calling pbump(1). It's always safe to write thecharater to the
		// buffer because we reserved an extra char at the end of our buffer
		// in the constructor.
        assert(std::less_equal<char *>()(pptr(), epptr()));
        *pptr() = character;
        pbump(1);

		// Flush the stream
        if (sync() == 0) {
            return character;
		}
    }

    return traits_type::eof();
}

/*!
	Synchronize stream buffer

	\return Returns 0 to indicates success, -1 to indicate failure.
*/
int LoggerBuffer::sync()
{
	return flush(false);
}

/*!
	Flushes stream buffer

	\param terminate if true ensures that a new line is added as last
	character
	\return Returns 0 to indicates success, -1 to indicate failure.
*/
int LoggerBuffer::flush(bool terminate)
{
	int status = 0;
	if (pptr() == pbase()) {
		return status;
	}

	// Identify line breakers
	std::vector<char *> linePointers;
	linePointers.push_back(pbase());
	for (char *p = pbase(), *e = pptr() - 1; p != e; ++p) {
		if (*p == '\n') {
			linePointers.push_back(p + 1);
		}
	}
	linePointers.push_back(pptr());

	// Detect if a new line must be added at the end
	terminate = (terminate && *(linePointers.back() - 1) != '\n');

	// Move the internal pointer
	std::ptrdiff_t nCharacters = pptr() - pbase();
	pbump(- nCharacters);

	// Write all lines
	for (unsigned int i = 0; i < linePointers.size() - 1; ++i) {
		char *firstCharacter = linePointers[i];
		char *lastCharacter  = linePointers[i + 1];
		std::ptrdiff_t lineSize = lastCharacter - firstCharacter;

		// Write to the console
		if (m_console && m_consoleEnabled) {
			if (!m_consolePrefix.empty()) {
				*m_console << m_consolePrefix << " :: ";
			}

			if (!m_context.empty()) {
				*m_console << m_context << " :: ";
			}

			if (!m_padding.empty()) {
				*m_console << m_padding;
			}

			bool success = m_console->write(firstCharacter, lineSize);
			if (!success) {
				status = -1;
			}

			if (terminate && lastCharacter == linePointers.back()) {
				*m_console << "\n";
			}
		}

		// Write to file
		if (m_file && m_fileEnabled && m_file->is_open()) {
			*m_file << "[" + getTimestamp() + "] ";

			if (!m_filePrefix.empty()) {
				*m_file << m_filePrefix << " :: ";
			}

			if (!m_context.empty()) {
				*m_file << m_context + " :: ";
			}

			if (!m_padding.empty()) {
				*m_file << m_padding;
			}

			bool success = m_file->write(firstCharacter, lineSize);
			if (!success) {
				status = -1;
			}

			if (terminate && lastCharacter == linePointers.back()) {
				*m_file << "\n";
			}
		}
	}

    return status;
}

/*!
	Enables the output on the console.

	\param enabled if set to true enables the console output.
*/
void LoggerBuffer::setConsoleEnabled(bool enabled)
{
	flush(true);

	m_consoleEnabled = enabled;
}

/*!
	Sets the stream to be used for the output on the console.

	\param console is the stream to be used for the output on the console
*/
void LoggerBuffer::setConsoleStream(std::ostream *console)
{
	flush(true);

	m_console = console;
}

/*!
	Gets the stream to be used for the output on the console.

	\result The stream to be used for the output on the console.
*/
std::ostream & LoggerBuffer::getConsoleStream()
{
	return *m_console;
}

/*!
	Sets the prefix for console output

	\param prefix is the prefix that will be prepended to every line of the
	console output
*/
void LoggerBuffer::setConsolePrefix(const std::string &prefix)
{
	flush(true);

	m_consolePrefix = prefix;
}

/*!
	Enables the output on the log file.

	The output on the log file can be enabled only id the log file path has
	been set. If there is no log file path specified, the output on the
	file will not be enabled.

	\param enabled if set to true enables the output on the log file.
*/
void LoggerBuffer::setFileEnabled(bool enabled)
{
	flush(true);

	m_fileEnabled = enabled;
}

/*!
	Sets the stream to be used for the output on the file.

	\param console is the stream to be used for the output on the file
*/
void LoggerBuffer::setFileStream(std::ofstream *file)
{
	flush(true);

	m_file = file;
}

/*!
	Gets the stream to be used for the output on the file.

	\result The stream to be used for the output on the file.
*/
std::ofstream & LoggerBuffer::getFileStream()
{
	return *m_file;
}

/*!
	Sets the prefix for file output

	\param prefix is the prefix that will be prepended to every line of the
	file output
*/
void LoggerBuffer::setFilePrefix(const std::string &prefix)
{
	flush(true);

	m_filePrefix = prefix;
}

/*!
	Sets the context of the output

	\param context is the context of the output
*/
void LoggerBuffer::setContext(const std::string &context)
{
	flush(true);

	m_context = context;
}

/*!
	Sets the padding to be prepended befor every message

	\param padding is the padding to be prepended befor every message
*/
void LoggerBuffer::setPadding(const std::string &padding)
{
	flush(true);

	m_padding = padding;
}

/*!
	Get current date/time, format is YYYY-MM-DD.HH:mm:ss
*/
const std::string LoggerBuffer::getTimestamp() const
{
	auto currentClock = std::chrono::system_clock::now();
	std::time_t currentTime = std::chrono::system_clock::to_time_t(currentClock);

	auto seconds = std::chrono::duration_cast<std::chrono::seconds>(currentClock.time_since_epoch());
	auto millisecs = std::chrono::duration_cast<std::chrono::milliseconds>(currentClock.time_since_epoch() - seconds);
	std::ostringstream millisecsConverter;
	millisecsConverter << (millisecs.count());

	char buffer[80];
	strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", localtime(&currentTime));
	std::string timestamp(buffer);
	timestamp += "." + millisecsConverter.str();

	return timestamp;
}

/*!
	@}
*/

/*!
	@ingroup Logger
	@{
*/

/*!
    \class Logger

    \brief Message logger.

    This class implements a message logger. The logger allows to write
    log messages both on the console and on a log file. The verbosity of
    the logger can be set for the console and the file independently.
*/

/*!
	Logger cache.
*/
Logger::Cache Logger::m_loggers = Logger::Cache();

/*!
	Default constructor.

	The constructor is private so that it can not be called.
*/
Logger::Logger(std::string name)
	: m_isMaster(true), m_buffer(256),
	  m_priority(log::NORMAL), m_visibility(log::MASTER),
	  m_consoleVerbosity(log::QUIET), m_fileVerbosity(log::QUIET)
{
	// Assigne the buffer to the stream
	rdbuf(&m_buffer);

	// Set buffer data
	m_buffer.setConsoleStream(&std::cout);
	m_buffer.setFileStream(&m_fileStream);

	// Set default file stream data
	m_fileHandler.setParallel(false);
	m_fileHandler.setAppendix("log");
	m_fileHandler.setName(name);
	m_fileHandler.setDirectory(".");

	// Set logger data
	setConsoleVerbosity(log::NORMAL);
	setFileVerbosity(log::QUIET);
	setPriority(log::NORMAL);
	setVisibility(log::MASTER);
}

/*!
	Destructor
*/
Logger::~Logger()
{
	if (m_fileStream.is_open()) {
		m_fileStream.close();
	}
}

/*!
	Returns an instance of the specified logger.

	This function returns an instance of the class. The firts time
	this function gets called it creates a new instance of the class.
	Calling the constructor publicly is not allowed. The constructor
	is private and is only called by this instance function.

	\param name is the name of the logger
	\result An instance of the specified logger.
*/
Logger & Logger::cout(std::string name)
{
	CacheIterator loggerItr = m_loggers.find(name);
	if (loggerItr != m_loggers.end()) {
		return *(loggerItr->second);
	}

	m_loggers[name] = std::unique_ptr<Logger>(new Logger(name));
	std::unique_ptr<Logger> &logger = m_loggers[name];

	return *logger;
}

/*!
	Sets the current priority of the messages.

	\param priority is the priority of the messags that will be printed
*/
void Logger::setPriority(log::Priority priority)
{
	m_priority = priority;

	setConsoleVerbosity(m_consoleVerbosity);
	setFileVerbosity(m_fileVerbosity);
}

/*!
	Gets the current priority of the messages.

	\result The gets the current priority of the messages.
*/
log::Priority Logger::getPriority()
{
	return m_priority;
}

/*!
	Sets the visibility of the messages.

	\param visibility is the visibility of the messags that will be printed
*/
void Logger::setVisibility(log::Visibility visibility)
{
	m_visibility = visibility;

	setConsoleVerbosity(m_consoleVerbosity);
	setFileVerbosity(m_fileVerbosity);
}

/*!
	Gets the visibility of the messages.

	\result The gets the visibility of the messages.
*/
log::Visibility Logger::getVisibility()
{
	return m_visibility;
}

/*!
	Sets the verbosity for the messages printed on the console.

	\param verbosity is the verbosity for the messages printed on the console
*/
void Logger::setConsoleVerbosity(log::Verbosity verbosity)
{
	m_consoleVerbosity = verbosity;
	if (m_consoleVerbosity == log::QUIET) {
		m_buffer.setConsoleEnabled(false);
	} else if (m_visibility == log::MASTER && !m_isMaster) {
		m_buffer.setConsoleEnabled(false);
	} else {
		m_buffer.setConsoleEnabled(m_consoleVerbosity >= m_priority);
	}
}

/*!
	Gets the verbosity for the messages printed on the console.

	\result The verbosity for the messages printed on the console.
*/
log::Verbosity Logger::getConsoleVerbosity()
{
	return m_consoleVerbosity;
}

/*!
	Sets the verbosity for the messages printed on the log file.

	\param verbosity is the verbosity for the messages printed on the log file
*/
void Logger::setFileVerbosity(log::Verbosity verbosity)
{
	m_fileVerbosity = verbosity;

	bool isFileEnabled;
	if (m_fileVerbosity == log::QUIET) {
		isFileEnabled = false;
	} else if (m_visibility == log::MASTER && !m_isMaster) {
		isFileEnabled = false;
	} else {
		isFileEnabled = (m_fileVerbosity >= m_priority);
	}
	m_buffer.setFileEnabled(isFileEnabled);

	if (isFileEnabled && !m_fileStream.is_open()) {
		openLogFile(false);
	}
}

/*!
	Gets the verbosity for the messages printed on the log file.

	\result The verbosity for the messages printed on the log file.
*/
log::Verbosity Logger::getFileVerbosity()
{
	return m_fileVerbosity;
}

/*!
	Sets the context used when printing the messages.

	\param context is the context used when printing the messages.
*/
void Logger::setContext(const std::string &context)
{
	m_context = context;

	m_buffer.setContext(m_context);
}

/*!
	Gets the context used when printing the messages.

	\result The the context used when printing the messages.
*/
std::string Logger::getContext()
{
	return m_context;
}

/*!
	Sets the indentation level used when printing the messages.

	\param delta is the relative indentation level
*/
void Logger::setIndentation(const int &delta)
{
	m_indentation = std::max(m_indentation + delta, 0);

	std::string padding = std::string(m_indentation, ' ');
	m_buffer.setPadding(padding);
}

/*!
	Gets the indentation level used when printing the messages.

	\result The the indentation level used when printing the messages.
*/
int Logger::getIndentation()
{
	return m_indentation;
}

/*!
	Resets the log file.
*/
void Logger::resetLogFile()
{
	openLogFile(true);
}

/*!
	Open the log file.
*/
void Logger::openLogFile(bool reset)
{
	std::string filePath = m_fileHandler.getPath();
	if (filePath.empty()) {
		return;
	}

	// Close previous stream
	if (m_fileStream.is_open()) {
		m_buffer.flush(true);
		m_fileStream.close();
	}

	// Open new stream
	std::ios_base::openmode mode;
	if (reset) {
		mode = std::ifstream::out;
	} else {
		mode = std::ifstream::app;
	}

	m_fileStream.rdbuf()->pubsetbuf(0, 0);
	m_fileStream.open(filePath, mode);
}

/*!
	Sets the parallel data of the logger.

	\param nProcessors is the total number of processors in the communicator
	\param rank is the parallel rank ofthe communicator
*/
void Logger::setParallel(int nProcessors, int rank)
{
	m_isMaster = (rank == 0);
	bool isParallel = (nProcessors > 1);

	// Set parallel data
	m_fileHandler.setParallel(isParallel);
	if (isParallel) {
		m_fileHandler.setBlock(rank);
	}

    // Rank prefix
	if (isParallel) {
		int nDigits = ceil(log10(nProcessors));
		std::ostringstream convert;
		convert << std::setw(nDigits) << rank;
		std::string rankPrefix = "#" + convert.str();

		m_buffer.setConsolePrefix(rankPrefix);
		m_buffer.setFilePrefix(rankPrefix);
	} else {
		m_buffer.setConsolePrefix("");
		m_buffer.setFilePrefix("");
	}

	// Reset verbosity
	setConsoleVerbosity(m_consoleVerbosity);
	setFileVerbosity(m_fileVerbosity);

	// Reopen the log file
	if (m_fileStream.is_open()) {
		openLogFile(false);
	}
}

/*!
	Prints a line in the log.

	\param line is the line to be printed
*/
void Logger::println(const std::string &line)
{
	print(line + '\n');
}

/*!
	Prints a line in the log.

	\param line is the line to be printed
	\param priority is the priority of the line that will be printed
*/
void Logger::println(const std::string &line, log::Priority priority)
{
	print(line + '\n', priority);
}

/*!
	Prints a line in the log.

	\param line is the line to be printed
	\param visibility is the visibility of the line that will be printed
*/
void Logger::println(const std::string &line, log::Visibility visibility)
{
	print(line + '\n', visibility);
}

/*!
	Prints a line in the log.

	\param line is the line to be printed
	\param priority is the priority of the line that will be printed
	\param visibility is the visibility of the line that will be printed
*/
void Logger::println(const std::string &line, log::Priority priority, log::Visibility visibility)
{
	print(line + '\n', priority, visibility);
}

/*!
	Prints a message in the log.

	\param message is the message to be printed
*/
void Logger::print(const std::string &message)
{
	(*this) << message;
}

/*!
	Prints a message in the log.

	\param message is the message to be printed
	\param priority is the priority of the message that will be printed
*/
void Logger::print(const std::string &message, log::Priority priority)
{
	setPriority(priority);

	(*this) << message;
}

/*!
	Prints a message in the log.

	\param message is the message to be printed
	\param visibility is the visibility of the message that will be printed
*/
void Logger::print(const std::string &message, log::Visibility visibility)
{
	setVisibility(visibility);

	(*this) << message;
}

/*!
	Prints a message in the log.

	\param message is the message to be printed
	\param priority is the priority of the message to be printed
	\param visibility is the visibility of the message that will be printed
*/
void Logger::print(const std::string &message, log::Priority priority, log::Visibility visibility)
{
	setPriority(priority);
	setVisibility(visibility);

	(*this) << message;
}

/*!
	@}
*/

/*!
	@ingroup Logger
	@{
*/

// Logger global functions
namespace log {

	/*!
		\enum Verbosity

		The verbosity of the message logger.

		\var Verbosity QUIET
		No messages wll be written to the output.

		\var Verbosity NORMAL
		Only messgaes with priority greater of equal "NORMAL" will be written
		to the output.

		\var Verbosity DEBUG
		Only messgaes with priority greater of equal "DEBUG" will be written
		to the output.
	*/

	/*!
		\typedef Priority

		The priority of a message.
	*/

	// Generic global functions

	/*!
		Returns an instance of the specified logger.

		\param name is the name of the logger
		\result An instance of the specified logger.
	*/
	Logger & cout(std::string name)
	{
		return Logger::cout(name);
	}

	// Manipulators global functions

	/*!
		Set the context of the output for the specified logger.

		\param logger is a reference pointing to the logger
		\param context is the context of the output
		\result A reference pointing to the logger received in input.
	*/
	Logger& setContext(Logger& logger, const std::string &context)
	{
		logger.setContext(context);

		return logger;
	}

	/*!
		Returns a logger manipulator that allows to change the context of the
		output.

		\param context is the context of the output
		\result A logger manipulator that allows to change the context of the
		output.
	*/
	LoggerManipulator<std::string> context(const std::string &context)
	{
		return LoggerManipulator<std::string>(setContext, context);
	}

	/*!
		Set the priority of the output for the specified logger.

		\param logger is a reference pointing to the logger
		\param priority is the priority of the output
		\result A reference pointing to the logger received in input.
	*/
	Logger& setPriority(Logger& logger, const log::Priority &priority)
	{
		logger.setPriority(priority);

		return logger;
	}

	/*!
		Returns a logger manipulator that allows to change the priority of the
		output.

		\param priority is the priority of the output
		\result A logger manipulator that allows to change the priority of the
		output.
	*/
	LoggerManipulator<log::Priority> priority(const log::Priority &priority)
	{
		return LoggerManipulator<log::Priority>(setPriority, priority);
	}

	/*!
		Set the visibility ofthe messages.

		\param logger is a reference pointing to the logger
		\param visibility is the visibility of the messages
		\result A reference pointing to the logger received in input.
	*/
	Logger& setVisibility(Logger& logger, const log::Visibility &visibility)
	{
		logger.setVisibility(visibility);

		return logger;
	}

	/*!
		Returns a logger manipulator that allows to change the visibility of the
		messages.

		\param visibility is the visibility of the messages
		\result A logger manipulator that allows to change the priority of the
		output.
	*/
	LoggerManipulator<log::Visibility> visibility(const log::Visibility &visibility)
	{
		return LoggerManipulator<log::Visibility>(setVisibility, visibility);
	}

	/*!
		Sets the verbosity for the messages printed on the console.

		\param logger is a reference pointing to the logger
		\param verbosity is the verbosity for the messages printed on the
		console
		\result A reference pointing to the logger received in input.
	*/
	Logger& setConsoleVerbosity(Logger& logger, const log::Verbosity &verbosity)
	{
		logger.setConsoleVerbosity(verbosity);

		return logger;
	}

	/*!
		Returns a logger manipulator that allows to change the verbosity for
		the messages printed on the console.

		\param verbosity is the verbosity for the messages printed on the
		console
		\result A logger manipulator that allows to change the verbosity for
		the messages printed on the console.
	*/
	LoggerManipulator<log::Verbosity> consoleVerbosity(const log::Verbosity &verbosity)
	{
		return LoggerManipulator<log::Verbosity>(setConsoleVerbosity, verbosity);
	}

	/*!
		Sets the verbosity for the messages printed on the file.

		\param logger is a reference pointing to the logger
		\param verbosity is the verbosity for the messages printed on the
		file
		\result A reference pointing to the logger received in input.
	*/
	Logger& setFileVerbosity(Logger& logger, const log::Verbosity &verbosity)
	{
		logger.setFileVerbosity(verbosity);

		return logger;
	}

	/*!
		Returns a logger manipulator that allows to change the verbosity for
		the messages printed on the file.

		\param verbosity is the verbosity for the messages printed on the
		file
		\result A logger manipulator that allows to change the verbosity for
		the messages printed on the file.
	*/
	LoggerManipulator<log::Verbosity> fileVerbosity(const log::Verbosity &verbosity)
	{
		return LoggerManipulator<log::Verbosity>(setFileVerbosity, verbosity);
	}

	/*!
		Sets the indentation level of the messages.

		\param logger is a reference pointing to the logger
		\param delta is the relative indentation level
		\result A reference pointing to the logger received in input.
	*/
	Logger& setIndentation(Logger& logger, const int &delta)
	{
		logger.setIndentation(delta);

		return logger;
	}

	/*!
		Returns a logger manipulator that allows to set the indentation level
		of the messages.

		\param delta is the relative indentation level
		\result A logger manipulator that allows to set the indentation level
		of the messages.
	*/
	LoggerManipulator<int> indent(const int &delta)
	{
		return LoggerManipulator<int>(setIndentation, delta);
	}
}

/*!
	@}
*/

}
