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
    Creates a new buffer.
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
    Synchronize stream buffer.

    \return Returns 0 to indicates success, -1 to indicate failure.
*/
int LoggerBuffer::sync()
{
    return flush(false);
}

/*!
    Flushes stream buffer.

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

            m_console->write(firstCharacter, lineSize);
            if ((m_console->rdstate() & std::ifstream::failbit ) != 0) {
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

            m_file->write(firstCharacter, lineSize);
            if ((m_file->rdstate() & std::ifstream::failbit ) != 0) {
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
    Sets the prefix for console output.

    \param prefix is the prefix that will be prepended to every line of the
    console output
*/
void LoggerBuffer::setConsolePrefix(const std::string &prefix)
{
    flush(true);

    m_consolePrefix = prefix;
}

/*!
    Gets the prefix for console output.

    \result The prefix for console output.
*/
std::string LoggerBuffer::getConsolePrefix() const
{
    return m_consolePrefix;
}

/*!
    Enables the output on the log file.

    The output on the log file can be enabled only id the log file path has
    been set. If there is no log file path specified, the output on the
    file will not be enabled.

    \param enabled if set to true enables the output on the log file
*/
void LoggerBuffer::setFileEnabled(bool enabled)
{
    flush(true);

    m_fileEnabled = enabled;
}

/*!
    Sets the stream to be used for the output on the file.

    \param file is the stream to be used for the output on the file
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
    Sets the prefix for file output.

    \param prefix is the prefix that will be prepended to every line of the
    file output
*/
void LoggerBuffer::setFilePrefix(const std::string &prefix)
{
    flush(true);

    m_filePrefix = prefix;
}

/*!
    Gets the prefix for file output.

    \result The prefix for file output.
*/
std::string LoggerBuffer::getFilePrefix() const
{
    return m_filePrefix;
}

/*!
    Sets the context of the output.

    \param context is the context of the output
*/
void LoggerBuffer::setContext(const std::string &context)
{
    flush(true);

    m_context = context;
}

/*!
    Sets the padding to be prepended befor every message.

    \param padding is the padding to be prepended befor every message
*/
void LoggerBuffer::setPadding(const std::string &padding)
{
    flush(true);

    m_padding = padding;
}

/*!
    Get current date/time, format is "YYYY-MM-DD.HH:mm:ss".
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
    Default constructor.

    The constructor is private so that it can not be called.
*/
Logger::Logger(std::ostream *consoleStream, std::ofstream *fileStream,
            const int &nProcessors, const int &rank)
    : m_nProcessors(nProcessors), m_rank(rank), m_buffer(256),
    m_priority(log::NORMAL), m_visibility(log::MASTER),
    m_consoleVerbosity(log::NORMAL), m_fileVerbosity(log::NORMAL)
{
    // Assigne the buffer to the stream
    rdbuf(&m_buffer);

    // Set buffer data
    setConsoleStream(consoleStream);
    setFileStream(fileStream);

    // Set parallel data
    if (m_nProcessors > 1) {
        int nDigits = ceil(log10(m_nProcessors));
        std::ostringstream convert;
        convert << std::setw(nDigits) << m_rank;
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

    // Set logger data
    setPriority(log::NORMAL);
    setVisibility(log::MASTER);
}

/*!
    Sets the stream to be used for the output on the console.

    \param console is the stream to be used for the output on the console
*/
void Logger::setConsoleStream(std::ostream *console)
{
    m_buffer.setConsoleStream(console);
}

/*!
    Gets the stream to be used for the output on the console.

    \result The stream to be used for the output on the console.
*/
std::ostream & Logger::getConsoleStream()
{
    return m_buffer.getConsoleStream();
}

/*!
    Sets the stream to be used for the output on the file.

    \param file is the stream to be used for the output on the file
*/
void Logger::setFileStream(std::ofstream *file)
{
    m_buffer.setFileStream(file);
}


/*!
    Gets the stream to be used for the output on the file.

    \result The stream to be used for the output on the file.
*/
std::ofstream & Logger::getFileStream()
{
    return m_buffer.getFileStream();
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
    Set the verbosity for the messages printed both on the console and on the
    log file.

    \param verbosity is the verbosity for the messages
*/
void Logger::setVerbosities(log::Verbosity verbosity)
{
    setConsoleVerbosity(verbosity);
    setFileVerbosity(verbosity);
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
    } else if (m_visibility == log::MASTER && (m_rank != 0)) {
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
    Gets the prefix for the messages printed on the console.

    \result The prefix for the messages printed on the console.
*/
std::string Logger::getConsolePrefix()
{
    return m_buffer.getConsolePrefix();
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
    } else if (m_visibility == log::MASTER && (m_rank != 0)) {
        isFileEnabled = false;
    } else {
        isFileEnabled = (m_fileVerbosity >= m_priority);
    }
    m_buffer.setFileEnabled(isFileEnabled);
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
    Gets the prefix for the messages printed on the log file.

    \result The prefix for the messages printed on the log file.
*/
std::string Logger::getFilePrefix()
{
    return m_buffer.getFilePrefix();
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
    Gets the processors in the communicator.

    \result The processors in the communicator.
*/
int Logger::getProcessorCount()
{
    return m_nProcessors;
}

/*!
    Gets the rank in the communicator.

    \result The rank in the communicator.
*/
int Logger::getRank()
{
    return m_rank;
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
    log::Priority prevPriority   = getPriority();
    log::Visibility prevVisibility = getVisibility();

    setPriority(priority);
    setVisibility(visibility);

    (*this) << message;

    setPriority(prevPriority);
    setVisibility(prevVisibility);
}

/*!
    @}
*/

/*!
    @ingroup Logger
    @{
*/

/*!
    \class LoggerManager

    \brief Manager for the loggers.

    This class implements a manager for the loggers. The manager allowes the
    different loggers to work together.
*/

/*
    Initialize logger manager instance.
*/
std::unique_ptr<LoggerManager> LoggerManager::m_manager = nullptr;

/*
    Default name for the logger.
*/
std::string LoggerManager::BITPIT_LOG_NAME = "bitpit";

/*
    Default directory for the logger.
*/
std::string LoggerManager::BITPIT_LOG_DIRECTORY = ".";

/*!
    Default constructor.
*/
LoggerManager::LoggerManager()
    : m_defaultName(BITPIT_LOG_NAME), m_defaultDirectory(BITPIT_LOG_DIRECTORY),
    m_mode(log::SEPARATE)
{

}

/*!
    Destructor.
*/
LoggerManager::~LoggerManager()
{
    // Get a list of the loggers
    std::vector<std::string> loggerNames;
    loggerNames.reserve(m_loggers.size());
    for (const auto &entry : m_loggers) {
        loggerNames.push_back(entry.first);
    }

    // Destry all the loggers
    for (const std::string &name : loggerNames) {
        destroy(name, true);
    }
}

/*!
    Returns an instance of the logger manager.

    \result An instance of the logger manager.
*/
LoggerManager & LoggerManager::manager()
{
    if (!m_manager) {
        m_manager = std::unique_ptr<LoggerManager>(new LoggerManager());
    }

    return *m_manager;
}

/*!
    Returns an instance of defualt logger.

    \result An instance of the default specified logger.
*/
Logger & LoggerManager::cout()
{
    return cout(m_defaultName);
}

/*!
    Returns an instance of the specified logger.

    This function returns an instance of the specified logger. If the logger
    does not exists a new instance will be created.

    \param name is the name of the logger
    \result An instance of the specified logger.
*/
Logger & LoggerManager::cout(const std::string &name)
{
    // The logger has to be created
    if (m_loggers.count(name) == 0) {
        if (!isInitialized()) {
            setMode(log::SEPARATE);
        }

        if (name == m_defaultName) {
            _create(name, false, m_defaultDirectory, 1, 0);
        } else {
            create(name, false, m_defaultDirectory, 1, 0);
        }
    }

    // Return the logger
    return *(m_loggers.at(name));
}

/*!
    Initializes the log manager.

    \param mode is the mode that will be set
    \param reset if true the log files will be reset
    \param nProcessors is the total number of processors in the communicator
    \param rank is the parallel rank in the communicator
*/
void LoggerManager::initialize(log::Mode mode, bool reset,
                            const int &nProcessors, const int &rank)
{
    initialize(mode, m_defaultName, reset, m_defaultDirectory, nProcessors, rank);
}

/*!
    Initializes the log manager.

    \param mode is the mode that will be set
    \param reset if true the log files will be reset
    \param directory is the defualt directory for saving the log files
    \param nProcessors is the total number of processors in the communicator
    \param rank is the parallel rank in the communicator
*/
void LoggerManager::initialize(log::Mode mode, bool reset, const std::string &directory,
                            const int &nProcessors, const int &rank)
{
    initialize(mode, m_defaultName, reset, directory, nProcessors, rank);
}

/*!
    Initializes the log manager.

    \param mode is the mode that will be set
    \param name is the name for the default logger
    \param reset if true the log files will be reset
    \param directory is the defualt directory for saving the log files
    \param nProcessors is the total number of processors in the communicator
    \param rank is the parallel rank in the communicator
*/
void LoggerManager::initialize(log::Mode mode, const std::string &name, bool reset,
                            const std::string &directory,
                            const int &nProcessors, const int &rank)
{
    if (isInitialized()) {
        log::cout().println("Logger initialization has to be called before creating the loggers.");
        return;
    }

    // Set mode
    setMode(mode);

    // Set the dafault data
    m_defaultName      = name;
    m_defaultDirectory = directory;

    // Create the logger
    _create(m_defaultName, reset, directory, nProcessors, rank);
}

/*!
    Creates a new logger.

    \param name is the name for the logger
    \param reset if true the log files will be reset
    \param nProcessors is the total number of processors in the communicator
    \param rank is the parallel rank in the communicator
*/
void LoggerManager::create(const std::string &name, bool reset,
                        const int &nProcessors, const int &rank)
{
    create(name, reset, m_defaultDirectory, nProcessors, rank);
}

/*!
    Creates a new logger.

    \param name is the name for the logger
    \param reset if true the log files will be reset
    \param directory is the directory for saving the log files
    \param nProcessors is the total number of processors in the communicator
    \param rank is the parallel rank in the communicator
*/
void LoggerManager::create(const std::string &name, bool reset,
                        const std::string &directory,
                        const int &nProcessors, const int &rank)
{
    // Its not possible to create a log with the default name nor a log
    // with the same name of an existent logger nor a log with an empty
    // name
    if (exists(name)) {
        m_loggerUsers[name]++;

        cout(name) << "Detected an attemp to create the logger \"" << name << "\" twice" << std::endl;
        cout(name) << "The previously created logger will be used." << std::endl;
        return;
    } else if (name == m_defaultName) {
        cout(m_defaultName) << "Detected an attemp to overwrite the default logger" << std::endl;
        cout(m_defaultName) << "The name of the default logger is \"" << name << "\"" << std::endl;
        return;
    } else if (name == "") {
        cout() << "Detected an attemp to create a logger with an empty name" << std::endl;
        return;
    }

    // Create the logger
    if (m_mode == log::SEPARATE) {
        _create(name, reset, directory, nProcessors, rank);
    } else {
        _create(name, cout(m_defaultName));
    }
}

/*!
    Destroys the specified logger.

    The logger can be shared among differen users, a logger is destoryed
    only if it has no users.

    \param name is the name of the logger
    \param force controls if the logger will be destory also if it still
    has users
    \result True if the logger has been destroyed, false otherwise.
*/
bool LoggerManager::destroy(const std::string &name, bool force)
{
    if (!exists(name)) {
        return false;
    }

    // Decrement the users
    int &nUsers = m_loggerUsers[name];
    nUsers--;

    // If the logger has no users we can delete it
    if (nUsers == 0 || force) {
        // Remove the logger from the manager
        m_loggers.erase(name);
        m_loggerUsers.erase(name);

        // If the logger has it own file stream, delete it
        if (m_fileStreams.count(name) > 0) {
            // Close the stream
            std::ofstream &fileStream = *(m_fileStreams.at(name));
            fileStream.close();

            // Remove the stream
            m_fileStreams.erase(name);
        }

        return true;
    } else {
        return false;
    }
}

/*!
    Check if the specified logger exists.

    \param name is the name of the logger
    \result True if the logger exists, false otherwise.
*/
bool LoggerManager::exists(const std::string &name) const
{
    return (m_loggers.count(name) != 0);
}

/*!
    Checks if the logger manager has been initialized.

    Explicit initialization of the manager is done using the function
    'initialize', implicit initialization happnes the first time a
    logger is created.

    \result Returns true if the logger manager has been initialized, false
    otherwise.
*/
bool LoggerManager::isInitialized() const
{
    return (m_loggers.size() > 0);
}

/*!
    Sets the policty for the loggers.

    The policy has to be set before creating the loggers.

    \param mode is the mode that will be set
    \result Returns true if the mode was successfully set, false otherwise.
*/
bool LoggerManager::setMode(log::Mode mode)
{
    if (isInitialized()) {
        cout().println("The policy has to be set before creating the loggers.");
        return false;
    }

    m_mode = mode;
    return true;
}

/*!
    Sets the policty for the loggers.

    \result The mode set for the loggers.
*/
log::Mode LoggerManager::getMode() const
{
    return m_mode;
}

/*!
    Internal function to create a logger.

    \param name is the name for the logger
    \param reset if true the log files will be reset
    \param directory is the directory for saving the log files
    \param nProcessors is the total number of processors in the communicator
    \param rank is the parallel rank in the communicator
*/
void LoggerManager::_create(const std::string &name, bool reset,
                            const std::string &directory,
                            const int &nProcessors, const int &rank)
{
    // Get the file path
    FileHandler fileHandler;
    fileHandler.setDirectory(directory);
    fileHandler.setName(name);
    fileHandler.setAppendix("log");
    fileHandler.setParallel(nProcessors > 1);
    if (nProcessors > 1) {
        fileHandler.setBlock(rank);
    }

    std::string filePath = fileHandler.getPath();

    // Create the file stream
    std::ios_base::openmode fileMode;
    if (reset) {
        fileMode = std::ofstream::out;
    } else {
        fileMode = std::ofstream::app;
    }

    m_fileStreams[name] = std::unique_ptr<std::ofstream>(new std::ofstream());
    std::ofstream &fileStream = *(m_fileStreams.at(name));
    fileStream.rdbuf()->pubsetbuf(0, 0);
    fileStream.open(filePath, fileMode);

    // Use cout as console stream
    std::ostream &consoleStream = std::cout;

    // Create the logger
    m_loggers[name]     = std::unique_ptr<Logger>(new Logger(&consoleStream, &fileStream, nProcessors, rank));
    m_loggerUsers[name] = 1;
}

/*!
    Internal function to create a logger.

    \param name is the name for the logger
    \param master is a reference to the logger from which the settings will
    be copied from
*/
void LoggerManager::_create(const std::string &name, Logger &master)
{
    std::ostream &consoleStream = master.getConsoleStream();
    std::ofstream &fileStream   = master.getFileStream();

    int rank        = master.getRank();
    int nProcessors = master.getProcessorCount();

    m_loggers[name]     = std::unique_ptr<Logger>(new Logger(&consoleStream, &fileStream, nProcessors, rank));
    m_loggerUsers[name] = 1;

    // Import logger settings
    Logger &logger = *(m_loggers.at(name));
    logger.setConsoleVerbosity(master.getConsoleVerbosity());
    logger.setFileVerbosity(master.getFileVerbosity());
    logger.setVisibility(master.getVisibility());
    logger.setPriority(master.getPriority());
}

/*!
    Sets the verbosity for the messages printed both on the console and on the
    log file.

    \param verbosity is the verbosity for the messages
*/
void LoggerManager::setVerbosities(log::Verbosity verbosity)
{
    for (auto &entry : m_loggers) {
        entry.second->setVerbosities(verbosity);
    }
}

/*!
    Sets the verbosity for the messages printed on the console.

    \param verbosity is the verbosity for the messages printed on the console
*/
void LoggerManager::setConsoleVerbosity(log::Verbosity verbosity)
{
    for (auto &entry : m_loggers) {
        entry.second->setConsoleVerbosity(verbosity);
    }
}

/*!
    Sets the verbosity for the messages printed on the file.

    \param verbosity is the verbosity for the messages printed on the file
*/
void LoggerManager::setFileVerbosity(log::Verbosity verbosity)
{
    for (auto &entry : m_loggers) {
        entry.second->setFileVerbosity(verbosity);
    }
}

/*!
    Gets the defualt logger name.

    \param The defualt logger name.
*/
std::string LoggerManager::getDefaultName() const
{
    return m_defaultName;
}

/*!
    Gets the defualt logger directory.

    \param The defualt logger directory.
*/
std::string LoggerManager::getDefaultDirectory() const
{
    return m_defaultDirectory;
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
        Returns the logger manager.

        \result The logger manager.
    */
    LoggerManager & manager()
    {
        return LoggerManager::manager();
    }

    /*!
        Returns an instance of the default logger.

        \result An instance of the default logger.
    */
    Logger & cout()
    {
        return manager().cout();
    }

    /*!
        Returns an instance of the specified logger.

        \param name is the name of the logger
        \result An instance of the specified logger.
    */
    Logger & cout(std::string name)
    {
        return manager().cout(name);
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
        Sets the verbosity for the messages printed both on the console and on
        the log file.

        \param logger is a reference pointing to the logger
        \param verbosity is the verbosity for the messages
        \result A reference pointing to the logger received in input.
    */
    Logger& setVerbosities(Logger& logger, const log::Verbosity &verbosity)
    {
        logger.setVerbosities(verbosity);

        return logger;
    }

    /*!
        Returns a logger manipulator that allows to change the verbosity for
        the messages printed both on the console and on the log file.

        \param verbosity is the verbosity for the messages
        \result A logger manipulator that allows to change the verbosity for
        the messages printed on the console.
    */
    LoggerManipulator<log::Verbosity> verbosities(const log::Verbosity &verbosity)
    {
        return LoggerManipulator<log::Verbosity>(setVerbosities, verbosity);
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
