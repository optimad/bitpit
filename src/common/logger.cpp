/*---------------------------------------------------------------------------*\
*
*  bitpit
*
*  Copyright (C) 2015-2021 OPTIMAD engineering Srl
*
*  -------------------------------------------------------------------------
*  License
*  This file is part of bitpit.
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

#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <functional>
#include <sys/types.h>
#include <sys/stat.h>
#include <type_traits>

#include "compiler.hpp"
#include "fileHandler.hpp"
#include "logger.hpp"

namespace bitpit{

/*!
    \class LoggerBuffer
    \ingroup common_logger
    \brief Stream buffer for the message logger.

    This class implements a stream buffer used by the message logger.
*/

/*!
    Creates a new buffer.

    \param nProcesses is the total number of processes in the communicator
    \param rank is the parallel rank in the communicator
    \param bufferSize is the size of the internal buffer
*/
LoggerBuffer::LoggerBuffer(int nProcesses, int rank, std::size_t bufferSize)
    : m_buffer(bufferSize + 1),
      m_nProcesses(nProcesses), m_rank(rank),
      m_consoleEnabled(false), m_console(&std::cout),
      m_fileEnabled(false), m_file(nullptr),
      m_settings(std::make_shared<Settings>())
{
    // Set the buffer
    char *bufferBegin = &m_buffer.front();
    setp(bufferBegin, bufferBegin + m_buffer.size() - 1);

    // Initialize console stream
    setConsoleEnabled(true);

    // Initialize file stream
    setFileEnabled(false);

    // Set parallel data
    if (m_nProcesses > 1) {
        int nDigits = ceil(log10(m_nProcesses));
        std::ostringstream convert;
        convert << std::setw(nDigits) << m_rank;
        m_rankPrefix = "#" + convert.str();
    } else {
        m_rankPrefix = "";
    }

    // Default settings
    m_settings->context = "";
    m_settings->indentation = 0;
    m_settings->consoleTimestampEnabled = false;
    m_settings->fileTimestampEnabled = false;
}

/*!
    Destructor
*/
LoggerBuffer::~LoggerBuffer()
{
    sync();
}

/*!
    Count the number of processes in the communicator.

    \result The number of processes in the communicator.
*/
int LoggerBuffer::getProcessCount() const
{
    return m_nProcesses;
}

/*!
    Gets the rank in the communicator.

    \result The rank in the communicator.
*/
int LoggerBuffer::getRank() const
{
    return m_rank;
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
        // Write the character to the buffer and then increment pptr() by
        // calling pbump(1). It's always safe to write the charater to the
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
    // Get buffer information
    const char *bufferBegin = pbase();
    const char *bufferEnd   = pptr();
    if (bufferBegin == bufferEnd) {
        return 0;
    }

    // Flush buffer a line at the time
    const char *lineBegin = nullptr;
    const char *lineEnd   = bufferBegin;
    while (lineEnd != bufferEnd) {
        // Update line information
        lineBegin = lineEnd;
        for (lineEnd = lineBegin; lineEnd != bufferEnd; ++lineEnd) {
            if (*lineEnd == '\n') {
                lineEnd += 1;
                break;
            }
        }

        // Detect if a new line must be added at the end of the line
        bool terminateLine = false;
        if (terminate && (lineEnd == bufferEnd)) {
            if (bufferBegin == bufferEnd) {
                terminateLine = true;
            } else if (*(bufferEnd - 1) != '\n') {
                terminateLine = true;
            }
        }

        // Get timestamp
        std::string consoleTimestamp;
        std::string fileTimestamp;
        if (m_settings->consoleTimestampEnabled || m_settings->fileTimestampEnabled) {
            std::string timestamp = getTimestamp();
            if (m_settings->consoleTimestampEnabled) {
                consoleTimestamp = timestamp;
            }
            if (m_settings->fileTimestampEnabled) {
                fileTimestamp = timestamp;
            }
        }

        // Flush line to console
        if (m_console && m_consoleEnabled) {
            int status = flushLine(*m_console, lineBegin, lineEnd, consoleTimestamp, terminateLine);
            if (status != 0) {
                return status;
            }
        }

        // Flush line to file
        if (m_file && m_fileEnabled && m_file->is_open()) {
            int status = flushLine(*m_file, lineBegin, lineEnd, fileTimestamp, terminateLine);
            if (status != 0) {
                return status;
            }
        }
    }

    // Reset the internal pointer
    pbump(bufferBegin - bufferEnd);

    return 0;
}

/*!
    Flushes stream buffer.

    \param stream is the steam the line will be flushed to
    \param begin refers to the first character of the line
    \param end refers to the past-the-last character of the line
    \param timestamp if the timestamp that will be printed at the beginning
    of the line, if an empty string is provided, no timestamp information
    will be printed
    \param terminate if set to true a new line character will be printed at
    \the end of the line
    \return Returns 0 to indicates success, -1 to indicate failure.
*/
int LoggerBuffer::flushLine(std::ostream &stream, const char *begin, const char *end,
                            const std::string &timestamp, bool terminate)
{
    if (!timestamp.empty()) {
        stream << "[" + timestamp + "] ";
        if ((stream.rdstate() & std::ifstream::failbit) != 0) {
            return -1;
        }
    }

    if (!m_rankPrefix.empty()) {
        stream << m_rankPrefix << " :: ";
        if ((stream.rdstate() & std::ifstream::failbit) != 0) {
            return -1;
        }
    }

    if (!m_settings->context.empty()) {
        stream << m_settings->context + " :: ";
        if ((stream.rdstate() & std::ifstream::failbit) != 0) {
            return -1;
        }
    }

    for (int i = 0; i < m_settings->indentation; ++i) {
        stream << " ";
        if ((stream.rdstate() & std::ifstream::failbit) != 0) {
            return -1;
        }
    }

    stream.write(begin, end - begin);
    if ((stream.rdstate() & std::ifstream::failbit) != 0) {
        return -1;
    }

    if (terminate) {
        stream << "\n";
        if ((stream.rdstate() & std::ifstream::failbit) != 0) {
            return -1;
        }
    }

    return 0;
}

/*!
    Enables the output on the console.

    \param enabled if set to true enables the console output.
*/
void LoggerBuffer::setConsoleEnabled(bool enabled)
{
    if (m_consoleEnabled == enabled) {
        return;
    }

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
    Enables the output on the log file.

    The output on the log file can be enabled only id the log file path has
    been set. If there is no log file path specified, the output on the
    file will not be enabled.

    \param enabled if set to true enables the output on the log file
*/
void LoggerBuffer::setFileEnabled(bool enabled)
{
    if (m_fileEnabled == enabled) {
        return;
    }

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
    Get a constant pointer to the buffer settings.

    \result A constant pointer to the buffer settings.
*/
const LoggerBuffer::Settings * LoggerBuffer::getSettings() const
{
    return m_settings.get();
}

/*!
    Sets buffer settings.

    \param settings are the buffer settings
*/
void LoggerBuffer::setSettings(const std::shared_ptr<Settings> &settings)
{
    if (settings.get() == m_settings.get()) {
        return;
    }

    flush(true);

    m_settings = settings;
}

/*!
    Get current date/time, format is "YYYY-MM-DD HH:mm:ss".
*/
std::string LoggerBuffer::getTimestamp() const
{
    static const std::string TIMESTAMP_FORMAT = "%Y-%m-%d %H:%M:%S";

    auto currentClock = std::chrono::system_clock::now();
    std::time_t currentTime = std::chrono::system_clock::to_time_t(currentClock);

    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(currentClock.time_since_epoch());
    auto millisecs = std::chrono::duration_cast<std::chrono::milliseconds>(currentClock.time_since_epoch() - seconds);
    std::array<char, 11> millisecsBuffer;
    std::snprintf(millisecsBuffer.data(), millisecsBuffer.size(), "%03u", static_cast<unsigned int>(millisecs.count()));

    std::string timestamp;
    timestamp.resize(23);
    std::strftime(&timestamp[0], timestamp.size(), TIMESTAMP_FORMAT.c_str(), localtime(&currentTime));
    timestamp[19] = '.';
    timestamp[20] = millisecsBuffer[0];
    timestamp[21] = millisecsBuffer[1];
    timestamp[22] = millisecsBuffer[2];

    return timestamp;
}

/*!
    \class Logger
    \ingroup common_logger
    \brief Message logger.

    This class implements a message logger. The logger allows to write
    log messages both on the console and on a log file. The verbosity of
    the logger can be set for the console and the file independently.
*/

/*!
    Default constructor.

    The constructor is private so that it can not be called.
*/
Logger::Logger(const std::string &name, const std::shared_ptr<LoggerBuffer> &buffer)
    : std::ios(nullptr), std::ostream(nullptr),
      m_name(name),
      m_buffer(buffer), m_bufferSettings(std::make_shared<LoggerBuffer::Settings>()),
      m_defaultSeverity(log::INFO), m_defaultVisibility(log::VISIBILITY_MASTER),
      m_consoleDisabledThreshold(log::NOTSET), m_consoleVerbosityThreshold(log::INFO),
      m_fileDisabledThreshold(log::NOTSET), m_fileVerbosityThreshold(log::INFO)
{
    // St stream buffer
    basic_ios<char>::rdbuf(m_buffer.get());

    // Default settings
    m_bufferSettings->context = "";
    m_bufferSettings->indentation = 0;
    m_bufferSettings->consoleTimestampEnabled = false;
    m_bufferSettings->fileTimestampEnabled = true;
}

/*!
    Provides an overriding verbosity level which takes precedence over the
    logger’s own level.

    When the need arises to temporarily throttle logging output down across
    the whole application, this function can be useful. Its effect is to
    disable all logging calls of severity level and below, so that if you
    call it with a value of INFO, then all INFO and DEBUG events would be
    discarded, whereas those of severity WARNING and above would be processed
    according to the logger’s effective level.

    If the a NOTSET level is specified, it effectively removes this overriding
    level, so that logging output again depends on the effective verbosity
    level of the logger.

    \param level is the overriding verbosity level
*/
void Logger::disable(log::Level level)
{
    disableConsole(level);
    disableFile(level);
}

/*!
    Provides an overriding console verbosity level which takes precedence over
    the logger’s own level.

    When the need arises to temporarily throttle logging output down across
    the whole application, this function can be useful. Its effect is to
    disable all logging calls of severity level and below, so that if you
    call it with a value of INFO, then all INFO and DEBUG events would be
    discarded, whereas those of severity WARNING and above would be processed
    according to the logger’s effective level.

    If the NOTSET level is specified, it effectively removes this overriding
    level, so that logging output again depends on the effective verbosity
    level of the logger.

    \param level is the overriding verbosity level
*/
void Logger::disableConsole(log::Level level)
{
    m_consoleDisabledThreshold = level;
}

/*!
    Provides an overriding file verbosity level which takes precedence over
    the logger’s own level.

    When the need arises to temporarily throttle logging output down across
    the whole application, this function can be useful. Its effect is to
    disable all logging calls of severity level and below, so that if you
    call it with a value of INFO, then all INFO and DEBUG events would be
    discarded, whereas those of severity WARNING and above would be processed
    according to the logger’s effective level.

    If the NOTSET level is specified, it effectively removes this overriding
    level, so that logging output again depends on the effective verbosity
    level of the logger.

    \param level is the overriding verbosity level
*/
void Logger::disableFile(log::Level level)
{
    m_fileDisabledThreshold = level;
}

/*!
    Sets the default severity of the messages.

    \param severity is the default severity of the messages
*/
void Logger::setDefaultSeverity(log::Level severity)
{
    // Set default severity
    assert(severity != log::Level::NOTSET);
    m_defaultSeverity = severity;
}

/*!
    Gets the default severity of the messages.

    \result The gets the default severity of the messages.
*/
log::Level Logger::getDefaultSeverity()
{
    return m_defaultSeverity;
}

/*!
    Sets the default priority of the messages.

    \param priority is the default priority of the messages
*/
void Logger::setPriority(log::Priority priority)
{
    setDefaultSeverity(priority);
}

/*!
    Gets the current priority of the messages.

    \result The gets the current priority of the messages.
*/
log::Priority Logger::getPriority()
{
    return getDefaultSeverity();
}

/*!
    Sets the default visibility of the messages.

    \param visibility is the default visibility of the messages
*/
void Logger::setDefaultVisibility(log::Visibility visibility)
{
    // Set default visibility
    assert(visibility != log::VISIBILITY_NOTSET);
    m_defaultVisibility = visibility;
}

/*!
    Gets the default visibility of the messages.

    \result The gets the visibility of the messages.
*/
log::Visibility Logger::getDefaultVisibility()
{
    return m_defaultVisibility;
}

/*!
    Sets the visibility of the messages.

    \param visibility is the default visibility of the messages
*/
void Logger::setVisibility(log::Visibility visibility)
{
    setDefaultVisibility(visibility);
}

/*!
    Gets the visibility of the messages.

    \result The gets the visibility of the messages.
*/
log::Visibility Logger::getVisibility()
{
    return getDefaultVisibility();
}

/*!
    Sets the verbosity threshold for messages printed on the console and on
    the log file.

    Logging messages which are less severe than the verbosity threshold
    will be ignored; logging messages which have severity equal to the
    verbosity threshold or higher will be printed.

    \param threshold is the verbosity threshold
*/
void Logger::setVerbosities(log::Level threshold)
{
    setConsoleVerbosity(threshold);
    setFileVerbosity(threshold);
}

/*!
    Enables the timestamp both on the console and on the log file.

    \param enabled if set to true enables the timestamp both on the console
    and on the log file.
*/
void Logger::setTimestampEnabled(bool enabled)
{
    setConsoleTimestampEnabled(enabled);
    setFileTimestampEnabled(enabled);
}

/*!
    Returns true if the timestamp is enabled on the console, false otherwise.

    \result Returns true if the timestamp is enabled on the console, false
    otherwise.
*/
bool Logger::isConsoleTimestampEnabled() const
{
    return m_bufferSettings->consoleTimestampEnabled;
}

/*!
    Enables the timestamp on the console.

    \param enabled if set to true enables the timestamp on the console.
*/
void Logger::setConsoleTimestampEnabled(bool enabled)
{
    if (m_bufferSettings->consoleTimestampEnabled == enabled) {
        return;
    }

    m_buffer->flush(true);

    m_bufferSettings->consoleTimestampEnabled = enabled;
}

/*!
    Sets the verbosity threshold for messages printed on the console.

    Logging messages which are less severe than the verbosity threshold
    will be ignored; logging messages which have severity equal to the
    verbosity threshold or higher will be printed.

    \param threshold is the verbosity threshold
*/
void Logger::setConsoleVerbosity(log::Level threshold)
{
    // Set verbosity threshold
    m_consoleVerbosityThreshold = threshold;
}

/*!
    Enable the buffer streams that will print a message with the specified properties.

    \param severity is the severity of the message to be printed
    \param visibility is the visibility of the message to be printed
*/
void Logger::enableBufferStreams(log::Level severity, log::Visibility visibility)
{
    // Enable console stream
    bool isConsoleEnabled = true;
    if (severity <= m_consoleDisabledThreshold) {
        isConsoleEnabled = false;
    } else if (visibility == log::VISIBILITY_MASTER && (m_buffer->getRank() != 0)) {
        isConsoleEnabled = false;
    } else {
        isConsoleEnabled = (severity >= m_consoleVerbosityThreshold);
    }
    m_buffer->setConsoleEnabled(isConsoleEnabled);

    // Enable file stream
    bool isFileEnabled = true;
    if (severity <= m_fileDisabledThreshold) {
        isFileEnabled = false;
    } else if (visibility == log::VISIBILITY_MASTER && (m_buffer->getRank() != 0)) {
        isFileEnabled = false;
    } else {
        isFileEnabled = (severity >= m_fileVerbosityThreshold);
    }
    m_buffer->setFileEnabled(isFileEnabled);
}

/*!
    Gets the verbosity threshold for messages printed on the console.

    Logging messages which are less severe than the verbosity threshold
    will be ignored; logging messages which have severity equal to the
    verbosity threshold or higher will be printed.

    \result The verbosity for the messages printed on the console.
*/
log::Level Logger::getConsoleVerbosity()
{
    return m_consoleVerbosityThreshold;
}

/*!
    Returns true if the timestamp is enabled on the log file, false otherwise.

    \result Returns true if the timestamp is enabled on the log file, false
    otherwise.
*/
bool Logger::isFileTimestampEnabled() const
{
    return m_bufferSettings->fileTimestampEnabled;
}

/*!
    Enables the timestamp on the log file.

    \param enabled if set to true enables the timestamp on the log file
*/
void Logger::setFileTimestampEnabled(bool enabled)
{
    if (m_bufferSettings->fileTimestampEnabled == enabled) {
        return;
    }

    m_buffer->flush(true);

    m_bufferSettings->fileTimestampEnabled = enabled;
}

/*!
    Sets the verbosity threshold for messages printed on the log file.

    Logging messages which are less severe than the verbosity threshold
    will be ignored; logging messages which have severity equal to the
    verbosity threshold or higher will be printed.

    \param threshold is the verbosity threshold
*/
void Logger::setFileVerbosity(log::Level threshold)
{
    // Set verbosity threshold
    m_fileVerbosityThreshold = threshold;
}

/*!
    Gets the verbosity threshold for messages printed on the log file.

    Logging messages which are less severe than the verbosity threshold
    will be ignored; logging messages which have severity equal to the
    verbosity threshold or higher will be printed.
*/
log::Level Logger::getFileVerbosity()
{
    return m_fileVerbosityThreshold;
}

/*!
    Sets the context used when printing the messages.

    \param context is the context used when printing the messages.
*/
void Logger::setContext(const std::string &context)
{
    if (m_bufferSettings->context == context) {
        return;
    }

    m_buffer->flush(true);

    m_bufferSettings->context = context;
}

/*!
    Gets the context used when printing the messages.

    \result The the context used when printing the messages.
*/
std::string Logger::getContext()
{
    return m_bufferSettings->context;
}

/*!
    Sets the indentation level used when printing the messages.

    \param delta is the relative indentation level
*/
void Logger::setIndentation(int delta)
{
    if (delta == 0) {
        return;
    }

    m_buffer->flush(true);

    m_bufferSettings->indentation += delta;
}

/*!
    Gets the indentation level used when printing the messages.

    \result The the indentation level used when printing the messages.
*/
int Logger::getIndentation()
{
    return m_bufferSettings->indentation;
}

/*!
    Gets the name of the logger

    \result The name of the logger.
*/
std::string Logger::getName() const
{
    return m_name;
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
    \param severity is the severity of the line that will be printed
*/
void Logger::println(const std::string &line, log::Level severity)
{
    print(line + '\n', severity);
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
    \param severity is the severity of the line that will be printed
    \param visibility is the visibility of the line that will be printed
*/
void Logger::println(const std::string &line, log::Level severity, log::Visibility visibility)
{
    print(line + '\n', severity, visibility);
}

/*!
    Prints a message in the log.

    \param message is the message to be printed
*/
void Logger::print(const std::string &message)
{
    print(message, getDefaultSeverity(), getDefaultVisibility());
}

/*!
    Prints a message in the log.

    \param message is the message to be printed
    \param severity is the severity of the line that will be printed
*/
void Logger::print(const std::string &message, log::Level severity)
{
    print(message, severity, getDefaultVisibility());
}

/*!
    Prints a message in the log.

    \param message is the message to be printed
    \param visibility is the visibility of the line that will be printed
*/
void Logger::print(const std::string &message, log::Visibility visibility)
{
    print(message, getDefaultSeverity(), visibility);
}

/*!
    Prints a message in the log.

    \param message is the message to be printed
    \param severity is the severity of the line to be printed
    \param visibility is the visibility of the line that will be printed
*/
void Logger::print(const std::string &message, log::Level severity, log::Visibility visibility)
{
    print<const std::string &>(message, severity, visibility);
}

/*!
    \class LoggerManager
    \ingroup common_logger
    \brief Manager for the loggers.

    This class implements a manager for the loggers. The manager allows the
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
    m_mode(log::MODE_SEPARATE)
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

    // Destroy all the loggers
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
    Returns an instance of default logger.

    The default severity of the messages will be set to the specified level,
    if no default severity is specified, the default severity will remain
    unaltered.

    The default visibility of the messages will be set to the specified level,
    if no default visibility is specified, the default visibility will remain
    unaltered.

    \param defaultSeverity is the default severity of the messages.
    \param defaultVisibility is the default visibility of the messages
    \result An instance of the default specified logger.
*/
Logger & LoggerManager::cout(log::Level defaultSeverity, log::Visibility defaultVisibility)
{
    return cout(m_defaultName, defaultSeverity, defaultVisibility);
}

/*!
    Returns an instance of the specified logger.

    This function returns an instance of the specified logger. If the logger
    does not exists a new instance will be created.

    The default severity of the messages will be set to the specified level,
    if no default severity is specified, the default severity will remain
    unaltered.

    The default visibility of the messages will be set to the specified level,
    if no default visibility is specified, the default visibility will remain
    unaltered.

    \param name is the name of the logger
    \param defaultSeverity is the default severity of the messages.
    \param defaultVisibility is the default visibility of the messages
    \result An instance of the specified logger.
*/
Logger & LoggerManager::cout(const std::string &name, log::Level defaultSeverity, log::Visibility defaultVisibility)
{
    // Get the logger
    if (m_loggers.count(name) == 0) {
        if (!isInitialized()) {
            setMode(log::MODE_SEPARATE);
        }

        if (name == m_defaultName) {
            _create(name, false, m_defaultDirectory, 1, 0);
        } else {
            create(name, false, m_defaultDirectory, 1, 0);
        }
    }

    Logger &logger = *(m_loggers.at(name));

    // Set default properties
    if (defaultSeverity != log::Level::NOTSET) {
        logger.setDefaultSeverity(defaultSeverity);
    }

    if (defaultVisibility != log::VISIBILITY_NOTSET) {
        logger.setDefaultVisibility(defaultVisibility);
    }

    // Return the logger
    return logger;
}

/*!
    Returns an instance of default logger.

    The default severity of the messages will be set to the CRITICAL level.

    The default visibility of the messages will be set to the specified level,
    if no default visibility is specified, the default visibility will remain
    unaltered.

    \param defaultVisibility is the default visibility of the messages
    \result An instance of the default specified logger.
*/
Logger & LoggerManager::critical(log::Visibility defaultVisibility)
{
    return critical(m_defaultName, defaultVisibility);
}

/*!
    Returns an instance of the specified logger.

    This function returns an instance of the specified logger. If the logger
    does not exists a new instance will be created.

    The default severity of the messages will be set to the CRITICAL level.

    The default visibility of the messages will be set to the specified level,
    if no default visibility is specified, the default visibility will remain
    unaltered.

    \param name is the name of the logger
    \param defaultVisibility is the default visibility of the messages
    \result An instance of the specified logger.
*/
Logger & LoggerManager::critical(const std::string &name, log::Visibility defaultVisibility)
{
    return cout(name, log::Level::CRITICAL, defaultVisibility);
}

/*!
    Returns an instance of default logger.

    The default severity of the messages will be set to the ERROR level.

    The default visibility of the messages will be set to the specified level,
    if no default visibility is specified, the default visibility will remain
    unaltered.

    \param defaultVisibility is the default visibility of the messages
    \result An instance of the default specified logger.
*/
Logger & LoggerManager::error(log::Visibility defaultVisibility)
{
    return error(m_defaultName, defaultVisibility);
}

/*!
    Returns an instance of the specified logger.

    This function returns an instance of the specified logger. If the logger
    does not exists a new instance will be created.

    The default severity of the messages will be set to the ERROR level.

    The default visibility of the messages will be set to the specified level,
    if no default visibility is specified, the default visibility will remain
    unaltered.

    \param name is the name of the logger
    \param defaultVisibility is the default visibility of the messages
    \result An instance of the specified logger.
*/
Logger & LoggerManager::error(const std::string &name, log::Visibility defaultVisibility)
{
    return cout(name, log::Level::ERROR, defaultVisibility);
}

/*!
    Returns an instance of default logger.

    The default severity of the messages will be set to the WARNING level.

    The default visibility of the messages will be set to the specified level,
    if no default visibility is specified, the default visibility will remain
    unaltered.

    \param defaultVisibility is the default visibility of the messages
    \result An instance of the default specified logger.
*/
Logger & LoggerManager::warning(log::Visibility defaultVisibility)
{
    return warning(m_defaultName, defaultVisibility);
}

/*!
    Returns an instance of the specified logger.

    This function returns an instance of the specified logger. If the logger
    does not exists a new instance will be created.

    The default severity of the messages will be set to the WARNING level.

    The default visibility of the messages will be set to the specified level,
    if no default visibility is specified, the default visibility will remain
    unaltered.

    \param name is the name of the logger
    \param defaultVisibility is the default visibility of the messages
    \result An instance of the specified logger.
*/
Logger & LoggerManager::warning(const std::string &name, log::Visibility defaultVisibility)
{
    return cout(name, log::Level::WARNING, defaultVisibility);
}

/*!
    Returns an instance of default logger.

    The default severity of the messages will be set to the INFO level.

    The default visibility of the messages will be set to the specified level,
    if no default visibility is specified, the default visibility will remain
    unaltered.

    \param defaultVisibility is the default visibility of the messages
    \result An instance of the default specified logger.
*/
Logger & LoggerManager::info(log::Visibility defaultVisibility)
{
    return info(m_defaultName, defaultVisibility);
}

/*!
    Returns an instance of the specified logger.

    This function returns an instance of the specified logger. If the logger
    does not exists a new instance will be created.

    The default severity of the messages will be set to the INFO level.

    The default visibility of the messages will be set to the specified level,
    if no default visibility is specified, the default visibility will remain
    unaltered.

    \param name is the name of the logger
    \param defaultVisibility is the default visibility of the messages
    \result An instance of the specified logger.
*/
Logger & LoggerManager::info(const std::string &name, log::Visibility defaultVisibility)
{
    return cout(name, log::Level::INFO, defaultVisibility);
}

/*!
    Returns an instance of default logger.

    The default severity of the messages will be set to the DEBUG level.

    The default visibility of the messages will be set to the specified level,
    if no default visibility is specified, the default visibility will remain
    unaltered.

    \param defaultVisibility is the default visibility of the messages
    \result An instance of the default specified logger.
*/
Logger & LoggerManager::debug(log::Visibility defaultVisibility)
{
    return debug(m_defaultName, defaultVisibility);
}

/*!
    Returns an instance of the specified logger.

    This function returns an instance of the specified logger. If the logger
    does not exists a new instance will be created.

    The default severity of the messages will be set to the DEBUG level.

    The default visibility of the messages will be set to the specified level,
    if no default visibility is specified, the default visibility will remain
    unaltered.

    \param name is the name of the logger
    \param defaultVisibility is the default visibility of the messages
    \result An instance of the specified logger.
*/
Logger & LoggerManager::debug(const std::string &name, log::Visibility defaultVisibility)
{
    return cout(name, log::Level::DEBUG, defaultVisibility);
}

/*!
    Initializes the log manager.

    \param mode is the mode that will be set
    \param reset if true the log files will be reset
    \param nProcesses is the total number of processes in the communicator
    \param rank is the parallel rank in the communicator
*/
void LoggerManager::initialize(log::Mode mode, bool reset,
                            int nProcesses, int rank)
{
    initialize(mode, m_defaultName, reset, m_defaultDirectory, nProcesses, rank);
}

/*!
    Initializes the log manager.

    \param mode is the mode that will be set
    \param reset if true the log files will be reset
    \param directory is the default directory for saving the log files
    \param nProcesses is the total number of processes in the communicator
    \param rank is the parallel rank in the communicator
*/
void LoggerManager::initialize(log::Mode mode, bool reset, const std::string &directory,
                            int nProcesses, int rank)
{
    initialize(mode, m_defaultName, reset, directory, nProcesses, rank);
}

/*!
    Initializes the log manager.

    \param mode is the mode that will be set
    \param name is the name for the default logger
    \param reset if true the log files will be reset
    \param directory is the default directory for saving the log files
    \param nProcesses is the total number of processes in the communicator
    \param rank is the parallel rank in the communicator
*/
void LoggerManager::initialize(log::Mode mode, const std::string &name, bool reset,
                            const std::string &directory,
                            int nProcesses, int rank)
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
    _create(m_defaultName, reset, directory, nProcesses, rank);
}

/*!
    Creates a new logger.

    \param name is the name for the logger
    \param reset if true the log files will be reset
    \param nProcesses is the total number of processes in the communicator
    \param rank is the parallel rank in the communicator
*/
void LoggerManager::create(const std::string &name, bool reset,
                           int nProcesses, int rank)
{
    create(name, reset, m_defaultDirectory, nProcesses, rank);
}

/*!
    Creates a new logger.

    \param name is the name for the logger
    \param reset if true the log files will be reset
    \param directory is the directory for saving the log files
    \param nProcesses is the total number of processes in the communicator
    \param rank is the parallel rank in the communicator
*/
void LoggerManager::create(const std::string &name, bool reset,
                           const std::string &directory,
                           int nProcesses, int rank)
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
    if (m_mode == log::MODE_SEPARATE) {
        _create(name, reset, directory, nProcesses, rank);
    } else {
        Logger &defaultLogger = cout(m_defaultName);
        _create(name, defaultLogger.m_buffer);
    }
}

/*!
    Destroys the specified logger.

    The logger can be shared among different users, a logger is destroyed
    only if it has no users.

    \param name is the name of the logger
    \param force controls if the logger will be destroyed also if it still
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
    'initialize', implicit initialization happens the first time a
    logger is created.

    \result Returns true if the logger manager has been initialized, false
    otherwise.
*/
bool LoggerManager::isInitialized() const
{
    return (m_loggers.size() > 0);
}

/*!
    Set the operational mode of the logger manager.

    There are two possible modes:
     - when the manager operates in SEPARATE mode, each logger will write
       the messages in its own file;
     - when the manager operates in COMBINED mode, all the loggers will write
       the messages on the same file.

    The mode has to be set before creating the first logger.

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
    Get the operational mode of the logger manager.

    \result The operational mode of the logger manager.
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
    \param nProcesses is the total number of processes in the communicator
    \param rank is the parallel rank in the communicator
*/
void LoggerManager::_create(const std::string &name, bool reset,
                            const std::string &directory,
                            int nProcesses, int rank)
{
    // Get the file path
    FileHandler fileHandler;
    fileHandler.setDirectory(directory);
    fileHandler.setName(name);
    fileHandler.setAppendix("log");
    fileHandler.setParallel(nProcesses > 1);
    if (nProcesses > 1) {
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

    // Create the buffer
    std::shared_ptr<LoggerBuffer> buffer = std::make_shared<LoggerBuffer>(nProcesses, rank, 256);
    buffer->setFileStream(&fileStream);
    buffer->setConsoleStream(&consoleStream);

    // Create the logger
    _create(name, buffer);
}

/*!
    Internal function to create a logger.

    \param name is the name for the logger
    \param buffer is the buffer that will be used for the logger
*/
void LoggerManager::_create(const std::string &name, std::shared_ptr<LoggerBuffer> &buffer)
{
    m_loggers[name]     = std::unique_ptr<Logger>(new Logger(name,  buffer));
    m_loggerUsers[name] = 1;
}

/*!
    Sets the verbosity threshold for messages printed on the console and on
    the log file.

    Logging messages which are less severe than the verbosity threshold
    will be ignored; logging messages which have severity equal to the
    verbosity threshold or higher will be printed.

    \param threshold is the verbosity threshold
*/
void LoggerManager::setVerbosities(log::Level threshold)
{
    for (auto &entry : m_loggers) {
        entry.second->setVerbosities(threshold);
    }
}

/*!
    Sets the verbosity threshold for messages printed on the console.

    Logging messages which are less severe than the verbosity threshold
    will be ignored; logging messages which have severity equal to the
    verbosity threshold or higher will be printed.

    \param threshold is the verbosity threshold
*/
void LoggerManager::setConsoleVerbosity(log::Level threshold)
{
    for (auto &entry : m_loggers) {
        entry.second->setConsoleVerbosity(threshold);
    }
}

/*!
    Sets the verbosity threshold for messages printed on the log file.

    Logging messages which are less severe than the verbosity threshold
    will be ignored; logging messages which have severity equal to the
    verbosity threshold or higher will be printed.

    \param threshold is the verbosity threshold
*/
void LoggerManager::setFileVerbosity(log::Level threshold)
{
    for (auto &entry : m_loggers) {
        entry.second->setFileVerbosity(threshold);
    }
}

/*!
    Gets the default logger name.

    \result The default logger name.
*/
std::string LoggerManager::getDefaultName() const
{
    return m_defaultName;
}

/*!
    Gets the default logger directory.

    \result The default logger directory.
*/
std::string LoggerManager::getDefaultDirectory() const
{
    return m_defaultDirectory;
}

// Logger global functions
namespace log {
    /*!
        \enum Mode

        Defines the operational mode of the logger manager.

        \var Mode COMBINED
        In this mode all the loggers will write the messages on the same file.

        \var Mode SEPARATE
        In this mode each logger will write the messages on its own file.
    */

    /*!
        \enum Level

        Define the logging levels.

        \var Level CRITICAL
        This level identifies messages categorized as critical errors.

        \var Level ERROR
        This level identifies messages categorized as errors.

        \var Level WARNING
        This level identifies messages categorized as warnings.

        \var Level INFO
        This level identifies messages categorized as information.

        \var Level NORMAL
        This level identifies messages categorized as information. The usage of
        this level is deprecated as it as been superseded by INFO.

        \var Level DEBUG
        This level identifies messages categorized as debug.

        \var Level QUIET
        This level identifies messages that will not be written. The usage of
        this level is deprecated, logging can be disabled using the functions
        disable, disableConsole and disableFile.

        \var Level NOTSET
        The level specified an undefined level.
    */

    /*!
        \typedef Severity

        Defines the severity levels that can be associated with a line.
        Severity is defined as a typedef of Level.
    */

    /*!
        \typedef Priority

        Defines the priority levels that can be associated with a line.
        The usage of this enum is deprecated as it as been superseded by
        Severity.
    */

    /*!
        \typedef Verbosity

        Defines the verbosity levels that can be associated with a logger.
        The usage of this enum is deprecated as it as been superseded by
        Level.
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

        The default severity of the messages will be set to the specified
        level, if no default severity is specified, the default severity
        will remain unaltered.

        The default visibility of the messages will be set to the specified
        level, if no default visibility is specified, the default visibility
        will remain unaltered.

        \param defaultSeverity is the default severity of the messages.
        \result An instance of the default logger.
    */
    Logger & cout(log::Level defaultSeverity, log::Visibility defaultVisibility)
    {
        return manager().cout(defaultSeverity, defaultVisibility);
    }

    /*!
        Returns an instance of the specified logger.

        The default severity of the messages will be set to the specified
        level, if no default severity is specified, the default severity
        will remain unaltered.

        The default visibility of the messages will be set to the specified
        level, if no default visibility is specified, the default visibility
        will remain unaltered.

        \param name is the name of the logger
        \param defaultSeverity is the default severity of the messages.
        \param defaultVisibility is the default visibility of the messages
        \result An instance of the specified logger.
    */
    Logger & cout(const std::string &name, log::Level defaultSeverity, log::Visibility defaultVisibility)
    {
        return manager().cout(name, defaultSeverity, defaultVisibility);
    }

    /*!
        Returns an instance of the default logger.

        The default severity of the messages will be set to the CRITICAL level.

        The default visibility of the messages will be set to the specified
        level, if no default visibility is specified, the default visibility
        will remain unaltered.

        \result An instance of the default logger.
    */
    Logger & critical(log::Visibility defaultVisibility)
    {
        return manager().critical(defaultVisibility);
    }

    /*!
        Returns an instance of the specified logger.

        The default severity of the messages will be set to the CRITICAL level.

        The default visibility of the messages will be set to the specified
        level, if no default visibility is specified, the default visibility
        will remain unaltered.

        \param name is the name of the logger
        \param defaultVisibility is the default visibility of the messages
        \result An instance of the specified logger.
    */
    Logger & critical(const std::string &name, log::Visibility defaultVisibility)
    {
        return manager().critical(name, defaultVisibility);
    }

    /*!
        Returns an instance of the default logger.

        The default severity of the messages will be set to the ERROR level.

        The default visibility of the messages will be set to the specified
        level, if no default visibility is specified, the default visibility
        will remain unaltered.

        \param defaultVisibility is the default visibility of the messages
        \result An instance of the default logger.
    */
    Logger & error(log::Visibility defaultVisibility)
    {
        return manager().error(defaultVisibility);
    }

    /*!
        Returns an instance of the specified logger.

        The default severity of the messages will be set to the ERROR level.

        The default visibility of the messages will be set to the specified
        level, if no default visibility is specified, the default visibility
        will remain unaltered.

        \param name is the name of the logger
        \param defaultVisibility is the default visibility of the messages
        \result An instance of the specified logger.
    */
    Logger & error(const std::string &name, log::Visibility defaultVisibility)
    {
        return manager().error(name, defaultVisibility);
    }

    /*!
        Returns an instance of the default logger.

        The default severity of the messages will be set to the WARNING level.

        The default visibility of the messages will be set to the specified
        level, if no default visibility is specified, the default visibility
        will remain unaltered.

        \param defaultVisibility is the default visibility of the messages
        \result An instance of the default logger.
    */
    Logger & warning(log::Visibility defaultVisibility)
    {
        return manager().warning(defaultVisibility);
    }

    /*!
        Returns an instance of the specified logger.

        The default severity of the messages will be set to the WARNING level.

        The default visibility of the messages will be set to the specified
        level, if no default visibility is specified, the default visibility
        will remain unaltered.

        \param name is the name of the logger
        \param defaultVisibility is the default visibility of the messages
        \result An instance of the specified logger.
    */
    Logger & warning(const std::string &name, log::Visibility defaultVisibility)
    {
        return manager().warning(name, defaultVisibility);
    }

    /*!
        Returns an instance of the default logger.

        The default severity of the messages will be set to the INFO level.

        The default visibility of the messages will be set to the specified
        level, if no default visibility is specified, the default visibility
        will remain unaltered.

        \param defaultVisibility is the default visibility of the messages
        \result An instance of the default logger.
    */
    Logger & info(log::Visibility defaultVisibility)
    {
        return manager().info(defaultVisibility);
    }

    /*!
        Returns an instance of the specified logger.

        The default severity of the messages will be set to the INFO level.

        The default visibility of the messages will be set to the specified
        level, if no default visibility is specified, the default visibility
        will remain unaltered.

        \param name is the name of the logger
        \param defaultVisibility is the default visibility of the messages
        \result An instance of the specified logger.
    */
    Logger & info(const std::string &name, log::Visibility defaultVisibility)
    {
        return manager().info(name, defaultVisibility);
    }

    /*!
        Returns an instance of the default logger.

        The default severity of the messages will be set to the DEBUG level.

        The default visibility of the messages will be set to the specified
        level, if no default visibility is specified, the default visibility
        will remain unaltered.

        \param defaultVisibility is the default visibility of the messages
        \result An instance of the default logger.
    */
    Logger & debug(log::Visibility defaultVisibility)
    {
        return manager().debug(defaultVisibility);
    }

    /*!
        Returns an instance of the specified logger.

        The default severity of the messages will be set to the DEBUG level.

        The default visibility of the messages will be set to the specified
        level, if no default visibility is specified, the default visibility
        will remain unaltered.

        \param name is the name of the logger
        \param defaultVisibility is the default visibility of the messages
        \result An instance of the specified logger.
    */
    Logger & debug(const std::string &name, log::Visibility defaultVisibility)
    {
        return manager().debug(name, defaultVisibility);
    }

    // Manipulators global functions

    /*!
        Set the context of the output for the specified logger.

        \param logger is a reference pointing to the logger
        \param context is the context of the output
        \result A reference pointing to the logger received in input.
    */
    Logger & setContext(Logger &logger, const std::string &context)
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
        Set the default severity of the messageg for the specified logger.

        \param logger is a reference pointing to the logger
        \param severity is the default severity of the messages
        \result A reference pointing to the logger received in input.
    */
    Logger & setDefaultSeverity(Logger &logger, const log::Level &severity)
    {
        logger.setDefaultSeverity(severity);

        return logger;
    }

    /*!
        Returns a logger manipulator that allows to change the default severity
        of the messages.

        \param severity is the default severity of the messages
        \result A logger manipulator that allows to change the severity of the
        output.
    */
    LoggerManipulator<log::Level> defaultSeverity(const log::Level &severity)
    {
        return LoggerManipulator<log::Level>(setDefaultSeverity, severity);
    }

    /*!
        Set the priority of the output for the specified logger.

        \param logger is a reference pointing to the logger
        \param priority is the priority of the output
        \result A reference pointing to the logger received in input.
    */
    Logger & setPriority(Logger &logger, const log::Priority &priority)
    {
        logger.setDefaultSeverity(priority);

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
        return LoggerManipulator<log::Priority>(setDefaultSeverity, priority);
    }

    /*!
        Set the default visibility of the messages.

        \param logger is a reference pointing to the logger
        \param visibility is the default visibility of the messages
        \result A reference pointing to the logger received in input.
    */
    Logger & setDefaultVisibility(Logger &logger, const log::Visibility &visibility)
    {
        logger.setDefaultVisibility(visibility);

        return logger;
    }

    /*!
        Returns a logger manipulator that allows to change the default
        visibility of the messages.

        \param visibility is the default visibility of the messages
        \result A logger manipulator that allows to change the priority of the
        output.
    */
    LoggerManipulator<log::Visibility> defaultVisibility(const log::Visibility &visibility)
    {
        return LoggerManipulator<log::Visibility>(setDefaultVisibility, visibility);
    }

    /*!
        Set the visibility of the messages.

        \param logger is a reference pointing to the logger
        \param visibility is the visibility of the messages
        \result A reference pointing to the logger received in input.
    */
    Logger & setVisibility(Logger &logger, const log::Visibility &visibility)
    {
        logger.setDefaultVisibility(visibility);

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
        return LoggerManipulator<log::Visibility>(setDefaultVisibility, visibility);
    }

    /*!
        Sets the verbosity thresholds for the messages printed both on the
        console and on the log file.

        \param logger is a reference pointing to the logger
        \param threshold is the verbosity threshold
        \result A reference pointing to the logger received in input.
    */
    Logger & setVerbosities(Logger &logger, const log::Level &threshold)
    {
        logger.setVerbosities(threshold);

        return logger;
    }

    /*!
        Returns a logger manipulator that allows to change the verbosity
        threshold for the messages printed both on the console and on the
        log file.

        \param threshold is the verbosity threshold
        \result A logger manipulator that allows to change the verbosity for
        the messages printed on the console.
    */
    LoggerManipulator<log::Level> verbosities(const log::Level &threshold)
    {
        return LoggerManipulator<log::Level>(setVerbosities, threshold);
    }

    /*!
        Sets the verbosity threshold for the messages printed on the console.

        \param logger is a reference pointing to the logger
        \param threshold is the verbosity threshold
        \result A reference pointing to the logger received in input.
    */
    Logger & setConsoleVerbosity(Logger &logger, const log::Level &threshold)
    {
        logger.setConsoleVerbosity(threshold);

        return logger;
    }

    /*!
        Returns a logger manipulator that allows to change the verbosity
        threshold for the messages printed on the console.

        \param threshold is the verbosity threshold
        \result A logger manipulator that allows to change the verbosity for
        the messages printed on the console.
    */
    LoggerManipulator<log::Level> consoleVerbosity(const log::Level &threshold)
    {
        return LoggerManipulator<log::Level>(setConsoleVerbosity, threshold);
    }

    /*!
        Sets the verbosity threshold for the messages printed on the file.

        \param logger is a reference pointing to the logger
        \param threshold is the verbosity threshold
        \result A reference pointing to the logger received in input.
    */
    Logger & setFileVerbosity(Logger &logger, const log::Level &threshold)
    {
        logger.setFileVerbosity(threshold);

        return logger;
    }

    /*!
        Returns a logger manipulator that allows to change the verbosity
        threshold for the messages printed on the file.

        \param threshold is the verbosity threshold
        \result A logger manipulator that allows to change the verbosity for
        the messages printed on the file.
    */
    LoggerManipulator<log::Level> fileVerbosity(const log::Level &threshold)
    {
        return LoggerManipulator<log::Level>(setFileVerbosity, threshold);
    }

    /*!
        Provides an overriding verbosity level which takes precedence over
        the logger’s own level. Its effect is to disable all logging calls
        of severity level and below.

        \param logger is a reference pointing to the logger
        \param level is the overriding verbosity level
        \result A reference pointing to the logger received in input.
    */
    Logger & disable(Logger &logger, const log::Level &level)
    {
        logger.disable(level);

        return logger;
    }

    /*!
        Returns a logger manipulator that allows to define an overriding
        verbosity level which takes precedence over the logger’s own level.
        The effect of this overriding verbosity level is to disable all
        logging calls of severity level and below.

        \param level is the overriding verbosity level
        \result A logger manipulator that allows to change the verbosity for
        the messages printed on the console.
    */
    LoggerManipulator<log::Level> disable(const log::Level &level)
    {
        return LoggerManipulator<log::Level>(disable, level);
    }

    /*!
        Provides an overriding verbosity level which takes precedence over
        the logger’s own console level. Its effect is to disable all logging
        calls of severity level and below.

        \param logger is a reference pointing to the logger
        \param level is the overriding verbosity level
        \result A reference pointing to the logger received in input.
    */
    Logger & disableConsole(Logger &logger, const log::Level &level)
    {
        logger.disableConsole(level);

        return logger;
    }

    /*!
        Returns a logger manipulator that allows to define an overriding
        verbosity level which takes precedence over the logger’s own console
        level. The effect of this overriding verbosity level is to disable
        all logging calls of severity level and below.

        \param level is the overriding verbosity level
        \result A logger manipulator that allows to change the verbosity for
        the messages printed on the console.
    */
    LoggerManipulator<log::Level> disableConsole(const log::Level &level)
    {
        return LoggerManipulator<log::Level>(disableConsole, level);
    }

    /*!
        Provides an overriding verbosity level which takes precedence over
        the logger’s own file level. Its effect is to disable all logging
        calls of severity level and below.

        \param logger is a reference pointing to the logger
        \param level is the overriding verbosity level
        \result A reference pointing to the logger received in input.
    */
    Logger & disableFile(Logger &logger, const log::Level &level)
    {
        logger.disableFile(level);

        return logger;
    }

    /*!
        Returns a logger manipulator that allows to define an overriding
        verbosity level which takes precedence over the logger’s own file
        level. The effect of this overriding verbosity level is to disable
        all logging calls of severity level and below.

        \param level is the overriding verbosity level
        \result A logger manipulator that allows to change the verbosity for
        the messages printed on the file.
    */
    LoggerManipulator<log::Level> disableFile(const log::Level &level)
    {
        return LoggerManipulator<log::Level>(disableFile, level);
    }

    /*!
        Sets the indentation level of the messages.

        \param logger is a reference pointing to the logger
        \param delta is the relative indentation level
        \result A reference pointing to the logger received in input.
    */
    Logger & setIndentation(Logger &logger, const int &delta)
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
    LoggerManipulator<int> indent(int delta)
    {
        return LoggerManipulator<int>(setIndentation, delta);
    }
}

}
