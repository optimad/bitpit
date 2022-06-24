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
#ifndef __BITPIT_LOGGER_HPP__
#define __BITPIT_LOGGER_HPP__

#include "compiler.hpp"

#include <memory>
#include <ostream>
#include <unordered_map>
#include <vector>

#define BITPIT_DEBUG_COUT(...) BITPIT_OVERLOAD_CALL(BITPIT_DEBUG_COUT, __VA_ARGS__)
#if BITPIT_ENABLE_DEBUG
#define BITPIT_DEBUG_COUT_0()     bitpit::log::cout()
#define BITPIT_DEBUG_COUT_1(NAME) bitpit::log::cout(NAME)
#else
#define BITPIT_DEBUG_COUT_0()     while(0) bitpit::log::cout()
#define BITPIT_DEBUG_COUT_1(NAME) while(0) bitpit::log::cout(NAME)
#endif

namespace bitpit{

namespace log {
    enum Mode {
        MODE_UNDEFINED = -1,
        MODE_SEPARATE = 0,
        MODE_COMBINE,
        COMBINED = MODE_COMBINE,
        SEPARATE = MODE_SEPARATE
    };

    enum Level {
        QUIET = 60,
        CRITICAL = 50,
        ERROR = 40,
        WARNING = 30,
        INFO = 20,
        NORMAL = INFO,
        DEBUG = 10,
        NOTSET = 0
    };

    enum Visibility {
        VISIBILITY_NOTSET = -1,
        VISIBILITY_MASTER = 0,
        VISIBILITY_GLOBAL,
        MASTER = VISIBILITY_MASTER,
        GLOBAL = VISIBILITY_GLOBAL,
    };

    typedef Level Severity;
    typedef Level Priority;
    typedef Level Verbosity;
}

// Logger buffer
class LoggerBuffer : public std::streambuf
{

public:
    LoggerBuffer(std::size_t bufferSize = 256);
    LoggerBuffer(LoggerBuffer const &other);
    ~LoggerBuffer();

    bool isConsoleTimestampEnabled() const;
    void setConsoleEnabled(bool enabled);
    void setConsoleTimestampEnabled(bool enabled);
    void setConsoleStream(std::ostream *console);
    std::ostream & getConsoleStream();
    void setConsolePrefix(const std::string &prefix);
    std::string getConsolePrefix()	const;

    bool isFileTimestampEnabled() const;
    void setFileEnabled(bool enabled);
    void setFileTimestampEnabled(bool enabled);
    void setFileStream(std::ofstream *file);
    std::ofstream & getFileStream();
    void setFilePrefix(const std::string &prefix);
    std::string getFilePrefix() const;

    void setContext(const std::string &context);
    void setPadding(const std::string &padding);

    int flush(bool terminate);

private:
    std::vector<char> m_buffer;

    std::string m_context;
    std::string m_padding;

    bool m_consoleEnabled;
    bool m_consoleTimestampEnabled;
    std::ostream *m_console;
    std::string m_consolePrefix;

    bool m_fileEnabled;
    bool m_fileTimestampEnabled;
    std::ofstream *m_file;
    std::string m_filePrefix;

    const std::string getTimestamp() const;

    int_type overflow(int_type ch) override;
    int sync() override;

};

// Logger
class Logger : public std::ostream
{

friend class LoggerManager;

public:
    Logger(const std::string &name,
           std::ostream *consoleStream, std::ofstream *fileStream,
           int nProcesses = 1, int rank = 0);

    int getProcessorCount();
    int getRank();

    void setContext(const std::string &context);
    std::string getContext();

    void setIndentation(int delta);
    int getIndentation();

    void disable(log::Level = log::Level::CRITICAL);
    void disableConsole(log::Level = log::Level::CRITICAL);
    void disableFile(log::Level = log::Level::CRITICAL);

    void setDefaultSeverity(log::Level severity);
    log::Level getDefaultSeverity();

    BITPIT_DEPRECATED(void setPriority(log::Priority priority));
    BITPIT_DEPRECATED(log::Priority getPriority());

    void setDefaultVisibility(log::Visibility visibility);
    log::Visibility getDefaultVisibility();

    BITPIT_DEPRECATED(void setVisibility(log::Visibility visibility));
    BITPIT_DEPRECATED(log::Visibility getVisibility());

    void setVerbosities(log::Level threshold);
    void setTimestampEnabled(bool enabled);

    bool isConsoleTimestampEnabled() const;
    void setConsoleTimestampEnabled(bool enabled);
    void setConsoleStream(std::ostream *console);
    std::ostream & getConsoleStream();
    std::string getConsolePrefix();
    void setConsoleVerbosity(log::Level threshold);
    log::Level getConsoleVerbosity();

    bool isFileTimestampEnabled() const;
    void setFileTimestampEnabled(bool enabled);
    void setFileStream(std::ofstream *file);
    std::ofstream & getFileStream();
    std::string getFilePrefix();
    void setFileVerbosity(log::Level threshold);
    log::Level getFileVerbosity();

    std::string getName() const;

    void println(const std::string &message);
    void println(const std::string &message, log::Level severity);
    void println(const std::string &message, log::Visibility visibility);
    void println(const std::string &message, const log::Level severity, log::Visibility visibility);

    void print(const std::string &message);
    void print(const std::string &message, log::Level severity);
    void print(const std::string &message, log::Visibility visibility);
    void print(const std::string &message, log::Level severity, log::Visibility visibility);

private:
    std::string m_name;
    int m_nProcesses;
    int m_rank;
    LoggerBuffer m_buffer;

    int m_indentation;
    std::string m_context;

    log::Level m_defaultSeverity;
    log::Visibility m_defaultVisibility;

    log::Level m_consoleDisabledThreshold;
    log::Level m_consoleVerbosityThreshold;

    log::Level m_fileDisabledThreshold;
    log::Level m_fileVerbosityThreshold;

    Logger(Logger const&);

    void setConsoleEnabled(log::Level severity, log::Visibility visibility);
    void setFileEnabled(log::Level severity, log::Visibility visibility);

};

// Logger manager
class LoggerManager
{

public:
    static std::string BITPIT_LOG_NAME;
    static std::string BITPIT_LOG_DIRECTORY;

    static LoggerManager & manager();

    ~LoggerManager();

    Logger & cout(log::Level defualtSeverity = log::Level::NOTSET, log::Visibility defualtVisibility = log::VISIBILITY_NOTSET);
    Logger & cout(const std::string &name, log::Level defualtSeverity = log::Level::NOTSET, log::Visibility defualtVisibility = log::VISIBILITY_NOTSET);

    Logger & critical(log::Visibility visibility = log::VISIBILITY_NOTSET);
    Logger & critical(const std::string &name,log::Visibility visibility = log::VISIBILITY_NOTSET);

    Logger & error(log::Visibility visibility = log::VISIBILITY_NOTSET);
    Logger & error(const std::string &name, log::Visibility visibility = log::VISIBILITY_NOTSET);

    Logger & warning(log::Visibility visibility = log::VISIBILITY_NOTSET);
    Logger & warning(const std::string &name, log::Visibility visibility = log::VISIBILITY_NOTSET);

    Logger & info(log::Visibility visibility = log::VISIBILITY_NOTSET);
    Logger & info(const std::string &name, log::Visibility visibility = log::VISIBILITY_NOTSET);

    Logger & debug(log::Visibility visibility = log::VISIBILITY_NOTSET);
    Logger & debug(const std::string &name, log::Visibility visibility = log::VISIBILITY_NOTSET);

    void initialize(log::Mode mode, bool reset,
                    int nProcesses, int rank);

    void initialize(log::Mode mode, bool reset = false,
                    const std::string &directory = BITPIT_LOG_DIRECTORY,
                    int nProcesses = 1, int rank = 0);

    void initialize(log::Mode mode, const std::string &name,
                    bool reset = false, const std::string &directory = BITPIT_LOG_DIRECTORY,
                    int nProcesses = 1, int rank = 0);

    void create(const std::string &name, bool reset = false,
                int nProcesses = 1, int rank = 0);

    void create(const std::string &name, bool reset, const std::string &directory,
                int nProcesses = 1, int rank = 0);

    bool destroy(const std::string &name, bool force = false);

    bool exists(const std::string &name) const;

    bool isInitialized() const;

    bool setMode(log::Mode mode);
    log::Mode getMode() const;

    void setVerbosities(log::Level threshold);
    void setConsoleVerbosity(log::Level threshold);
    void setFileVerbosity(log::Level threshold);

    std::string getDefaultName() const;
    std::string getDefaultDirectory() const;

private:
    static std::unique_ptr<LoggerManager> m_manager;

    std::string m_defaultName;
    std::string m_defaultDirectory;
    log::Mode m_mode;

    std::unordered_map<std::string, std::unique_ptr<Logger>> m_loggers;
    std::unordered_map<std::string, int> m_loggerUsers;
    std::unordered_map<std::string, std::unique_ptr<std::ofstream>> m_fileStreams;

    LoggerManager();

    LoggerManager(LoggerManager const&) = delete;
    LoggerManager& operator=(LoggerManager const&) = delete;

    void _create(const std::string &name, bool reset, const std::string &directory, int nProcesses, int rank);
    void _create(const std::string &name, Logger &master);

};

/*!
    \ingroup common_logger
    \brief The namespace 'log' contains routines for interacting with the
    message logger.
*/
namespace log {

    // Generic global functions
    LoggerManager & manager();

    Logger & cout(log::Level defualtSeverity = log::Level::NOTSET, log::Visibility defualtVisibility = log::VISIBILITY_NOTSET);
    Logger & cout(const std::string &name, log::Level defualtSeverity = log::Level::NOTSET, log::Visibility defualtVisibility = log::VISIBILITY_NOTSET);

    Logger & critical(log::Visibility visibility = log::VISIBILITY_NOTSET);
    Logger & critical(const std::string &name, log::Visibility visibility = log::VISIBILITY_NOTSET);

    Logger & error(log::Visibility visibility = log::VISIBILITY_NOTSET);
    Logger & error(const std::string &name, log::Visibility visibility = log::VISIBILITY_NOTSET);

    Logger & warning(log::Visibility visibility = log::VISIBILITY_NOTSET);
    Logger & warning(const std::string &name, log::Visibility visibility = log::VISIBILITY_NOTSET);

    Logger & info(log::Visibility visibility = log::VISIBILITY_NOTSET);
    Logger & info(const std::string &name, log::Visibility visibility = log::VISIBILITY_NOTSET);

    Logger & debug(log::Visibility visibility = log::VISIBILITY_NOTSET);
    Logger & debug(const std::string &name, log::Visibility visibility = log::VISIBILITY_NOTSET);

    // Manipulators global functions

    /*!
        Struct that allows to manipulate a logger.
    */
    template<typename T>
    struct LoggerManipulator {
        Logger & (*f) (Logger &, const T &);
        T value;

        /*!
            Creates a new logger manipulator
        */
        LoggerManipulator(Logger & (*ff)(Logger &, const T &), const T & ss)
            : f(ff), value(ss) {

        }
    };

    /*!
        Applies a manipulator to a logger.
    */
    template<typename T>
    Logger & operator<<(Logger &logger, LoggerManipulator<T>&& m)
    {
        return m.f(logger, m.value);
    }

    Logger & setContext(Logger &logger, const std::string &context);
    LoggerManipulator<std::string> context(const std::string &context);

    Logger & setDefaultSeverity(Logger &logger, const log::Level &severity);
    LoggerManipulator<log::Level> defaultSeverity(const log::Level &severity);

    BITPIT_DEPRECATED(Logger & setPriority(Logger &logger, const log::Priority &priority));
    BITPIT_DEPRECATED(LoggerManipulator<log::Level> priority(const log::Priority &priority));

    BITPIT_DEPRECATED(Logger & setVisibility(Logger &logger, const log::Visibility &visibility));
    BITPIT_DEPRECATED(LoggerManipulator<log::Visibility> visibility(const log::Visibility &visibility));

    Logger & setDefaultVisibility(Logger &logger, const log::Visibility &visibility);
    LoggerManipulator<log::Visibility> defaultVisibility(const log::Visibility &visibility);

    Logger & setVerbosities(Logger &logger, const log::Level &threshold);
    LoggerManipulator<log::Level> verbosities(const log::Level &threshold);

    Logger & setConsoleVerbosity(Logger &logger, const log::Level &threshold);
    LoggerManipulator<log::Level> consoleVerbosity(const log::Level &threshold);

    Logger & setFileVerbosity(Logger &logger, const log::Level &threshold);
    LoggerManipulator<log::Level> fileVerbosity(const log::Level &threshold);

    Logger & disable(Logger &logger, const log::Level &verbosity = log::Level::CRITICAL);
    LoggerManipulator<log::Level> disable(const log::Level &verbosity = log::Level::CRITICAL);

    Logger & disableConsole(Logger &logger, const log::Level &verbosity = log::Level::CRITICAL);
    LoggerManipulator<log::Level> disableConsole(const log::Level &verbosity = log::Level::CRITICAL);

    Logger & disableFile(Logger &logger, const log::Level &verbosity = log::Level::CRITICAL);
    LoggerManipulator<log::Level> disableFile(const log::Level &verbosity = log::Level::CRITICAL);

    Logger & setIndentation(Logger &logger, const int &delta);
    LoggerManipulator<int> indent(int delta);

}

}

#endif
