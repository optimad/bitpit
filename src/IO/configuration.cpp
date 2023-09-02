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


#include "configuration.hpp"
#include "configuration_XML.hpp"

#if BITPIT_ENABLE_RAPIDJSON
#include "configuration_JSON.hpp"
#endif

#include <stringUtils.hpp>

namespace bitpit {

/*!
    \class ConfigParser
    \ingroup Configuration
    \brief Configuration file parser

    This class implements a configuration file parser.
*/

/*
    Undefined versione
*/
const int ConfigParser::VERSION_UNDEFINED = -1;


/*!
    Construct a new configuration parser.

    \param root is the name of the root element
*/
ConfigParser::ConfigParser(const std::string &root)
    : Config(false)
{
    reset(root);
}

/*!
    Construct a new configuration parser.

    \param root is the name of the root element
    \param multiSections if set to true the configuration parser will allow
    multiple sections with the same name
*/
ConfigParser::ConfigParser(const std::string &root, bool multiSections)
    : Config(multiSections)
{
    reset(root);
}

/*!
    Construct a new configuration parser.

    \param root is the name of the root element
    \param version is the required version
*/
ConfigParser::ConfigParser(const std::string &root, int version)
    : Config(false)
{
    reset(root, version);
}

/*!
    Construct a new configuration parser.

    \param root is the name of the root element
    \param version is the required version
    \param multiSections if set to true the configuration parser will allow
    multiple sections with the same name
*/
ConfigParser::ConfigParser(const std::string &root, int version, bool multiSections)
    : Config(multiSections)
{
    reset(root, version);
}

/*!
    Resets the configuration parser. MultiSection property set in construction
    is not affected.

    \param root is the name of the root element
*/
void ConfigParser::reset(const std::string &root)
{
    m_root         = root;
    m_checkVersion = false;
    m_version      = VERSION_UNDEFINED;

    clear();
}

/*!
    Resets the configuration parser. MultiSections property set in construction
    is not affected.

    \param root is the name of the root element
    \param version is the version
*/
void ConfigParser::reset(const std::string &root, int version)
{
    m_root         = root;
    m_checkVersion = true;
    m_version      = version;

    clear();
}

/*!
    Resets the configuration parser. MultiSections property set in construction
    is forced to the new value given by the parameter multiSection

    \param root is the name of the root element
    \param version is the version
    \param multiSections if set to true the configuration parser will allow
    multiple sections with the same name
*/
void ConfigParser::reset(const std::string &root, int version, bool multiSections)
{
    m_root         = root;
    m_checkVersion = true;
    m_version      = version;

    m_multiSections = multiSections;

    clear();
}

/*!
    Read the specified configuration file.

    Configuration file can be either XML or JSON files (if JSON support was
    found at compile time).

    \param filename is the filename of the configuration file
    \param append controls if the configuration file will be appended to the
    current configuration or if the current configuration will be overwritten
    with the contents of the configuration file.
*/
void ConfigParser::read(const std::string &filename, bool append)
{
    // Processing not-append requests
    if (!append) {
        clear();
    }

    // Read the configuration
    std::string extension = "";
    std::size_t dotPosition = filename.find_last_of(".");
    if (dotPosition != std::string::npos) {
        extension = filename.substr(dotPosition + 1);
        extension = utils::string::trim(extension);
    }

    if (extension == "xml" || extension == "XML") {
        config::XML::readConfiguration(filename, m_root, m_checkVersion, m_version, this);
#if BITPIT_ENABLE_RAPIDJSON
    } else if (extension == "json" || extension == "JSON") {
        config::JSON::readConfiguration(filename, this);
#endif
    } else {
        throw std::runtime_error("ConfigParser::read - Unsupported file format");
    }
}

/*!
    Write the configuration to the specified file.

    Configuration file can be either XML or JSON files (if JSON support was
    found at compile time).


    \param filename is the filename where the configuration will be written to
*/
void ConfigParser::write(const std::string &filename) const
{
    // Write the configuration
    std::string extension = "";
    std::size_t dotPosition = filename.find_last_of(".");
    if (dotPosition != std::string::npos) {
        extension = filename.substr(dotPosition + 1);
        extension = utils::string::trim(extension);
    }

    if (extension == "xml" || extension == "XML") {
        config::XML::writeConfiguration(filename, m_root, m_version, this);
#if BITPIT_ENABLE_RAPIDJSON
    } else if (extension == "json" || extension == "JSON") {
        bool prettify = true;
        config::JSON::writeConfiguration(filename, this, prettify);
#endif
    } else {
        throw std::runtime_error("ConfigParser::write - Unsupported file format");
    }
}

/*!
    \class GlobalConfigParser
    \ingroup Configuration
    \brief Global configuration file parser

    This class implements a global configuration file parser.
*/

/*
    Initialize the global instance of the configuration file parser.
*/
std::unique_ptr<GlobalConfigParser> GlobalConfigParser::m_parser = nullptr;

/*
    Initialize the defualt name of the root element.
*/
const std::string GlobalConfigParser::DEFAULT_ROOT_NAME = "bitpit";

/*
    Initialize the default version of the configuration file.
*/
const int GlobalConfigParser::DEFAULT_VERSION = 1;

/*!
    Default constructor.
*/
GlobalConfigParser::GlobalConfigParser()
    : ConfigParser(DEFAULT_ROOT_NAME, DEFAULT_VERSION)
{
}

/*!
    Constructor a new parser.
*/
GlobalConfigParser::GlobalConfigParser(const std::string &name, int version)
    : ConfigParser(name, version)
{
}

/*!
    Constructor a new parser.
*/
GlobalConfigParser::GlobalConfigParser(const std::string &name, bool multiSections)
    : ConfigParser(name, DEFAULT_VERSION, multiSections)
{
}

/*!
    Constructor a new parser.
*/
GlobalConfigParser::GlobalConfigParser(const std::string &name, int version, bool multiSections)
    : ConfigParser(name, version, multiSections)
{
}

/*!
    Returns a global configuration with default options.
    \result The global instance of the configuration file parser.
*/
GlobalConfigParser & GlobalConfigParser::parser()
{
    if (!m_parser) {
        m_parser = std::unique_ptr<GlobalConfigParser>(new GlobalConfigParser());
    }

    return *m_parser;
}

// Routines for interacting with the global configuration file parser
namespace config {

    /*!
        Global instance of the configuration file parser.
    */
    GlobalConfigParser &root = GlobalConfigParser::parser();

    /*!
        Resets to custom root element name, all other options are forcefully
        reset to default. MultiSection property remains unaffected by the reset

        \param name is the name of the root element
    */
    void reset(const std::string &name)
    {
        root.reset(name);
    }

    /*!
        Resets to custom root element name and version. , all other options are
        forcefully reset to default. MultiSection property remains unaffected by
        the reset.

        \param name is the name of the root element
        \param version is the version
    */
    void reset(const std::string &name, int version)
    {
        root.reset(name, version);
    }

    /*!
        Resets to custom root element name, version and multiSections property

        \param name is the name of the root element
        \param version is the version
        \param multiSections boolean to control multiSections property of Config tree
    */
    void reset(const std::string &name, int version, bool multiSections)
    {
        root.reset(name, version, multiSections);
    }

    /*!
        Read the specified configuration file.

        Configuration file can be either XML or JSON files (if JSON support was
        found at compile time).

        \param filename is the filename of the configuration file
        \param append controls if the configuration file will be appended to the
        current configuration or if the current configuration will be overwritten
        with the contents of the configuration file.
    */
    void read(const std::string &filename, bool append)
    {
        root.read(filename, append);
    }

    /*!
        Write the configuration to the specified file.

        Configuration file can be either XML or JSON files (if JSON support was
        found at compile time).

        \param filename is the filename where the configuration will be written to
    */
    void write(const std::string &filename)
    {
        root.write(filename);
    }

}

}
