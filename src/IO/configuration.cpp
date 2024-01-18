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
#include "configuration_tree.hpp"

#include <stringUtils.hpp>

#include <sstream>

namespace bitpit {

/*!
    \class ConfigParser
    \ingroup Configuration
    \brief Configuration file parser

    This class implements a configuration file parser.
*/

/*
    Undefined version.
*/
const int ConfigParser::VERSION_UNDEFINED = -1;


/*!
    Construct a new configuration parser.

    \param rootName is the name of the root element
*/
ConfigParser::ConfigParser(const std::string &rootName)
    : Config(false)
{
    reset(rootName);
}

/*!
    Construct a new configuration parser.

    \param rootName is the name of the root element
    \param multiSections if set to true the configuration parser will allow
    multiple sections with the same name
*/
ConfigParser::ConfigParser(const std::string &rootName, bool multiSections)
    : Config(multiSections)
{
    reset(rootName);
}

/*!
    Construct a new configuration parser.

    \param rootName is the name of the root element
    \param version is the required version
*/
ConfigParser::ConfigParser(const std::string &rootName, int version)
    : Config(false)
{
    reset(rootName, version);
}

/*!
    Construct a new configuration parser.

    \param rootName is the name of the root element
    \param version is the required version
    \param multiSections if set to true the configuration parser will allow
    multiple sections with the same name
*/
ConfigParser::ConfigParser(const std::string &rootName, int version, bool multiSections)
    : Config(multiSections)
{
    reset(rootName, version);
}

/*!
    Resets the configuration parser.

    MultiSection property set in construction is not affected.

    \param rootName is the name of the root element
*/
void ConfigParser::reset(const std::string &rootName)
{
    m_rootName     = rootName;
    m_checkVersion = false;
    m_version      = VERSION_UNDEFINED;

    clear();
}

/*!
    Resets the configuration parser. MultiSections property set in construction
    is not affected.

    \param rootName is the name of the root element
    \param version is the version
*/
void ConfigParser::reset(const std::string &rootName, int version)
{
    m_rootName     = rootName;
    m_checkVersion = true;
    m_version      = version;

    clear();
}

/*!
    Resets the configuration parser. MultiSections property set in construction
    is forced to the new value given by the parameter multiSection

    \param rootName is the name of the root element
    \param version is the version
    \param multiSections if set to true the configuration parser will allow
    multiple sections with the same name
*/
void ConfigParser::reset(const std::string &rootName, int version, bool multiSections)
{
    m_rootName         = rootName;
    m_checkVersion = true;
    m_version      = version;

    m_multiSections = multiSections;

    clear();
}

/*!
    Read the specified configuration file.

    Configuration can be read either from XML or JSON files.

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

    config::SourceFormat fileFormat;
    if (extension == "xml" || extension == "XML") {
        fileFormat = config::SOURCE_FORMAT_XML;
    } else if (extension == "json" || extension == "JSON") {
        fileFormat = config::SOURCE_FORMAT_JSON;
    } else {
        throw std::runtime_error("ConfigParser::read - Unsupported file format");
    }

    bool checkFileVersion = false;
    if (fileFormat == config::SOURCE_FORMAT_XML) {
        checkFileVersion = m_checkVersion;
    }

    config::tree::readConfiguration(filename, fileFormat, m_rootName, checkFileVersion, m_version, this);
}

/*!
    Read the specified content.

    Configuration can be read from a content either in XML or JSON format.

    \param format is the format of the content
    \param content is the content that will be read
    \param append controls if the configuration file will be appended to the
    current configuration or if the current configuration will be overwritten
    with the contents of the configuration file.
*/
void ConfigParser::read(config::SourceFormat format, const std::string &content, bool append)
{
    // Processing not-append requests
    if (!append) {
        clear();
    }

    // Read the configuration
    bool checkFileVersion = false;
    if (format == config::SOURCE_FORMAT_XML) {
        checkFileVersion = m_checkVersion;
    }

    std::stringstream contentStream(content, std::ios_base::in | std::ios_base::out | std::ios_base::binary);

    config::tree::readConfiguration(contentStream, format, m_rootName, checkFileVersion, m_version, this);
}

/*!
    Write the configuration to the specified file.

    Configuration can be written either to XML or JSON files.

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

    config::SourceFormat fileFormat;
    if (extension == "xml" || extension == "XML") {
        fileFormat = config::SOURCE_FORMAT_XML;
    } else if (extension == "json" || extension == "JSON") {
        fileFormat = config::SOURCE_FORMAT_JSON;
    } else {
        throw std::runtime_error("ConfigParser::read - Unsupported file format");
    }

    config::tree::writeConfiguration(filename, fileFormat, m_rootName, m_version, this);
}

/*!
    Write the configuration to the specified string.

    Configuration can be written either in XML or JSON format.

    \param format is the format that will be used to write the configuration
    \param content on output will contain the configuration in the specified format
*/
void ConfigParser::write(config::SourceFormat format, std::string *content) const
{
    std::stringstream contentStream;
    config::tree::writeConfiguration(contentStream, format, m_rootName, m_version, this);

    *content = contentStream.str();
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
    Initialize the default name of the root element.
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

    \param rootName is the name of the root element
    \param version is the required version
*/
GlobalConfigParser::GlobalConfigParser(const std::string &rootName, int version)
    : ConfigParser(rootName, version)
{
}

/*!
    Constructor a new parser.

    \param rootName is the name of the root element
    \param multiSections if set to true the configuration parser will allow
    multiple sections with the same name
*/
GlobalConfigParser::GlobalConfigParser(const std::string &rootName, bool multiSections)
    : ConfigParser(rootName, DEFAULT_VERSION, multiSections)
{
}

/*!
    Constructor a new parser.

    \param rootName is the name of the root element
    \param version is the required version
    \param multiSections if set to true the configuration parser will allow
    multiple sections with the same name
*/
GlobalConfigParser::GlobalConfigParser(const std::string &rootName, int version, bool multiSections)
    : ConfigParser(rootName, version, multiSections)
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

        \param rootName is the name of the root element
    */
    void reset(const std::string &rootName)
    {
        root.reset(rootName);
    }

    /*!
        Resets to custom root element name and version. , all other options are
        forcefully reset to default. MultiSection property remains unaffected by
        the reset.

        \param rootName is the name of the root element
        \param version is the version
    */
    void reset(const std::string &rootName, int version)
    {
        root.reset(rootName, version);
    }

    /*!
        Resets to custom root element name, version and multiSections property

        \param rootName is the name of the root element
        \param version is the version
        \param multiSections boolean to control multiSections property of Config tree
    */
    void reset(const std::string &rootName, int version, bool multiSections)
    {
        root.reset(rootName, version, multiSections);
    }

    /*!
        Read the specified configuration file.

        Configuration can be read either from XML or JSON files.

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
        Read the specified configuration file.

        Configuration file can be either XML or JSON files (if JSON support was
        found at compile time).

        \param format is the format of the content
        \param content is the content that will be read
        \param append controls if the configuration file will be appended to the
        current configuration or if the current configuration will be overwritten
        with the contents of the configuration file.
    */
    void read(config::SourceFormat format, const std::string &content, bool append)
    {
        root.read(format, content, append);
    }

    /*!
        Write the configuration to the specified file.

        Configuration can be written either to XML or JSON files.

        \param filename is the filename where the configuration will be written to
    */
    void write(const std::string &filename)
    {
        root.write(filename);
    }

    /*!
        Write the configuration to the specified string.

        Configuration can be written either in XML or JSON format.

        \param format is the format that will be used to write the configuration
        \param content on output will contain the configuration in the specified format
    */
    void write(config::SourceFormat format, std::string *content)
    {
        root.write(format, content);
    }

}

}
