/*---------------------------------------------------------------------------*\
*
*  bitpit
*
*  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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

#include <stdexcept>
#include <libxml/encoding.h>
#include <libxml/parser.h>

#include "configuration.hpp"
#include "configuration_XML.hpp"

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
ConfigParser::ConfigParser(const std::string &root, const int &version)
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
ConfigParser::ConfigParser(const std::string &root, const int &version, bool multiSections)
    : Config(multiSections)
{
    reset(root, version);
}

/*!
    Resets the configuration parser

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
    Resets the configuration parser

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
    Read the specified configuration file.

    \param filename is the filename of the configuration file
    \param append controls if the configuration file will be appended to the
    current configuration or if the current configuration will be overwritten
    with the contents of the configuration file.
*/
void ConfigParser::read(const std::string &filename, bool append)
{
    // Macro to check API for match with the DLL we are using
    LIBXML_TEST_VERSION

    // Read the XML file
    xmlDoc *doc = xmlReadFile(filename.c_str(), NULL, 0);
    if (doc == nullptr) {
        throw std::runtime_error("Could not parse configuration file \"" + filename + "\"");
    }

    // Get the root element
    xmlNode * rootElement = xmlDocGetRootElement(doc);

    // check if the root name is the requeste one
    std::string rootName(reinterpret_cast<const char*>(rootElement->name));
    if (rootName != m_root) {
        throw std::runtime_error("The name of the root element is not \"" + m_root + "\"");
    }

    // Check if the version is supported
    const xmlChar *versionAttributeName = reinterpret_cast<const xmlChar *>("version");
    if (m_checkVersion && xmlHasProp(rootElement, versionAttributeName)) {
        xmlChar *versionValue = xmlGetProp(rootElement, versionAttributeName);
        std::string versionString((char *) versionValue);
        xmlFree(versionValue);

        int version;
        std::istringstream(versionString) >> version;

        if (version != m_version) {
            throw std::runtime_error("The version ofthe file is not not \"" + std::to_string(m_version) + "\"");
        }
    }

    // Load the options in the configuration file
    if (!append) {
        clear();
    }
    config::XML::readNode(rootElement->children, this);

    // Clean-up
    xmlFreeDoc(doc);
    xmlCleanupParser();
}

/*!
    Write the configuration to the specified file.

    \param filename is the filename where the configuration will be written to
*/
void ConfigParser::write(const std::string &filename) const
{
    int status;

    // Create a new XmlWriter for DOM tree, with no compression
    xmlTextWriterPtr writer = xmlNewTextWriterFilename(filename.c_str(), 0);
    if (writer == NULL) {
        throw std::runtime_error("Error creating the xml writer");
    }

    xmlTextWriterSetIndent(writer, 1);

    // Start the document
    status = xmlTextWriterStartDocument(writer, NULL, config::XML::DEFAULT_ENCODING.c_str(), NULL);
    if (status < 0) {
        throw std::runtime_error("Error at xmlTextWriterStartDocument");
    }

    // Start the root element
    xmlChar *elementName = config::XML::encodeString(m_root, config::XML::DEFAULT_ENCODING);
    status = xmlTextWriterStartElement(writer, BAD_CAST elementName);
    if (status < 0) {
        throw std::runtime_error("Error at xmlTextWriterStartElement");
    }

    // Add an attribute with version
    std::ostringstream versionStream;
    versionStream << m_version;

    xmlChar *versionAttr = config::XML::encodeString(versionStream.str(), config::XML::DEFAULT_ENCODING);
    status = xmlTextWriterWriteAttribute(writer, BAD_CAST "version", BAD_CAST versionAttr);
    if (status < 0) {
        throw std::runtime_error("Error at xmlTextWriterWriteAttribute");
    }

    // Write the configuration
    config::XML::writeNode(writer, this, config::XML::DEFAULT_ENCODING);

    // End section
    status = xmlTextWriterEndElement(writer);
    if (status < 0) {
        throw std::runtime_error("Error at xmlTextWriterEndElement");
    }

    // Close the document
    status = xmlTextWriterEndDocument(writer);
    if (status < 0) {
        throw std::runtime_error("Error at xmlTextWriterEndDocument");
    }

    // Write the XML
    xmlFreeTextWriter(writer);
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
    Returns the global configuration.

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
        Resets the root element name.

        \param name is the name of the root element
    */
    void reset(const std::string &name)
    {
        root.reset(name);
    }

    /*!
        Resets the root element name and version.

        \param name is the name of the root element
        \param version is the version
    */
    void reset(const std::string &name, int version)
    {
        root.reset(name, version);
    }

    /*!
        Read the specified configuration file.

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

        \param filename is the filename where the configuration will be written to
    */
    void write(const std::string &filename)
    {
        root.write(filename);
    }

}

}
