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

#include <libxml/encoding.h>
#include <libxml/parser.h>
#include <stdexcept>
#include <string>

#include "configuration_XML.hpp"

namespace bitpit {

namespace config {

namespace XML {

/*!
    Default encoding.
*/
const std::string DEFAULT_ENCODING = "ISO-8859-1";

/*!
    Read an XML node.

    \param root is the root node to be read
    \param config is the configuration where the options and section will be
    stored
*/
void readNode(xmlNodePtr root, Config *config)
{
    for (xmlNode *node = root; node; node = node->next) {
        if (node->type != XML_ELEMENT_NODE) {
            readNode(node->children, config);
            continue;
        }

        bool isSection = (xmlChildElementCount(node) != 0);
        std::string key(reinterpret_cast<const char*>(node->name));
        if (isSection) {
            Config::Section *section;
            if (!config->isMultiSectionsEnabled() && config->hasSection(key)) {
                section = &(config->getSection(key));
            } else {
                section = &(config->addSection(key));
            }
            readNode(node->children, section);
        } else {
            xmlChar *nodeContent = xmlNodeGetContent(node);
            std::string value(reinterpret_cast<const char*>(nodeContent));
            config->set(key, value);
            xmlFree(nodeContent);
        }
    }
}

/*!
    Write an XML node.

    \param rootNode is the root node to be written
    \param config is the configuration where the options and section will be
    read
*/
void writeNode(xmlTextWriterPtr writer, const Config *config, const std::string &encoding)
{
    int status;

    // Write the options
    for (const auto &entry : config->getOptions()) {
        const std::string &key   = entry.first;
        const std::string &value = entry.second;

        xmlChar *elementName = encodeString(key, encoding);
        xmlChar *elementText = encodeString(value, encoding);
        int status = xmlTextWriterWriteFormatElement(writer, BAD_CAST elementName, "%s", elementText);
        if (elementText) {
            xmlFree(elementText);
        }
        if (elementName) {
            xmlFree(elementName);
        }
        if (status < 0) {
            throw std::runtime_error("Error at xmlTextWriterWriteFormatElement");
        }
    }

    // Write the sections
    for (auto &entry : config->getSections()) {
        const std::string &key = entry.first;
        const Config::Section *section = entry.second.get();

        // Start the section
        xmlChar *elementName = encodeString(key, encoding);
        status = xmlTextWriterStartElement(writer, BAD_CAST elementName);
        if (elementName) {
            xmlFree(elementName);
        }
        if (status < 0) {
            throw std::runtime_error("Error at xmlTextWriterStartElement");
        }

        // Write the section
        writeNode(writer, section, encoding);

        // End section
        status = xmlTextWriterEndElement(writer);
        if (status < 0) {
            throw std::runtime_error("Error at xmlTextWriterEndElement");
        }
    }
}

/*!
    Encodes the specified string into UTF-8 for processing with libxml2 APIs.

    \param int is the string in a given encoding
    \param encoding isthe encoding used
    \result The converted UTF-8 string, or NULL in case of error.
*/
xmlChar * encodeString(const std::string &in, const std::string &encoding)
{
    if (in.empty()) {
        return nullptr;
    }

    xmlCharEncodingHandlerPtr handler = xmlFindCharEncodingHandler(encoding.c_str());
    if (!handler) {
        return nullptr;
    }

    int size = (int) in.length() + 1;
    int out_size = size * 2 - 1;
    xmlChar *out = (unsigned char *) xmlMalloc((size_t) out_size);
    if (out == nullptr) {
        return nullptr;
    }

    int temp = size - 1;
    int ret = handler->input(out, &out_size, (const xmlChar *) in.c_str(), &temp);
    if ((ret < 0) || (temp - size + 1)) {
        xmlFree(out);
        out = 0;
    } else {
        out = (unsigned char *) xmlRealloc(out, out_size + 1);
        out[out_size] = 0;  /*null terminating out */
    }

    return out;
}

/*!
    Read the specified configuration file.

    \param filename is the filename of the configuration file
    \param rootname name of the root. If this name does not match xml doc root
        name return an error
    \param checkVersion boolean to enable the control of XML version
    \param version number of version to check. If checkVersion is true and
           version does not match the xml one, return an error
    \param rootConfig pointer to Config tree to store data parsed from document.
*/
void readConfiguration(const std::string &filename, const std::string &rootname, bool checkVersion,
                       int version, Config *rootConfig)
{
    if (!rootConfig) {
        throw std::runtime_error("XML::readConfiguration Null Config tree structure passed");
    }

    // Macro to check API for match with the DLL we are using
    LIBXML_TEST_VERSION

    // Read the XML file
    xmlDoc *doc = xmlReadFile(filename.c_str(), NULL, 0);
    if (doc == nullptr) {
        throw std::runtime_error("Could not parse XML configuration file \"" + filename + "\"");
    }

    // Get the root element
    xmlNode * rootElement = xmlDocGetRootElement(doc);

    // check if the root name is the requeste one
    std::string rootXMLName(reinterpret_cast<const char*>(rootElement->name));
    if (rootXMLName != rootname) {
        throw std::runtime_error("The name of the root XML element is not \"" + rootname + "\"");
    }

    // Check if the version is supported
    const xmlChar *versionAttributeName = reinterpret_cast<const xmlChar *>("version");
    if (checkVersion && xmlHasProp(rootElement, versionAttributeName)) {
        xmlChar *versionValue = xmlGetProp(rootElement, versionAttributeName);
        std::string versionString((char *) versionValue);
        xmlFree(versionValue);

        int versionXML;
        std::istringstream(versionString) >> versionXML;

        if (versionXML != version) {
            throw std::runtime_error("The version of the XML file is not not \"" + std::to_string(version) + "\"");
        }
    }

    readNode(rootElement->children, rootConfig);

    // Clean-up
    xmlFreeDoc(doc);
    xmlCleanupParser();
}

/*!
    Write the configuration to the specified file.

    \param filename is the filename where the configuration will be written to
    \param rootname name of the root to assign to the XML file.
    \param version number of version to assing to the XML file.
    \param rootConfig pointer to the Config tree to be written on file.
*/
void writeConfiguration(const std::string &filename, const std::string &rootname, int version,
                        const Config *rootConfig)
{
    if (!rootConfig) {
        throw std::runtime_error("XML::writeConfiguration Null Config tree structure passed");
    }

    int status;

    // Create a new XmlWriter for DOM tree, with no compression
    xmlTextWriterPtr writer = xmlNewTextWriterFilename(filename.c_str(), 0);
    if (writer == NULL) {
        throw std::runtime_error("Error creating the xml writer");
    }

    xmlTextWriterSetIndent(writer, 1);

    // Start the document
    status = xmlTextWriterStartDocument(writer, NULL, DEFAULT_ENCODING.c_str(), NULL);
    if (status < 0) {
        throw std::runtime_error("Error at xmlTextWriterStartDocument");
    }

    // Start the root element
    xmlChar *elementName = encodeString(rootname, DEFAULT_ENCODING);
    status = xmlTextWriterStartElement(writer, BAD_CAST elementName);
    if (elementName) {
        xmlFree(elementName);
    }
    if (status < 0) {
        throw std::runtime_error("Error at xmlTextWriterStartElement");
    }

    // Add an attribute with version
    std::ostringstream versionStream;
    versionStream << version;

    xmlChar *versionAttr = encodeString(versionStream.str(), DEFAULT_ENCODING);
    status = xmlTextWriterWriteAttribute(writer, BAD_CAST "version", BAD_CAST versionAttr);
    if (versionAttr) {
        xmlFree(versionAttr);
    }
    if (status < 0) {
        throw std::runtime_error("Error at xmlTextWriterWriteAttribute");
    }

    // Write the configuration
    writeNode(writer, rootConfig, DEFAULT_ENCODING);

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

}

}

}
