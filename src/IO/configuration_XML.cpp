/*---------------------------------------------------------------------------*\
*
*  bitpit
*
*  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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
    for (auto &entry : config->getOptions()) {
        std::string key   = entry.first;
        std::string value = entry.second;

        xmlChar *elementName = encodeString(key, encoding);
        xmlChar *elementText = encodeString(value, encoding);
        status = xmlTextWriterWriteFormatElement(writer, BAD_CAST elementName, "%s", elementText);
        if (status < 0) {
            throw std::runtime_error("Error at xmlTextWriterWriteFormatElement");
        }
    }

    // Write the sections
    for (auto &entry : config->getSections()) {
        std::string key = entry.first;
        const Config::Section *section = entry.second.get();

        // Start the section
        xmlChar *elementName = encodeString(key, encoding);
        int status = xmlTextWriterStartElement(writer, BAD_CAST elementName);
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

}

}

}
