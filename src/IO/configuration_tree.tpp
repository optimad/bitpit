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
#ifndef __BITPIT_CONFIGURATION_TREE_TPP__
#define __BITPIT_CONFIGURATION_TREE_TPP__

namespace bitpit {

namespace config {

namespace tree {

/*!
    Read the specified configuration from the given source.

    \param source is the source of the configuration file
    \param format is the format of the file
    \param rootName name of the root to assign to the configuration file
    \param checkVersion boolean to enable the control of the version
    \param version number of version to assign to the configuration file
    \param rootConfig pointer to the configuration tree to be written on file
    \param[in,out] rootConfig pointer to configuration tree to store data parsed from document.
*/
template<typename Source>
void readConfiguration(Source &source, SourceFormat format, const std::string &rootName,
                       bool checkVersion, int version, Config *rootConfig)
{
    if (!rootConfig) {
        throw std::runtime_error("config::tree::readConfiguration Null Config tree structure passed");
    }

    // Read property tree
    boost::property_tree::ptree propertyTree;
    if (format == SOURCE_FORMAT_XML) {
        boost::property_tree::read_xml(source, propertyTree);
    } else if (format == SOURCE_FORMAT_JSON) {
        boost::property_tree::read_json(source, propertyTree);
    }

    // Check if a root node exists
    bool hasRootNode = hasRootSection(format);

    // Validate root node
    if (hasRootNode) {
        // Get the root node
        auto rootNode = propertyTree.begin();

        // Check if the root name matches the requested one
        if (rootNode->first != rootName) {
            throw std::runtime_error("The name of the root element is not \"" + rootName + "\"");
        }

        // Check if the version matches the requested one
        if (checkVersion) {
            try {
                int fileVersion = rootNode->second.get<int>("<xmlattr>.version");
                if (fileVersion != version) {
                    throw std::runtime_error("The version of the configuration file is not not \"" + std::to_string(version) + "\"");
                }
            } catch (boost::property_tree::ptree_bad_path &error) {
                BITPIT_UNUSED(error);
                throw std::runtime_error("Unable to identify the version of the configuration");
            }
        }
    }

    // Read configuration
    const boost::property_tree::ptree *configTree;
    if (hasRootNode) {
        configTree = &(propertyTree.begin()->second);
    } else {
        configTree = &propertyTree;
    }

    importTree(*configTree, rootConfig);
}

/*!
    Write the specified configuration to the given destination.

    \param source is the destination of the configuration file
    \param format is the format of the file
    \param rootName name of the root to assign to the configuration file
    \param checkVersion boolean to enable the control of the version
    \param version number of version to assign to the configuration file
    \param rootConfig pointer to the configuration tree to be written on file
    \param[in,out] rootConfig pointer to configuration tree to store data parsed from document.
*/
template<typename Destination>
void writeConfiguration(Destination &destination, SourceFormat format, const std::string &rootName,
                        int version, const Config *rootConfig)
{
    if (!rootConfig) {
        throw std::runtime_error("config::tree::writeConfiguration Null Config tree structure passed");
    }

    // Check if a root node exists
    bool hasRootNode = hasRootSection(format);

    // Create an empty property tree
    boost::property_tree::ptree propertyTree;

    // Create the root node
    if (hasRootNode) {
        // Create the root
        propertyTree.put<std::string>(rootName, "");

        // Get the root node
        auto rootNode = propertyTree.begin();

        // Add an attribute with version
        rootNode->second.put<int>("<xmlattr>.version", version);
    }

    // Export the configuration into the property tree
    boost::property_tree::ptree *configTree;
    if (hasRootNode) {
        configTree = &(propertyTree.begin()->second);
    } else {
        configTree = &propertyTree;
    }

    bool areArraysAllowed = hasArraySupport(format);

    exportTree(*rootConfig, areArraysAllowed, configTree);

    // Write the tree
    writeTree(destination, format, propertyTree);
}

/*!
    Write the specified tree to the given destination.

    \param destination is the destination of the configuration file
    \param format is the format that will be used for writing the configuration
    \param propertyTree is the property tree that will be written
*/
template<typename Destination>
void writeTree(Destination &destination, SourceFormat format, boost::property_tree::ptree &propertyTree)
{
    if (format == SOURCE_FORMAT_XML) {
        boost::property_tree::xml_writer_settings<std::string> settings(' ', 4);
        write_xml(destination, propertyTree, settings);
    } else if (format == SOURCE_FORMAT_JSON) {
        write_json(destination, propertyTree, true);
    }
}

}

}

}

#endif
