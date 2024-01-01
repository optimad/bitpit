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

#include "configuration_tree.hpp"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <unordered_set>

namespace bitpit {

namespace config {

namespace tree {

/*!
    Read the specified configuration file.

    \param filename is the filename of the configuration file
    \param format is the format of the file
    \param rootname name of the root to assign to the configuration file
    \param checkVersion boolean to enable the control of the version
    \param version number of version to assign to the configuration file
    \param rootConfig pointer to the configuration tree to be written on file
    \param[in,out] rootConfig pointer to configuration tree to store data parsed from document.
*/
void readConfiguration(const std::string &filename, FileFormat format,
                       const std::string &rootname, bool checkVersion,
                       int version, Config *rootConfig)
{
    if (!rootConfig) {
        throw std::runtime_error("config::tree::readConfiguration Null Config tree structure passed");
    }

    // Read property tree
    boost::property_tree::ptree propertyTree;
    if (format == FILE_FORMAT_XML) {
        boost::property_tree::read_xml(filename, propertyTree);
    } else if (format == FILE_FORMAT_JSON) {
        boost::property_tree::read_json(filename, propertyTree);
    }

    // Check if a root node exists
    bool hasRootNode = hasRootSection(format);

    // Validate root node
    if (hasRootNode) {
        // Get the root node
        auto rootNode = propertyTree.begin();

        // Check if the root name matches the requested one
        if (rootNode->first != rootname) {
            throw std::runtime_error("The name of the root element is not \"" + rootname + "\"");
        }

        // Check if the version matches the requested one
        if (checkVersion) {
            try {
                int fileVersion = rootNode->second.get<int>("<xmlattr>.version");
                if (fileVersion != version) {
                    throw std::runtime_error("The version of the configuration file is not not \"" + std::to_string(version) + "\"");
                }
            } catch (boost::property_tree::ptree_bad_path &error) {
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
    Write the configuration to the specified file.

    \param filename is the filename where the configuration will be written to
    \param format is the format of the file
    \param rootname name of the root to assign to the configuration file
    \param version number of version to assign to the configuration file
    \param rootConfig pointer to the configuration tree to be written on file
*/
void writeConfiguration(const std::string &filename, FileFormat format,
                        const std::string &rootname, int version,
                        const Config *rootConfig)
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
        propertyTree.put<std::string>(rootname, "");

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

    // Write the file
    if (format == FILE_FORMAT_XML) {
        boost::property_tree::xml_writer_settings<std::string> settings(' ', 4);
        write_xml(filename, propertyTree, std::locale(), settings);
    } else if (format == FILE_FORMAT_JSON) {
        write_json(filename, propertyTree, std::locale(), true);
    }
}

/*!
    Import the specified tree.

    \param root is the root node to be imported
    \param config is the configuration where the options and section will be stored
*/
void importTree(boost::property_tree::ptree const &root, Config *config)
{
    for (boost::property_tree::ptree::value_type const &node : root) {
        std::string key = node.first;
        if (key == "<xmlattr>") {
            continue;
        } else if (key == "") {
            continue;
        }

        bool isSection = (node.second.size() != 0);
        if (isSection) {
            // Check if the section defines an array
            //
            // A node can define a single section or an array of sections. Arrays of sections are
            // stored as key/value pairs where the key is an empty string.
            bool isArray = true;
            for (boost::property_tree::ptree::value_type const &itemCandidate : node.second) {
                if (itemCandidate.first != "") {
                    isArray = false;
                    break;
                }
            }

            // Read section
            if (isArray) {
                if (!config->isMultiSectionsEnabled()) {
                    throw std::runtime_error("config::tree:readNode reading array but Config MultiSection disabled");
                }

                for (boost::property_tree::ptree::value_type const &item : node.second) {
                    Config::Section *section = &(config->addSection(key));
                    importTree(item.second, section);
                }
            } else {
                Config::Section *section;
                if (!config->isMultiSectionsEnabled() && config->hasSection(key)) {
                    section = &(config->getSection(key));
                } else {
                    section = &(config->addSection(key));
                }

                importTree(node.second, section);
            }
        } else {
            std::string value = readNodeValue(node.second);
            config->set(key, value);
        }
    }
}

/*!
    Export the specified tree.

    \param config is the configuration where the options and section will be read
    \param areArraysAllowed controls if sections with the same name can be grouped into arrays
    \param root is the root node to be written
*/
void exportTree(const Config &config, bool areArraysAllowed, boost::property_tree::ptree *root)
{
    // Write the options
    for (const auto &entry : config.getOptions()) {
        const std::string &key   = entry.first;
        const std::string &value = entry.second;
        root->add(key, value);
    }

    // Write the sections
    std::unordered_set<std::string> processedKeys;
    for (const auto &entry : config.getSections()) {
        const std::string &key = entry.first;
        if (areArraysAllowed && processedKeys.count(key) != 0) {
            continue;
        }

        // Create the section tree
        root->push_back(boost::property_tree::ptree::value_type(key, boost::property_tree::ptree()));
        auto &sectionChild = root->back();
        auto &sectionRoot = sectionChild.second;

        // Check if the section defines an array
        bool isArray = areArraysAllowed && (config.getSectionCount(key) > 1);
        if (isArray) {
            // Fill the section with the elements of the array
            //
            // Arrays of sections should be stored as key/value pairs where the key is an
            // empty string.
            for (const auto &itemEntry : config.getSections(key)) {
                // Create the item tree
                sectionRoot.push_back(boost::property_tree::ptree::value_type("", boost::property_tree::ptree()));
                auto &itemChild = sectionRoot.back();

                // Fill the section
                const Config::Section &itemSection = *itemEntry;
                exportTree(itemSection, areArraysAllowed, &(itemChild.second));
            }
        } else {
            // Fill the section
            const Config::Section &section = *(entry.second);
            exportTree(section, areArraysAllowed, &sectionRoot);
        }

        // The key has been processed
        if (areArraysAllowed) {
            processedKeys.insert(key);
        }
    }
}

/*!
    Read the value of the specified node.

    \param config is the configuration where the options and section will be read
    \result The value of the specified node.
*/
std::string readNodeValue(const boost::property_tree::ptree &node)
{
    std::string value = node.get_value("");

    std::string lowercaseValue = value;
    std::transform(lowercaseValue.begin(), lowercaseValue.end(), lowercaseValue.begin(), ::tolower);
    if (value == "true") {
        return std::to_string(1);
    } else if (value == "false") {
        return std::to_string(0);
    } else {
        return value;
    }
}

}

}

}
