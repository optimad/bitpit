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
#ifndef __BITPIT_CONFIGURATION_TREE_HPP__
#define __BITPIT_CONFIGURATION_TREE_HPP__

#include "configuration_config.hpp"
#include "configuration_common.hpp"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace bitpit {

namespace config {

namespace tree {

template<typename Source>
void readConfiguration(Source &source, SourceFormat format, const std::string &rootName,
                       bool checkVersion, int version, Config *rootConfig);
template<typename Destination>
void writeConfiguration(Destination &destination, SourceFormat format, const std::string &rootName,
                        int version, const Config *rootConfig);

void importTree(boost::property_tree::ptree const &root, Config *config);
void exportTree(const Config &config, bool areArraysAllowed, boost::property_tree::ptree *root);

template<typename Destination>
void writeTree(Destination &destination, SourceFormat format, boost::property_tree::ptree &propertyTree);

void writeTree(const std::string &filename, SourceFormat format, boost::property_tree::ptree &propertyTree);

std::string readNodeValue(const boost::property_tree::ptree &node);

}

}

}

// Include template implementation
#include "configuration_tree.tpp"

#endif
