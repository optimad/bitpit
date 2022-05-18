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
#ifndef __BITPIT_CONFIGURATION_JSON_HPP__
#define __BITPIT_CONFIGURATION_JSON_HPP__

#if HAS_RAPIDJSON_LIB
#include <rapidjson/document.h>

#include "configuration_config.hpp"

namespace bitpit {

namespace config {

namespace JSON {

void readConfiguration(const std::string &filename, Config *rootConfig);
void readBufferConfiguration(const std::string &source, Config *rootConfig);
void readNode(const std::string &key, const rapidjson::Value &value, Config *config);

void writeConfiguration(const std::string &filename, const Config *rootConfig, bool prettify = true);
void writeBufferConfiguration(std::string &source, const Config *rootConfig);
void writeNode(const Config *config, rapidjson::Value &rootJSONData, rapidjson::Document::AllocatorType &allocator);

std::string decodeValue(const rapidjson::Value &value);
rapidjson::Value encodeValue(const std::string &stringValue, rapidjson::Document::AllocatorType &allocator);

}

}

}
#endif

#endif
