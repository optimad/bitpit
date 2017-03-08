/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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
#ifndef __BITPIT_CONFIGURATION_XML_HPP__
#define __BITPIT_CONFIGURATION_XML_HPP__

#include <libxml/tree.h>
#include <libxml/xmlwriter.h>

#include "configuration_config.hpp"

namespace bitpit {

namespace config {

namespace XML {

extern const std::string DEFAULT_ENCODING;

void readNode(xmlNodePtr root, Config *config);
void writeNode(xmlTextWriterPtr writer, const Config *config,
               const std::string &encoding = DEFAULT_ENCODING);

xmlChar * encodeString(const std::string &in, const std::string &encoding);

}

}

}

#endif
