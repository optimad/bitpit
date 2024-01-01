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

#include "configuration_common.hpp"

namespace bitpit {

namespace config {

/*!
    Check if the specified file format has a root section.

    \param format is the format of the file
    \result Return true if the specified file format has a root section, false otherwise.
*/
bool hasRootSection(FileFormat format)
{
    if (format == FILE_FORMAT_XML) {
        return true;
    }

    return false;
}

/*!
    Check if the specified file format supports arrays.

    \param format is the format of the file
    \result Return true if the specified file format support arrays, false otherwise.
*/
bool hasArraySupport(FileFormat format)
{
    if (format == FILE_FORMAT_JSON) {
        return true;
    }

    return false;
}

}

}
