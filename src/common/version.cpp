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

#include "version.hpp"

#include "version_verbose.hpp"

namespace bitpit {

/*!
 * Get the bitpit version.
 *
 * \param verbose if set to true, verbose version will be returned.
 */
const std::string & getVersion(bool verbose)
{
    const static std::string VERSION_STRING         = BITPIT_VERSION;
    const static std::string VERBOSE_VERSION_STRING = BITPIT_VERBOSE_VERSION;

    if (verbose) {
        return VERBOSE_VERSION_STRING;
    } else {
        return VERSION_STRING;
    }
}

}
