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
#ifndef __BITPIT_CONFIGURATION_TPP__
#define __BITPIT_CONFIGURATION_TPP__

#include <sstream>
#include <stdexcept>

namespace bitpit {

/*!
    Gets the specified option.

    If the option does not exists an exception is thrown.

    \param key is the name of the option
    \result The specified option.
*/
template<typename T>
T Config::get(const std::string &key) const
{
    T value;
    if (std::istringstream(get(key)) >> value) {
        return value;
    } else {
        throw std::runtime_error("Unable to convert the option \"" + key + "\"");
    }
}

/*!
    Gets the specified option.

    If the option does not exists the fallback value is returned.

    \param key is the name of the option
    \param fallback is the value that will be returned if the specified
    options does not exist
    \result The specified option or the fallback value if the specified
    options does not exist.
*/
template<typename T>
T Config::get(const std::string &key, const T &fallback) const
{
    if (hasOption(key)) {
        return get<T>(key);
    } else {
        return fallback;
    }
}

/*!
    Set the given option to the specified value

    \param key is the name of the option
    \param value is the value of the option
*/
template<typename T>
void Config::set(const std::string &key, const T &value)
{
    std::ostringstream valueStream;
    if (valueStream << value) {
        set(key, valueStream.str());
    } else {
        throw std::runtime_error("Unable to convert the option \"" + key + "\"");
    }
}

}

#endif
