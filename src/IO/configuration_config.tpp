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
    if (std::istringstream(getOption(key).value) >> value) {
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
        getOption(key).value = valueStream.str();
    } else {
        throw std::runtime_error("Unable to convert the option \"" + key + "\"");
    }
}

/*!
    Gets the value of the specified option attribute.

    If the option or the attribute does not exists, an exception is thrown.

    \param key is the name of the option
    \param name is the name of the attribute
    \result The value of the specified attribute.
*/
template<typename T>
T Config::getAttribute(const std::string &key, const std::string &name) const
{
    T value;
    if (std::istringstream(getOption(key).attributes.at(name)) >> value) {
        return value;
    } else {
        throw std::runtime_error("Unable to convert the option \"" + key + "\"");
    }
}

/*!
    Gets the value of the specified option attribute.

    If the option does not exists, an exception will be thrown. However, if
    the attribute do not exists, the fallback walue will be returned

    \param key is the name of the option
    \param name is the name of the attribute
    \param fallback is the value that will be returned if the specified
    attribute does not exist
    \result The value of the specified attribute or the fallback value if
    the options or the attribute does not exist.
*/
template<typename T>
T Config::getAttribute(const std::string &key, const std::string &name, const T &fallback) const
{
    const Option &option = getOption(key);
    if (option.attributes.count(name) > 0) {
        return getAttribute<T>(key);
    }

    return fallback;
}

/*!
    Set the value of the specified option attribute.

    If the option does not exists, an exception will be thrown. However,
    if the attribute does not exists, a new attribute will be added.

    \param key is the name of the option
    \param name is the name of the attribute
    \param value is the value of the attribute
*/
template<typename T>
void Config::setAttribute(const std::string &key, const std::string &name, const T &value)
{
    std::ostringstream valueStream;
    if (valueStream << value) {
        getOption(key).attributes[name] = valueStream.str();
    } else {
        throw std::runtime_error("Unable to convert the option \"" + key + "\"");
    }
}

}

#endif
