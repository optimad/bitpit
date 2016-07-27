/*---------------------------------------------------------------------------*\
*
*  bitpit
*
*  Copyright (C) 2015-2016 OPTIMAD engineering Srl
*
*  -------------------------------------------------------------------------
*  License
*  This file is part of bitbit.
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

#include <stdexcept>

#include "configuration_config.hpp"

namespace bitpit {

/*!
    @ingroup Configuration
    @{
*/

/*!
    \class Config

    \brief Configuration storage

    This class implements a configuration storage.
*/

/*!
    Create a new configuration.
*/
Config::Config()
{
    m_options  = std::unique_ptr<Options>(new Options());
    m_sections = std::unique_ptr<Sections>(new Sections());
}

/*!
    Count the number of options.

    \result The number of options stored.
*/
int Config::getOptionCount() const
{
    return m_options->size();
}

/*!
    Gets a reference to the stored options.

    \result A reference to the stored options.
*/
Config::Options & Config::getOptions()
{
    return const_cast<Options &>(static_cast<const Config &>(*this).getOptions());
}

/*!
    Gets a constant reference to the stored options.

    \result A constant reference to the stored options.
*/
const Config::Options & Config::getOptions() const
{
    return *m_options;
}

/*!
    Checks if the specified option exists.

    \param key is the name of the option
    \result True is the option exists, false otherwise.
*/
bool Config::hasOption(const std::string &key) const
{
    return (m_options->count(key) > 0);
}

/*!
    Gets the specified option.

    If the option does not exists an exception is thrown.

    \param key is the name of the option
    \result The specified option.
*/
std::string Config::get(const std::string &key) const
{
    return m_options->at(key);
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
std::string Config::get(const std::string &key, const std::string &fallback) const
{
    if (hasOption(key)) {
        return get(key);
    } else {
        return fallback;
    }
}

/*!
    Set the given option to the specified value

    \param key is the name of the option
    \param value is the value of the option
*/
void Config::set(const std::string &key, const std::string &value)
{
    (*m_options)[key] = value;
}

/*!
    Remove the specified option.

    \param key is the name of the option
    \result Returns true if the option existed, otherwise it returns false.
*/
bool Config::removeOption(const std::string &key)
{
    return (m_options->erase(key) != 0);
}

/*!
    Count the number of sections.

    \result The number of sections stored.
*/
int Config::getSectionCount() const
{
    return m_sections->size();
}

/*!
    Gets a reference to the stored sections.

    \result A reference to the stored sections.
*/
Config::Sections & Config::getSections()
{
    return const_cast<Sections &>(static_cast<const Config &>(*this).getSections());
}

/*!
    Gets a constant reference to the stored options.

    \result A constant reference to the stored options.
*/
const Config::Sections & Config::getSections() const
{
    return *m_sections;
}

/*!
    Checks if the specified section exists.

    \param key is the name of the section
    \result True is the section exists, false otherwise.
*/
bool Config::hasSection(const std::string &key) const
{
    return (m_sections->count(key) > 0);
}

/*!
    Gets a reference to the specified section.

    If the section does not exists an exception is thrown.

    \param key is the name of the section
    \result A reference to the specified section.
*/
Config::Section & Config::getSection(const std::string &key)
{
    return const_cast<Section &>(static_cast<const Config &>(*this).getSection(key));
}

/*!
    Gets a constant reference to the specified section.

    If the section does not exists an exception is thrown.

    \param key is the name of the section
    \result A constant reference to the specified section.
*/
const Config::Section & Config::getSection(const std::string &key) const
{
    return *(m_sections->at(key));
}

/*!
    Add a section with the specified name.

    If a section with the given name already exists, an exception is raised.

    \param key is the name of the section
    \return A reference to the newly added section.
*/
Config::Section & Config::addSection(const std::string &key)
{
    if (hasSection(key)) {
        throw std::runtime_error("A section named \"" + key + "\" already esists");
    }

    (*m_sections)[key] = std::unique_ptr<Section>(new Section());

    return *((*m_sections)[key]);
}

/*!
    Remove the specified section.

    \param key is the name of the section
    \result Returns true if the section existed, otherwise it returns false.
*/
bool Config::removeSection(const std::string &key)
{
    return (m_sections->erase(key) != 0);
}

/*!
    Clear the configuration.
*/
void Config::clear()
{
    m_options->clear();
    m_sections->clear();
}

/*!
    Get a constant reference of the specified section.

    If the section does not exists an exception is thrown.

    \param key is the name of the section
    \result A constant reference to the specified section.
*/
const Config::Section & Config::operator[](const std::string &key) const
{
    return getSection(key);
}

/*!
    Get a reference of the specified section.

    If the section does not exists an exception is thrown.

    \param key is the name of the section
    \result A reference to the specified section.
*/
Config::Section & Config::operator[](const std::string &key)
{
    return getSection(key);
}

/*!
    Write the specified configuration to screen.

    \param configuration is the configuration to dump
    \param indentLevel is the indentation level
*/
void Config::dump(std::ostream &out, int indentLevel) const
{
    const int INDENT_SIZE = 2;

    std::string padding(INDENT_SIZE, ' ');
    std::string indent(INDENT_SIZE * indentLevel, ' ');

    out << std::endl;
    out << indent << "Options..." << std::endl;
    if (getOptionCount() > 0) {
        for (auto &entry : getOptions()) {
            out << indent << padding << entry.first << " = " << entry.second << std::endl;
        }
    } else {
        out << indent << padding << "No options." << std::endl;
    }

    ++indentLevel;
    for (auto &entry : getSections()) {
        out << std::endl;
        out << indent << padding << "::: Section " << entry.first << " :::" << std::endl;
        entry.second->dump(out, indentLevel);
    }
}

/*!
    @}
*/

}
