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
#ifndef __BITPIT_CONFIGURATION_HPP__
#define __BITPIT_CONFIGURATION_HPP__

#include <memory>
#include <unordered_map>

#include <configuration_config.hpp>

namespace bitpit {

// Configuration file parser
class ConfigParser : public Config
{

public:
    ConfigParser(const std::string &root);
    ConfigParser(const std::string &root, const int &version);

    void reset(const std::string &root);
    void reset(const std::string &root, int version);

    void read(const std::string &filename, bool append = true);
    void write(const std::string &filename) const;

private:
    static const int VERSION_UNDEFINED;

    std::string m_root;

    bool m_checkVersion;
    int m_version;

};

// Global configuration file parser
class GlobalConfigParser : public ConfigParser
{

public:
    static GlobalConfigParser & parser();

private:
    static const std::string DEFAULT_ROOT_NAME;
    static const int DEFAULT_VERSION;

    static std::unique_ptr<GlobalConfigParser> m_parser;

    GlobalConfigParser();
    GlobalConfigParser(const std::string &name, int version);

    GlobalConfigParser(GlobalConfigParser const&) = delete;
    GlobalConfigParser& operator=(GlobalConfigParser const&) = delete;

};


/*!
    \brief The namespace 'config' contains routines for interacting with the
    global configuration file parser.
*/
namespace config {

    extern GlobalConfigParser & root;

    void reset();
    void reset(const std::string &name);
    void reset(const std::string &name, int version);

    void read(const std::string &filename, bool append = true);
    void write(const std::string &filename);

};

}

#endif
