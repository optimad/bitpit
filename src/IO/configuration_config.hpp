/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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
#ifndef __BITPIT_CONFIGURATION_CONFIG_HPP__
#define __BITPIT_CONFIGURATION_CONFIG_HPP__

#include <memory>
#include <map>
#include <vector>

namespace bitpit {

// Configuration storage
class Config
{

public:
    typedef Config Section;
    typedef std::vector<Section *> MultiSection;
    typedef std::vector<const Section *> ConstMultiSection;
    typedef std::map<std::string, std::string> Options;
    typedef std::multimap<std::string, std::unique_ptr<Config>> Sections;

    Config(bool multiSections);
    Config(const Config &other);
    Config(Config &&other) = default;

    Config & operator=(Config other);
    Config & operator=(Config &&other) = default;

    virtual ~Config();

    void swap(Config &other);

    bool isMultiSectionsEnabled() const;

    int getOptionCount() const;
    Options & getOptions();
    const Options & getOptions() const;
    bool hasOption(const std::string &key) const;
    std::string get(const std::string &key) const;
    std::string get(const std::string &key, const std::string &fallback) const;
    void set(const std::string &key, const std::string &value);
    bool removeOption(const std::string &key);

    int getSectionCount() const;
    int getSectionCount(const std::string &key) const;
    Sections & getSections();
    const Sections & getSections() const;
    MultiSection getSections(const std::string &key);
    ConstMultiSection getSections(const std::string &key) const;
    bool hasSection(const std::string &key) const;
    Section & getSection(const std::string &key);
    const Section & getSection(const std::string &key) const;
    Section & addSection(const std::string &key);
    bool removeSection(const std::string &key);

    void clear();

    void dump(std::ostream &out, int indentLevel = 0) const;

    Section & operator[](const std::string &key);
    const Section & operator[](const std::string &key) const;

    template<typename T>
    T get(const std::string &key) const;

    template<typename T>
    T get(const std::string &key, const T &fallback) const;

    template<typename T>
    void set(const std::string &key, const T &value);

private:
    bool m_multiSections;
    std::unique_ptr<Options> m_options;
    std::unique_ptr<Sections> m_sections;

};

}

#include "configuration_config.tpp"

#endif
