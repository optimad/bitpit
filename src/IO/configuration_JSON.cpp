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

#include "configuration_JSON.hpp"

#if HAS_RAPIDJSON_LIB

#include "rapidjson/encodedstream.h"
#include <rapidjson/error/en.h>
#include "rapidjson/filereadstream.h"
#include <rapidjson/filewritestream.h>
#include <rapidjson/prettywriter.h>

#include <string>
#include <cstdio>
#include <cassert>
#include <set>
#include <iomanip>

namespace bitpit {

namespace config {

namespace JSON {

/*!
    Reading a json dictionary from file and fill a bitpit::Config tree
    accordingly.

    Inside, the method creates a DOM structure from json file and scan it
    using readNode method.

    JSON root document object does not have a key name and a version unlike
    XML.

    Encoding is always of UTF type (8,16,32). Here it is adopted the standard
    UTF-8 encoding. JSON dictionary can be written in other UTF format, the
    conversion is automatically handled by rapidjson parsing.

    \param[in] filename is the path of the JSON file
    \param[in,out] rootConfig is a pointer to Config tree that will be used
    to store the data read from the document
*/
void readConfiguration(const std::string &filename, Config *rootConfig)
{
    if (!rootConfig) {
        throw std::runtime_error("JSON::readDoc Null Config tree structure passed");
    }

    // Open the file
    //
    // To be compatible with Windows platforms, the file is open as a binary
    // file (this works fine also on Posix platforms).
    std::FILE *fp = std::fopen(filename.c_str(), "rb");

    if (!fp) {
        std::string message = "JSON:readDoc impossible to write on file : ";
        message += filename.c_str();
        throw std::runtime_error(message.c_str());
    }

    // Create a buffer for reading stream purposes.
    std::vector<char> readBuffer(65536);
    rapidjson::FileReadStream bis(fp, readBuffer.data(), readBuffer.size());

    // Wrap the readstream into an auto encoding handle. This will let
    // absorb any UTF encoded file into a regular UTF-8 document.
    rapidjson::AutoUTFInputStream<unsigned, rapidjson::FileReadStream> eis(bis);

    // Parse UTF file into UTF-8 in memory
    rapidjson::Document jsonRoot;
    jsonRoot.ParseStream<0, rapidjson::AutoUTF<unsigned> >(eis);

    // Close the file
    std::fclose(fp);

    // Handle parse errors
    if (jsonRoot.HasParseError()) {
        std::string message = "JSON:readDoc error of type : ";
        message += std::string(rapidjson::GetParseError_En(jsonRoot.GetParseError()));
        throw std::runtime_error(message.c_str());
    }

    // Root should be an object, bitpit doens't support root arrays.
    assert(jsonRoot.IsObject() && "JSON:readDoc parsed document root is not an object");

    // Fill the configuration
    readNode("", jsonRoot, rootConfig);
}

/*!
    Read recursively a json object content and fill accordingly the Config tree
    branch.

    The current method visit internal members of JSON object and parse tehm in
    Config::Option or if subnested JSON objects or arrays are present in other
    Config::Section through a recursive call to itself. If arrays are present,
    MultiSection feature of the bitpit::Config tree must be enabled

    The method can be used directly on the root content of the JSON document,
    specifying as key an empty string.

    BEWARE: JSON array support is limited to arrays of JSON objects. Each object
    is converted into a Config::Section and stored with the same key name of the
    array (using MultiSection property of Config). Pod JSON arrays (bools, numbers,
    strings) are not supported since Config tree does not support MultiOptions.
    An attempt to read an array of pod elements will throw an error.


    \param[in] key key name associated of the JSON data.
    \param[in] value JSON data, can be a pod type (bool, number, null,
    str ing) or another JSON object/JSON array of contents. In that case
    automatic recursion will occur.
    \param[in] config pointer to configuration tree to be filled
*/
void readNode(const std::string &key, const rapidjson::Value &value, Config *config)
{
    if (value.IsObject()) {
        // Get the section
        Config::Section *section;
        if (key.empty()) {
            // If the key is empty we are reading the root.
            section = config;
        } else if (!config->isMultiSectionsEnabled() &&config->hasSection(key)) {
            // Get the requestd section
            section = &(config->getSection(key));
        } else {
            // Create a new section.
            //
            // If the name already exists, it will be handled as a multi-section
            // in the tree map.
            section = &(config->addSection(key));
        }

        // Read the contents of the section
        for (auto itr = value.MemberBegin(); itr != value.MemberEnd(); ++itr) {
            readNode(std::string(itr->name.GetString()), itr->value, section);
        }

    } else if (value.IsArray()) {
        // Check if multi-section configurations are supported
        bool multiSectionsEnabled = config->isMultiSectionsEnabled();
        if (!multiSectionsEnabled) {
            throw std::runtime_error("JSON:readNode reading array but Config MultiSection disabled");
        }

        // Read the array contents
        for (auto itr = value.Begin(); itr != value.End(); ++itr) {
            if (itr->IsObject()) {
                readNode(key, *itr, config);
            } else {
                throw std::runtime_error("JSON:readNode only arrays of JSON Objects are supported!");
            }
        }

    } else {
        // Read a single JSON POD type (bool, numbers, strings null).
        config->set(key, decodeValue(value));
    }
}

/*!
    Write a bitpit::Config tree contents to json formatted dictionary.

    NOTE: JSON root document object does not have a key name and a version
    unlike XML.

    Inside the method creates an empty rapidjson Document, set the root object
    and fill it recursively from the Config tree using the writeNode method.
    Then, it writes the Document on file with standard encoding of UTF-8 type.

    Option prettify controls the json dictionary formatting. If false write it
    in a unique compact string, otherwise write it down in a more readable tree
    structure.

    \param[in] filename is the path to a valid json file.
    \param[in] rootConfig pointer to the Config tree to be written on file
    \param[in] prettify if set to true, a prettyfied human-readable JSON will
    be written, otherwise a compatc JSON will be written
*/
void writeConfiguration(const std::string &filename, const Config *rootConfig, bool prettify)
{
    if (!rootConfig) {
        throw std::runtime_error("JSON::writeConfiguration Null Config tree structure passed");
    }

    // DOM Document is GenericDocument<UTF8<>>
    rapidjson::Document jsonRoot;

    // Get the allocator
    rapidjson::Document::AllocatorType &allocator = jsonRoot.GetAllocator();

    // Create a root node and recursively fill it
    jsonRoot.SetObject();
    writeNode(rootConfig, jsonRoot, allocator);

    // Write on file
    //
    // To be compatible with Windows platforms, the file is open as a binary
    // file (this works fine also on Posix platforms).
    std::FILE *fp = std::fopen(filename.c_str(), "wb");
    if (!fp) {
        std::string message = "JSON:writeConfiguration impossible to write on file : ";
        message += filename.c_str();
        throw std::runtime_error(message.c_str());
    }

    // Create a buffer for writing purposes.
    std::vector<char> writeBuffer(65536);
    rapidjson::FileWriteStream bos(fp, writeBuffer.data(), writeBuffer.size());

    // Write the file
    if (prettify) {
         // Use prettyfied UTF-8 JSON format with default indentantion
         rapidjson::PrettyWriter<rapidjson::FileWriteStream> writer(bos);
         jsonRoot.Accept(writer);

    } else {
        // Use single-ultracompact UTF-8 JSON format
        rapidjson::Writer<rapidjson::FileWriteStream> writer(bos);
        jsonRoot.Accept(writer);
    }

    // Close the file.
    std::fclose(fp);
}

/*!
    Write recursively Config tree contents to JSON objects

    The current method visits the config Tree and write Option and nested
    sections as internal members of JSON object. MultiSections will be
    treated as array of JSON objects.

    The method can be used directly on the root content of the JSON document,
    specifying as key an empty string.

    \param[in] config pointer to local Config tree section to be parsed.
    \param[in] rootJSONdata, reference JSON object value to fill
    \param[in] allocator  Allocator of reference rapidjson Document
*/
void writeNode(const Config *config, rapidjson::Value &rootJSONData, rapidjson::Document::AllocatorType &allocator)
{
    // Write the options
    for (const auto &entry : config->getOptions()) {
        const std::string &configKey   = entry.first;
        const std::string &configValue = entry.second;

        rapidjson::Value jsonKey;
        jsonKey.SetString(configKey, allocator);

        rapidjson::Value jsonValue = encodeValue(configValue, allocator);

        rootJSONData.AddMember(jsonKey, jsonValue, allocator);
    }

    // Get list of unique section keys
    std::set<std::string> listOfSectionsKeys;
    for (const auto &entry : config->getSections()) {
        listOfSectionsKeys.insert(entry.first);
    }

    // Write the sections
    for (const std::string &configKey : listOfSectionsKeys) {
        const Config::ConstMultiSection &keySections = config->getSections(configKey);
        if (keySections.empty()) {
                continue;
        }

        std::size_t nEntries = keySections.size();

        rapidjson::Value object(rapidjson::kObjectType);
        rapidjson::Value arrayOfObjects(rapidjson::kArrayType);
        if (nEntries > 1) {
            arrayOfObjects.Reserve(nEntries, allocator);
        }

        // Process section content
        for (std::size_t i = 0; i< nEntries; ++i) {
            writeNode(keySections[i], object, allocator);
            if (nEntries > 1) {
                arrayOfObjects.PushBack(object, allocator);
                rapidjson::Value swapDummy(rapidjson::kObjectType);
                object.Swap(swapDummy);
            }
        }

        rapidjson::Value jsonKey;
        jsonKey.SetString(configKey, allocator);

        if (nEntries > 1) {
            rootJSONData.AddMember(jsonKey, arrayOfObjects, allocator);
        } else {
            rootJSONData.AddMember(jsonKey, object, allocator);
        }
    }
}

/*!
    Convert a JSON value into a string.

    Supported JSON values are:
    - bool
    - number (int/int64, uint/uint64 or float/double)
    - string
    - null

    Objects and array structures are not stringfyied. Unsupported values and
    null will be converted in an empty string.

    \param[in] value is the JSON value
    \result The string corresponding to the specified JSON value.
*/
std::string decodeValue(const rapidjson::Value &value)
{
    if (value.IsBool()) {
        return std::to_string(int(value.GetBool()));
    } else if (value.IsString()) {
        return std::string(value.GetString());
    } else if (value.IsNumber()) {
        // To avoid problems, numbers are always casted to the maximum allowed
        // datatype.
        std::ostringstream stringStream;
        if (value.IsUint64()) {
            // Try uint casting first
            stringStream << value.GetUint64();

        } else if (value.IsInt64()) {
            // Then the normal int casting
            stringStream << value.GetInt64();
        } else {
            // Otherwise try it as double
            stringStream << std::scientific << std::setprecision(8) << value.GetDouble();
        }

        return stringStream.str();
    } else {
        return std::string("");
    }
}

/*!
    Convert a string into a JSON.

    Supported types are:
    - unsigned long int (e.g., 123456789);
    - signed long int (e.g. -123456789);
    - double (+/-123456789.123456789 or in scientific notation);
    - null (emtpy strings).

    Strings that does't match any know data types are converted to JSON strings.

    \param[in] stringValue is the string that will be converted
    \param[in] allocator is the JSON allocator
    \result The JSON value corresponding to the specified string.
*/
rapidjson::Value encodeValue(const std::string &stringValue, rapidjson::Document::AllocatorType &allocator)
{
    rapidjson::Value value;

    // Empty string are null values
    if (stringValue.empty()) {
        return value;
    }

    // Strings that contains non-numeric charactes are strings
    std::size_t nonNumericCharacterPos = stringValue.find_first_not_of("Ee+-.0123456789");
    if (nonNumericCharacterPos != std::string::npos) {
        value.SetString(stringValue, allocator);
        return value;
    }

    // Convert the string into a number
    std::istringstream stringStream(stringValue);
    if (stringValue.find_first_of("Ee.") != std::string::npos) {
        // Convert the string into a double
        double number;
        stringStream >> number;
        value.SetDouble(number);
    } else if (stringValue.find_first_of("-") != std::string::npos) {
        // Convert the string into a long
        int64_t number;
        stringStream >> number;
        value.SetInt64(number);

    } else {
        // Convert the string into an unsigned long
        uint64_t number;
        stringStream >> number;
        value.SetUint64(number);
    }

    return value;
}

}

}

}

#endif
