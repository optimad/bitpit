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

#ifndef __BITPIT_BINARY_STREAM_HPP__
#define __BITPIT_BINARY_STREAM_HPP__

#include <array>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace bitpit {

// Stream operators
class IBinaryStream;
class OBinaryStream;

template<typename T, std::size_t d>
IBinaryStream &operator>>(IBinaryStream &stream, std::array<T, d> &data);
template<typename T, std::size_t d>
OBinaryStream &operator<<(OBinaryStream &stream, const std::array<T, d> &data);

template<typename T>
IBinaryStream& operator>>(IBinaryStream &stream, std::vector<T> &data);
template<typename T>
OBinaryStream& operator<<(OBinaryStream &stream, const std::vector<T> &data);

template<typename K, typename T>
IBinaryStream &operator>>(IBinaryStream &stream, std::pair<K, T> &data);
template<typename K, typename T>
OBinaryStream &operator<<(OBinaryStream &stream, const std::pair<K, T> &data);

template<typename K, typename T>
IBinaryStream &operator>>(IBinaryStream &stream, std::map<K, T> &data);
template<typename K, typename T>
OBinaryStream &operator<<(OBinaryStream &stream, const std::map<K, T> &data);

template<typename K, typename T>
IBinaryStream &operator>>(IBinaryStream &stream, std::unordered_map<K, T> &data);
template<typename K, typename T>
OBinaryStream &operator<<(OBinaryStream &stream, const std::unordered_map<K, T> &data);

template<typename T>
IBinaryStream &operator>>(IBinaryStream &stream, std::unordered_set<T> &data);
template<typename T>
OBinaryStream &operator<<(OBinaryStream &stream, const std::unordered_set<T> &data);

template<typename T>
IBinaryStream & operator>>(IBinaryStream &stream, T &data);
template<typename T>
OBinaryStream & operator<<(OBinaryStream &stream, const T &data);

template<>
IBinaryStream & operator>>(IBinaryStream &stream, std::string &data);
template<>
OBinaryStream & operator<<(OBinaryStream &stream, const std::string &data);

// Binary stream
class BinaryStream {

public:
    virtual ~BinaryStream() = default;

    bool eof() const;

    std::streampos tellg() const;
    bool seekg(std::size_t pos);
    bool seekg(std::streamoff offset, std::ios_base::seekdir way);

    char * data();
    const char * data() const;

    virtual void setSize(std::size_t size);
    std::size_t getSize() const;

    std::size_t getCapacity() const;

    int getChunkSize();
    int getChunkCount();

protected:
    std::vector<char> m_buffer;
    std::size_t m_size;
    std::size_t m_pos;

    BinaryStream();
    BinaryStream(std::size_t capacity);
    BinaryStream(const char *buffer, std::size_t capacity);
    BinaryStream(const std::vector<char> &buffer);

    void open(const char *buffer, std::size_t capacity);
    void open(std::size_t capacity);

    void setCapacity(std::size_t capacity);

private:
    int m_chunkSize;

};

class IBinaryStream : public BinaryStream {

template<typename T>
friend IBinaryStream & (operator>>)(IBinaryStream &stream, T &data);

public:
    IBinaryStream(void);
    IBinaryStream(std::size_t size);
    IBinaryStream(const char *buffer, std::size_t size);
    IBinaryStream(const std::vector<char> &buffer);

    void open(const char *buffer, std::size_t size);
    void open(std::size_t size);

    void read(char *data, std::size_t size);

};

class OBinaryStream : public BinaryStream {

template<typename T>
friend OBinaryStream & (operator<<)(OBinaryStream &stream, const T &data);

public:
    OBinaryStream();
    OBinaryStream(std::size_t size);

    void open(std::size_t size);

    void setSize(std::size_t size) override;
    void squeeze();

    void write(const char *data, std::size_t size);

private:
    bool m_expandable;

};

}

// Template implementation
#include "binary_stream.tpp"

#endif
