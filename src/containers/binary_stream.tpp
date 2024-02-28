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

namespace bitpit {

/*!
* Read the specified array from the stream.
*
* \param[in] stream is the input stream
* \param[in] data is the data to be streamed
* \result Returns the same input stream received in input.
*/
template<typename T, std::size_t SIZE>
IBinaryStream & operator>>(IBinaryStream &stream, std::array<T, SIZE> &data)
{
    for (T &value : data) {
        stream >> value;
    }

    return stream;
}

/*!
* Write the specified array into the stream.
*
* \param[in] stream is the output stream
* \param[in] data is the data to be streamed
* \result Returns the same output stream received in input.
*/
template<typename T, std::size_t SIZE>
OBinaryStream & operator<<(OBinaryStream &stream, const std::array<T, SIZE> &data)
{
    for (const T &value : data) {
        stream << value;
    }

    return stream;
}

/*!
* Read the specified vector from the stream.
*
* Input vector will be resized to match the size of the vector stored in the
* stream.
*
* \param[in] stream is the input stream
* \param[in] data is the data to be streamed
* \result Returns the updated input stream.
*/
template<typename T>
IBinaryStream & operator>>(IBinaryStream &stream, std::vector<T> &data)
{
    std::size_t size;
    stream.read(reinterpret_cast<char *>(&size), sizeof(size));
    data.resize(size);

    for (T &value : data) {
        stream >> value;
    }

    return stream;
}

/*!
* Write the specified vector into the stream.
*
* Along with vector data, also the size of the vector is stored into the
* stream.
*
* \param[in] stream is the output stream
* \param[in] data is the vector to be streamed
* \result Returns the updated output stream.
*/
template<typename T>
OBinaryStream & operator<<(OBinaryStream &stream, const std::vector<T> &data)
{
    std::size_t size = data.size();
    stream.write(reinterpret_cast<const char *>(&size), sizeof(size));

    for (const T &value : data) {
        stream << value;
    }

    return stream;
}

/*!
* Read the specified pair from the stream.
*
* \param[in] stream is the input stream
* \param[in] data is the pair to be streamed
* \result Returns the same input stream received in input.
*/
template<typename T1, typename T2>
IBinaryStream & operator>>(IBinaryStream &stream, std::pair<T1, T2> &data)
{
    T1 first;
    T2 second;
    stream >> first;
    stream >> second;
    data = std::make_pair(std::move(first), std::move(second));

    return stream;
}

/*!
* Write the specified pair into the stream.
*
* \param[in] stream is the input stream
* \param[in] data is the data to be streamed
* \result Returns the same input stream received in input.
*/
template<typename T1, typename T2>
OBinaryStream & operator<<(OBinaryStream &stream, const std::pair<T1, T2> &data)
{
    stream << data.first;
    stream << data.second;

    return stream;
}

/*!
* Read the specified map from the stream.
*
* Input map will be clear and its contents will be replaced with the ones stored into
* the stream.
*
* \param[in] stream is the input stream
* \param[in] data is the data to be streamed
* \result Returns the same input stream received in input.
*/
template<typename K, typename T>
IBinaryStream & operator>>(IBinaryStream &stream, std::unordered_map<K, T> &data)
{
    std::size_t size;
    stream.read(reinterpret_cast<char *>(&size), sizeof(size));

    data.clear();
    data.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
        K key;
        T value;
        stream >> key;
        stream >> value;
        data[key] = std::move(value);
    }

    return stream;
}

/*!
* Write the specified map into the stream.
*
* Along with map data, also the size of the map is stored into the
* stream.
*
* \param[in] stream is the output stream
* \param[in] data is the data to be streamed
* \result Returns the same output stream received in input.
*/
template<typename K, typename T>
OBinaryStream & operator<<(OBinaryStream &stream, const std::unordered_map<K, T> &data)
{
    std::size_t size = data.size();
    stream.write(reinterpret_cast<const char *>(&size), sizeof(size));

    for (const auto &entry : data) {
        stream << entry.first;
        stream << entry.second;
    }

    return stream;
}

/*!
* Read the specified map from the stream.
*
* Input map will be clear and its contents will be replaced with the ones stored into
* the stream.
*
* \param[in] stream is the input stream
* \param[in] data is the data to be streamed
* \result Returns the same input stream received in input.
*/
template<typename K, typename T>
IBinaryStream & operator>>(IBinaryStream &stream, std::map<K, T> &data)
{
    std::size_t size;
    stream.read(reinterpret_cast<char *>(&size), sizeof(size));

    data.clear();
    for (std::size_t i = 0; i < size; ++i) {
        K key;
        T value;
        stream >> key;
        stream >> value;
        data[key] = std::move(value);
    }

    return stream;
}

/*!
* Write the specified map into the stream.
*
* Along with map data, also the size of the map is stored into the
* stream.
*
* \param[in] stream is the output stream
* \param[in] data is the data to be streamed
* \result Returns the same output stream received in input.
*/
template<typename K, typename T>
OBinaryStream & operator<<(OBinaryStream &stream, const std::map<K, T> &data)
{
    std::size_t size = data.size();
    stream.write(reinterpret_cast<const char *>(&size), sizeof(size));

    for (const auto &entry : data) {
        stream << entry.first;
        stream << entry.second;
    }

    return stream;
}

/*!
* Read the specified set from the stream.
*
* Input set will be clear and its contents will be replaced with the ones stored into
* the stream.
*
* \param[in] stream is the input stream
* \param[in] data is the data to be streamed
* \result Returns the same input stream received in input.
*/
template<typename T>
IBinaryStream & operator>>(IBinaryStream &stream, std::unordered_set<T> &data)
{
    std::size_t size;
    stream.read(reinterpret_cast<char *>(&size), sizeof(size));

    data.clear();
    data.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
        T value;
        stream >> value;
        data.insert(std::move(value));
    }

    return stream;
}

/*!
* Write the specified set into the stream.
*
* Along with set data, also the size of the map is stored into the
* stream.
*
* \param[in] stream is the output stream
* \param[in] data is the data to be streamed
* \result Returns the same output stream received in input.
*/
template<typename T>
OBinaryStream & operator<<(OBinaryStream &stream, const std::unordered_set<T> &data)
{
    std::size_t size = data.size();
    stream.write(reinterpret_cast<const char *>(&size), sizeof(size));

    for (const auto &value : data) {
        stream << value;
    }

    return stream;
}

/*!
* Read the specified value from the stream.
*
* \param[in] stream is the input stream
* \param[in] data is the value to be streamed
* \result Returns the updated input stream.
*/
template<typename T>
IBinaryStream & operator>>(IBinaryStream &stream, T &data)
{
    stream.read(reinterpret_cast<char *>(&data), sizeof(T));

    return stream;
}

/*!
* Write the specified value into the stream.
*
* \param[in] stream is the output stream
* \param[in] data is the data to be streamed
* \result Returns the updated output stream.
*/
template<typename T>
OBinaryStream & operator<<(OBinaryStream &stream, const T &data)
{
    stream.write(reinterpret_cast<const char *>(&data), sizeof(T));

    return stream;
}

}
