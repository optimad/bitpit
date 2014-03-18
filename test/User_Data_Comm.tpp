/*
 * Userdatacomm.cpp
 *
 *  Created on: 18/mar/2014
 *      Author: Marco Cisternino
 */

#include "Userdatacomm.hpp"

template<class Data>
User_data_comm<Data>::User_data_comm(Data & data_) : data(data_){};

template<class Data>
User_data_comm<Data>::~User_data_comm() {};

template<class Data>
inline size_t User_data_comm<Data>::fixedSize() const {
}

template<class Data>
inline size_t User_data_comm<Data>::size(const uint32_t e) const {
}

template<class Data>
template<class Buffer>
inline void User_data_comm<Data>::gather(Buffer& buff, const uint32_t constUnsignedInt) {
}

template<class Data>
template<class Buffer>
inline void User_data_comm<Data>::scatter(Buffer& buff,	const uint32_t constUnsignedInt) {
}

