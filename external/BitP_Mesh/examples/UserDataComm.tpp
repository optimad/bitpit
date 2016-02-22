/*
 * UserDataComm.tpp
 *
 *  Created on: 18/mar/2014
 *      Author: Marco Cisternino
 */

template<class Data>
UserDataComm<Data>::UserDataComm(Data & data_, Data & ghostData_) : data(data_), ghostData(ghostData_){};

template<class Data>
UserDataComm<Data>::~UserDataComm() {};

template<class Data>
inline size_t UserDataComm<Data>::fixedSize() const {
	return 0;
};

template<class Data>
inline size_t UserDataComm<Data>::size(const uint32_t e) const {
	return sizeof(double);
};

template<class Data>
template<class Buffer>
inline void UserDataComm<Data>::gather(Buffer& buff, const uint32_t e) {
	buff.write(data[e]);
};

template<class Data>
template<class Buffer>
inline void UserDataComm<Data>::scatter(Buffer& buff,	const uint32_t e) {
	buff.read(ghostData[e]);
};

