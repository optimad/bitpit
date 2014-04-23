/*
 * User_Data_Comm.tpp
 *
 *  Created on: 18/mar/2014
 *      Author: Marco Cisternino
 */

template<class Data>
User_Data_Comm<Data>::User_Data_Comm(Data & data_, Data & ghostData_) : data(data_), ghostData(ghostData_){};

template<class Data>
User_Data_Comm<Data>::~User_Data_Comm() {};

template<class Data>
inline size_t User_Data_Comm<Data>::fixedSize() const {
	return 0;
};

template<class Data>
inline size_t User_Data_Comm<Data>::size(const uint32_t e) const {
	return sizeof(double);
};

template<class Data>
template<class Buffer>
inline void User_Data_Comm<Data>::gather(Buffer& buff, const uint32_t e) {
	buff.write(data[e]);
};

template<class Data>
template<class Buffer>
inline void User_Data_Comm<Data>::scatter(Buffer& buff,	const uint32_t e) {
	buff.read(ghostData[e]);
};

