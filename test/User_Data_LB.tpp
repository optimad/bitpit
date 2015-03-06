/*
 * User_Data_LB.tpp
 *
 *  Created on: 27/mar/2014
 *      Author: Marco Cisternino
 */

template<class D>
inline size_t User_Data_LB<D>::fixedSize() const {
	return 0;
}

template<class D>
inline size_t User_Data_LB<D>::size(const uint32_t e) const {
	return sizeof(double);
}

template<class D>
inline void User_Data_LB<D>::move(const uint32_t from, const uint32_t to) {
	data[to] = data[from];
}

template<class D>
template<class Buffer>
inline void User_Data_LB<D>::gather(Buffer& buff, const uint32_t e) {
	buff.write(data[e]);
}

template<class D>
template<class Buffer>
inline void User_Data_LB<D>::scatter(Buffer& buff, const uint32_t e) {
	buff.read(data[e]);
}

template<class D>
inline void User_Data_LB<D>::assign(uint32_t stride, uint32_t length) {
	Data dataCopy = data;
	typename Data::iterator first = dataCopy.begin() + stride;
	typename Data::iterator last = first + length;
	data.assign(first,last);
	data.shrink_to_fit();
	first = dataCopy.end();
	last = dataCopy.end();
};

template<class D>
inline void User_Data_LB<D>::resize(uint32_t newSize) {
	data.resize(newSize);
}

template<class D>
inline void User_Data_LB<D>::shrink() {
	data.shrink_to_fit();
}

template<class D>
inline User_Data_LB<D>::User_Data_LB(Data& data_) : data(data_){}

template<class D>
inline User_Data_LB<D>::~User_Data_LB() {}
