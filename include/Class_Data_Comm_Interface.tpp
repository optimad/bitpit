/*
 * Class_Data_Comm_Interface.tpp
 *
 *  Created on: 18/mar/2014
 *      Author: Marco Cisternino
 */

template <class Impl>
Class_Data_Comm_Interface<Impl>::Class_Data_Comm_Interface(){}

template <class Impl>
size_t Class_Data_Comm_Interface<Impl>::size(const uint32_t e) const {
	return getImpl().size(e);
};

template <class Impl>
size_t Class_Data_Comm_Interface<Impl>::fixedSize() const {
	return getImpl().fixedSize();
};


template<class Impl>
template<class Buffer>
void Class_Data_Comm_Interface<Impl>::gather(Buffer& buff, const uint32_t e) {
	return getImpl().gather(buff,e);
}

template<class Impl>
template<class Buffer>
void Class_Data_Comm_Interface<Impl>::scatter(Buffer& buff,	const uint32_t e) {
	return getImpl().scatter(buff,e);
}


template <class Impl>
Impl& Class_Data_Comm_Interface<Impl>::getImpl() {
	return static_cast<Impl &>(*this);
}

template <class Impl>
const Impl& Class_Data_Comm_Interface<Impl>::getImpl() const{
	return static_cast<const Impl &>(*this);
}
