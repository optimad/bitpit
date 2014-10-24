template<class Impl>
inline Class_Data_LB_Interface<Impl>::Class_Data_LB_Interface() {};

template<class Impl>
inline size_t Class_Data_LB_Interface<Impl>::size(const uint32_t e) const {
	return getImpl().size(e);
}

template<class Impl>
inline size_t Class_Data_LB_Interface<Impl>::fixedSize() const {
	return getImpl().fixedSize();
}

template<class Impl>
inline void Class_Data_LB_Interface<Impl>::move(const uint32_t from, const uint32_t to) {
	return getImpl().move(from,to);
}

template<class Impl>
template<class Buffer>
inline void Class_Data_LB_Interface<Impl>::gather(Buffer& buff,	const uint32_t e) {
	return getImpl().gather(buff,e);
}

template<class Impl>
template<class Buffer>
inline void Class_Data_LB_Interface<Impl>::scatter(Buffer& buff, const uint32_t e) {
	return getImpl().scatter(buff,e);
}

template<class Impl>
inline Impl& Class_Data_LB_Interface<Impl>::getImpl() {
	return static_cast<Impl &>(*this);
}

template<class Impl>
inline const Impl& Class_Data_LB_Interface<Impl>::getImpl() const {
	return static_cast<const Impl &>(*this);
}

