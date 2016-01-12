template<class Impl>
inline ClassDataLBInterface<Impl>::ClassDataLBInterface() {};

template<class Impl>
inline size_t ClassDataLBInterface<Impl>::size(const uint32_t e) const {
	return getImpl().size(e);
}

template<class Impl>
inline size_t ClassDataLBInterface<Impl>::fixedSize() const {
	return getImpl().fixedSize();
}

template<class Impl>
inline void ClassDataLBInterface<Impl>::move(const uint32_t from, const uint32_t to) {
	return getImpl().move(from,to);
}

template<class Impl>
template<class Buffer>
inline void ClassDataLBInterface<Impl>::gather(Buffer& buff,	const uint32_t e) {
	return getImpl().gather(buff,e);
}

template<class Impl>
template<class Buffer>
inline void ClassDataLBInterface<Impl>::scatter(Buffer& buff, const uint32_t e) {
	return getImpl().scatter(buff,e);
}

template<class Impl>
inline void ClassDataLBInterface<Impl>::assign(uint32_t stride, uint32_t length){
	return getImpl().assign(stride, length);
}

template<class Impl>
inline void ClassDataLBInterface<Impl>::resize(uint32_t newSize){
	return getImpl().resize(newSize);
}

template<class Impl>
inline void ClassDataLBInterface<Impl>::resizeGhost(uint32_t newSize){
	return getImpl().resizeGhost(newSize);
}

template<class Impl>
inline void ClassDataLBInterface<Impl>::shrink(){
	return getImpl().shrink();
}


template<class Impl>
inline Impl& ClassDataLBInterface<Impl>::getImpl() {
	return static_cast<Impl &>(*this);
}

template<class Impl>
inline const Impl& ClassDataLBInterface<Impl>::getImpl() const {
	return static_cast<const Impl &>(*this);
}

