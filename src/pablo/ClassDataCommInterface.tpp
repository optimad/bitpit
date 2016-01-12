template <class Impl>
ClassDataCommInterface<Impl>::ClassDataCommInterface(){}

template <class Impl>
size_t ClassDataCommInterface<Impl>::size(const uint32_t e) const {
	return getImpl().size(e);
};

template <class Impl>
size_t ClassDataCommInterface<Impl>::fixedSize() const {
	return getImpl().fixedSize();
};


template<class Impl>
template<class Buffer>
void ClassDataCommInterface<Impl>::gather(Buffer& buff, const uint32_t e) {
	return getImpl().gather(buff,e);
}

template<class Impl>
template<class Buffer>
void ClassDataCommInterface<Impl>::scatter(Buffer& buff,	const uint32_t e) {
	return getImpl().scatter(buff,e);
}


template <class Impl>
Impl& ClassDataCommInterface<Impl>::getImpl() {
	return static_cast<Impl &>(*this);
}

template <class Impl>
const Impl& ClassDataCommInterface<Impl>::getImpl() const{
	return static_cast<const Impl &>(*this);
}
