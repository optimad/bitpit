template <class Impl>
DataCommInterface<Impl>::DataCommInterface(){}

/*! Its user specification computes the specific size of data for an element.
 * \param[in] e Element local index.
 * \return the size of the data for the e element
 */
template <class Impl>
size_t DataCommInterface<Impl>::size(const uint32_t e) const {
	return getImpl().size(e);
};

/*! Its user specification computes the same data size for every element in the grid.
 * \return the size of the data for every element
 */
template <class Impl>
size_t DataCommInterface<Impl>::fixedSize() const {
	return getImpl().fixedSize();
};

/*! Its user specification writes the e element data to be communicated in the buffer.
 *
 * The user has not to care about the buffer but a char buffer is available in PABLO,
 * Class_Comm_Buffer. This class has an important method, Class_Comm_Buffer#write.
 *
 * This method has to be used to allocate any single element datum in the
 * communication buffer, as follow
 * ~~~~~~~~~~~~~~~~~~~{.c}
 * buff.write(userdatum)
 * ~~~~~~~~~~~~~~~~~~~
 * where userdatum can be any MPI compatible POD variable associated to the e element.
 *
 * In case of a vector of double, called userdata, to store data,
 * ~~~~~~~~~~~~~~~~~~~{.c}
 * buff.write(userdata[e])
 * ~~~~~~~~~~~~~~~~~~~
 * \param[in] buff Template communication buffer
 * \param[in] e The element local index
 */
template<class Impl>
template<class Buffer>
void DataCommInterface<Impl>::gather(Buffer& buff, const uint32_t e) {
	return getImpl().gather(buff,e);
}

/*! Its user specification reads the e element data from the communication buffer
 * and store them in the ghost user data container.
 *
 * The user has not to care about the buffer but a char buffer is available in PABLO,
 * Class_Comm_Buffer. This class has an important method, Class_Comm_Buffer#read.
 *
 * This method has to be used to read any single element datum from the
 * communication buffer, as follow
 * ~~~~~~~~~~~~~~~~~~~{.c}
 * buff.read(userdatum)
 * ~~~~~~~~~~~~~~~~~~~
 * where userdatum can be any MPI compatible POD variable associated to the e element.
 *
 * In case of a vector of double, called userdata, to store data,
 * ~~~~~~~~~~~~~~~~~~~{.c}
 * buff.read(userdata[e])
 * ~~~~~~~~~~~~~~~~~~~
 * \param[in] buff Template communication buffer
 * \param[in] e The element local index
 */
template<class Impl>
template<class Buffer>
void DataCommInterface<Impl>::scatter(Buffer& buff,	const uint32_t e) {
	return getImpl().scatter(buff,e);
}


template <class Impl>
Impl& DataCommInterface<Impl>::getImpl() {
	return static_cast<Impl &>(*this);
}

template <class Impl>
const Impl& DataCommInterface<Impl>::getImpl() const{
	return static_cast<const Impl &>(*this);
}
