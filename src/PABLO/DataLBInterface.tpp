/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

template<class Impl>
inline DataLBInterface<Impl>::DataLBInterface() {};

/*! Its user specification computes the specific size of data for an element.
 * \param[in] e The element local index.
 * \return the size of the data for the e element
 */
template<class Impl>
inline size_t DataLBInterface<Impl>::size(const uint32_t e) const {
	return getImpl().size(e);
}

/*! Its user specification computes the same data size for every element in the grid.
 * \return The size of the data for every element
 */
template<class Impl>
inline size_t DataLBInterface<Impl>::fixedSize() const {
	return getImpl().fixedSize();
}

/*! Its user specification moves the data from the "from" element to the "to" element
 * \param[in] from Element local index from.
 * \param[in] to Element local index to.
 */
template<class Impl>
inline void DataLBInterface<Impl>::move(const uint32_t from, const uint32_t to) {
	return getImpl().move(from,to);
}

/*! Its user specification writes the e element data to be communicated
 * during the load balance in the buffer.
 *
 * The user has not to care about the buffer but a char buffer is available in PABLO,
 * Class_Comm_Buffer. This class has an important method, Class_Comm_Buffer#write.
 * This method has to be used to allocate any single element datum
 * in the communication buffer, as follow
 * ~~~~~~~~~~~~~~~~~~~{.c}
 * buff.write(userdatum)
 * ~~~~~~~~~~~~~~~~~~~
 * where userdatum can be any MPI compatible POD variable associated to the e element.
 * In case of a vector of double, called userdata, to store data,
 * ~~~~~~~~~~~~~~~~~~~{.c}
 * buff.write(userdata[e])
 * ~~~~~~~~~~~~~~~~~~~
 * \param[in] buff Template communication buffer
 * \param[in] e The element local index
 */
template<class Impl>
template<class Buffer>
inline void DataLBInterface<Impl>::gather(Buffer& buff,	const uint32_t e) {
	return getImpl().gather(buff,e);
}

/*! Its user specification reads the e element data from the communication buffer
 * and store them in the user data container.
 *
 * The user has not to care about the buffer but a char buffer is available in PABLO,
 * Class_Comm_Buffer. This class has an important method, Class_Comm_Buffer#read().
 * This method has to be used to read any single element datum from
 * the communication buffer, as follow
 * ~~~~~~~~~~~~~~~~~~~{.c}
 * buff.read(userdatum)
 * ~~~~~~~~~~~~~~~~~~~
 * where userdatum can be any MPI compatible POD variable associated to the e element.
 * In case of a vector of double, called userdata, to store data,
 * ~~~~~~~~~~~~~~~~~~~{.c}
 * buff.read(userdata[e])
 * ~~~~~~~~~~~~~~~~~~~
 * \param[in] buff Template communication buffer
 * \param[in] e The element local index
 */
template<class Impl>
template<class Buffer>
inline void DataLBInterface<Impl>::scatter(Buffer& buff, const uint32_t e) {
	return getImpl().scatter(buff,e);
}

/*!  Its user specification extracts contiguous element from a data container.
 * 	This method is used during the very first load balance where the grid from being serial becomes parallel.
 * 	From the initial serial container, which is the same on every process, this method takes the data relatives to the process sub-domain, storing them in the same resized parallel local container.
 * \param[in] stride The location where the process has to start reading from
 * \param[in] length The length of the data container the process has to read
 */
template<class Impl>
inline void DataLBInterface<Impl>::assign(uint32_t stride, uint32_t length){
	return getImpl().assign(stride, length);
}

/*!  Its user specification resizes the element data container to the newSize value.
 * \param[in] newSize The new size of the user data container
 */
template<class Impl>
inline void DataLBInterface<Impl>::resize(uint32_t newSize){
	return getImpl().resize(newSize);
}

/*!  Its user specification resizes the ghost element data container to the newSize value.
 * \param[in] newSize The new size of the user ghost data container
 */
template<class Impl>
inline void DataLBInterface<Impl>::resizeGhost(uint32_t newSize){
	return getImpl().resizeGhost(newSize);
}

/*!  Its user specification reduces the capacity of the container to its size for those containers which has different values for capacity and size.
 * 	If the user container has no capacity concept, nothing has to be done.
 */
template<class Impl>
inline void DataLBInterface<Impl>::shrink(){
	return getImpl().shrink();
}

template<class Impl>
inline Impl& DataLBInterface<Impl>::getImpl() {
	return static_cast<Impl &>(*this);
}

template<class Impl>
inline const Impl& DataLBInterface<Impl>::getImpl() const {
	return static_cast<const Impl &>(*this);
}

}
