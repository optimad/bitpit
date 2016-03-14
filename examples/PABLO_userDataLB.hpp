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

/*
 * User_data_LB.hpp
 *
 *  Created on: 27/mar/2014
 *      Author: Marco Cisternino
 */

#ifndef __BITPIT_USER_DATA_LB_HPP__
#define __BITPIT_USER_DATA_LB_HPP__

#include "DataCommInterface.hpp"

/*!  \cond  EXAMPLE_CLASSES */
/*! @ingroup PABLO 
 * @{ 
 *
 *  \brief User implementation of CRTP interface for load balance data communications in Message Passing Paradigm framework
 *
 *  UserDataLB is an example of user specification of the data communication interface based on the Curiously Recurrent Template Pattern.
 *  The user has to implement his interface class(es) in order to define how his data have to be written and read in the communication buffer.
 *  These classes have to be derived from the template base class bitpit::DataLBInterface using as template argument value the derived class.
 *
 *  The choice of the members of the class is completely up to the user and they have to be useful to access the internal data container.
 *  In any case, the user must at least implement all the methods reported in this example.
 *  This specific example is intended for communication during the load balance of POD data stored in a container provided with random access operator.
 *  See the examples to better understand how to use it.
 */
template <class D>
class UserDataLB : public bitpit::DataLBInterface<UserDataLB<D> >{
public:

	typedef D Data;

	Data& data;
	Data& ghostdata;

	size_t fixedSize() const;
	size_t size(const uint32_t e) const;
	void move(const uint32_t from, const uint32_t to);

	template<class Buffer>
	void gather(Buffer & buff, const uint32_t e);

	template<class Buffer>
	void scatter(Buffer & buff, const uint32_t e);

	void assign(uint32_t stride, uint32_t length);
	void resize(uint32_t newSize);
	void resizeGhost(uint32_t newSize);
	void shrink();

	UserDataLB(Data& data_, Data& ghostdata_);
	~UserDataLB();
};
/*!  \endcond  */

/*!This method is automatically called by the manager and it is intended to define the constant size of the data to be communicated per grid element.
 * If all the element of the grid communicate the same size of data, this method must return a value different from zero and equal to the number of byte to be communicated for every element.
 * Otherwise, it must return zero and the different data size for each element must be specified by the size method
 * \param[out] the constant size in bytes of the data to be communicated per element.
 */
template<class D>
inline size_t UserDataLB<D>::fixedSize() const {
	return 0;
}

/*!This method is automatically called by the manager and it is intended to define the variable size of the data to be communicated of every grid element.
 * In order to make the manager use this method, the fixedsize method has to return zero.
 * Implementing this method, the user can pass to the manager the specific data size to be communicated for the element e.
 * \param[in] index of the internal element to be communicated.
 * \param[out] the size in bytes of the data tobe communicated for the element e
 */
template<class D>
inline size_t UserDataLB<D>::size(const uint32_t e) const {
	return sizeof(double);
}

template<class D>
inline void UserDataLB<D>::move(const uint32_t from, const uint32_t to) {
	data[to] = data[from];
}

/*!This method is automatically called by the manager and it is intended to write user data to the char communication buffer.
 * The user has to specify in its implementation the way data can be read from the char buffer.
 * The manager provide the user with a buffer and its read method in order to simply read POD data from the buffer.
 * The user has to define the way his data can be read from the buffer by decomposing them in POD data and by using the buffer read method to take them from the buffer and to store them in the ghost data container.
 * In this example we suppose that data is a container of POD data having the random access operator.
 * \param[in] e index of the ghost element to be communicated.
 * \param[in] buff is the char buffer used to communicate user data. The user has not to take care of the buffer, but its method write and read. These methods are intended to write/read POD data to/fromthe buffer.
 */
template<class D>
template<class Buffer>
inline void UserDataLB<D>::gather(Buffer& buff, const uint32_t e) {
	buff.write(data[e]);
}

/*!This method is automatically called by the manager and it is intended to read user data from the char communication buffer and store them in the data container.
 * The user has to specify in its implementation the way data can be read from the char buffer.
 * The manager provide the user with a buffer and its read method in order to simply read POD data from the buffer.
 * The user has to define the way his data can be read from the buffer by decomposing them in POD data and by using the buffer read method to take them from the buffer and to store them in the ghost data container.
 * In this example we suppose that data is a container of POD data having the random access operator.
 * \param[in] e index of the element to be communicated.
 * \param[in] buff is the char buffer used to communicate user data. The user has not to take care of the buffer, but its method write and read. These methods are intended to write/read POD data to/fromthe buffer.
 */
template<class D>
template<class Buffer>
inline void UserDataLB<D>::scatter(Buffer& buff, const uint32_t e) {
	buff.read(data[e]);
}

/*!This method is automatically called by the manager and it intended to be used during the static load balance.
 * At the beginning of any application of the this manager, the entire grid is on every MPI process. At the first load balance the grid is partitioned and distributed to the processes. Data are distributed using this method, by assigning the right length to the the process starting from the right stride.
 * The user has to tell the manager how to start reading the length of data from the stride position and how to put only these data in the same container, deleting the rest. It is a restriction operation of the container to a contiguous part of it.
 * In this examples data is supposed to be a container with the assign operator, as std::vector.
 * \param[in] stride the initial position of the data to be 
 */
template<class D>
inline void UserDataLB<D>::assign(uint32_t stride, uint32_t length) {
	Data dataCopy = data;
	typename Data::iterator first = dataCopy.begin() + stride;
	typename Data::iterator last = first + length;
	data.assign(first,last);
#if defined(__INTEL_COMPILER)
#else
	data.shrink_to_fit();
#endif
	first = dataCopy.end();
	last = dataCopy.end();
};


/*!This method is automatically called by the manager and it intended to be used during the static and the dynamic load balance.
 * During the load balance process the user data container has to change its size to adapt itself to the new partition of the domain.
 * The user has to tell the manager how to change the size of his data containers.
 * In this examples data is supposed to be a container with the resize operator, as std::vector.
 * \param[in] newSize is the new size of the container
 */
template<class D>
inline void UserDataLB<D>::resize(uint32_t newSize) {
	data.resize(newSize);
}

/*!This method is automatically called by the manager and it intended to be used during the static and the dynamic load balance.
 * During the load balance process the ghost user data container has to change its size to adapt itself to the new partition of the domain.
 * The user has to tell the manager how to change the size of his ghost data containers.
 * In this examples data is supposed to be a container with the resize operator, as std::vector.
 * \param[in] newSize is the new size of the container
 */
template<class D>
inline void UserDataLB<D>::resizeGhost(uint32_t newSize) {
	ghostdata.resize(newSize);
}

/*!This method is automatically called by the manager and it intended to be used during the static and the dynamic load balance.
 * During the load balance process the user data container changes its size to adapt itself to the new partition of the domain. If the container can change its size and its capacity separately, this method is intended to get the container capacity equal to its size
 * The user has to tell the manager how to change the capacity of his data containers to its size
 * In this examples data is supposed to be a container with the shrink_to_fit operator, as std::vector.
 */
template<class D>
inline void UserDataLB<D>::shrink() {
#if defined(__INTEL_COMPILER)
#else
	data.shrink_to_fit();
#endif
}

/*!Custom constructor of the user data communication class for load balance
 * \param[in] internal data structure.
 * \param[in] ghost data structure.
 */
template<class D>
inline UserDataLB<D>::UserDataLB(Data& data_, Data& ghostdata_) : data(data_), ghostdata(ghostdata_){}

/*!Destructor of the user data communication class for load balance.
 *
 */
template<class D>
inline UserDataLB<D>::~UserDataLB() {}


#endif /* __BITPIT_USER_DATA_LB_HPP__ */
