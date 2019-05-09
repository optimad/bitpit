/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
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

/*! Its user specification writes the e element data to be communicated
 * during the load balance in the buffer.
 *
 * The user has to use the buffer buff as an output binary stream.
 * A stream operator (<<) is provided to allocate any single element datum
 * in the communication buffer, as follow
 * ~~~~~~~~~~~~~~~~~~~{.c}
 * buff << userdatum
 * ~~~~~~~~~~~~~~~~~~~
 * where userdatum can be any MPI compatible POD variable associated to the e element.
 * In case of a vector of double, called userdata, to store data,
 * ~~~~~~~~~~~~~~~~~~~{.c}
 * buff << userdata[e]
 * ~~~~~~~~~~~~~~~~~~~
 * \param[in] buff Output communication buffer
 * \param[in] e The element local index
 */
template<class Impl>
template<class Buffer>
void DataCommInterface<Impl>::gather(Buffer& buff, const uint32_t e) {
	return getImpl().gather(buff,e);
}

/*! Its user specification reads the e element data from the communication buffer
 * and store them in the user data container.
 *
 * The user has to use the buffer buff as an input binary stream.
 * A stream operator (>>) is provided to read any single element datum from
 * the communication buffer, as follow
 * ~~~~~~~~~~~~~~~~~~~{.c}
 * buff >> userdatum
 * ~~~~~~~~~~~~~~~~~~~
 * where userdatum can be any MPI compatible POD variable associated to the e element.
 * In case of a vector of double, called userdata, to store data,
 * ~~~~~~~~~~~~~~~~~~~{.c}
 * buff >> userdata[e]
 * ~~~~~~~~~~~~~~~~~~~
 * \param[in] buff Input communication buffer
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

}
