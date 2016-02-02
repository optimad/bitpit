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
 * UserDataLB.tpp
 *
 *  Created on: 27/mar/2014
 *      Author: Marco Cisternino
 */

template<class D>
inline size_t UserDataLB<D>::fixedSize() const {
	return 0;
}

template<class D>
inline size_t UserDataLB<D>::size(const uint32_t e) const {
	return sizeof(double);
}

template<class D>
inline void UserDataLB<D>::move(const uint32_t from, const uint32_t to) {
	data[to] = data[from];
}

template<class D>
template<class Buffer>
inline void UserDataLB<D>::gather(Buffer& buff, const uint32_t e) {
	buff.write(data[e]);
}

template<class D>
template<class Buffer>
inline void UserDataLB<D>::scatter(Buffer& buff, const uint32_t e) {
	buff.read(data[e]);
}

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

template<class D>
inline void UserDataLB<D>::resize(uint32_t newSize) {
	data.resize(newSize);
}

template<class D>
inline void UserDataLB<D>::resizeGhost(uint32_t newSize) {
	ghostdata.resize(newSize);
}

template<class D>
inline void UserDataLB<D>::shrink() {
#if defined(__INTEL_COMPILER)
#else
	data.shrink_to_fit();
#endif
}

template<class D>
inline UserDataLB<D>::UserDataLB(Data& data_, Data& ghostdata_) : data(data_), ghostdata(ghostdata_){}

template<class D>
inline UserDataLB<D>::~UserDataLB() {}
