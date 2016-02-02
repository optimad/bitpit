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
 * Userdatacomm.hpp
 *
 *  Created on: 18/mar/2014
 *      Author: Marco Cisternino
 */

#ifndef __BITPIT_USERDATACOMM_HPP__
#define __BITPIT_USERDATACOMM_HPP__

#include "DataCommInterface.hpp"

/*!  \cond  EXAMPLE_CLASSES */
template <class D>
class UserDataComm : public bitpit::DataCommInterface< UserDataComm<D> > {
public:

	typedef D Data;

	Data & data;
	Data & ghostData;

	size_t fixedSize() const;
	size_t size(const uint32_t e) const;

	template<class Buffer>
	void gather(Buffer & buff, const uint32_t e);

	template<class Buffer>
	void scatter(Buffer & buff, const uint32_t e);

	UserDataComm(Data & data_, Data & ghostData_);
	~UserDataComm();
};
/*!  \endcond  */

#include "PABLO_userDataComm.tpp"

#endif /* __BITPIT_USERDATACOMM_HPP__ */
