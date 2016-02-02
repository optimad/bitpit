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

// Communications Class for Ghosts Update

#ifndef __BITPIT_DATA_COMM_INTERFACE_HPP__
#define __BITPIT_DATA_COMM_INTERFACE_HPP__

namespace bitpit {

/*!
 *  \ingroup PABLO
 *  @{
 *
 *	\date			09/sep/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Base class for data communications
 *
 *	This class is the base class used to implement the user interface to data communications.
 *
 *	The Curiously Recurrent Template Pattern is exploited to achieve the interface.
 *	By this way the interface is based on static polymorphism with no extra cost at runtime.
 *
 *	The user has to implement his communication classes by deriving them from this class.
 *	The mechanism implies that the derived class derives from a template
 *	base class and that the template parameter is the derived class itself, as follow
 *	class Derived : public Base<Derived>{...}
 *
 *	The user has to implement all the methods of the base class in his derived class.
 *	These user's methods will really do the job.
 *
 *	Easily speaking, only the user knows his data and through the interface
 *	specialization he states the size of element data, how to write/read them
 *	 in a communication buffer.
 *	Any MPI compatible POD datum can be written and read in the communication buffer.
 */

template <class Impl>
class DataCommInterface {
public:
	size_t size(const uint32_t e) const;
	size_t fixedSize() const;

	template<class Buffer>
	void gather(Buffer & buff,const uint32_t e);

	template<class Buffer>
	void scatter(Buffer & buff,const uint32_t e);

protected:
	DataCommInterface();

private:
	//BartonHackman trick
	Impl& getImpl();
	const Impl& getImpl() const;
};

/*  @{ */

}

#include "DataCommInterface.tpp"

#endif /* __BITPIT_DATA_COMM_INTERFACE_HPP__ */
