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

#ifndef __BITPIT_INTERFACE_HPP__
#define __BITPIT_INTERFACE_HPP__

#include <array>
#include <memory>

#include "bitpit_containers.hpp"

#include "element.hpp"

namespace bitpit {

class Interface : public Element {

public:
	Interface();
	Interface(const long &id, ElementInfo::Type type = ElementInfo::UNDEFINED);

	Interface(Interface&& other) = default;
	Interface& operator=(Interface&& other) = default;

	bool isBorder() const;

	std::array<std::array<double, 3>, 3> evalRotationFromCartesian();
	static std::array<std::array<double, 3>, 3> evalRotationFromCartesian(std::array<double, 3> &versor);
	std::array<std::array<double, 3>, 3> evalRotationToCartesian();
	static std::array<std::array<double, 3>, 3> evalRotationToCartesian(std::array<double, 3> &versor);
	static std::array<std::array<double, 3>, 3> evalRotationTranspose(const std::array<std::array<double, 3>, 3> &R);

	void setOwner(const long &owner, const int &onwerFace);
	void unsetOwner();
	long getOwner() const;
	int getOwnerFace() const;

	void setNeigh(const long &neigh, const int &onwerFace);
	void unsetNeigh();
	long getNeigh() const;
	int getNeighFace() const;

	std::array<long, 2> getOwnerNeigh() const;

protected:

private:
	long m_owner;
	int m_ownerFace;

	long m_neigh;
	int m_neighFace;

	Interface(const Interface &other) = delete;
	Interface& operator = (const Interface &other) = delete;

};

extern template class bitpit::PiercedVector<Interface>;
extern template class bitpit::PositionalPiercedVector<Interface>;

}

#endif
