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

#ifndef __BITP_MESH_INTERFACE_HPP__
#define __BITP_MESH_INTERFACE_HPP__

/*! \file */

#include <array>
#include <memory>

#include "element.hpp"

/*!
	\ingroup patch
	@{
*/

class Interface : public Element {

public:
	Interface();
	Interface(const long &id, ElementInfo::Type type = ElementInfo::UNDEFINED);

	Interface(Interface&& other) = default;
	Interface& operator=(Interface&& other) = default;

	bool is_border() const;

	std::array<std::array<double, 3>, 3> eval_rotation_from_cartesian();
	static std::array<std::array<double, 3>, 3> eval_rotation_from_cartesian(std::array<double, 3> &versor);
	std::array<std::array<double, 3>, 3> eval_rotation_to_cartesian();
	static std::array<std::array<double, 3>, 3> eval_rotation_to_cartesian(std::array<double, 3> &versor);
	static std::array<std::array<double, 3>, 3> eval_rotation_transpose(const std::array<std::array<double, 3>, 3> &R);

	void set_owner(const long &owner, const int &onwerFace);
	void unset_owner();
	long get_owner() const;
	int get_owner_face() const;

	void set_neigh(const long &neigh, const int &onwerFace);
	void unset_neigh();
	long get_neigh() const;
	int get_neigh_face() const;

	std::array<long, 2> get_owner_neigh() const;

protected:

private:
	long m_owner;
	int m_ownerFace;

	long m_neigh;
	int m_neighFace;

	Interface(const Interface &other) = delete;
	Interface& operator = (const Interface &other) = delete;

};

/*!
	@}
*/

#endif
