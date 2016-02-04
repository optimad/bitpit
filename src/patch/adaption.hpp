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

#ifndef __PATCHMAN_TYPE_HPP__
#define __PATCHMAN_TYPE_HPP__

/*! \file */

#include <vector>

namespace bitpit {

struct Adaption
{
	enum Type {
		TYPE_UNKNOWN = -1,
		TYPE_CREATION = 0,
		TYPE_DELETION,
		TYPE_REFINEMENT,
		TYPE_COARSENING,
		TYPE_RENUMBERING
	};

	enum Entity {
		ENTITY_UNKNOWN = -1,
		ENTITY_CELL,
		ENTITY_INTERFACE
	};

	struct Info
	{
		Info()
			: type(TYPE_UNKNOWN), entity(ENTITY_UNKNOWN)
		{
		}

		Type type;
		Entity entity;
		std::vector<unsigned long> previous;
		std::vector<unsigned long> current;
	};
};

class Patch;

class FlatMapping
{

public:
	FlatMapping();
	FlatMapping(Patch *patch);

	virtual ~FlatMapping();

	virtual void update(const std::vector<Adaption::Info> adaptionData) = 0;

	const std::vector<long> & get_numbering() const;
	const std::vector<long> & get_mapping() const;

protected:
	Patch *m_patch;
	std::vector<long> m_numbering;
	std::vector<long> m_mapping;

};


class CellFlatMapping : public FlatMapping
{

public:
	CellFlatMapping();
	CellFlatMapping(Patch *patch);

	~CellFlatMapping();

	void update(const std::vector<Adaption::Info> adaptionData);

};

}

#endif
