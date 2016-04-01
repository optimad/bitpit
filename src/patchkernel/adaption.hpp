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

#ifndef __BITPIT_ADAPTION_HPP__
#define __BITPIT_ADAPTION_HPP__

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
		TYPE_RENUMBERING,
		TYPE_PARTITION_SEND,
		TYPE_PARTITION_RECV,
		TYPE_PARTITION_NOTICE
	};

	enum Entity {
		ENTITY_UNKNOWN = -1,
		ENTITY_CELL,
		ENTITY_INTERFACE
	};

	struct Info
	{
		Info()
			: type(TYPE_UNKNOWN), entity(ENTITY_UNKNOWN), rank(-1)
		{
		}

		Type type;
		Entity entity;
		int rank;
		std::vector<unsigned long> previous;
		std::vector<unsigned long> current;
	};
};

class PatchKernel;

class FlatMapping
{

public:
	FlatMapping();
	FlatMapping(PatchKernel *patch);

	virtual ~FlatMapping();

	virtual void update(const std::vector<Adaption::Info> adaptionData) = 0;

	const std::vector<long> & getNumbering() const;
	const std::vector<long> & getMapping() const;

protected:
	PatchKernel *m_patch;
	std::vector<long> m_numbering;
	std::vector<long> m_mapping;

};


class CellFlatMapping : public FlatMapping
{

public:
	CellFlatMapping();
	CellFlatMapping(PatchKernel *patch);

	~CellFlatMapping();

	void update(const std::vector<Adaption::Info> adaptionData);

};

}

#endif
