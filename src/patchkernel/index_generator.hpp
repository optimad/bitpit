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

#ifndef __BITPIT_INDEX_GENERATOR_HPP__
#define __BITPIT_INDEX_GENERATOR_HPP__

#include <deque>
#include <iostream>

namespace bitpit {

class IndexGenerator {

public:
	static const long NULL_ID;

	IndexGenerator();

	long generateId();
	long getLatestId();
	long getHighestId();
	bool isIdAssigned(long id);
	void setAssignedId(long id);
	void trashId(const long &id);
	void reset();

	void dump(std::ostream &stream);
	void restore(std::istream &stream);

private:
	long m_latest;
	long m_highest;
	std::deque<long> m_trash;

	int getBinaryArchiveVersion();

};

}

#endif
