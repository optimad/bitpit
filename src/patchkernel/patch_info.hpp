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

#ifndef __BITPIT_PATCH_INFO_HPP__
#define __BITPIT_PATCH_INFO_HPP__

#include <unordered_map>
#include <vector>

namespace bitpit {

class PatchKernel;

class PatchInfo {

public:
	virtual void extract(PatchKernel const *patch) = 0;

protected:
	PatchKernel const *m_patch;

};

#if BITPIT_ENABLE_MPI==1
class PatchGlobalInfo : public PatchInfo {

public:
	PatchGlobalInfo(PatchKernel const *patch);

	void extract(PatchKernel const *patch);

	int getCellRankFromLocal(long id);
	int getCellRankFromGlobal(long id);
	long getCellGlobalId(long id);
	const std::unordered_map<long, long> & getCellGlobalMap();

private:
	std::unordered_map<long, long> m_cellLocalToGlobalMap;
	std::vector<long> m_nGlobalInternals;

};
#endif

}

#endif
