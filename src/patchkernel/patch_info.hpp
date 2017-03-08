/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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

#ifndef __BITPIT_PATCH_INFO_HPP__
#define __BITPIT_PATCH_INFO_HPP__

#include <unordered_map>
#include <vector>

namespace bitpit {

class PatchKernel;

class PatchInfo {

public:
	virtual ~PatchInfo();

	void setPatch(PatchKernel const *patch);

	void reset();
	void extract();
	void update();

protected:
	PatchKernel const *m_patch;

	PatchInfo(PatchKernel const *patch);

	void setPatch(PatchKernel const *patch, bool initialize);

	virtual void _init() = 0;
	virtual void _reset() = 0;
	virtual void _extract() = 0;

};

#if BITPIT_ENABLE_MPI==1
class PatchGlobalInfo : public PatchInfo {

public:
	PatchGlobalInfo(PatchKernel const *patch = nullptr);

	long getCellGlobalCount() const;

	int getCellRankFromLocal(long id) const;
	int getCellRankFromGlobal(long id) const;
	long getCellGlobalId(long id) const;

	long getCellGlobalOffset() const;
	long getCellGlobalOffset(int rank) const;
	const std::unordered_map<long, long> & getCellGlobalMap() const;

protected:
	void _init();
	void _reset();
	void _extract();

private:
	std::unordered_map<long, long> m_cellLocalToGlobalMap;
	std::vector<long> m_nGlobalInternals;

};
#endif

}

#endif
