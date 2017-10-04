/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

#ifndef __BITPIT_PATCH_MANAGER_HPP__
#define __BITPIT_PATCH_MANAGER_HPP__

#include <memory>
#include <iostream>
#include <unordered_map>

#include <bitpit_IO.hpp>

namespace bitpit {

class PatchManager {

friend class PatchKernel;

public:
	static int const AUTOMATIC_ID;

    static PatchManager & manager();

	PatchKernel * get(int id);

	void dump(std::ostream &stream);
	void restore(std::istream &stream);

	void dumpAll(std::ostream &stream);
	void restoreAll(std::istream &stream);

private:
    static std::unique_ptr<PatchManager> m_manager;

	IndexGenerator<int> m_idGenerator;
	std::vector<PatchKernel *> m_patchOrder;
	std::unordered_map<PatchKernel *, int> m_patchIds;

	PatchManager();

    PatchManager(PatchManager const&) = delete;
    PatchManager& operator=(PatchManager const&) = delete;

	int registerPatch(PatchKernel *patch, int id = AUTOMATIC_ID);
	void unregisterPatch(PatchKernel *patch);

};

/*!
    \ingroup patchkernel

    \brief The namespace 'patch' contains routines for interacting with the
    patch manager.
*/
namespace patch {

    // Generic global functions
    PatchManager & manager();

}

}

#endif
