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

#ifndef __BITPIT_SURFUNSTRUCTURED_HPP__
#define __BITPIT_SURFUNSTRUCTURED_HPP__

#include <array>
#include <vector>

#include "bitpit_IO.hpp"
#include "bitpit_patchkernel.hpp"

namespace bitpit {

class SurfUnstructured : public SurfaceKernel {

public:
	using PatchKernel::locatePoint;

        // Constructors
	SurfUnstructured();
	SurfUnstructured(const int &id, int patch_dim = 2, int space_dim = 3);
	SurfUnstructured(std::istream &stream);

        // Clone
        std::unique_ptr<PatchKernel> clone() const;

        // Setters
	void setExpert(bool expert);

        // Search algorithms
        long locatePoint(const std::array<double, 3> &point);

        // Evaluations
        void extractEdgeNetwork(SurfUnstructured &);

        // I/O routines
        unsigned short importSTL(const std::string &, int PIDOffset = 0, bool PIDSquah = false);
        unsigned short importSTL(const std::string &, const bool &, int PIDOffset = 0, bool PIDSquah = false);
        unsigned short exportSTL(const std::string &, const bool &, bool flag = true);
        unsigned short exportSTL(const std::string &, const bool &, const bool &, bool flag);
        unsigned short importDGF(const std::string &, int PIDOffset = 0, bool PIDSquah = false);
        unsigned short exportDGF(const std::string &);

protected:
	SurfUnstructured(const SurfUnstructured &other) = default;

	int _getDumpVersion() const;
	void _dump(std::ostream &stream);
	void _restore(std::istream &stream);

	static ElementInfo::Type getSTLFacetType(int nFacetVertices);
	static ElementInfo::Type getDGFFacetType(int nFacetVertices);

	unsigned short importSTL(STLObj &STL, int PIDOffset, bool PIDSquash);

	unsigned short exportSTLSingle(const std::string &, const bool &, bool flag = true);
	unsigned short exportSTLMulti(const std::string &, bool flag = true);

private:

};

}

#endif
