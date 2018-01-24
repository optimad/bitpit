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
	SurfUnstructured(int patch_dim, int space_dim);
	SurfUnstructured(const int &id, int patch_dim, int space_dim);
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
        unsigned short importSTL(const std::string &, int PIDOffset = 0, bool PIDSquash = false);
        unsigned short importSTL(const std::string &, const bool &, int PIDOffset = 0, bool PIDSquash = false, std::unordered_map<int, std::string> *PIDNames = nullptr);
        unsigned short exportSTL(const std::string &, const bool &, bool flag = true);
        unsigned short exportSTL(const std::string &, const bool &, const bool &, bool flag, std::unordered_map<int, std::string> *PIDNames = nullptr);
        unsigned short importDGF(const std::string &, int PIDOffset = 0, bool PIDSquash = false);
        unsigned short exportDGF(const std::string &);

protected:
	SurfUnstructured(const SurfUnstructured &other) = default;

	int _getDumpVersion() const;
	void _dump(std::ostream &stream) const;
	void _restore(std::istream &stream);

	static ElementType getSTLFacetType(int nFacetVertices);
	static ElementType getDGFFacetType(int nFacetVertices);

	unsigned short importSTL(STLObj &STL, int PIDOffset, bool PIDSquash, std::unordered_map<int, std::string> *PIDNames = nullptr);

	unsigned short exportSTLSingle(const std::string &, const bool &, bool flag = true);
	unsigned short exportSTLMulti(const std::string &, bool flag = true, std::unordered_map<int, std::string> *PIDNames = nullptr);

private:

};

}

#endif
