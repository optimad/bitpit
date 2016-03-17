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

#ifndef __BITPIT_SURFUNSTRUCTURED_HPP__
#define __BITPIT_SURFUNSTRUCTURED_HPP__

#include <array>
#include <vector>

#include "bitpit_IO.hpp"
#include "bitpit_patchkernel.hpp"

namespace bitpit {

class SurfUnstructured : public SurfaceKernel {

public:
	using PatchKernel::isPointInside;
	using PatchKernel::locatePoint;

        // Constructors
	SurfUnstructured(const int &id, int patch_dim = 2, int space_dim = 3);

        // Setters
	void setExpert(bool expert);

        // Modifiers
        void buildAdjacencies(void);
        void updateAdjacencies(const std::vector<long>&);

        // Search algorithms
        bool isPointInside(const std::array<double, 3> &point);
        long locatePoint(const std::array<double, 3> &point);

        // Evaluations
        void extractEdgeNetwork(SurfUnstructured &);

        // I/O routines
        unsigned short importSTL(const std::string &, const bool &);
        unsigned short exportSTL(const std::string &, const bool &, bool flag = true);
        unsigned short importDGF(const std::string &);
        unsigned short exportDGF(const std::string &);

protected:
	const std::vector<Adaption::Info> _updateAdaption(bool trackAdaption);
	bool _markCellForRefinement(const long &id);
	bool _markCellForCoarsening(const long &id);
	bool _enableCellBalancing(const long &id, bool enabled);

private:

};

}

#endif
