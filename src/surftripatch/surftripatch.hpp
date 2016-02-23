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

#ifndef __BITPIT_SURFTRIPATCH_HPP__
#define __BITPIT_SURFTRIPATCH_HPP__

#include <array>
#include <vector>

#include "bitpit_patch.hpp"

namespace bitpit {

class SurfTriPatch : public Patch {

public:

	using Patch::isPointInside;
	using Patch::locatePoint;

        // Types definitions
        typedef double (SurfTriPatch::*eval_f_)(const long&, int&);

        // Static constant
        static const unsigned short SELECT_TRIANGLE;
        static const unsigned short SELECT_QUAD;
        static const unsigned short SELECT_ALL;

	SurfTriPatch(const int &id);

	~SurfTriPatch();

	void setExpert(bool expert);

	double evalCellVolume(const long &id);
	double evalCellSize(const long &id);

	double evalInterfaceArea(const long &id);
	std::array<double, 3> evalInterfaceNormal(const long &id);

	bool isPointInside(const std::array<double, 3> &point);
	long locatePoint(const std::array<double, 3> &point);
        void buildAdjacencies(void);
        void updateAdjacencies(const std::vector<long>&);

        //TODO: double evalCellArea(const long &);
        double evalEdgeLength(const long&, const int&);
        double evalMinEdgeLength(const long &, int &);
        double evalMaxEdgeLength(const long &, int &);
        double evalAngleAtVertex(const long&, const int&);
        double evalMinAngleAtVertex(const long&, int &);
        double evalMaxAngleAtVertex(const long&, int &);
        array<double, 3> evalFacetNormal(const long&);
        double evalAspectRatio(const long&, int&);
        double evalFacetArea(const long&);
        vector<double> computeHistogram(
            eval_f_                     ,
            vector<double>              &,
            long                        &,
            int                          n_int = 8,
            unsigned short               mask = SELECT_ALL
        );
        void displayQualityStats(
            ostream                     &,
            unsigned int                 padding = 0
        );

protected:
	const std::vector<Adaption::Info> _update(bool trackAdaption);
	bool _markCellForRefinement(const long &id);
	bool _markCellForCoarsening(const long &id);
	bool _enableCellBalancing(const long &id, bool enabled);

private:
        void displayHistogram(
            const long                  &,
            const vector<double>        &,
            const vector<double>        &,
            const std::string           &,
            ostream                     &,
            unsigned int                 padding = 0
        );
        bool compareSelectedTypes(const unsigned short &, const ElementInfo::Type &);
        static const std::map<ElementInfo::Type, unsigned short>     m_selectionTypes;
};

}

#endif
