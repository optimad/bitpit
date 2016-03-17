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

#ifndef __BITPIT_SURFACE_KERNEL_HPP__
#define __BITPIT_SURFACE_KERNEL_HPP__

#include "patch_kernel.hpp"

namespace bitpit {

class SurfaceKernel : public PatchKernel {

public:
        // Static constant
        static const unsigned short SELECT_TRIANGLE;
        static const unsigned short SELECT_QUAD;
        static const unsigned short SELECT_ALL;
        static const std::map<ElementInfo::Type, unsigned short>     m_selectionTypes;

        // Types definitions
        typedef double (SurfaceKernel::*eval_f_)(const long&, int&);

	SurfaceKernel(const int &id, const int &patch_dim, const int &space_dim, bool expert);

        int getSpaceDimension(void) const;

	virtual ~SurfaceKernel();
        virtual std::array<double, 3> evalCellCentroid(const long &);
        virtual double evalCellArea(const long &);
        virtual double evalEdgeLength(const long&, const int&);
        virtual double evalMinEdgeLength(const long &, int &);
        virtual double evalMaxEdgeLength(const long &, int &);
        virtual double evalAngleAtVertex(const long&, const int&);
        virtual double evalMinAngleAtVertex(const long&, int &);
        virtual double evalMaxAngleAtVertex(const long&, int &);
        virtual double evalAspectRatio(const long&, int&);
        virtual std::array<double, 3> evalFacetNormal(const long&);
        std::array<double, 3> evalEdgeNormal(const long&, const int&);
        virtual std::array<double, 3> evalVertexNormal(const long&, const int&);
        double evalCellSize(const long &id);

        void displayQualityStats(ostream&, unsigned int padding = 0);
        vector<double> computeHistogram(eval_f_, std::vector<double>&, long&, int n_int = 8, unsigned short mask = SELECT_ALL);

private:
        bool compareSelectedTypes(const unsigned short &, const ElementInfo::Type &);
        void displayHistogram(const long&, const std::vector<double>&, const std::vector<double>&, const std::string&, std::ostream&, unsigned int padding = 0);

protected:
        int                     m_spaceDim;
        
};

}

#endif
