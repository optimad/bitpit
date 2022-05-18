/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

#include "bitpit_patchkernel.hpp"
#include "bitpit_lineunstructured.hpp"

namespace bitpit {

class SurfUnstructured : public SurfaceKernel {

public:
    using PatchKernel::locatePoint;

    // Constructors
#if BITPIT_ENABLE_MPI==1
    SurfUnstructured(MPI_Comm communicator, std::size_t haloSize = 1);
    SurfUnstructured(int dimension, MPI_Comm communicator, std::size_t haloSize = 1);
    SurfUnstructured(int id, int dimension, MPI_Comm communicator, std::size_t haloSize = 1);
    SurfUnstructured(std::istream &stream, MPI_Comm communicator, std::size_t haloSize = 1);
#else
    SurfUnstructured();
    SurfUnstructured(int dimension);
    SurfUnstructured(int id, int dimension);
    SurfUnstructured(std::istream &stream);
#endif

    // Clone
    std::unique_ptr<PatchKernel> clone() const override;

    // Setters
    void setExpert(bool expert);

    // Search algorithms
    long locatePoint(const std::array<double, 3> &point) const override;

    // Evaluations
    void extractEdgeNetwork(LineUnstructured &net);

    // I/O routines
    int importSTL(const std::string &filename, int PIDOffset = 0, bool PIDSquash = false);
    int importSTL(const std::string &filename, bool isBinary, int PIDOffset = 0, bool PIDSquash = false, std::unordered_map<int, std::string> *PIDNames = nullptr);
    int importSTL(const std::string &filename, STLReader::Format format, bool joinFactes, int PIDOffset = 0, bool PIDSquash = false, std::unordered_map<int, std::string> *PIDNames = nullptr);
    int exportSTL(const std::string &filename, bool isBinary);
    int exportSTL(const std::string &filename, bool isBinary, bool isMulti, std::unordered_map<int, std::string> *PIDNames = nullptr);
    int importDGF(const std::string &filename, int PIDOffset = 0, bool PIDSquash = false);
    int importDGF(const std::string &filename, bool joinFactes, int PIDOffset = 0, bool PIDSquash = false);
    int exportDGF(const std::string &filename);

protected:
    int _getDumpVersion() const override;
    void _dump(std::ostream &stream) const override;
    void _restore(std::istream &stream) override;

    static ElementType getDGFFacetType(int nFacetVertices);

    int exportSTLSingle(const std::string &name, bool isBinary);
    int exportSTLMulti(const std::string &name, std::unordered_map<int, std::string> *PIDNames = nullptr);

};

}

#endif
