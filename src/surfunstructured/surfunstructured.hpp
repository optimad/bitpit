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

#include "bitpit_IO.hpp"
#include "bitpit_patchkernel.hpp"

namespace bitpit {

class SurfUnstructured : public SurfaceKernel {

public:
    using PatchKernel::locatePoint;

    // Constructors
#if BITPIT_ENABLE_MPI==1
    SurfUnstructured(MPI_Comm communicator);
    SurfUnstructured(int patch_dim, int space_dim, MPI_Comm communicator);
    SurfUnstructured(int id, int patch_dim, int space_dim, MPI_Comm communicator);
    SurfUnstructured(std::istream &stream, MPI_Comm communicator);
#else
    SurfUnstructured();
    SurfUnstructured(int patch_dim, int space_dim);
    SurfUnstructured(int id, int patch_dim, int space_dim);
    SurfUnstructured(std::istream &stream);
#endif

    // Clone
    std::unique_ptr<PatchKernel> clone() const override;

    // Setters
    void setExpert(bool expert);

    // Search algorithms
    long locatePoint(const std::array<double, 3> &point) const override;

    // Evaluations
    void extractEdgeNetwork(SurfUnstructured &net);

    // I/O routines
    int importSTL(const std::string &filename, int PIDOffset = 0, bool PIDSquash = false);
    int importSTL(const std::string &filename, bool isBinary, int PIDOffset = 0, bool PIDSquash = false, std::unordered_map<int, std::string> *PIDNames = nullptr);
    int exportSTL(const std::string &filename, bool isBinary);
    int exportSTL(const std::string &filename, bool isBinary, bool isMulti, std::unordered_map<int, std::string> *PIDNames = nullptr);
    int importDGF(const std::string &filename, int PIDOffset = 0, bool PIDSquash = false);
    int exportDGF(const std::string &filename);

protected:
    SurfUnstructured(const SurfUnstructured &other) = default;

    int _getDumpVersion() const override;
    void _dump(std::ostream &stream) const override;
    void _restore(std::istream &stream) override;

    static ElementType getDGFFacetType(int nFacetVertices);

    int importSTL(const std::string &name, STLReader::Format format, int PIDOffset, bool PIDSquash, std::unordered_map<int, std::string> *PIDNames = nullptr);

    int exportSTLSingle(const std::string &name, bool isBinary);
    int exportSTLMulti(const std::string &name, std::unordered_map<int, std::string> *PIDNames = nullptr);

};

}

#endif
