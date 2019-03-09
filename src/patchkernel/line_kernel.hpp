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

#ifndef __BITPIT_LINE_KERNEL_HPP__
#define __BITPIT_LINE_KERNEL_HPP__

#include "patch_kernel.hpp"

namespace bitpit {

class LineKernel : public PatchKernel {

public:
    int getVolumeCodimension() const override;
    int getSurfaceCodimension() const override;
    int getLineCodimension() const override;
    int getPointCodimension() const override;

    virtual double evalCellLength(long id) const;
    double evalCellSize(long id) const override;
    virtual std::array<double, 3> evalCellNormal(long id, const std::array<double, 3> &orientation = {{0., 0., 1.}}) const;

private:
    void initialize();

protected:
#if BITPIT_ENABLE_MPI==1
    LineKernel(MPI_Comm communicator, std::size_t haloSize, bool expert);
    LineKernel(int dimension, MPI_Comm communicator, std::size_t haloSize, bool expert);
    LineKernel(int id, int dimension, MPI_Comm communicator, std::size_t haloSize, bool expert);
#else
    LineKernel(bool expert);
    LineKernel(int dimension, bool expert);
    LineKernel(int id, int dimension, bool expert);
#endif

};

}

#endif
