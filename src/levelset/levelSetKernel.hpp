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

# ifndef __BITPIT_LEVELSET_KERNEL_HPP__
# define __BITPIT_LEVELSET_KERNEL_HPP__

// Standard Template Library
# include <array>
# include <memory>
# include <vector>
# include <unordered_map>

# if BITPIT_ENABLE_MPI
# include <mpi.h>
# include "bitpit_communications.hpp"
# endif
# include "bitpit_common.hpp"

namespace bitpit{

namespace adaption{
    struct Info;
}

class VolumeKernel ;

class LevelSetObject ;
class LevelSetSignPropagator;

class LevelSetKernel{

    private:
    mutable std::unordered_map<long, std::array<double,3>> m_cellCentroids; /**< Cached cell center coordinates*/

    protected:
    VolumeKernel*                               m_mesh;        /**< Pointer to underlying mesh*/
# if BITPIT_ENABLE_MPI
    MPI_Comm                                    m_communicator; /**< MPI communicator */
# endif

    public:
    virtual ~LevelSetKernel() ;
    LevelSetKernel() ;
    LevelSetKernel( VolumeKernel *) ;

    virtual VolumeKernel *                      getMesh() const;

    const std::array<double,3> &                computeCellCentroid(long) const;
    virtual double                              computeCellIncircle(long) const;
    virtual double                              computeCellCircumcircle(long) const;

    virtual void                                clearGeometryCache();
    virtual void                                updateGeometryCache(const std::vector<adaption::Info> &);

    virtual bool                                intersectCellPlane(long, const std::array<double,3> &, const std::array<double,3> &, double);

    BITPIT_DEPRECATED(bool                      isPointInCell(long, const std::array<double,3> &) const);
    BITPIT_DEPRECATED(double                    isCellInsideBoundingBox(long, const std::array<double,3> &, const std::array<double,3> & ) const);

# if BITPIT_ENABLE_MPI
    void                                        initializeCommunicator();
    MPI_Comm                                    getCommunicator() const;
    void                                        freeCommunicator();
    bool                                        isCommunicatorSet() const;

    std::unique_ptr<DataCommunicator>           createDataCommunicator() const;
# endif

    virtual std::unique_ptr<LevelSetSignPropagator>     createSignPropagator() const ;

};

}

#endif
