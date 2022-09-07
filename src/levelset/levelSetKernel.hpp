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

# include "bitpit_CG.hpp"
# if BITPIT_ENABLE_MPI
# include <mpi.h>
# include "bitpit_communications.hpp"
# endif
# include "bitpit_common.hpp"
# include "bitpit_operators.hpp"
# include "bitpit_patchkernel.hpp"

# include "levelSetObject.hpp"
# include "levelSetSignPropagator.hpp"

namespace bitpit{

namespace adaption{
    struct Info;
}

class VolumeKernel ;

class LevelSetObject ;
class LevelSetSignPropagator;

class LevelSetKernel {

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

    virtual std::array<double, 3>               computeCellCentroid(long) const = 0;
    virtual double                              computeCellTangentRadius(long) const = 0;
    virtual double                              computeCellBoundingRadius(long) const = 0;

    virtual void                                update(const std::vector<adaption::Info> &);

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

template<typename CellCacheEntry>
class LevelSetCachedKernel : public LevelSetKernel {

    private:
    mutable std::unordered_map<long, CellCacheEntry>  m_cellCache;  /**< Cell cache */

    protected:
    const CellCacheEntry &                      computeCellCacheEntry( long id ) const;

    virtual void                                pruneCellCache(const std::vector<adaption::Info> &);

    public:
    LevelSetCachedKernel() ;
    LevelSetCachedKernel( VolumeKernel *) ;

    void                                        update(const std::vector<adaption::Info> &) override;

};

}

// Include template implementations
#include "levelSetKernel.tpp"

#endif
