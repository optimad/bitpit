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

# ifndef __BITPIT_LEVELSET_KERNEL_HPP__
# define __BITPIT_LEVELSET_KERNEL_HPP__

// Standard Template Library
# include <array>
# include <vector>
# include <unordered_map>

# if BITPIT_ENABLE_MPI
# include <mpi.h>
# endif

namespace bitpit{

namespace adaption{
    class Info ;
}

class VolumeKernel ;

class LevelSetObject ;

class LevelSetKernel{

    private:
    std::unordered_map<long, std::array<double,3>> m_cellCentroids; /**< Cached cell center coordinates*/

    protected:
    VolumeKernel*                               m_mesh;        /**< Pointer to underlying mesh*/
# if BITPIT_ENABLE_MPI
    MPI_Comm                                    m_commMPI;     /**< MPI communicator */
# endif

    public:
    virtual ~LevelSetKernel() ;
    LevelSetKernel() ;
    LevelSetKernel( VolumeKernel *) ;

    VolumeKernel*                               getMesh() const;

    virtual double                              computeSizeNarrowBandFromLS(LevelSetObject*, const bool &);
    virtual double                              computeRSearchFromCell(long) = 0;

    void                                        clearGeometryCache();
    void                                        updateGeometryCache(const std::vector<adaption::Info> &);

    const std::array<double,3> &                computeCellCentroid(long);
    double                                      isCellInsideBoundingBox(long, std::array<double,3> , std::array<double,3> );

# if BITPIT_ENABLE_MPI
    MPI_Comm                                    getCommunicator() const;
    void                                        freeCommunicator();
    bool                                        isCommunicatorSet() const;
    bool                                        assureMPI() ;
# endif

};

}

#endif
