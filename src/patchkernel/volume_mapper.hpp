/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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

#ifndef __BITPIT_VOLUME_MAPPER_HPP__
#define __BITPIT_VOLUME_MAPPER_HPP__

#if BITPIT_ENABLE_MPI
#    include <mpi.h>
#endif
#include <string>
#include <vector>
#include <unordered_map>

#include "bitpit_patchkernel.hpp"

namespace bitpit {

/*!
    \ingroup mapping
    \brief The namespace 'mapping' contains structures for working with the MapperVolOctree class.
 */
namespace mapping
{

/*!
* Type of mapping relationship between elements.
*/
enum Type {
    TYPE_UNKNOWN = 0,
    TYPE_REFINEMENT,
    TYPE_COARSENING,
    TYPE_RENUMBERING
};

/*!
* Entity involved in a mapping relationship.
*/
enum Entity {
    ENTITY_UNKNOWN = -1,
    ENTITY_CELL,
};

/**
 * \class Info
 * \ingroup volumepatches
 *
 * \brief The Info is the structure to store info about mapping between elements.
 *
 */
struct Info
{
    /**
     * Default constructor.
     */
    Info()
    : type(TYPE_UNKNOWN), entity(ENTITY_UNKNOWN)
    {
    }

    /**
     * Custom constructor.
     * \param[in] user_type Type of mapping item
     * \param[in] user_entity Mapped entity
     */
    Info(Type user_type, Entity user_entity)
    : type(user_type), entity(user_entity)
    {
    }

    Type type;                  /**< Type of mapping item. */
    Entity entity;              /**< Mapped entity. */
    std::vector<long> ids;      /**< Local ids of the mapped elements. */
#if BITPIT_ENABLE_MPI
    std::vector<int> ranks;     /**< Rank owners of the mapped elements. */
#endif

};

}

class VolumeMapper {

public:
    virtual ~VolumeMapper() = default;

    void clear();
    void clearMapping();
    void clearInverseMapping();

    const bitpit::PiercedStorage<mapping::Info> & getMapping();
    const bitpit::PiercedStorage<mapping::Info> & getInverseMapping();

    void initialize(bool fillInv = false);

    virtual void adaptionPrepare(const std::vector<adaption::Info> &infoAdapt, bool reference = true) = 0;
    virtual void adaptionAlter(const std::vector<adaption::Info> &infoAdapt, bool reference = true, bool fillInv = false) = 0;
    virtual void adaptionCleanup() = 0;

#if BITPIT_ENABLE_MPI
    virtual bool checkPartition() = 0;
    virtual std::map<int, std::vector<long> > getReceivedMappedIds() = 0;
    virtual std::map<int, std::vector<long> > getSentMappedIds() = 0;
    virtual std::map<int, std::vector<long> > getSentReferenceIds() = 0;
#endif

protected:
    VolumeKernel *m_referencePatch;    /**< Pointer to reference mesh.*/
    VolumeKernel *m_mappedPatch;       /**< Pointer to mapped mesh.*/

    PiercedStorage<mapping::Info> m_mapping;                     /**< Mapping info for each cell of reference mesh.
                                                                     The mapping info is treated as a set of adaption info related to
                                                                     an adaption of the mapped mesh toward the reference mesh. */

    PiercedStorage<mapping::Info> m_inverseMapping;              /**< Inverse mapping info for each cell of mapped mesh. */

    std::unordered_map<long, mapping::Info> m_previousMapping;   /**< Temporary mapping used during a mesh adaptation. */

#if BITPIT_ENABLE_MPI
    MPI_Comm m_communicator; /**< MPI communicator */
    int m_rank;              /**< Local rank of process. */
    int m_nProcs;            /**< Number of processes. */
#endif

#if BITPIT_ENABLE_MPI
    VolumeMapper(VolumeKernel *referencePatch, VolumeKernel *mappedPatch, MPI_Comm communicator = MPI_COMM_WORLD);
#else
    VolumeMapper(VolumeKernel *referencePatch, VolumeKernel *mappedPatch);
#endif

    virtual void _mapMeshes(bool fillInverse) = 0;

#if BITPIT_ENABLE_MPI
    void initializeCommunicator(MPI_Comm communicator);
    MPI_Comm getCommunicator() const;
    bool isCommunicatorSet() const;
    void freeCommunicator();
#endif

};

}
#endif
