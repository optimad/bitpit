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

#ifndef __BITPIT_VOLOCTREE_MAPPER_HPP__
#define __BITPIT_VOLOCTREE_MAPPER_HPP__

#if BITPIT_ENABLE_MPI
#    include <mpi.h>
#endif
#include <string>
#include <vector>
#include <unordered_map>

#include "volume_mapper.hpp"
#include "voloctree.hpp"
#include "voloctree_mapper.hpp"

namespace bitpit {

class VolOctreeMapper: public VolumeMapper {

public:
#if BITPIT_ENABLE_MPI
    VolOctreeMapper(VolOctree *referencePatch, VolOctree *mappedPatch, MPI_Comm communicator = MPI_COMM_WORLD);
#else
    VolOctreeMapper(VolOctree *referencePatch, VolOctree *mappedPatch);
#endif

    void adaptionPrepare(const std::vector<adaption::Info> &adaptionInfo, bool reference = true);
    void adaptionAlter(const std::vector<adaption::Info> &adaptionInfo, bool reference = true, bool inverseFilled = false);
    void adaptionCleanup();

#if BITPIT_ENABLE_MPI
    bool checkPartition();
    void clearPartitionMapping();
    std::map<int, std::vector<long>> getReceivedMappedIds();
    std::map<int, std::vector<long>> getSentMappedIds();
    std::map<int, std::vector<long>> getSentReferenceIds();
#endif

private:
    /**
     * \class OctantIR
     * \ingroup volumepatches
     *
     * \brief The OctantIR is an internal structure used by a mapper object to
     * store info about an octant.
     */
    struct OctantIR {
        Octant octant;  /**< Octant object. */
        long id;        /**< Octant local id. */
        long globalId;  /**< Octant global id. */
        int rank;       /**< Octant rank. */

        /**
         * Constructor.
         */
        OctantIR()
        {
        };

        /**
         * Constructor.
         *
         * \param[in] _octant Octant object
         * \param[in] _id Local octant id
         * \param[in] _globalId Global octant id
         * \param[in] _rank Octant rank
         */
        OctantIR(Octant _octant, long _id, long _globalId = -1, int _rank = -1)
        {
            octant = _octant;
            id = _id;
            globalId = _globalId;
            rank = _rank;
        };
    };

#if BITPIT_ENABLE_MPI
    /**
     * \class PartitionMapper
     * \ingroup volumepatches
     *
     * \brief The PartitionMapper is an internal structure of a mapper object
     * used to store info about the partitioning of two mapped patches.
     */
    struct PartitionMapper {
        std::vector<OctantIR> list_rec_octantIR_before; /**< List of received Octants from lower rank processes. */
        std::vector<OctantIR> list_rec_octantIR_after;  /**< List of received Octants from higher rank processes. */

        std::unordered_map<int, std::unordered_map<long, OctantIR *>> map_rank_rec_octantIR;       /**< List of received Octants from other processes. */
        std::unordered_map<int, std::unordered_map<long, mapping::Info>> map_rank_inverseMapping;  /**< Inverse mapper terms related to processes different from the local rank. */
        std::unordered_map<int, std::unordered_map<long, mapping::Info>> map_rank_previousMapping; /**< Temporary inverse mapper for each processes used during a patch adaptation. */

        std::vector<OctantIR> list_sent_octantIR;   /**< List of sent Octants to other processes. */

        std::vector<uint64_t> partitionFDReference;  /**< First descendant partitioning structure of reference patch. */
        std::vector<uint64_t> partitionLDReference;  /**< Last descendant partitioning structure of reference patch. */
        std::vector<uint64_t> partitionFDMapped;     /**< First descendant partitioning structure of mapped patch. */
        std::vector<uint64_t> partitionLDMapped;     /**< Last descendant partitioning structure of mapped patch. */

        /**
         * Default constructor.
         */
        PartitionMapper(){};

        /**
         * Clear members.
         */
        void clear(){
            list_rec_octantIR_before.clear();
            list_rec_octantIR_after.clear();
            list_sent_octantIR.clear();
            map_rank_rec_octantIR.clear();
            map_rank_inverseMapping.clear();
            map_rank_previousMapping.clear();
            partitionFDReference.clear();
            partitionLDReference.clear();
            partitionFDMapped.clear();
            partitionLDMapped.clear();
        }

        /**
         * Clear only list members.
         */
        void clearLists(){
            list_rec_octantIR_before.clear();
            list_rec_octantIR_after.clear();
            list_sent_octantIR.clear();
        }
    };
#endif

#if BITPIT_ENABLE_MPI
    PartitionMapper m_partitionIR;  /**< Partitioning info structure. */
#endif

    void _mapMeshes(bool fillInverse);
    void _mapMeshesSamePartition(const std::vector<OctantIR> * octantsReference, const std::vector<OctantIR> * octantsMapped,
                                 bool fillInverse, long *indRef);

    void _mappingAdaptionReferenceUpdate(const std::vector<adaption::Info> &adaptionInfo, bool fillInverse = false);
    void _mappingAdaptionMappedUpdate(const std::vector<adaption::Info> &adaptionInfo);

#if BITPIT_ENABLE_MPI
    void clearPartitionMappingLists();

    void _mapMeshPartitioned(bool fillInverse);

    bool _recoverPartition();

    void _communicateInverseMapper(const std::vector<OctantIR> *octantsIRMapped,
                                   const std::unordered_map<long, mapping::Info> &inverseGlobalMapping,
                                   std::unordered_map<long, mapping::Info> *inverseLocalMapping);

    void _communicateInverseMapperBack();

    void _communicateMappedAdaptionInfo(const std::vector<adaption::Info> &adaptionInfoMap, std::vector<adaption::Info> *adaptionInfoRef);
#endif

};

}

#endif
