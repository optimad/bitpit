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

#include "voloctree_mapper.hpp"

namespace bitpit {

/**
 * \class VolOctreeMapper
 * \ingroup volumepatches
 *
 * \brief The VolOctreeMapper is the class to map two meshes of class VolOctree.
 *
 * The VolOctreeMapper allows to map meshes of class VolOctree. The meshes are
 * defined as a reference mesh and a mapped mesh.
 *
 * The object can provide a direct mapper and a inverse mapper between only the
 * Cells of the meshes.
 *
 * The information given by a mapper object is analogous to an adaptation info
 * to adapt the mapped mesh to the reference one (or vice versa in case of
 * inverse mapper).
 *
 * The two meshes have to be imperatively linked at declaration of the mapper
 * object.
 *
 * To compute the mapper the first time call initialize method. Then, if the
 * reference mesh OR the mapped mesh (one at a time) is adapted, the mapper
 * can be adapted together by passing the adaptation information to the prepare
 * and alter methods. A load-balancing procedure is not allowed, i.e. the mapper
 * has to be entirely recomputed after a load-balance of a mesh.
 *
 * Note. The domain of the two meshes has to be the same, i.e. the meshes must
 * have same origin and length.
 */

/**
 * Constructor.
 *
 * \param[in] referencePatch is the reference mesh
 * \param[in] mappedPatch is the mapped mesh
 */
#if BITPIT_ENABLE_MPI
/**
 * \param[in] communicator is the MPI communicator
 */
VolOctreeMapper::VolOctreeMapper(bitpit::VolOctree *referencePatch, bitpit::VolOctree *mappedPatch, MPI_Comm communicator)
    : VolumeMapper(referencePatch, mappedPatch, communicator)
#else
VolOctreeMapper::VolOctreeMapper(bitpit::VolOctree *referencePatch, bitpit::VolOctree *mappedPatch)
    : VolumeMapper(referencePatch, mappedPatch)
#endif
{
}

#if BITPIT_ENABLE_MPI
/**
 * Clear partition mapping members
 */
void VolOctreeMapper::clearPartitionMapping()
{
    m_partitionIR.clear();
}

/**
 * Clear only partition mapping lists
 */
void VolOctreeMapper::clearPartitionMappingLists()
{
    m_partitionIR.clearLists();
}
#endif

/**
 * Prepare the mapper for an adaption of the reference OR the mapped mesh.
 *
 * The adaptation of only one mesh at time is allowed. Calling this method
 * is mandatory to perform an update of the mapper after an adaptation.
 *
 * \param[in] adaptionInfo are the adaptation info that describe the changes
 * of the patch
 * \param[in] reference if set to true the reference mesh will be adapted,
 * is set to false the mapped one will be adapted
 */
void VolOctreeMapper::adaptionPrepare(const std::vector<adaption::Info> &adaptionInfo, bool reference)
{
    m_previousMapping.clear();
    PiercedStorage<mapping::Info> *pmapper;

    if (reference) {
        pmapper = &m_mapping;
    } else {
        pmapper = &m_inverseMapping;
    }

    m_previousMapping.clear();
#if BITPIT_ENABLE_MPI
    m_partitionIR.map_rank_previousMapping.clear();
#endif
    for (const adaption::Info &info : adaptionInfo) {
        if (info.type == adaption::Type::TYPE_PARTITION_SEND ||
                info.type == adaption::Type::TYPE_PARTITION_RECV ||
                info.type == adaption::Type::TYPE_PARTITION_NOTICE) {
            throw std::runtime_error ("Mapper: type of adation not supported : " + std::to_string(info.type));
        } else {
            if (info.type != adaption::Type::TYPE_DELETION && info.type != adaption::Type::TYPE_CREATION) {
                for (long id : info.previous) {
                    m_previousMapping[id] = (*pmapper).at(id);
                }
            }
        }
    }
}

/**
 * Update the mapper after an adaption of the reference OR the mapped mesh.
 *
 * The adaptationPrepare method of the mapper has to called before the
 * adaptation of the mesh.
 *
 * \param[in] adaptionInfo are the adaptation info that describe the changes of
 * the patch
 * \param[in] reference if set to true the reference mesh will be adapted,
 * is set to false the mapped one will be adapted
 * \param[in] inverseFilled if set to true the inverse mapped was filled
 * during the mesh mapper computing. If the adapted mesh is the mapped one,
 * the inverse mapper has to be necessarily filled during the first computation
 * in initialization)
 */
void VolOctreeMapper::adaptionAlter(const std::vector<adaption::Info> &adaptionInfo, bool reference, bool inverseFilled)
{
    if (reference) {
        _mappingAdaptionReferenceUpdate(adaptionInfo, inverseFilled);
    } else {
        _mappingAdaptionMappedUpdate(adaptionInfo);
    }
}

/**
 * Clear the mapper updating internal structures.
 */
void VolOctreeMapper::adaptionCleanup()
{
    m_previousMapping.clear();
}

/**
 * Update the mapper after an adaption of the reference mesh.
 *
 * \param[in] adaptionInfo are the adaptation info that describe the changes of
 * the patch
 * \param[in] inverseFilled if set to true the inverse mapped was filled
 * during the mesh mapper computing. If the adapted mesh is the mapped one,
 * the inverse mapper has to be necessarily filled during the first computation
 * in initialization)
 */
void VolOctreeMapper::_mappingAdaptionReferenceUpdate(const std::vector<adaption::Info> &adaptionInfo, bool inverseFilled)
{
    VolOctree *adaptedPatch = static_cast<VolOctree*>(m_referencePatch);
    VolOctree *mappedPatch  = static_cast<VolOctree*>(m_mappedPatch);

    PiercedStorage<mapping::Info> *mappingAdapted = &m_mapping;
    PiercedStorage<mapping::Info> *mappingMapped  = &m_inverseMapping;

    bool changedPartition = false;
#if BITPIT_ENABLE_MPI
    bool checkPart = true;
    checkPart = checkPartition();
    if (!checkPart) {
        changedPartition = _recoverPartition();
    }
#endif

    if (!changedPartition) {
        for (const adaption::Info &info : adaptionInfo) {
            if (info.type == adaption::Type::TYPE_PARTITION_SEND ||
                    info.type == adaption::Type::TYPE_PARTITION_RECV ||
                    info.type == adaption::Type::TYPE_PARTITION_NOTICE) {
                throw std::runtime_error ("Mapper: type of adation not supported : " + std::to_string(info.type));
            } else if (info.type == adaption::Type::TYPE_DELETION ||
                            info.type == adaption::Type::TYPE_CREATION) {
                // Do nothing
            } else {
                if (inverseFilled) {
                    // Erase previous id in mappingMapped elements
                    for (long idprevious : info.previous) {
                        const std::vector<long> &mappedIds = m_previousMapping[idprevious].ids;

                        std::size_t imapped = 0;
                        for (long idp : mappedIds) {
#if BITPIT_ENABLE_MPI
                            if (checkPart || m_previousMapping[idprevious].ranks[imapped] == m_rank) {
#endif
                                std::vector<long>::iterator it = std::find((*mappingMapped)[idp].ids.begin(), (*mappingMapped)[idp].ids.end(), idprevious);
                                if (it != (*mappingMapped)[idp].ids.end()) {
                                    (*mappingMapped)[idp].ids.erase(it);
#if BITPIT_ENABLE_MPI
                                    int dist = std::distance((*mappingMapped)[idp].ids.begin(), it);
                                    (*mappingMapped)[idp].ranks.erase((*mappingMapped)[idp].ranks.begin()+dist);
#endif
                                }
#if BITPIT_ENABLE_MPI
                            } else {
                                mapping::Info &inverseMappingInfo = m_partitionIR.map_rank_inverseMapping[m_previousMapping[idprevious].ranks[imapped]][idp];
                                std::vector<long>::iterator it = std::find(inverseMappingInfo.ids.begin(), inverseMappingInfo.ids.end(), idprevious);
                                if (it != inverseMappingInfo.ids.end()) {
                                    int dist = std::distance(inverseMappingInfo.ids.begin(), it);
                                    inverseMappingInfo.ids.erase(it);
                                    inverseMappingInfo.ranks.erase(inverseMappingInfo.ranks.begin() + dist);
                                }
                            }
#endif
                            imapped++;
                        }
                    }
                }

                switch(info.type) {

                // Renumbering
                case adaption::Type::TYPE_RENUMBERING:
                {
                    long id = info.current[0];
                    mapping::Info &mappingInfo = (*mappingAdapted)[id];

                    mappingInfo.ids.clear();
                    mappingInfo.type = m_previousMapping[info.previous[0]].type;
                    mappingInfo.entity = mapping::Entity::ENTITY_CELL;
                    mappingInfo.ids = m_previousMapping[info.previous[0]].ids;
#if BITPIT_ENABLE_MPI
                    mappingInfo.ranks = m_previousMapping[info.previous[0]].ranks;
#endif

                    if (inverseFilled) {
                        const std::vector<long> &mappedIds = mappingInfo.ids;

                        std::size_t imapped = 0;
                        for (long idp : mappedIds) {
#if BITPIT_ENABLE_MPI
                            if (checkPart || mappingInfo.ranks[imapped] == m_rank) {
#endif
                                (*mappingMapped)[idp].ids.push_back(id);
#if BITPIT_ENABLE_MPI
                                (*mappingMapped)[idp].ranks.push_back(info.rank);
                            } else {
                                mapping::Info &inverseMappingInfo = m_partitionIR.map_rank_inverseMapping[mappingInfo.ranks[imapped]][idp];
                                inverseMappingInfo.ids.push_back(id);
                                inverseMappingInfo.ranks.push_back(info.rank);
                            }
#endif
                            imapped++;
                        }
                    }
                }
                break;

                // Refinement and coarsening
                case adaption::Type::TYPE_REFINEMENT:
                case adaption::Type::TYPE_COARSENING:

                    // Clear current mapper
                    for (long id : info.current) {
                        mapping::Info &mappingInfo = (*mappingAdapted)[id];

                        mappingInfo.ids.clear();
#if BITPIT_ENABLE_MPI
                        mappingInfo.ranks.clear();
#endif

                        mappingInfo.entity = mapping::Entity::ENTITY_CELL;
                    }

                    for (long idprevious : info.previous) {
                        const std::vector<long> &mappedIds = m_previousMapping[idprevious].ids;
#if BITPIT_ENABLE_MPI
                        const std::vector<int> &mappedRanks = m_previousMapping[idprevious].ranks;
#endif

                        std::size_t imapped = 0;
                        for (long mappedId : mappedIds) {
                            VolOctree::OctantInfo oinfoprev;
                            Octant * octm ;
                            uint64_t mortonmapped;
                            uint64_t mortonlastdescmapped;
                            uint8_t levelmapped;
#if BITPIT_ENABLE_MPI
                            int mappedRank;
                            if (checkPart || mappedRanks[imapped] == m_rank) {
#endif
                                oinfoprev = mappedPatch->getCellOctant(mappedId);
                                octm = mappedPatch->getOctantPointer(oinfoprev);
                                mortonmapped = mappedPatch->getTree().getMorton(octm);
                                mortonlastdescmapped = mappedPatch->getTree().getLastDescMorton(octm);
                                levelmapped = mappedPatch->getCellLevel(mappedId);
#if BITPIT_ENABLE_MPI
                                mappedRank = mappedRanks[imapped];
                            } else {
                                OctantIR *poct = m_partitionIR.map_rank_rec_octantIR[mappedRanks[imapped]][mappedId];
                                octm = &poct->octant;
                                mortonmapped = mappedPatch->getTree().getMorton(&poct->octant);
                                mortonlastdescmapped = mappedPatch->getTree().getLastDescMorton(&poct->octant);
                                levelmapped = mappedPatch->getTree().getLevel(&poct->octant);
                                mappedRank = mappedRanks[imapped];
                            }
#endif

                            for (long id : info.current) {
                                // Retrieve level of current cell
                                uint8_t level;
                                uint64_t morton, mortonlastdesc;
                                level = adaptedPatch->getCellLevel(id);
                                VolOctree::OctantInfo octantIfo = adaptedPatch->getCellOctant(id);
                                Octant *octant = adaptedPatch->getOctantPointer(octantIfo);
                                morton = adaptedPatch->getTree().getMorton(octant);
                                mortonlastdesc = adaptedPatch->getTree().getLastDescMorton(octant);

                                bool checkmorton = false;

                                checkmorton |= (morton >= mortonmapped && mortonlastdesc <= mortonlastdescmapped);
                                checkmorton |= (morton <= mortonmapped && mortonlastdesc >= mortonlastdescmapped);

                                if (checkmorton) {
                                    if (level == levelmapped) {
                                        mapping::Info &mappingInfo = (*mappingAdapted)[id];
                                        if (mappingInfo.ids.size() == 0) {
                                            mappingInfo.type = mapping::Type::TYPE_RENUMBERING;
                                            mappingInfo.ids.push_back(mappedId);
#if BITPIT_ENABLE_MPI
                                            mappingInfo.ranks.push_back(mappedRank);
#endif
                                        }
#if BITPIT_ENABLE_MPI
                                        if (inverseFilled) {
                                            if (checkPart || mappedRank == m_rank) {
#endif
                                                (*mappingMapped)[mappedId].type = mapping::Type::TYPE_RENUMBERING;
                                                (*mappingMapped)[mappedId].ids.clear();
                                                (*mappingMapped)[mappedId].ids.push_back(id);
#if BITPIT_ENABLE_MPI
                                                (*mappingMapped)[mappedId].ranks.push_back(info.rank);
                                            } else {
                                                mapping::Info &inverseMappingInfo = m_partitionIR.map_rank_inverseMapping[mappedRank][mappedId];
                                                inverseMappingInfo.type = mapping::Type::TYPE_RENUMBERING;
                                                inverseMappingInfo.ids.push_back(id);
                                                inverseMappingInfo.ranks.push_back(info.rank);
                                            }
                                        }
#endif
                                    }
                                    else if (level > levelmapped) {
                                        mapping::Info &mappingInfo = (*mappingAdapted)[id];
                                        if (mappingInfo.ids.size() == 0) {
                                            mappingInfo.type = mapping::Type::TYPE_REFINEMENT;
                                            mappingInfo.ids.push_back(mappedId);
#if BITPIT_ENABLE_MPI
                                            mappingInfo.ranks.push_back(mappedRank);
#endif
                                        }
#if BITPIT_ENABLE_MPI
                                        if (inverseFilled) {
                                            if (checkPart || mappedRank == m_rank) {
#endif
                                                if (std::find((*mappingMapped)[mappedId].ids.begin(), (*mappingMapped)[mappedId].ids.end(), id) == (*mappingMapped)[mappedId].ids.end()) {
                                                (*mappingMapped)[mappedId].ids.push_back(id);
                                                (*mappingMapped)[mappedId].type = mapping::Type::TYPE_COARSENING;
#if BITPIT_ENABLE_MPI
                                                (*mappingMapped)[mappedId].ranks.push_back(info.rank);
#endif
                                                }
#if BITPIT_ENABLE_MPI
                                            } else {
                                                mapping::Info &inverseMappingInfo = m_partitionIR.map_rank_inverseMapping[mappedRank][mappedId];
                                                inverseMappingInfo.type = mapping::Type::TYPE_COARSENING;
                                                inverseMappingInfo.ids.push_back(id);
                                                inverseMappingInfo.ranks.push_back(info.rank);
                                            }
                                        }
#endif
                                    } else if (level < levelmapped) {
                                        mapping::Info &mappingInfo = (*mappingAdapted)[id];
                                        if (std::find(mappingInfo.ids.begin(), mappingInfo.ids.end(), mappedId) == mappingInfo.ids.end()) {
                                            mappingInfo.type = mapping::Type::TYPE_COARSENING;
                                            mappingInfo.ids.push_back(mappedId);
#if BITPIT_ENABLE_MPI
                                            mappingInfo.ranks.push_back(mappedRank);
#endif
                                        }
#if BITPIT_ENABLE_MPI
                                        if (inverseFilled) {
                                            if (checkPart || mappedRank == m_rank ) {
#endif
                                                if ((*mappingMapped)[mappedId].ids.size() == 0) {
                                                (*mappingMapped)[mappedId].type = mapping::Type::TYPE_REFINEMENT;
                                                (*mappingMapped)[mappedId].ids.push_back(id);
#if BITPIT_ENABLE_MPI
                                                (*mappingMapped)[mappedId].ranks.push_back(info.rank);
#endif
                                                }
#if BITPIT_ENABLE_MPI
                                            } else {
                                                mapping::Info &inverseMappingInfo = m_partitionIR.map_rank_inverseMapping[mappedRank][mappedId];
                                                if (inverseMappingInfo.ids.size() == 0) {
                                                    inverseMappingInfo.type = mapping::Type::TYPE_REFINEMENT;
                                                    inverseMappingInfo.ids.push_back(id);
                                                    inverseMappingInfo.ranks.push_back(info.rank);
                                                }
                                            }
                                        }
#endif
                                    } // End level < levelreference
                                } // End checkmorton
                            } // End id current

                            imapped++;
                        } // End idmapped
                    } // End idprevious
                    break;

                // Default
                default:
                    break;

                } // End switch
            } // End if deletion
        } // End end info

#if BITPIT_ENABLE_MPI
        if (inverseFilled) {
            _communicateInverseMapperBack();
        }
#endif
    } else {
        initialize(inverseFilled);
    }
}

/**
 * Update the mapping after an adaption of the mapped mesh.
 *
 * \param[in] adaptionInfo are the adaptation info that describe the changes of
 * the patch
 */
void VolOctreeMapper::_mappingAdaptionMappedUpdate(const std::vector<adaption::Info> &adaptionInfo)
{
    VolOctree *adaptedPatch   = static_cast<VolOctree*>(m_mappedPatch);
    VolOctree *referencePatch = static_cast<VolOctree*>(m_referencePatch);

    PiercedStorage<mapping::Info> *mappingAdapted   = &m_inverseMapping;
    PiercedStorage<mapping::Info> *mappingReference = &m_mapping;

    bool changedPartition = false;
#if BITPIT_ENABLE_MPI
    bool checkPart = true;
    checkPart = checkPartition();
    if (!checkPart) {
        changedPartition = _recoverPartition();
    }
#endif

    if (!m_inverseMapping.getKernel()) {
        throw std::runtime_error("Updating mapping with adaption of mapped mesh possible only if inverse mapper is filled");
    }

    if (!changedPartition) {
        std::vector<adaption::Info> adaptionInfoRef;
#if BITPIT_ENABLE_MPI
        if (!checkPart) {
            _communicateMappedAdaptionInfo(adaptionInfo, &adaptionInfoRef);
        } else {
            adaptionInfoRef = adaptionInfo;
        }
#else
        adaptionInfoRef = adaptionInfo;
#endif

        for (const adaption::Info &info : adaptionInfoRef) {
            if (info.type == adaption::Type::TYPE_PARTITION_SEND ||
                    info.type == adaption::Type::TYPE_PARTITION_RECV ||
                    info.type == adaption::Type::TYPE_PARTITION_NOTICE) {
                throw std::runtime_error ("Mapper: type of adation not supported : " + std::to_string(info.type));
            } else if (info.type == adaption::Type::TYPE_DELETION ||
                    info.type == adaption::Type::TYPE_CREATION) {
                // Do nothing
            } else {
                // Erase previous id in mappingReference elements
                for (long idprevious : info.previous) {
                    std::vector<long> *mappedIds;
#if BITPIT_ENABLE_MPI
                    if (checkPart || info.rank == m_rank) {
#endif
                        mappedIds = &(m_previousMapping[idprevious].ids);
#if BITPIT_ENABLE_MPI
                    } else {
                        mappedIds = &(m_partitionIR.map_rank_previousMapping[info.rank][idprevious].ids);
                    }
#endif
                    for (long idp : *mappedIds) {
                        std::vector<long>::iterator it = std::find((*mappingReference)[idp].ids.begin(), (*mappingReference)[idp].ids.end(), idprevious);
                        if (it != (*mappingReference)[idp].ids.end()) {
                            (*mappingReference)[idp].ids.erase(it);
#if BITPIT_ENABLE_MPI
                            int dist = std::distance((*mappingReference)[idp].ids.begin(), it);
                            (*mappingReference)[idp].ranks.erase((*mappingReference)[idp].ranks.begin()+dist);
#endif
                        }
                    }
                }

                // Renumbering
                if (info.type == adaption::Type::TYPE_RENUMBERING) {
                    long id = info.current[0];
                    mapping::Info &mappingInfo = (*mappingAdapted)[id];

                    std::vector<long> *mappedIds;
#if BITPIT_ENABLE_MPI
                    if (checkPart || info.rank == m_rank) {
#endif
                        mappingInfo.ids.clear();
                        mappingInfo.type = m_previousMapping[info.previous[0]].type;
                        mappingInfo.entity = mapping::Entity::ENTITY_CELL;
                        mappingInfo.ids = m_previousMapping[info.previous[0]].ids;
                        mappedIds = &(mappingInfo.ids);
#if BITPIT_ENABLE_MPI
                        mappingInfo.ranks = m_previousMapping[info.previous[0]].ranks;
                    } else {
                        for (long _idp : info.previous) {
                            if (m_partitionIR.map_rank_inverseMapping[info.rank].count(_idp)) {
                                m_partitionIR.map_rank_inverseMapping[info.rank].erase(_idp);
                            }
                        }

                        m_partitionIR.map_rank_inverseMapping[info.rank][id].ids.clear();
                        m_partitionIR.map_rank_inverseMapping[info.rank][id].type = m_partitionIR.map_rank_previousMapping[info.rank][info.previous[0]].type;
                        m_partitionIR.map_rank_inverseMapping[info.rank][id].entity = mapping::Entity::ENTITY_CELL;
                        m_partitionIR.map_rank_inverseMapping[info.rank][id].ids = m_partitionIR.map_rank_previousMapping[info.rank][info.previous[0]].ids;
                        m_partitionIR.map_rank_inverseMapping[info.rank][id].ranks = m_partitionIR.map_rank_previousMapping[info.rank][info.previous[0]].ranks;
                        mappedIds = &(m_partitionIR.map_rank_inverseMapping[info.rank][id].ids);
                    }
#endif
                    for (long idp : *mappedIds) {
                        (*mappingReference)[idp].ids.push_back(id);
#if BITPIT_ENABLE_MPI
                        (*mappingReference)[idp].ranks.push_back(info.rank);
#endif
                    }
                }

                // Refinement and coarsening
                if (info.type == adaption::Type::TYPE_REFINEMENT || info.type == adaption::Type::TYPE_COARSENING) {
                    // Clear current mapper
                    for (long id : info.current) {
#if BITPIT_ENABLE_MPI
                        if (checkPart || info.rank == m_rank) {
#endif
                            mapping::Info &mappingInfo = (*mappingAdapted)[id];
                            mappingInfo.ids.clear();
                            mappingInfo.entity = mapping::Entity::ENTITY_CELL;
#if BITPIT_ENABLE_MPI
                            mappingInfo.ranks.clear();
                        } else {
                            m_partitionIR.map_rank_inverseMapping[info.rank][id].ids.clear();
                            m_partitionIR.map_rank_inverseMapping[info.rank][id].ranks.clear();
                            m_partitionIR.map_rank_inverseMapping[info.rank][id].entity = mapping::Entity::ENTITY_CELL;
                        }
#endif
                    }

                    for (long idprevious : info.previous) {
                        std::vector<long> *mappedIds;
                        std::vector<int> *mappedRanks;
#if BITPIT_ENABLE_MPI
                        if (checkPart || info.rank == m_rank) {
#endif
                            mappedIds   = &(m_previousMapping[idprevious].ids);
#if BITPIT_ENABLE_MPI
                            mappedRanks = &(m_previousMapping[idprevious].ranks);
                        } else {
                            mappedIds   = &(m_partitionIR.map_rank_previousMapping[info.rank][idprevious].ids);
                            mappedRanks = &(m_partitionIR.map_rank_previousMapping[info.rank][idprevious].ranks);
                        }
#endif
                        std::size_t icount = 0;
                        for (long mappedId : *mappedIds) {
                            VolOctree::OctantInfo oinfoprev = referencePatch->getCellOctant(mappedId);
                            const Octant *octant = referencePatch->getOctantPointer(oinfoprev);
                            uint64_t mortonreference = referencePatch->getTree().getMorton(octant);
                            uint64_t mortonlastdescreference = referencePatch->getTree().getLastDescMorton(octant);
                            uint8_t levelreference = referencePatch->getCellLevel(mappedId);
#if BITPIT_ENABLE_MPI
                            // IT'S ALWAYS = m_rank?!...
                            int mappedRank = (*mappedRanks)[icount];
#endif
                            for (long id : info.current) {
                                // Retrieve level of current cell
                                uint8_t level;
                                uint64_t morton, mortonlastdesc;
#if BITPIT_ENABLE_MPI
                                if (checkPart || info.rank == m_rank) {
#endif
                                    level = adaptedPatch->getCellLevel(id);
                                    VolOctree::OctantInfo octantIfo = adaptedPatch->getCellOctant(id);
                                    const Octant *octant = adaptedPatch->getOctantPointer(octantIfo);
                                    morton = adaptedPatch->getTree().getMorton(octant);
                                    mortonlastdesc = adaptedPatch->getTree().getLastDescMorton(octant);
#if BITPIT_ENABLE_MPI
                                } else {
                                    OctantIR *octantIR = m_partitionIR.map_rank_rec_octantIR[info.rank][id];
                                    level = octantIR->octant.getLevel();
                                    morton = adaptedPatch->getTree().getMorton(&octantIR->octant);
                                    mortonlastdesc = adaptedPatch->getTree().getLastDescMorton(&octantIR->octant);
                                }
#endif

                                bool checkmorton = false;

                                checkmorton |= (morton >= mortonreference && mortonlastdesc <= mortonlastdescreference);
                                checkmorton |= (morton <= mortonreference && mortonlastdesc >= mortonlastdescreference);

                                if (checkmorton) {

                                    if (level == levelreference) {
#if BITPIT_ENABLE_MPI
                                        if (checkPart || info.rank == m_rank) {
#endif
                                            mapping::Info &mappingInfo = (*mappingAdapted)[id];
                                            if (mappingInfo.ids.size() == 0) {
                                                mappingInfo.type = mapping::Type::TYPE_RENUMBERING;
                                                mappingInfo.ids.push_back(mappedId);
#if BITPIT_ENABLE_MPI
                                                mappingInfo.ranks.push_back(info.rank);
#endif
                                            }
#if BITPIT_ENABLE_MPI
                                        } else {
                                            if (m_partitionIR.map_rank_inverseMapping[info.rank][id].ids.size() == 0) {
                                                m_partitionIR.map_rank_inverseMapping[info.rank][id].type = mapping::Type::TYPE_RENUMBERING;
                                                m_partitionIR.map_rank_inverseMapping[info.rank][id].ids.push_back(mappedId);
                                                m_partitionIR.map_rank_inverseMapping[info.rank][id].ranks.push_back(mappedRank);
                                            }
                                        }
#endif
                                        (*mappingReference)[mappedId].type = mapping::Type::TYPE_RENUMBERING;
                                        (*mappingReference)[mappedId].ids.clear();
                                        (*mappingReference)[mappedId].ids.push_back(id);
#if BITPIT_ENABLE_MPI
                                        (*mappingReference)[mappedId].ranks.push_back(info.rank);
#endif
                                    } else if (level > levelreference) {
#if BITPIT_ENABLE_MPI
                                        if (checkPart || info.rank == m_rank) {
#endif
                                            mapping::Info &mappingInfo = (*mappingAdapted)[id];
                                            if (mappingInfo.ids.size() == 0) {
                                                mappingInfo.type = mapping::Type::TYPE_REFINEMENT;
                                                mappingInfo.ids.push_back(mappedId);
#if BITPIT_ENABLE_MPI
                                                mappingInfo.ranks.push_back(mappedRank);
#endif
                                            }
#if BITPIT_ENABLE_MPI
                                        } else {
                                            if (m_partitionIR.map_rank_inverseMapping[info.rank][id].ids.size() == 0) {
                                                m_partitionIR.map_rank_inverseMapping[info.rank][id].type = mapping::Type::TYPE_REFINEMENT;
                                                m_partitionIR.map_rank_inverseMapping[info.rank][id].ids.push_back(mappedId);
                                                m_partitionIR.map_rank_inverseMapping[info.rank][id].ranks.push_back(mappedRank);
                                            }
                                        }
#endif
                                        if (std::find((*mappingReference)[mappedId].ids.begin(), (*mappingReference)[mappedId].ids.end(), id) == (*mappingReference)[mappedId].ids.end()) {
                                            (*mappingReference)[mappedId].ids.push_back(id);
                                            (*mappingReference)[mappedId].type = mapping::Type::TYPE_COARSENING;
#if BITPIT_ENABLE_MPI
                                            (*mappingReference)[mappedId].ranks.push_back(info.rank);
#endif
                                        }
                                    } else if (level < levelreference) {
#if BITPIT_ENABLE_MPI
                                        if (checkPart || info.rank == m_rank) {
#endif
                                            mapping::Info &mappingInfo = (*mappingAdapted)[id];
                                            if (std::find(mappingInfo.ids.begin(), mappingInfo.ids.end(), mappedId) == mappingInfo.ids.end()) {
                                                mappingInfo.type = mapping::Type::TYPE_COARSENING;
                                                mappingInfo.ids.push_back(mappedId);
#if BITPIT_ENABLE_MPI
                                                mappingInfo.ranks.push_back(mappedRank);
#endif
                                            }
#if BITPIT_ENABLE_MPI
                                        } else {
                                            if (std::find(m_partitionIR.map_rank_inverseMapping[info.rank][id].ids.begin(), m_partitionIR.map_rank_inverseMapping[info.rank][id].ids.end(), mappedId) == m_partitionIR.map_rank_inverseMapping[info.rank][id].ids.end()) {
                                                m_partitionIR.map_rank_inverseMapping[info.rank][id].type = mapping::Type::TYPE_COARSENING;
                                                m_partitionIR.map_rank_inverseMapping[info.rank][id].ids.push_back(mappedId);
                                                m_partitionIR.map_rank_inverseMapping[info.rank][id].ranks.push_back(mappedRank);
                                            }
                                        }
#endif
                                        if ((*mappingReference)[mappedId].ids.size() == 0) {
                                            (*mappingReference)[mappedId].type = mapping::Type::TYPE_REFINEMENT;
                                            (*mappingReference)[mappedId].ids.push_back(id);
#if BITPIT_ENABLE_MPI
                                            (*mappingReference)[mappedId].ranks.push_back(info.rank);
#endif
                                        }
                                    } // End level < levelreference
                                } // End checkmorton
                            } // End id current

                            icount++;
                        } // Enf idmapped
                    } // End idprevious
                }
            } // End if deletion
        } // End info

#if BITPIT_ENABLE_MPI
        _communicateInverseMapperBack();
#endif

    } else {
        initialize(true);
    }
}

#if BITPIT_ENABLE_MPI
/**
 * Get the list of the octants of the partitions of the mapped mesh, different
 * from the local rank, overlapped to the local partition of the reference mesh.
 *
 * \result Returns the map with for each rank (key) the list of the octants
 * (argument) of the mapped mesh overlapped with the local partition of the
 * reference mesh.
 */
std::map<int, std::vector<long>> VolOctreeMapper::getReceivedMappedIds()
{
    std::map<int, std::vector<long>> received;

    for (OctantIR octir : m_partitionIR.list_rec_octantIR_before) {
        received[octir.rank].push_back(octir.id);
    }

    for (OctantIR octir : m_partitionIR.list_rec_octantIR_after) {
        received[octir.rank].push_back(octir.id);
    }

    return received;
}

/**
 * Get the list of the octants of the local partition of the reference mesh
 * overlapped to a different partition of the mapped mesh.
 *
 * \result Returns the map with for each rank (key) the list of the octants
 * (argument) of the reference mesh overlapped with the partition (rank) of
 * the mapped mesh.
 */
std::map<int, std::vector<long>> VolOctreeMapper::getSentReferenceIds()
{
    std::map<int, std::vector<long>> sent;

    // Recover id to be recv/send
    std::map<int, std::set<long>> rankIdSend;

    for (Cell &cell : m_referencePatch->getCells()) {
        long id = cell.getId();
        auto info = m_mapping[id];
        for (int rank : info.ranks) {
            if (rank != m_referencePatch->getRank()) {
                rankIdSend[rank].insert(id);
            }
        }
    }

    for (int rank = 0; rank < m_nProcs; rank++) {
        sent[rank].reserve(rankIdSend[rank].size());
    }

    for (std::map<int, std::set<long>>::iterator it=rankIdSend.begin(); it!=rankIdSend.end(); ++it) {
        for (long id : it->second) {
            sent[it->first].push_back(id);
        }
    }

    return sent;
}

/**
 * Get the list of the octants of the local partition of the mapped mesh
 * overlapped to a different partition of the reference mesh.
 *
 * \result Returns the map with for each rank (key) the list of the octants
 * (argument) of the mapped mesh overlapped with the partition (rank) of the
 * reference mesh.
 */
std::map<int, std::vector<long>> VolOctreeMapper::getSentMappedIds()
{
    std::map<int, std::vector<long>> sent;

    for (OctantIR octir : m_partitionIR.list_sent_octantIR) {
        sent[octir.rank].push_back(octir.id);
    }

    return sent;
}

#endif

/**
 * Map an input VolOctree mesh on a VolOctree reference mesh already set in
 * constructor.
 *
 * Requirement: the meshes have to be defined on the same identical domain.
 *
 * \param[in] fillInverse if set to true the inverse mapped (reference mesh
 * to input mesh) will be filled
 */
void VolOctreeMapper::_mapMeshes(bool fillInverse)
{
    std::array<double,3> originR = static_cast<VolOctree*>(m_referencePatch)->getOrigin();
    std::array<double,3> originM = static_cast<VolOctree*>(m_mappedPatch)->getOrigin();

    if (!(utils::DoubleFloatingEqual()(static_cast<VolOctree*>(m_referencePatch)->getLength(),static_cast<VolOctree*>(m_mappedPatch)->getLength()))
            || !(utils::DoubleFloatingEqual()(originR[0],originM[0]))
            || !(utils::DoubleFloatingEqual()(originR[1],originM[1]))
            || !(utils::DoubleFloatingEqual()(originR[2],originM[2]))) {
        throw std::runtime_error ("mesh mapper: different domain of VolOctree meshes not allowed.");
    }

#if BITPIT_ENABLE_MPI
    clearPartitionMapping();
#endif

    m_mapping.setDynamicKernel(&m_referencePatch->getCells(), PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED);
    if (fillInverse) {
        m_inverseMapping.setDynamicKernel(&m_mappedPatch->getCells(), PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED);
    } else {
        clearInverseMapping();
    }

#if BITPIT_ENABLE_MPI
    if (!(m_referencePatch->isPartitioned())) {
#endif
        long indRef = 0;
        _mapMeshesSamePartition(nullptr, nullptr, fillInverse, &indRef);

#if BITPIT_ENABLE_MPI
    } else {
        bool checkPart = checkPartition();
        if (checkPart) {
            long indRef = 0;
            _mapMeshesSamePartition(nullptr, nullptr, fillInverse, &indRef);
        } else {
            _mapMeshPartitioned(fillInverse);
        }
    }
#endif
}

/**
 * Map an input list of octants on a list of octants of a reference mesh.
 *
 * Requirements: the meshes have to be defined on the same identical domain
 * and the lists are considered as identically parallel partitioned (i.e.,
 * the first and the last descendant octants of the input list are contained
 * in the covered domain of the reference list).
 *
 * \param[in] octantsIRReference is the list of reference octants, if set to
 * null the whole local mesh is used
 * \param[in] octantsIRMapped is the list of mapped octants, is set to null
 * the whole local mesh is used
 * \param[in] fillInverse if set to true the inverse mapped (reference mesh
 * to input mesh) will be filled
 * \param[in,out] indRef is the dtarting octant (ending octant in output) of
 * reference mesh to compute the mapping (default = 0).
 */
void VolOctreeMapper::_mapMeshesSamePartition(const std::vector<OctantIR> *octantsIRReference, const std::vector<OctantIR> *octantsIRMapped,
                                              bool fillInverse, long *indRef)
{
    bitpit::VolOctree *referencePatch = static_cast<VolOctree*>(m_referencePatch);
    bitpit::VolOctree *mappedPatch    = static_cast<VolOctree*>(m_mappedPatch);

    // Fill IR with meshes if list pointer is null
#if BITPIT_ENABLE_MPI
    bool localMapped = false;
#endif

    std::vector<OctantIR> tempOctantsIRReference;
    if (!octantsIRReference) {
        long n = referencePatch->getInternalCount();
        tempOctantsIRReference.reserve(n);
        for (long i = 0; i < n; i++) {
            VolOctree::OctantInfo octantIfoRef(i, true);
            const Octant *octRef = static_cast<VolOctree*>(m_referencePatch)->getOctantPointer(octantIfoRef);
            long idRef = referencePatch->getOctantId(octantIfoRef);
#if BITPIT_ENABLE_MPI
            tempOctantsIRReference.emplace_back(*octRef, idRef, idRef, m_rank);
#else
            tempOctantsIRReference.emplace_back(*octRef, idRef, idRef);
#endif
        }
        octantsIRReference = &tempOctantsIRReference;
    }

    std::vector<OctantIR> tempOctantsIRMapped;
    if (!octantsIRMapped) {
        long n = mappedPatch->getInternalCount();
        tempOctantsIRMapped.reserve(n);
        for (long i = 0; i < n; i++) {
            VolOctree::OctantInfo octantIfoMap(i, true);
            const Octant *octMap = mappedPatch->getOctantPointer(octantIfoMap);
            long idMap = mappedPatch->getOctantId(octantIfoMap);
#if BITPIT_ENABLE_MPI
            tempOctantsIRMapped.emplace_back(*octMap, idMap, idMap, m_rank);
#else
            tempOctantsIRMapped.emplace_back(*octMap, idMap, idMap);
#endif
        }
        octantsIRMapped = &tempOctantsIRMapped;
#if BITPIT_ENABLE_MPI
        localMapped = true;
#endif
    }

    // Define a map for inverse mapper structure
    //
    // If serial the map is transfer to inverse mapper directly.
    //
    // If mapped mesh is parallel and partitioned the map object is
    // communicated at the end of the procedure.
    long nRef  = octantsIRReference->size();
    long nMap  = octantsIRMapped->size();

    long indMap = 0;
    if (*indRef < nRef && nMap > 0) {
        // Place indMap at first last desc morton greater than first morton of
        // reference octants
        uint64_t morton = referencePatch->getTree().getMorton(&octantsIRReference->at(*indRef).octant);

        uint64_t mortonMapLastdesc;
        for (indMap = 0; indMap < nMap; ++indMap) {
            mortonMapLastdesc = mappedPatch->getTree().getLastDescMorton(&octantsIRMapped->at(indMap).octant);
            if (mortonMapLastdesc >= morton) {
                break;
            }
        }
    }

    std::unordered_map<long, mapping::Info> inverseGlobalMapping;
    while (*indRef < nRef && indMap < nMap) {
        long idRef = octantsIRReference->at(*indRef).id;
        long idMap = octantsIRMapped->at(indMap).id;
        long gidMap = octantsIRMapped->at(indMap).globalId;
#if BITPIT_ENABLE_MPI
        int rank = octantsIRMapped->at(indMap).rank;
#endif

        m_mapping[idRef].entity = mapping::Entity::ENTITY_CELL;
        if (fillInverse) {
            inverseGlobalMapping[gidMap].entity = mapping::Entity::ENTITY_CELL;
        }

        if (octantsIRMapped->at(indMap).octant.getLevel() == octantsIRReference->at(*indRef).octant.getLevel()) {
            m_mapping[idRef].ids.push_back(idMap);
            m_mapping[idRef].type = mapping::Type::TYPE_RENUMBERING;
#if BITPIT_ENABLE_MPI
            m_mapping[idRef].ranks.push_back(rank);
#endif
            if (fillInverse) {
                inverseGlobalMapping[gidMap].ids.push_back(idRef);
                inverseGlobalMapping[gidMap].type = mapping::Type::TYPE_RENUMBERING;
#if BITPIT_ENABLE_MPI
                inverseGlobalMapping[gidMap].ranks.push_back(m_rank);
#endif
            }
            (*indRef)++;
            indMap++;
        }
        else if (octantsIRMapped->at(indMap).octant.getLevel() > octantsIRReference->at(*indRef).octant.getLevel()) {
            m_mapping[idRef].type = mapping::Type::TYPE_COARSENING;

            uint64_t mortonlastdesc = referencePatch->getTree().getLastDescMorton(&octantsIRReference->at(*indRef).octant);
            uint64_t mortonMap = mappedPatch->getTree().getMorton(&octantsIRMapped->at(indMap).octant);
            uint64_t mortonMapLastDesc = mappedPatch->getTree().getLastDescMorton(&octantsIRMapped->at(indMap).octant);
            while(mortonMap <= mortonlastdesc) {
                m_mapping[idRef].ids.push_back(idMap);
#if BITPIT_ENABLE_MPI
                m_mapping[idRef].ranks.push_back(rank);
#endif
                if (fillInverse) {
                    inverseGlobalMapping[gidMap].type = mapping::Type::TYPE_REFINEMENT;
                    inverseGlobalMapping[gidMap].ids.push_back(idRef);
#if BITPIT_ENABLE_MPI
                    inverseGlobalMapping[gidMap].ranks.push_back(m_rank);
#endif
                }
                indMap++;
                if (indMap == nMap) {
                    break;
                }
                idMap = octantsIRMapped->at(indMap).id;
                gidMap = octantsIRMapped->at(indMap).globalId;
#if BITPIT_ENABLE_MPI
                rank = octantsIRMapped->at(indMap).rank;
#endif
                mortonMap = mappedPatch->getTree().getMorton(&octantsIRMapped->at(indMap).octant);
                mortonMapLastDesc = mappedPatch->getTree().getLastDescMorton(&octantsIRMapped->at(indMap).octant);
            }
            if (mortonMapLastDesc >= mortonlastdesc) {
                (*indRef)++;
            }
        }
        else if (octantsIRMapped->at(indMap).octant.getLevel() < octantsIRReference->at(*indRef).octant.getLevel()) {

            uint64_t morton = referencePatch->getTree().getMorton(&octantsIRReference->at(*indRef).octant);
            uint64_t mortonlastdescmesh = mappedPatch->getTree().getLastDescMorton(&octantsIRMapped->at(indMap).octant);

            if (fillInverse) {
                inverseGlobalMapping[gidMap].type = mapping::Type::TYPE_COARSENING;
            }

            while (morton <= mortonlastdescmesh) {
                m_mapping[idRef].type = mapping::Type::TYPE_REFINEMENT;
                m_mapping[idRef].ids.push_back(idMap);
#if BITPIT_ENABLE_MPI
                m_mapping[idRef].ranks.push_back(rank);
#endif
                if (fillInverse) {
                    inverseGlobalMapping[gidMap].ids.push_back(idRef);
#if BITPIT_ENABLE_MPI
                    inverseGlobalMapping[gidMap].ranks.push_back(m_rank);
#endif
                }
                (*indRef)++;
                if (*indRef == nRef) {
                    break;
                }
                morton = referencePatch->getTree().getMorton(&octantsIRReference->at(*indRef).octant);
                idRef = octantsIRReference->at(*indRef).id;
            }
            indMap++;
        }
    }

    // Fill inverse mapping
    if (fillInverse) {
        std::unordered_map<long, mapping::Info> inverseLocalMapping;
#if BITPIT_ENABLE_MPI
        if (mappedPatch->isPartitioned() && !localMapped && !checkPartition()) {
            _communicateInverseMapper(octantsIRMapped, inverseGlobalMapping, &inverseLocalMapping);
        } else {
            inverseGlobalMapping.swap(inverseLocalMapping);
        }
#else
        inverseGlobalMapping.swap(inverseLocalMapping);
#endif

        for (const auto &mappingEntry : inverseLocalMapping) {
            m_inverseMapping[mappingEntry.first].type = mappingEntry.second.type;
            for (long inverseMappedId : mappingEntry.second.ids){
                m_inverseMapping[mappingEntry.first].ids.push_back(inverseMappedId);
            }
#if BITPIT_ENABLE_MPI
            for (int inverseMappedRank : mappingEntry.second.ranks){
                m_inverseMapping[mappingEntry.first].ranks.push_back(inverseMappedRank);
            }
#endif
        }
    }
}

#if BITPIT_ENABLE_MPI
/**
 * Check if the reference and the mapped meshes have the same geometric
 * partitioning.
 *
 * \result Returns true if the meshes have identical geometric partitioning,
 * false otherwise.
 */
bool VolOctreeMapper::checkPartition()
{
    std::vector<uint64_t> partitionMapped    = static_cast<VolOctree*>(m_mappedPatch)->getTree().getPartitionLastDesc();
    std::vector<uint64_t> partitionReference = static_cast<VolOctree*>(m_referencePatch)->getTree().getPartitionLastDesc();
    for (int rank = 0; rank < m_nProcs; ++rank) {
        if (partitionReference[rank] != partitionMapped[rank]) {
            return false;
        }
    }

    return true;
}

/**
 * Mapping computing between two partitioned meshes (stored in class members).
 *
 * \param[in] fillInverse if set to true the inverse mapped (reference mesh
 * to input mesh) will be filled
 */
void VolOctreeMapper::_mapMeshPartitioned(bool fillInverse)
{
    _recoverPartition();

    // Fill IR with reference mesh
    //
    // TODO: make a method to do that
    long n = static_cast<VolOctree*>(m_referencePatch)->getInternalCount();
    std::vector<OctantIR> octantsIRReference;
    octantsIRReference.reserve(n);
    for (long i = 0; i < n; i++) {
        VolOctree::OctantInfo octantIfoRef(i, true);
        const Octant *octRef = static_cast<VolOctree*>(m_referencePatch)->getOctantPointer(octantIfoRef);
        long idRef = static_cast<VolOctree*>(m_referencePatch)->getOctantId(octantIfoRef);
        octantsIRReference.emplace_back(*octRef, idRef, idRef, m_rank);
    }

    long indRef = 0;
    _mapMeshesSamePartition(&octantsIRReference, &m_partitionIR.list_rec_octantIR_before, fillInverse, &indRef);
    _mapMeshesSamePartition(&octantsIRReference, nullptr, fillInverse, &indRef);
    _mapMeshesSamePartition(&octantsIRReference, &m_partitionIR.list_rec_octantIR_after, fillInverse, &indRef);
}

/**
 * Recover the overlapped partitions between the mapped meshes.
 *
 * \result Return true if the partition structures are different from the ones
 * stored in class (e.g. after an adaptation or after a load balancing).
 */
bool VolOctreeMapper::_recoverPartition()
{
    // Now the mapped mesh is repartitioned (is copied...) over the processes
    // to match the partitioning between the two meshes

    // Find owner of reference partition over mapped mesh partition
    VolOctree* referencePatch = static_cast<VolOctree*>(m_referencePatch);
    VolOctree* mappedPatch = static_cast<VolOctree*>(m_mappedPatch);
    std::vector<uint64_t> partitionLDReference = referencePatch->getTree().getPartitionLastDesc();
    std::vector<uint64_t> partitionLDMapped = mappedPatch->getTree().getPartitionLastDesc();
    std::vector<uint64_t> partitionFDReference = referencePatch->getTree().getPartitionFirstDesc();
    std::vector<uint64_t> partitionFDMapped = mappedPatch->getTree().getPartitionFirstDesc();

    if (m_partitionIR.partitionFDReference.size() != 0) {
        if (m_partitionIR.partitionFDMapped != partitionFDMapped ||
                m_partitionIR.partitionLDMapped != partitionLDMapped ||
                m_partitionIR.partitionFDReference != partitionFDReference ||
                m_partitionIR.partitionLDReference != partitionLDReference) {
            m_partitionIR.partitionFDReference.clear();
            m_partitionIR.partitionLDReference.clear();
            m_partitionIR.partitionFDMapped.clear();
            m_partitionIR.partitionLDMapped.clear();
            return true;
        }
    }

    std::map<int, std::vector<int>> frommapped_rank;
    std::map<int,std::vector<int>> toreference_rank;
    for (int ref_rank = 0; ref_rank < m_nProcs; ref_rank++) {
        uint64_t local_first_morton = partitionFDReference[ref_rank];
        uint64_t local_last_morton = partitionLDReference[ref_rank];

        for (int irank = 0; irank < m_nProcs; irank++) {
            if (irank != ref_rank) {
                if (partitionFDMapped[irank] <= local_first_morton && partitionLDMapped[irank] >= local_first_morton) {
                    frommapped_rank[ref_rank].push_back(irank);
                    toreference_rank[irank].push_back(ref_rank);
                }
                else if (partitionFDMapped[irank] <= local_last_morton && partitionLDMapped[irank] >= local_last_morton) {
                    frommapped_rank[ref_rank].push_back(irank);
                    toreference_rank[irank].push_back(ref_rank);
                }
                else if (partitionFDMapped[irank] <= local_first_morton && partitionLDMapped[irank] >= local_last_morton) {
                    frommapped_rank[ref_rank].push_back(irank);
                    toreference_rank[irank].push_back(ref_rank);
                }
                else if (partitionFDMapped[irank] >= local_first_morton && partitionLDMapped[irank] <= local_last_morton) {
                    frommapped_rank[ref_rank].push_back(irank);
                    toreference_rank[irank].push_back(ref_rank);
                }
            }
        }
    }

    clearPartitionMappingLists();

    // Locally the mapped mesh build the lists of octants to send
    std::map<int, std::vector<Octant>> list_octant;
    std::map<int, std::vector<long>> list_id;
    std::map<int, std::vector<long>> list_globalId;
    uint32_t idx = 0;
    uint64_t morton = mappedPatch->getTree().getLastDescMorton(idx);
    for (int reference_rank : toreference_rank[m_rank]) {
        uint64_t reference_first_morton = partitionFDReference[reference_rank];
        uint64_t reference_last_morton = partitionLDReference[reference_rank];
        while (morton < reference_first_morton) {
            idx++;
            if (idx == mappedPatch->getTree().getNumOctants()) {
                break;
            }
            morton = mappedPatch->getTree().getLastDescMorton(idx);
        }
        while (morton < reference_last_morton) {
            Octant oct = *mappedPatch->getTree().getOctant(idx);
            list_octant[reference_rank].push_back(oct);
            VolOctree::OctantInfo octantIfo(idx, true);
            long id = mappedPatch->getOctantId(octantIfo);
            list_id[reference_rank].push_back(id);
            long globalId = mappedPatch->getTree().getGlobalIdx(idx);
            list_globalId[reference_rank].push_back(globalId);
            m_partitionIR.list_sent_octantIR.emplace_back(oct, id, globalId, reference_rank);
            idx++;
            if (idx == mappedPatch->getTree().getNumOctants()) {
                break;
            }
            morton = mappedPatch->getTree().getMorton(idx);
        }
    }

    // Communications
    //
    // Send local mapped octants to reference rank and receive mapped octants
    // from other processes
    // Note. The communicated octants have to be ordered by Morton index.
    // For this reason the communications must be performed ordered by rank.
    std::size_t octantBinarySize = Octant::getBinarySize() + 2 * sizeof(long);

    // Build send buffers
    DataCommunicator octCommunicator(m_communicator);

    // Set size
    //
    // TODO: make accessible global variables in ParaTree
    for (int reference_rank : toreference_rank[m_rank]) {
        const std::vector<Octant> &rankOctants = list_octant[reference_rank];
        std::size_t nRankOctants = rankOctants.size();

        // Set buffer size
        std::size_t buffSize = sizeof(std::size_t) + nRankOctants * octantBinarySize;
        octCommunicator.setSend(reference_rank,buffSize);

        // Fill buffer with octants
        SendBuffer &sendBuffer = octCommunicator.getSendBuffer(reference_rank);
        sendBuffer << nRankOctants;
        for (std::size_t n = 0; n < nRankOctants; ++n) {
            const Octant &octant = rankOctants[n];
            sendBuffer << octant;

            long id = list_id[reference_rank][n];
            sendBuffer << id;

            long globalId = list_globalId[reference_rank][n];
            sendBuffer << globalId;
        }
    }

    octCommunicator.discoverRecvs();
    octCommunicator.startAllRecvs();
    octCommunicator.startAllSends();

    std::vector<int> recvRanks = octCommunicator.getRecvRanks();
    std::sort(recvRanks.begin(),recvRanks.end());

    std::vector<OctantIR> &list_rec_octantIR_before = m_partitionIR.list_rec_octantIR_before;
    std::vector<OctantIR> &list_rec_octantIR_after = m_partitionIR.list_rec_octantIR_after;
    std::vector<OctantIR> *_list_rec_octantIR;

    list_rec_octantIR_before.clear();
    list_rec_octantIR_after.clear();

    for (int rank : recvRanks) {
        octCommunicator.waitRecv(rank);

        if (rank < m_rank) {
            _list_rec_octantIR = &list_rec_octantIR_before;
        } else if (rank > m_rank) {
            _list_rec_octantIR = &list_rec_octantIR_after;
        } else {
            break;
        }

        RecvBuffer &recvBuffer = octCommunicator.getRecvBuffer(rank);

        std::size_t nRecvOctants;
        recvBuffer >> nRecvOctants;
        for (std::size_t n = 0; n < nRecvOctants; ++n) {
            Octant octant;
            recvBuffer >> octant;

            long id;
            recvBuffer >> id;

            long globalId;
            recvBuffer >> globalId;

            _list_rec_octantIR->emplace_back(octant, id, globalId, rank);
        }

        for (OctantIR &octantIR : *_list_rec_octantIR) {
            long id = octantIR.id;
            m_partitionIR.map_rank_rec_octantIR[rank][id] = &octantIR;
        }
    }

    octCommunicator.waitAllSends();

    m_partitionIR.partitionLDReference = partitionLDReference;
    m_partitionIR.partitionLDMapped = partitionLDMapped;
    m_partitionIR.partitionFDReference = partitionFDReference;
    m_partitionIR.partitionFDMapped = partitionFDMapped;

    return false;
}

/**
 * Communicate inverse mapping info of overlapped partitions between processes.
 *
 * \param[in] octantsIRMapped is the lisst of overlapped octants of mapped mesh
 * \param[in] inverseGlobalMapping is the inverse global mapping
 * \param[out] inverseLocalMapping on output contains the inverse mapping
 */
void VolOctreeMapper::_communicateInverseMapper(const std::vector<OctantIR> *octantsIRMapped,
                                                const std::unordered_map<long, mapping::Info> &inverseGlobalMapping,
                                                std::unordered_map<long, mapping::Info> *inverseLocalMapping)
{
    // Recover mapping elements to send (to partitions of mapped mesh)
    std::set<int> toRanks;
    std::map<int, std::vector<long>> toRankGlobalId;
    std::map<long, long > globalIdToId;
    for (const OctantIR &octantIR : *octantsIRMapped) {
        if (!inverseGlobalMapping.count(octantIR.globalId)){
            continue;
        }
        toRankGlobalId[octantIR.rank].push_back(octantIR.globalId);
        globalIdToId[octantIR.globalId] = octantIR.id;
        toRanks.insert(octantIR.rank);
    }

    // Build send buffers
    DataCommunicator mapCommunicator(m_communicator);

    // Set size
    for (int rank : toRanks) {
        // Get buffer size
        std::size_t buffSize = 0;
        for (long globalId : toRankGlobalId[rank]) {
            const mapping::Info &info = inverseGlobalMapping.at(globalId);
            int mappedSize = info.ids.size();
            std::size_t infoBytes = std::size_t(sizeof(long) + sizeof(int) + sizeof(int) + sizeof(int) + (sizeof(long))*mappedSize);
            buffSize += infoBytes;
        }
        buffSize += std::size_t(sizeof(int));
        mapCommunicator.setSend(rank, buffSize);

        // Fill buffer with octants and local map for inverse mapping
        SendBuffer &sendBuffer = mapCommunicator.getSendBuffer(rank);
        sendBuffer << int(toRankGlobalId[rank].size());
        for (long globalId : toRankGlobalId[rank]) {
            const mapping::Info &info = inverseGlobalMapping.at(globalId);
            sendBuffer << globalId;
            sendBuffer << int(info.type);
            sendBuffer << int(info.entity);
            int mappedSize = info.ids.size();
            sendBuffer << mappedSize;
            for (long refId : info.ids) {
                sendBuffer << refId;
            }

            m_partitionIR.map_rank_inverseMapping[rank][globalIdToId[globalId]] = info;

        }
    }

    mapCommunicator.discoverRecvs();
    mapCommunicator.startAllRecvs();
    mapCommunicator.startAllSends();

    int nCompletedRecvs = 0;
    while (nCompletedRecvs < mapCommunicator.getRecvCount()) {
        int rank = mapCommunicator.waitAnyRecv();
        RecvBuffer &recvBuffer = mapCommunicator.getRecvBuffer(rank);

        int nof;
        recvBuffer >> nof;
        for (int i = 0; i < nof; i++) {
            long globalId;
            recvBuffer >> globalId;

            uint32_t idx = static_cast<VolOctree*>(m_mappedPatch)->getTree().getLocalIdx(globalId);
            VolOctree::OctantInfo octantIfo(idx, true);

            int type;
            recvBuffer >> type;

            int entity;
            recvBuffer >> entity;

            long id = static_cast<VolOctree*>(m_mappedPatch)->getOctantId(octantIfo);
            mapping::Info &info = (*inverseLocalMapping)[id];
            info.type = mapping::Type(type);
            info.entity = mapping::Entity(entity);
            int nmap;
            recvBuffer >> nmap;
            for (int j=0; j<nmap; j++) {
                long idref;
                recvBuffer >> idref;
                info.ids.push_back(idref);
                info.ranks.push_back(rank);
            }
        }

        ++nCompletedRecvs;
    }

    mapCommunicator.waitAllSends();
}

/**
 * Communicate inverse mapping info of overlapped partitions back to the
 * processes of mapped partitions.
 */
void VolOctreeMapper::_communicateInverseMapperBack()
{
    // Recover mapping elements to send (to partitions of mapped mesh)
    std::set<int> toRanks;
    std::map<int, std::vector<long>> toRankId;
    for (const auto &mappingEntry : m_partitionIR.map_rank_inverseMapping) {
        toRanks.insert(mappingEntry.first);

        const std::unordered_map<long, mapping::Info> &inverseMappingInfo = m_partitionIR.map_rank_inverseMapping[mappingEntry.first];
        for (const auto &mappingEntry : inverseMappingInfo) {
            toRankId[mappingEntry.first].push_back(mappingEntry.first);
        }
    }

    // Build send buffers
    DataCommunicator mapCommunicator(m_communicator);

    // Set size
    for (int rank : toRanks) {
        const std::unordered_map<long, mapping::Info> &inverseMappingInfo = m_partitionIR.map_rank_inverseMapping[rank];

        // Get buffer size
        std::size_t buffSize = 0;
        for (long id : toRankId[rank]) {
            const mapping::Info &info = inverseMappingInfo.at(id);
            int mappedSize = info.ids.size();
            std::size_t infoBytes = std::size_t(sizeof(long) + sizeof(int) + sizeof(int) + sizeof(int) + (sizeof(long))*mappedSize);
            buffSize += infoBytes;
        }
        buffSize += std::size_t(sizeof(int));
        mapCommunicator.setSend(rank,buffSize);

        // Fill buffer with octants and local map for inverse mapping
        SendBuffer &sendBuffer = mapCommunicator.getSendBuffer(rank);
        sendBuffer << int(toRankId[rank].size());
        for (long id : toRankId[rank]) {
            const mapping::Info &info = inverseMappingInfo.at(id);
            sendBuffer << id;
            sendBuffer << int(info.type);
            sendBuffer << int(info.entity);
            int mappedSize = info.ids.size();
            sendBuffer << mappedSize;
            for (long refId : info.ids) {
                sendBuffer << refId;
            }
        }
    }

    mapCommunicator.discoverRecvs();
    mapCommunicator.startAllRecvs();
    mapCommunicator.startAllSends();

    int nCompletedRecvs = 0;
    while (nCompletedRecvs < mapCommunicator.getRecvCount()) {
        int rank = mapCommunicator.waitAnyRecv();
        RecvBuffer &recvBuffer = mapCommunicator.getRecvBuffer(rank);

        int nof;
        recvBuffer >> nof;
        for (int i = 0; i < nof; i++) {
            long id;
            recvBuffer >> id;

            int type;
            recvBuffer >> type;

            int entity;
            recvBuffer >> entity;

            mapping::Info &info = m_inverseMapping[id];
            info.type = mapping::Type(type);
            info.entity = mapping::Entity(entity);
            int nmap;
            recvBuffer >> nmap;
            for (int j=0; j<nmap; j++) {
                long idref;
                recvBuffer >> idref;
                info.ids.push_back(idref);
                info.ranks.push_back(rank);
            }
        }

        ++nCompletedRecvs;
    }

    mapCommunicator.waitAllSends();

}

/**
 * Communicate adaption info of overlapped partitions of the mapped mesh to
 * the processes of reference partitions.
 *
 * \param[in] adaptionInfo are the adaptation info that describe the changes
 * of the local mapped partitions
 * \param[out] adaptionInfoRef are the adaptation info that describe the
 * changes of the local reference partitions
 */
void VolOctreeMapper::_communicateMappedAdaptionInfo(const std::vector<adaption::Info> &adaptionInfoMap, std::vector<adaption::Info> *adaptionInfoRef)
{
    // Recover mapping elements to send (to partitions of mapped mesh)
    std::set<int> toRanks;
    std::map<int, std::set<long>> toRankInd;
    std::map<int, std::vector<long>> toRankId;

    std::size_t n = 0;
    for (auto &info : adaptionInfoMap) {
        for (long id : info.previous) {
            auto itPreviousMapping = m_previousMapping.find(id);
            if (itPreviousMapping == m_previousMapping.end()){
                continue;
            }
            for (int rank : itPreviousMapping->second.ranks) {
                toRanks.insert(rank);
                toRankInd[rank].insert(n);
                toRankId[rank].push_back(id);
            }
        }
        ++n;
    }

    // Build send buffers
    DataCommunicator mapCommunicator(m_communicator);

    // Initialize communications
    for (int rank : toRanks) {
        std::size_t buffSize = 0;
        for (long id : toRankId[rank]) {
            const mapping::Info &info = m_previousMapping.at(id);
            int mappedSize = info.ids.size();
            std::size_t infoBytes = std::size_t(sizeof(long) + 3*sizeof(int) + (sizeof(long))*mappedSize);
            buffSize += infoBytes;
        }
        buffSize += std::size_t(sizeof(int));

        for (long ind : toRankInd[rank]) {
            const adaption::Info &info = adaptionInfoMap[ind];
            int currentSize = info.current.size();
            int previousSize = info.previous.size();
            std::size_t infoBytes = std::size_t(5*sizeof(int) + (sizeof(long))*(currentSize+previousSize));
            buffSize += infoBytes;
        }
        buffSize += std::size_t(sizeof(int));

        mapCommunicator.setSend(rank,buffSize);

        // Fill buffer with octants and local map for inverse mapper
        SendBuffer &sendBuffer = mapCommunicator.getSendBuffer(rank);
        sendBuffer << int(toRankId[rank].size());
        for (long id : toRankId[rank]) {
            auto itPreviousMapping = m_previousMapping.find(id);
            if (itPreviousMapping == m_previousMapping.end()){
                continue;
            }
            const mapping::Info &info = itPreviousMapping->second;
            sendBuffer << id;
            sendBuffer << int(info.type);
            sendBuffer << int(info.entity);
            int mappedSize = info.ids.size();
            sendBuffer << mappedSize;
            for (long refId : info.ids) {
                sendBuffer << refId;
            }
        }

        sendBuffer << int(toRankInd[rank].size());
        for (long ind : toRankInd[rank]) {
            const adaption::Info &info = adaptionInfoMap[ind];
            sendBuffer << int(info.type);
            sendBuffer << int(info.entity);
            int currentSize = info.current.size();
            int previousSize = info.previous.size();
            sendBuffer << currentSize;
            for (long Id : info.current) {
                sendBuffer << Id;
            }
            sendBuffer << previousSize;
            for (long Id : info.previous) {
                sendBuffer << Id;
            }

            // Rank of adaption is used to identify the rank where the adaption
            // occurs
            sendBuffer << m_rank;
        }

    }

    // Execute communications
    mapCommunicator.discoverRecvs();
    mapCommunicator.startAllRecvs();
    mapCommunicator.startAllSends();

    // Retrieve data
    int nCompletedRecvs = 0;
    while (nCompletedRecvs < mapCommunicator.getRecvCount()) {
        int rank = mapCommunicator.waitAnyRecv();

        RecvBuffer &recvBuffer = mapCommunicator.getRecvBuffer(rank);
        int nofMapper;
        recvBuffer >> nofMapper;
        for (int i = 0; i < nofMapper; i++) {
            long id;
            recvBuffer >> id;

            int mtype;
            recvBuffer >> mtype;

            int mentity;
            recvBuffer >> mentity;

            mapping::Info info;
            info.type = mapping::Type(mtype);
            info.entity = mapping::Entity(mentity);
            int nmap;
            recvBuffer >> nmap;
            for (int j=0; j<nmap; j++) {
                long idref;
                recvBuffer >> idref;
                info.ids.push_back(idref);
                info.ranks.push_back(m_rank);
            }
            m_partitionIR.map_rank_previousMapping[rank][id] = info;
            m_partitionIR.map_rank_inverseMapping[rank].erase(id);
        }

        int nofInfo;
        recvBuffer >> nofInfo;
        for (int i = 0; i < nofInfo; i++) {
            int type;
            recvBuffer >> type;

            int entity;
            recvBuffer >> entity;

            adaption::Info info;
            info.type = adaption::Type(type);
            info.entity = adaption::Entity(entity);
            int ncurrent;
            recvBuffer >> ncurrent;
            for (int j=0; j<ncurrent; j++) {
                long id;
                recvBuffer >> id;
                info.current.push_back(id);
            }
            int nprevious;
            recvBuffer >> nprevious;
            for (int j=0; j<nprevious; j++) {
                long id;
                recvBuffer >> id;
                info.previous.push_back(id);
            }
            int arank;
            recvBuffer >> arank;
            info.rank = arank;
            adaptionInfoRef->push_back(info);
        }

        ++nCompletedRecvs;
    }

    mapCommunicator.waitAllSends();
}

#endif

}
