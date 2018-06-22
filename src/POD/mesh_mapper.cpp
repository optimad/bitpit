/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

#include <cassert>

#include "mesh_mapper.hpp"

namespace bitpit {

/**
 * \class MeshMapper
 * \ingroup Ref
 *
 * \brief The MeshMapper is the class to map two meshes.
 *
 * The MeshMapper allows to map meshes of class: VolOctree.
 *
 */

/**
 * Creates a new MeshMapper object.
 */
# if BITPIT_ENABLE_MPI
/**
 * \param[in] comm The MPI communicator used by the Ref object. MPI_COMM_WORLD is the default value.
 */
MeshMapper::MeshMapper(MPI_Comm comm)
# else
MeshMapper::MeshMapper()
# endif
: m_mapper(1)
{

#if BITPIT_ENABLE_MPI
    m_communicator = MPI_COMM_NULL;
    initializeCommunicator(comm);
    MPI_Comm_size(m_communicator, &m_nProcs);
    MPI_Comm_rank(m_communicator, &m_rank);
#else
    m_rank = 0;
    m_nProcs = 1;
#endif

}

/**
 * Destructor of MeshMapper
 */
MeshMapper::~MeshMapper()
{
}

/**
 * Clear mapping members
 */
void MeshMapper::clear()
{
    clearMapping();
    clearInverseMapping();
}

/**
 * Clear direct mapping
 */
void MeshMapper::clearMapping()
{
    if (m_mapper.getKernel() != nullptr)
        m_mapper.unsetKernel(true);
}

/**
 * Clear inverse mapping
 */
void MeshMapper::clearInverseMapping()
{
    if (m_invmapper.getKernel() != nullptr)
        m_invmapper.unsetKernel(false);
}

/**
 * Get direct mapping
 */
const bitpit::PiercedStorage<mapping::Info> & MeshMapper::getMapping()
{
    return m_mapper;
}

/**
 * Get inverse mapping
 */
const bitpit::PiercedStorage<mapping::Info> & MeshMapper::getInverseMapping()
{
    return m_invmapper;
}

/**
 * Map an input mesh on a reference mesh. The specialization of the meshes are automatically recovered
 * by the method. VolumeKernel is the Base class from which the two meshes have to be derived.
 * Allowed mesh classes:
 * 1. VolOctree.
 *
 * \param[in] meshReference Pointer to VolumeKernel reference mesh
 * \param[in] meshMapped Pointer to VolumeKernel input mesh to map
 * \param[in] fillInv If true even the inverse mapping (reference mesh to input mesh) is filled.
 */
void MeshMapper::mapMeshes(bitpit::VolumeKernel * meshReference, bitpit::VolumeKernel * meshMapped, bool fillInv)
{
    clear();

    {
        VolOctree* _meshReference = dynamic_cast<VolOctree*>(meshReference);
        VolOctree* _meshMapped = dynamic_cast<VolOctree*>(meshMapped);
        if (_meshReference && _meshMapped){
            m_referenceMesh = meshReference;
            m_mappedMesh = meshMapped;
            _mapMeshes(_meshReference, _meshMapped, fillInv);
        }
    }
}


void MeshMapper::mappingAdaptionPreparare(const std::vector<adaption::Info> & infoAdapt, bool reference)
{
    m_previousmapper.clear();
    PiercedStorage<adaption::Info>* pmapper;

    if (reference)
        pmapper = &m_mapper;
    else
        pmapper = &m_invmapper;

    for (const adaption::Info & info : infoAdapt){
        for (const long & id : info.previous)
            m_previousmapper[id] = pmapper->at(id);
    }
}

void MeshMapper::mappingAdaptionUpdate(const std::vector<adaption::Info> & infoAdapt, bool reference, bool fillInv)
{

    VolOctree* meshAdapted;
    PiercedStorage<adaption::Info>* mapperAdapted;
    VolOctree* meshMapped;
    PiercedStorage<adaption::Info>* mapperMapped;
    if (reference){
        meshAdapted = dynamic_cast<VolOctree*>(m_referenceMesh);
        mapperAdapted = &m_mapper;
        meshMapped = dynamic_cast<VolOctree*>(m_mappedMesh);
        mapperMapped = &m_invmapper;
    }
    else{
        meshMapped = dynamic_cast<VolOctree*>(m_referenceMesh);
        mapperAdapted = &m_invmapper;
        meshAdapted = dynamic_cast<VolOctree*>(m_mappedMesh);
        mapperMapped = &m_mapper;
    }

    assert(meshAdapted != nullptr);
    assert(meshMapped != nullptr);

    for (const adaption::Info & info : infoAdapt){
        if (info.type == adaption::Type::TYPE_RENUMBERING){
            long id = info.current[0];
            (*mapperAdapted)[id].current.clear();
            (*mapperAdapted)[id].previous.clear();

            (*mapperAdapted)[id].current.push_back(id);
            (*mapperAdapted)[id].type = m_previousmapper[info.previous[0]].type;
            (*mapperAdapted)[id].entity = adaption::Entity::ENTITY_CELL;
            (*mapperAdapted)[id].previous = m_previousmapper[info.previous[0]].previous;

            if (fillInv){
                for (long idp : (*mapperAdapted)[id].previous){
                    std::vector<long>::iterator it = std::find((*mapperMapped)[idp].previous.begin(), (*mapperMapped)[idp].previous.end(), info.previous[0]);
                    assert(it != (*mapperMapped)[idp].previous.end());
                    *it = id;
                }
            }
        }
        if (info.type == adaption::Type::TYPE_REFINEMENT){
            std::size_t iprevious = 0;
            long idprevious = m_previousmapper[info.previous[0]].previous[iprevious];
            VolOctree::OctantInfo oinfoprev = meshMapped->getCellOctant(idprevious);
            uint64_t morton = meshMapped->getTree().getMorton(oinfoprev.id);
            for (const long & id : info.current){
                (*mapperAdapted)[id].current.clear();
                (*mapperAdapted)[id].previous.clear();
                (*mapperAdapted)[id].entity = adaption::Entity::ENTITY_CELL;
                (*mapperAdapted)[id].current.push_back(id);
                uint8_t level =  meshAdapted->getCellLevel(id);
                if (level == meshMapped->getCellLevel(idprevious)){
                    (*mapperAdapted)[id].type = adaption::Type::TYPE_RENUMBERING;
                    (*mapperAdapted)[id].previous.push_back(idprevious);

                    if (fillInv){
                        (*mapperMapped)[idprevious].type = adaption::Type::TYPE_RENUMBERING;
                        (*mapperMapped)[idprevious].previous.clear();
                        (*mapperMapped)[idprevious].previous.push_back(id);
                    }
                }
                else if (level > meshMapped->getCellLevel(idprevious)){
                    (*mapperAdapted)[id].type = adaption::Type::TYPE_REFINEMENT;
                    (*mapperAdapted)[id].previous.push_back(idprevious);

                    if (fillInv){
                        std::vector<long>::iterator it = std::find((*mapperMapped)[idprevious].previous. begin(), (*mapperMapped)[idprevious].previous.end(), info.previous[0]);
                        if (it != (*mapperMapped)[idprevious].previous.end()){
                            *it = id;
                        }
                        else{
                            (*mapperMapped)[idprevious].previous.push_back(id);
                        }
                        (*mapperMapped)[idprevious].type = adaption::Type::TYPE_COARSENING;
                    }
                }
                else if (level < meshMapped->getCellLevel(idprevious)){
                    (*mapperAdapted)[id].type = adaption::Type::TYPE_COARSENING;
                    VolOctree::OctantInfo oinfo = meshAdapted->getCellOctant(id);
                    uint64_t mortonlastdesc = meshAdapted->getTree().getLastDescMorton(oinfo.id);
                    while (morton <= mortonlastdesc && iprevious < m_previousmapper[info.previous[0]].previous.size()){

                        idprevious = m_previousmapper[info.previous[0]].previous[iprevious];
                        oinfoprev = meshMapped->getCellOctant(idprevious);
                        morton = meshMapped->getTree().getMorton(oinfoprev.id);
                        (*mapperAdapted)[id].previous.push_back(idprevious);

                        if (fillInv){
                            (*mapperMapped)[idprevious].type = adaption::Type::TYPE_REFINEMENT;
                            (*mapperMapped)[idprevious].previous[0] = id;
                        }

                        iprevious++;
                    }
                }
            }
        }
        if (info.type == adaption::Type::TYPE_COARSENING){
            long id = info.current[0];
            (*mapperAdapted)[id].current.clear();
            (*mapperAdapted)[id].previous.clear();
            (*mapperAdapted)[id].entity = adaption::Entity::ENTITY_CELL;
            (*mapperAdapted)[id].current.push_back(id);
            uint8_t level =  meshAdapted->getCellLevel(id);
            std::unordered_set<long> idsprev;
            for (const long & id_ : info.previous){
                for (const long & idp : m_previousmapper[id_].previous){
                    idsprev.insert(idp);
                }
            }
            if (idsprev.size() == 1){
                long idprevious = *idsprev.begin();
                if (level == meshMapped->getCellLevel(idprevious)){
                    (*mapperAdapted)[id].type = adaption::Type::TYPE_RENUMBERING;

                    if (fillInv){
                        (*mapperMapped)[idprevious].type = adaption::Type::TYPE_RENUMBERING;
                        (*mapperMapped)[idprevious].previous.clear();
                        (*mapperMapped)[idprevious].previous.push_back(id);
                    }
                }
                else{
                    (*mapperAdapted)[id].type = adaption::Type::TYPE_REFINEMENT;

                    if (fillInv){
                        std::vector<long>::iterator it = std::find((*mapperMapped)[idprevious].previous. begin(), (*mapperMapped)[idprevious].previous.end(), info.previous[0]);
                        assert(it != (*mapperMapped)[idprevious].previous.end());
                        *it = id;
                        (*mapperMapped)[idprevious].type = adaption::Type::TYPE_COARSENING;
                    }

                }
            }
            else{
                (*mapperAdapted)[id].type = adaption::Type::TYPE_COARSENING;

                if (fillInv){
                    for (const long & idprevious : idsprev){
                        (*mapperMapped)[idprevious].type = adaption::Type::TYPE_REFINEMENT;
                        (*mapperMapped)[idprevious].previous[0] = id;
                    }
                }

            }
            for (const long & id_ : idsprev){
                (*mapperAdapted)[id].previous.push_back(id_);
            }
        }
    }

}



/**
 * Map an input VolOctree mesh on a VolOctree reference mesh.
 * Requirement : the meshes have to be defined on the same identical domain.
 *
 * \param[in] meshReference Pointer to reference mesh
 * \param[in] meshMapped Pointer to input mesh to map
 * \param[in] fillInv If true even the inverse mapping (reference mesh to input mesh) is filled.
 */
void MeshMapper::_mapMeshes(bitpit::VolOctree * meshReference, bitpit::VolOctree * meshMapped, bool fillInv)
{

    if ( (meshReference->getLength() != meshMapped->getLength()) || (meshReference->getOrigin() != meshMapped->getOrigin()) )
        throw std::runtime_error ("mesh mapper: different domain of VolOctree meshes not allowed.");

    m_mapper.setDynamicKernel(&meshReference->getCells(), PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED);
    if (fillInv)
        m_invmapper.setDynamicKernel(&meshMapped->getCells(), PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED);
    else
        clearInverseMapping();

#if BITPIT_ENABLE_MPI==1

    if (!(meshReference->isPartitioned())){
#endif

        _mapMeshesSamePartition(meshReference, meshMapped, fillInv);

#if BITPIT_ENABLE_MPI==1
    }
    else{

        std::vector<uint64_t> partitionReference = meshReference->getTree().getPartitionLastDesc();
        std::vector<uint64_t> partitionMapped = meshMapped->getTree().getPartitionLastDesc();

        bool checkPartition = (partitionReference[m_rank] == partitionMapped[m_rank]);
        if (m_rank != 0)
            checkPartition = checkPartition && (partitionReference[m_rank-1] == partitionMapped[m_rank-1]);

        if (checkPartition){
            _mapMeshesSamePartition(meshReference, meshMapped, fillInv);
        }
        else{
            throw std::runtime_error("Mapping with different partitions not supported.");
        }

    }
#endif

}

/**
 * Map an input VolOctree mesh on a VolOctree reference mesh.
 * Requirements : the meshes have to be defined on the same identical domain;
 *                the meshes must have the same parallel partitioning in terms of last descendant octant.
 *
 * \param[in] meshReference Pointer to reference mesh
 * \param[in] meshMapped Pointer to input mesh to map
 * \param[in] fillInv If true even the inverse mapping (reference mesh to input mesh) is filled.
 */
void MeshMapper::_mapMeshesSamePartition(bitpit::VolOctree * meshReference, bitpit::VolOctree * meshMapped, bool fillInv)
{

    long nRef  = meshReference->getInternalCount();
    long nMap  = meshMapped->getInternalCount();

    long indRef  = 0;
    long indMap = 0;

    while (indRef < nRef && indMap < nMap){
        VolOctree::OctantInfo octinfoRef(indRef, true);
        long idRef = meshReference->getOctantId(octinfoRef);
        VolOctree::OctantInfo octinfoMap(indMap, true);
        long idMap = meshMapped->getOctantId(octinfoMap);

        m_mapper[idRef].entity = bitpit::adaption::Entity::ENTITY_CELL;
        m_mapper[idRef].rank = m_rank;
        if (fillInv){
            m_invmapper[idMap].entity = bitpit::adaption::Entity::ENTITY_CELL;
            m_invmapper[idMap].rank = m_rank;
        }

        if (meshMapped->getCellLevel(idMap) == meshReference->getCellLevel(idRef)){
            m_mapper[idRef].current.push_back(idRef);
            m_mapper[idRef].previous.push_back(idMap);
            m_mapper[idRef].type = bitpit::adaption::Type::TYPE_RENUMBERING;
            if (fillInv){
                m_invmapper[idMap].current.push_back(idMap);
                m_invmapper[idMap].previous.push_back(idRef);
                m_invmapper[idMap].type = bitpit::adaption::Type::TYPE_RENUMBERING;
            }
            indRef++;
            indMap++;
        }
        else if (meshMapped->getCellLevel(idMap) > meshReference->getCellLevel(idRef)){
            m_mapper[idRef].current.push_back(idRef);
            m_mapper[idRef].type = bitpit::adaption::Type::TYPE_COARSENING;

            uint64_t mortonlastdesc = meshReference->getTree().getLastDescMorton(indRef);
            uint64_t mortonMap = meshMapped->getTree().getMorton(indMap);

            while(mortonMap <= mortonlastdesc && indMap < nMap){
                m_mapper[idRef].previous.push_back(idMap);
                if (fillInv){
                    m_invmapper[idMap].current.push_back(idMap);
                    m_invmapper[idMap].type = bitpit::adaption::Type::TYPE_REFINEMENT;
                    m_invmapper[idMap].previous.push_back(idRef);
                }
                indMap++;
                octinfoMap = VolOctree::OctantInfo(indMap, true);
                idMap = meshMapped->getOctantId(octinfoMap);
                mortonMap = meshMapped->getTree().getMorton(indMap);
            }
            indRef++;
        }
        else if (meshMapped->getCellLevel(idMap) < meshReference->getCellLevel(idRef)){

            uint64_t morton= meshReference->getTree().getMorton(indRef);
            uint64_t mortonlastdescmesh = meshMapped->getTree().getLastDescMorton(indMap);

            if (fillInv){
                m_invmapper[idMap].current.push_back(idMap);
                m_invmapper[idMap].type = bitpit::adaption::Type::TYPE_COARSENING;
            }

            while (morton <= mortonlastdescmesh && indRef < nRef){
                m_mapper[idRef].current.push_back(idRef);
                m_mapper[idRef].type = bitpit::adaption::Type::TYPE_REFINEMENT;
                m_mapper[idRef].previous.push_back(idMap);
                if (fillInv){
                    m_invmapper[idMap].previous.push_back(idRef);
                }
                indRef++;
                morton = meshReference->getTree().getMorton(indRef);
                octinfoRef = VolOctree::OctantInfo(indRef, true);
                idRef = meshReference->getOctantId(octinfoRef);
            }
            indMap++;
        }
    }

}


#if BITPIT_ENABLE_MPI
/**
 * Initializes the MPI communicator to be used for parallel communications.
 *
 * \param communicator is the communicator.
 */
void MeshMapper::initializeCommunicator(MPI_Comm communicator)
{
    // Communication can be set just once
    if (isCommunicatorSet())
        throw std::runtime_error ("MeshMapper communicator can be set just once");

    // The communicator has to be valid
    if (communicator == MPI_COMM_NULL)
        throw std::runtime_error ("MeshMapper communicator is not valid");

    // Create a copy of the user-specified communicator
    //
    // No library routine should use MPI_COMM_WORLD as the communicator;
    // instead, a duplicate of a user-specified communicator should always
    // be used.
    MPI_Comm_dup(communicator, &m_communicator);
}

/**
 * Returns the MPI communicator stored within LevelSetKernel.
 * @return MPI communicator.
 */
MPI_Comm MeshMapper::getCommunicator() const
{
    return m_communicator;
}

/**
 * Checks if the communicator to be used for parallel communications has
 * already been set.
 *
 * \result Returns true if the communicator has been set, false otherwise.
 */
bool MeshMapper::isCommunicatorSet() const
{

    return (getCommunicator() != MPI_COMM_NULL);
}

/**
 * Frees the MPI communicator associated to the patch
 */
void MeshMapper::freeCommunicator()
{
    if (!isCommunicatorSet())
        return;

    int finalizedCalled;
    MPI_Finalized(&finalizedCalled);
    if (finalizedCalled)
        return;

    MPI_Comm_free(&m_communicator);
}
#endif


}
