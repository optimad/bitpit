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

/**
 * \example mapper_voloctree_example_00004.cpp
 *
 * \brief Mesh mapper updating computing between partitioned voloctree meshes.
 * <b>To run</b>: ./mapper_voloctree_example_00004 \n
 */

#include <ctime>
#include <array>
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_POD.hpp"
#include "bitpit_patchkernel.hpp"
#include "bitpit_voloctree.hpp"

using namespace bitpit;
using namespace pod;

/**
 * Run the examples.
 */
void runReferenceAdaptation()
{

    log::cout() << "  >> 2D octree patch" << "\n";

    /**
     * Create the tree
     */
    double x_0 = 10.;
    double y_0 = 20.;
    double z_0 = 30.;
    double l   = 1.5;

    std::unique_ptr<PabloUniform> treePointer = std::unique_ptr<PabloUniform>(new PabloUniform(x_0, y_0, z_0, l, 2));
    PabloUniform &octree = *treePointer;

    std::cout << " Origin : ( " << octree.getX0() << ", " << octree.getY0() << ", " << octree.getZ0() << " )" << std::endl;
    std::cout << " Length : " << octree.getL() << std::endl;

    /** Refine and write the octree*/
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();

    /**
     * Create the patch from the existing tree
     */
    /** Create the original patch */
    VolOctree *patch_2D_original = new VolOctree(std::move(treePointer), &treePointer);

    patch_2D_original->buildAdjacencies();
    patch_2D_original->buildInterfaces();

    patch_2D_original->update();

#if BITPIT_ENABLE_MPI==1
    /** Partition the patch */
    patch_2D_original->partition(true);
#endif

    std::vector<uint64_t> refineList;
    refineList.push_back(  7);
    refineList.push_back( 13);
    refineList.push_back( 15);
    refineList.push_back( 26);
    refineList.push_back( 27);
    refineList.push_back( 31);
    refineList.push_back( 37);
    refineList.push_back( 39);
    refineList.push_back( 49);
    refineList.push_back( 50);
    refineList.push_back( 51);
    refineList.push_back( 53);
    refineList.push_back( 55);
    refineList.push_back( 63);
    refineList.push_back( 78);
    refineList.push_back(100);
    refineList.push_back(102);
    refineList.push_back(105);
    refineList.push_back(108);
    refineList.push_back(109);
    refineList.push_back(110);
    refineList.push_back(135);
    refineList.push_back(141);
    refineList.push_back(143);
    refineList.push_back(146);
    refineList.push_back(147);
    refineList.push_back(151);
    refineList.push_back(153);
    refineList.push_back(154);
    refineList.push_back(155);
    refineList.push_back(157);
    refineList.push_back(159);
    refineList.push_back(165);
    refineList.push_back(167);
    refineList.push_back(183);
    refineList.push_back(198);
    refineList.push_back(204);
    refineList.push_back(206);
    refineList.push_back(225);
    refineList.push_back(228);
    refineList.push_back(229);
    refineList.push_back(230);


    for (uint64_t ind : refineList) {
#if BITPIT_ENABLE_MPI==1
        int owner = patch_2D_original->getTree().getOwnerRank(ind);
        if (patch_2D_original->getRank() == owner){
            uint32_t lind = patch_2D_original->getTree().getLocalIdx(ind, owner);
#else
            uint32_t lind = patch_2D_original->getTree().getLocalIdx(ind);
#endif
            VolOctree::OctantInfo octinfo(lind, true);
            long id = patch_2D_original->getOctantId(octinfo);
            patch_2D_original->markCellForRefinement(id);
#if BITPIT_ENABLE_MPI==1
        }
#endif
    }

    patch_2D_original->update(true);

    /**
     * Create the new tree
     */
    std::unique_ptr<PabloUniform> treePointer2 = std::unique_ptr<PabloUniform>(new PabloUniform(x_0, y_0, z_0, l, 2));
    PabloUniform &octree2 = *treePointer2;

    /** Refine and write the octree */
    octree2.adaptGlobalRefine();
    octree2.adaptGlobalRefine();
    octree2.adaptGlobalRefine();
    octree2.adaptGlobalRefine();

    /** Create a new patch */
    VolOctree *patch_2D = new VolOctree(std::move(treePointer2), &treePointer2);
#if BITPIT_ENABLE_MPI==1
    patch_2D->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);
#endif

    /** Adapt the patch */
    for (int k = 0; k < 10; ++k) {
        long nCells = patch_2D->getCellCount();
        log::cout() << std::endl;
        log::cout() << ">> Marking the cells to adapt... " << std::endl;

        for (int i = 0; i < 100; ++i) {
            long cellId;
#if BITPIT_ENABLE_MPI==1
            if (patch_2D->getRank() == 0)
#endif
                cellId = rand() % nCells * 2;
#if BITPIT_ENABLE_MPI==1
            MPI_Bcast(&cellId, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

            if (!patch_2D->getCells().exists(cellId)) {
                continue;
            }

            patch_2D->markCellForRefinement(cellId);
        }

        for (int i = 0; i < 55; ++i) {
            long cellId;
#if BITPIT_ENABLE_MPI==1
            if (patch_2D->getRank() == 0)
#endif
                cellId = rand() % nCells * 2;
#if BITPIT_ENABLE_MPI==1
            MPI_Bcast(&cellId, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

            if (!patch_2D->getCells().exists(cellId)) {
                continue;
            }

            patch_2D->markCellForCoarsening(cellId);
            for (auto neighId : patch_2D->findCellNeighs(cellId)) {
                patch_2D->markCellForCoarsening(neighId);
            }
        }

        log::cout() << std::endl;
        log::cout() << ">> Initial number of cells... " << nCells << std::endl;

        patch_2D->buildAdjacencies();
        patch_2D->buildInterfaces();

        patch_2D->update(true);

        nCells = patch_2D->getCellCount();
        log::cout() << ">> Final number of cells... " << nCells << std::endl;
    }

#if BITPIT_ENABLE_MPI==1
    /** Partition the patch */
    patch_2D->partition(true);
#endif

    /** Show patch info */
    log::cout() << "Cell count:   " << patch_2D->getCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch_2D->getVertexCount() << std::endl;

    patch_2D->getVTK().setName("mesh_random_testR.0");
    patch_2D->write();


    /** Create mapper object */
    VolOctreeMapper mapobject(patch_2D, patch_2D_original);

    /** Map the two meshes */
    mapobject.initialize(true);

    /** Define data on original mesh and write */
    PiercedStorage<double> data(1, &patch_2D_original->getCells());
    std::vector<double> vdata(patch_2D_original->getInternalCount());
    int count = 0;
    for (Cell & cell : patch_2D_original->getCells()){
        if (cell.isInterior()){
            long id = cell.getId();
            data[id] = double(patch_2D_original->getCellLevel(id));
            vdata[count] = data[id];
            count++;
        }
    }

    patch_2D_original->getVTK().setName("mesh_original_testR.0");
    patch_2D_original->getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::CELL, vdata);
#if BITPIT_ENABLE_MPI==1
    patch_2D_original->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);
#endif
    patch_2D_original->write();


    /** Map data on second mesh and write */
    {
        PiercedStorage<double> data2(1, &patch_2D->getCells());
        data2.fill(0.0);

        const PiercedStorage<mapping::Info> & mapper = mapobject.getMapping();

#if BITPIT_ENABLE_MPI==1
        //Communicate needed data for mapping
        std::map<int, std::map<long, double> > datarec;
        std::map<int, std::map<long, double> > volrec;
        {
            std::map<int, std::vector<long> > rankIDrec = mapobject.getReceivedMappedIds();
            std::map<int, std::vector<long> > rankIDsend = mapobject.getSentMappedIds();

            //build send buffers
            MPI_Comm comm = MPI_COMM_WORLD;
            DataCommunicator dataCommunicator(comm);
            MPI_Barrier(comm);
            std::size_t bytes = uint8_t(2*sizeof(double));
            for (std::pair<const int, std::vector<long int> > val : rankIDsend){
                int rank = val.first;
                //set size
                std::size_t buffSize = val.second.size() * bytes;
                dataCommunicator.setSend(rank,buffSize);
                //fill buffer with octants
                SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(rank);
                for (long & ID : val.second){
                    sendBuffer << data[ID];
                    sendBuffer << patch_2D_original->evalCellVolume(ID);
                }
            }

            dataCommunicator.discoverRecvs();
            dataCommunicator.startAllRecvs();
            dataCommunicator.startAllSends();

            for (std::pair<const int, std::vector<long int> > val : rankIDrec){
                int rank = val.first;
                dataCommunicator.waitRecv(rank);
                RecvBuffer & recvBuffer = dataCommunicator.getRecvBuffer(rank);
                for (long & ID : val.second){
                    recvBuffer >> datarec[rank][ID];
                    recvBuffer >> volrec[rank][ID];
                }
            }
            dataCommunicator.waitAllSends();

        }
#endif

        std::vector<double> vdata2(patch_2D->getInternalCount());
        count = 0;
        for (Cell & cell : patch_2D->getCells()){
            if (cell.isInterior()){
                long id = cell.getId();
                if (mapper[id].type == mapping::Type::TYPE_RENUMBERING){
#if BITPIT_ENABLE_MPI==1
                    if (mapper[id].ranks[0] == patch_2D->getRank()){
#endif
                        data2[id] = data[mapper[id].ids[0]];
                        vdata2[count] = data2[id];
#if BITPIT_ENABLE_MPI==1
                    }
                    else{
                        data2[id] = datarec[mapper[id].ranks[0]][mapper[id].ids[0]];
                        vdata2[count] = data2[id];
                    }
#endif
                }
                else if (mapper[id].type == mapping::Type::TYPE_COARSENING){
                    data2[id] = 0.0;
                    double Volume = 0.0;
                    int i = 0;
                    for (long idd : mapper[id].ids){
#if BITPIT_ENABLE_MPI==1
                        if (mapper[id].ranks[i] == patch_2D->getRank()){
#endif
                            data2[id] += data[idd] * patch_2D_original->evalCellVolume(idd);
                            Volume += patch_2D_original->evalCellVolume(idd);
#if BITPIT_ENABLE_MPI==1
                        }
                        else{
                            data2[id] += datarec[mapper[id].ranks[i]][idd] * volrec[mapper[id].ranks[i]][idd];
                            Volume += volrec[mapper[id].ranks[i]][idd];
                        }
#endif
                        i++;
                    }
                    data2[id] /= Volume;
                    vdata2[count] = data2[id];
                }
                else if (mapper[id].type == mapping::Type::TYPE_REFINEMENT){
#if BITPIT_ENABLE_MPI==1
                    if (mapper[id].ranks[0] == patch_2D->getRank()){
#endif
                        data2[id] = data[mapper[id].ids[0]];
                        vdata2[count] = data2[id];
#if BITPIT_ENABLE_MPI==1
                    }
                    else{
                        data2[id] = datarec[mapper[id].ranks[0]][mapper[id].ids[0]];
                        vdata2[count] = data2[id];
                    }
#endif
                }
                count++;
            }
        }

        patch_2D->getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::CELL, vdata2);
        patch_2D->getVTK().setName("mesh_random_testR.0");
        patch_2D->write();

    }

    /** Adapt again and update the mapper */
    long nCells = patch_2D->getCellCount();
    log::cout() << std::endl;
    log::cout() << ">> Marking the cells to adapt... " << std::endl;

    srand(1517497123);

    for (int i = 0; i < 30; ++i) {
        long cellId;
#if BITPIT_ENABLE_MPI==1
        if (patch_2D->getRank() == 0)
#endif
            cellId = rand() % nCells * 2;
#if BITPIT_ENABLE_MPI==1
        MPI_Bcast(&cellId, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

        if (!patch_2D->getCells().exists(cellId)) {
            continue;
        }

        patch_2D->markCellForRefinement(cellId);

    }

    for (int i = 0; i < 60; ++i) {
        long cellId;
#if BITPIT_ENABLE_MPI==1
        if (patch_2D->getRank() == 0)
#endif
            cellId = rand() % nCells * 2;
#if BITPIT_ENABLE_MPI==1
        MPI_Bcast(&cellId, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

        if (!patch_2D->getCells().exists(cellId)) {
            continue;
        }

        patch_2D->markCellForCoarsening(cellId);
        for (auto neighId : patch_2D->findCellNeighs(cellId)) {
            patch_2D->markCellForCoarsening(neighId);
        }

    }

    log::cout() << std::endl;
    log::cout() << ">> Initial number of cells... " << nCells << std::endl;

    std::vector<adaption::Info> infoAdapt = patch_2D->adaptionPrepare(true);

    mapobject.adaptionPrepare(infoAdapt, true);

    infoAdapt = patch_2D->adaptionAlter(true);

    mapobject.adaptionAlter(infoAdapt, true, true);

    patch_2D->adaptionCleanup();

    mapobject.adaptionCleanup();

    nCells = patch_2D->getCellCount();
    log::cout() << ">> Final number of cells... " << nCells << std::endl;

    /** Re-Map data on second mesh and write */
    PiercedStorage<double> data2(1, &patch_2D->getCells());
    data2.fill(0.0);
    {
        const PiercedStorage<mapping::Info> & mapper = mapobject.getMapping();

#if BITPIT_ENABLE_MPI==1
        //Communicate needed data for mapping
        std::map<int, std::map<long, double> > datarec;
        std::map<int, std::map<long, double> > volrec;
        {
            std::map<int, std::vector<long> > rankIDrec = mapobject.getReceivedMappedIds();
            std::map<int, std::vector<long> > rankIDsend = mapobject.getSentMappedIds();

            //build send buffers
            MPI_Comm comm = MPI_COMM_WORLD;
            DataCommunicator dataCommunicator(comm);
            MPI_Barrier(comm);
            std::size_t bytes = uint8_t(2*sizeof(double));
            for (std::pair<const int, std::vector<long int> > val : rankIDsend){
                int rank = val.first;
                //set size
                std::size_t buffSize = val.second.size() * bytes;
                dataCommunicator.setSend(rank,buffSize);
                //fill buffer with octants
                SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(rank);
                for (long & ID : val.second){
                    sendBuffer << data[ID];
                    sendBuffer << patch_2D_original->evalCellVolume(ID);
                }
            }

            dataCommunicator.discoverRecvs();
            dataCommunicator.startAllRecvs();
            dataCommunicator.startAllSends();

            for (std::pair<const int, std::vector<long int> > val : rankIDrec){
                int rank = val.first;
                dataCommunicator.waitRecv(rank);
                RecvBuffer & recvBuffer = dataCommunicator.getRecvBuffer(rank);
                for (long & ID : val.second){
                    recvBuffer >> datarec[rank][ID];
                    recvBuffer >> volrec[rank][ID];
                }
            }
            dataCommunicator.waitAllSends();
        }
#endif

        std::vector<double> vdata2(patch_2D->getInternalCount());
        count = 0;
        for (Cell & cell : patch_2D->getCells()){
            if (cell.isInterior()){
                long id = cell.getId();
                if (mapper[id].type == mapping::Type::TYPE_RENUMBERING){
#if BITPIT_ENABLE_MPI==1
                    if (mapper[id].ranks[0] == patch_2D->getRank()){
#endif
                        data2[id] = data[mapper[id].ids[0]];
                        vdata2[count] = data2[id];
#if BITPIT_ENABLE_MPI==1
                    }
                    else{
                        data2[id] = datarec[mapper[id].ranks[0]][mapper[id].ids[0]];
                        vdata2[count] = data2[id];
                    }
#endif
                }
                else if (mapper[id].type == mapping::Type::TYPE_COARSENING){
                    data2[id] = 0.0;
                    double Volume = 0.0;
                    int i = 0;
                    for (long idd : mapper[id].ids){
#if BITPIT_ENABLE_MPI==1
                        if (mapper[id].ranks[i] == patch_2D->getRank()){
#endif
                            data2[id] += data[idd] * patch_2D_original->evalCellVolume(idd);
                            Volume += patch_2D_original->evalCellVolume(idd);
#if BITPIT_ENABLE_MPI==1
                        }
                        else{
                            data2[id] += datarec[mapper[id].ranks[i]][idd] * volrec[mapper[id].ranks[i]][idd];
                            Volume += volrec[mapper[id].ranks[i]][idd];
                        }
#endif
                        i++;
                    }
                    data2[id] /= Volume;
                    vdata2[count] = data2[id];
                }
                else if (mapper[id].type == mapping::Type::TYPE_REFINEMENT){
#if BITPIT_ENABLE_MPI==1
                    if (mapper[id].ranks[0] == patch_2D->getRank()){
#endif
                        data2[id] = data[mapper[id].ids[0]];
                        vdata2[count] = data2[id];
#if BITPIT_ENABLE_MPI==1
                    }
                    else{
                        data2[id] = datarec[mapper[id].ranks[0]][mapper[id].ids[0]];
                        vdata2[count] = data2[id];
                    }
#endif
                }
                count++;
            }
        }

        patch_2D->getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::CELL, vdata2);
        patch_2D->getVTK().setName("mesh_random_testR.1");
        patch_2D->write();

    }

    /** Inverse Mapping data on first mesh with inverse mapping and write */
    PiercedStorage<double> data2inv(1, &patch_2D_original->getCells());
    data2inv.fill(0.0);
    {
        const PiercedStorage<mapping::Info> & mapper = mapobject.getInverseMapping();

#if BITPIT_ENABLE_MPI==1
        //Communicate needed data for mapping
        std::map<int, std::map<long, double> > datarec;
        std::map<int, std::map<long, double> > volrec;
        {

            // recover id to be rec/sent
            std::map<int, std::set<long> > rankIDsend;

            const PiercedStorage<mapping::Info> & _mapper = mapobject.getMapping();

            for (Cell & cell : patch_2D->getCells()){
                long ID = cell.getId();
                auto info = _mapper[ID];
                int i = 0;
                for (int rank : info.ranks){
                    if (rank != patch_2D->getRank()){
                        rankIDsend[rank].insert(ID);
                    }
                    i++;
                }
            }

            //build send buffers
            MPI_Comm comm = MPI_COMM_WORLD;
            DataCommunicator dataCommunicator(comm);
            MPI_Barrier(comm);
            std::size_t bytes = uint8_t(2*sizeof(double) + sizeof(long));
            for (std::pair<const int, std::set<long int> > val : rankIDsend){
                int rank = val.first;
                //set size
                std::size_t buffSize = val.second.size() * bytes;
                dataCommunicator.setSend(rank,buffSize);
                //fill buffer with octants
                SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(rank);
                for (long ID : val.second){
                    sendBuffer << ID;
                    sendBuffer << data2[ID];
                    sendBuffer << patch_2D->evalCellVolume(ID);
                }
            }

            dataCommunicator.discoverRecvs();
            dataCommunicator.startAllRecvs();
            dataCommunicator.startAllSends();

            std::vector<int> recvRanks = dataCommunicator.getRecvRanks();
            std::sort(recvRanks.begin(),recvRanks.end());

            for (int rank : recvRanks){
                dataCommunicator.waitRecv(rank);
                RecvBuffer & recvBuffer = dataCommunicator.getRecvBuffer(rank);
                long nof = recvBuffer.getSize() / bytes;
                for (long ii = 0; ii < nof; ii++){
                    long ID;
                    recvBuffer >> ID;
                    recvBuffer >> datarec[rank][ID];
                    recvBuffer >> volrec[rank][ID];
                }
            }
            dataCommunicator.waitAllSends();

        }
#endif

        std::vector<double> vdata2inv(patch_2D_original->getInternalCount());
        count = 0;
        for (Cell & cell : patch_2D_original->getCells()){
            if (cell.isInterior()){
                long id = cell.getId();
                if (mapper[id].type == mapping::Type::TYPE_RENUMBERING){
#if BITPIT_ENABLE_MPI==1
                    if (mapper[id].ranks[0] == patch_2D_original->getRank()){
#endif
                        data2inv[id] = data2[mapper[id].ids[0]];
                        vdata2inv[count] = data2inv[id];
#if BITPIT_ENABLE_MPI==1
                    }
                    else{
                        data2inv[id] = datarec[mapper[id].ranks[0]][mapper[id].ids[0]];
                        vdata2inv[count] = data2inv[id];
                    }
#endif
                }
                else if (mapper[id].type == mapping::Type::TYPE_COARSENING){
                    data2inv[id] = 0.0;
                    double Volume = 0.0;
                    int i = 0;
                    for (long idd : mapper[id].ids){
#if BITPIT_ENABLE_MPI==1
                        if (mapper[id].ranks[i] == patch_2D_original->getRank()){
#endif
                            data2inv[id] += data2[idd] * patch_2D->evalCellVolume(idd);
                            Volume += patch_2D->evalCellVolume(idd);
#if BITPIT_ENABLE_MPI==1
                        }
                        else{
                            data2inv[id] += datarec[mapper[id].ranks[i]][idd] * volrec[mapper[id].ranks[i]][idd];
                            Volume += volrec[mapper[id].ranks[i]][idd];
                        }
#endif
                        i++;
                    }
                    data2inv[id] /= Volume;
                    vdata2inv[count] = data2inv[id];
                }
                else if (mapper[id].type == mapping::Type::TYPE_REFINEMENT){
#if BITPIT_ENABLE_MPI==1
                    if (mapper[id].ranks[0] == patch_2D_original->getRank()){
#endif
                        data2inv[id] = data2[mapper[id].ids[0]];
                        vdata2inv[count] = data2inv[id];
#if BITPIT_ENABLE_MPI==1
                    }
                    else{
                        data2inv[id] = datarec[mapper[id].ranks[0]][mapper[id].ids[0]];
                        vdata2inv[count] = data2inv[id];
                    }
#endif
                }
                count++;
            }
        }

        patch_2D_original->getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::CELL, vdata2inv);
        patch_2D_original->getVTK().setName("mesh_original_testR.1");
        patch_2D_original->write();

    }

}

void runMappedAdaptation()
{

    log::cout() << "  >> 2D octree patch" << "\n";

    /**
     * Create the tree
     */
    double x_0 = 10.;
    double y_0 = 20.;
    double z_0 = 30.;
    double l   = 1.5;

    std::unique_ptr<PabloUniform> treePointer = std::unique_ptr<PabloUniform>(new PabloUniform(x_0, y_0, z_0, l, 2));
    PabloUniform &octree = *treePointer;

    std::cout << " Origin : ( " << octree.getX0() << ", " << octree.getY0() << ", " << octree.getZ0() << " )" << std::endl;
    std::cout << " Length : " << octree.getL() << std::endl;

    /** Refine and write the octree*/
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();

    /**
     * Create the patch from the existing tree
     */
    /** Create the original patch */
    VolOctree *patch_2D_original = new VolOctree(std::move(treePointer), &treePointer);

    patch_2D_original->buildAdjacencies();
    patch_2D_original->buildInterfaces();

    patch_2D_original->update();

#if BITPIT_ENABLE_MPI==1
    /** Partition the patch */
    patch_2D_original->partition(true);
#endif

    std::vector<uint64_t> refineList;
    refineList.push_back(  7);
    refineList.push_back( 13);
    refineList.push_back( 15);
    refineList.push_back( 26);
    refineList.push_back( 27);
    refineList.push_back( 31);
    refineList.push_back( 37);
    refineList.push_back( 39);
    refineList.push_back( 49);
    refineList.push_back( 50);
    refineList.push_back( 51);
    refineList.push_back( 53);
    refineList.push_back( 55);
    refineList.push_back( 63);
    refineList.push_back( 78);
    refineList.push_back(100);
    refineList.push_back(102);
    refineList.push_back(105);
    refineList.push_back(108);
    refineList.push_back(109);
    refineList.push_back(110);
    refineList.push_back(135);
    refineList.push_back(141);
    refineList.push_back(143);
    refineList.push_back(146);
    refineList.push_back(147);
    refineList.push_back(151);
    refineList.push_back(153);
    refineList.push_back(154);
    refineList.push_back(155);
    refineList.push_back(157);
    refineList.push_back(159);
    refineList.push_back(165);
    refineList.push_back(167);
    refineList.push_back(183);
    refineList.push_back(198);
    refineList.push_back(204);
    refineList.push_back(206);
    refineList.push_back(225);
    refineList.push_back(228);
    refineList.push_back(229);
    refineList.push_back(230);


    for (uint64_t ind : refineList) {
#if BITPIT_ENABLE_MPI==1
        int owner = patch_2D_original->getTree().getOwnerRank(ind);
        if (patch_2D_original->getRank() == owner){
            uint32_t lind = patch_2D_original->getTree().getLocalIdx(ind, owner);
#else
            uint32_t lind = patch_2D_original->getTree().getLocalIdx(ind);
#endif
            VolOctree::OctantInfo octinfo(lind, true);
            long id = patch_2D_original->getOctantId(octinfo);
            patch_2D_original->markCellForRefinement(id);
#if BITPIT_ENABLE_MPI==1
        }
#endif
    }

    patch_2D_original->update(true);

    /**
     * Create the new tree
     */
    std::unique_ptr<PabloUniform> treePointer2 = std::unique_ptr<PabloUniform>(new PabloUniform(x_0, y_0, z_0, l, 2));
    PabloUniform &octree2 = *treePointer2;

    /** Refine and write the octree */
    octree2.adaptGlobalRefine();
    octree2.adaptGlobalRefine();
    octree2.adaptGlobalRefine();
    octree2.adaptGlobalRefine();

    /** Create a new patch */
    VolOctree *patch_2D = new VolOctree(std::move(treePointer2), &treePointer2);
#if BITPIT_ENABLE_MPI==1
    patch_2D->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);
#endif

    /** Adapt the patch */
    for (int k = 0; k < 10; ++k) {
        long nCells = patch_2D->getCellCount();
        log::cout() << std::endl;
        log::cout() << ">> Marking the cells to adapt... " << std::endl;

        for (int i = 0; i < 100; ++i) {
            long cellId;
#if BITPIT_ENABLE_MPI==1
            if (patch_2D->getRank() == 0)
#endif
                cellId = rand() % nCells * 2;
#if BITPIT_ENABLE_MPI==1
            MPI_Bcast(&cellId, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

            if (!patch_2D->getCells().exists(cellId)) {
                continue;
            }

            patch_2D->markCellForRefinement(cellId);
        }

        for (int i = 0; i < 55; ++i) {
            long cellId;
#if BITPIT_ENABLE_MPI==1
            if (patch_2D->getRank() == 0)
#endif
                cellId = rand() % nCells * 2;
#if BITPIT_ENABLE_MPI==1
            MPI_Bcast(&cellId, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

            if (!patch_2D->getCells().exists(cellId)) {
                continue;
            }

            patch_2D->markCellForCoarsening(cellId);
            for (auto neighId : patch_2D->findCellNeighs(cellId)) {
                patch_2D->markCellForCoarsening(neighId);
            }
        }

        log::cout() << std::endl;
        log::cout() << ">> Initial number of cells... " << nCells << std::endl;

        patch_2D->buildAdjacencies();
        patch_2D->buildInterfaces();

        patch_2D->update(true);

        nCells = patch_2D->getCellCount();
        log::cout() << ">> Final number of cells... " << nCells << std::endl;
    }


#if BITPIT_ENABLE_MPI==1
    /** Partition the patch */
    patch_2D->partition(true);
#endif

    /** Show patch info */
    log::cout() << "Cell count:   " << patch_2D->getCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch_2D->getVertexCount() << std::endl;


    patch_2D->getVTK().setName("mesh_random_testM.0");
    patch_2D->write();

    /** Define data on original mesh and write */
    PiercedStorage<double> data(1, &patch_2D_original->getCells());
    std::vector<double> vdata(patch_2D_original->getInternalCount());
    int count = 0;
    for (Cell & cell : patch_2D_original->getCells()){
        if (cell.isInterior()){
            long id = cell.getId();
            data[id] = double(patch_2D_original->getCellLevel(id));
            vdata[count] = data[id];
            count++;
        }
    }

    patch_2D_original->getVTK().setName("mesh_original_testM.0");
    patch_2D_original->getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::CELL, vdata);
#if BITPIT_ENABLE_MPI==1
    patch_2D_original->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);
#endif
    patch_2D_original->write();

    /** Create mapper object */
    VolOctreeMapper mapobject(patch_2D_original, patch_2D);

    /** Map the two meshes */
    mapobject.initialize(true);

    /** Inverse Mapping data on second mesh with inverse mapping and write */
    {
        PiercedStorage<double> data2(1, &patch_2D->getCells());
        data2.fill(0.0);

        const PiercedStorage<mapping::Info> & mapper = mapobject.getInverseMapping();

#if BITPIT_ENABLE_MPI==1
        //Communicate needed data for mapping
        std::map<int, std::map<long, double> > datarec;
        std::map<int, std::map<long, double> > volrec;
        {

            // recover id to be rec/sent
            std::map<int, std::set<long> > rankIDsend;

            const PiercedStorage<mapping::Info> & _mapper = mapobject.getMapping();

            for (Cell & cell : patch_2D_original->getCells()){
                long ID = cell.getId();
                auto info = _mapper[ID];
                int i = 0;
                for (int rank : info.ranks){
                    if (rank != patch_2D_original->getRank()){
                        rankIDsend[rank].insert(ID);
                    }
                    i++;
                }
            }

            //build send buffers
            MPI_Comm comm = MPI_COMM_WORLD;
            DataCommunicator dataCommunicator(comm);
            MPI_Barrier(comm);
            std::size_t bytes = uint8_t(2*sizeof(double) + sizeof(long));
            for (std::pair<const int, std::set<long int> > val : rankIDsend){
                int rank = val.first;
                //set size
                std::size_t buffSize = val.second.size() * bytes;
                dataCommunicator.setSend(rank,buffSize);
                //fill buffer with octants
                SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(rank);
                for (long ID : val.second){
                    sendBuffer << ID;
                    sendBuffer << data[ID];
                    sendBuffer << patch_2D_original->evalCellVolume(ID);
                }
            }

            dataCommunicator.discoverRecvs();
            dataCommunicator.startAllRecvs();
            dataCommunicator.startAllSends();

            std::vector<int> recvRanks = dataCommunicator.getRecvRanks();
            std::sort(recvRanks.begin(),recvRanks.end());

            for (int rank : recvRanks){
                dataCommunicator.waitRecv(rank);
                RecvBuffer & recvBuffer = dataCommunicator.getRecvBuffer(rank);
                long nof = recvBuffer.getSize() / bytes;
                for (long ii = 0; ii < nof; ii++){
                    long ID;
                    recvBuffer >> ID;
                    recvBuffer >> datarec[rank][ID];
                    recvBuffer >> volrec[rank][ID];
                }
            }
            dataCommunicator.waitAllSends();

        }
#endif

        std::vector<double> vdata2(patch_2D->getInternalCount());
        count = 0;
        for (Cell & cell : patch_2D->getCells()){
            if (cell.isInterior()){
                long id = cell.getId();
                if (mapper[id].type == mapping::Type::TYPE_RENUMBERING){
#if BITPIT_ENABLE_MPI==1
                    if (mapper[id].ranks[0] == patch_2D->getRank()){
#endif
                        data2[id] = data[mapper[id].ids[0]];
                        vdata2[count] = data2[id];
#if BITPIT_ENABLE_MPI==1
                    }
                    else{
                        data2[id] = datarec[mapper[id].ranks[0]][mapper[id].ids[0]];
                        vdata2[count] = data2[id];
                    }
#endif
                }
                else if (mapper[id].type == mapping::Type::TYPE_COARSENING){
                    data2[id] = 0.0;
                    double Volume = 0.0;
                    int i = 0;
                    for (long idd : mapper[id].ids){
#if BITPIT_ENABLE_MPI==1
                        if (mapper[id].ranks[i] == patch_2D->getRank()){
#endif
                            data2[id] += data[idd] * patch_2D_original->evalCellVolume(idd);
                            Volume += patch_2D_original->evalCellVolume(idd);
#if BITPIT_ENABLE_MPI==1
                        }
                        else{
                            data2[id] += datarec[mapper[id].ranks[i]][idd] * volrec[mapper[id].ranks[i]][idd];
                            Volume += volrec[mapper[id].ranks[i]][idd];
                        }
#endif
                        i++;
                    }
                    data2[id] /= Volume;
                    vdata2[count] = data2[id];
                }
                else if (mapper[id].type == mapping::Type::TYPE_REFINEMENT){
#if BITPIT_ENABLE_MPI==1
                    if (mapper[id].ranks[0] == patch_2D->getRank()){
#endif
                        data2[id] = data[mapper[id].ids[0]];
                        vdata2[count] = data2[id];
#if BITPIT_ENABLE_MPI==1
                    }
                    else{
                        data2[id] = datarec[mapper[id].ranks[0]][mapper[id].ids[0]];
                        vdata2[count] = data2[id];
                    }
#endif
                }
                count++;
            }
        }

        patch_2D->getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::CELL, vdata2);
        patch_2D->getVTK().setName("mesh_random_testM.0");
        patch_2D->write();

    }


    /** Adapt again and update the mapper */
    long nCells = patch_2D->getCellCount();
    log::cout() << std::endl;
    log::cout() << ">> Marking the cells to adapt... " << std::endl;

    srand(1517497123);

    for (int i = 0; i < 300; ++i) {
        long cellId;
#if BITPIT_ENABLE_MPI==1
        if (patch_2D->getRank() == 0)
#endif
            cellId = rand() % nCells * 2;
#if BITPIT_ENABLE_MPI==1
        MPI_Bcast(&cellId, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

        if (!patch_2D->getCells().exists(cellId)) {
            continue;
        }

        patch_2D->markCellForRefinement(cellId);

    }

    for (int i = 0; i < 200; ++i) {
        long cellId;
#if BITPIT_ENABLE_MPI==1
        if (patch_2D->getRank() == 0)
#endif
            cellId = rand() % nCells * 2;
#if BITPIT_ENABLE_MPI==1
        MPI_Bcast(&cellId, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

        if (!patch_2D->getCells().exists(cellId)) {
            continue;
        }

        patch_2D->markCellForCoarsening(cellId);
        for (auto neighId : patch_2D->findCellNeighs(cellId)) {
            patch_2D->markCellForCoarsening(neighId);
        }

    }

    log::cout() << std::endl;
    log::cout() << ">> Initial number of cells... " << nCells << std::endl;

    std::vector<adaption::Info> infoAdapt = patch_2D->adaptionPrepare(true);

    mapobject.adaptionPrepare(infoAdapt, false);

    infoAdapt = patch_2D->adaptionAlter(true);

    mapobject.adaptionAlter(infoAdapt, false, true);

    patch_2D->adaptionCleanup();

    mapobject.adaptionCleanup();

    nCells = patch_2D->getCellCount();
    log::cout() << ">> Final number of cells... " << nCells << std::endl;


    /** Re-Map data on second mesh with inverse mapping and write */
    PiercedStorage<double> data2(1, &patch_2D->getCells());
    data2.fill(0.0);
    {

        const PiercedStorage<mapping::Info> & mapper = mapobject.getInverseMapping();

#if BITPIT_ENABLE_MPI==1
        //Communicate needed data for mapping
        std::map<int, std::map<long, double> > datarec;
        std::map<int, std::map<long, double> > volrec;
        {

            // recover id to be rec/sent
            std::map<int, std::set<long> > rankIDsend;

            const PiercedStorage<mapping::Info> & _mapper = mapobject.getMapping();

            for (Cell & cell : patch_2D_original->getCells()){
                long ID = cell.getId();
                auto info = _mapper[ID];
                int i = 0;
                for (int rank : info.ranks){
                    if (rank != patch_2D_original->getRank()){
                        rankIDsend[rank].insert(ID);
                    }
                    i++;
                }
            }

            //build send buffers
            MPI_Comm comm;
            MPI_Comm_dup(MPI_COMM_WORLD, &comm);
            DataCommunicator dataCommunicator(comm);
            MPI_Barrier(comm);
            std::size_t bytes = uint8_t(2*sizeof(double) + sizeof(long));
            for (std::pair<const int, std::set<long int> > val : rankIDsend){
                int rank = val.first;
                //set size
                std::size_t buffSize = val.second.size() * bytes;
                dataCommunicator.setSend(rank,buffSize);
                //fill buffer with octants
                SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(rank);
                for (long ID : val.second){
                    sendBuffer << ID;
                    sendBuffer << data[ID];
                    sendBuffer << patch_2D_original->evalCellVolume(ID);
                }
            }

            dataCommunicator.discoverRecvs();
            dataCommunicator.startAllRecvs();
            dataCommunicator.startAllSends();

            std::vector<int> recvRanks = dataCommunicator.getRecvRanks();
            std::sort(recvRanks.begin(),recvRanks.end());

            for (int rank : recvRanks){
                dataCommunicator.waitRecv(rank);
                RecvBuffer & recvBuffer = dataCommunicator.getRecvBuffer(rank);
                long nof = recvBuffer.getSize() / bytes;
                for (long ii = 0; ii < nof; ii++){
                    long ID;
                    recvBuffer >> ID;
                    recvBuffer >> datarec[rank][ID];
                    recvBuffer >> volrec[rank][ID];
                }
            }
            dataCommunicator.waitAllSends();

        }
#endif

        std::vector<double> vdata2(patch_2D->getInternalCount());
        count = 0;
        for (Cell & cell : patch_2D->getCells()){
            if (cell.isInterior()){
                long id = cell.getId();
                if (mapper[id].type == mapping::Type::TYPE_RENUMBERING){
#if BITPIT_ENABLE_MPI==1
                    if (mapper[id].ranks[0] == patch_2D->getRank()){
#endif
                        data2[id] = data[mapper[id].ids[0]];
                        vdata2[count] = data2[id];
#if BITPIT_ENABLE_MPI==1
                    }
                    else{
                        data2[id] = datarec[mapper[id].ranks[0]][mapper[id].ids[0]];
                        vdata2[count] = data2[id];
                    }
#endif
                }
                else if (mapper[id].type == mapping::Type::TYPE_COARSENING){
                    data2[id] = 0.0;
                    double Volume = 0.0;
                    int i = 0;
                    for (long idd : mapper[id].ids){

                        if (id == 348){
                            std::cout << "mapped id : " << idd << std::endl;
                        }


#if BITPIT_ENABLE_MPI==1
                        if (mapper[id].ranks[i] == patch_2D->getRank()){
#endif
                            data2[id] += data[idd] * patch_2D_original->evalCellVolume(idd);
                            Volume += patch_2D_original->evalCellVolume(idd);
#if BITPIT_ENABLE_MPI==1
                        }
                        else{
                            data2[id] += datarec[mapper[id].ranks[i]][idd] * volrec[mapper[id].ranks[i]][idd];
                            Volume += volrec[mapper[id].ranks[i]][idd];
                        }
#endif
                        i++;
                    }
                    data2[id] /= Volume;
                    vdata2[count] = data2[id];
                }
                else if (mapper[id].type == mapping::Type::TYPE_REFINEMENT){
#if BITPIT_ENABLE_MPI==1
                    if (mapper[id].ranks[0] == patch_2D->getRank()){
#endif
                        data2[id] = data[mapper[id].ids[0]];
                        vdata2[count] = data2[id];
#if BITPIT_ENABLE_MPI==1
                    }
                    else{
                        data2[id] = datarec[mapper[id].ranks[0]][mapper[id].ids[0]];
                        vdata2[count] = data2[id];
                    }
#endif
                }
                count++;
            }
        }

        patch_2D->getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::CELL, vdata2);
        patch_2D->getVTK().setName("mesh_random_testM.1");
        patch_2D->write();

    }


    /** Final Re-Map data on first mesh and write */
    {
        PiercedStorage<double> data2inv(1, &patch_2D_original->getCells());
        data2inv.fill(0.0);

        const PiercedStorage<mapping::Info> & mapper = mapobject.getMapping();

#if BITPIT_ENABLE_MPI==1
        //Communicate needed data for mapping
        std::map<int, std::map<long, double> > datarec;
        std::map<int, std::map<long, double> > volrec;
        {
            std::map<int, std::vector<long> > rankIDrec = mapobject.getReceivedMappedIds();
            std::map<int, std::vector<long> > rankIDsend = mapobject.getSentMappedIds();

            //build send buffers
            MPI_Comm comm;
            MPI_Comm_dup(MPI_COMM_WORLD, &comm);
            DataCommunicator dataCommunicator(comm);
            MPI_Barrier(comm);
            std::size_t bytes = uint8_t(2*sizeof(double));
            for (std::pair<const int, std::vector<long int> > val : rankIDsend){
                int rank = val.first;
                //set size
                std::size_t buffSize = val.second.size() * bytes;
                dataCommunicator.setSend(rank,buffSize);
                //fill buffer with octants
                SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(rank);
                for (long & ID : val.second){
                    sendBuffer << data2[ID];
                    sendBuffer << patch_2D->evalCellVolume(ID);
                }
            }

            dataCommunicator.discoverRecvs();
            dataCommunicator.startAllRecvs();
            dataCommunicator.startAllSends();

            for (std::pair<const int, std::vector<long int> > val : rankIDrec){
                int rank = val.first;
                dataCommunicator.waitRecv(rank);
                RecvBuffer & recvBuffer = dataCommunicator.getRecvBuffer(rank);
                for (long & ID : val.second){
                    recvBuffer >> datarec[rank][ID];
                    recvBuffer >> volrec[rank][ID];
                }
            }
            dataCommunicator.waitAllSends();

        }
#endif

        std::vector<double> vdata2inv(patch_2D_original->getInternalCount());
        count = 0;
        for (Cell & cell : patch_2D_original->getCells()){
            if (cell.isInterior()){
                long id = cell.getId();
                if (mapper[id].type == mapping::Type::TYPE_RENUMBERING){
#if BITPIT_ENABLE_MPI==1
                    if (mapper[id].ranks[0] == patch_2D_original->getRank()){
#endif
                        data2inv[id] = data2[mapper[id].ids[0]];
                        vdata2inv[count] = data2inv[id];
#if BITPIT_ENABLE_MPI==1
                    }
                    else{
                        data2inv[id] = datarec[mapper[id].ranks[0]][mapper[id].ids[0]];
                        vdata2inv[count] = data2inv[id];
                    }
#endif
                }
                else if (mapper[id].type == mapping::Type::TYPE_COARSENING){
                    data2inv[id] = 0.0;
                    double Volume = 0.0;
                    int i = 0;
                    for (long idd : mapper[id].ids){
#if BITPIT_ENABLE_MPI==1
                        if (mapper[id].ranks[i] == patch_2D_original->getRank()){
#endif
                            data2inv[id] += data2[idd] * patch_2D->evalCellVolume(idd);
                            Volume += patch_2D->evalCellVolume(idd);
#if BITPIT_ENABLE_MPI==1
                        }
                        else{
                            data2inv[id] += datarec[mapper[id].ranks[i]][idd] * volrec[mapper[id].ranks[i]][idd];
                            Volume += volrec[mapper[id].ranks[i]][idd];
                        }
#endif
                        i++;
                    }
                    data2inv[id] /= Volume;
                    vdata2inv[count] = data2inv[id];
                }
                else if (mapper[id].type == mapping::Type::TYPE_REFINEMENT){
#if BITPIT_ENABLE_MPI==1
                    if (mapper[id].ranks[0] == patch_2D_original->getRank()){
#endif
                        data2inv[id] = data2[mapper[id].ids[0]];
                        vdata2inv[count] = data2inv[id];
#if BITPIT_ENABLE_MPI==1
                    }
                    else{
                        data2inv[id] = datarec[mapper[id].ranks[0]][mapper[id].ids[0]];
                        vdata2inv[count] = data2inv[id];
                    }
#endif
                }
                count++;
            }
        }

#if BITPIT_ENABLE_MPI==1
        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Barrier(comm);
#endif
        patch_2D_original->getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::CELL, vdata2inv);
        patch_2D_original->getVTK().setName("mesh_original_testM.1");
        patch_2D_original->write();

    }

}

/**
 * Main program.
 */

int main(int argc, char *argv[])
{

#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#endif


    /** Run the examples **/
    try {
        runReferenceAdaptation();
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    try {
        runMappedAdaptation();
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

}
