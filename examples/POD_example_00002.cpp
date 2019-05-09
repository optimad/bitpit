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
 * \example POD_example_00002.cpp
 * 
 * \brief POD basis computation using voloctree.
 * This example computes the POD basis starting from a database of simulations
 * defined on the same mesh and evaluate the reconstruction of a snapshot
 * not included in the database testing the hybrid interface.
 * The reconstruction is performed by methods based on POD field structures and
 * on PiercedStorage containers.
 * Finally, the reconstruction of a PiercedStorage field
 * is performed by using an adapted mesh field different from the POD mesh.
 * <b>To run</b>: ./POD_example_00002 \n
 */ 

#include <array>
#if BITPIT_ENABLE_MPI
#include <mpi.h>
#endif

#include "pod.hpp"
#include "bitpit_patchkernel.hpp"
#include "bitpit_voloctree.hpp"

using namespace bitpit;

/**
 * Run the example.
 */ 
void run(int rank, int nProcs)
{
    if (rank==0) 
        std::cout<< "0. Compute POD ..." << std::endl;
    /**<Create POD object.*/
    POD pod;

    /**<Add snapshots to database.*/   
    for (int i=1; i<11; i++)
        pod.addSnapshot("./data", "test."+std::to_string(i));

    /**<Set POD.*/    
    pod.setMeshType(POD::MeshType::VOLOCTREE);
    pod.setStaticMesh(true);
    pod.setWriteMode(POD::WriteMode::DEBUG);
    pod.setMemoryMode(POD::MemoryMode::MEMORY_NORMAL);
    pod.setModeCount(9);

    pod.setDirectory("pod");
    pod.setName("pod.test.solver");

    /**<Compute the POD basis.*/ 
    pod.run();

    if (rank==0) 
        std::cout<< "1. Restore POD ..." << std::endl;

    /**<Create restored POD object.*/
    POD podr;

    /**<Set restored POD.*/     
    podr.setDirectory("pod");
    podr.setName("pod.test.solver");
    podr.restore();

    if (rank==0) 
        std::cout<< "2. Read snapshot for reconstruction ..." << std::endl;

    VolumeKernel* meshr = new VolOctree();
    {
        int dumpBlock = (nProcs > 1) ? rank : -1;
        std::string filename = "./data/test.0.mesh";
        IBinaryArchive binaryReader(filename, dumpBlock);
#if BITPIT_ENABLE_MPI	
        meshr->setCommunicator(MPI_COMM_WORLD);
#endif
        meshr->restore(binaryReader.getStream());
        binaryReader.close();
    }

    pod::PODField fieldr, reconr;
    std::size_t nf;

    int dumpBlock = (nProcs > 1) ? rank : -1;
    std::string filename = "./data/test.0.data";
    IBinaryArchive binaryReader(filename, dumpBlock);
    std::istream &dataStream = binaryReader.getStream();

    /**<Restore solved cells.*/
    fieldr.mask = std::unique_ptr<PiercedStorage<bool>>(new PiercedStorage<bool>(1, &meshr->getCells()));
    fieldr.mask->restore(dataStream);

    /**<Restore scalar fields.*/
    std::size_t nsf;
    utils::binary::read(dataStream, nsf);
    std::vector<std::string> namesf;
    namesf.resize(nsf);
    for (std::size_t i = 0; i <nsf; ++i){
        utils::binary::read(dataStream, namesf[i]);
    }

    fieldr.scalar = std::unique_ptr<pod::ScalarStorage>(new pod::ScalarStorage(nsf, &meshr->getCells()));
    fieldr.scalar->restore(dataStream);

    /**<Restore vector fields.*/
    std::size_t nvf;
    utils::binary::read(dataStream, nvf);
    std::vector<std::array<std::string,3>> namevf;
    namevf.resize(nvf);
    for (std::size_t i = 0; i < nvf; ++i){
        utils::binary::read(dataStream, namevf[i]);
    }

    fieldr.vector = std::unique_ptr<pod::VectorStorage>(new pod::VectorStorage(nvf, &meshr->getCells()));
    fieldr.vector->restore(dataStream);
    binaryReader.close();

    nf=nsf+nvf;

    /**<Create sensors mask.*/
    pod.setSensorMask(*(fieldr.mask));

    /**<Reconstruct fields.*/
    pod.reconstructFields(fieldr,reconr);

    if (rank==0){
        std::cout<< ">> N modes:" << pod.getModeCount() << std::endl;

        std::cout<< ">> Reconstruction coeffs:" << std::endl;
        for (std::size_t i = 0; i < nf; ++i)
            std::cout<< pod.getReconstructionCoeffs()[i] << std::endl;
    }


    /**<Test hybrid interface using optimad solver format.*/
    PiercedStorage<double> gfield0(nsf+3*nvf+1, &meshr->getCells(), PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED);
    PiercedStorage<double> gfieldAdapt(nsf+3*nvf+1, &meshr->getCells(), PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED);
    std::unordered_set<long> targetCells;

    std::map<std::string, std::size_t> fields;

    for (Cell & cell : meshr->getCells()){
        long id = cell.getId();
        for (std::size_t ifield=0; ifield<nsf; ifield++){
            gfield0.at(id, ifield) = fieldr.scalar->at(id, ifield);
            gfieldAdapt.at(id, ifield) = fieldr.scalar->at(id, ifield);
        }

        for (std::size_t ifield=0; ifield<nvf; ifield++){
            for (std::size_t j=0; j<3; j++){
                gfield0.at(id, ifield*3+nsf+j) = fieldr.vector->at(id, ifield)[j];
                gfieldAdapt.at(id, ifield*3+nsf+j) = fieldr.vector->at(id, ifield)[j];
            }

        }
        if (fieldr.mask->at(id)){
            targetCells.insert(id);
            gfield0.at(id, nsf+3*nvf) = fieldr.mask->at(id);
            gfieldAdapt.at(id, nsf+3*nvf) = fieldr.mask->at(id);
        }
    }

    for (std::size_t ifield=0; ifield<nsf; ifield++){
        fields[namesf[ifield]] = ifield;
    }
    for (std::size_t ifield=0; ifield<nvf; ifield++){
        for (std::size_t j=0; j<3; j++){
            fields[namevf[ifield][j]] = ifield*3+nsf+j;
        }
    }

    pod.reconstructFields(gfield0, meshr, fields, &targetCells);

    if (rank==0){
        std::cout<< " " << std::endl;
        std::cout<< ">> Reconstruction coeffs hybrid interface:" << std::endl;
        for (std::size_t i = 0; i < nf; ++i)
            std::cout<< pod.getReconstructionCoeffs()[i] << std::endl;

        std::cout<< " " << std::endl;
        std::cout<< ">> Reconstruction error hybrid interface" << std::endl;
        double maxerr = 0.0;
        for (Cell & cell : meshr->getCells()){
            long id = cell.getId();
            for (std::size_t i = 0; i < nsf; ++i)
                maxerr = std::max(maxerr, std::abs(gfield0.at(id, i) - fieldr.scalar->at(id, i)));
            for (std::size_t i = 0; i < nvf; ++i){
                for (std::size_t j=0; j<3; j++)
                    maxerr = std::max(maxerr, std::abs(gfield0.at(id, nsf+3*i+j) - fieldr.vector->at(id, i)[j]));
            }
        }
        std::cout<< ">> Max error :" << maxerr << std::endl;

    }


    /**<Test hybrid interface using optimad solver format on an adapted mesh.*/
    /**<Adapt the patch with random markers.*/
    long nCells = meshr->getCellCount();
    log::cout() << std::endl;
    log::cout() << ">> Marking the cells to adapt... " << std::endl;

    for (int i = 0; i < 100; ++i) {
        long cellId = rand() % nCells * 2;
        if (!meshr->getCells().exists(cellId)) {
            continue;
        }

        for (auto neighId : meshr->findCellNeighs(cellId)) {
            meshr->markCellForRefinement(neighId);
        }
    }

    for (int i = 0; i < 50; ++i) {
        long cellId = rand() % nCells * 2;
        if (!meshr->getCells().exists(cellId)) {
            continue;
        }

        if (fieldr.mask->at(cellId)){
            meshr->markCellForCoarsening(cellId);
            for (auto neighId : meshr->findCellNeighs(cellId)) {
                if (fieldr.mask->at(neighId))
                    meshr->markCellForCoarsening(neighId);
            }
        }
    }

    log::cout() << std::endl;
    log::cout() << ">> Initial number of cells... " << nCells << std::endl;

    /**<Preadapt and adapt to update data.*/
    {

        std::vector<adaption::Info> preadaptInfo = meshr->adaptionPrepare(true);

        PiercedVector<std::vector<double>> oldData(meshr->getCellCount());
        for (adaption::Info & info : preadaptInfo){
            for (long & id : info.previous){
                oldData.insert(id, std::vector<double>(nsf+3*nvf+1, 0.0));
                for (std::size_t i=0; i<nsf+3*nvf+1; i++)
                    oldData[id][i] = gfield0.at(id,i);
            }
        }

        std::vector<adaption::Info> adaptInfo = meshr->adaptionAlter(true);

        nCells = meshr->getCellCount();
        log::cout() << ">> Final number of cells... " << nCells << std::endl;

        for (adaption::Info & info : adaptInfo){
            if (info.type == adaption::Type::TYPE_RENUMBERING){
                for (std::size_t i=0; i<nsf+3*nvf+1; i++){
                    gfield0.at(info.current[0],i) = oldData[info.previous[0]][i];
                    gfieldAdapt.at(info.current[0],i) = oldData[info.previous[0]][i];
                }
            }
            if (info.type == adaption::Type::TYPE_REFINEMENT){
                for (long & id : info.current){
                    for (std::size_t i=0; i<nsf+3*nvf+1; i++){
                        gfield0.at(id,i) = oldData[info.previous[0]][i];
                        gfieldAdapt.at(id,i) = oldData[info.previous[0]][i];
                    }
                }
            }
            if (info.type == adaption::Type::TYPE_COARSENING){
                for (std::size_t i=0; i<nsf+3*nvf+1; i++){
                    gfield0.at(info.current[0], i) = 0.0;
                    gfieldAdapt.at(info.current[0], i) = 0.0;
                }
                for (long & id : info.previous){
                    for (std::size_t i=0; i<nsf+3*nvf; i++){
                        gfield0.at(info.current[0],i) += oldData[id][i] / info.previous.size();
                        gfieldAdapt.at(info.current[0],i) += oldData[id][i] / info.previous.size();
                    }
                    gfield0.at(info.current[0], nsf+3*nvf) = std::max(gfield0.at(info.current[0], nsf+3*nvf), oldData[id][nsf+3*nvf]);
                    gfieldAdapt.at(info.current[0], nsf+3*nvf) = std::max(gfieldAdapt.at(info.current[0], nsf+3*nvf), oldData[id][nsf+3*nvf]);
                }
            }
        }

        meshr->adaptionCleanup();

    }

    /**<Reconstruct fields with pre-computed POD and evaluate maximum error.
     * The reconstruction overwrite the fields value- > use a gfieldAdapt
     * with original data to evaluate errors.*/
    targetCells.clear();
    for (Cell & cell : meshr->getCells()){
        long id = cell.getId();
        if (gfield0.at(id, nsf+3*nvf))
            targetCells.insert(id);
    }

    /**<Test pod reconstruction by using dynamic mesh mode.
     * The field is now defined on a different mesh of POD basis.*/
    pod.setStaticMesh(false);
    pod.reconstructFields(gfield0, meshr, fields, &targetCells);

    if (rank==0){
        std::cout<< " " << std::endl;
        std::cout<< ">> Reconstruction coeffs hybrid interface:" << std::endl;
        for (std::size_t i = 0; i < nf; ++i)
            std::cout<< pod.getReconstructionCoeffs()[i] << std::endl;

        std::cout<< " " << std::endl;
        std::cout<< ">> Reconstruction error hybrid interface" << std::endl;
        double maxerr = 0.0;
        for (Cell & cell : meshr->getCells()){
            long id = cell.getId();
            for (std::size_t i = 0; i < nsf+3*nvf; ++i)
                maxerr = std::max(maxerr, std::abs(gfield0.at(id, i) - gfieldAdapt.at(id, i)));
        }
        std::cout<< ">> Max error :" << maxerr << std::endl;

    }

}

/**
 * Main program.
 */

int main(int argc, char *argv[]) 
{
#if BITPIT_ENABLE_MPI
    MPI_Init(&argc,&argv);
#endif 

    int nProcs;
    int rank;

#if BITPIT_ENABLE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    nProcs = 1;
    rank   = 0;   
#endif 

    /** Run the example */
    try {
        run(rank, nProcs);
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI
    MPI_Finalize();
#endif

}
