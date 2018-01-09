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

/**
 * \example POD_example_reconstruction_00002.cpp
 * 
 * \brief POD basis computation using voloctree.
 * This example computes the POD basis starting from a database of simulations
 * defined on the same mesh and evaluate the reconstruction of a snapshot
 * not included in the database testing the hybrid interface.
 * <b>To run</b>: ./voloctree_example_00002 \n
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
    for (int i=0; i<11; i++)
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
        std::cout<< ">> Reconstruction coeffs:" << std::endl;
        for (std::size_t i = 0; i < nf; ++i)
            std::cout<< pod.getReconstructionCoeffs()[i] << std::endl;
    }

    /**<Test hybrid interface using optimad solver format.*/
    PiercedStorage<double> gfield0(nsf+3*nvf, &meshr->getCells());
    std::map<std::string, std::size_t> scalarFields;
    std::map<std::array<std::string, 3>, std::array<std::size_t,3>> vectorFields;
    std::unordered_set<long> targetCells;

    std::map<std::string, std::size_t> fields;

    for (Cell & cell : meshr->getCells()){
        long id = cell.getId();
        for (std::size_t ifield=0; ifield<nsf; ifield++)
            gfield0.at(id, ifield) = fieldr.scalar->at(id, ifield);

        for (std::size_t ifield=0; ifield<nvf; ifield++){
            for (std::size_t j=0; j<3; j++)
                gfield0.at(id, ifield*3+nsf+j) = fieldr.vector->at(id, ifield)[j];

        }
        if (fieldr.mask->at(id))
            targetCells.insert(id);
    }

    for (std::size_t ifield=0; ifield<nsf; ifield++){
        scalarFields[namesf[ifield]] = ifield;
        fields[namesf[ifield]] = ifield;
    }
    for (std::size_t ifield=0; ifield<nvf; ifield++){
        for (std::size_t j=0; j<3; j++){
            vectorFields[namevf[ifield]][j] = ifield*3+nsf+j;
            fields[namevf[ifield][j]] = ifield*3+nsf+j;
        }
    }

    pod.reconstructFields(gfield0, meshr, scalarFields, vectorFields, &targetCells);

    if (rank==0){
        std::cout<< " " << std::endl;
        std::cout<< ">> Reconstruction coeffs hybrid interface:" << std::endl;
        for (std::size_t i = 0; i < nf; ++i)
            std::cout<< pod.getReconstructionCoeffs()[i] << std::endl;
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

    // Run the example
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
