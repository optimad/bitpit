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

#include <array>
#if BITPIT_ENABLE_MPI
#include <mpi.h>
#endif

#include "bitpit_voloctree.hpp"
#include "pod.hpp"

using namespace bitpit;

/** Subtest 001
 * 
 * Testing basic POD features on a 2D octree patch
 * 
 * \param rank is the rank of the process
 * \param nProcs is the number of the processes
 */
int subtest_001(int rank, int nProcs)
{
    if (rank == 0)
        log::cout() << "Creating 2D patch..." << std::endl;

    /**<Set coordinates of the origin and size of the 2D octree patch.*/
    std::array<double, 3> origin = {{0., 0., 0.}};
    double length = 2*BITPIT_PI;
    double dh = length/100;

    /**<Create the patch.*/ 
    VolumeKernel * mesh = new VolOctree(2, origin, length, dh);
#if BITPIT_ENABLE_MPI    
    mesh->setCommunicator(MPI_COMM_WORLD);
#endif
    mesh->buildAdjacencies();
    mesh->buildInterfaces();
    mesh->update();

#if BITPIT_ENABLE_MPI
    mesh->partition(true);
#endif 

    int archiveVersion = 1;
    int dumpBlock = (nProcs > 1) ? rank : -1;

    /**<Create synthetic fields.*/
    if (rank == 0)
        log::cout() << ">> Creating synthetic fields... " << std::endl;
    bitpit::PiercedStorage<double> fields(1, &mesh->getCells());
    bitpit::PiercedStorage<std::array<double, 3>> fieldv(1, &mesh->getCells());
    bitpit::PiercedStorage<bool> mask(1, &mesh->getCells());  
    mask.fill(true);

    for (int k=0; k<2; k++){   
        for (bitpit::Cell & cell : mesh->getCells()){
            long id = cell.getId();
            if (k%2){
                fields.at(id) = std::cos((k)*mesh->evalCellCentroid(id)[0]);
                fieldv.at(id) = {{1.,0.,0.}};
            }else{
                fields.at(id) = std::sin((k+1)*mesh->evalCellCentroid(id)[1]);
                fieldv.at(id) = {{0.,1.,0.}};
            }
        }
        /**<Dump the snapshots.*/
        {
            if (rank == 0)
                log::cout() << "Dumping snapshot " << k << std::endl;
            std::string header = "octree snapshot";
            OBinaryArchive binaryWriter2D("snapshot."+std::to_string(k)+".data", archiveVersion, header, dumpBlock);
            std::ostream &dataStream = binaryWriter2D.getStream();
            mask.dump(dataStream);
            utils::binary::write(dataStream,std::size_t(1));
            utils::binary::write(dataStream,std::string("scalar"));
            fields.dump(dataStream);
            std::array<std::string,3> namevf = {{"vector_x","vector_y","vector_z"}};
            utils::binary::write(dataStream,std::size_t(1));
            utils::binary::write(dataStream,namevf);
            fieldv.dump(dataStream);
            binaryWriter2D.close();
        }

        /**<Dump the mesh.*/
        {
            if (rank == 0)
                log::cout() << "Dumping mesh " << k << std::endl;
            std::string header = "octree patch";
            OBinaryArchive binaryWriter2D("snapshot."+std::to_string(k)+".mesh", archiveVersion, header,dumpBlock);
            mesh->dump(binaryWriter2D.getStream());
            binaryWriter2D.close();
        }
    }

    /**<Create POD object.*/
    POD pod;

    /**<Add snapshots to database.*/   
    pod.addSnapshot(".", "snapshot.0");
    pod.addSnapshot(".", "snapshot.1");

    /**<Set POD.*/    
    pod.setMeshType(POD::MeshType::VOLOCTREE);
    pod.setStaticMesh(true);
    pod.setWriteMode(POD::WriteMode::DEBUG);
    pod.setMemoryMode(POD::MemoryMode::MEMORY_NORMAL);
    pod.setReconstructionMode(POD::ReconstructionMode::PROJECTION);
    pod.setModeCount(2);
    pod.setDirectory(".");
    pod.setName("s001_pod");

    /**<Compute the POD basis.*/ 
    pod.evalMeanMesh();
    pod.fillListActiveIDs(mask);
    pod.evalCorrelation();
    pod.evalEigen();
    pod.evalModes();
    pod.setSensorMask(mask);

    /**<Reconstruct first snapshot.*/    
    pod.addReconstructionSnapshot(".", "snapshot.0");
    pod.evalReconstruction();
    std::vector<std::vector<double> > rcoeffs0 = pod.getReconstructionCoeffs();

    /**<Reconstruct second snapshot.*/  
    pod.addReconstructionSnapshot(".", "snapshot.1");
    pod.evalReconstruction();
    std::vector<std::vector<double> > rcoeffs1 = pod.getReconstructionCoeffs();

    std::vector<std::vector<double> > sum = rcoeffs0+rcoeffs1;

    if (rank==0){ 
        std::cout<< ">> Reconstruction coeffs:" << std::endl;
        std::cout<< rcoeffs0 << endl;
        std::cout<< rcoeffs1 << endl;
    }

    for (int i=0; i<2; i++)
        if (bitpit::abs(*sum[i].data()) > 1e-8 ){
            if (rank==0)
                std::cout<< "\ntest failed" << std::endl;
            return 1;
        }
    return 0;
}

/*!
 * Main program.
 */
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
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

    int status = 1;
    try {
        status = subtest_001(rank, nProcs);

    } catch (const std::exception &exception) {
        log::cout() << "test_podvoloctree_00001 exited with an error of type :" << exception.what() << std::endl;
        exit(1);
    }

#if BITPIT_ENABLE_MPI
    MPI_Finalize();
#endif 

    return status;
}

