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
 * \example POD_example_00003.cpp
 * 
 * \brief POD leave-1-out error map computation using voloctree.
 * This example uses the leave-1-out cross-validation method to compute the reconstruction
 * error map starting from a database of simulations defined on the same mesh.
 * It evaluates also the bounding box containing all those cells whose error is equal or greater
 * than an assigned threshold.
 * <b>To run</b>: ./POD_example_00003 \n
 */ 

#include <array>
#if BITPIT_ENABLE_MPI
#include <mpi.h>
#endif

#include "pod.hpp"

using namespace bitpit;

/**
 * Run the example.
 */ 
void run()
{
    /**<Create POD object.*/
    POD pod;

    /**<Add snapshots to database.*/   
    for (int i=0; i<10; i++)
        pod.addSnapshot("./data", "test."+to_string(i));

    /**<Set POD.*/    
    pod.setMeshType(POD::MeshType::VOLOCTREE);
    pod.setStaticMesh(true); 
    pod.setUseMean(false);
    pod.setWriteMode(POD::WriteMode::DEBUG);
    pod.setMemoryMode(POD::MemoryMode::MEMORY_NORMAL);
    pod.setEnergyLevel(99);
    
    pod.setDirectory("pod");
    pod.setName("pod.test.solver");    

    /**<Remove snapshots from the leave-1-out method.
     * These snapshots are always used in the POD bases computation and
     * the corresponding reconstruction error is never evaluated. */   
    for (int i=0; i<5; i++)
        pod.removeLeave1outSnapshot("./data", "test."+to_string(2*i)); 
    
    /**<Compute the error map through the leave-1-out method.*/ 
    pod.leave1out();
    
    /**<Set target error fields used in the bounding box evaluation.*/
    std::vector<std::string> namesf {"p","a"};
    std::vector<std::array<std::string,3>> namevf {}; //{{"u_x", "u_y","u_z"}};  
    pod.setTargetErrorFields(namesf,namevf);

    /**<Set error threshold for the bounding box evaluation.*/
    pod.setErrorThreshold(0.001);
    
    /**<Evaluate the bounding box of the target error fields.*/
    pod.evalErrorBoundingBox();
}

/**
 * Main program.
 */

int main(int argc, char *argv[]) 
{
#if BITPIT_ENABLE_MPI
    MPI_Init(&argc,&argv);
#endif    

    /** Run the example */
    try {
        run();
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI
    MPI_Finalize();
#endif

}
