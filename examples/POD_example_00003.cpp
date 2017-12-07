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
 * \brief leave-1-out computation using voloctree.
 * This example computes the leave-1-out cross-validation starting from a 
 * database of simulations defined on the same mesh and evaluate the error maps
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

    /**<Remove snapshots from the leave-1-out method.*/   
    for (int i=0; i<5; i++)
        pod.removeLeave1outSnapshot("./data", "test."+to_string(2*i)); 
    
    /**<Compute the error map through the leave-1-out method.*/ 
    pod.leave1out();
}

/**
 * Main program.
 */

int main(int argc, char *argv[]) 
{
#if BITPIT_ENABLE_MPI
    MPI_Init(&argc,&argv);
#endif    

    // Run the example
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
