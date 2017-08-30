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

#if BITPIT_ENABLE_MPI==1
#   include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_PABLO.hpp"
#include "bitpit_IO.hpp"

using namespace std;
using namespace bitpit;

/*!
* Subtest 001
*
* Testing basic features of a 2D octree.
*
* \param rank is the rank of the process
* \param nProcs is the number of processes
*/
int subtest_001(int rank, int nProcs)
{
   /**<Instantation of a 2D para_tree object with default constructor.*/
    ParaTree ptreedefault(2);
    /**<Write the para_tree in physical domain.*/
    ptreedefault.write("Pablo001_default");

    /**<Instantation and setup of a custom (named custom) logfile.*/
    log::manager().create("custom", false, nProcs, rank);
    log::cout("custom") << fileVerbosity(log::NORMAL);
    log::cout("custom") << consoleVerbosity(log::QUIET);

    /**<Set coordinates of the origin and size of a 2D custom para_tree object.*/
    double X, Y, Z, L;
    X = 10.0; Y = 20.0; Z = 0.0; L = 250.0;
    int dim;
    dim = 2;
    /**<Instantation of a 2D para_tree object with custom constructor.*/
    PabloUniform ptreecustom(X,Y,Z,L,dim,"custom");
    /**<Write the para_tree in physical domain.*/
    ptreecustom.write("Pablo001_custom");

    return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    int nProcs;
    int rank;
#if BITPIT_ENABLE_MPI==1
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    nProcs = 1;
    rank   = 0;
#endif

    // Initialize the logger
    log::manager().initialize(log::SEPARATE, false, nProcs, rank);
    log::cout() << fileVerbosity(log::NORMAL);
    log::cout() << consoleVerbosity(log::QUIET);

    // Run the subtests
    log::cout() << "Testing basic octree features" << std::endl;

    int status;
    try {
        status = subtest_001(rank, nProcs);
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
