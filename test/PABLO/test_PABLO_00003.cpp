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
* Testing refinement and coarsening.
*/
int subtest_001()
{
	/**<Instantation of a 2D para_tree object.*/
	ParaTree pablo(2);

    /**<Set 2:1 balance only through faces.*/
    pablo.setBalanceCodimension(1);
    uint32_t idx=0;
    pablo.setBalance(idx,true);

    /**<Compute the connectivity and write the para_tree.*/
    pablo.computeConnectivity();
    pablo.write("Pablo003_iter0");

    /**<Refine globally two level and write the para_tree.*/
    for (int iter=1; iter<3; iter++){
        pablo.adaptGlobalRefine();
        pablo.updateConnectivity();
        pablo.write("Pablo003_iter"+to_string(static_cast<unsigned long long>(iter)));
    }

    /**<Define a center point and a radius.*/
    double xc, yc;
    xc = yc = 0.5;
    double radius = 0.4;

    /**<Simple adapt() [refine] 6 times the octants with at least one node inside the circle.*/
    int fiter;
    for (int iter=3; iter<9; iter++){
        uint32_t nocts = pablo.getNumOctants();
        for (unsigned int i=0; i<nocts; i++){
            /**<Compute the nodes of the octant.*/
            vector<array<double,3> > nodes = pablo.getNodes(i);
            for (int j=0; j<4; j++){
                double x = nodes[j][0];
                double y = nodes[j][1];
                if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
                    pablo.setMarker(i, 1);
                }
            }
        }
        /**<Adapt octree.*/
        pablo.adapt();

        /**<Update the connectivity and write the para_tree.*/
        pablo.updateConnectivity();
        pablo.write("Pablo003_iter"+to_string(static_cast<unsigned long long>(iter)));
    }

    /**<Simple adapt() [coarse] 3 times the octants with at least one node inside the 2nd circle.*/
    /**<Define a center point and a radius.*/
    double xc2, yc2;
    xc2 = yc2 = 0.5;
    double radius2 = 0.2;
    for (int iter=9; iter<12; iter++){
        uint32_t nocts = pablo.getNumOctants();
        for (unsigned int i=0; i<nocts; i++){
            /**<Compute the nodes of the octant.*/
            vector<array<double,3> > nodes = pablo.getNodes(i);
            for (int j=0; j<4; j++){
                double x = nodes[j][0];
                double y = nodes[j][1];
                if ((pow((x-xc2),2.0)+pow((y-yc2),2.0) <= pow(radius2,2.0))){
                    pablo.setMarker(i, -1);
                }
            }
            if (iter == 11) pablo.setBalance(i,false);
        }
        /**<Adapt octree.*/
        pablo.adapt();

        /**<Update the connectivity and write the para_tree.*/
        pablo.updateConnectivity();
        pablo.write("Pablo003_iter"+to_string(static_cast<unsigned long long>(iter)));
        fiter = iter;
    }


    uint32_t nocts = pablo.getNumOctants();
    for (unsigned int i=0; i<nocts; i++){
        pablo.setBalance(i,true);
    }

    /**<Adapt octree.*/
    pablo.adapt();

    /**<Update the connectivity and write the para_tree.*/
    pablo.updateConnectivity();
    pablo.write("Pablo003_iter"+to_string(static_cast<unsigned long long>(fiter+1)));

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
        status = subtest_001();
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
