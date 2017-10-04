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

#include <mpi.h>

#include "bitpit_common.hpp"
#include "bitpit_PABLO.hpp"
#include "bitpit_IO.hpp"

using namespace std;
using namespace bitpit;

/*!
* Subtest 001
*
* Testing basic features of a 2D octree.
*/
int subtest_001()
{
    /**<Instantation of a 2D para_tree object.*/
    ParaTree pablo(2);

    /**<Set NO 2:1 balance for the octree.*/
    uint32_t idx=0;
    pablo.setBalance(idx,false);

    /**<Compute the connectivity and write the para_tree.*/
    pablo.computeConnectivity();
    pablo.write("PabloParallel001_iter0");

    /**<Refine globally two level and write the para_tree.*/
    for (int iter=1; iter<3; iter++){
        pablo.adaptGlobalRefine();
        pablo.updateConnectivity();
        pablo.write("PabloParallel001_iter"+to_string(static_cast<unsigned long long>(iter)));
    }

    /**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
    pablo.loadBalance();

    /**<Define a center point and a radius.*/
    double xc, yc;
    xc = yc = 0.5;
    double radius = 0.4;

    /**<Simple adapt() (refine) 6 times the octants with at least one node inside the circle.*/
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

        /**<(Load)Balance the octree over the processes.*/
        pablo.loadBalance();

        /**<Update the connectivity and write the para_tree.*/
        pablo.updateConnectivity();
        pablo.write("PabloParallel001_iter"+to_string(static_cast<unsigned long long>(iter)));
    }

    return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);

	// Initialize the logger
	int nProcs;
	int	rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	log::manager().initialize(log::COMBINED, true, nProcs, rank);
	log::cout().setVisibility(log::GLOBAL);

	// Run the subtests
    log::cout() << "Testing octree basic features" << std::endl;

	int status;
	try {
		status = subtest_001();
		if (status != 0) {
			return status;
		}
	} catch (const std::exception &exception) {
		log::cout() << exception.what();
		exit(1);
	}

	MPI_Finalize();
}
