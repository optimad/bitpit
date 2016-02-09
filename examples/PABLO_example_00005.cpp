/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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

#include "bitpit_PABLO.hpp"

using namespace std;
using namespace bitpit;

// =================================================================================== //
/*!
	\example PABLO_example_00005.cpp

	\brief Parallel 2D adaptive mesh refinement (AMR) with data using PABLO

	This example is the parallel version of PABLO_example_00002.

	Here one quadtree is instantiated in the same way on every process.
	Therefore, each process refines globally the initial quadrant 3 times.
	Till this moment, there's no parallel paradigm in action, actually the code
	is replicated on every process.

	At this point a load-balance method call is performed and the grid is
	partitioned and distributed among the processes of the world communicator.
	From now on, each process owns one portion of the grid.

	The test continues as in example 00002, but after every refinement/coarsening
	a load-balance call is introduced to keep the computational burden well distributed.

	<b>To run</b>: ./PABLO_example_00005 \n

	<b>To see the result visit</b>: <a href="http://optimad.github.io/PABLO/">PABLO website</a> \n

*/
// =================================================================================== //

int main(int argc, char *argv[]) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		int iter = 0;

		/**<Instantation of a 2D pablo uniform object.*/
		PabloUniform pablo5;

		/**<Set NO 2:1 balance for the octree.*/
		uint32_t idx=0;
		pablo5.setBalance(idx,false);

		/**<Compute the connectivity and write the octree.*/
		pablo5.computeConnectivity();
		pablo5.write("pablo00005_iter"+to_string(static_cast<unsigned long long>(iter)));

		/**<Refine globally two level and write the octree.*/
		for (iter=1; iter<3; iter++){
			pablo5.adaptGlobalRefine();
			pablo5.updateConnectivity();
			pablo5.write("pablo00005_iter"+to_string(static_cast<unsigned long long>(iter)));
		}

#if ENABLE_MPI==1
		/**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
		pablo5.loadBalance();
#endif

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.4;

		/**<Simple adapt() (refine) 6 times the octants with at least one node inside the circle.*/
		for (iter=3; iter<9; iter++){
			uint32_t nocts = pablo5.getNumOctants();
			for (unsigned int i=0; i<nocts; i++){
				/**<Compute the nodes of the octant.*/
				vector<array<double,3> > nodes = pablo5.getNodes(i);
				for (int j=0; j<4; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						pablo5.setMarker(i, 1);
					}
				}
			}

			/**<Adapt octree.*/
			pablo5.adapt();

#if ENABLE_MPI==1
			/**<(Load)Balance the octree over the processes.*/
			pablo5.loadBalance();
#endif

			/**<Update the connectivity and write the octree.*/
			pablo5.updateConnectivity();
			pablo5.write("pablo00005_iter"+to_string(static_cast<unsigned long long>(iter)));
		}

		/**<Coarse globally one level and write the octree.*/
		pablo5.adaptGlobalCoarse();
		pablo5.updateConnectivity();
		pablo5.write("pablo00005_iter"+to_string(static_cast<unsigned long long>(iter)));


		/**<Define a center point and a radius.*/
		xc = yc = 0.35;
		radius = 0.15;

		/**<Simple adapt() 5 times the octants with at least one node inside the circle.*/
		for (iter=10; iter<15; iter++){
			uint32_t nocts = pablo5.getNumOctants();
			for (unsigned int i=0; i<nocts; i++){
				/**<Compute the nodes of the octant.*/
				vector<array<double,3> > nodes = pablo5.getNodes(i);
				for (int j=0; j<4; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					/**<Set refinement marker=-1 (coarse it one time) for octants inside a circle.*/
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						pablo5.setMarker(i, -1);
					}
				}
			}
			/**<Adapt octree.*/
			pablo5.adapt();

#if ENABLE_MPI==1
			/**<(Load)Balance the octree over the processes.*/
			pablo5.loadBalance();
#endif

			/**<Update the connectivity and write the octree.*/
			pablo5.updateConnectivity();
			pablo5.write("pablo00005_iter"+to_string(static_cast<unsigned long long>(iter)));
		}

		/**<Coarse globally one level and write the octree.*/
		pablo5.adaptGlobalCoarse();
		pablo5.updateConnectivity();
		pablo5.write("pablo00005_iter"+to_string(static_cast<unsigned long long>(iter)));

#if ENABLE_MPI==1
	}
	MPI::Finalize();
#endif
}
