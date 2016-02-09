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
	\example PABLO_example_00002.cpp

	\brief 2D adaptive mesh refinement (AMR) using PABLO

	This example creates a 2D Octree mesh on the square domain [0,1]x[0,1].

	In this example the domain is iteratively refined using a slight modification
	of the first refinement criterion of PABLO_example_00001.

	Here, each octant is marked for further refinement if at least one of its vertices
	lies within a circle of specified radius. Again the iterative refinement is performed
	 for a specified number of iterations.

	 Using the present criterion the geometry is captured better than in example PABLO_example_00001.

	Then a global coarsening iteration is performed,
	followed by an iterative coarsening of each quadrant within a smaller circle inside
	the original geometry.
	The coarsening procedure is repeated 5 times and the 2:1 balancing option is
	de-activated.

	<b>To run</b>: ./PABLO_example_00002 \n

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
		PabloUniform pablo2;

		/**<Set NO 2:1 balance for ancestor octant.*/
		uint32_t idx=0;
		pablo2.setBalance(idx,true);

		/**<Compute the connectivity and write the octree.*/
		pablo2.computeConnectivity();
		pablo2.write("pablo00002_iter"+to_string(static_cast<unsigned long long>(iter)));

		/**<Refine globally two level and write the octree.*/
		for (iter=1; iter<3; iter++){
			pablo2.adaptGlobalRefine();
			pablo2.updateConnectivity();
			pablo2.write("pablo00002_iter"+to_string(static_cast<unsigned long long>(iter)));
		}

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.4;

		/**<Simple adapt() 6 times the octants with at least one node inside the circle.*/
		for (iter=3; iter<9; iter++){
			uint32_t nocts = pablo2.getNumOctants();
			for (unsigned int i=0; i<nocts; i++){
				/**<Compute the nodes of the octant.*/
				vector<array<double,3> > nodes = pablo2.getNodes(i);
				for (int j=0; j<4; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						pablo2.setMarker(i, 1);
					}
				}
			}
			/**<Adapt octree.*/
			pablo2.adapt();

			/**<Update the connectivity and write the octree.*/
			pablo2.updateConnectivity();
			pablo2.write("pablo00002_iter"+to_string(static_cast<unsigned long long>(iter)));
		}

		/**<Coarse globally one level and write the octree.*/
		pablo2.adaptGlobalCoarse();
		pablo2.updateConnectivity();
		pablo2.write("pablo00002_iter"+to_string(static_cast<unsigned long long>(iter)));


		/**<Define a center point and a radius.*/
		xc = yc = 0.35;
		radius = 0.15;

		/**<Simple adapt() 5 times the octants with at least one node inside the circle.*/
		for (iter=10; iter<15; iter++){
			uint32_t nocts = pablo2.getNumOctants();
			for (unsigned int i=0; i<nocts; i++){
				/**<Compute the nodes of the octant.*/
				vector<array<double,3> > nodes = pablo2.getNodes(i);
				for (int j=0; j<4; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					/**<Set refinement marker=-1 (coarse it one time) for octants inside a circle.*/
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						pablo2.setMarker(i, -1);
					}
				}
			}
			/**<Adapt octree, update connectivity and write.*/
			pablo2.adapt();
			pablo2.updateConnectivity();
			pablo2.write("pablo00002_iter"+to_string(static_cast<unsigned long long>(iter)));
		}
#if ENABLE_MPI==1
	}
	MPI::Finalize();
#endif
}

