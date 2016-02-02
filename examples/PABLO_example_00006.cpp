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

#include "ParaTree.hpp"

using namespace std;
using namespace bitpit;

// =================================================================================== //

int main(int argc, char *argv[]) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		int iter = 0;

		/**<Instantation of a 2D para_tree object.*/
		ParaTree pablo13;

		/**<Set NO 2:1 balance for the octree.*/
		uint32_t idx=0;
		pablo13.setBalance(idx,false);

		/**<Compute the connectivity and write the para_tree.*/
		pablo13.computeConnectivity();
		pablo13.write("Pablo13_iter"+to_string(static_cast<unsigned long long>(iter)));

		/**<Refine globally two level and write the para_tree.*/
		for (iter=1; iter<3; iter++){
			pablo13.adaptGlobalRefine();
			pablo13.updateConnectivity();
			pablo13.write("Pablo13_iter"+to_string(static_cast<unsigned long long>(iter)));
		}

#if ENABLE_MPI==1
		/**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
		pablo13.loadBalance();
#endif

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.4;

		/**<Simple adapt() (refine) 6 times the octants with at least one node inside the circle.*/
		for (iter=3; iter<9; iter++){
			uint32_t nocts = pablo13.getNumOctants();
			for (int i=0; i<nocts; i++){
				/**<Compute the nodes of the octant.*/
				vector<array<double,3> > nodes = pablo13.getNodes(i);
				for (int j=0; j<4; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						pablo13.setMarker(i, 1);
					}
				}
			}

			/**<Adapt octree.*/
			pablo13.adapt();

#if ENABLE_MPI==1
			/**<(Load)Balance the octree over the processes.*/
			pablo13.loadBalance();
#endif

			/**<Update the connectivity and write the para_tree.*/
			pablo13.updateConnectivity();
			pablo13.write("Pablo13_iter"+to_string(static_cast<unsigned long long>(iter)));
		}

		/**<Coarse globally one level and write the para_tree.*/
		pablo13.adaptGlobalCoarse();
		pablo13.updateConnectivity();
		pablo13.write("Pablo13_iter"+to_string(static_cast<unsigned long long>(iter)));


		/**<Define a center point and a radius.*/
		xc = yc = 0.35;
		radius = 0.15;

		/**<Simple adapt() 5 times the octants with at least one node inside the circle.*/
		for (iter=10; iter<15; iter++){
			uint32_t nocts = pablo13.getNumOctants();
			for (int i=0; i<nocts; i++){
				/**<Compute the nodes of the octant.*/
				vector<array<double,3> > nodes = pablo13.getNodes(i);
				for (int j=0; j<4; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					/**<Set refinement marker=-1 (coarse it one time) for octants inside a circle.*/
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						pablo13.setMarker(i, -1);
					}
				}
			}
			/**<Adapt octree.*/
			pablo13.adapt();

#if ENABLE_MPI==1
			/**<(Load)Balance the octree over the processes.*/
			pablo13.loadBalance();
#endif

			/**<Update the connectivity and write the para_tree.*/
			pablo13.updateConnectivity();
			pablo13.write("Pablo13_iter"+to_string(static_cast<unsigned long long>(iter)));
		}
#if ENABLE_MPI==1
	}
	MPI::Finalize();
#endif
}
