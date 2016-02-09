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

#if ENABLE_MPI==1
#include "PABLO_userDataComm.hpp"
#include "PABLO_userDataLB.hpp"
#endif

using namespace std;
using namespace bitpit;

// =================================================================================== //
/*!
	\example PABLO_example_00007.cpp

	\brief Parallel 2D adaptive mesh refinement (AMR) with data using PABLO

	The example is the parallel version of PABLO_example_00004.

	Here the main focus is on the load-balance of both grid and data.
	The grid is refined several times together with the data and their inheritance
	follows the same rules like in example 00004.
	Until the last refinement no parallel paradigm is in action: every process owns
	the entire grid.

	After this refinement, the load-balance with data is introduced, giving as result
	a parallel distribution of grid and data.

	Even a constant vector of weight is used, in order to show that the loadBalance
	can be performed by using a weight function for each cell.

	<b>To run</b>: ./PABLO_example_00007 \n

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
		PabloUniform pablo7;

		/**<Set NO 2:1 balance for the octree.*/
		int idx = 0;
		pablo7.setBalance(idx,false);

		/**<Refine globally five level and write the octree.*/
		for (iter=1; iter<6; iter++){
			pablo7.adaptGlobalRefine();
		}

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo7.getNumOctants();
		uint32_t nghosts = pablo7.getNumGhosts();
		vector<double> oct_data(nocts, 0.0), ghost_data(nghosts, 0.0);

		/**<Assign a data (distance from center of a circle) to the octants with at least one node inside the circle.*/
		for (unsigned int i=0; i<nocts; i++){
			/**<Compute the nodes of the octant.*/
			vector<array<double,3> > nodes = pablo7.getNodes(i);
			/**<Compute the center of the octant.*/
			array<double,3> center = pablo7.getCenter(i);
			for (int j=0; j<4; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					oct_data[i] = (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0));
				}
			}
		}

		/**<Update the connectivity and write the octree.*/
		iter = 0;
		pablo7.updateConnectivity();
		pablo7.writeTest("pablo00007_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

		/**<Adapt two times with data injection on new octants.*/
		int start = 1;
		/**<Weight.*/
		vector<double> weight(nocts, 1.0),weightGhost;
		for (iter=start; iter<start+2; iter++){
			for (unsigned int i=0; i<nocts; i++){
				/**<Compute the nodes of the octant.*/
				vector<array<double,3> > nodes = pablo7.getNodes(i);
				/**<Compute the center of the octant.*/
				array<double,3> center = pablo7.getCenter(i);
				for (int j=0; j<4; j++){
					weight[i] = 2.0;
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						if (center[0]<=xc){

							/**<Set to refine to the octants in the left side of the domain inside a circle.*/
							pablo7.setMarker(i,1);
							weight[i] = 1.0;
						}
						else{

							/**<Set to coarse to the octants in the right side of the domain inside a circle.*/
							pablo7.setMarker(i,-1);
							weight[i] = 1.0;
						}
					}
				}
			}

			/**<Adapt the octree and map the data in the new octants.*/
			vector<double> oct_data_new;
			vector<double> weight_new;
			vector<uint32_t> mapper;
			vector<bool> isghost;
			pablo7.adapt(true);
			nocts = pablo7.getNumOctants();
			oct_data_new.resize(nocts, 0.0);
			weight_new.resize(nocts,0.0);

			/**<Assign to the new octant the average of the old children if it is new after a coarsening;
			 * while assign to the new octant the data of the old father if it is new after a refinement.
			 */
			for (uint32_t i=0; i<nocts; i++){
				pablo7.getMapping(i, mapper, isghost);
				if (pablo7.getIsNewC(i)){
					for (int j=0; j<4; j++){
						oct_data_new[i] += oct_data[mapper[j]]/4;
						weight_new[i] += weight[mapper[j]];
					}
				}
				else if (pablo7.getIsNewR(i)){
					oct_data_new[i] += oct_data[mapper[0]];
					weight_new[i] += weight[mapper[0]];
				}
				else{
					oct_data_new[i] += oct_data[mapper[0]];
					weight_new[i] += weight[mapper[0]];
				}
			}

			/**<Update the connectivity and write the octree.*/
			pablo7.updateConnectivity();
			pablo7.writeTest("pablo00007_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data_new);

			oct_data = oct_data_new;
			weight = weight_new;
		}

#if ENABLE_MPI==1
		/**<(Load)Balance the octree over the processes with communicating the data.*/
		UserDataLB<vector<double> > data_lb(weight,weightGhost);
		pablo7.loadBalance(data_lb, &weight);
#endif

		double tot = 0.0;
		for (unsigned int i=0; i<weight.size(); i++){
			tot += weight[i];
		}

		/**<Update the connectivity and write the octree.*/
		pablo7.updateConnectivity();
		pablo7.writeTest("pablo00007_iter"+to_string(static_cast<unsigned long long>(iter)), weight);

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}



