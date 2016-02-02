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

#include "PabloUniform.hpp"
#if ENABLE_MPI==1
#include "PABLO_userDataComm.hpp"
#include "PABLO_userDataLB.hpp"
#endif

using namespace std;

// =================================================================================== //

int main(int argc, char *argv[]) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		int iter = 0;

		/**<Instantation of a 3D para_tree object.*/
		PabloUniform pablo116(0,0,0,2,3);

		/**<Set NO 2:1 balance for the octree.*/
		int idx = 0;
		pablo116.setBalance(idx,false);

		/**<Refine globally five level and write the para_tree.*/
		for (iter=1; iter<6; iter++){
			pablo116.adaptGlobalRefine();
		}

#if ENABLE_MPI==1
		/**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
		pablo116.loadBalance();
#endif

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo116.getNumOctants();
		uint32_t nghosts = pablo116.getNumGhosts();
		vector<double> oct_data(nocts, 0.0), ghost_data(nghosts, 0.0);

		/**<Assign a data (distance from center of a cylinder) to the octants with at least one node inside the cylinder.*/
		for (int i=0; i<nocts; i++){
			/**<Compute the nodes of the octant.*/
			vector<array<double,3> > nodes = pablo116.getNodes(i);
			/**<Compute the center of the octant.*/
			array<double,3> center = pablo116.getCenter(i);
			for (int j=0; j<8; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					oct_data[i] = (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0));
				}
			}
		}

		/**<Update the connectivity and write the para_tree.*/
		iter = 0;
		pablo116.updateConnectivity();
		pablo116.writeTest("Pablo116_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

		/**<Adapt two times with data injection on new octants.*/
		int start = 1;
		for (iter=start; iter<start+2; iter++){
			for (int i=0; i<nocts; i++){
				/**<Compute the nodes of the octant.*/
				vector<array<double,3> > nodes = pablo116.getNodes(i);
				/**<Compute the center of the octant.*/
				array<double,3> center = pablo116.getCenter(i);
				for (int j=0; j<8; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						if (center[0]<=xc){

							/**<Set to refine to the octants in the left side of the domain inside a circle.*/
							pablo116.setMarker(i,1);
						}
						else{

							/**<Set to coarse to the octants in the right side of the domain inside a circle.*/
							pablo116.setMarker(i,-1);
						}
					}
				}
			}

			/**<Adapt the octree and map the data in the new octants.*/
			vector<double> oct_data_new;
			vector<uint32_t> mapper;
			vector<bool> isghost;
			pablo116.adapt(true);
			nocts = pablo116.getNumOctants();
			oct_data_new.resize(nocts, 0.0);

			/**<Assign to the new octant the average of the old children if it is new after a coarsening;
			 * while assign to the new octant the data of the old father if it is new after a refinement.
			 */
			for (uint32_t i=0; i<nocts; i++){
				pablo116.getMapping(i, mapper, isghost);
				if (pablo116.getIsNewC(i)){
					for (int j=0; j<8; j++){
						if (isghost[j]){
							oct_data_new[i] += ghost_data[mapper[j]]/8;
						}
						else{
							oct_data_new[i] += oct_data[mapper[j]]/8;
						}
					}
				}
				else if (pablo116.getIsNewR(i)){
					oct_data_new[i] += oct_data[mapper[0]];
				}
				else{
					oct_data_new[i] += oct_data[mapper[0]];
				}
			}

			oct_data = oct_data_new;
			vector<double>().swap(oct_data_new);

			/**<Update the connectivity and write the para_tree.*/
			pablo116.updateConnectivity();
			pablo116.writeTest("Pablo116_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

		}

#if ENABLE_MPI==1
		/**<PARALLEL TEST: (Load)Balance the octree over the processes with communicating the data.
		 * Preserve the family compact up to 4 levels over the max deep reached in the octree.*/
		uint8_t levels = 4;
		UserDataLB<vector<double> > data_lb(oct_data,ghost_data);
		pablo116.loadBalance(data_lb, levels);
#endif

		/**<Update the connectivity and write the para_tree.*/
		pablo116.updateConnectivity();
		pablo116.writeTest("Pablo116_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
