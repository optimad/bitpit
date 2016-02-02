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

int main(int argc, char *argv[]) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		int iter = 0;
		/**<Instantation of a 2D para_tree object.*/
		ParaTree pablo17a;
		ParaTree pablo17b;

		/**<Set NO 2:1 balance for the octree.*/
		int idx = 0;
		pablo17a.setBalance(idx,false);

		/**<Refine globally five level and write the para_tree.*/
		for (iter=1; iter<6; iter++){
			pablo17a.adaptGlobalRefine();
			pablo17b.adaptGlobalRefine();
		}

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo17a.getNumOctants();
		uint32_t nghosts = pablo17a.getNumGhosts();
		vector<double> oct_data(nocts, 0.0);
		vector<double> ghost_data(nghosts, 0.0);

		/**<Assign a data (distance from center of a circle) to the octants with at least one node inside the circle.*/
		for (int i=0; i<nocts; i++){
			/**<Compute the center of the octant.*/
			array<double,3> center = pablo17a.getCenter(i);
			oct_data[i] = (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0));
			/**<Compute the nodes of the octant.*/
			vector<array<double,3> > nodes = pablo17a.getNodes(i);
			for (int j=0; j<4; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					if (center[0]<=xc){

						/**<Set to refine to the octants in the left side of the domain.*/
						pablo17a.setMarker(i,1);
					}
					else{

						/**<Set to coarse to the octants in the right side of the domain.*/
						pablo17a.setMarker(i,-1);
					}
				}
			}
		}

		/**<Update the connectivity and write the para_tree.*/
		iter = 0;
		pablo17a.updateConnectivity();
		pablo17b.updateConnectivity();
		pablo17a.write("Pablo17a_iter"+to_string(static_cast<unsigned long long>(iter)));

		/**<Adapt two times with data injection on new octants.*/
		int start = 1;
		for (iter=start; iter<start+2; iter++){
			for (int i=0; i<nocts; i++){
				/**<Compute the center of the octant.*/
				array<double,3> center = pablo17a.getCenter(i);
				/**<Compute the nodes of the octant.*/
				vector<array<double,3> > nodes = pablo17a.getNodes(i);
				for (int j=0; j<4; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						if (center[0]<xc){

							/**<Set to refine to the octants in the left side of the domain inside a circle.*/
							pablo17a.setMarker(i,1);
						}
						else{

							/**<Set to coarse to the octants in the right side of the domain inside a circle.*/
							pablo17a.setMarker(i,-1);
						}
					}
				}
			}

			/**<Adapt the octree and map the data in the new octants.*/
			vector<double> oct_data_new;
			vector<uint32_t> mapper;
			vector<bool> isghost;
			pablo17a.adapt(true);
			nocts = pablo17a.getNumOctants();
			oct_data_new.resize(nocts, 0.0);

			/**<Assign to the new octant the average of the old children if it is new after a coarsening;
			 * while assign to the new octant the data of the old father if it is new after a refinement.
			 */
			for (uint32_t i=0; i<nocts; i++){
				pablo17a.getMapping(i, mapper, isghost);
				if (pablo17a.getIsNewC(i)){
					for (int j=0; j<4; j++){
						if (isghost[j]){
							oct_data_new[i] += ghost_data[mapper[j]]/4;
						}
						else{
							oct_data_new[i] += oct_data[mapper[j]]/4;
						}
					}
				}
				else if (pablo17a.getIsNewR(i)){
					oct_data_new[i] += oct_data[mapper[0]];
				}
				else{
					oct_data_new[i] += oct_data[mapper[0]];
				}
			}
			oct_data = oct_data_new;
			oct_data_new.clear();

			/**<Update the connectivity and write the para_tree.*/
			pablo17a.updateConnectivity();
			pablo17a.writeTest("Pablo17a_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

		}

#if ENABLE_MPI==1
		/**<PARALLEL TEST: (Load)Balance the octree over the processes with communicating the data.
		 * Preserve the family compact up to 4 levels over the max deep reached in the octree.*/
		UserDataLB<vector<double> > data_lb(oct_data,ghost_data);
		pablo17a.loadBalance(data_lb);
		pablo17b.loadBalance();

		/**<Update the connectivity and write the para_tree.*/
		pablo17a.updateConnectivity();
		pablo17a.writeTest("Pablo17a_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

#endif


#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
