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
#endif

using namespace std;
using namespace bitpit;

// =================================================================================== //
/*!
	\example PABLO_example_00006.cpp

	\brief Parallel 2D smoothing data using PABLO

	The example is the parallel version of PABLO_example_00003.

	In order to perform the smoothing procedure in parallel,
	ghost elements and data communications towards them are needed.
	Now ghost quadrants and ghost data are not only instantiated, as in example 00003,
	but actually used to perform smoothing across processor borders.

	<b>To run</b>: ./PABLO_example_00006 \n

	<b>To see the result visit</b>: <a href="http://optimad.github.io/PABLO/">PABLO website</a> \n

*/
// =================================================================================== //

int main(int argc, char *argv[]) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		int iter = 0;
		int dim = 2;

		/**<Instantation of a 2D pablo uniform object.*/
		PabloUniform pablo6;

		/**<Refine globally four level and write the octree.*/
		for (iter=1; iter<5; iter++){
			pablo6.adaptGlobalRefine();
		}

#if ENABLE_MPI==1
		/**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
		pablo6.loadBalance();
#endif

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo6.getNumOctants();
		uint32_t nghosts = pablo6.getNumGhosts();
		vector<double> oct_data(nocts, 0.0), ghost_data(nghosts, 0.0);

		/**<Assign a data to the octants with at least one node inside the circle.*/
		for (int i=0; i<nocts; i++){
			vector<array<double,3> > nodes = pablo6.getNodes(i);
			for (int j=0; j<4; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					oct_data[i] = 1.0;
				}
			}
		}

		/**<Assign a data to the ghost octants (PARALLEL TEST) with at least one node inside the circle.*/
		for (int i=0; i<nghosts; i++){
			/**<Compute the nodes of the octant (Use pointer for ghost).*/
			Octant *oct = pablo6.getGhostOctant(i);
			vector<array<double,3> > nodes = pablo6.getNodes(oct);
			for (int j=0; j<4; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					ghost_data[i] = 1.0;
				}
			}
		}

		/**<Update the connectivity and write the para_tree.*/
		iter = 0;
		pablo6.updateConnectivity();
		pablo6.writeTest("pablo00006_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

		/**<Smoothing iterations on initial data*/
		int start = iter + 1;
		for (iter=start; iter<start+25; iter++){
			vector<double> oct_data_smooth(nocts, 0.0);
			vector<uint32_t> neigh, neigh_t;
			vector<bool> isghost, isghost_t;
			uint8_t iface, nfaces, codim;
			for (int i=0; i<nocts; i++){
				neigh.clear();
				isghost.clear();

				/**<Find neighbours through faces (codim=1) and edges (codim=2) of the octants*/
				for (codim=1; codim<dim+1; codim++){
					if (codim == 1){
						nfaces = 4;
					}
					else if (codim == 2){
						nfaces = 4;
					}
					for (iface=0; iface<nfaces; iface++){
						pablo6.findNeighbours(i,iface,codim,neigh_t,isghost_t);
						neigh.insert(neigh.end(), neigh_t.begin(), neigh_t.end());
						isghost.insert(isghost.end(), isghost_t.begin(), isghost_t.end());
					}
				}

				/**<Smoothing data with the average over the one ring neighbours of octants*/
				oct_data_smooth[i] = oct_data[i]/(neigh.size()+1);
				for (int j=0; j<neigh.size(); j++){
					if (isghost[j]){
						oct_data_smooth[i] += ghost_data[neigh[j]]/(neigh.size()+1);
					}
					else{
						oct_data_smooth[i] += oct_data[neigh[j]]/(neigh.size()+1);
					}
				}
			}

			/**<Update the connectivity and write the para_tree.*/
			pablo6.updateConnectivity();
			pablo6.writeTest("pablo00006_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data_smooth);

#if ENABLE_MPI==1
			/**<Communicate the data of the octants and the ghost octants between the processes.*/
			UserDataComm<vector<double> > data_comm(oct_data_smooth, ghost_data);
			pablo6.communicate(data_comm);

#endif
			oct_data = oct_data_smooth;

		}
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}


