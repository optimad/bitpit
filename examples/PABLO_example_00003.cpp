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
	\example PABLO_example_00003.cpp

	\brief 2D smoothing data using PABLO

	This example shows how to use PABLO's methods to find neighboring quadrants of a
	specified element.
	Neighbors search can be performed through faces, edges and nodes
	(in 2D case, only faces and nodes).

	In this example a 2D octree is refined four times and then a set of data is
	assigned to the mesh using STL vectors.
	More specifically, an integer equal to 1 is assigned uniformly to all
	quadrants within a circle, and 0 to the remaining quadrants.
	In this example neighbor-search is used to perform a simple moving-average-smoothing
	procedure of the data.

	The obtained results show the time-evolution of data over 25 smoothing iterations.

	<b>To run</b>: ./PABLO_example_00003 \n

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
		PabloUniform pablo3;

		/**<Refine globally four level and write the octree.*/
		for (iter=1; iter<5; iter++){
			pablo3.adaptGlobalRefine();
		}

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo3.getNumOctants();
		vector<double> oct_data(nocts, 0.0);

		/**<Assign a data to the octants with at least one node inside the circle.*/
		for (unsigned int i=0; i<nocts; i++){
			/**<Compute the nodes of the octant.*/
			vector<array<double,3> > nodes = pablo3.getNodes(i);
			for (int j=0; j<4; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					oct_data[i] = 1.0;
				}
			}
		}

		/**<Update the connectivity and write the octree.*/
		iter = 0;
		pablo3.updateConnectivity();
		pablo3.writeTest("pablo00003_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

		/**<Smoothing iterations on initial data*/
		int start = 1;
		for (iter=start; iter<start+25; iter++){
			vector<double> oct_data_smooth(nocts, 0.0);
			vector<uint32_t> neigh, neigh_t;
			vector<bool> isghost, isghost_t;
			uint8_t iface, nfaces;
			int codim;
			for (unsigned int i=0; i<nocts; i++){
				neigh.clear();
				isghost.clear();

				/**<Find neighbours through faces (codim=1) and edges (codim=2) of the octants*/
				for (codim=1; codim<3; codim++){
					if (codim == 1){
						nfaces = 4;
					}
					else if (codim == 2){
						nfaces = 4;
					}
					for (iface=0; iface<nfaces; iface++){
						pablo3.findNeighbours(i,iface,codim,neigh_t,isghost_t);
						neigh.insert(neigh.end(), neigh_t.begin(), neigh_t.end());
						isghost.insert(isghost.end(), isghost_t.begin(), isghost_t.end());
					}
				}

				/**<Smoothing data with the average over the one ring neighbours of octants*/
				oct_data_smooth[i] = oct_data[i]/(neigh.size()+1);
				for (unsigned int j=0; j<neigh.size(); j++){
					if (isghost[j]){
						/**< Do nothing - No ghosts: is a serial test.*/
					}
					else{
						oct_data_smooth[i] += oct_data[neigh[j]]/(neigh.size()+1);
					}
				}
			}

			/**<Update the connectivity and write the octree.*/
			pablo3.updateConnectivity();
			pablo3.writeTest("pablo00003_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data_smooth);

			oct_data = oct_data_smooth;
		}
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
