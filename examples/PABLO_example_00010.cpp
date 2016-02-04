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
	\example PABLO_example_00010.cpp

	\brief Parallel 3D dynamic adaptive mesh refinement (AMR) with data using PABLO

	In this example an Adaptive Mesh Refinement (AMR) around a moving
	sphere is performed.

	<b>To run</b>: ./PABLO_example_00010 \n

	<b>To see the result visit</b>: <a href="http://optimad.github.io/PABLO/">PABLO website</a> \n

*/
// =================================================================================== //

int main(int argc, char *argv[]) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		int iter = 0;

		/**<Instantation of a 3D para_tree object.*/
		PabloUniform pablo10(3);

		/**<Set NO 2:1 balance for the octree.*/
		int idx = 0;
		pablo10.setBalance(idx,false);

		/**<Refine globally five level and write the para_tree.*/
		for (iter=1; iter<6; iter++){
			pablo10.adaptGlobalRefine();
		}

#if ENABLE_MPI==1
		/**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
		pablo10.loadBalance();
#endif

		/**<Define the center point and the radius of the sphere.*/
		double xc, yc, zc;
		double radius = 0.15;
		double t0 = 0;
		double Dt = 0.000025;
		double omega = 2.0*3.14/0.001;

		/**<Define the center point and the trajectory of the sphere.*/
		xc = 0.25*cos(omega* Dt) + 0.5 ;
		yc = 0.25*sin(omega* Dt) + 0.5;
		zc = 100*Dt;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo10.getNumOctants();
		vector<double> oct_data(nocts, 0.0);

		/**<Adapt itend times eith data re-computing at each iteration.*/
		int itstart = 1;
		int itend = 460;

		for (iter=itstart; iter<itend; iter++){

			/**<Update the position of the sphere on the trajectory.*/
			xc = 0.25*cos(omega* Dt* (double) iter) + 0.5 ;
			yc = 0.25*sin(omega* Dt* (double) iter) + 0.5;
			zc = 100*Dt*iter;

			for (int i=0; i<nocts; i++){
				bool inside = false;
				/**<Compute the nodes of the octant.*/
				vector<array<double,3> > nodes = pablo10.getNodes(i);
				/**<Compute the center of the octant.*/
				array<double,3> center = pablo10.getCenter(i);
				oct_data[i] = sqrt((pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0)+pow((center[2]-zc),2.0)));
				for (int j=0; j<8; j++){
					double x = nodes[j][0];
					double y = nodes[j][1];
					double z = nodes[j][2];
					if (pow((x-xc),2.0)+pow((y-yc),2.0)+pow((z-zc),2.0) <= pow(radius,2.0)){
						if (pablo10.getLevel(i) < 6){
							/**<Set to refine inside the sphere.*/
							pablo10.setMarker(i,1);
						}
						else{
							pablo10.setMarker(i,0);
						}
						inside = true;
					}
				}
				if (pablo10.getLevel(i) > 5 && !inside){
					/**<Set to coarse if the octant has a level higher than 5.*/
					pablo10.setMarker(i,-1);
				}
			}

			/**<Adapt the octree.*/
			bool adapt = pablo10.adapt();

#if ENABLE_MPI==1
			/**<PARALLEL TEST: (Load)Balance the octree over the processes with communicating the data.*/
			pablo10.loadBalance();
#endif

			/**<Re-Assign to the new octants the data after an adaption.*/
			nocts = pablo10.getNumOctants();
			vector<double> oct_data_new(nocts, 0.0);
			for (int i=0; i<nocts; i++){
				array<double,3> center = pablo10.getCenter(i);
				oct_data_new[i] = sqrt((pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0)+pow((center[2]-zc),2.0)));
			}

			oct_data.resize(nocts);
			oct_data = oct_data_new;
			oct_data_new.clear();

			/**<Update the connectivity and write the para_tree.*/
			pablo10.updateConnectivity();
			pablo10.writeTest("pablo00010_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);
		}
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
