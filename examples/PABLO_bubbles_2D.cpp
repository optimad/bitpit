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
	\example PABLO_bubbles_2D.cpp

	\brief 2D dynamic adaptive mesh refinement (AMR) using PABLO

	This example creates a 2D Octree mesh on the square domain [0,1]x[0,1].

	The domain is adapted to track a set of moving bubbles with random initialization.

	<b>To run</b>: ./PABLO_bubbles_2D \n

	<b>To see the result visit</b>: <a href="http://optimad.github.io/PABLO/">PABLO website</a> \n
*/


/*!  \cond  EXAMPLE_CLASSES */

/**<Declaration of a class bubble with center and radius.*/
class bubble{
public:
	double c[2];
	double r;
};

/*!  \endcond   */

int main(int argc, char *argv[]) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		int iter = 0;

		/**<Instantation of a 2D para_tree object.*/
		PabloUniform pabloBB;

		/**<Set 2:1 balance for the octree.*/
		pabloBB.setBalanceCodimension(1);
		int idx = 0;
		pabloBB.setBalance(idx,true);

		/**<Set periodic condition on face 0 (and 1) (test 2:1 balance for periodic boundaries).*/
		pabloBB.setPeriodic(0);

		/**<Refine globally four level and write the para_tree.*/
		for (iter=1; iter<4; iter++){
			pabloBB.adaptGlobalRefine();
		}

#if ENABLE_MPI==1
		/**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
		pabloBB.loadBalance();
#endif

		/**<Define and initialize a set of bubbles and their trajectories.*/
		time_t Time = time(NULL);
		srand(Time);
		if(pabloBB.getRank() == 0)
			cout << "the seed = " << Time << endl;

#if ENABLE_MPI==1
		int nb = 50;
#else
		int nb = 10;
#endif
		vector<bubble> BB;
		vector<bubble> BB0;
		vector<double> DY;
		vector<double> OM;
		vector<double> AA;

		for (int i=0; i<nb; i++){
			double randc[2];
			randc[0] = 0.8 * (double) (rand()) /  RAND_MAX + 0.1;
			randc[1] = (double) (rand()) /  RAND_MAX - 0.5;
			double randr = 0.1 * (double) (rand()) / RAND_MAX + 0.02;
			double dy = 0.005 + 0.05 * (double) (rand()) / RAND_MAX;
			double omega = 0.5 * (double) (rand()) / RAND_MAX;
			double aa = 0.15 * (double) (rand()) / RAND_MAX;
			bubble bb;
			bb.c[0] = randc[0];
			bb.c[1] = randc[1];
			bb.r = randr;
			BB.push_back(bb);
			BB0.push_back(bb);
			DY.push_back(dy);
			OM.push_back(omega);
			AA.push_back(aa);
		}
		/**<Initialize time and timestep.*/
		double t0 = 0;
		double t = t0;
		double Dt = 0.5;

		/**<Adapt itend times with refinement on the interface of the bubbles.*/
		int itstart = 1;
		int itend = 200;

		/**<Perform time iterations.*/
		for (iter=itstart; iter<itend; iter++){
			if(pabloBB.getRank()==0) cout << "iter " << iter << endl;
			t += Dt;

			/**<Update bubbles position.*/
			for (int i=0; i<nb; i++){
				BB[i].c[0] = BB0[i].c[0] + AA[i]*cos(OM[i]*t);
				BB[i].c[1] = BB[i].c[1]+ Dt*DY[i];
			}

			/**<Adapting (refinement and coarsening).*/
			bool adapt = true;
			while (adapt){
				octantIterator it, itend = pabloBB.getInternalOctantsEnd();
				for (it=pabloBB.getInternalOctantsBegin(); it!=itend; ++it){
					bool inside = false;
					/**<Compute the nodes of the octant.*/
					vector<array<double,3> > nodes = pabloBB.getNodes((*it));
					/**<Compute the center of the octant.*/
					array<double,3> center = pabloBB.getCenter((*it));
					int ib = 0;
					while (!inside && ib<nb){
						double xc = BB[ib].c[0];
						double yc = BB[ib].c[1];
						double radius = BB[ib].r;
						/**<Set marker with condition on center or nodes of the octant.*/
						for (int j=0; j<4; j++){
							double x = nodes[j][0];
							double y = nodes[j][1];
							if ( ((!inside) &&
									(pow((x-xc),2.0)+pow((y-yc),2.0) <= 1.25*pow(radius,2.0) &&
											pow((x-xc),2.0)+pow((y-yc),2.0) >= 0.75*pow(radius,2.0)))
											|| ((!inside) && (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0) <= 1.25*pow(radius,2.0) &&
													pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0) >= 0.75*pow(radius,2.0)))){
								if (int(pabloBB.getLevel((*it))) < 9){
									/**<Set to refine inside a band around the interface of the bubbles.*/
									pabloBB.setMarker((*it),1);
								}
								else{
									pabloBB.setMarker((*it),0);
								}
								inside = true;
							}
						}
						ib++;
					}
					if (int(pabloBB.getLevel((*it))) > 0 && !inside){
						/**<Set to coarse outside the band if the octant has a level higher than 6.*/
						pabloBB.setMarker((*it),5-pabloBB.getLevel((*it)));
					}
				}

				itend = pabloBB.getPboundOctantsEnd();
				for (it=pabloBB.getPboundOctantsBegin(); it!=itend; ++it){
					bool inside = false;
					/**<Compute the nodes of the octant.*/
					vector<array<double,3> > nodes = pabloBB.getNodes((*it));
					/**<Compute the center of the octant.*/
					array<double,3> center = pabloBB.getCenter((*it));
					int ib = 0;
					while (!inside && ib<nb){
						double xc = BB[ib].c[0];
						double yc = BB[ib].c[1];
						double radius = BB[ib].r;
						/**<Set marker with condition on center or nodes of the octant.*/
						for (int j=0; j<4; j++){
							double x = nodes[j][0];
							double y = nodes[j][1];
							if ( ((!inside) &&
									(pow((x-xc),2.0)+pow((y-yc),2.0) <= 1.25*pow(radius,2.0) &&
											pow((x-xc),2.0)+pow((y-yc),2.0) >= 0.75*pow(radius,2.0)))
											|| ((!inside) && (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0) <= 1.25*pow(radius,2.0) &&
													pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0) >= 0.75*pow(radius,2.0)))){
								if (pabloBB.getLevel((*it)) < 9){
									/**<Set to refine inside a band around the interface of the bubbles.*/
									pabloBB.setMarker((*it),1);
								}
								else{
									pabloBB.setMarker((*it),0);
								}
								inside = true;
							}
						}
						ib++;
					}
					if (pabloBB.getLevel((*it)) > 0 && !inside){
						/**<Set to coarse outside the band if the octant has a level higher than 6.*/
						pabloBB.setMarker((*it),5-pabloBB.getLevel((*it)));
					}
				}

				/**<Adapt the octree.*/
				adapt = pabloBB.adapt();

			}

#if ENABLE_MPI==1
				/**<PARALLEL TEST: (Load)Balance the octree over the processes with communicating the data.*/
				pabloBB.loadBalance();
#endif

			/**<Update the connectivity and write the para_tree.*/
			pabloBB.updateConnectivity();
			pabloBB.write("PabloBubble_iter"+to_string(static_cast<unsigned long long>(iter)));

		}
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
