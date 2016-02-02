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
void test004() {

    int iter = 0;
    int dim = 3;

    /**<Instantation of a 3D para_tree object.*/
    ParaTree pablo104(3);

    /**<Refine globally four level and write the para_tree.*/
    for (iter=1; iter<4; iter++){
        pablo104.adaptGlobalRefine();
    }

    /**<Define a center point and a radius.*/
    double xc, yc;
    xc = yc = 0.5;
    double radius = 0.25;

    /**<Define vectors of data.*/
    uint32_t nocts = pablo104.getNumOctants();
    vector<double> oct_data(nocts, 0.0);

    /**<Assign a data to the octants with at least one node inside the cylinder.*/
    for (int i=0; i<nocts; i++){
        /**<Compute the nodes of the octant.*/
        vector<array<double,3> > nodes = pablo104.getNodes(i);
        for (int j=0; j<8; j++){
            double x = nodes[j][0];
            double y = nodes[j][1];
            if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
                oct_data[i] = 1.0;
            }
        }
    }

    /**<Update the connectivity and write the para_tree.*/
    iter = 0;
    pablo104.updateConnectivity();
    pablo104.writeTest("Pablo004_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

    /**<Smoothing iterations on initial data*/
    int start = 1;
    for (iter=start; iter<start+25; iter++){
        vector<double> oct_data_smooth(nocts, 0.0);
        vector<uint32_t> neigh, neigh_t;
        vector<bool> isghost, isghost_t;
        uint8_t iface, codim, nfaces;
        for (int i=0; i<nocts; i++){
            neigh.clear();
            isghost.clear();
            /**<Find neighbours through faces (codim=1), edges (codim=2) and nodes (codim=3) of the octants*/
            for (codim=1; codim<dim+1; codim++){
                if (codim == 1){
                    nfaces = 6;
                }
                else if (codim == 2){
                    nfaces = 12;
                }
                else if (codim == 3){
                    nfaces = 8;
                }
                for (iface=0; iface<nfaces; iface++){
                    pablo104.findNeighbours(i,iface,codim,neigh_t,isghost_t);
                    neigh.insert(neigh.end(), neigh_t.begin(), neigh_t.end());
                    isghost.insert(isghost.end(), isghost_t.begin(), isghost_t.end());
                }
            }

            /**<Smoothing data with the average over the one ring neighbours of octants*/
            oct_data_smooth[i] = oct_data[i]/(neigh.size()+1);
            for (int j=0; j<neigh.size(); j++){
                if (isghost[j]){
                    /**< Do nothing - No ghosts: is a serial test.*/
                }
                else{
                    oct_data_smooth[i] += oct_data[neigh[j]]/(neigh.size()+1);
                }
            }
        }

        /**<Update the connectivity and write the para_tree.*/
        pablo104.updateConnectivity();
        pablo104.writeTest("Pablo004_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data_smooth);

        oct_data = oct_data_smooth;
    };

    return ;
} ;

// =================================================================================== //
int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling Pablo Test routines*/
        test004() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
