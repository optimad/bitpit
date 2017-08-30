/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
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

#if BITPIT_ENABLE_MPI==1
#   include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_PABLO.hpp"
#include "bitpit_IO.hpp"

using namespace std;
using namespace bitpit;

/*!
* Subtest 001
*
* Testing refinement.
*/
int subtest_001()
{
    /**<Instantation of a 2D para_tree object.*/
    ParaTree pablo(2);

    /**<Compute the connectivity and write the para_tree.*/
    pablo.computeConnectivity();
    pablo.write("Pablo002_iter0");

    /**<Refine globally one level and write the para_tree.*/
    pablo.adaptGlobalRefine();
    pablo.updateConnectivity();
    pablo.write("Pablo002_iter1");

    /**<Define a center point.*/
    double xc, yc;
    xc = yc = 0.5;

    /**<Set 2:1 balance only through faces.*/
    pablo.setBalanceCodimension(1);

    /**<Set NO 2:1 balance in the right side of domain.*/
    uint32_t nocts = pablo.getNumOctants();
    for (unsigned int i=0; i<nocts; i++){
        array<double,3> center = pablo.getCenter(i);
        double x = center[0];
        if (x>xc)
            pablo.setBalance(i,false);
    }

    /**<Define a radius.*/
    double radius = 0.4;

    /**<Simple adapt() nref1 times in the lower area of domain.*/
    int nref1 = 6;
    for (int iter=0; iter<nref1; iter++){
        nocts = pablo.getNumOctants();
        for (unsigned int i=0; i<nocts; i++){
            /**<Extract Octant (pointer use).*/
            Octant *oct = pablo.getOctant(i);
            /**<Compute center of the octant.*/
            array<double,3> center = pablo.getCenter(oct);
            double x = center[0];
            double y = center[1];
            /**<Set refinement marker=1 for octants inside a circle.*/
            if ((pow((x-xc),2.0)+pow((y-yc),2.0) < pow(radius,2.0)) &&
                    (y<yc)){
                pablo.setMarker(oct, 1);
            }
        }
        /**<Adapt octree, update connectivity and write.*/
        pablo.adapt();
        pablo.updateConnectivity();
        pablo.write("Pablo002_iter"+to_string(static_cast<unsigned long long>(iter+2)));
    }


    /**<While adapt() nref2 times in the upper area of domain.
     * (Useful if you work with center of octants) */
    int nref2 = 5;
    int iter = 0;
    bool done = true;
    while(iter<=nref2){
        done = true;
        while(done)
        {
            nocts = pablo.getNumOctants();
            for (unsigned int i=0; i<nocts; i++){
                /**<Compute center of the octant (index use).*/
                array<double,3> center = pablo.getCenter(i);
                double x = center[0];
                double y = center[1];
                if ((pow((x-xc),2.0)+pow((y-yc),2.0) < pow(radius,2.0)) &&
                        (y>yc) && iter<=nref2 && pablo.getLevel(i)<=iter+1){

                    /**<Set refinement marker=1 for octants inside a circle.*/
                    pablo.setMarker(i, 1);
                }
            }
            done = pablo.adapt();
            pablo.updateConnectivity();
            pablo.write("Pablo002_iter"+to_string(static_cast<unsigned long long>(iter+nref1+2)));
        }
        iter++;
    }
    /**<Globally refine one level, update the connectivity and write the para_tree.*/
    pablo.adaptGlobalRefine();
    pablo.updateConnectivity();
    pablo.write("Pablo002_iter"+to_string(static_cast<unsigned long long>(iter+nref1+3)));

    return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    int nProcs;
    int rank;
#if BITPIT_ENABLE_MPI==1
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    nProcs = 1;
    rank   = 0;
#endif

    // Initialize the logger
    log::manager().initialize(log::SEPARATE, false, nProcs, rank);
    log::cout() << fileVerbosity(log::NORMAL);
    log::cout() << consoleVerbosity(log::QUIET);

    // Run the subtests
    log::cout() << "Testing refinement" << std::endl;

    int status;
    try {
        status = subtest_001();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
