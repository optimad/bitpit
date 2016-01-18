#include "ParaTree.hpp"

using namespace std;

// =================================================================================== //
void testParallel001() {

    /**<Instantation of a 2D para_tree object.*/
    ParaTree pablo12;

    /**<Set NO 2:1 balance for the octree.*/
    uint32_t idx=0;
    pablo12.setBalance(idx,false);

    /**<Compute the connectivity and write the para_tree.*/
    pablo12.computeConnectivity();
    pablo12.write("PabloParallel001_iter0");

    /**<Refine globally two level and write the para_tree.*/
    for (int iter=1; iter<3; iter++){
        pablo12.adaptGlobalRefine();
        pablo12.updateConnectivity();
        pablo12.write("PabloParallel001_iter"+to_string(static_cast<unsigned long long>(iter)));
    }

#if ENABLE_MPI==1
    /**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
    pablo12.loadBalance();
#endif

    /**<Define a center point and a radius.*/
    double xc, yc;
    xc = yc = 0.5;
    double radius = 0.4;

    /**<Simple adapt() (refine) 6 times the octants with at least one node inside the circle.*/
    for (int iter=3; iter<9; iter++){
        uint32_t nocts = pablo12.getNumOctants();
        for (int i=0; i<nocts; i++){
            /**<Compute the nodes of the octant.*/
            vector<array<double,3> > nodes = pablo12.getNodes(i);
            for (int j=0; j<4; j++){
                double x = nodes[j][0];
                double y = nodes[j][1];
                if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
                    pablo12.setMarker(i, 1);
                }
            }
        }

        /**<Adapt octree.*/
        pablo12.adapt();

#if ENABLE_MPI==1
        /**<(Load)Balance the octree over the processes.*/
        pablo12.loadBalance();
#endif

        /**<Update the connectivity and write the para_tree.*/
        pablo12.updateConnectivity();
        pablo12.write("PabloParallel001_iter"+to_string(static_cast<unsigned long long>(iter)));
    }

    return ;
}

// =================================================================================== //
int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling Pablo Test routines*/

        testParallel001() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
