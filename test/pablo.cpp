#include "preprocessor_defines.dat"
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"

using namespace std;

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if NOMPI==0
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling Pablo Test routines*/

        test0() ;
        test1() ;
        test2() ;
        test12() ;
        test104() ;

#if NOMPI==0
	}

	MPI::Finalize();
#endif
}

// =================================================================================== //

void test0() {

    /**<Instantation of a 2D para_tree object with default constructor.*/
    Class_Para_Tree<2> ptreedefault;
    /**<Write the para_tree in physical domain.*/
    ptreedefault.write("Pablo0_default");
    /**<Write the para_tree in logical domain.*/
    ptreedefault.writeLogical("Pablo0_default_logical");

    /**<Set coordinates of the origin and size of a 2D custom para_tree object.*/
    double X, Y, Z, L;
    X = 10.0; Y = 20.0; Z = 0.0; L = 250.0;
    /**<Instantation of a 2D para_tree object with custom constructor.*/
    Class_Para_Tree<2> ptreecustom(X, Y, Z, L);
    /**<Write the para_tree in physical domain.*/
    ptreecustom.write("Pablo0_custom");
    /**<Write the para_tree in logical domain.*/
    ptreecustom.writeLogical("Pablo0_custom_logical");

    return 0;

}

// =================================================================================== //

void test1() {

    /**<Instantation of a 2D para_tree object.*/
    Class_Para_Tree<2> pablo1;

    /**<Compute the connectivity and write the para_tree.*/
    pablo1.computeConnectivity();
    pablo1.write("Pablo1_iter0");

    /**<Refine globally one level and write the para_tree.*/
    pablo1.adaptGlobalRefine();
    pablo1.updateConnectivity();
    pablo1.write("Pablo1_iter1");

    /**<Define a center point.*/
    double xc, yc;
    xc = yc = 0.5;

    /**<Set 2:1 balance only through faces.*/
    pablo1.setBalanceCodimension(1);

    /**<Set NO 2:1 balance in the right side of domain.*/
    uint32_t nocts = pablo1.getNumOctants();
    for (int i=0; i<nocts; i++){
        vector<double> center = pablo1.getCenter(i);
        double x = center[0];
        double y = center[1];
        if (x>xc)
            pablo1.setBalance(i,false);
    }

    /**<Define a radius.*/
    double radius = 0.4;

    /**<Simple adapt() nref1 times in the lower area of domain.*/
    int nref1 = 6;
    for (int iter=0; iter<nref1; iter++){
        nocts = pablo1.getNumOctants();
        for (int i=0; i<nocts; i++){
            /**<Extract Octant (pointer use).*/
            Class_Octant<2> *oct = pablo1.getOctant(i);
            /**<Compute center of the octant.*/
            vector<double> center = pablo1.getCenter(oct);
            double x = center[0];
            double y = center[1];

            /**<Set refinement marker=1 for octants inside a circle.*/
            if ((pow((x-xc),2.0)+pow((y-yc),2.0) < pow(radius,2.0)) &&
                    (y<yc)){
                pablo1.setMarker(oct, 1);
            }
        }
        /**<Adapt octree, update connectivity and write.*/
        pablo1.adapt();
        pablo1.updateConnectivity();
        pablo1.write("Pablo1_iter"+to_string(static_cast<unsigned long long>(iter+2)));
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
            nocts = pablo1.getNumOctants();
            for (int i=0; i<nocts; i++){
                /**<Compute center of the octant (index use).*/
                vector<double> center = pablo1.getCenter(i);
                double x = center[0];
                double y = center[1];
                if ((pow((x-xc),2.0)+pow((y-yc),2.0) < pow(radius,2.0)) &&
                        (y>yc) && iter<=nref2 && pablo1.getLevel(i)<=iter+1){

                    /**<Set refinement marker=1 for octants inside a circle.*/
                    pablo1.setMarker(i, 1);
                }
            }
            done = pablo1.adapt();
            pablo1.updateConnectivity();
            pablo1.write("Pablo1_iter"+to_string(static_cast<unsigned long long>(iter+nref1+2)));
        }
        iter++;
    }
    /**<Globally refine one level, update the connectivity and write the para_tree.*/
    pablo1.adaptGlobalRefine();
    pablo1.updateConnectivity();
    pablo1.write("Pablo1_iter"+to_string(static_cast<unsigned long long>(iter+nref1+3)));

    return 0;
}

// =================================================================================== //

void test2() {

    /**<Instantation of a 2D para_tree object.*/
    Class_Para_Tree<2> pablo2;

    /**<Set 2:1 balance only through faces.*/
    pablo2.setBalanceCodimension(1);
    uint32_t idx=0;
    pablo2.setBalance(idx,true);

    /**<Compute the connectivity and write the para_tree.*/
    pablo2.computeConnectivity();
    pablo2.write("Pablo2_iter0");

    /**<Refine globally two level and write the para_tree.*/
    for (int iter=1; iter<3; iter++){
        pablo2.adaptGlobalRefine();
        pablo2.updateConnectivity();
        pablo2.write("Pablo2_iter"+to_string(static_cast<unsigned long long>(iter)));
    }

    /**<Define a center point and a radius.*/
    double xc, yc;
    xc = yc = 0.5;
    double radius = 0.4;

    /**<Simple adapt() 6 times the octants with at least one node inside the circle.*/
    for (int iter=3; iter<9; iter++){
        uint32_t nocts = pablo2.getNumOctants();
        for (int i=0; i<nocts; i++){
            /**<Compute the nodes of the octant.*/
            vector<vector<double> > nodes = pablo2.getNodes(i);
            for (int j=0; j<global2D.nnodes; j++){
                double x = nodes[j][0];
                double y = nodes[j][1];
                if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
                    pablo2.setMarker(i, 1);
                }
            }
        }
        /**<Adapt octree.*/
        pablo2.adapt();

        /**<Update the connectivity and write the para_tree.*/
        pablo2.updateConnectivity();
        pablo2.write("Pablo2_iter"+to_string(static_cast<unsigned long long>(iter)));
    }

    return 0;
}

// =================================================================================== //

void test12() {

    /**<Instantation of a 2D para_tree object.*/
    Class_Para_Tree<2> pablo12;

    /**<Set NO 2:1 balance for the octree.*/
    uint32_t idx=0;
    pablo12.setBalance(idx,false);

    /**<Compute the connectivity and write the para_tree.*/
    pablo12.computeConnectivity();
    pablo12.write("Pablo12_iter0");

    /**<Refine globally two level and write the para_tree.*/
    for (int iter=1; iter<3; iter++){
        pablo12.adaptGlobalRefine();

        pablo12.updateConnectivity();
        pablo12.write("Pablo12_iter"+to_string(static_cast<unsigned long long>(iter)));
    }

#if NOMPI==0
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
            vector<vector<double> > nodes = pablo12.getNodes(i);
            for (int j=0; j<global2D.nnodes; j++){
                double x = nodes[j][0];
                double y = nodes[j][1];
                if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
                    pablo12.setMarker(i, 1);
                }
            }
        }

        /**<Adapt octree.*/
        pablo12.adapt();

#if NOMPI==0
        /**<(Load)Balance the octree over the processes.*/
        pablo12.loadBalance();
#endif

        /**<Update the connectivity and write the para_tree.*/
        pablo12.updateConnectivity();
        pablo12.write("Pablo12_iter"+to_string(static_cast<unsigned long long>(iter)));
    }

    return 0;
}


// =================================================================================== //

void test104() {

    int iter = 0;
    int dim = 3;

    /**<Instantation of a 3D para_tree object.*/
    Class_Para_Tree<3> pablo104;

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
        vector<vector<double> > nodes = pablo104.getNodes(i);
        for (int j=0; j<global3D.nnodes; j++){
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
    pablo104.writeTest("Pablo104_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

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
                    nfaces = global3D.nfaces;
                }
                else if (codim == 2){
                    nfaces = global3D.nedges;
                }
                else if (codim == 3){
                    nfaces = global3D.nnodes;
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
        pablo104.writeTest("Pablo104_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data_smooth);

        oct_data = oct_data_smooth;
    };

    return 0;
} ;

