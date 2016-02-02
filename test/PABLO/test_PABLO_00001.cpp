#include "ParaTree.hpp"
#include "PabloUniform.hpp"

using namespace std;

// =================================================================================== //

void test001() {

    /**<Instantation of a 2D para_tree object with default constructor.*/
    ParaTree ptreedefault;
    /**<Write the para_tree in physical domain.*/
    ptreedefault.write("Pablo001_default");

    /**<Set coordinates of the origin and size of a 2D custom para_tree object.*/
    double X, Y, Z, L;
    X = 10.0; Y = 20.0; Z = 0.0; L = 250.0;
    /**<Instantation of a 2D para_tree object with custom constructor.*/
    PabloUniform ptreecustom(X,Y,Z,L);
    /**<Write the para_tree in physical domain.*/
    ptreecustom.write("Pablo001_custom");

    return ;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling Pablo Test routines*/

        test001() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
