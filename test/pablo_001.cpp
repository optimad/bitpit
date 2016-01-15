#include "ClassParaTree.hpp"

using namespace std;

// =================================================================================== //

void test0() {

    /**<Instantation of a 2D para_tree object with default constructor.*/
    ClassParaTree ptreedefault;
    /**<Write the para_tree in physical domain.*/
    ptreedefault.write("Pablo0_default");

    /**<Set coordinates of the origin and size of a 2D custom para_tree object.*/
    double X, Y, Z, L;
    X = 10.0; Y = 20.0; Z = 0.0; L = 250.0;
    /**<Instantation of a 2D para_tree object with custom constructor.*/
    ClassParaTree ptreecustom(X, Y, Z, L);
    /**<Write the para_tree in physical domain.*/
    ptreecustom.write("Pablo0_custom");

    return ;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling Pablo Test routines*/

        test0() ;

#if ENABLE_MPI
	}

	MPI::Finalize();
#endif
}
