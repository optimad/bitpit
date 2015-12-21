#include "classParaTree.hpp"

using namespace std;

// =================================================================================== //

void test0() {

    /**<Instantation of a 2D para_tree object with default constructor.*/
    classParaTree ptreedefault;
    /**<Write the para_tree in physical domain.*/
    ptreedefault.write("Pablo0_default");
//    /**<Write the para_tree in logical domain.*/
//    ptreedefault.writeLogical("Pablo0_default_logical");

    /**<Set coordinates of the origin and size of a 2D custom para_tree object.*/
    double X, Y, Z, L;
    X = 10.0; Y = 20.0; Z = 0.0; L = 250.0;
    /**<Instantation of a 2D para_tree object with custom constructor.*/
    classParaTree ptreecustom(X, Y, Z, L);
    /**<Write the para_tree in physical domain.*/
    ptreecustom.write("Pablo0_custom");
//    /**<Write the para_tree in logical domain.*/
//    ptreecustom.writeLogical("Pablo0_custom_logical");

    return ;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if NOMPI==0
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling Pablo Test routines*/

        test0() ;

#if NOMPI==0
	}

	MPI::Finalize();
#endif
}
