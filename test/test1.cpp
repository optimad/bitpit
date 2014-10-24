#include "preprocessor_defines.dat"
#include <mpi.h>
#include "Class_Global.hpp"
#include "Class_Para_Tree.hpp"
#include "ioFunct.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	{

		Class_Para_Tree<2> pablo1;
		pablo1.computeConnectivity();
		pablo1.write("Pablo1_iter0");
		pablo1.adaptGlobalRefine();
		pablo1.updateConnectivity();
		pablo1.write("Pablo1_iter1");


		double xc, yc;
		xc = yc = 0.5;

		//Set NO 2:1 balance right area of domain
		uint32_t nocts = pablo1.getNumOctants();
		for (int i=0; i<nocts; i++){
			vector<double> center = pablo1.getCenter(i);
			double x = center[0];
			double y = center[1];
			if (x>xc)
				pablo1.setBalance(i,false);
		}

		double radius = 0.4;

		// Simple adapt() in upper area of domain
		int nref1 = 6;
		for (int iter=0; iter<nref1; iter++){
			nocts = pablo1.getNumOctants();
			for (int i=0; i<nocts; i++){
				Class_Octant<2> *oct = pablo1.getOctant(i);
				vector<double> center = pablo1.getCenter(oct);
				double x = center[0];
				double y = center[1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) < pow(radius,2.0)) &&
						(y<yc)){
					pablo1.setMarker(oct, 1);
				}
			}
			pablo1.adapt();
			pablo1.updateConnectivity();
			pablo1.write("Pablo1_iter"+to_string(iter+2));
		}

		// While adapt() in downer area of domain
		int nref2 = 5;
		int iter = 0;
		bool done = true;
		while(iter<=nref2){
			done = true;
		while(done)
		{
			nocts = pablo1.getNumOctants();
			for (int i=0; i<nocts; i++){
				Class_Octant<2> *oct = pablo1.getOctant(i);
				vector<double> center = pablo1.getCenter(oct);
				double x = center[0];
				double y = center[1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) < pow(radius,2.0)) &&
						(y>yc) && iter<=nref2 && oct->getLevel()<=iter+1){
					pablo1.setMarker(oct, 1);
				}
			}
			done = pablo1.adapt();
			pablo1.updateConnectivity();
			pablo1.write("Pablo1_iter"+to_string(iter+nref1+2));
		}
		iter++;
		}
		pablo1.adaptGlobalRefine();
		pablo1.updateConnectivity();
		pablo1.write("Pablo1_iter"+to_string(iter+nref1+3));
	}

	MPI::Finalize();

}

