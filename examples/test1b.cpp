#include "ParaTree.hpp"

using namespace std;

// =================================================================================== //

int main(int argc, char *argv[]) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Instantation of a 2D para_tree object.*/
		ParaTree pablo1;

		/**<Compute the connectivity and write the para_tree.*/
		pablo1.computeConnectivity();
		pablo1.write("Pablo1b_iter0");

		/**<Refine globally one level and write the para_tree.*/
		pablo1.adaptGlobalRefine();
		pablo1.updateConnectivity();
		pablo1.write("Pablo1b_iter1");

		/**<Define a center point.*/
		double xc, yc;
		xc = yc = 0.5;

		/**<Set 2:1 balance through faces and edges.*/
		pablo1.setBalanceCodimension(2);

		/**<Set NO 2:1 balance in the right side of domain.*/
		uint32_t nocts = pablo1.getNumOctants();
		for (int i=0; i<nocts; i++){
			array<double,3> center = pablo1.getCenter(i);
			double x = center[0];
			double y = center[1];
			if (x>xc)
				pablo1.setBalance(i,false);
		}

		/**<Define a radius.*/
		double radius = 0.4;

		/**<Simple adapt() nref1 times in lower area of domain.*/
		int nref1 = 6;
		for (int iter=0; iter<nref1; iter++){
			nocts = pablo1.getNumOctants();
			for (int i=0; i<nocts; i++){
				/**<Extract Octant (pointer use).*/
				Octant *oct = pablo1.getOctant(i);
				/**<Compute center of the octant.*/
				array<double,3> center = pablo1.getCenter(oct);
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
			pablo1.write("Pablo1b_iter"+to_string(static_cast<unsigned long long>(iter+2)));
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
					array<double,3> center = pablo1.getCenter(i);
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
				pablo1.write("Pablo1b_iter"+to_string(static_cast<unsigned long long>(iter+nref1+2)));
			}
			iter++;
		}
		/**<Globally refine one level, update the connectivity and write the para_tree.*/
		pablo1.adaptGlobalRefine();
		pablo1.updateConnectivity();
		pablo1.write("Pablo1b_iter"+to_string(static_cast<unsigned long long>(iter+nref1+3)));
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
