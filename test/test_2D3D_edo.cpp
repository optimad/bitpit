/*
 * test_2D3D.cpp
 *
 *  Created on: 14/apr/2014
 *      Author: Marco Cisternino
 */

#include "preprocessor_defines.dat"
#include <mpi.h>
#include "global.hpp"
#include "Class_Para_Tree.hpp"
#include "ioFunct.hpp"

using namespace std;

int main(int argc, char *argv[]) {

/*
	MPI::Init(argc, argv);

		{
			double X, Y, Z, L;
			uint8_t level0 = MAX_LEVEL_2D;
			X = 0.0; Y = 0.0; Z = 0.0; L = 100.0;
			Class_Para_Tree<2> ptree(X, Y, Z, L);

			clock_t start = clock();
			clock_t end = clock();

			ptree.octree.setBalance(0,false);
			uint8_t levell = 4;
			ptree.octree.setMarker(0, levell);
			bool done = ptree.adapt();
			ptree.loadBalance();
			uint64_t nocts = ptree.octree.getNumOctants();


			// LOGO TEST
			vector<double> C1 = {double(global2D.max_length)*0.5, double(global2D.max_length)*0.5, 0.0};
			double R1 = double(global2D.max_length)*0.4;
			vector<double> C2 = {double(global2D.max_length)*0.35, double(global2D.max_length)*0.65, 0.0};
			double R2 = double(global2D.max_length)*0.1;
			vector<double> C3 = {double(global2D.max_length)*0.65, double(global2D.max_length)*0.65, 0.0};
			double R3 = double(global2D.max_length)*0.1;
			vector<double> C4 = {double(global2D.max_length)*0.325, double(global2D.max_length)*0.675, 0.0};
			double R4 = double(global2D.max_length)*0.025;
			vector<double> C5 = {double(global2D.max_length)*0.675, double(global2D.max_length)*0.625, 0.0};
			double R5 = double(global2D.max_length)*0.025;

 			int nref = 5;
			for (int l=0; l<nref; l++){
				for (int i=0; i<nocts; i++){
					double* center;
					Class_Octant<2> oct = ptree.octree.extractOctant(i);
					center = oct.getCenter();
					//				if (sqrt(pow((center[0]-double(global2D.max_length)*0.5),2.0)+pow((center[1]-double(global2D.max_length)*0.5),2.0)+pow((center[2]-double(global2D.max_length)*0.5),2.0)) <= double(global2D.max_length)*0.4){
					if((pow((center[0]-C1[0]),2.0)+pow((center[1]-C1[1]),2.0)+pow((center[2]-C1[2]),2.0) <= pow(R1+oct.getSize(),2.0))
							&& (pow((center[0]-C1[0]),2.0)+pow((center[1]-C1[1]),2.0)+pow((center[2]-C1[2]),2.0) >= pow(R1-oct.getSize(),2.0)) ){
						ptree.octree.setMarker(i,1);
					}
					if((pow((center[0]-C2[0]),2.0)+pow((center[1]-C2[1]),2.0)+pow((center[2]-C2[2]),2.0) <= pow(R2+oct.getSize(),2.0))
							&& (pow((center[0]-C2[0]),2.0)+pow((center[1]-C2[1]),2.0)+pow((center[2]-C2[2]),2.0) >= pow(R2-oct.getSize(),2.0)) ){
						ptree.octree.setMarker(i,1);
					}
					if((pow((center[0]-C3[0]),2.0)+pow((center[1]-C3[1]),2.0)+pow((center[2]-C3[2]),2.0) <= pow(R3+oct.getSize(),2.0))
							&& (pow((center[0]-C3[0]),2.0)+pow((center[1]-C3[1]),2.0)+pow((center[2]-C3[2]),2.0) >= pow(R3-oct.getSize(),2.0)) ){
						ptree.octree.setMarker(i,1);
					}
					if((global2D.max_length*0.1*sin(4.0*M_PI*center[0]/global2D.max_length + M_PI*0.5) + global2D.max_length*0.4 <= center[1]+7.5*oct.getSize())
						 && (global2D.max_length*0.1*sin(4.0*M_PI*center[0]/global2D.max_length + M_PI*0.5) + global2D.max_length*0.4 >= center[1]-7.5*oct.getSize())){
	//					 && center[0] >= global2D.max_length*0.05 && center[0] <= global2D.max_length*0.95){
						ptree.octree.setMarker(i,1);
					}
					if((global2D.max_length*0.025*sin(2.0*4.0*M_PI*center[0]/global2D.max_length + M_PI*1.5) + global2D.max_length*0.75 <= center[1]+1.0*oct.getSize())
						 && (global2D.max_length*0.025*sin(2.0*4.0*M_PI*center[0]/global2D.max_length + M_PI*1.5) + global2D.max_length*0.75 >= center[1]-1.0*oct.getSize())
						 && (center[0] >= global2D.max_length*0.35) && (center[0] <= global2D.max_length*0.65)){
						//ptree.octree.setMarker(i,1);
					}
					//MexHat
					if((center[1]-global2D.max_length*0.75 <= sin(pow(center[0],2.0)+pow(center[2],2.0))/(pow(center[0],2.0)+pow(center[2],2.0))+oct.getSize())
							&& (center[1]-global2D.max_length*0.75 >= sin(pow(center[0],2.0)+pow(center[2],2.0))/(pow(center[0],2.0)+pow(center[2],2.0))-oct.getSize())){
	//					ptree.octree.setMarker(i,1);
					}
					delete[] center;
				}
				typename Class_Para_Tree<2>::u32vector mapper,internal,ghosts,borders;
				ptree.octree.computeIntersections();
				bool done = ptree.adapt(mapper,internal,ghosts,borders);
				ptree.loadBalance();
				nocts = ptree.octree.getNumOctants();


				ptree.octree.updateConnectivity();
//				ptree.octree.updateGhostsConnectivity();
				writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,("Pablo_"+to_string(l+3)));

				end = clock();
				float seconds = (float)(end - start) / CLOCKS_PER_SEC;
				writeLog(" ");
				writeLog("---------------------------------------------");
				writeLog(" ");
				writeLog(" CPU time (sec)			:	"+to_string(seconds));
				writeLog(" ");
				writeLog("---------------------------------------------");
			}


//			int nref = 4;
//			for (int l=nref+1; l<nref+8; l++){
//				for (int i=0; i<nocts; i++){
//					double* center;
//					Class_Octant<2> oct = ptree.octree.extractOctant(i);
//					center = oct.getCenter();
//					//				if (sqrt(pow((center[0]-double(global2D.max_length)*0.5),2.0)+pow((center[1]-double(global2D.max_length)*0.5),2.0)+pow((center[2]-double(global2D.max_length)*0.5),2.0)) <= double(global2D.max_length)*0.4){
//					if((pow((center[0]-C2[0]),2.0)+pow((center[1]-C2[1]),2.0)+pow((center[2]-C2[2]),2.0) <= pow(R2,2.0))){
//						ptree.octree.setMarker(i,-1);
//						ptree.octree.setBalance(i,false);
//					}
//					if((pow((center[0]-C3[0]),2.0)+pow((center[1]-C3[1]),2.0) <= pow(R3,2.0))){
//						ptree.octree.setMarker(i,-1);
//						ptree.octree.setBalance(i,false);
//					}
//					if((global2D.max_length*0.1*sin(4.0*M_PI*center[0]/global2D.max_length + M_PI*0.5) + global2D.max_length*0.4 <= center[1]+6.5*oct.getSize())
//						 && (global2D.max_length*0.1*sin(4.0*M_PI*center[0]/global2D.max_length + M_PI*0.5) + global2D.max_length*0.4 >= center[1]-6.5*oct.getSize())){
//						ptree.octree.setMarker(i,-1);
//						ptree.octree.setBalance(i,false);
//					}
//					delete[] center;
//				}
//				//TEST ADAPT WITH MAPPER
//				typename Class_Para_Tree<2>::u32vector mapper,internal,ghosts,borders;
//				ptree.octree.computeIntersections();
//				bool done = ptree.adapt(mapper,internal,ghosts,borders);
//				ptree.loadBalance();
//				nocts = ptree.octree.getNumOctants();
//
//				ptree.octree.updateConnectivity();
////				ptree.octree.updateGhostsConnectivity();
//				writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,("Pablo_"+to_string(l+3)));
//
//				end = clock();
//				float seconds = (float)(end - start) / CLOCKS_PER_SEC;
//				writeLog(" ");
//				writeLog("---------------------------------------------");
//				writeLog(" ");
//				writeLog(" CPU time (sec)			:	"+to_string(seconds));
//				writeLog(" ");
//				writeLog("---------------------------------------------");
//			}

//			ptree.octree.updateConnectivity();
////				ptree.octree.updateGhostsConnectivity();
//			writeLocalTree(ptree.octree.nodes,ptree.octree.connectivity,ptree.octree.ghostsnodes,ptree.octree.ghostsconnectivity,ptree,("Pablo_"+to_string(0)));

		}

	MPI::Finalize();
*/

	return 0;

}
