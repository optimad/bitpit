/*
 ============================================================================
 Name        : PABLO.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Compute Pi in MPI C++
 ============================================================================
 */
#include "preprocessor_defines.dat"
#include <math.h> 
#include "mpi.h" 
#include <iostream>
#include "Class_Para_Tree.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	MPI::Init(argc, argv);

	Class_Para_Tree ptree;


	ptree.octree.setMarker(0,1);
	ptree.octree.refine();
	ptree.octree.setMarker(2,1);
	ptree.octree.refine();

	Class_Octant oct_test;
	uint8_t sizehf=0;
	uint8_t sizem=0;
	uint8_t sizen=0;
	bool isghost;
	for (int i=0; i<6; i++){//nface; i++){
		oct_test =	ptree.octree.extractOctant(6);
		cout << "Bound " << i << "-face : " << oct_test.getBound(i) << endl;
/*		uint64_t *hfneigh = oct_test.computeHalfSizeMorton(i,sizehf);
		for (int j=0; j<sizehf; j++){
			cout << "Morton half-size idx=1 iface " << i << " : " << hfneigh[j] << endl;
		}
		uint64_t *mneigh = oct_test.computeMinSizeMorton(i,3, sizem);
		for (int j=0; j<sizem; j++){
			cout << "Morton min-size idx=1 iface " << i << " : " << mneigh[j] << endl;
		}*/
		uint64_t *idxneigh = ptree.octree.findNeighbours(6, i, sizen, isghost);
		for (int j=0; j<sizen; j++){
			cout << "Idx neigh idx=6 iface " << i << " : " << idxneigh[j] << endl;
		}
	}
/*
	uint32_t (*nodes)[DIM] = oct_test.getNodes();
	for(int i=0; i<nnodes; i++){
		for (int j=0; j<DIM; j++){
			cout << "node " << i << "  coord " << j << " : " << nodes[i][j] << endl;
		}
	}
	uint64_t numoctants = ptree.octree.getNumOctants();
	cout << "Num Octants : " << numoctants << endl;
	ptree.octree.computeConnectivity();
	for (int i=0; i<ptree.octree.nodes.size(); i++){
		cout << " x " << ptree.octree.nodes[i][0] << "  y " << ptree.octree.nodes[i][1] << "  z " << ptree.octree.nodes[i][2] << " morton " << mortonEncode_magicbits(ptree.octree.nodes[i][0], ptree.octree.nodes[i][1], ptree.octree.nodes[i][2]) << endl;
	}
	for (int i=0; i<numoctants; i++){
		cout << ptree.octree.connectivity[i][0] << endl;
	}
*/

	ptree.octree.clearConnectivity();
	writeLog("PABLO");

	MPI::Finalize();

	return 0;
}

