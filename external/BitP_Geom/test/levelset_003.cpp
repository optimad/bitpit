/*!
 *	\date			10/jul/2014
 *	\authors		Alessandro Alaia
 *	\authors		Haysam Telib
 *	\authors		Edoardo Lombardi
 *	\version		0.1
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *
 *	\brief Level Set Class Demos
 */

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

//Standard Template Library
# include <ctime>
# include <chrono>

// bitpit
# include "BitP_Geom_BASE.hpp"
# include "BitP_Geom_LEVELSET.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace bitpit;


/*!Demo for 2D level set of complex geometries on a Pablo octree mesh.
 */
int main( int argc, char *argv[]){

#if NOMPI==0
        MPI::Init(argc,argv);
#endif
	// ========================================================================== //
	// VARIABLES DECLARATION                                                      //
	// ========================================================================== //

	// Input geometry
	int						nSTL = 1;
	vector<bitpit::SurfTriPatch>	STL(nSTL);

    int                     dimensions(2) ;
	PabloUniform            Mesh(dimensions);


	// Counters
	// none


	// ========================================================================== //
	// LOAD STL GEOMETRY                                                          //
	// ========================================================================== //
	{
		// Scope variables ------------------------------------------------------ //
		// none

		// Output message ------------------------------------------------------- //
		if (Mesh.getRank() == 0) cout << " - Loading stl geometry" << endl;

		// Load stl geometry ---------------------------------------------------- //
		STL[0].Import_dgf("./data/naca0012.dgf");
// 		STL[0].Export_vtu("naca0012.vtu");
//		STL[0].Import_dgf("./data/2Dcircle.dgf");

		for(int i=0; i<nSTL; i++){


            STL[i].SetDimensions(dimensions) ;
            STL[i].Clean() ;
            STL[i].GenerateVNormals(3-dimensions) ;
            STL[i].DestroyAdjacency() ;
            STL[i].DestroyEdge() ;

			// Display info --------------------------------------------------------- //
			if (Mesh.getRank() == 0) cout << "n. vertex: " << STL[i].nVertex << endl;
			if (Mesh.getRank() == 0) cout << "n. simplex: " << STL[i].nSimplex << endl;

		}
	}

    // ========================================================================== //
    // CREATE MESH                                                                //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        int             i,j;
        double          dx, dy, dz;
        array<double,3> meshMin, meshMax ;

        meshMin.fill(0.0) ;
        meshMax.fill(0.0) ;

        vector< array<double,3> > stlMin(nSTL), stlMax(nSTL) ;
        
        
        // Output message ------------------------------------------------------- //
        cout << " - Setting mesh" << endl;

        // Bounding box for input geometry -------------------------------------- //
        for (j=0; j<nSTL; j++){
            STL[j].BoundingBox( stlMin[j], stlMax[j] ) ;
        }//next j

        CGElem::UnionAABB( stlMin, stlMax, meshMin, meshMax ) ;

        meshMin = meshMin - 0.1*(meshMax -meshMin) ;
        meshMax = meshMax + 0.1*(meshMax -meshMin) ;

        // Set mesh ------------------------------------------------------------- //
        Mesh.setL( max(meshMax[0]-meshMin[0],meshMax[1]-meshMin[1]) ) ;
        Mesh.setOrigin( meshMin ) ;

		if (Mesh.getRank() == 0) cout << " - Setting pablo mesh" << endl;
		for (i=0; i<5; i++){
			Mesh.adaptGlobalRefine();
		}

#if NOMPI==0
		//if (Mesh.getRank() == 0) cout << " - Balancing pablo mesh" << endl;
		//Mesh.loadBalance();
#endif

		for (i=0; i<2; i++){
			Mesh.adaptGlobalRefine();
		}


    }

	

	// ========================================================================== //
	// COMPUTE LEVEL SET                                                          //
	// ========================================================================== //

	{
		vector< LevelSetTriOctr > LSP;
		LSP.resize(nSTL);

		for (int istl=0; istl<nSTL; istl++){
			LSP[istl].setLog("ls"+to_string(istl)+".log");
			LSP[istl].setMesh(&Mesh);
			LSP[istl].setTri(&STL[istl]);
		}

		int i = 0;

		// Compute level set  in narrow band ------------------------------------- //
		for (int istl=0; istl<nSTL; istl++){
			LSP[istl].compute();
		}

		// Export level set ------------------------------------------------------- //
		if (Mesh.getRank() == 0) cout << " - Exporting data" << endl;
		for (int istl=0; istl<nSTL; istl++){
			Mesh.writeTest("ls"+to_string(istl)+"_it0", LSP[istl].getLS());//, LSGH);
		}
		if (Mesh.getRank() == 0) cout << " - Exported data" << endl;


		// Adapt mesh by level set on STL ----------------------------------------- //
//ht		i = 0;
//ht		vector<int> targetstep(nSTL,4);
//ht
//ht		while(targetstep != levSTL){
//ht
//ht			for (int istl=0; istl<nSTL; istl++){
//ht
//ht				targetstep[istl] = min(targetstep[istl]+1,levSTL[istl]);
//ht
//ht				int target = targetstep[istl];
//ht				bool done = true;
//ht				while (done){
//ht					start = std::chrono::system_clock::now();
//ht					done = adaptOnStl( LSP[istl], target);
//ht					end = chrono::system_clock::now();
//ht					elapsed_seconds = chrono::duration_cast<chrono::milliseconds>(end-start).count();
//ht					if (Mesh.rank == 0) cout << "elapsed time: " << elapsed_seconds << " ms" << endl;
//ht
//ht					if (done){
//ht						i++;
//ht						// Export level set ----------------------------------------------------- //
//ht						if (Mesh.rank == 0) cout << " - Exporting data" << endl;
//ht						Mesh.computeConnectivity();
//ht						for (int jstl=0; jstl<nSTL; jstl++)
//ht							Mesh.writeTest("ls"+to_string(jstl)+"_it"+to_string(i), LSP[jstl].getSdf());//, LSGH);
//ht						Mesh.clearConnectivity();
//ht						if (Mesh.rank == 0) cout << " - Exported data" << endl;
//ht					}
//ht
//ht				}//end while done
//ht			}//next istl
//ht
//ht#if NOMPI==0
//ht			//Load Balance
//ht			if (Mesh.rank == 0) cout << " - Load Balance" << endl;
//ht			//User_VLSData_LB<vector<Class_LevelSet_Stl<Class_Para_Tree<2> > > > data_vls(LSP);
//ht			User_LSData_LB<Class_LevelSet_Stl<Class_Para_Tree<2> > > data_ls(LSP[0]);
//ht#endif
//ht
//ht			int nocts = Mesh.getNumOctants();
//ht			vector<double> weight(nocts, 1.0);
//ht			for (int istl=0; istl<nSTL; istl++){
//ht				for (int j=0; j<nocts; j++){
//ht					weight[j] += LSP[istl].getSimplexList(j).size();
//ht				}
//ht			}
//ht
//ht#if NOMPI==0
//ht			Mesh.loadBalance(data_ls, &weight);
//ht			postLoadBalance(LSP);
//ht#endif
//ht			i++;
//ht			// Export level set ----------------------------------------------------- //
//ht			if (Mesh.rank == 0) cout << " - Exporting data" << endl;
//ht			Mesh.computeConnectivity();
//ht			for (int jstl=0; jstl<nSTL; jstl++)
//ht				Mesh.writeTest("ls"+to_string(jstl)+"_it"+to_string(i), LSP[jstl].getSdf());//, LSGH);
//ht			Mesh.clearConnectivity();
//ht			if (Mesh.rank == 0) cout << " - Exported data" << endl;
//ht
//ht		}//end while targetstep
//ht
	}
	
#if NOMPI==0
    MPI::Finalize();
#endif

    return 0;

};


