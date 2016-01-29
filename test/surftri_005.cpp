// ========================================================================== //
//                - Class_SurfTri - Examples of usage                         //
//                                                                            //
// Example of usage for Class_SurfTri (grid manager for unstructured meshes)  //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author   : Alessandro Alaia                                                //
// Version  : v3.0                                                            //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDE                                                                    //
// ========================================================================== //

// Standard Template Library
# include <iostream>
# include <chrono>
# include <ctime>

# include <Class_SurfTri.hpp>

// -------------------------------------------------------------------------- //
void Demo_CleaningTools3(
        void
        ) {

    // ========================================================================== //
    // void Demo_CleaningTools3(                                                  //
    //     void)                                                                  //
    //                                                                            //
    // Demo for Class_SurfTri cleaning tools. check on TrueDoubleSimplex and      //
    // Ring with safe failure are performed
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    Class_SurfTri           Tri;
    ivector2D		        E, S2E;
    // Counters
    // none

    // ========================================================================== //
    // OUTPUT MESSAGE                                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << "============= Class_SurfTri: cleaning tools 3 demo ==============" << endl;
    }

    // ========================================================================== //
    // LOAD TASSELATION FROM STL FILE                                             //
    // ========================================================================== //
    {
        // Import from dgf file ------------------------------------------------- //
        Tri.Import_stl("./data/buddha.stl", true);

        cout<<"number of vertices  "<<Tri.nVertex<<endl;  
        cout<<"number of simplicies  "<<Tri.nSimplex<<endl;  
        // Output stats --------------------------------------------------------- //

        Tri.SetTolerance();				

    }

    // ========================================================================== //
    // REMOVE DOUBLE VERTICES                                                     //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Finding duplicated vertices" << endl;

        // Output stats --------------------------------------------------------- //
        ivector1D list = Tri.FindDoubleVertex();  
        Tri.RemoveDoubleVertex();
        cout << endl << " - Going to remove  " << list.size() <<"   vertices"<<endl<<endl;				
        Tri.BuildAdjacency();	
    }

    // ========================================================================== //
    // REMOVE 0 Area SIMPLICIES                                                 //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Finding 0 area simplicies" << endl;

        ivector1D list;
        // Output stats --------------------------------------------------------- //
        list = Tri.Find0AreaSimplex();
        cout << endl << " - Going to remove  " << list.size() << " 0 area simplicies"<<endl<<endl;				
        Tri.RemoveSimplex(list);
        Tri.UpdateAdjacency(list);
    }

    // ========================================================================== //
    // REMOVE ISOLATED SIMPLICIES                                                 //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Finding isolated simplicies" << endl;

        ivector1D list;
        // Output stats --------------------------------------------------------- //
        list = Tri.FindIsolatedSimplex();
        cout << endl << " - Going to remove  " << list.size() << "  simplicies"<<endl<<endl;				
        Tri.RemoveSimplex(list);
        Tri.UpdateAdjacency(list);
    }
    // ========================================================================== //
    // REMOVE ISOLATED VERTICES                                                   //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Finding isolated vertices" << endl;

        ivector1D list;
        // Output stats --------------------------------------------------------- //
        list = Tri.FindIsolatedVertex();
        cout << endl << " - Going to remove  " << list.size() << "  vertices"<<endl<<endl;				
        Tri.RemoveVertex(list);
        Tri.UpdateAdjacency(list);
    }


    // ========================================================================== //
    // REMOVE DOUBLE SIMPLICIES                                                   //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        ivector1D           list;

        // Output message ------------------------------------------------------- //
        cout << endl << " - Removing true double simplicies" << endl << endl;

        // Output stats --------------------------------------------------------- //
        list = Tri.FindTrueDoubleSimplex();
        cout << endl << " - Going to remove  " << list.size() <<"  simplicies"<<endl<<endl;	
        Tri.RemoveSimplex(list);
        Tri.UpdateAdjacency(list);
    }

    Tri.FixNodeNumb();
    Tri.ResizeVertex();
    Tri.ResizeSimplex();
    Tri.ResizeNormal();
    Tri.ResizeAdjacency();

    Tri.Stats(cout);
    // TODO VTK Output

    // ========================================================================== //
    // CLOSING MESSAGE                                                            //
    // ========================================================================== //
    {
        //     // Scope variables ------------------------------------------------------ //
        //     // none
        // 
        //     ivector2D edges, edgesadj;
        //     Tri.BuildEdges(edges,edgesadj);
        // 
        //     int counter = 0;
        //     cout<<"Check every 1-ring in tassellation vertices...."<<endl;
        // 
        //     bvector1D checkV(Tri.nVertex, false);
        // 				
        //    for(int T=0; T<Tri.nSimplex; ++T)
        //    {
        //     for(int j=0; j<Tri.Simplex[0].size(); ++j)
        //     {
        //      bool check,isRing=true;
        //      if(!checkV[Tri.Simplex[T][j]])
        // 	{										
        //          ivector1D list = Tri.Ring_1(T,j,check, isRing);
        // 	 if(!isRing) {	cout<<"Failed Ring Computation for vertex "<<Tri.Simplex[T][j]<<" ring size was"<<list.size()<<endl;
        // 			counter++;								     
        // 	             }	
        //      	checkV[Tri.Simplex[T][j]] = true;
        //   	}	
        //     }
        //   }

        //  cout<<"Total failed Rings caught:   "<<counter; 	
        //  if(counter>0) {cout<<"   Check your triangulation."<<endl;}
        // Output message ------------------------------------------------------- //
        cout << "======================== DEMO: done!! =========================" << endl;
    }

    return; 
};

// ========================================================================== //
int main(
        void
        ) {

    // ========================================================================== //
    // int main(                                                                  //
    //     void)                                                                  //
    //                                                                            //
    // Demo with example of usage for Class_SurfTri variables. Main.              //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //


    // Counters
    // none

    // Output message ------------------------------------------------------- //
    cout << "===================== Class_SurfTri DEMO ===================== " << endl;
    // Scope variables ------------------------------------------------------ //
    // none

    // Run demo ------------------------------------------------------------- //
    Demo_CleaningTools3();

    return 0;
}


