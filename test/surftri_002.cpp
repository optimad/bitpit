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
void Demo_CleaningTools2(
        void
        ) {

    // ========================================================================== //
    // void Demo_CleaningTools2(                                                  //
    //     void)                                                                  //
    //                                                                            //
    // Demo for Class_SurfTri cleaning tools.                                     //
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

    // Counters
    // none

    // ========================================================================== //
    // OUTPUT MESSAGE                                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << "============ Class_SurfTri: cleaning tools demo 2 =============" << endl;
    }

    // ========================================================================== //
    // LOAD TASSELATION FROM DGF FILE                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl << " - Loading surface tasselation" << endl << endl;

        // Import from dgf file ------------------------------------------------- //
        Tri.Import_stl("./data/buddha.stl",true);

        // Output stats --------------------------------------------------------- //
        //     Tri.Stats(cout);

        // Export to vtk format ------------------------------------------------- //
        //     Tri.Export_vtu("tri_0.vtu");

    }

    // ========================================================================== //
    // REMOVE DUPLICATED VERTICES                                                 //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        Class_SurfTri                                Tri2 = Tri;
        chrono::time_point<chrono::system_clock>     start, end;
        int                                          elapsed_seconds;

        // Output message ------------------------------------------------------- //
        cout << endl << " - Removing duplicated vertices" << endl << endl;

        // Output stats --------------------------------------------------------- //
        start = std::chrono::system_clock::now();
        Tri2.RemoveDoubleVertex();
        Tri2.ResizeVertex();
        end = chrono::system_clock::now();
        elapsed_seconds = chrono::duration_cast<chrono::milliseconds>(end-start).count();
        cout << "elapsed time: " << elapsed_seconds << " ms" << endl;
        cout << "# of vertices:" << endl;
        cout << " before cleaning: " << Tri.nVertex << ", after cleaning: " << Tri2.nVertex << endl;
        cout << "# of simplicies:" << endl;
        cout << " before cleaning: " << Tri.nSimplex << ", after cleaning: " << Tri2.nSimplex << endl;;
        cout << " double vertices found after cleaning: " << Tri2.CountDoubleVertex() << endl;

        // Measure error w.r.t original surface --------------------------------- //
        double err = 0.0;
        int    m;
        for (int T = 0; T < Tri.nSimplex; T++) {
            m = Tri.Simplex[T].size();
            for (int i = 0; i < m; i++) {
                err += norm2(Tri.Vertex[Tri.Simplex[T][i]] - Tri2.Vertex[Tri2.Simplex[T][i]]);
            } //next i
        } //next T
        cout << "error: " << err << endl;

        // Output stats --------------------------------------------------------- //
        //     Tri.Stats(cout);

        // Export to vtk format ------------------------------------------------- //
        //     Tri.Export_vtu("tri_1.vtu");
        //      Tri2.Export_vtu("tri_1b.vtu");
    }

    // ========================================================================== //
    // BUILD ADJACENCIES                                                          //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        Class_SurfTri                                Tri2 = Tri;
        chrono::time_point<chrono::system_clock>     start, end;
        int                                          elapsed_seconds;

        // Output message ------------------------------------------------------- //
        cout << endl << " - Building adjacencies" << endl << endl;

        // Output stats --------------------------------------------------------- //
        start = std::chrono::system_clock::now();
        Tri2.BuildAdjacency();
        end = chrono::system_clock::now();
        elapsed_seconds = chrono::duration_cast<chrono::milliseconds>(end-start).count();
        cout << "elapsed time: " << elapsed_seconds << " ms" << endl;

    }

    // ========================================================================== //
    // CLOSING MESSAGE                                                            //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

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
    Demo_CleaningTools2();

    return 0;
}


