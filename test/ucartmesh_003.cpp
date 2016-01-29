// ========================================================================== //
//                  - DEMO FOR CARTESIAN MESH MANAGER -                       //
//                                                                            //
// Main                                                                       //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author      :   Alessandro Alaia                                           //
// Version     :   v2.0                                                       //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <iostream>

# include <Class_UCartMesh.hpp>
# include <UCartMesh.hpp>

using namespace std;
// -------------------------------------------------------------------------- //
void Demo_Class_UCartMesh2D(
        void
        ) {

    // ========================================================================== //
    // void Demo_Class_UCartMesh2D(                                               //
    //     void)                                                                  //
    //                                                                            //
    // Example of usage of Class_UCartMesh2D.                                     //
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
    Class_UCartMesh2D           Mesh1, Mesh2;

    // Counters
    // none

    // ========================================================================== //
    // OUTPUT MESSAGE                                                             //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        // none

        // Output message ------------------------------------------------------- //
        cout << endl;
        cout << "================= Class_UCartMesh2D DEMO =================" << endl;
    }

    // ========================================================================== //
    // SETTING MESH & MESH MANIPULATION                                           //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        dvector1D       xlim(2, 0.0), ylim(2, 0.0);
        int             nx = 100, ny = 100;

        // Output message ------------------------------------------------------- //
        cout << " - Setting mesh and mesh manipulation" << endl;

        // Set mesh ------------------------------------------------------------- //
        xlim[0] = 0.0;  xlim[1] = 1.0;
        ylim[0] = 0.0;  ylim[1] = 1.0;
        Mesh1.SetMesh(xlim, ylim, nx, ny);
        Mesh2 = Mesh1;
        // TODO VTK output

        // Scale mesh ----------------------------------------------------------- //
        dvector1D       scale(2, 2.0);
        Mesh2.Scale(scale);
        // TODO VTK output

        // Translate mesh ------------------------------------------------------- //
        dvector1D       transl(2, 0.0);
        transl[0] = -0.507;    transl[1] = -0.523;
        Mesh2.Translate(transl);
        // TODO VTK output

    }

    // ========================================================================== //
    // CELL DATA INTERPOLATION                                                    //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        int                 J, i, j;
        dvector1D           Scalar1(Mesh1.nx * Mesh1.ny, 1.0);
        dvector1D           Scalar2(Mesh2.nx * Mesh2.ny, 0.0);
        dvector2D           Vectorial1(Mesh1.nx * Mesh1.ny, dvector1D(3, -1.0));
        dvector2D           Vectorial2(Mesh2.nx * Mesh2.ny, dvector1D(3, 0.0));
        dvector1D           P(2, 0.0);

        // Output message ------------------------------------------------------- //
        cout << " - Cell data interpolation" << endl;

        // Intialize cell data on finer mesh ------------------------------------ //
        for (i = 0; i < Mesh1.nx; i++) {
            for (j = 0; j < Mesh1.ny; j++) {
                J = Mesh1.AccessCellData(i,j);
                Scalar1[J] = Mesh1.xnode[i] + Mesh1.ynode[j] + 1.0;
                Vectorial1[J][0] = Mesh1.xnode[i] + Mesh1.ynode[j] + 1.0;
                Vectorial1[J][1] = Mesh1.xnode[i] + Mesh1.ynode[j] + 1.0;
                Vectorial1[J][2] = Mesh1.xnode[i] + Mesh1.ynode[j] + 1.0;
            } //next j
        } //next i

        // Interpolate cell data on the coarser mesh ---------------------------- //
        for (i = 0; i < Mesh2.nx; i++) {
            for (j = 0; j < Mesh2.ny; j++) {
                J = Mesh2.AccessCellData(i,j);
                P[0] = Mesh2.xnode[i];
                P[1] = Mesh2.ynode[j];
                Mesh1.interpolateCellData(P, Scalar1, Scalar2[J], Vectorial1, Vectorial2[J]);
            } //next j
        } //next i

        // Export vectorial data ------------------------------------------------ //
        // TODO VTK output

    }

    // ========================================================================== //
    // POINT DATA INTERPOLATION                                                   //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        int                 J, i, j;
        dvector1D           Scalar1((Mesh1.nx+1) * (Mesh1.ny+1), 1.0);
        dvector1D           Scalar2((Mesh2.nx+1) * (Mesh2.ny+1), 0.0);
        dvector2D           Vectorial1((Mesh1.nx+1) * (Mesh1.ny+1), dvector1D(3, -1.0));
        dvector2D           Vectorial2((Mesh2.nx+1) * (Mesh2.ny+1), dvector1D(3, 0.0));
        dvector1D           P(2, 0.0);

        // Output message ------------------------------------------------------- //
        cout << " - Point data interpolation" << endl;

        // Intialize cell data on finer mesh ------------------------------------ //
        for (i = 0; i < Mesh1.nx+1; i++) {
            for (j = 0; j < Mesh1.ny+1; j++) {
                J = Mesh1.AccessPointData(i,j);
                Scalar1[J] = Mesh1.xedge[i] + Mesh1.yedge[j] + 1.0;
                Vectorial1[J][0] = Mesh1.xedge[i] + Mesh1.yedge[j] + 1.0;
                Vectorial1[J][1] = Mesh1.xedge[i] + Mesh1.yedge[j] + 1.0;
                Vectorial1[J][2] = Mesh1.xedge[i] + Mesh1.yedge[j] + 1.0;
            } //next j
        } //next i

        // Interpolate cell data on the coarser mesh ---------------------------- //
        for (i = 0; i < Mesh2.nx+1; i++) {
            for (j = 0; j < Mesh2.ny+1; j++) {
                J = Mesh2.AccessPointData(i,j);
                P[0] = Mesh2.xedge[i];
                P[1] = Mesh2.yedge[j];
                Mesh1.interpolatePointData(P, Scalar1, Scalar2[J], Vectorial1, Vectorial2[J]);
            } //next j
        } //next i

        // Export vectorial data ------------------------------------------------ //
        // TODO VTK output

    }

    return; 
};

// ========================================================================== //
// MAIN                                                                       //
int main(
        void
        ) {

    Demo_Class_UCartMesh2D();

    return(0); 
};


