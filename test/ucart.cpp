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

// CC_lib
// none


// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;

// ========================================================================== //
// MAIN                                                                       //
// ========================================================================== //
int main(
        void
        ) {

    // ========================================================================== //
    // int main(                                                                  //
    //     void)                                                                  //
    //                                                                            //
    // Main function for class_UCartMesh2D demo.                                  //
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
    int         selection;

    // Counters
    // none

    // ========================================================================== //
    // OUTPUT MESSAGE                                                             //
    // ========================================================================== //
    cout << "DEMO for Class_UCartMesh"  << endl;
    cout << "Select demo:"              << endl;
    cout << "0. 2D case"                << endl;
    cout << "1. 3D case"                << endl;
    cout << "2. UCartMesh 2D"                << endl;
    cout << "3. UCartMesh 3D"                << endl;
    //cin >> selection;
    selection =2 ;

    // ========================================================================== //
    //  CALL DEMO FUNCTION                                                        //
    // ========================================================================== //
    Demo_Class_UCartMesh2D(); 
    Demo_Class_UCartMesh3D(); 
    Demo2D_UCartMesh(); 
    Demo3D_UCartMesh(); 

    return(0); 
};
// ========================================================================== //
//                  - DEMO FOR CARTESIAN MESH MANAGER -                       //
//                                                                            //
// Demo functions for Class_UCartMesh usage.                                  //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author      :   Alessandro Alaia                                           //
// Version     :   v2.0                                                       //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //
// -------------------------------------------------------------------------- //

void Demo2D_UCartMesh(
        void
        ){

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
    UCartMesh           Mesh1, Mesh2;

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
        cout << "================= UCartMesh 2D DEMO =================" << endl;
    }

    // ========================================================================== //
    // SETTING MESH & MESH MANIPULATION                                           //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        darray3E        B0, B1;
        iarray3E        nc ;

        // Output message ------------------------------------------------------- //
        cout << " - Setting mesh and mesh manipulation" << endl;

        // Set mesh ------------------------------------------------------------- //
        B0.fill(0.) ; 
        B1.fill(1.) ; 
        nc.fill(100) ; 

        Mesh1.SetMesh(B0, B1, nc, 2);
        Mesh2 = Mesh1;
        Mesh2.Export_vtr("original.vtr");

        // Scale mesh ----------------------------------------------------------- //
        darray3E        scale;
        scale.fill(2.0) ;

        Mesh2.Scale(scale);
        Mesh2.Export_vtr("scaled.vtr");

        // Translate mesh ------------------------------------------------------- //
        darray3E        transl;
        transl[0] = -0.507;    transl[1] = -0.523;    transl[2] = 0.0 ;
        Mesh2.Translate(transl);
        Mesh2.Export_vtr("translated.vtr");    

    }

    // ========================================================================== //
    // CELL DATA INTERPOLATION                                                    //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        int                 J;
        dvector1D           Scalar1( Mesh1.getNCells(), 1.0);
        dvector1D           Scalar2( Mesh2.getNCells(), 0.0);
        dvecarr3E           Vectorial1( Mesh1.getNCells(), darray3E() );
        dvecarr3E           Vectorial2( Mesh2.getNCells(), darray3E() );
        darray3E            P;

        P[2] = 0.0 ;

        // Output message ------------------------------------------------------- //
        cout << " - Cell data interpolation" << endl;

        // Intialize cell data on finer mesh ------------------------------------ //
        for( J=0 ; J<Mesh1.getNCells(); ++J){
            P = Mesh1.getCenter(J) ;
            Scalar1[J] =  P[0] + P[1] +1.0;
            Vectorial1[J][0] = P[0] + P[1] + 1.0;
            Vectorial1[J][1] = P[0] + P[1] + 1.0;
            Vectorial1[J][2] = 0.0 ;
        } //next j

        // Interpolate cell data on the coarser mesh ---------------------------- //
        for (J = 0; J < Mesh2.getNCells(); J++) {
            P = Mesh2.getCenter( J ) ;
            Mesh1.interpolateCellData(P, Scalar1, Scalar2[J], Vectorial1, Vectorial2[J]);
        } //next j

        // Export vectorial data ------------------------------------------------ //
        Mesh1.Export_CellData_vtr("2Doriginal_CSdata.vtr" , "scalar_field" , Scalar1);
        Mesh2.Export_CellData_vtr("2Dmapped_CSdata.vtr" , "scalar_field" , Scalar2);
        Mesh1.Export_CellData_vtr("2Doriginal_CVdata.vtr" , "vectorial_field" , Vectorial1);
        Mesh2.Export_CellData_vtr("2Dmapped_CVdata.vtr" , "vectorial_field" , Vectorial2);
    }



    // ========================================================================== //
    // POINT DATA INTERPOLATION                                                   //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        int                 J ;
        dvector1D           Scalar1( Mesh1.getNPoints(), 1.0);
        dvector1D           Scalar2( Mesh2.getNPoints(), 0.0);
        dvecarr3E           Vectorial1( Mesh1.getNPoints(), darray3E() );
        dvecarr3E           Vectorial2( Mesh1.getNPoints(), darray3E() );
        darray3E            P;

        // Output message ------------------------------------------------------- //
        cout << " - Point data interpolation" << endl;

        // Intialize cell data on finer mesh ------------------------------------ //
        for (J = 0; J < Mesh1.getNPoints(); J++) {

            P = Mesh1.getPoint(J);

            Scalar1[J]       = P[0] + P[1] + 1.0;

            Vectorial1[J][0] = P[0] + P[1] + 1.0;
            Vectorial1[J][1] = P[0] + P[1] + 1.0;
            Vectorial1[J][2] = 0.0 ;
        } //next j

        // Interpolate cell data on the coarser mesh ---------------------------- //
        for (J = 0; J < Mesh2.getNPoints(); J++) {

            P = Mesh2.getPoint(J);

            Mesh1.interpolatePointData(P, Scalar1, Scalar2[J], Vectorial1, Vectorial2[J]);
        } //next j

        // Export vectorial data ------------------------------------------------ //
        Mesh1.Export_PointData_vtr("2Doriginal_PSdata.vtr" , "scalar_field" , Scalar1);
        Mesh2.Export_PointData_vtr("2Dmapped_PSdata.vtr" , "scalar_field" , Scalar2);
        Mesh1.Export_PointData_vtr("2Doriginal_PVdata.vtr" , "vectorial_field" , Vectorial1);
        Mesh2.Export_PointData_vtr("2Dmapped_PVdata.vtr" , "vectorial_field" , Vectorial2);

    }

    return; 
}; 

void Demo3D_UCartMesh(
        void
        ){

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
    UCartMesh           Mesh1, Mesh2;

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
        cout << "================= UCartMesh 3D DEMO =================" << endl;
    }

    // ========================================================================== //
    // SETTING MESH & MESH MANIPULATION                                           //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        darray3E        B0, B1;
        iarray3E        nc ;

        // Output message ------------------------------------------------------- //
        cout << " - Setting mesh and mesh manipulation" << endl;

        // Set mesh ------------------------------------------------------------- //
        B0.fill(0.) ; 
        B1.fill(1.) ; 
        nc.fill(100) ; 

        Mesh1.SetMesh(B0, B1, nc,3);
        Mesh2 = Mesh1;
        Mesh2.Export_vtr("original.vtr");

        // Scale mesh ----------------------------------------------------------- //
        darray3E        scale;
        scale.fill(2.0) ;

        Mesh2.Scale(scale);
        Mesh2.Export_vtr("scaled.vtr");

        // Translate mesh ------------------------------------------------------- //
        darray3E        transl;
        transl[0] = -0.507;    transl[1] = -0.523;    transl[2] = -0.497;
        Mesh2.Translate(transl);
        Mesh2.Export_vtr("translated.vtr");    

    }

    // ========================================================================== //
    // CELL DATA INTERPOLATION                                                    //
    // ========================================================================== //
    {

        cout << Mesh1.getNCells() << endl ;
        cout << Mesh2.getNCells() << endl ;

        // Scope variables ------------------------------------------------------ //
        int                 J, i, j, k;
        dvector1D           Scalar1( Mesh1.getNCells(), 1.0);
        dvector1D           Scalar2( Mesh2.getNCells(), 0.0);
        dvecarr3E           Vectorial1( Mesh1.getNCells(), darray3E() );
        dvecarr3E           Vectorial2( Mesh2.getNCells(), darray3E() );
        darray3E            P;

        // Output message ------------------------------------------------------- //
        cout << " - Cell data interpolation" << endl;

        // Intialize cell data on finer mesh ------------------------------------ //
        for( J=0 ; J<Mesh1.getNCells(); ++J){
            P = Mesh1.getCenter(J);
            Scalar1[J] =  P[0] + P[1] + P[2] + 1.0 ;
            Vectorial1[J][0] = P[0] + P[1] + P[2] + 1.0;
            Vectorial1[J][1] = P[0] + P[1] + P[2] + 1.0;
            Vectorial1[J][2] = P[0] + P[1] + P[2] + 1.0;
        } //next j

        // Interpolate cell data on the coarser mesh ---------------------------- //
        for (k = 0; k < Mesh2.getNCells(2); k++) {
            for (j = 0; j < Mesh2.getNCells(1); j++) {
                for (i = 0; i < Mesh2.getNCells(0); i++) {
                    P = Mesh2.getCenter(i,j,k);
                    Mesh1.interpolateCellData(P, Scalar1, Scalar2[J], Vectorial1, Vectorial2[J]);
                } //next k
            } //next j
        } //next i

        // Export vectorial data ------------------------------------------------ //
        Mesh1.Export_CellData_vtr("3Doriginal_CSdata.vtr" , "scalar_field" , Scalar1);
        Mesh2.Export_CellData_vtr("3Dmapped_CSdata.vtr" , "scalar_field" , Scalar2);
        Mesh1.Export_CellData_vtr("3Doriginal_CVdata.vtr" , "vectorial_field" , Vectorial1);
        Mesh2.Export_CellData_vtr("3Dmapped_CVdata.vtr" , "vectorial_field" , Vectorial2);
    }



    // ========================================================================== //
    // POINT DATA INTERPOLATION                                                   //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        int                 J, i, j, k;
        dvector1D           Scalar1( Mesh1.getNPoints(), 1.0);
        dvector1D           Scalar2( Mesh2.getNPoints(), 0.0);
        dvecarr3E           Vectorial1( Mesh1.getNPoints(), darray3E() );
        dvecarr3E           Vectorial2( Mesh1.getNPoints(), darray3E() );
        darray3E            P;

        // Output message ------------------------------------------------------- //
        cout << " - Point data interpolation" << endl;

        // Intialize cell data on finer mesh ------------------------------------ //
        for (k = 0; k < Mesh1.getNPoints(2); k++) {
            for (j = 0; j < Mesh1.getNPoints(1); j++) {
                for (i = 0; i < Mesh1.getNPoints(0); i++) {
                    P = Mesh1.getPoint(i,j,k);
                    J   =   Mesh1.AccessPoint( i, j, k ) ;

                    Scalar1[J]       = P[0] + P[1] + P[2] + 1.0;
                    Vectorial1[J][0] = P[0] + P[1] + P[2] + 1.0;
                    Vectorial1[J][1] = P[0] + P[1] + P[2] + 1.0;
                    Vectorial1[J][2] = P[0] + P[1] + P[2] + 1.0;
                } //next k
            } //next j
        } //next i

        // Interpolate cell data on the coarser mesh ---------------------------- //
        for (k = 0; k < Mesh2.getNPoints(2); k++) {
            for (j = 0; j < Mesh2.getNPoints(1); j++) {
                for (i = 0; i < Mesh2.getNPoints(0); i++) {
                    P = Mesh2.getPoint(i,j,k);
                    Mesh1.interpolatePointData(P, Scalar1, Scalar2[J], Vectorial1, Vectorial2[J]);
                } //next k
            } //next j
        } //next i

        // Export vectorial data ------------------------------------------------ //
        Mesh1.Export_PointData_vtr("3Doriginal_PSdata.vtr" , "scalar_field" , Scalar1);
        Mesh2.Export_PointData_vtr("3Dmapped_PSdata.vtr" , "scalar_field" , Scalar2);
        Mesh1.Export_PointData_vtr("3Doriginal_PVdata.vtr" , "vectorial_field" , Vectorial1);
        Mesh2.Export_PointData_vtr("3Dmapped_PVdata.vtr" , "vectorial_field" , Vectorial2);

    }

    return; 
};
// ========================================================================== //
//                  - DEMO FOR CARTESIAN MESH MANAGER -                       //
//                                                                            //
// Demo functions for Class_UCartMesh usage.                                  //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author      :   Alessandro Alaia                                           //
// Version     :   v2.0                                                       //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

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
        Mesh2.Export_vtr("original.vtr");

        // Scale mesh ----------------------------------------------------------- //
        dvector1D       scale(2, 2.0);
        Mesh2.Scale(scale);
        Mesh2.Export_vtr("scaled.vtr");

        // Translate mesh ------------------------------------------------------- //
        dvector1D       transl(2, 0.0);
        transl[0] = -0.507;    transl[1] = -0.523;
        Mesh2.Translate(transl);
        Mesh2.Export_vtr("translated.vtr");    

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
        Mesh1.Export_CellData_vtr("original_CSdata.vtr" , "scalar_field" , Scalar1);
        Mesh2.Export_CellData_vtr("mapped_CSdata.vtr" , "scalar_field" , Scalar2);
        Mesh1.Export_CellData_vtr("original_CVdata.vtr" , "vectorial_field" , Vectorial1);
        Mesh2.Export_CellData_vtr("mapped_CVdata.vtr" , "vectorial_field" , Vectorial2);

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
        Mesh1.Export_PointData_vtr("original_PSdata.vtr" , "scalar_field" , Scalar1);
        Mesh2.Export_PointData_vtr("mapped_PSdata.vtr" , "scalar_field" , Scalar2);
        Mesh1.Export_PointData_vtr("original_PVdata.vtr" , "vectorial_field" , Vectorial1);
        Mesh2.Export_PointData_vtr("mapped_PVdata.vtr" , "vectorial_field" , Vectorial2);

    }

    return; 
};

// -------------------------------------------------------------------------- //
void Demo_Class_UCartMesh3D(
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
    Class_UCartMesh3D           Mesh1, Mesh2;

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
        cout << "================= Class_UCartMesh3D DEMO =================" << endl;
    }

    // ========================================================================== //
    // SETTING MESH & MESH MANIPULATION                                           //
    // ========================================================================== //
    {
        // Scope variables ------------------------------------------------------ //
        dvector1D       xlim(2, 0.0), ylim(2, 0.0), zlim(2, 0.0);
        int             nx = 100, ny = 100, nz = 100;

        // Output message ------------------------------------------------------- //
        cout << " - Setting mesh and mesh manipulation" << endl;

        // Set mesh ------------------------------------------------------------- //
        xlim[0] = 0.0;  xlim[1] = 1.0;
        ylim[0] = 0.0;  ylim[1] = 1.0;
        zlim[0] = 0.0;  zlim[1] = 1.0;
        Mesh1.SetMesh(xlim, ylim, zlim, nx, ny, nz);
        Mesh2 = Mesh1;
        Mesh2.Export_vtr("original.vtr");

        // Scale mesh ----------------------------------------------------------- //
        dvector1D       scale(3, 2.0);
        Mesh2.Scale(scale);
        Mesh2.Export_vtr("scaled.vtr");

        // Translate mesh ------------------------------------------------------- //
        dvector1D       transl(3, 0.0);
        transl[0] = -0.507;    transl[1] = -0.523;    transl[2] = -0.497;
        Mesh2.Translate(transl);
        Mesh2.Export_vtr("translated.vtr");    

    }

    // ========================================================================== //
    // CELL DATA INTERPOLATION                                                    //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        int                 J, i, j, k;
        dvector1D           Scalar1(Mesh1.nx * Mesh1.ny * Mesh1.nz, 1.0);
        dvector1D           Scalar2(Mesh2.nx * Mesh2.ny * Mesh2.nz, 0.0);
        dvector2D           Vectorial1(Mesh1.nx * Mesh1.ny * Mesh1.nz, dvector1D(3, -1.0));
        dvector2D           Vectorial2(Mesh2.nx * Mesh2.ny * Mesh2.nz, dvector1D(3, 0.0));
        dvector1D           P(3, 0.0);

        // Output message ------------------------------------------------------- //
        cout << " - Cell data interpolation" << endl;

        // Intialize cell data on finer mesh ------------------------------------ //
        for (i = 0; i < Mesh1.nx; i++) {
            for (j = 0; j < Mesh1.ny; j++) {
                for (k = 0; k < Mesh1.nz; k++) {
                    J = Mesh1.AccessCellData(i,j,k);
                    Scalar1[J] = Mesh1.znode[k] + Mesh1.xnode[i] + Mesh1.ynode[j] + 1.0;
                    Vectorial1[J][0] = Mesh1.znode[k] + Mesh1.xnode[i] + Mesh1.ynode[j] + 1.0;
                    Vectorial1[J][1] = Mesh1.znode[k] + Mesh1.xnode[i] + Mesh1.ynode[j] + 1.0;
                    Vectorial1[J][2] = Mesh1.znode[k] + Mesh1.xnode[i] + Mesh1.ynode[j] + 1.0;
                } //next k
            } //next j
        } //next i

        // Interpolate cell data on the coarser mesh ---------------------------- //
        for (i = 0; i < Mesh2.nx; i++) {
            for (j = 0; j < Mesh2.ny; j++) {
                for (k = 0; k < Mesh2.nz; k++) {
                    J = Mesh2.AccessCellData(i,j,k);
                    P[0] = Mesh2.xnode[i];
                    P[1] = Mesh2.ynode[j];
                    P[2] = Mesh2.znode[k];
                    Mesh1.interpolateCellData(P, Scalar1, Scalar2[J], Vectorial1, Vectorial2[J]);
                } //next k
            } //next j
        } //next i

        // Export vectorial data ------------------------------------------------ //
        // Mesh1.Export_CellData_vtr("original_CSdata.vtr" , "scalar_field" , Scalar1);
        // Mesh2.Export_CellData_vtr("mapped_CSdata.vtr" , "scalar_field" , Scalar2);
        // Mesh1.Export_CellData_vtr("original_CVdata.vtr" , "vectorial_field" , Vectorial1);
        // Mesh2.Export_CellData_vtr("mapped_CVdata.vtr" , "vectorial_field" , Vectorial2);

    }

    // ========================================================================== //
    // POINT DATA INTERPOLATION                                                   //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        int                 J, i, j, k;
        dvector1D           Scalar1((Mesh1.nx+1) * (Mesh1.ny+1) * (Mesh1.nz+1), 1.0);
        dvector1D           Scalar2((Mesh2.nx+1) * (Mesh2.ny+1) * (Mesh1.nz+1), 0.0);
        dvector2D           Vectorial1((Mesh1.nx+1) * (Mesh1.ny+1) * (Mesh1.nz+1), dvector1D(3, -1.0));
        dvector2D           Vectorial2((Mesh2.nx+1) * (Mesh2.ny+1) * (Mesh1.nz+1), dvector1D(3, 0.0));
        dvector1D           P(3, 0.0);

        // Output message ------------------------------------------------------- //
        cout << " - Point data interpolation" << endl;

        // Intialize cell data on finer mesh ------------------------------------ //
        for (i = 0; i < Mesh1.nx+1; i++) {
            for (j = 0; j < Mesh1.ny+1; j++) {
                for (k = 0; k < Mesh1.nz+1; k++) {
                    J = Mesh1.AccessPointData(i,j,k);
                    Scalar1[J] = Mesh1.zedge[k] + Mesh1.xedge[i] + Mesh1.yedge[j] + 1.0;
                    Vectorial1[J][0] = Mesh1.zedge[k] + Mesh1.xedge[i] + Mesh1.yedge[j] + 1.0;
                    Vectorial1[J][1] = Mesh1.zedge[k] + Mesh1.xedge[i] + Mesh1.yedge[j] + 1.0;
                    Vectorial1[J][2] = Mesh1.zedge[k] + Mesh1.xedge[i] + Mesh1.yedge[j] + 1.0;
                } //next k
            } //next j
        } //next i

        // Interpolate cell data on the coarser mesh ---------------------------- //
        for (i = 0; i < Mesh2.nx+1; i++) {
            for (j = 0; j < Mesh2.ny+1; j++) {
                for (k = 0; k < Mesh2.nz+1; k++) {
                    J = Mesh2.AccessPointData(i,j,k);
                    P[0] = Mesh2.xedge[i];
                    P[1] = Mesh2.yedge[j];
                    P[2] = Mesh2.zedge[k];
                    Mesh1.interpolatePointData(P, Scalar1, Scalar2[J], Vectorial1, Vectorial2[J]);
                } //next k
            } //next j
        } //next i

        // Export vectorial data ------------------------------------------------ //
        // Mesh1.Export_PointData_vtr("original_PSdata.vtr" , "scalar_field" , Scalar1);
        // Mesh2.Export_PointData_vtr("mapped_PSdata.vtr" , "scalar_field" , Scalar2);
        // Mesh1.Export_PointData_vtr("original_PVdata.vtr" , "vectorial_field" , Vectorial1);
        // Mesh2.Export_PointData_vtr("mapped_PVdata.vtr" , "vectorial_field" , Vectorial2);

    }

    return; 
};
