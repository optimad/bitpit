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

// ========================================================================== //
// MAIN                                                                       //
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
    Demo2D_UCartMesh(); 

    return(0); 
};


