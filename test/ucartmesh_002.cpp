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

        Mesh1.setMesh(B0, B1, nc,3);
        Mesh2 = Mesh1;
        Mesh2.ExportVtr("./","original");

        // Scale mesh ----------------------------------------------------------- //
        darray3E        scale;
        scale.fill(2.0) ;

        Mesh2.Scale(scale);
        Mesh2.ExportVtr("./","scaled");

        // Translate mesh ------------------------------------------------------- //
        darray3E        transl;
        transl[0] = -0.507;    transl[1] = -0.523;    transl[2] = -0.497;
        Mesh2.Translate(transl);
        Mesh2.ExportVtr("./","translated");

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
        int                 n, N ;
        vector<int>         intStencil ;
        vector<double>      intWeights ;
        // Output message ------------------------------------------------------- //
        cout << " - Cell data interpolation" << endl;

        // Intialize cell data on finer mesh ------------------------------------ //
        for( J=0 ; J<Mesh1.getNCells(); ++J){
            P = Mesh1.getCellCenter(J);
            Scalar1[J] =  P[0] + P[1] + P[2] + 1.0 ;
            Vectorial1[J][0] = P[0] + P[1] + P[2] + 1.0;
            Vectorial1[J][1] = P[0] + P[1] + P[2] + 1.0;
            Vectorial1[J][2] = P[0] + P[1] + P[2] + 1.0;
        } //next j

        // Interpolate cell data on the coarser mesh ---------------------------- //
        for (k = 0; k < Mesh2.getNCells(2); k++) {
            for (j = 0; j < Mesh2.getNCells(1); j++) {
                for (i = 0; i < Mesh2.getNCells(0); i++) {
                    P = Mesh2.getCellCenter(i,j,k);

                    N = Mesh1.linearCellInterpolation( P, intStencil, intWeights );

                    Scalar2[J] = 0. ;
                    Vectorial2[J] = {0.,0.,0.} ;
                    for( n=0; n<N; ++n){
                        Scalar2[J] += intWeights[n] * Scalar1[ intStencil[n] ] ;
                        Vectorial2[J] = Vectorial2[J] + intWeights[n] * Vectorial1[ intStencil[n] ] ;
                    };

                } //next k
            } //next j
        } //next i

        // Export vectorial data ------------------------------------------------ //
        Mesh1.ExportVtr("./", "3Doriginal_CSdata", "scalar_field", "Cell", Scalar1);
        Mesh2.ExportVtr("./", "3Dmapped_CSdata", "scalar_field", "Cell", Scalar2);
        Mesh1.ExportVtr("./", "3Doriginal_CVdata", "vectorial_field", "Cell", Vectorial1);
        Mesh2.ExportVtr("./", "3Dmapped_CVdata", "vectorial_field", "Cell", Vectorial2);
    }



    // ========================================================================== //
    // POINT DATA INTERPOLATION                                                   //
    // ========================================================================== //
    {

        // Scope variables ------------------------------------------------------ //
        int                 J, i, j, k;
        dvector1D           Scalar1( Mesh1.getNNodes(), 1.0);
        dvector1D           Scalar2( Mesh2.getNNodes(), 0.0);
        dvecarr3E           Vectorial1( Mesh1.getNNodes(), darray3E() );
        dvecarr3E           Vectorial2( Mesh1.getNNodes(), darray3E() );

        darray3E            P;
        int                 n, N ;
        vector<int>         intStencil ;
        vector<double>      intWeights ;
        // Output message ------------------------------------------------------- //
        cout << " - Point data interpolation" << endl;

        // Intialize cell data on finer mesh ------------------------------------ //
        for (k = 0; k < Mesh1.getNNodes(2); k++) {
            for (j = 0; j < Mesh1.getNNodes(1); j++) {
                for (i = 0; i < Mesh1.getNNodes(0); i++) {
                    P = Mesh1.getNodeCoordinates(i,j,k);
                    J   =   Mesh1.NodeLinearId( i, j, k ) ;

                    Scalar1[J]       = P[0] + P[1] + P[2] + 1.0;
                    Vectorial1[J][0] = P[0] + P[1] + P[2] + 1.0;
                    Vectorial1[J][1] = P[0] + P[1] + P[2] + 1.0;
                    Vectorial1[J][2] = P[0] + P[1] + P[2] + 1.0;
                } //next k
            } //next j
        } //next i

        // Interpolate cell data on the coarser mesh ---------------------------- //
        for (k = 0; k < Mesh2.getNNodes(2); k++) {
            for (j = 0; j < Mesh2.getNNodes(1); j++) {
                for (i = 0; i < Mesh2.getNNodes(0); i++) {
                    P = Mesh2.getNodeCoordinates(i,j,k);

                    N = Mesh1.linearNodeInterpolation( P, intStencil, intWeights );

                    Scalar2[J] = 0. ;
                    Vectorial2[J] = {0.,0.,0.} ;
                    for( n=0; n<N; ++n){
                        Scalar2[J] += intWeights[n] * Scalar1[ intStencil[n] ] ;
                        Vectorial2[J] = Vectorial2[J] + intWeights[n] * Vectorial1[ intStencil[n] ] ;
                    };

                } //next k
            } //next j
        } //next i

        // Export vectorial data ------------------------------------------------ //
        Mesh1.ExportVtr("./", "3Doriginal_PSdata", "scalar_field", "Point", Scalar1);
        Mesh2.ExportVtr("./", "3Dmapped_PSdata", "scalar_field", "Point", Scalar2);
        Mesh1.ExportVtr("./", "3Doriginal_PVdata", "vectorial_field", "Point", Vectorial1);
        Mesh2.ExportVtr("./", "3Dmapped_PVdata", "vectorial_field", "Point", Vectorial2);

    }

    return; 
};
// ========================================================================== //
// MAIN                                                                       //
int main(
        void
        ) {

    Demo3D_UCartMesh(); 

    return(0); 
};

