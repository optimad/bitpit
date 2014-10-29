// ========================================================================== //
//                         - Class_VolTri -                                   //
//                                                                            //
// Grid manager for unstructured volume meshes.                               //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author   : Alessandro Alaia                                                //
// Version  : v2.0                                                            //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
# include "Class_VolTri.hpp"

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //

// TOPOLOGY ================================================================= //

// -------------------------------------------------------------------------- //
bool Class_VolTri::SameFace(
    int             A,
    int             i,
    int             B,
    int             j
) {

// ========================================================================== //
// bool Class_VolTri::SameFace(                                               //
//     int             A,                                                     //
//     int             i,                                                     //
//     int             B,                                                     //
//     int             j)                                                     //
//                                                                            //
// Check for coincident faces.                                                //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - A     : int, 1st simplex global index                                    //
// - i     : int, face local index on simplex A                               //
// - B     : int, 2nd simplex global index                                    //
// - j     : int, face local index on simplex B                               //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - check : bool, returns true if face (A, i) and (B, j) are coincident      //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                check = false;
ivector1D           face_vlistA, face_vlistB;

// ========================================================================== //
// CHECK FOR COINCIDENT FACES                                                 //
// ========================================================================== //
face_vlistA = FaceVertices(A, i);
face_vlistB = FaceVertices(B, j);
if (face_vlistA.size() == face_vlistB.size()) {
    sort(face_vlistA.begin(), face_vlistA.end());
    sort(face_vlistB.begin(), face_vlistB.end());
    check = (face_vlistA == face_vlistB);
}

return(check); };

// -------------------------------------------------------------------------- //
bool Class_VolTri::SameSimplex(
    int             A,
    int             B
) {

// ========================================================================== //
// bool Class_VolTri::SameSimplex(                                            //
//     int             A,                                                     //
//     int             B)                                                     //
//                                                                            //
// Check if two simplicies are coincident.                                    //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - A     : int, 1st simplex global index                                    //
// - B     : int, 2nd simplex global index                                    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - check : bool, true if A and B are coincident, false otherwise            //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool        check = false;

// Counters
// none

// ========================================================================== //
// CHECK FOR COINCIDENT SIMPLICIES                                            //
// ========================================================================== //
if (e_type[A] == e_type[B]) {

    // Scope variables
    ivector1D       vlist_A(infos[e_type[A]].n_vert);
    ivector1D       vlist_B(infos[e_type[B]].n_vert);

    // Sort vertex list in ascending order
    vlist_A = Simplex[A];
    sort(vlist_A.begin(), vlist_A.end());
    vlist_B = Simplex[B];
    sort(vlist_A.begin(), vlist_B.end());
    check = (vlist_A == vlist_B);
    
}

return(check); };

// -------------------------------------------------------------------------- //
ivector1D Class_VolTri::FaceVertices(
    int             S,
    int             i
) {

// ========================================================================== //
// ivector1D Class_VolTri::FaceVertices(                                      //
//     int             S,                                                     //
//     int             i)                                                     //
//                                                                            //
// Return list of global indices of vertices for the i-th face of simplex S.  //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - S     : int, simplex global index                                        //
// - i     : int, face local index                                            //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - vlist : ivector1D, list of global indices of vertices of the i-th face   //
//           of simplex S.                                                    //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                 n = infos[e_type[S]].faces[i].size();
ivector1D           vlist(n, -1);

// Counters
int                 j;

// ========================================================================== //
// RETURN GLOBAL INDEX OF VERTICES                                            //
// ========================================================================== //
for (j = 0; j < n; ++j) {
    vlist[j] = Simplex[S][infos[e_type[S]].faces[i][j]];
} //next i

return(vlist); };

// GEOMETRY ================================================================= //

// -------------------------------------------------------------------------- //
dvector1D Class_VolTri::Baricenter(
    int             T
) {

// ========================================================================== //
// dvector1D Class_VolTri::Baricenter(                                        //
//     int             T)                                                     //
//                                                                            //
// Compute simplex baricenter.                                                //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T       : int, simplex global index                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - P       : dvector1D, with simplex baricenter coordinates.                //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                     dim = Vertex[0].size(), n = infos[e_type[T]].n_vert;
dvector1D               P(dim, 0.0);

// Counters
int                     i, j;

// ========================================================================== //
// FIND ISOLATED VERTEXES                                                     //
// ========================================================================== //
for (i = 0; i < n; ++i) {
    for (j = 0; j < dim; ++j) {
        P[j] += Vertex[Simplex[T][i]][j];
    } //next j
} //next i
for (j = 0; j < dim; ++j) {
    P[j] = P[j]/((double) n);
} //next j

return(P); };

// -------------------------------------------------------------------------- //
dvector1D Class_VolTri::FaceCenter(
    int             T,
    int             i
) {

// ========================================================================== //
// dvector1D Class_VolTri::FaceCenter(                                        //
//     int             T,                                                     //
//     int             i)                                                     //
//                                                                            //
// Returns baricenter of the i-th face of the T-th simplex.                   //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex globa index                                        //
// - i      : int, face local index                                           //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - C      : dvector1D, with face center coordinates                         //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int              dim = Vertex[0].size();
int              n;
ivector1D        f_vert;
dvector1D        C(dim, 0.0);

// Counters
int              V;
int              j, k;

// ========================================================================== //
// COMPUTE FACE CENTER                                                        //
// ========================================================================== //

// Get face vertex list
f_vert = FaceVertices(T, i);
n = f_vert.size();
for (j = 0; j < n; ++j) {
    V = f_vert[j];
    for (k = 0; k < dim; ++k) {
        C[k] += Vertex[V][k];
    } //next k
} //next j
for (j = 0; j < dim; ++j) {
    C[j] = C[j]/((double) n);
} //next j

return(C); };

// -------------------------------------------------------------------------- //
dvector1D Class_VolTri::CircumCenter(
    int             T
) {

// ========================================================================== //
// dvector1D Class_VolTri::CircumCenter(                                      //
//     int             T)                                                     //
//                                                                            //
// Compute circum center of simplex T.                                        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex global index                                       //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
int     const   iter_max = 10;
double  const   abs_toll = 1.0e-6;

// Local variables
bool            check;
int             dim = Vertex[0].size();
dvector1D       C(dim, 0.0), dC;
dvector1D       d, r;
dvector1D       dG;
dvector2D       dr;
dvector3D       d2r;
dvector2D       d2G;

// Counters
int             i, j, k;
int             n, m;
int             U, V;
int             iter;

// ========================================================================== //
// SET PARAMETERS                                                             //
// ========================================================================== //

// Linear system dimensions
if ((e_type[T] == 5) || (e_type[T] == 9)) {
    n = 2;
}
else {
    n = 3;
}
dG.resize(n, 0.0);
d2G.resize(n, dvector1D(n, 0.0));
dC.resize(n, 0.0);

// number of simplex vertices
m = infos[e_type[T]].n_vert;
d.resize(m, 0.0);
r.resize(m, 0.0);
dr.resize(m, dvector1D(n, 0.0));
d2r.resize(m, dvector2D(n, dvector1D(n, 0.0)));

// ========================================================================== //
// ITERATIVE PROCEDURE                                                        //
// ========================================================================== //

// Initial estimate for circum center --------------------------------------- //
C = Baricenter(T);

// Newton-Raphson iterations ------------------------------------------------ //
iter = 0;
check = true;
U = Simplex[T][0];
while (check && (iter < iter_max)) {

    // Reset r.h.s. and coeffs. matrix
    for (i = 0; i < n; ++i) {
        dG[i] = 0.0;
        for (j = 0; j < n; ++j) {
            d2G[i][j] = 0.0;
        } //next j
    } //next i

    // Update distances
    for (i = 0; i < m; ++i) {
        V = Simplex[T][i];
        d[i] = norm_2(C - Vertex[V]);
    } //next i

    // Update residuals
    for (i = 1; i < m; ++i) {
        r[i] = d[i] - d[0];
    } //next i

    // Update gradient of residuals
    for (i = 1; i < m; ++i) {
        V = Simplex[T][i];
        dr[i] = (C - Vertex[V])/d[i] - (C - Vertex[U])/d[0];
    } //next i

    // Update Hessians of residuals
    for (i = 1; i < m; ++i) {
        V = Simplex[T][i];
        for (j = 0; j < n; ++j) {
            for (k = 0; k < n; ++k) {
                d2r[i][j][k] = -(C[j] - Vertex[V][j]) * (C[k] - Vertex[V][k]) / (d[i]*d[i]*d[i])
                               +(C[j] - Vertex[U][j]) * (C[k] - Vertex[U][k]) / (d[0]*d[0]*d[0]);
            } //next k
        } //next j
        for (j = 0; j < n; ++j) {
            d2r[i][j][j] += 1.0/d[i] - 1.0/d[0];
        } //next j
    }

    // r.h.s
    for (i = 1; i < m; ++i) {
        V = Simplex[T][i];
        dG = dG + 2.0 * r[i] * dr[i];
    } //next i
    dG = -1.0*dG;

    // coeff. matrix
    for (i = 0; i < n; ++i) {    
        for (j = 0; j < n; ++j) {
            for (k = 1; k < m; ++k) {
                d2G[i][j] += 2.0 * (dr[k][i]*dr[k][j] + r[i] * d2r[k][i][j]);
            } //next k
        } //next j
    } //next i

    // Solve linear system
    Cramer(d2G, dG, dC);

    // Convergence criterion
    check = (norm_2(dC) > abs_toll);
    iter++;

    // Update solution
    for (i = 0; i < n; ++i) {
        C[i] += dC[i];
    } //next i

} //next iter

return(C); };

        // ----------------------------------------------------------------------------------- //
        double Class_VolTri::FaceArea(int T, int i) {

        // =================================================================================== //
        // double Class_VolTri::FaceArea(int T, int i)                                         //
        //                                                                                     //
        // Compute area of the i-th face of the T-th simplex.                                  //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - T    : int, simplex global index                                                  //
        // - i    : int, face local index.                                                     //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - A    : double, face area                                                          //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        int              dim = Vertex[0].size();
        double           A = 0.0;

        // Counters
        int              j;

        // =================================================================================== //
        // COMPUTE FACE AREA                                                                   //
        // =================================================================================== //
        switch(e_type[T]) {
            case 6  : goto label_tria;  break;
            case 8  : goto label_quad;  break;
            case 12 : goto label_tetra; break;
            case 16 : goto label_pyram; break;
            case 18 : goto label_prism; break;
            case 22 : goto label_dhexa; break;
            case 24 : goto label_rhexa; break;

        };

label_tria: {
        j = (i + 1) % 3;
        A = norm_2(Vertex[Simplex[T][i]] - Vertex[Simplex[T][j]]);
        goto label_quit;
        }
label_quad: {
        j = (i + 1) % 4;
        A = norm_2(Vertex[Simplex[T][i]] - Vertex[Simplex[T][j]]);
        goto label_quit;
        }
label_tetra: {
        goto label_quit;
        }
label_pyram: {
        goto label_quit;
        }
label_prism: {
        goto label_quit;
        }
label_dhexa: {
        goto label_quit;
        }
label_rhexa: {
        goto label_quit;
        }

label_quit:
        return(A); };

        // ----------------------------------------------------------------------------------- //
        double Class_VolTri::minFaceArea(int T, int &i) {

        // =================================================================================== //
        // double Class_VolTri::minFaceArea(int T, int i)                                      //
        //                                                                                     //
        // Find the simplex face with the minimal area.                                        //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - T      : int, simplex globa index                                                 //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - A      : double, min face area                                                    //
        // - i      : int, local index of min area face                                        //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        double      A, B;

        // Counters
        int         j;

        // =================================================================================== //
        // FIND MIN FACE AREA                                                                  //
        // =================================================================================== //
        switch(e_type[T]) {
            case 6  : goto label_tria;  break;
            case 8  : goto label_quad;  break;
            case 12 : goto label_tetra; break;
            case 16 : goto label_pyram; break;
            case 18 : goto label_prism; break;
            case 22 : goto label_dhexa; break;
            case 24 : goto label_rhexa; break;
        };
label_tria: {
        A = FaceArea(T, 0);
        i = 0;
        for (j = 1; j < 3; j++) {
            B = FaceArea(T, j);
            if (B < A) {
                A = B;
                i = j;
            }
        } //next i
        goto label_quit;
}
label_quad: {
        A = FaceArea(T, 0);
        i = 0;
        for (j = 1; j < 4; j++) {
            B = FaceArea(T, j);
            if (B < A) {
                A = B;
                i = j;
            }
        } //next i
        goto label_quit;
}
label_tetra: {
        goto label_quit;
}
label_pyram: {
        goto label_quit;
}
label_prism: {
        goto label_quit;
}
label_dhexa: {
        goto label_quit;
}
label_rhexa: {
        goto label_quit;
}

label_quit:
        return(A); };

        // ----------------------------------------------------------------------------------- //
        double Class_VolTri::maxFaceArea(int T, int &i) {

        // =================================================================================== //
        // double Class_VolTri::maxFaceArea(int T, int i)                                      //
        //                                                                                     //
        // Find the simplex face with the maximal area.                                        //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - T      : int, simplex globa index                                                 //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - A      : double, min face area                                                    //
        // - i      : int, local index of max area face                                        //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        double      A, B;

        // Counters
        int         j;

        // =================================================================================== //
        // FIND MIN FACE AREA                                                                  //
        // =================================================================================== //
        switch(e_type[T]) {
            case 6  : goto label_tria;  break;
            case 8  : goto label_quad;  break;
            case 12 : goto label_tetra; break;
            case 16 : goto label_pyram; break;
            case 18 : goto label_prism; break;
            case 22 : goto label_dhexa; break;
            case 24 : goto label_rhexa; break;
        };
label_tria: {
        A = FaceArea(T, 0);
        i = 0;
        for (j = 1; j < 3; j++) {
            B = FaceArea(T, j);
            if (B > A) {
                A = B;
                i = j;
            }
        } //next i
        goto label_quit;
}
label_quad: {
        A = FaceArea(T, 0);
        i = 0;
        for (j = 1; j < 4; j++) {
            B = FaceArea(T, j);
            if (B > A) {
                A = B;
                i = j;
            }
        } //next i
        goto label_quit;
}
label_tetra: {
        goto label_quit;
}
label_pyram: {
        goto label_quit;
}
label_prism: {
        goto label_quit;
}
label_dhexa: {
        goto label_quit;
}
label_rhexa: {
        goto label_quit;
}

label_quit:
        return(A); };

        // ----------------------------------------------------------------------------------- //
        double Class_VolTri::Volume(int T) {

        // =================================================================================== //
        // double Class_VolTri::Volume(int T)                                                  //
        //                                                                                     //
        // Compute volume of a given simplex.                                                  //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - T        : int, simplex global index                                              //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - V        : double, simplex volume                                                 //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        int                dim = Vertex[0].size();
        double             V;

        // Counters


        // =================================================================================== //
        // COMPUTE SIMPLEX VOLUME                                                              //
        // =================================================================================== //
        switch(e_type[T]) {
            case 6  : goto label_tria; break;
            case 8  : goto label_quad; break;
            case 12 : goto label_tetra; break;
            case 16 : goto label_pyram; break;
            case 18 : goto label_prism; break;
            case 22 : goto label_dhexa; break;
            case 24 : goto label_rhexa; break;
        };

label_tria: {
//         if (dim == 2) {
//             V = 0.5*abs(Cross_Product(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                       Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]));
//         }
//         else if (dim == 3) {
//             V = 0.5*norm_2(Cross_Product(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                          Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]), dim);
//         }
        }
label_quad: {
//         if (dim == 2) {
//             V = 0.5*abs(Cross_Product(Vertex[Simplex[T][3]] - Vertex[Simplex[T][2]],
//                                       Vertex[Simplex[T][2]] - Vertex[Simplex[T][0]]))
//               + 0.5*abs(Cross_Product(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                       Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]));
//         }
//         else if (dim == 3) {
//             V = 0.5*norm_2(Cross_Product(Vertex[Simplex[T][3]] - Vertex[Simplex[T][2]],
//                                          Vertex[Simplex[T][2]] - Vertex[Simplex[T][0]]), 3)
//               + 0.5*norm_2(Cross_Product(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                          Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]), 3);
//         }
        }
label_tetra: {
//         if (dim == 2) {
//             V = 0.0;
//         }
//         else if (dim == 3) {
//             V = 1.0/6.0 * abs(Dot_Product(Vertex[Simplex[T][0]] - Vertex[Simplex[T][3]],
//                                           Cross_Product(Vertex[Simplex[T][1]] - Vertex[Simplex[T][3]],
//                                                         Vertex[Simplex[T][2]] - Vertex[Simplex[T][3]]), dim));
//         }
        }
label_pyram: {
        if (dim == 2) {
            V = 0.0; // to be done
        }
        else if (dim == 3) {
            V = 0.0; // to be done
        }
        }
label_prism: {
        if (dim == 2) {
            V = 0.0; // to be done
        }
        else if (dim == 3) {
            V = 0.0; // to be done
        }
        }
label_dhexa: {
        if (dim == 2) {
            V = 0.0; // to be done
        }
        else if (dim == 3) {
            V = 0.0; // to be done
        }
        }
label_rhexa: {
        if (dim == 2) {
            V = 0.0; // to be done
        }
        else if (dim == 3) {
            V = 0.0; // to be done
        }
        }

        return(V); };

                // ----------------------------------------------------------------------------------- //
        double Class_VolTri::Volume(dvector2D &X, int T) {

        // =================================================================================== //
        // double Class_VolTri::Volume(int T)                                                  //
        //                                                                                     //
        // Compute volume of a given simplex.                                                  //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - T        : int, simplex global index                                              //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - V        : double, simplex volume                                                 //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        int                dim = X[0].size();
        double             V;

        // Counters


        // =================================================================================== //
        // COMPUTE SIMPLEX VOLUME                                                              //
        // =================================================================================== //
        switch(e_type[T]) {
            case 6  : goto label_tria; break;
            case 8  : goto label_quad; break;
            case 12 : goto label_tetra; break;
            case 16 : goto label_pyram; break;
            case 18 : goto label_prism; break;
            case 22 : goto label_dhexa; break;
            case 24 : goto label_rhexa; break;
        };

label_tria: {
//         if (dim == 2) {
//             V = 0.5*abs(Cross_Product(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                       Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]));
//         }
//         else if (dim == 3) {
//             V = 0.5*norm_2(Cross_Product(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                          Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]), dim);
//         }
        }
label_quad: {
//         if (dim == 2) {
//             V = 0.5*abs(Cross_Product(Vertex[Simplex[T][3]] - Vertex[Simplex[T][2]],
//                                       Vertex[Simplex[T][2]] - Vertex[Simplex[T][0]]))
//               + 0.5*abs(Cross_Product(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                       Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]));
//         }
//         else if (dim == 3) {
//             V = 0.5*norm_2(Cross_Product(Vertex[Simplex[T][3]] - Vertex[Simplex[T][2]],
//                                          Vertex[Simplex[T][2]] - Vertex[Simplex[T][0]]), 3)
//               + 0.5*norm_2(Cross_Product(Vertex[Simplex[T][2]] - Vertex[Simplex[T][1]],
//                                          Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]]), 3);
//         }
        }
label_tetra: {
//         if (dim == 2) {
//             V = 0.0;
//         }
//         else if (dim == 3) {
//             V = 1.0/6.0 * abs(Dot_Product(Vertex[Simplex[T][0]] - Vertex[Simplex[T][3]],
//                                           Cross_Product(Vertex[Simplex[T][1]] - Vertex[Simplex[T][3]],
//                                                         Vertex[Simplex[T][2]] - Vertex[Simplex[T][3]]), dim));
//         }
        }
label_pyram: {
        if (dim == 2) {
            V = 0.0; // to be done
        }
        else if (dim == 3) {
            V = 0.0; // to be done
        }
        }
label_prism: {
        if (dim == 2) {
            V = 0.0; // to be done
        }
        else if (dim == 3) {
            V = 0.0; // to be done
        }
        }
label_dhexa: {
        if (dim == 2) {
            V = 0.0; // to be done
        }
        else if (dim == 3) {
            V = 0.0; // to be done
        }
        }
label_rhexa: {
        if (dim == 2) {
            V = 0.0; // to be done
        }
        else if (dim == 3) {
            V = 0.0; // to be done
        }
        }

        return(V); };
