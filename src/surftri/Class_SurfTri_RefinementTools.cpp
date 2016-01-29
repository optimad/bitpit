// ========================================================================== //
//                         - Class_SurfTri -                                  //
//                                                                            //
// Grid manager for unstructured meshes.                                      //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author   : Alessandro Alaia                                                //
// Version  : v3.0                                                            //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
# include "Class_SurfTri.hpp"

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void Class_SurfTri::invert_loc_num(
    int          I
) {

// ========================================================================== //
// void Class_SurfTri::invert_loc_num(                                        //
//     int    I)                                                              //
//                                                                            //
// Invert local numbering on simplex I                                        //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - I    : int, global index of simplex                                      //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool             flag_a, flag_n;
ivector1D        idummy1D;
ivector2D        idummy2D;

// Counters
int              i, m = Simplex[I].size();

// ========================================================================== //
// CHECK OPTIONAL DATA STRUCTURES                                             //
// ========================================================================== //
flag_n =  ((Normal.size() > 0) && (Normal.size() >= nSimplex));
flag_a =  ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex));

// ========================================================================== //
// INVERT LOCAL NUMBERING                                                     //
// ========================================================================== //

// Adjust vertex-simplex connectivity
idummy1D = Simplex[I];
for (i = 0; i < m; i++) {
    Simplex[I][i] = idummy1D[m-i-1];
} //next i

// Adjust simplex-simplex adjacency
if (flag_a) {
    idummy2D = Adjacency[I];
    for (i = 0; i < m; i++) {
        Adjacency[I][i] = idummy2D[m-i-1];
    } //next i
}

// Reverse normal
if (flag_n) {
    Normal[I] = -1.0*Normal[I];
}

return; }

// Simplex refinement ======================================================= //

// -------------------------------------------------------------------------- //
void Class_SurfTri::BinaryRefinement(
    double       h
) {

// ========================================================================== //
// void Class_SurfTri::BinaryRefinement(                                      //
//     double       h)                                                        //
//                                                                            //
// Refine tasselation using the mid-point rule. Each simplex is splitted      //
// into various simplicies by joining simplex vertices with its baricenter.   //
// The number of newly created simplicies depends on simplex type.            //
// Refinement is iteratively performed until simplex are above the given      //
// threshold.                                                                 //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - h      : double, simplex area                                            //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                            flag_n, flag_a;
int                             nV_new, nS_new;
int                             dim = Vertex[0].size();
int                             n, m;
double                          A;
pair<int, double>               dummy;
darray3E                        P, temp;

// Counters
int                             i;
int                             T;

P.fill(0.) ; temp.fill(0.) ;

// ========================================================================== //
// PARAMETERS                                                                 //
// ========================================================================== //
flag_a = ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex));
flag_n = ((Normal.size() > 0) && (Normal.size() >= nSimplex));

// ========================================================================== //
// RESIZE DATA STRUCTURE                                                      //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    // none

    // Number of simplicies which will be generated ------------------------- //
    nS_new = 0;
    nV_new = 0;
    for (T = 0; T < nSimplex; ++T) {
        m = Simplex[T].size();
        A = Area(T);
        if (A > h) {
            n = (int) (log10(A/h)/log10(m)) + 1;
            nS_new += pow(m, n) - 1;
            for (i = 0; i < n; ++i) {
                nV_new += pow(m, i);
            } //next i
        }
    } //next T

    // Resize data structure ------------------------------------------------ //
    Vertex.resize(nVertex + nV_new, temp);
    Simplex.resize(nSimplex + nS_new);
    if (flag_a) { Adjacency.resize(nSimplex + nS_new); };
    if (flag_n) { Normal.resize(nSimplex + nS_new, temp); }
    
}

// ========================================================================== //
// REFINEMENT LOOP                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    bitpit::LIFOStack< pair<int, double> >  stack(nSimplex + nS_new);

    // Initialize stack ----------------------------------------------------- //
    for (T = 0; T < nSimplex; T++) {
        dummy.first = T;
        dummy.second = Area(T);
        stack.push(dummy);
    } //next T

    // Iterative refinement ------------------------------------------------- //
    while (stack.TOPSTK > 0) {

        // Pop item from stack
        dummy = stack.pop();

        // Retrieve item info
        T = dummy.first;
        A = dummy.second;
        n = Simplex[T].size();

        if (A > h) {
            P = Baricenter(T);
            if (n == 0) {}
            else if (n == 1) {}
            else if (n == 2) {
                split_1segm2segm(T, P);
                if (0.5*A > h) {
                    dummy.first = T;
                    dummy.second = 0.5*A;
                    stack.push(dummy);
                    dummy.first = nSimplex-1;
                    dummy.second = 0.5*A;
                    stack.push(dummy);
                }
            }
            else if (n == 3) {
                //split_1tri3tri(T, P);
                if (A/3.0 > h) {
                }
            }
            else if (n == 4) {
                //split_1quad4quad(T, P);
                if (A/4.0 > h) {
                }
            }
            else {
                P = Baricenter(T);
                //split_1polyNpoly(T, P);                
                if (A/((double) n) > h) {
                }
            }
        }
    } //next item
}

// ========================================================================== //
// RESIZE DATA STRUCTURE                                                      //
// ========================================================================== //
ResizeVertex();
ResizeSimplex();
if (flag_a) { ResizeAdjacency(); }
if (flag_n) { ResizeNormal(); }

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::BinaryRefinement(
    dvecarr3E   &V,
    double       h
) {

// ========================================================================== //
// void Class_SurfTri::BinaryRefinement(                                      //
//     dvecarr3E   &V,                                                        //
//     double       h)                                                        //
//                                                                            //
// Refine tasselation using the mid-point rule. Each simplex is splitted      //
// into various simplicies by joining simplex vertices with its baricenter.   //
// The number of newly created simplicies depends on simplex type.            //
// Refinement is iteratively performed until simplex are above the given      //
// threshold. Vertex coordinate list is provided externally.                  //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - V      : dvecarr3E, vertex coordinate list. V[i][0], V[i][1] are the     //
//            x, y, coordinates of the i-th vertex                            //
// - h      : double, simplex area                                            //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                            flag_n, flag_a;
int                             nV = V.size();
int                             nV_new, nS_new;
int                             dim = V[0].size();
int                             n, m;
double                          A;
pair<int, double>               dummy;
darray3E                        P, temp;

// Counters
int                             i;
int                             T;

P.fill(0.) ; temp.fill(0.) ;

// ========================================================================== //
// PARAMETERS                                                                 //
// ========================================================================== //
flag_a = ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex));
flag_n = ((Normal.size() > 0) && (Normal.size() >= nSimplex));

// ========================================================================== //
// RESIZE DATA STRUCTURE                                                      //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    // none

    // Number of simplicies which will be generated ------------------------- //
    nS_new = 0;
    nV_new = 0;
    for (T = 0; T < nSimplex; ++T) {
        m = Simplex[T].size();
        A = Area(T, V);
        if (A > h) {
            n = (int) (log2(A) - log2(h)) - 1;
            nS_new += pow(m, n);
            for (i = 0; i < n; ++i) {
                nV_new += pow(m, i);
            } //next i
        }
    } //next T

    // Resize data structure ------------------------------------------------ //
    cout << "resizing" << endl;
    V.resize(nV + nV_new, temp);
    Simplex.resize(nSimplex + nS_new);
    if (flag_a) { Adjacency.resize(nSimplex + nS_new); };
    if (flag_n) { Normal.resize(nSimplex + nS_new, temp); }
    
}

// ========================================================================== //
// REFINEMENT LOOP                                                            //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    bitpit::LIFOStack< pair<int, double> >  stack(nSimplex + nS_new);

    // Initialize stack ----------------------------------------------------- //
    for (T = 0; T < nSimplex; T++) {
        dummy.first = T;
        dummy.second = Area(T, V);
        stack.push(dummy);
    } //next T

    // Iterative refinement ------------------------------------------------- //
    while (stack.TOPSTK > 0) {

        // Pop item from stack
        dummy = stack.pop();

        // Retrieve item info
        T = dummy.first;
        A = dummy.second;
        n = Simplex[T].size();

        if (A > h) {
            P = Baricenter(T, V);
            if (n == 0) {}
            else if (n == 1) {}
            else if (n == 2) {
                V.push_back(P);
                split_1segm2segm(T, (int) V.size()-1, V);
                if (0.5*A > h) {
                    dummy.first = T;
                    dummy.second = 0.5*A;
                    stack.push(dummy);
                    dummy.first = nSimplex-1;
                    dummy.second = 0.5*A;
                    stack.push(dummy);
                }
            }
            else if (n == 3) {
                //split_1tri3tri(T, P);
                if (A/3.0 > h) {
                }
            }
            else if (n == 4) {
                //split_1quad4quad(T, P);
                if (A/4.0 > h) {
                }
            }
            else {
                //split_1polyNpoly(T, P);                
                if (A/((double) n) > h) {
                }
            }
        }
    } //next item
}

return; };


// ----------------------------------------------------------------------------------- //
void Class_SurfTri::Collapse_2Simplex(int T, int rule) {

// =================================================================================== //
// void Class_SurfTri::Collapse_2Simplex(int T, int rule)                              //
//                                                                                     //
// Collapse a 2-simplex using the given rule:                                          //
//    rule = 0  --> vertex 1 (local index) is collapsed onto vertex 0                  //
//    rule = 1  --> vertex 0 (local index) is collapsed onto vertex 1                  //
//    rule = -1 --> vertex 0 and 1 (local index) are collapsed at midpoint             //
// =================================================================================== //
// INPUT                                                                               //
// =================================================================================== //
// - T       : int, simplex global index                                               //
// - rule    : int, rule used to collapse simplex.                                     //
// =================================================================================== //
// OUTPUT                                                                              //
// =================================================================================== //
// - none                                                                              //
// =================================================================================== //

// =================================================================================== //
// VARIABLES DECLARATION                                                               //
// =================================================================================== //

// Local variables
int                    A, B, V;
bvector1D              flag(2, false);
ivector1D              adj;
vector<int>::iterator  it;

// Counters
int                i, j, k;

darray3E               P, temp;
P.fill(0.) ; temp.fill(0.) ;

// =================================================================================== //
// COLLAPSE SIMPLEX                                                                    //
// =================================================================================== //

// Vertex to be retained ------------------------------------------------------------- //
switch (rule) {
    case -1:
        flag[0] = flag[1] = true;
        P = Baricenter(T);
        AddVertex(P);
        V = nVertex-1;
    break;
    case 0:
        flag[1] = true;
        V = Simplex[T][0];
    break;
    case 1:
        flag[0] = true;
        V = Simplex[T][1];
    break;
}

// Update simplex-vertex connectivity ------------------------------------------------ //
for (i = 0; i < 2; i++) {
    if (flag[i]) {
        for (j = 0; j < Adjacency[T][i].size(); j++) {
            A = Adjacency[T][i][j];
            if (A >= 0) {
                k = vertex(A, Simplex[T][i]);
                Simplex[A][k] = V;
            }
        } //next j
    }
} //next i

// Update adjacencies ---------------------------------------------------------------- //

// Update neighbors
for (i = 0; i < 2; i++) {
    adj.resize(0);
    k = (i + 1) % Simplex[T].size();
    if (Adjacency[T][k][0] >= 0) {
        adj.resize(Adjacency[T][k].size());
        adj = Adjacency[T][k];
    }
    for (j = 0; j < Adjacency[T][i].size(); j++) {
        A = Adjacency[T][i][j];
        if (A >= 0) {
            k = edge(A, T);
            it = find(Adjacency[A][k].begin(), Adjacency[A][k].end(), T);
            Adjacency[A][k].erase(it);
            Adjacency[A][k].insert(Adjacency[A][k].end(), adj.begin(), adj.end());
        }
    } //next j
} //next i

// Update simplex T
Adjacency[T].resize(2, ivector1D(1, -1));
Adjacency[T][0][0] = -1;
Adjacency[T][1][0] = -1;


return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::split_1segm2segm(
    int          T,
    darray3E    &P
) {

// ========================================================================== //
// void Class_SurfTri::split_1segm2segm(                                      //
//     int          T,                                                        //
//     dvector1D   &P)                                                        //
//                                                                            //
// Split a given 2-simplex in two 2-simplicies by inserting point P.          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex global index                                       //
// - P      : dvector1D, point coordinates                                    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                flag_n, flag_a;
int                 s_type ;
ivector1D           idummy1D(2, -1);
darray3E            ddummy1D;
ivector2D           idummy2D(2, ivector1D(1, -1));

// Counters
int                 S = nSimplex;

// ========================================================================== //
// CHECK SIMPLEX TYPE                                                         //
// ========================================================================== //

// Simplex type
s_type = Simplex[T].size();
if (s_type != 2) { return; }

// ========================================================================== //
// SET PARAMETERS                                                             //
// ========================================================================== //
flag_n = ((Normal.size() > 0) && (Normal.size() >= nSimplex));
flag_a = ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex));

// ========================================================================== //
// SPLIT SIMPLEX                                                              //
// ========================================================================== //

// Modify vertex list ------------------------------------------------------- //
AddVertex(P);

// Modify simplex-vertex connectivity --------------------------------------- //
idummy1D[0] = nVertex-1;
idummy1D[1] = Simplex[T][1];
AddSimplex(idummy1D);
S = nSimplex-1;
Simplex[T][1] = nVertex-1;

// Modify simplex-simplex adjacencies --------------------------------------- //
if (flag_a) {
    idummy2D[0][0] = T;
    idummy2D[1] = Adjacency[T][1];
    SetAdjacency(S, idummy2D);
    Adjacency[T][1].resize(1);
    Adjacency[T][1][0] = S;
}

// Generate normals for new simplex ----------------------------------------- //
if (flag_n) {
    ddummy1D = Vertex[Simplex[S][1]] - Vertex[Simplex[S][0]];
    ddummy1D = ddummy1D/norm2(ddummy1D);
    SetNormal(S, ddummy1D);
    Normal[T] = Vertex[Simplex[T][1]] - Vertex[Simplex[T][0]];
    Normal[T] = Normal[T]/norm2(Normal[T]);
}

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::split_1segm2segm(
    int          T,
    int          V,
    dvecarr3E   &X
) {

// ========================================================================== //
// void Class_SurfTri::split_1segm2segm(                                      //
//     int          T,                                                        //
//     int          V,                                                        //
//     dvecarr3E   &X)                                                        //
//                                                                            //
// Split a given 2-simplex in two 2-simplicies by inserting point P.          //
// Vertex coordinate list is provided externally.                             //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T      : int, simplex global index                                       //
// - V      : int, global index of vertex used to split segment               //
// - X      : dvecarr3E, external vertex coorindate list. X[i][0], X[i][1],   //
//            ... are the x, y, ... coordinates of the i-th vertex            //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                flag_n, flag_a;
int                 s_type;
ivector1D           idummy1D(2, -1);
darray3E            ddummy1D;
ivector2D           idummy2D(2, ivector1D(1, -1));

// Counters
int                 S = nSimplex;

// ========================================================================== //
// CHECK SIMPLEX TYPE                                                         //
// ========================================================================== //

// Simplex type
s_type = Simplex[T].size();
if (s_type != 2) { return; }

// ========================================================================== //
// SET PARAMETERS                                                             //
// ========================================================================== //
flag_n = ((Normal.size() > 0) && (Normal.size() >= nSimplex));
flag_a = ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex));

// ========================================================================== //
// SPLIT SIMPLEX                                                              //
// ========================================================================== //

// Modify simplex-vertex connectivity --------------------------------------- //
idummy1D[0] = V;
idummy1D[1] = Simplex[T][1];
AddSimplex(idummy1D);
S = nSimplex-1;
Simplex[T][1] = V;

// Modify simplex-simplex adjacencies --------------------------------------- //
if (flag_a) {
    idummy2D[0][0] = T;
    idummy2D[1] = Adjacency[T][1];
    SetAdjacency(S, idummy2D);
    Adjacency[T][1].resize(1);
    Adjacency[T][1][0] = S;
}

// Generate normals for new simplex ----------------------------------------- //
if (flag_n) {
    ddummy1D = X[Simplex[S][1]] - X[Simplex[S][0]];
    ddummy1D = ddummy1D/norm2(ddummy1D);
    SetNormal(S, ddummy1D);
    Normal[T] = X[Simplex[T][1]] - X[Simplex[T][0]];
    Normal[T] = Normal[T]/norm2(Normal[T]);
}

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::Split_2Simplex(
    int          T
) {

// ========================================================================== //
// void Class_SurfTri::Split_2Simplex(                                        //
//     int T)                                                                 //
//                                                                            //
// Split a given 2-simplex at mid-point.                                      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - T        : int, simplex global index.                                    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                   flag_a, flag_n;
int                    V;
vector<int>::iterator  it;
ivector1D              idummy1D(2, -1);
darray3E               P;

// Counters
int                    A;
int                    i, j;

// ========================================================================== //
// INITIALIZE VARIABLES                                                       //
// ========================================================================== //
flag_a = ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex));
flag_n = ((Simplex.size() > 0) && (Simplex.size() >= nSimplex));

// ========================================================================== //
// SPLIT SIMPLEX                                                              //
// ========================================================================== //

// Add new vertex to the tasselation ---------------------------------------- //
P = Baricenter(T);
AddVertex(P);
V = nVertex - 1;

// Update simplex-vertex connectivity --------------------------------------- //

// Add new simplex to the tasselation
idummy1D[0] = V;
idummy1D[1] = Simplex[T][1];
AddSimplex(idummy1D);

// Update simplex T
Simplex[T][1] = V;

// Update adjacency --------------------------------------------------------- //
if (flag_a) {

    // Update adjacency for the newly created simplex
    Adjacency.push_back(ivector2D(2, ivector1D(1, -1)));
    Adjacency[nSimplex-1][0][0] = T;
    Adjacency[nSimplex-1][1] = Adjacency[T][1];
    
    // Update adjacency for simplex T
    Adjacency[T][1].resize(1);
    Adjacency[T][1][0] = nSimplex-1;
    
    // Update adjacency for neighboring simplicies
    if (Adjacency[nSimplex-1][1][0] >= 0) {
        for (i = 0; i < Adjacency[nSimplex-1][1].size(); i++) {
            A = Adjacency[nSimplex-1][1][i];
            j = vertex(A, Simplex[nSimplex-1][1]);
            it = find(Adjacency[A][j].begin(), Adjacency[A][j].end(), T);
            *it = nSimplex-1;
        } //next i
    }
}

// Update normals ----------------------------------------------------------- //
if (flag_n) {
    ResizeNormal();
    Normal[nSimplex-1] = Normal[nSimplex-2];
}

return; };

// -------------------------------------------------------------------------- //
void Class_SurfTri::Split_2Simplex(
    dvecarr3E   &X,
    int          T
) {

// ========================================================================== //
// void Class_SurfTri::Split_2Simplex(                                        //
//     dvecarr3E   &X,                                                        //
//     int          T)                                                        //
//                                                                            //
// Split a given 2-simplex at mid-point. Vertex list is provided externally.  //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - X        : dvecarr3E, vertex coordinate list. X[i][0], X[i][1], ...      //
//              are the x, y, ... coordinates of the i-th vertex              //
// - T        : int, simplex global index.                                    //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                   flag_a, flag_n;
int                    V, nV = X.size();
vector<int>::iterator  it;
ivector1D              idummy1D(2, -1);
darray3E               P;

// Counters
int                    A;
int                    i, j;

// ========================================================================== //
// INITIALIZE VARIABLES                                                       //
// ========================================================================== //
flag_a = ((Adjacency.size() > 0) && (Adjacency.size() >= nSimplex));
flag_n = ((Normal.size() > 0) && (Normal.size() >= nSimplex));

// ========================================================================== //
// SPLIT SIMPLEX                                                              //
// ========================================================================== //

// Add new vertex to the tasselation ---------------------------------------- //
P = Baricenter(T, X);
X.push_back(P);
V = nV;

// Update simplex-vertex connectivity --------------------------------------- //

// Add new simplex to the tasselation
idummy1D[0] = V;
idummy1D[1] = Simplex[T][1];
AddSimplex(idummy1D);

// Update simplex T
Simplex[T][1] = V;

// Update adjacency --------------------------------------------------------- //
if (flag_a) {

    // Update adjacency for the newly created simplex
    Adjacency.push_back(ivector2D(2, ivector1D(1, -1)));
    Adjacency[nSimplex-1][0][0] = T;
    Adjacency[nSimplex-1][1] = Adjacency[T][1];
    
    // Update adjacency for simplex T
    Adjacency[T][1].resize(1);
    Adjacency[T][1][0] = nSimplex-1;
    
    // Update adjacency for neighboring simplicies
    if (Adjacency[nSimplex-1][1][0] >= 0) {
        for (i = 0; i < Adjacency[nSimplex-1][1].size(); i++) {
            A = Adjacency[nSimplex-1][1][i];
            j = vertex(A, Simplex[nSimplex-1][1]);
            it = find(Adjacency[A][j].begin(), Adjacency[A][j].end(), T);
            *it = nSimplex-1;
        } //next i
    }
}

// Update normals ----------------------------------------------------------- //
if (flag_n) {
    ResizeNormal();
    Normal[nSimplex-1] = Normal[nSimplex-2];
}

return; };
    
// ----------------------------------------------------------------------------------- //
void Class_SurfTri::SplitEdge(int T, int i) {

// =================================================================================== //
// void Class_SurfTri::SplitEdge(int T, int i)                                         //
//                                                                                     //
// Split edge at mid-point. (2-SIMPLICIES ONLY)                                        //
// =================================================================================== //
// INPUT                                                                               //
// =================================================================================== //
// - T     : int, simplex global index                                                 //
// - i     : int, edge local index.                                                    //
// =================================================================================== //
// OUTPUT                                                                              //
// =================================================================================== //
// - none                                                                              //
// =================================================================================== //

// =================================================================================== //
// VARIABLES DECLARATION                                                               //
// =================================================================================== //

// Local variables
int              V1, V2, V3;
ivector1D        iummy(3, -1);
darray3E         P;


// Counters
int              A, B, C, S;
int              j, k;
int              ii;

// =================================================================================== //
// UPDATE VERTEX LIST                                                                  //
// =================================================================================== //

// Simplex vertexes ------------------------------------------------------------------ //

// Vertex local numbering
j = (i + 1) % 3;
k = (j + 1) % 3;

// Vertex global index
V1 = Simplex[T][i];
V2 = Simplex[T][j];
V3 = Simplex[T][k];

// Add the new vertex ---------------------------------------------------------------- //
P = 0.5 * (Vertex[Simplex[T][i]] + Vertex[Simplex[T][j]]);
AddVertex(P);

// =================================================================================== //
// UPDATE SIMPLEX - ADJACENCY LIST                                                     //
// =================================================================================== //

// Update simplex list --------------------------------------------------------------- //

// Overwrite simplex T
Simplex[T][0] = V1;
Simplex[T][1] = nVertex-1;
Simplex[T][2] = V3;

// Add new simplex
iummy[0] = nVertex-1;
iummy[1] = V2;
iummy[2] = V3;
AddSimplex(iummy);

// Adjust adjacencies ---------------------------------------------------------------- //

// Old adjacencies
A = Adjacency[T][i][0];
B = Adjacency[T][j][0];
C = Adjacency[T][k][0];

// Overwrite adjacencies for simplex T
Adjacency[T][0][0] = A;
Adjacency[T][1][0] = nSimplex-1;
Adjacency[T][2][0] = C;

// Adjacencies for the new simplex
Adjacency.push_back(ivector2D(3, ivector1D(1,-1)));
if (A >= 0) {
    Adjacency[nSimplex-1][0][0] = nSimplex;
}
Adjacency[nSimplex-1][1][0] = B;
Adjacency[nSimplex-1][2][0] = T;

// Adjacencies for the new simplex - neighbors
if (B >= 0) {
    j = edge(B, T);
    Adjacency[B][j][0] = nSimplex - 1;
}

// =================================================================================== //
// NEIGHBOR SIMPLEX                                                                    //
// =================================================================================== //
if (A >= 0) {

    // Update simplex list ----------------------------------------------------------- //

    // Vertex local numbering
    i = edge(A, T);
    j = (i + 1) % 3;
    k = (j + 1) % 3;

    // Vertex global index
    S = A;
    V1 = Simplex[S][i];
    V2 = Simplex[S][j];
    V3 = Simplex[S][k];

    // Overwrite simplex S
    Simplex[S][0] = nVertex - 1;
    Simplex[S][1] = V2;
    Simplex[S][2] = V3;

    // Add new simplex
    iummy[0] = nVertex - 1;
    iummy[1] = V3;
    iummy[2] = V1;
    AddSimplex(iummy);

    // Adjust adjacencies ------------------------------------------------------------ //

    // Old adjacencies
    A = Adjacency[S][i][0];
    B = Adjacency[S][j][0];
    C = Adjacency[S][k][0];

    // Overwrite adjacencies for simplex S
    Adjacency[S][0][0] = A;
    Adjacency[S][1][0] = B;
    Adjacency[S][2][0] = nSimplex - 1;

    // Adjacencies for the new simplex
    Adjacency.push_back(ivector2D(3, ivector1D(1, -1)));
    Adjacency[nSimplex-1][0][0] = S;
    Adjacency[nSimplex-1][1][0] = C;
    Adjacency[nSimplex-1][2][0] = nSimplex - 2;

    // Adjacencies for the new simplex - neighbors
    if (C >= 0) {
        j = edge(C, S);
        Adjacency[C][j][0] = nSimplex - 1;
    }
}

return; };

// ----------------------------------------------------------------------------------- //
void Class_SurfTri::CollapseEdge(int T, int i, int m) {

// =================================================================================== //
// void Class_SurfTri::CollapseEdge(int T, int i, int m)                               //
//                                                                                     //
// Collapse a given edge in a simplex (2-simplicies only)                              //
// =================================================================================== //
// INPUT                                                                               //
// =================================================================================== //
// - T     : int, simplex global index                                                 //
// - i     : int, edge local index                                                     //
// - m     : int, vertex to be retained (local index)                                  //
// =================================================================================== //
// OUTPUT                                                                              //
// =================================================================================== //
// - none                                                                              //
// =================================================================================== //

// =================================================================================== //
// VARIABLES DECLARATION                                                               //
// =================================================================================== //

// Local variables
bool             flag, ring_flag = true;
int              retained, lost;
darray3E         P;
ivector1D        ring1;

// Counters
int              A, B, C;
int              j, k;
int              ii, jj;

// =================================================================================== //
// COLLAPSE EDGE AT MIDPOINT                                                           //
// =================================================================================== //

// Vertex local index ---------------------------------------------------------------- //
j = (i + 1) % 3;
k = (j + 1) % 3;

// Set coordinates for target point -------------------------------------------------- //
if (m == -1) {

    // Edge mid_point coordinates
    P = Edge_midPoint(T, i);

    // Overwrite vertex coordinates
    Vertex[Simplex[T][i]] = P;

    // Retained point
    retained = i;
    lost = j;
}
else if (m == i) {

    // Retained
    retained = i;
    lost = j;

}
else {

    // Retained point
    retained = j;
    lost = i;

}

// =================================================================================== //
// UPDATE SIMPLEX LIST                                                                 //
// =================================================================================== //

// 1-ring around lost vertex
ring1 = Ring_1(T, lost, ring_flag);

// Update simplicies in the 1-ring of lost vertex
for (jj = 0; jj < ring1.size(); jj++) {
    ii = 0;
    flag = true;
    while (flag && (ii < 3)) {
        if (Simplex[ring1[jj]][ii] == Simplex[T][lost]) {
            Simplex[ring1[jj]][ii] = Simplex[T][retained];
            flag = false;
        }
        ii++;
    } //next ii
} //next jj

// =================================================================================== //
// ADJUST ADJACENCY LIST FOR SIMPLEX T                                                 //
// =================================================================================== //

// Old adjacencies
A = Adjacency[T][i][0];
B = Adjacency[T][j][0];
C = Adjacency[T][k][0];

// Update adjacencies for simplex T
Adjacency[T][0][0] = -1;
Adjacency[T][1][0] = -1;
Adjacency[T][2][0] = -1;

// Update adjacencies for neighboring elements
if (C >= 0) {
    j = edge(C, T);
    Adjacency[C][j][0] = B;
}
if (B >= 0) {
    j = edge(B, T);
    Adjacency[B][j][0] = C;
}

// =================================================================================== //
// ADJUST ADJACENCY LIST FOR T-NEIGHBOR                                                //
// =================================================================================== //
if (A >= 0) {

    // Vertex local index
    i = edge(A, T);
    j = (i + 1) % 3;
    k = (j + 1) % 3;

    // Old adjacencies
    T = A;
    A = Adjacency[T][i][0];
    B = Adjacency[T][j][0];
    C = Adjacency[T][k][0];

    // Update adjacencies for simplex T
    Adjacency[T][0][0] = -1;
    Adjacency[T][1][0] = -1;
    Adjacency[T][2][0] = -1;

    // Update adjacencies for neighboring elements
    if (C >= 0) {
        j = edge(C, T);
        Adjacency[C][j][0] = B;
    }
    if (B >= 0) {
        j = edge(B, T);
        Adjacency[B][j][0] = C;
    }

}

return; };

// Voronoi diagrams ================================================================== //

// -----------------(Ricontrollare)--------------------------------------------------- //
void Class_SurfTri::Voronoi(Class_SurfTri &Voronoi) {

// =================================================================================== //
// void Class_SurfTri::Voronoi(Class_SurfTri &Voronoi)                                 //
//                                                                                     //
// Compute approximated Voronoi tasselation from the surface triangulation.            //
// =================================================================================== //
// INPUT                                                                               //
// =================================================================================== //
// - Voronoi   : Class_SurfTri, with Voronoi diagram                                   //
// =================================================================================== //
// OUTPUT                                                                              //
// =================================================================================== //
// - none                                                                              //
// =================================================================================== //

// =================================================================================== //
// VARIABLES DECLARATION                                                               //
// =================================================================================== //

// Local variables
bool             inside, s1, s2, s3, ring_flag = true;
darray3E         P;
darray3E         xP, V1, V2, V3;
ivector1D        S2CC(nSimplex, -1);
ivector2D        S2MP(nSimplex, ivector1D(3, -1));

// Counters
int              I, J;
int              i, j, k, e, d;

// =================================================================================== //
// COMPUTE EDGES MIDPOINTS                                                             //
// =================================================================================== //
for (i = 0; i < nSimplex; i++) {
    for (j = 0; j < 3; j++) {
        if (S2MP[i][j] == -1) {

            // Compute edge midpoint
            P = Edge_midPoint(i, j);

            // Add vertex to the Voronoi vertex list
            Voronoi.AddVertex(P);

            // Update {Simplex->MidPoint} map
            S2MP[i][j] = Voronoi.nVertex-1;
            I = Adjacency[i][j][0];
            J = edge(I, i);
            S2MP[I][J] = Voronoi.nVertex-1;

        }
    } //next j
} //next i

// =================================================================================== //
// COMPUTE TRIANGLES CIRCUMCENTERS                                                     //
// =================================================================================== //
for (i = 0; i < nSimplex; i++) {

    // Tringle vertexes
    V1 = Vertex[Simplex[i][0]];
    V2 = Vertex[Simplex[i][1]];
    V3 = Vertex[Simplex[i][2]];

    // Compute triangle circumcenter
    xP = CircumCenter(i);

    // Add circumcenter to Voronoi vertex list and update {Simplex->CircumCenters} map
    s1 = PointsOnSameSide(xP,V3,V1,V2);
    s2 = PointsOnSameSide(xP,V1,V2,V3);
    s3 = PointsOnSameSide(xP,V2,V3,V1);
    inside = (s1 && s2 && s3);
    if (inside) {
        Voronoi.AddVertex(P);
        S2CC[i] = Voronoi.nVertex - 1;
    }
    else {
        if (!s1) {
            S2CC[i] = S2MP[i][0];
        }
        if (!s2) {
            S2CC[i] = S2MP[i][1];
        }
        if (!s3) {
            S2CC[i] = S2MP[i][2];
        }
    }

} //next i

// =================================================================================== //
// COMPUTE VORONOI CELLS                                                               //
// =================================================================================== //

// Resize data structure
Voronoi.Simplex.resize(nVertex);

// Create Voronoi cells
for (i = 0; i < nSimplex; i++) {
    for (j = 0; j < 3; j++) {
        if (Voronoi.Simplex[Simplex[i][j]].size() == 0) {

            // Scope variables
            bool           flag;
            ivector1D      Ring1, iummy;

            // Select triangles in the 1-ring of vertex
            Ring1 = Ring_1(i, j, ring_flag);

            // Add vertex to the Voronoi cell
            for (k = 0; k < Ring1.size(); k++) {

                // Add edge midpoint
                e = 0;
                d = -1;
                flag = true;
                while(flag && (e < 3)) {
                    if (Simplex[Ring1[k]][e] == Simplex[i][j]) { 
                        d = e;
                        flag = false;
                    }
                    e++;
                } //next e
                iummy.push_back(S2MP[Ring1[k]][d]);

                // Add simplex circumcenter
                if (S2CC[Ring1[k]] >= 0) {
                    iummy.push_back(S2CC[Ring1[k]]);
                }

            } //next k

            // Update Voronoi cell list
            Voronoi.Simplex[Simplex[i][j]] = iummy;
            Voronoi.nSimplex++;
        }
    } //next j
} //next i


return; };

// -----------------(implementare)---------------------------------------------------- //
void Class_SurfTri::Voronoi(Class_SurfTri &Voronoi, dvector2D &V) {

return; };
