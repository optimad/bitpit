// ========================================================================== //
//                 - COMPUTATIONAL GEOMETRY PACKAGE -                         //
//                                                                            //
// Routines for computational geometry.                                       //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author      :   Alessandro Alaia                                           //
// Version     :   v1.0                                                       //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
# include "Delaunay.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace CGPolygon2D;

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //

// -------------------------------------------------------------------------- //
void QDelaunay::QDelaunay2D::PlaceVertices(
    Class_VolTri                    &Tri,
    Class_CG_Polygon2D              &Dpoly,
    vector<Class_CG_Polygon2D>      &Hpoly,
    ivector2D                       *domain,
    vector<ivector2D*>              &holes,
    vector<ivector2D*>              &segments
) {

// ========================================================================== //
// void QDelaunay::QDelaunay2D::PlaceVertices(                             //
//     Class_VolTri                    &Tri,                                  //
//     Class_CG_Polygon2D              &Dpoly,                                //
//     vector<Class_CG_Polygon2D>      &Hpoly,                                //
//     ivector2D                       *domain,                               //
//     vector<ivector2D*>              &holes,                                //
//     vector<ivector2D*>              &segments)                             //
//                                                                            //
// Plce vertices inside domain for quad mesh generation.                      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - Tri        : Class_VolTri, unstructured volume mesh manager, containing  //
//                vertex of surface mesh.                                     //
// - Dpoly      : Class_CG_Polygon, polygonal shape describing domain         //
//                boundaries                                                  //
// - Hpoly      : vector<Class_CG_Polygon>, polygonal shapes describing each  //
//                domain hole.                                                //
// - domain     : ivector2D*, pointer to vertex-simplex connectivity for      //
//                domain boundaries.                                          //
// - holes      : vector<vector2D*>, vector of pointers to vertex-simplex     //
//                connectivity for each domain hole.                          //
// - segments   : vector<ivector2D*>, vector of pointers to vecrtex-simplex   //
//                connectivity for each domain segment                        //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Parameters
double const                            pi = 3.14159265358979;

// Local variables
int                                     nbins = 1;
int                                     nH = holes.size(), nS = segments.size();
double                                  dx, dy;
array<double, 2>                        xlim, ylim;
ivector2D                               ring1;
a3vector2D                              B_normal;
a3vector3D                              H_normal(nH);
a3vector3D                              S_normal(nS);
queue< pair< int, array<double, 2> > >  Vqueue;                       
vector<bitpit::KdTree<2, a3vector1D> >  Ktree(nbins*nbins, bitpit::KdTree<2, a3vector1D>(nbins));

// Counters
int                                     H, S;


Tri.Vertex.resize(10*Tri.nVertex*Tri.nVertex);

// ========================================================================== //
// CHECK INPUT                                                                //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    // none

    // Check for empty surface mesh ----------------------------------------- //
    if (Tri.nVertex == 0) { return; }
}

// ========================================================================== //
// COMPUTE SIMPLEX NORMALS                                                    //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int                 n;
    int                 T, A, B;
    array<double, 2>    v;

    // Compute simplex normals for each edge on domain boundaries ----------- //
    n = (*domain).size();
    B_normal.resize(n);
    for (T = 0; T < n; ++T) {
        A = (*domain)[T][0];
        B = (*domain)[T][1];
        B_normal[T][0] = Tri.Vertex[B][0] - Tri.Vertex[A][0];
        B_normal[T][1] = Tri.Vertex[B][1] - Tri.Vertex[A][1];
        B_normal[T] = B_normal[T]/norm2(B_normal[T]);
        swap(B_normal[T][0], B_normal[T][1]);
        B_normal[T][1] = -B_normal[T][1];
    } //next T

    // Compute simplex normals for each edge on domain holes ---------------- //
    for (H = 0; H < nH; ++H) {
        n = (*(holes[H])).size();
        H_normal[H].resize(n);
        for (T = 0; T < n; ++T) {
            A = (*(holes[H]))[T][0];
            B = (*(holes[H]))[T][1];
            H_normal[H][T][0] = Tri.Vertex[B][0] - Tri.Vertex[A][0];
            H_normal[H][T][1] = Tri.Vertex[B][1] - Tri.Vertex[A][1];
            H_normal[H][T] = H_normal[H][T]/norm2(H_normal[H][T]);
            swap(H_normal[H][T][0], H_normal[H][T][1]);
            H_normal[H][T][1] = -H_normal[H][T][1];
        } //next T
    } //next H

    // Compute simplex normals for each edge on domain segments ------------- //
    for (S = 0; S < nS; ++S) {
        n = (*(segments[S])).size();
        S_normal[S].resize(n);
        for (T = 0; T < n; ++T) {
            A = (*(segments[S]))[T][0];
            B = (*(segments[S]))[T][1];
            S_normal[S][T][0] = Tri.Vertex[B][0] - Tri.Vertex[A][0];
            S_normal[S][T][1] = Tri.Vertex[B][1] - Tri.Vertex[A][1];
            S_normal[S][T] = S_normal[S][T]/norm2(S_normal[S][T]);
            swap(S_normal[S][T][0], S_normal[S][T][1]);
            S_normal[S][T][1] = -S_normal[S][T][1];
        } //next T
    } //next H
}

// ========================================================================== //
// COMPUTE 1-RING OF VERTICES ON SURFACE MESH                                 //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int                     i;
    int                     T, V;
    int                     n;

    // Initialize variables ------------------------------------------------- //
    ring1.resize(Tri.nVertex);

    // Build vertex-simplex connectivity ------------------------------------ //

    // Domain boundaries
    n = (*domain).size();
    for (T = 0; T < n; ++T) {
        for (i = 0; i < 2; ++i) {
            V = (*domain)[T][i];
            ring1[V].push_back(T);
        } //next i
    } //next T

    // Domain holes
    for (H = 0; H < nH; ++H) {
        n = (*(holes[H])).size();
        for (T = 0; T < n; ++T) {
            for (i = 0; i < 2; ++i) {
                V = (*(holes[H]))[T][i];
                ring1[V].push_back(T);
            } //next i
        } //next T
    } //next H

    // Domain segments
    for (S = 0; S < nS; ++S) {
        n = (*(segments[S])).size();
        for (T = 0; T < n; ++T) {
            for (i = 0; i < 2; ++i) {
                V = (*(segments[S]))[T][i];
                ring1[V].push_back(T);
            } //next i
        } //next T
    } //next H
}

// ========================================================================== //
// INIITALIZE BOUNDING BOX                                                    //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    int                         T, V;
    int                         n;
    int                         i;

    // Compute vertex bounding box ------------------------------------------ //

    // Initialize bounding box limits
    V = (*domain)[0][0];
    xlim[0] = xlim[1] = Tri.Vertex[V][0];
    ylim[0] = ylim[1] = Tri.Vertex[V][1];

    // Update bounding box limits for domain boundaries
    n = (*domain).size();
    for (T = 0; T < n; ++T) {
        for (i = 0; i < 2; ++i) {
            V = (*domain)[T][i];
            xlim[0] = min(xlim[0], Tri.Vertex[V][0]);
            xlim[1] = max(xlim[1], Tri.Vertex[V][0]);
            ylim[0] = min(ylim[0], Tri.Vertex[V][1]);
            ylim[1] = max(ylim[1], Tri.Vertex[V][1]);
        } //next i
    } //next T

    // Compute bins width --------------------------------------------------- //
    dx = (xlim[1] - xlim[0])/((double) nbins);
    dy = (ylim[1] - ylim[0])/((double) nbins);
    xlim[0] -= 0.05*dx;
    xlim[1] += 0.05*dx;
    ylim[0] -= 0.05*dy;
    ylim[1] += 0.05*dy;
    dx = (xlim[1] - xlim[0])/((double) nbins);
    dy = (ylim[1] - ylim[0])/((double) nbins);    

}

// ========================================================================== //
// INIITALIZE INSERTION PROCEDURE                                             //
// ========================================================================== //
{
    // Scope variables ------------------------------------------------------ //
    int                                                 n;
    int                                                 i, j, ix, iy, I;
    int                                                 V, T, A, B;
    array<double, 2>                                    v, P;
    a3vector1D*                                         X_;
    pair<int, array<double, 2> >                        dummy;
    dvector1D                                           xlim(2, 0.0), ylim(2, 0.0);

    // Sort vertices on bins (domain boundaries) ---------------------------- //
    n = (*domain).size();
    for (T = 0; T < n; ++T) {

        // Vertex global index
        V = (*domain)[T][0];

        // Simplex incident on vertex V
        A = ring1[V][0];
        B = ring1[V][1];

        // Compute vertex normal
        v[0] = - (B_normal[A][0] + B_normal[B][0]);
        v[1] = - (B_normal[A][1] + B_normal[B][1]);
        v = v/norm2(v);

        // Insert vertex into FIFO queue
        dummy.first = V;
        dummy.second = v;
        Vqueue.push(dummy);

        // Locate bins for vertex V
        ix = (int) ((Tri.Vertex[V][0] - xlim[0])/dx);
        iy = (int) ((Tri.Vertex[V][1] - ylim[0])/dy);

        // Insert vertices into kd-tree
        if ((ix >= 0) && (ix <= nbins-1)
         && (iy >= 0) && (iy <= nbins-1)) {
            I = nbins*iy + ix;
            X_ = &Tri.Vertex[V];
            Ktree[I].insert(X_);
        }
    } //next T

}


//Tri.Export_vtu("vertex_cloud_0.vtu");

// ========================================================================== //
// ITERATIVELY INSERT VERTICES                                                //
// ========================================================================== //
{

    // Scope variables ------------------------------------------------------ //
    bool                                check;
    int                                 v_counter = Tri.nVertex, iter;
    int                                 W;
    int                                 V, I;
    int                                 di, dj, i, j, k, ix, iy;
    double                              h = 1.0/16.0;
    a3vector1D                          X, *X_ = &X;
    array<double, 3>                    P;
    array<double, 2>                    v;
    pair<int, array<double, 2> >        dummy;

    // Iteratively insert new vertices -------------------------------------- //
    iter = 0;
    v.fill(0.0);
    P.fill(0.0);
    while (!Vqueue.empty()) {
    // for (int ii = 0; ii < 4; ii++) {
        iter++;
        //cout << Vqueue.size() << endl;
        //cout << "------------------------------------" << endl;

        // Extract next item from FIFO queue
        dummy = Vqueue.front();

        // Retrieve item's info
        V = dummy.first;
        v = dummy.second;
        //cout << "vertex: " << V << " normal: " << v << endl;

        // Remove extracted element from queue
        Vqueue.pop();

        // Compute insertion depth
        // h = interp(h_field, Tri.Vertex[V]);

        // Check if candidate vertices are close to accepted vertices
        CGElem::RotateVector(v, -0.5*pi);
        for (k = 0; k < 3; k++) {

            // Compute coordinates for the new vertex
            //cout << "inserting candidate vertex " << k << " along " << v << endl;
            P[0] = Tri.Vertex[V][0] + h * v[0];
            X[0] = P[0];
            P[1] = Tri.Vertex[V][1] + h * v[1];
            X[1] = P[1];
            //cout << "  candidate point is: " << X << endl;

            // Bins which new vertex belongs to
            ix = (int) ((P[0] - xlim[0])/dx);
            iy = (int) ((P[1] - ylim[0])/dy);
            //cout << "  new vertex is included in bin: (" << ix << ", " << iy << ")" << endl;

            // Proximity check
            check = true;
            //cout << "h: " << h << ", dx: " << dx << ", dy: " << dy << endl;
            di = nbins; //((int) (h/dx)) + 1;
            dj = nbins; //((int) (h/dy)) + 1;
            // di = 5;
            // dj = 5;
            //cout << "di: " << di << ", dj: " << dj << endl;
            i = -di;
            while ((i <= di) && check) {
                j = -dj;
                while ((j <= dj) && check) {
                    if ((ix+i >= 0) && (ix+i <= nbins-1)
                     && (iy+j >= 0) && (iy+j <= nbins-1)) {
                        I = nbins*(iy+j) + ix+i;
                        W = Ktree[I].hNeighbor(X_, 0.5*h, (iter == 69));
                        check = (check && (W < 0));
                    }
                    j++;
                } //next j
                i++;
            } //next i
            if (iter == 69) {
                cout << "found: " << W << endl;
                cout << "  candidate point is: " << X << endl;
                cout << "  pointer points to: " << (*X_) << endl;
                for (int ii = 0; ii < Ktree[I].n_nodes; ++ii) {
                    if (norm2((*(Ktree[I].nodes[ii].object_)) - X) <= 0.75*h) {
                        cout << "dummy check node: " << ii << endl;
                    }
                }
                cout << "proximity check is: " << check << endl;
            }
            //if (!check) {
            //    cout << (*(Ktree[I].nodes[W].object_)) << " " << P << endl;
            //    cout << norm2((*(Ktree[I].nodes[W].object_)) - X) << " "  << 0.5*h << endl;
            //}

            // Check if candidate point is inside domain
            if (check) {
                //cout << "checking inclusion in domain" << endl;
                check = (check && Dpoly.IncludePoint(P));
                //if (!check) { cout << "  point is outside domain" << endl; }
            }

            // Check if candidate point is outside holes
            if (check) {
                H = 0;
                while (check && (H < nH)) {
                    //cout << "checking inclusion in hole: " << H << endl;
                    check = (check && (!Hpoly[H].IncludePoint(P)));
                    H++;
                } //next H
                //if (!check) { cout << "  point is inside hole: " << H-1 << endl; }
            }

            if (check) {

                // Insert candidate vertex into kd-tree
                //cout << "adding vertex to mesh data structure" << endl;
                Tri.AddVertex(X);
                // if ((ix >= 0) && (ix <= nbins-1)
                 // && (iy >= 0) && (iy <= nbins-1)) {
                    //cout << "inserting vertex into kd-tree" << endl;
                    I = nbins*(iy) + ix;
                    Ktree[I].insert(&(Tri.Vertex[Tri.nVertex-1]));
                // }

                // Insert candidate vertex into FIFO queue
                //cout << "inserting vertex into FIFO queue" << endl;
                dummy.first = Tri.nVertex-1;
                dummy.second = v;
                Vqueue.push(dummy);

                // Update counters
                v_counter++;
            }

            // Update insertion direction
            CGElem::RotateVector(v, 0.5*pi);
            
        } //next k

//         {
//             stringstream        nome;
//             nome << "vertex_cloud_" << iter << ".vtu";
//             Tri.Export_vtu(nome.str());
//         }

    } //next items
    
}

// Debug only (export vertex cloud into .vtk format) ------------------------ //
//Tri.Export_vtu("vertex_cloud.vtu");

return; };
