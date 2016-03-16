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
# include "Operators.hpp"
# include "CGBase.hpp"

// ========================================================================== //
// IMPLEMENTATIONS                                                            //
// ========================================================================== //

// Constructors ============================================================= //

// -------------------------------------------------------------------------- //
CGPolygon2D::Class_CG_Polygon2D::Class_CG_Polygon2D(
    void
) {

// ========================================================================== //
// CGPolygon2D::Class_CG_Polygon2D::Class_CG_Polygon2D(                      //
//     void)                                                                  //
//                                                                            //
// Default constructor for Class_CG_Polygon2D. Initialize an empty polygon.   //
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
// none

// Counters
// none

// ========================================================================== //
// INITIALIZE VARIABLES                                                       //
// ========================================================================== //

// Number of faces
n_faces = 0;

// infos
info.status = -2;
info.is_convex = false;
info.is_clockwise = false;
info.xlim.fill(0.0);
info.ylim.fill(0.0);

return; }

// -------------------------------------------------------------------------- //
CGPolygon2D::Class_CG_Polygon2D::Class_CG_Polygon2D(
    std::vector<std::array<double,3>>           &X
) {

// ========================================================================== //
// CGPolygon2D::Class_CG_Polygon2D::Class_CG_Polygon2D(                      //
//     dvector2D           &X)                                                //
//                                                                            //
// Custom constructor #1 for Class_CG_Polygon2D. Initialize a polygon with    //
// n faces (n = X.size()), whose vertices are specified in X. Vertex must     //
// be given in consecutive order. No repeated vertex are allowed.             //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - X          : dvector2D, vertex coordinate list. X[i][0], X[i][1] are the //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                 n = X.size();

// ========================================================================== //
// INITIALIZE DATA STRUCTURE                                                  //
// ========================================================================== //

// Number of faces
n_faces = n;

// Vertex list
v_list.resize(n_faces);
v_list = X ;

// Infos
if (n == 0) {
}
else {
    if (n < 3) {
        info.status = 0;
    }
    else {
        info.status = 1;
    }
}
IsClockWise();
IsConvex();
BoundingBox();

return; }

// Destructors ============================================================== //

// -------------------------------------------------------------------------- //
CGPolygon2D::Class_CG_Polygon2D::~Class_CG_Polygon2D(
    void
) {

// ========================================================================== //
// CGPolygon2D::Class_CG_Polygon2D::~Class_CG_Polygon2D(                     //
//     void)                                                                  //
//                                                                            //
// Default destructor for Class_CG_Polygon2D variables.                       //
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
// none

// Counters
// none

// ========================================================================== //
// CLEAN POLYGON DATA STRUCTURE                                               //
// ========================================================================== //

// General info
n_faces = 0;

// Vertex list
v_list.clear();

// Infos
info.status = -2;
info.is_convex = false;
info.is_clockwise = false;

return; }

// Assignement operators ==================================================== //
CGPolygon2D::Class_CG_Polygon2D& CGPolygon2D::Class_CG_Polygon2D::operator=(
    const CGPolygon2D::Class_CG_Polygon2D    &S
) {

// ========================================================================== //
// CG_Polygon2D::Class_CG_Polygon2D&                                          //
//      CG_Polygon2D::Class_CG_Polygon2D::operator=(                          //
//     const CG_Polygon2D::Class_CG_Polygon2D    &S)                          //
//                                                                            //
// Assignament operator for Class_CG_Polygon2D variables.                     //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - S    : Class_CG_Polygon2D, source polygon                                //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - this : &Class_CG_Polygon2D, copy of source variables                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// COPY S INTO *THIS                                                          //
// ========================================================================== //

// General info
n_faces = S.n_faces;

// Vertex list
v_list.resize(n_faces);
v_list = S.v_list;

// Infos
info.is_clockwise = S.info.is_clockwise;
info.is_convex = S.info.is_convex;
info.xlim = S.info.xlim;
info.ylim = S.info.ylim;
info.status = S.info.status;

return(*this); }

// Methods ================================================================== //

// -------------------------------------------------------------------------- //
void CGPolygon2D::Class_CG_Polygon2D::AddVertex(
    std::array<double,3>             &P
) {

// ========================================================================== //
// void CGPolygon2D::Class_CG_Polygon2D::AddVertex(                          //
//     dvector1D           &P)                                                //
//                                                                            //
// Add vertex to polygon vertex list. The new vertex is added at the end of   //
// the vertex list. Check for convexity is re-run.                            //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P      : dvector1D, new vertex coordinates                               //
// ========================================================================== //
// - OUTPUT                                                                   //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
//none

// Counters
// none

// ========================================================================== //
// ADD VERTEX TO VERTEX LIST                                                  //
// ========================================================================== //

// Update v_list
v_list.push_back(P);
n_faces++;

// Check status
if (info.status != -1) {
    if (n_faces < 3) {
        info.status = 0;
    }
    else {
        info.status = 1;
    }
}

// Update bounding box
info.xlim[0] = std::min(info.xlim[0], P[0]);
info.xlim[1] = std::max(info.xlim[1], P[0]);
info.ylim[0] = std::min(info.ylim[0], P[1]);
info.ylim[1] = std::max(info.ylim[1], P[1]);

// Check for convexity
IsClockWise();
IsConvex();

return; };

// -------------------------------------------------------------------------- //
void CGPolygon2D::Class_CG_Polygon2D::AddVertices(
    std::vector<std::array<double,3>>           &X
) {

// ========================================================================== //
// void CGPolygon2D::Class_CG_Polygon2D::AddVertices(                        //
//     dvector2D           &X)                                                //
//                                                                            //
// Add multiple vertices to polygon vertex list. New vertices are added at    //
// the end of the vertex list. Check for convexity is re-run.                 //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - X      : dvector2D, new vertex coordinates                               //
// ========================================================================== //
// - OUTPUT                                                                   //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
int                         n = X.size();

// Counters
int                         i;

// ========================================================================== //
// ADD VERTEX TO VERTEX LIST                                                  //
// ========================================================================== //

// Update v_list
//2array v_list.resize(n_faces+n);
//2array for (i = n_faces; i < n; ++i) {
//2array     V[0] = X[i][0];
//2array     V[1] = X[i][1];
//2array     v_list[i] = V;
//2array }

v_list.insert(v_list.end(), X.begin(), X.end()) ;
n_faces += n ;

// Check status
if (info.status != -1) {
    if (n_faces < 3) {
        info.status = 0;
    }
    else {
        info.status = 1;
    }
}

// Update bounding box
for (i = 0; i < n; ++i) {
    info.xlim[0] = std::min(info.xlim[0], X[i][0]);
    info.xlim[1] = std::max(info.xlim[1], X[i][0]);
    info.ylim[0] = std::min(info.ylim[0], X[i][1]);
    info.ylim[1] = std::max(info.ylim[1], X[i][1]);
} //next i

// Check for convexity
IsClockWise();
IsConvex();

return; };

// -------------------------------------------------------------------------- //
void CGPolygon2D::Class_CG_Polygon2D::BoundingBox(
    void
) {

// ========================================================================== //
// void CGPolygon2D::Class_CG_Polygon2D::BoundingBox(                        //
//     void)                                                                  //
//                                                                            //
// Set polygon bounding box.                                                  //
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
// none

// Counters
int                 i;

// ========================================================================== //
// CHECK POLYGON STATUS                                                       //
// ========================================================================== //
if (info.status < 0) { return; }

// ========================================================================== //
// COMPUTE LIMITS OF POLYGON BOUNDING BOX                                     //
// ========================================================================== //
info.xlim[0] = info.xlim[1] = v_list[0][0];
info.ylim[0] = info.ylim[1] = v_list[0][1];
for (i = 1; i < n_faces; ++i) {
    info.xlim[0] = std::min(info.xlim[0], v_list[i][0]);
    info.xlim[1] = std::max(info.xlim[1], v_list[i][0]);
    info.ylim[0] = std::min(info.ylim[0], v_list[i][1]);
    info.ylim[1] = std::max(info.ylim[1], v_list[i][1]);
} //next i

return; }

// -------------------------------------------------------------------------- //
void CGPolygon2D::Class_CG_Polygon2D::IsClockWise(
    void
) {

// ========================================================================== //
// void CGPolygon2D::Class_CG_Polygon2D::IsClockWise(                        //
//     void)                                                                  //
//                                                                            //
// Set flag 'is_clockwise' to 'true' if polygon is clockwise, to 'false'      //
// otherwise.                                                                 //
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
double              A = 0.0;
std::array<double, 3>    u, v;

// Counters
int                 i, j;

// ========================================================================== //
// CHECK POLYGON STATUS                                                       //
// ========================================================================== //
if (info.status <= 0) { return; }

// ========================================================================== //
// CHECK IF POLYGON IS CLOCKWISE                                              //
// ========================================================================== //
for (i = 1; i < n_faces; ++i) {
    j = (i + 1) % n_faces;
    u = v_list[i] - v_list[0];
    v = v_list[j] - v_list[0];
    A += norm2(crossProduct(u, v));
} //next i

info.is_clockwise = (A <= 0.0);

return; }

// -------------------------------------------------------------------------- //
void CGPolygon2D::Class_CG_Polygon2D::IsConvex(
    void
) {

// ========================================================================== //
// void CGPolygon2D::Class_CG_Polygon2D::IsConvex(                           //
//     void)                                                                  //
//                                                                            //
// Set flag 'is_convex' to true if polygon is convex, to 'false' otherwise.   //
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
double              s, f;
std::array<double, 3>    u, v;

// Counters
int                 i, j, k;

// ========================================================================== //
// CHECK POLYGON STATUS                                                       //
// ========================================================================== //
if (info.status <= 0) { return; }

// ========================================================================== //
// CHECK CONVEXITY                                                            //
// ========================================================================== //
if (info.is_clockwise) { f = -1.0; }
else                   { f = 1.0; }
info.is_convex = true;
j = 0;
while (info.is_convex && (j < n_faces)) {
    i = (j + n_faces - 1) % n_faces;
    k = (j + 1) % n_faces;
    u = v_list[j] - v_list[i];
    u = u/norm2(u);
    v = v_list[k] - v_list[j];
    v = v/norm2(v);
    s = norm2(crossProduct(v, u));
    info.is_convex = (f*s <= 0.0);
    j++;
} //next vertex

return; }

// -------------------------------------------------------------------------- //
void CGPolygon2D::Class_CG_Polygon2D::display(
    std::ostream             &out
) {

// ========================================================================== //
// void CGPolygon2D::Class_CG_Polygon2D::display(                            //
//     ostream             &out)                                              //
//                                                                            //
// Display polygon infos.                                                     //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - out      : ostream, output stream                                        //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - none                                                                     //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
int                     i;

// ========================================================================== //
// DISPLAY POLYGON INFOS                                                      //
// ========================================================================== //

// Header
out << "2D Polygon: " << std::endl;

// General infos
out << " n faces:        " << n_faces << std::endl;
out << " status:         ";
switch (info.status) {
    case(-2): { out << "polygon is empty" << std::endl; break;}
    case(-1): { out << "vertex list is empty" << std::endl; break;}
    case( 0): { out << "polygon is ill-defined" << std::endl; break;}
    case( 1): { out << "ok" << std::endl; break;}
};
if (info.status >= 0) {
    out << " polygon is:     ";
    if (info.is_convex) { out << "convex" << std::endl; }
    else           { out << "concave" << std::endl; }
    out << " polygon dir.:   ";
    if (info.is_clockwise) { out << "clockwise" << std::endl; }
    else              { out << "counter-clockwise" << std::endl; }
}
// bounding box
if (info.status == 1) {
    out << " polygon bounds: " << "x: " << info.xlim << ", y: " << info.ylim << std::endl;
}

// vertex list
if (info.status == 1) {
    out << " vertex list:    {";
    for (i = 0; i < n_faces-1; ++i) {
        out << "(" << v_list[i][0] << ", " << v_list[i][1] << ", " << v_list[i][2] <<"); ";
    }
    i = n_faces-1;
    out << "(" << v_list[i][0] << ", " << v_list[i][1] << ", " << v_list[i][2] << ")";
    out << "}" << std::endl;
}

return; }

// -------------------------------------------------------------------------- //
bool CGPolygon2D::Class_CG_Polygon2D::ConvexIncludePoint(
    std::array<double, 3>    &P
) {

// ========================================================================== //
// bool CGPolygon2D::Class_CG_Polygon2D::ConvexIncludePoint(                 //
//     std::array<double, 2>    &P)                                                //
//                                                                            //
// Check wheter the polygon includes a point P (Convex case). P must not lie  //
// on polygon edges.                                                          //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P     : dvector1D, point coordinates                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - check : bool, 'true' if P is inside the polygon, 'false' otherwise       //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                check = true;
double              s, f;
std::array<double, 3>    u, v;

// Counters
int                 i, j;

// ========================================================================== //
// CHECK IF POLYGON INCLUDES POINT P                                          //
// ========================================================================== //

// Check with bounding box
if ((P[0] < info.xlim[0])
 || (P[0] > info.xlim[1])
 || (P[1] < info.ylim[0])
 || (P[1] > info.ylim[1])) { return(false); }

// Check using cross productu rule
if (info.is_clockwise) { f = -1.0; }
else                { f = 1.0; }
i = 0;
while (check && (i < n_faces)) {
    j = (i + 1) % n_faces;
    u = v_list[j] - v_list[i];
    u = u/norm2(u);
    v = P - v_list[i];
    v = v/norm2(v);
    s = norm2(crossProduct(u, v));
    check = (f*s >= 0.0);
    i++;
}

return(check); };

// -------------------------------------------------------------------------- //
bool CGPolygon2D::Class_CG_Polygon2D::NonConvexIncludePoint(
    std::array<double, 3>    &P
) {

// ========================================================================== //
// bool CGPolygon2D::Class_CG_Polygon2D::NonConvexIncludePoint(              //
//     std::array<double, 2>    &P)                                                //
//                                                                            //
// Check wheter the polygon includes a point P (non-Convex case). P must not  //
// lie on polygon edges.                                                      //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P     : dvector1D, point coordinates                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - check : bool, 'true' if P is inside the polygon, 'false' otherwise       //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                check;
int                 wn = 0;
double              s;
std::array<double, 3>    u, v;
// Counters
int                 i, j;

// ========================================================================== //
// CHECK IF POLYGON INCLUDES POINT P                                          //
// ========================================================================== //

// Check with bounding box
if ((P[0] < info.xlim[0])
 || (P[0] > info.xlim[1])
 || (P[1] < info.ylim[0])
 || (P[1] > info.ylim[1])) { return(false); }
 
// Check with winding number rule
for (i = 0; i < n_faces; ++i) {
    j = (i + 1) % n_faces;
    if (abs(v_list[i][1] - v_list[j][1]) >= 1.0e-12) {
        u = v_list[j] - v_list[i];
        u = u/norm2(u);
        v = P - v_list[i];
        v = v/norm2(v);
        s = norm2(crossProduct(u, v) );
        if ((v_list[i][1] < P[1]) && (v_list[j][1] >= P[1])) {
            if (s >= 0.0) {
                wn++;
            }
        }
        else if ((v_list[i][1] > P[1]) && (v_list[j][1] <= P[1])) {
            if (s <= 0.0) {
                wn--;
            }
        }
    }
} //next i
check = (wn != 0);

return(check); };

// -------------------------------------------------------------------------- //
bool CGPolygon2D::Class_CG_Polygon2D::IncludePoint(
    std::array<double, 3>    &P
) {

// ========================================================================== //
// bool CGPolygon2D::Class_CG_Polygon2D::IncludePoint(                       //
//     std::array<double, 2>    &P)                                                //
//                                                                            //
// Check wheter the polygon includes a point P. P must not lie on polygon     //
// edges.                                                                     //
// ========================================================================== //
// INPUT                                                                      //
// ========================================================================== //
// - P     : dvector1D, point coordinates                                     //
// ========================================================================== //
// OUTPUT                                                                     //
// ========================================================================== //
// - check : bool, 'true' if P is inside the polygon, 'false' otherwise       //
// ========================================================================== //

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
bool                check;

// Counters
// none

// ========================================================================== //
// CHECK INCLUSION                                                            //
// ========================================================================== //
if (info.is_convex) {
    check = ConvexIncludePoint(P);
}
else {
    check = NonConvexIncludePoint(P);
}


return(check); }
