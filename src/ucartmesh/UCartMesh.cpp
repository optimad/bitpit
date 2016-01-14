/*!
 *  \ingroup    UCartMesh
 *  @{
 */

#include "UCartMesh.hpp"

# include "Operators.hpp"
# include "Class_VTK.hpp"

using namespace std;

// ========================================================================== //
// METHODS IMPLEMENTATIONS FOR UCartMesh                              //
// ========================================================================== //

// Constructors ------------------------------------------------------------- //

/* -------------------------------------------------------------------------- */
/*!
 *   Default constructor creates an empty 3D mesh 
 */
UCartMesh::UCartMesh( ){

    dim = 3 ;

    // Mesh extent
    B0.fill(0.0) ;
    B1.fill(0.0) ;

    // Mesh size
    nc.fill(0);

    center.resize(3) ;
    edge.resize(3) ;

    // Mesh spacing
    h.fill(0.0);


    whichDirection[0] = 0 ;
    whichDirection[1] = 0 ;
    whichDirection[2] = 1 ;
    whichDirection[3] = 1 ;
    whichDirection[4] = 2 ;
    whichDirection[5] = 2 ;

    whichStep[0] = -1 ;
    whichStep[1] = +1 ;
    whichStep[2] = -1 ;
    whichStep[3] = +1 ;
    whichStep[4] = -1 ;
    whichStep[5] = +1 ;

    status = 0;

};

/* -------------------------------------------------------------------------- */
/*!
 *   Generic constructor 
 *   \param[in] P0 min coordinate of mesh
 *   \param[in] P1 max coordinate of mesh
 *   \param[in] N number of cell in each direction
 *   \param[in] dimension number of space dimensions [2/3]
 */
UCartMesh::UCartMesh( darray3E const &P0, darray3E const &P1, iarray3E const &N, int const dimension) :UCartMesh(){

    setMesh( P0, P1, N, dimension);
    status = 0 ;

};

/* -------------------------------------------------------------------------- */
/*!
 *   2D mesh constructor 
 *   \param[in] P0 min coordinate of mesh
 *   \param[in] P1 max coordinate of mesh
 *   \param[in] I number of cell in first direction
 *   \param[in] J number of cell in second direction
 */
UCartMesh::UCartMesh( darray3E const &P0, darray3E const &P1, int const &I, int const &J) :UCartMesh(){

    iarray3E    N={I,J,1} ;     

    setMesh( P0, P1, N, 2);
    status = 0 ;

};

/* -------------------------------------------------------------------------- */
/*!
 *   3D mesh constructor 
 *   \param[in] P0 min coordinate of mesh
 *   \param[in] P1 max coordinate of mesh
 *   \param[in] I number of cell in first direction
 *   \param[in] J number of cell in second direction
 *   \param[in] K number of cell in third direction
 */
UCartMesh::UCartMesh( darray3E const &P0, darray3E const &P1, int const &I, int const &J, int const &K) :UCartMesh(){

    iarray3E    N={I,J,K} ;     

    setMesh( P0, P1, N, 3);
    status = 0 ;

};

/* -------------------------------------------------------------------------- */
/*!
 *   Destructor 
 */
UCartMesh::~UCartMesh( ){

    for( int d=0; d<3; ++d){
        dvector1D().swap(center[d]);
        dvector1D().swap(edge[d]);
    };


};

/* -------------------------------------------------------------------------- */
/*!
 *   Assignment operator
 */
UCartMesh& UCartMesh::operator=(
        const UCartMesh &B
        ) {


    // Number of cells
    nc = B.nc;
    np = B.np;

    h = B.h;

    CellsInIJPlane = B.CellsInIJPlane ;
    NodesInIJPlane = B.NodesInIJPlane ;

    nCells = B.nCells ;
    nNodes = B.nNodes ;

    // Mesh limits
    B0 = B.B0;
    B1 = B.B1;

    // Resize mesh data structure ----------------------------------------------- //
    ResizeMesh();

    // Copy cell edges and cell centers ----------------------------------------- //
    edge = B.edge;
    center = B.center;

    return(*this); 
};

/* -------------------------------------------------------------------------- */
/*!  Set new mesh
 *   \param[in]     A0      mesh min point
 *   \param[in]     A1      mesh max point
 *   \param[in]     N       number of cells in each direction
 *   \param[in]     dims    number of cells in each direction
 */
void UCartMesh::setMesh( darray3E const & A0, darray3E const & A1, iarray3E const & N, int const &dims ){

    // Counters
    int       i, d;

    dim             = dims ;
    nc              = N ;
    np              = N + 1;

    if( dim ==2 ){
        nc[2] = 1 ;
        np[2] = 1 ;
    };

    nCells          = nc[0]*nc[1]*nc[2] ;
    nNodes         = np[0]*np[1]*np[2] ;
    CellsInIJPlane  = nc[0]*nc[1] ;
    NodesInIJPlane = np[0]*np[1] ;


    // Mesh limits
    B0 = A0 ;
    B1 = A1 ;

    // Number of mesh cells

    // Resize mesh data structure ----------------------------------------------- //
    ResizeMesh();

    // Create mesh -------------------------------------------------------------- //

    // Mesh spacing
    for( d=0; d<dim; ++d){
        h[d] = (B1[d] - B0[d])/((double) nc[d]);
    };

    // vetices
    for( d=0; d<dim; ++d){

        for( i=0; i<np[d]; ++i){
            edge[d][i] = B0[d] + ((double) i) * h[d];
        };

    };


    // Cells centers
    for( d=0; d<dim; ++d){

        for( i=0; i<nc[d]; ++i){
            center[d][i] = edge[d][i] + 0.5 *h[d] ;
        };
    };

    status++ ;

    return; 

};

/* -------------------------------------------------------------------------- */
/*!  Get status of mesh; staus is increased each time mesh is modified
 *   \return void
 */
void UCartMesh::ClearMesh( ){

    // Mesh limits
    B0.fill(0.0) ;
    B1.fill(0.0) ;

    h.fill(0.0) ;

    // Number of cells
    nc.fill(0) ;
    np.fill(0) ;

    // Resize mesh data structure
    ResizeMesh();

    status++ ;

    return; 
}

/* -------------------------------------------------------------------------- */
/*!  Resize of data structures of UCartMesh
 */
void UCartMesh::ResizeMesh( ){

    int d ;  

    for( d=0; d<3; ++d){
        center[d].resize( nc[d], 0.0 ) ;
        edge[d].resize( np[d], 0.0);
    };

    return; 
};

/* -------------------------------------------------------------------------- */
/*!  Get number of cells 
 *   \return number of cells in mesh
 */
int UCartMesh::getNCells(){

    return nCells ;
};

/* -------------------------------------------------------------------------- */
/*!  Get number of cells in one direction
 *   \param[in] d cartesian direction
 *   \return number of cells in direction
 */
int UCartMesh::getNCells( int d){

    return nc[d] ;
};

/* -------------------------------------------------------------------------- */
/*!  Get number of nodes  
 *   \return number of points in mesh
 */
int UCartMesh::getNNodes(){

    return nNodes  ;
};

/* -------------------------------------------------------------------------- */
/*!  Get number of points in one direction
 *   \param[in] d cartesian direction
 *   \return number of points in direction
 */
int UCartMesh::getNNodes( int d){

    return np[d]  ;
};

/* -------------------------------------------------------------------------- */
/*! Get mesh spacing
 *  \return mesh spacing in all directions 
 */
darray3E    UCartMesh::getSpacing(){
    return h ;
};

/* -------------------------------------------------------------------------- */
/*! Get mesh spacing in one direction
 *  \param[in]  d   direction
 *  \return mesh spacing in direction 
 */
double      UCartMesh::getSpacing( int d){
    return h[d] ;
};

/* -------------------------------------------------------------------------- */
/*! Get mesh dimension
 *  \return mesh dimension
 */
int      UCartMesh::getDimension( ){
    return dim ;
};

/* -------------------------------------------------------------------------- */
/*!  Get axis aligned bounding box of mesh
 *   \param[out] A0 min point
 *   \param[out] A1 max point
 */
void UCartMesh::getBoundingBox( darray3E &A0, darray3E &A1){

    A0 = B0; 
    A1 = B1; 

    return ;
}

/* -------------------------------------------------------------------------- */
/*!  Get status of mesh; staus is increased each time mesh is modified
 *   \return status
 */
int UCartMesh::getStatus(  ){

    return  status ;

};

/* -------------------------------------------------------------------------- */
/*! Tanslate mesh by given offset
 *  \param[in]  ds      offset
 */
void UCartMesh::Translate( darray3E const &ds ){ 

    // Counters
    int     d;

    // Mesh limits
    B0 = B0 + ds;
    B1 = B1 + ds;

    // Cells edges
    for( d=0; d<dim; ++d){
        edge[d] = edge[d] + ds[d] ;
        center[d] = center[d] + ds[d] ;
    };

    status++ ;

    return; 
};

/* -------------------------------------------------------------------------- */
/*! Scale mesh with respect to given origin by given factor 
 *  \param[in]  s       scale factor
 *  \param[in]  origin  scale origin
 */
void UCartMesh::Scale( darray3E const &s, darray3E const &origin ){

    // Counters
    int     d;


    h = h *s;

    // Mesh limits
    for( d=0; d<dim; ++d){
        B0[d] = origin[d] +s[d] *( B0[d] - origin[d]) ;
        B1[d] = origin[d] +s[d] *( B1[d] - origin[d]) ;

        edge[d] = origin[d] +s[d] *( edge[d] - origin[d]) ;

        center[d] = origin[d] +s[d] *( center[d] - origin[d]) ;
    };

    status++ ;

    return; 

};

/* -------------------------------------------------------------------------- */
/*! Scale mesh with respect to its min coordinate by given factor 
 *  \param[in]  s       scale factor
 */
void UCartMesh::Scale( darray3E const &s ){

    Scale( s, B0);

    return; 

};

/* -------------------------------------------------------------------------- */
/*! Compute the cartesian indices of cell which contains given point.
 *  If point is outside mesh closest cell is returned
 *  \param[in]  P       point coordinates
 *  \return             cartesian cell indices
 */
iarray3E UCartMesh::CellCartesianId( darray3E const &P ){

    int         d ;
    iarray3E    id;

    id.fill(0) ;

    for( d=0; d<dim; ++d){
        id[d] = min( nc[d]-1, max(0, (int) floor( (P[d] - B0[d])/h[d] )) );
    };

    return id; 
};

/* -------------------------------------------------------------------------- */
/*! Compute the cartesian indices of cell which contains given point in a 2D Mesh.
 *  If point is outside mesh closest cell is returned
 *  \param[in]  P       point coordinates
 *  \param[out]  i      first cartesian coordinate
 *  \param[out]  j      second cartesian coordinate
 */
void UCartMesh::CellCartesianId( darray3E const &P, int &i, int &j ){

    i = min( nc[0]-1, max(0, (int) floor( (P[0] - B0[0])/h[0] )) );
    j = min( nc[1]-1, max(0, (int) floor( (P[1] - B0[1])/h[1] )) );

};

/* -------------------------------------------------------------------------- */
/*! Compute the cartesian indices of cell which contains given point in a 3D Mesh.
 *  If point is outside mesh closest cell is returned
 *  \param[in]  P       point coordinates
 *  \param[out]  i       first cartesian coordinate
 *  \param[out]  j       second cartesian coordinate
 *  \param[out]  k       third cartesian coordinate
 */
void UCartMesh::CellCartesianId( darray3E const &P, int &i, int &j, int &k ){

    i = min( nc[0]-1, max(0, (int) floor( (P[0] - B0[0])/h[0] )) );
    j = min( nc[1]-1, max(0, (int) floor( (P[1] - B0[1])/h[1] )) );
    k = min( nc[2]-1, max(0, (int) floor( (P[2] - B0[2])/h[2] )) );

};

/* -------------------------------------------------------------------------- */
/*! Transformation from linear index to cartesian indices for generic meshes.
 * No check on bounds is performed
 *  \param[in]  J       linear cell index
 *  \return             cell cartesien indices
 */
iarray3E UCartMesh::CellCartesianId( int const &J ){

    iarray3E    id ;

    id[0] = J % nc[0] ;
    id[2] = J / CellsInIJPlane;
    id[1] =  (J - id[2] *CellsInIJPlane ) /nc[0]  ;


    return id; 
};

/* -------------------------------------------------------------------------- */
/*! Transformation from linear index to cartesian indices for 2D meshes.
 * No check on bounds is performed
 *  \param[in]  J       linear cell index
 *  \param[out]  i       first cartesian coordinate
 *  \param[out]  j       second cartesian coordinate
 */
void UCartMesh::CellCartesianId( int const &J, int &i, int &j ){

    i = J % nc[0] ;
    j = J /nc[0]  ;


    return ; 
};

/* -------------------------------------------------------------------------- */
/*! Transformation from linear index to cartesian indices for 3D meshes
 * No check on bounds is performed
 *  \param[in]  J       linear cell index
 *  \param[out]  i       first cartesian coordinate
 *  \param[out]  j       second cartesian coordinate
 *  \param[out]  k       third cartesian coordinate
 */
void UCartMesh::CellCartesianId( int const &J, int &i, int &j, int &k ){

    i = J % nc[0] ;
    k = J / CellsInIJPlane;
    j =  (J - k *CellsInIJPlane ) /nc[0]  ;

    return ; 
};

/* -------------------------------------------------------------------------- */
/*! Compute the linear index of cell which contains given point.
 *  If point is outside mesh closest cell is returned
 *  \param[in]  P       point coordinates
 *  \return             linear cell index
 */
int UCartMesh::CellLinearId( darray3E const &P ){

    // Counters
    int         n ;
    iarray3E    id;


    id = CellCartesianId(P) ;
    n  = CellLinearId(id) ;

    return n; 

};

/* -------------------------------------------------------------------------- */
/*! Transformation from cartesian indices to linear indices.
 * No check on bounds is performed
 *  \param[in]  i       first cartesian index
 *  \param[in]  j       second cartesian index
 *  \param[in]  k       third cartesian index
 *  \return             linear cell index
 */
int UCartMesh::CellLinearId( int const &i, int const &j ,int const &k ){

    return( CellsInIJPlane*k + nc[0]*j +i );

};

/* -------------------------------------------------------------------------- */
/*! Transformation from cartesian indices to linear indices.
 * No check on bounds is performed
 *  \param[in]  id      cartesian indices
 *  \return             linear cell index
 */
int UCartMesh::CellLinearId( iarray3E const &id ){

    return( CellLinearId(id[0], id[1], id[2] ) ); 
};

/* -------------------------------------------------------------------------- */
/*!  Get cell center coordinates through cartesian indices
 *   \param[in] i first index
 *   \param[in] j second index
 *   \param[in] k third index (0 for 2D mesh)
 *   \return cell center 
 */
darray3E UCartMesh::getCellCenter( int i, int j, int k ){

    darray3E    P ;

    P[0]= center[0][i] ;
    P[1]= center[1][j] ;
    P[2]= center[2][k] ;

    return P;

} ;

/* -------------------------------------------------------------------------- */
/*!  Get cell center coordinates through cartesian indices
 *   \param[in] id array with cartesian indices
 *   \return cell center 
 */
darray3E UCartMesh::getCellCenter( iarray3E id ){

    darray3E    P ;

    P[0]= center[0][id[0]] ;
    P[1]= center[1][id[1]] ;
    P[2]= center[2][id[2]] ;

    return P;

} ;

/* -------------------------------------------------------------------------- */
/*!  Get cell center coordinates through linear index
 *   \param[in] J linear index
 *   \return cell center 
 */
darray3E UCartMesh::getCellCenter( int J ){

    return  getCellCenter( CellCartesianId(J) );

};

/* -------------------------------------------------------------------------- */
/*!  Get min and max node coordinatess of cell through cartesian indices
 *   \param[in] i first cartesian index
 *   \param[in] j secon cartesian index
 *   \param[out] C0 min node coordinates
 *   \param[out] C1 max node coordinates
 */
void UCartMesh::getCellBoundingBox( int const &i, int const &j, darray3E &C0, darray3E &C1 ){

    C0 =  getNodeCoordinates( i, j );
    C1 =  getNodeCoordinates( i+1, j+1 );

};

/* -------------------------------------------------------------------------- */
/*!  Get min and max node coordinatess of cell through cartesian indices
 *   \param[in] i first cartesian index
 *   \param[in] j secon cartesian index
 *   \param[in] k third cartesian index 
 *   \param[out] C0 min node coordinates
 *   \param[out] C1 max node coordinates
 */
void UCartMesh::getCellBoundingBox( int const &i, int const &j, int const &k, darray3E &C0, darray3E &C1 ){

    iarray3E    id({ i, j, k}) ;

    getCellBoundingBox( id, C0, C1) ;

};

/* -------------------------------------------------------------------------- */
/*!  Get min and max node coordinates of cell through cartesian indices
 *   \param[in] id array of cartesian indices
 *   \param[out] C0 min node coordinates
 *   \param[out] C1 max node coordinates
 */
void UCartMesh::getCellBoundingBox( iarray3E const &id0, darray3E &C0, darray3E &C1 ){


    int         d;
    iarray3E    id1(id0) ;

    for( d=0; d<dim; ++d){
        id1[d]++ ;
    };

    C0 =  getNodeCoordinates( id0 );
    C1 =  getNodeCoordinates( id1 );

    return  ;


};

/* -------------------------------------------------------------------------- */
/*!  Get min and max node coordinates of cell through linear index
 *   \param[in] J linear index
 *   \param[out] C0 min node coordinates
 *   \param[out] C1 max node coordinates
 */
void UCartMesh::getCellBoundingBox( int const &J, darray3E &C0, darray3E &C1 ){

    int         d ;
    iarray3E    id0, id1 ;

    id0 = CellCartesianId(J) ;

    for( d=0; d<dim; ++d){
        id1[d] = id0[d] +1 ;
    };

    C0 =  getNodeCoordinates( id0 );
    C1 =  getNodeCoordinates( id1 );

    return  ;

};

/* -------------------------------------------------------------------------- */
/*!  Get inode coordinates through cartesian indices
 *   \param[in] i first index
 *   \param[in] j second index
 *   \param[in] k third index (0 for 2D mesh)
 *   \return node coordinates
 */
darray3E UCartMesh::getNodeCoordinates( int i, int j, int k ){

    darray3E    P ;

    P[0]= edge[0][i] ;
    P[1]= edge[1][j] ;
    P[2]= edge[2][k] ;

    return P;

} ;

/* -------------------------------------------------------------------------- */
/*!  Get node coordinates through cartesian indices
 *   \param[in] id array with cartesian indices
 *   \return node coordinates
 */
darray3E UCartMesh::getNodeCoordinates( iarray3E id ){

    darray3E    P ;

    P[0]= edge[0][id[0]] ;
    P[1]= edge[1][id[1]] ;
    P[2]= edge[2][id[2]] ;

    return P;

} ;

/* -------------------------------------------------------------------------- */
/*!  Get node coordinates through linear index
 *   \param[in] J linear index
 *   \return node coordinates
 */
darray3E UCartMesh::getNodeCoordinates( int J ){

    return  getNodeCoordinates( NodeCartesianId(J) );

};

/* -------------------------------------------------------------------------- */
/*! Get linear index of closest point neighbour in given direction.
 *  If neighbour does not exist boundary point is returned
 *  \param[in]  I       reference point index
 *  \param[in]  dir     direction [0/1/2/3/4/5] 
 *  \return neighbour index
 */
int      UCartMesh::getNodeNeighbour( int const &I, int const &dir){

    int     d = whichDirection[dir], step = whichStep[dir];

    return  getNodeNeighbour(I,d,step) ;
};

/* -------------------------------------------------------------------------- */
/*! Get linear index of point neighbour in given direction and distance.
 *  If neighbour does not exist boundary point is returned
 *  \param[in]  I       reference point index
 *  \param[in]  d       direction [0/1/2] 
 *  \param[in]  step    distance
 *  \return neighbour index
 */
int      UCartMesh::getNodeNeighbour( int const &I, int const &d, int const &step){

    iarray3E    i;

    i       = NodeCartesianId(I) ;
    i[d]    += step ; 
    i[d]    = max( min( i[d], nc[d] ), 0 ) ;

    return  NodeLinearId(i) ;
};

/* -------------------------------------------------------------------------- */
/*! Get linear index of closest cell neighbour in given direction.
 *  If neighbour does not exist boundary cell is returned
 *  \param[in]  I       reference point index
 *  \param[in]  dir     direction [0/1/2/3/4/5] 
 *  \return neighbour index
 */
int      UCartMesh::getCellNeighbour( int const &I, int const &dir){

    int         d = whichDirection[dir], step = whichStep[dir];

    return  getCellNeighbour(I,d,step) ;
};

/* -------------------------------------------------------------------------- */
/*! Get linear index of cell neighbour in given direction and distance.
 *  If neighbour does not exist boundary cell is returned
 *  \param[in]  I       reference point index
 *  \param[in]  d       direction [0/1/2] 
 *  \param[in]  step    distance
 *  \return neighbour index
 */
int      UCartMesh::getCellNeighbour( int const &I, int const &d, int const &step){

    iarray3E    i;

    i       = CellCartesianId(I) ;
    i[d]    += step ; 
    i[d]    = max( min( i[d], nc[d]-1 ), 0 ) ;

    return  CellLinearId(i) ;
};

/* -------------------------------------------------------------------------- */
/*! Compute the cartesian indices of closest node to given point.
 *  \param[in]  P       point coordinates
 *  \return             cartesian indices
 */
iarray3E UCartMesh::NodeCartesianId( darray3E const &P ){

    int         d ;
    iarray3E    id;

    id.fill(0) ;

    for( d=0; d<dim; ++d){
        id[d] = min( np[d]-1, max(0, (int) round( (P[d] - B0[d])/h[d] )) );
    };

    return id; 
};

/* -------------------------------------------------------------------------- */
/*! Compute the cartesian indices of closest node to given point in a 2D Mesh
 *  \param[in]  P       point coordinates
 *  \param[out]  i      first cartesian index
 *  \param[out]  j      second cartesian index
 */
void UCartMesh::NodeCartesianId( darray3E const &P, int &i, int &j ){


    i = min( np[0]-1, max(0, (int) round( (P[0] - B0[0])/h[0] )) );
    j = min( np[1]-1, max(0, (int) round( (P[1] - B0[1])/h[1] )) );

    return ;

};

/* -------------------------------------------------------------------------- */
/*! Compute the cartesian indices of closest node to given point in a 3D Mesh
 *  \param[in]  P       point coordinates
 *  \param[out]  i      first cartesian index
 *  \param[out]  j      second cartesian index
 *  \param[out]  k      third cartesian index
 */
void UCartMesh::NodeCartesianId( darray3E const &P, int &i, int &j, int &k ){


    i = min( np[0]-1, max(0, (int) round( (P[0] - B0[0])/h[0] )) );
    j = min( np[1]-1, max(0, (int) round( (P[1] - B0[1])/h[1] )) );
    k = min( np[2]-1, max(0, (int) round( (P[2] - B0[2])/h[2] )) );

    return ;

};

/* -------------------------------------------------------------------------- */
/*! Transformation from linear index to cartesian indices.
 * No check on bounds is performed
 *  \param[in]  J       node linear index
 *  \return             node cartesian indices
 */
iarray3E UCartMesh::NodeCartesianId( int const &J ){

    // Local variables
    iarray3E    id ;

    id[0] = J % np[0] ;
    id[2] = J / NodesInIJPlane;
    id[1] =  (J - id[2] *NodesInIJPlane ) /np[0]  ;

    return id; 

};

/* -------------------------------------------------------------------------- */
/*! Transformation from linear index to cartesian indices in a 3D mesh
 * No check on bounds is performed
 *  \param[in]  J       node linear index
 *  \param[out]  i      first cartesian  index
 *  \param[out]  j      second cartesian  index
 */
void UCartMesh::NodeCartesianId( int const &J, int &i, int &j ){

    i = J %np[0] ;
    j = J /np[0]  ;

    return ; 

};

/* -------------------------------------------------------------------------- */
/*! Transformation from linear index to cartesian indices in a 3D mesh
 * No check on bounds is performed
 *  \param[in]  J       node linear index
 *  \param[out]  i      first cartesian  index
 *  \param[out]  j      second cartesian  index
 *  \param[out]  k      third cartesian  index
 */
void UCartMesh::NodeCartesianId( int const &J, int &i, int &j, int &k ){

    i = J % np[0] ;
    k = J / NodesInIJPlane;
    j =  (J - k *NodesInIJPlane ) /np[0]  ;

    return ; 

};

/* -------------------------------------------------------------------------- */
/*! Transformation from cartesian indices to linear indices.
 * No check on bounds is performed
 *  \param[in]  i       first cartesian index
 *  \param[in]  j       second cartesian index
 *  \param[in]  k       third cartesian index
 *  \return             linear index of node
 */
int UCartMesh::NodeLinearId( int const &i, int const &j, int const &k ){

    return (NodesInIJPlane *k + np[0]*j + i); 

};

/* -------------------------------------------------------------------------- */
/*! Transformation from cartesian indices to linear indices.
 * No check on bounds is performed
 *  \param[in]  id      node cartesian indices
 *  \return             node linear index
 */
int UCartMesh::NodeLinearId( darray3E const &P ){

    // Counters
    iarray3E    id( NodeCartesianId(P) );

    return NodeLinearId(id); 

};

/* -------------------------------------------------------------------------- */
/*! Transformation from cartesian indices to linear indices.
 * No check on bounds is performed
 *  \param[in]  id      node cartesian indices
 *  \return             node linear index
 */
int UCartMesh::NodeLinearId( iarray3E const &id){

    return ( NodeLinearId(id[0], id[1], id[2] ) ); 

};

/* -------------------------------------------------------------------------- */
/*! Calculate cell subset indices form cartesian indices
 *  \param[in]  i0      min cartesian indices
 *  \param[in]  i1      max cartesian indices
 *  \return             cell linear indices of subset mesh
 */
ivector1D UCartMesh::CellSubSet( iarray3E const &i0, iarray3E const &i1 ){

    int                     i, j, k; 
    ivector1D               ids;
    ivector1D::iterator     it;

    i  =  i1[0]-i0[0]+1  ;
    j  =  i1[1]-i0[1]+1  ;
    k  =  i1[2]-i0[2]+1  ;

    i  =  i *j *k ;
    ids.resize(i) ;

    it = ids.begin() ;

    for( k=i0[2]; k<=i1[2]; ++k){
        for( j=i0[1]; j<=i1[1]; ++j){
            for( i=i0[0]; i<=i1[0]; ++i){

                *it = CellLinearId( i, j, k) ;            
                ++it ;

            };
        };
    };

    return ids; 

};

/* -------------------------------------------------------------------------- */
/*! Calculate cell subset indices form linear indices
 *  \param[in]  I0      min linear indices
 *  \param[in]  I1      max linear indices
 *  \return             cell linear indices of subset mesh
 */
ivector1D UCartMesh::CellSubSet( int const &I0, int const &I1 ){

    return CellSubSet( CellCartesianId(I0), CellCartesianId(I1) ); 

};

/* -------------------------------------------------------------------------- */
/*! Calculate cell subset indices form min and max point.
 *  The cell conataining the points (or closest to them) are maintained
 *  \param[in]  P0      min point
 *  \param[in]  P1      max point
 *  \return             cell linear indices of subset mesh
 */
ivector1D UCartMesh::CellSubSet( darray3E const &P0, darray3E const &P1 ){

    return CellSubSet( CellCartesianId(P0), CellCartesianId(P1) ); 

};

/* -------------------------------------------------------------------------- */
/*! Calculate subset indices form cartesian indices
 *  \param[in]  i0      min cartesian indices
 *  \param[in]  i1      max cartesian indices
 *  \return             node linear indices of subset mesh
 */
ivector1D UCartMesh::NodeSubSet( iarray3E const &i0, iarray3E const &i1 ){

    int                     i, j, k; 
    ivector1D               ids;
    ivector1D::iterator     it;

    i  =  i1[0]-i0[0]+1  ;
    j  =  i1[1]-i0[1]+1  ;
    k  =  i1[2]-i0[2]+1  ;

    i  =  i *j *k ;
    ids.resize(i) ;

    it = ids.begin() ;

    for( k=i0[2]; k<=i1[2]; ++k){
        for( j=i0[1]; j<=i1[1]; ++j){
            for( i=i0[0]; i<=i1[0]; ++i){

                *it = NodeLinearId( i, j, k) ;            
                ++it ;

            };
        };
    };

    return ids; 

};

/* -------------------------------------------------------------------------- */
/*! Calculate node subset indices form linear indices
 *  \param[in]  I0      min linear indices
 *  \param[in]  I1      max linear indices
 *  \return             node linear indices of subset mesh
 */
ivector1D UCartMesh::NodeSubSet( int const &I0, int const &I1 ){

    return NodeSubSet( NodeCartesianId(I0), NodeCartesianId(I1) ); 

};

/* -------------------------------------------------------------------------- */
/*! Calculate node subset indices form min and max point.
 *  The nodes closest to the points are used as limites 
 *  \param[in]  P0      min point
 *  \param[in]  P1      max point
 *  \return             cell linear indices of subset mesh
 */
ivector1D UCartMesh::NodeSubSet( darray3E const &P0, darray3E const &P1 ){

    return NodeSubSet( NodeCartesianId(P0), NodeCartesianId(P1) ); 

};

/* -------------------------------------------------------------------------- */
/*! Check if point lies within the mesh
 *  \param[in]  P      
 *  \return             true if point lies within grid
 */
bool UCartMesh::PointInGrid(  darray3E const &P ){

    int     d;

    for( d=0; d<dim; ++d){

        if( P[d]< B0[d] || P[d] > B1[d] ){
            return false;
        };
    };

    return true ;
};

/* -------------------------------------------------------------------------- */
/*! Check if point lies within the mesh
 *  \param[in]  P      
 *  \param[out] I       cartesian indices of cell containing the point      
 *  \return             true if point lies within grid
 */
bool UCartMesh::PointInGrid( darray3E const &P, iarray3E &I){

    int     d;

    for( d=0; d<dim; ++d){

        if( P[d]< B0[d] || P[d] > B1[d] ){
            return false;
        };
    };

    I = CellCartesianId(P) ;

    return true ;
};

/* -------------------------------------------------------------------------- */
/*! Check if point lies within the mesh
 *  \param[in]  P      
 *  \param[out] i       first cartesian index of cell containing the point      
 *  \param[out] j       second cartesian index of cell containing the point      
 *  \param[out] k       third cartesian index of cell containing the point      
 *  \return             true if point lies within grid
 */
bool UCartMesh::PointInGrid( darray3E const &P, int &i, int &j, int &k ){

    bool        inGrid ;
    iarray3E    I;

    if(PointInGrid(P,I) ){

        i = I[0] ;
        j = I[1] ;
        k = I[2] ;

        return true ;
    } ;

    return false ;
};

/* -------------------------------------------------------------------------- */
/*! Check if point lies within the mesh
 *  \param[in]  P      
 *  \param[out] I       linear index of cell containing the point      
 *  \return             true if point lies within grid
 */
bool UCartMesh::PointInGrid( darray3E const &P, int &I ){

    int     d;

    for( d=0; d<dim; ++d){

        if( P[d]< B0[d] || P[d] > B1[d] ){
            return false;
        };
    };

    I = CellLinearId(P) ;

    return true ;
};

/* -------------------------------------------------------------------------- */
/*!  Save cartesian mesh as unstructured grid
 *  \param[out]     nV  number of vertices
 *  \param[out]     nS  number of cells
 *  \param[out]     V   vertex coordinates
 *  \param[out]     S   cell-vertex connectivity
 *  \param[out]     A   cell-cell adjacencies
 */
void UCartMesh::Cart2Unstr( int &nV, int &nS, vector<darray3E> &V, ivector2D &S, ivector3D &A ){

    // Local variables
    int         nv, ns;

    // Counters
    int         i, j, k, J;

    // Number of new vertices/simplicies
    nv = nNodes ;
    ns = nCells ;

    // Resize vertex list
    V.resize(nV + nv);

    // Resize simplex list
    S.resize(nS + ns, ivector1D(pow(2,dim), -1));

    // Resize adjacency
    A.resize(nS + ns, ivector2D(2*dim, ivector1D(1, -1)));

    // ========================================================================== //
    // GENERATE VERTEX LIST                                                       //
    // ========================================================================== //
    for (k = 0; k < np[2]; k++) {
        for (j = 0; j < np[1]; j++) {
            for (i = 0; i < np[0]; i++) {
                J = NodeLinearId(i,j,k);
                V[J][0] = edge[0][i] ;
                V[J][1] = edge[1][j] ;
                V[J][2] = edge[2][k] ; 
                nV++;
            };
        }
    } 

    // ========================================================================== //
    // SIMPLEX-VERTEX CONNECTIVITY                                                //
    // ========================================================================== //
    for (k = 0; k < nc[2]; k++) {
        for (j = 0; j < nc[1]; j++) {
            for (i = 0; i < nc[0]; i++) {
                J = CellLinearId(i,j,k);

                S[J][0] = NodeLinearId(i,j,k);
                S[J][1] = NodeLinearId(i+1,j,k);
                S[J][2] = NodeLinearId(i,j+1,k);
                S[J][3] = NodeLinearId(i+1,j+1,k);

                if(dim==3){
                    S[J][4] = NodeLinearId(i,j,k+1);
                    S[J][5] = NodeLinearId(i+1,j,k+1);
                    S[J][6] = NodeLinearId(i,j+1,k+1);
                    S[J][7] = NodeLinearId(i+1,j+1,k+1);
                }

            }
        }
    } 

    // ========================================================================== //
    // SIMPLEX-VERTEX ADJACENCY                                                   //
    // ========================================================================== //
    for (k = 0; k < nc[2]; k++) {
        for (j = 0; j < nc[1]; j++) {
            for (i = 0; i < nc[0]; i++) {
                J = CellLinearId(i,j);

                if (i != 0)     { A[J][0][0] = CellLinearId(i-1,j,k); }
                if (i != nc[0]) { A[J][1][0] = CellLinearId(i+1,j,k); }

                if (j != 0)     { A[J][2][0] = CellLinearId(i,j-1,k); }
                if (j != nc[1]) { A[J][3][0] = CellLinearId(i,j+1,k); }

                if (k != 0)     { A[J][4][0] = CellLinearId(i,j,k-1); }
                if (k != nc[2]) { A[J][5][0] = CellLinearId(i,j,k+1); }

            }
        }
    }

    return; 
}

/* -------------------------------------------------------------------------- */
/*!  Save mesh in .vtr fie format (paraview)
 *  \param[in]     dir  relative or absolute path to directory with final "/"
 *  \param[in]     filename  name of output file without suffix
 */
void UCartMesh::ExportVtr(string dir, string filename) {

    // Local variables

    class VTKOut : public VTK_RectilinearGrid<VTKOut>{
        private:
            dvector2D   *ptredge ;

        public:
            VTKOut(string dir, string filename, dvector2D &edges) :VTK_RectilinearGrid<VTKOut>( ){ 
                SetNames( dir, filename) ;
                SetCodex( "appended") ;
                ptredge = &edges ;
                SetDimensions( 0, (*ptredge)[0].size()-1, 0, (*ptredge)[1].size()-1, 0, (*ptredge)[2].size()-1 ) ;

            } ;

            void Flush(  fstream &str, string codex_, string name  ){

                if( name == "x_Coord"){
                    flush_binary( str, (*ptredge)[0] ) ;
                }

                else if( name == "y_Coord"){
                    flush_binary( str, (*ptredge)[1] ) ;
                }

                else if( name == "z_Coord"){
                    flush_binary( str, (*ptredge)[2] ) ;
                };

                return ;
            };

    };

    VTKOut      myVTK( dir, filename, edge) ;
    myVTK.Write() ;

    return; 
};

//* -------------------------------------------------------------------------- */
/*!  Save mesh and floating data in .vtr file format (paraview)
 *  \param[in]     dir  relative or absolute path to directory with final "/"
 *  \param[in]     filename  name of output file without suffix
 *  \param[in]     dataname  name of data field
 *  \param[in]     location  Cell or point data ["Cell"/"Point"]
 *  \param[in]     data      data field
 */
void UCartMesh::ExportVtr(string dir, string filename, string dataname, string location, vector<double> &data ) {


    class VTKOut : public VTK_RectilinearGrid<VTKOut>{
        private:
            dvector2D   *ptredge ;
            dvector1D   *ptrdata ;

        public:
            VTKOut(string dir, string filename, dvector2D &edges, string dataname, string location, vector<double> &data) :VTK_RectilinearGrid<VTKOut>( ){ 
                SetNames( dir, filename) ;
                SetCodex( "appended") ;
                ptredge = &edges ;
                SetDimensions( 0, (*ptredge)[0].size()-1, 0, (*ptredge)[1].size()-1, 0, (*ptredge)[2].size()-1 ) ;

                ptrdata = &data ;
                AddData( dataname, 1, "Float64", location ) ;

            } ;

            void Flush(  fstream &str, string codex_, string name  ){

                if( name == "x_Coord"){
                    flush_binary( str, (*ptredge)[0] ) ;
                }

                else if( name == "y_Coord"){
                    flush_binary( str, (*ptredge)[1] ) ;
                }

                else if( name == "z_Coord"){
                    flush_binary( str, (*ptredge)[2] ) ;
                }
                
                else{
                    flush_binary( str, (*ptrdata) ) ;
                };

                return ;
            };

    };

    VTKOut      myVTK( dir, filename, edge, dataname, location, data) ;
    myVTK.Write() ;

    return; 

};

//* -------------------------------------------------------------------------- */
/*!  Save mesh and integer data in .vtr file format (paraview)
 *  \param[in]     dir  relative or absolute path to directory with final "/"
 *  \param[in]     filename  name of output file without suffix
 *  \param[in]     dataname  name of data field
 *  \param[in]     location  Cell or point data ["Cell"/"Point"]
 *  \param[in]     data      data field
 */
void UCartMesh::ExportVtr(string dir, string filename, string dataname, string location, ivector1D &data ) {


    class VTKOut : public VTK_RectilinearGrid<VTKOut>{
        private:
            dvector2D   *ptredge ;
            ivector1D   *ptrdata ;

        public:
            VTKOut(string dir, string filename, dvector2D &edges, string dataname, string location, ivector1D &data) :VTK_RectilinearGrid<VTKOut>( ){ 
                SetNames( dir, filename) ;
                SetCodex( "appended") ;
                ptredge = &edges ;
                SetDimensions( 0, (*ptredge)[0].size()-1, 0, (*ptredge)[1].size()-1, 0, (*ptredge)[2].size()-1 ) ;

                ptrdata = &data ;
                AddData( dataname, 1, "UInt32", location, "appended" ) ;

            } ;

            void Flush(  fstream &str, string codex_, string name  ){

                if( name == "x_Coord"){
                    flush_binary( str, (*ptredge)[0] ) ;
                }

                else if( name == "y_Coord"){
                    flush_binary( str, (*ptredge)[1] ) ;
                }

                else if( name == "z_Coord"){
                    flush_binary( str, (*ptredge)[2] ) ;
                }
                
                else{
                    flush_binary( str, (*ptrdata) ) ;
                };

                return ;
            };

    };

    VTKOut      myVTK( dir, filename, edge, dataname, location, data) ;
    myVTK.Write() ;

    return; 

};

//* -------------------------------------------------------------------------- */
/*!  Save mesh and vector field data in .vtr file format (paraview)
 *  \param[in]     dir  relative or absolute path to directory with final "/"
 *  \param[in]     filename  name of output file without suffix
 *  \param[in]     dataname  name of data field
 *  \param[in]     location  Cell or point data ["Cell"/"Point"]
 *  \param[in]     data      data field
 */
void UCartMesh::ExportVtr(string dir, string filename, string dataname, string location, vector<array<double,3>> &data ) {


    class VTKOut : public VTK_RectilinearGrid<VTKOut>{
        private:
            dvector2D                   *ptredge ;
             vector<array<double,3>>    *ptrdata ;

        public:
            VTKOut(string dir, string filename, dvector2D &edges, string dataname, string location, vector<array<double,3>> &data) :VTK_RectilinearGrid<VTKOut>( ){ 
                SetNames( dir, filename) ;
                SetCodex( "appended") ;
                ptredge = &edges ;
                SetDimensions( 0, (*ptredge)[0].size()-1, 0, (*ptredge)[1].size()-1, 0, (*ptredge)[2].size()-1 ) ;

                ptrdata = &data ;
                AddData( dataname, 3, "Float64", location, "appended" ) ;

            } ;

            void Flush(  fstream &str, string codex_, string name  ){

                if( name == "x_Coord"){
                    flush_binary( str, (*ptredge)[0] ) ;
                }

                else if( name == "y_Coord"){
                    flush_binary( str, (*ptredge)[1] ) ;
                }

                else if( name == "z_Coord"){
                    flush_binary( str, (*ptredge)[2] ) ;
                }
                
                else{
                    flush_binary( str, (*ptrdata) ) ;
                };

                return ;
            };

    };

    VTKOut      myVTK( dir, filename, edge, dataname, location, data) ;
    myVTK.Write() ;

    return; 

};


/* -------------------------------------------------------------------------- */
/*!  Transform cell data to point data by calculating the mean of incident cells in each vertex
 *  \param[in]     CellData  Data on cells
 *  \param[out]    NodeData  Data on nodes
 */
void UCartMesh::CellData2NodeData( vector<double> &CellData, vector<double> &NodeData ){
    // ========================================================================== //

    // Local variables
    int            K, J;
    int            ip, jp, kp;
    vector<int>     NodeIter ;

    // Counters
    int            i, j, k, l, m, n;

    // ========================================================================== //
    // RESIZE OUTPUT VARIABLES                                                    //
    // ========================================================================== //
    NodeData.resize( getNNodes() );
    NodeIter.resize( getNNodes() );

    // ========================================================================== //
    // CONVERT CELL DATA INTO POINT DATA                                          //
    // ========================================================================== //
    for (k = 0; k < nc[2]; k++) {
        for (j = 0; j < nc[1]; j++) {
            for (i = 0; i < nc[0]; i++) {

                // Cell index
                K = CellLinearId(i, j, k);

                for (n = 0; n < dim-1; n++) {
                    for (m = 0; m < 2; m++) {
                        for (l = 0; l < 2; l++) {

                            // Point index
                            ip = i + l;
                            jp = j + m;
                            kp = j + n;
                            J = NodeLinearId(ip,jp,kp);

                            NodeData[J] = NodeData[J] + CellData[K]; 
                            NodeIter[J]++ ;
                        } //next n
                    } //next m
                } //next l

            } //next k
        } //next j
    } //next i


    for( J=0; J<nNodes; ++J){
        NodeData[J] = NodeData[J] / ((float) NodeIter[J]) ;
    };

    return; 

};

/* -------------------------------------------------------------------------- */
/*!  Transform node data to cell data by calculating the mean of incident vertices of each cell
 *  \param[in]    PointData  Data on nodes
 *  \param[out]   CellData  Data on cells
 */
void UCartMesh::NodeData2CellData( vector<double> &NodeData, vector<double> &CellData ){


    // Local variables
    int     K, J;
    int     ip, jp, kp;
    double  factor ;

    // Counters
    int     i, j, k, l, m, n;

    factor  =   pow(0.5,dim) ;

    // ========================================================================== //
    // RESIZE OUTPUT VARIABLES                                                    //
    // ========================================================================== //
    CellData.resize(getNCells(), 0.0);

    // ========================================================================== //
    // CONVERT POINT DATA TO CELL DATA                                            //
    // ========================================================================== //
    for (k = 0; k < nc[2]; k++) {
        for (j = 0; j < nc[1]; j++) {
            for (i = 0; i < nc[0]; i++) {

                K = CellLinearId(i,j,k);
                for (n = 0; n < 2; n++) {
                    for (m = 0; m < 2; m++) {
                        for (l = 0; l < 2; l++) {
                            ip = i + l;
                            jp = j + m;
                            kp = k + n;
                            J = NodeLinearId(ip,jp,kp);
                            CellData[K] = CellData[K] + factor * NodeData[J];
                        } //next n
                    } //next m
                } //next l
            } //next k
        } //next j
    } //next i

    return; 
};

/* -------------------------------------------------------------------------- */
/*! Calculates bi-/ tri- linear interpolation stencil on cells for a given point.
 * At boundaries stencil is reduced to assure positive weights.
 * If the point is outside a null-stencil is returned
 *  \param[in]      P       Point coordinates
 *  \param[out]     stencil linear indices of the interpolation stencil
 *  \param[out]     weights weights associated to stencil
 *  \return         number of cells used in the interpolation stencil. 0 id point outside the grid
 */
int UCartMesh::linearCellInterpolation( darray3E &P, vector<int> &stencil, vector<double> &weights ){

    int                         nStencil(0) ;
    int                         d, i, j, k;
    iarray3E                    i0, i1 ;
    
    array< array<int,2>, 3>     cStencil ;
    array< array<double,2>, 3>  cWeights ;

    array<int,3>                nS({1,1,1}) ;

    vector<int>::iterator       itrStencil ;
    vector<double>::iterator    itrWeights ;

    if( PointInGrid(P, i0) ){

        nStencil = 1 ;

        for( d=0; d<dim; ++d){

            // Find cell index
            if( P[d] < center[d][i0[d] ] ){
                i0[d] = i0[d]-1 ;
            };

            i1[d] = i0[d] +1 ;

            if( i0[d] < 0 ){
                cStencil[d][0] = 0 ; 
                cWeights[d][0] = 1. ; 

            } 

            else if( i1[d] > nc[d]-1 ){
                cStencil[d][0] = nc[d]-1 ; 
                cWeights[d][0] = 1. ; 

            } 

            else{
                nS[d] = 2 ;

                cStencil[d][0] = i0[d] ;
                cStencil[d][1] = i1[d] ;

                cWeights[d][1] = (P[d] - center[d][i0[d]]) /h[d]   ;
                cWeights[d][0] = 1.0 - cWeights[d][1] ;  
            }


        };

        for( d=dim; d<3; ++d){
            cStencil[d][0] = 0 ;
            cWeights[d][0] = 1. ;
        };

        nStencil = nS[0] *nS[1] *nS[2] ;

        stencil.resize( nStencil ) ;
        weights.resize( nStencil ) ;

        itrStencil = stencil.begin() ;
        itrWeights = weights.begin() ;

        for( k=0; k<nS[2]; ++k){
            for( j=0; j<nS[1]; ++j){
                for( i=0; i<nS[0]; ++i){

                    int &is = cStencil[0][i] ;
                    int &js = cStencil[1][j] ;
                    int &ks = cStencil[2][k] ;

                    double &iw = cWeights[0][i] ;
                    double &jw = cWeights[1][j] ;
                    double &kw = cWeights[2][k] ;

                    *itrStencil =  CellLinearId(is,js,ks) ;
                    *itrWeights =  iw *jw *kw ;

                    ++itrStencil ;
                    ++itrWeights ;


                }
            }
        }

    }

    return nStencil; 

};

/* -------------------------------------------------------------------------- */
/*! Calculates bi-/ tri- linear interpolation stencil on nodes for a given point.
 * If the point is outside grid a null-stencil is returned
 *  \param[in]      P       Point coordinates
 *  \param[out]     stencil linear indices of the interpolation stencil
 *  \param[out]     weights weights associated to stencil
 *  \return         number of cells used in the interpolation stencil. 0 id point outside the grid
 */
int UCartMesh::linearNodeInterpolation( darray3E &P, vector<int> &stencil, vector<double> &weights ){


    int                         nStencil(0) ;
    int                         d, i, j, k;
    iarray3E                    i0, i1 ;

    array< array<int,2>, 3>     cStencil ;
    array< array<double,2>, 3>  cWeights ;

    vector<int>::iterator       itrStencil ;
    vector<double>::iterator    itrWeights ;

    if( PointInGrid(P,i0) ){

        nStencil = pow(2,dim) ;

        for(d=0; d<dim; ++d){
            i1[d] = i0[d] +1 ;

            cStencil[d][0] = i0[d] ;
            cStencil[d][1] = i1[d] ;

            cWeights[d][1] = ( P[d] - edge[d][i0[d]]) /h[d]  ;
            cWeights[d][0] = 1.0 - cWeights[d][1] ;  
        };

        for( d=dim; d<3; ++d){
            cStencil[d][0] = 0 ;
            cWeights[d][0] = 1. ;
        };

        stencil.resize(nStencil) ;
        weights.resize(nStencil) ;

        itrStencil = stencil.begin() ;
        itrWeights = weights.begin() ;

        for( k=0; k<dim-1; ++k){
            for( j=0; j<2; ++j){
                for( i=0; i<2; ++i){

                    int &is = cStencil[0][i] ;
                    int &js = cStencil[1][j] ;
                    int &ks = cStencil[2][k] ;

                    double &iw = cWeights[0][i] ;
                    double &jw = cWeights[1][j] ;
                    double &kw = cWeights[2][k] ;

                    *itrStencil =  NodeLinearId(is,js,ks) ;
                    *itrWeights =  iw *jw *kw ;

                    ++itrStencil ;
                    ++itrWeights ;


                }
            }
        }

    };


    return nStencil; 

};


/* @} */
