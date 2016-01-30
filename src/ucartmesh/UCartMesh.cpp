/*!
 *  \ingroup    UCartMesh
 *  @{
 */

# include <bitpit_operators.hpp>
# include <bitpit_IO.hpp>

#include "UCartMesh.hpp"

//using namespace std;

// ========================================================================== //
// METHODS IMPLEMENTATIONS FOR UCartMesh                              //
// ========================================================================== //

// Constructors ------------------------------------------------------------- //

/* -------------------------------------------------------------------------- */
/*!
 *   Default constructor creates an empty 3D mesh 
 */
UCartMesh::UCartMesh( ){

    m_dim = 3 ;

    // Mesh extent
    m_B0.fill(0.0) ;
    m_B1.fill(0.0) ;

    // Mesh size
    m_nc.fill(0);
    m_np.fill(0);

    m_center.resize(3) ;
    m_edge.resize(3) ;

    // Mesh spacing
    m_h.fill(0.0);


    m_whichDirection[0] = 0 ;
    m_whichDirection[1] = 0 ;
    m_whichDirection[2] = 1 ;
    m_whichDirection[3] = 1 ;
    m_whichDirection[4] = 2 ;
    m_whichDirection[5] = 2 ;

    m_whichStep[0] = -1 ;
    m_whichStep[1] = +1 ;
    m_whichStep[2] = -1 ;
    m_whichStep[3] = +1 ;
    m_whichStep[4] = -1 ;
    m_whichStep[5] = +1 ;

    m_status = 0;

};

/* -------------------------------------------------------------------------- */
/*!
 *   Generic constructor 
 *   \param[in] P0 min coordinate of mesh
 *   \param[in] P1 max coordinate of mesh
 *   \param[in] N number of cell in each direction
 *   \param[in] dimension number of space dimensions [2/3]
 */
UCartMesh::UCartMesh( std::array<double,3> const &P0, std::array<double,3> const &P1, std::array<int,3> const &N, int const dimension) :UCartMesh(){

    setMesh( P0, P1, N, dimension);
    m_status = 0 ;

};

/* -------------------------------------------------------------------------- */
/*!
 *   2D mesh constructor 
 *   \param[in] P0 min coordinate of mesh
 *   \param[in] P1 max coordinate of mesh
 *   \param[in] I number of cell in first direction
 *   \param[in] J number of cell in second direction
 */
UCartMesh::UCartMesh( std::array<double,3> const &P0, std::array<double,3> const &P1, int const &I, int const &J) :UCartMesh(){

    std::array<int,3>    N={I,J,1} ;     

    setMesh( P0, P1, N, 2);
    m_status = 0 ;

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
UCartMesh::UCartMesh( std::array<double,3> const &P0, std::array<double,3> const &P1, int const &I, int const &J, int const &K) :UCartMesh(){

    std::array<int,3>    N={I,J,K} ;     

    setMesh( P0, P1, N, 3);
    m_status = 0 ;

};

/* -------------------------------------------------------------------------- */
/*!
 *   Destructor 
 */
UCartMesh::~UCartMesh( ){

    for( int d=0; d<3; ++d){
        std::vector<double>().swap(m_center[d]);
        std::vector<double>().swap(m_edge[d]);
    };


};

/* -------------------------------------------------------------------------- */
/*!
 *   Assignment operator
 */
UCartMesh& UCartMesh::operator=(
        const UCartMesh &B
        ) {

    setMesh( B.m_B0, m_B1, B.m_nc, B.m_dim) ;

    // Copy cell m_edges and cell m_centers ----------------------------------------- //
    m_edge = B.m_edge;
    m_center = B.m_center;

    return(*this); 
};

/* -------------------------------------------------------------------------- */
/*!  Set new mesh
 *   \param[in]     A0      mesh min point
 *   \param[in]     A1      mesh max point
 *   \param[in]     N       number of cells in each direction
 *   \param[in]     dims    number of cells in each direction
 */
void UCartMesh::setMesh( std::array<double,3> const & A0, std::array<double,3> const & A1, std::array<int,3> const & N, int const &dims ){

    // Counters
    int       i, d;

    m_dim             = dims ;
    m_nc              = N ;
    m_np              = N + 1;

    if( m_dim ==2 ){
        m_nc[2] = 1 ;
        m_np[2] = 1 ;
    };

    // Number of mesh cells
    m_nCells          = m_nc[0]*m_nc[1]*m_nc[2] ;
    m_nNodes         = m_np[0]*m_np[1]*m_np[2] ;
    m_CellsInIJPlane  = m_nc[0]*m_nc[1] ;
    m_NodesInIJPlane = m_np[0]*m_np[1] ;


    // Mesh limits
    m_B0 = A0 ;
    m_B1 = A1 ;


    // Resize mesh data structure ----------------------------------------------- //
    for( d=0; d<3; ++d){
        m_center[d].resize( m_nc[d], 0.0 ) ;
        m_edge[d].resize( m_np[d], 0.0);
    };

    // Create mesh -------------------------------------------------------------- //

    // Mesh spacing
    for( d=0; d<m_dim; ++d){
        m_h[d] = (m_B1[d] - m_B0[d])/((double) m_nc[d]);
    };

    // vetices
    for( d=0; d<m_dim; ++d){

        for( i=0; i<m_np[d]; ++i){
            m_edge[d][i] = m_B0[d] + ((double) i) * m_h[d];
        };

    };

    // Cells m_centers
    for( d=0; d<m_dim; ++d){

        for( i=0; i<m_nc[d]; ++i){
            m_center[d][i] = m_edge[d][i] + 0.5 *m_h[d] ;
        };
    };

    m_status++ ;

    return; 

};

/* -------------------------------------------------------------------------- */
/*!  Clears the mesh data
 *   \return void
 */
void UCartMesh::clearMesh( ){

    // Mesh limits
    m_B0.fill(0.0) ;
    m_B1.fill(0.0) ;

    m_h.fill(0.0) ;

    // Number of cells
    m_nc.fill(0) ;
    m_np.fill(0) ;

    // Resize mesh data structure
    for( int d=0; d<3; ++d){
        m_center[d].clear() ;
        m_edge[d].clear() ;
    };

    m_status++ ;

    return; 
}

/* -------------------------------------------------------------------------- */
/*!  Get number of cells 
 *   \return number of cells in mesh
 */
int UCartMesh::getNCells(){

    return m_nCells ;
};

/* -------------------------------------------------------------------------- */
/*!  Get number of cells in one direction
 *   \param[in] d cartesian direction
 *   \return number of cells in direction
 */
int UCartMesh::getNCells( int d){

    return m_nc[d] ;
};

/* -------------------------------------------------------------------------- */
/*!  Get number of nodes  
 *   \return number of points in mesh
 */
int UCartMesh::getNNodes(){

    return m_nNodes  ;
};

/* -------------------------------------------------------------------------- */
/*!  Get number of points in one direction
 *   \param[in] d cartesian direction
 *   \return number of points in direction
 */
int UCartMesh::getNNodes( int d){

    return m_np[d]  ;
};

/* -------------------------------------------------------------------------- */
/*! Get mesh spacing
 *  \return mesh spacing in all directions 
 */
std::array<double,3>    UCartMesh::getSpacing(){
    return m_h ;
};

/* -------------------------------------------------------------------------- */
/*! Get mesh spacing in one direction
 *  \param[in]  d   direction
 *  \return mesh spacing in direction 
 */
double      UCartMesh::getSpacing( int d){
    return m_h[d] ;
};

/* -------------------------------------------------------------------------- */
/*! Get mesh dimension
 *  \return mesh dimension
 */
int      UCartMesh::getDimension( ){
    return m_dim ;
};

/* -------------------------------------------------------------------------- */
/*!  Get axis aligned bounding box of mesh
 *   \param[out] A0 min point
 *   \param[out] A1 max point
 */
void UCartMesh::getBoundingBox( std::array<double,3> &A0, std::array<double,3> &A1){

    A0 = m_B0; 
    A1 = m_B1; 

    return ;
}

/* -------------------------------------------------------------------------- */
/*!  Get status of mesh; staus is increased each time mesh is modified
 *   \return status
 */
int UCartMesh::getStatus(  ){

    return  m_status ;

};

/* -------------------------------------------------------------------------- */
/*! Tanslate mesh by given offset
 *  \param[in]  ds      offset
 */
void UCartMesh::translate( std::array<double,3> const &ds ){ 

    // Counters
    int     d;

    // Mesh limits
    m_B0 = m_B0 + ds;
    m_B1 = m_B1 + ds;

    // Cells m_edges
    for( d=0; d<m_dim; ++d){
        m_edge[d] = m_edge[d] + ds[d] ;
        m_center[d] = m_center[d] + ds[d] ;
    };

    m_status++ ;

    return; 
};

/* -------------------------------------------------------------------------- */
/*! scale mesh with respect to given origin by given factor 
 *  \param[in]  s       scale factor
 *  \param[in]  origin  scale origin
 */
void UCartMesh::scale( std::array<double,3> const &s, std::array<double,3> const &origin ){

    // Counters
    int     d;


    m_h = m_h *s;

    // Mesh limits
    for( d=0; d<m_dim; ++d){
        m_B0[d] = origin[d] +s[d] *( m_B0[d] - origin[d]) ;
        m_B1[d] = origin[d] +s[d] *( m_B1[d] - origin[d]) ;

        m_edge[d] = origin[d] +s[d] *( m_edge[d] - origin[d]) ;

        m_center[d] = origin[d] +s[d] *( m_center[d] - origin[d]) ;
    };

    m_status++ ;

    return; 

};

/* -------------------------------------------------------------------------- */
/*! scale mesh with respect to its min coordinate by given factor 
 *  \param[in]  s       scale factor
 */
void UCartMesh::scale( std::array<double,3> const &s ){

    scale( s, m_B0);

    return; 

};

/* -------------------------------------------------------------------------- */
/*! Compute the cartesian indices of cell which contains given point.
 *  If point is outside mesh closest cell is returned
 *  \param[in]  P       point coordinates
 *  \return             cartesian cell indices
 */
std::array<int,3> UCartMesh::getCellCartesianId( std::array<double,3> const &P ){

    int         d ;
    std::array<int,3>    id;

    id.fill(0) ;

    for( d=0; d<m_dim; ++d){
        id[d] = std::min( m_nc[d]-1, std::max(0, (int) floor( (P[d] - m_B0[d])/m_h[d] )) );
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
void UCartMesh::getCellCartesianId( std::array<double,3> const &P, int &i, int &j ){

    i = std::min( m_nc[0]-1, std::max(0, (int) floor( (P[0] - m_B0[0])/m_h[0] )) );
    j = std::min( m_nc[1]-1, std::max(0, (int) floor( (P[1] - m_B0[1])/m_h[1] )) );

};

/* -------------------------------------------------------------------------- */
/*! Compute the cartesian indices of cell which contains given point in a 3D Mesh.
 *  If point is outside mesh closest cell is returned
 *  \param[in]  P       point coordinates
 *  \param[out]  i       first cartesian coordinate
 *  \param[out]  j       second cartesian coordinate
 *  \param[out]  k       third cartesian coordinate
 */
void UCartMesh::getCellCartesianId( std::array<double,3> const &P, int &i, int &j, int &k ){

    i = std::min( m_nc[0]-1, std::max(0, (int) floor( (P[0] - m_B0[0])/m_h[0] )) );
    j = std::min( m_nc[1]-1, std::max(0, (int) floor( (P[1] - m_B0[1])/m_h[1] )) );
    k = std::min( m_nc[2]-1, std::max(0, (int) floor( (P[2] - m_B0[2])/m_h[2] )) );

};

/* -------------------------------------------------------------------------- */
/*! Transformation from linear index to cartesian indices for generic meshes.
 * No check on bounds is performed
 *  \param[in]  J       linear cell index
 *  \return             cell cartesien indices
 */
std::array<int,3> UCartMesh::getCellCartesianId( int const &J ){

    std::array<int,3>    id ;

    id[0] = J % m_nc[0] ;
    id[2] = J / m_CellsInIJPlane;
    id[1] =  (J - id[2] *m_CellsInIJPlane ) /m_nc[0]  ;


    return id; 
};

/* -------------------------------------------------------------------------- */
/*! Transformation from linear index to cartesian indices for 2D meshes.
 * No check on bounds is performed
 *  \param[in]  J       linear cell index
 *  \param[out]  i       first cartesian coordinate
 *  \param[out]  j       second cartesian coordinate
 */
void UCartMesh::getCellCartesianId( int const &J, int &i, int &j ){

    i = J % m_nc[0] ;
    j = J /m_nc[0]  ;


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
void UCartMesh::getCellCartesianId( int const &J, int &i, int &j, int &k ){

    i = J % m_nc[0] ;
    k = J / m_CellsInIJPlane;
    j =  (J - k *m_CellsInIJPlane ) /m_nc[0]  ;

    return ; 
};

/* -------------------------------------------------------------------------- */
/*! Compute the linear index of cell which contains given point.
 *  If point is outside mesh closest cell is returned
 *  \param[in]  P       point coordinates
 *  \return             linear cell index
 */
int UCartMesh::getCellLinearId( std::array<double,3> const &P ){

    // Counters
    int         n ;
    std::array<int,3>    id;


    id = getCellCartesianId(P) ;
    n  = getCellLinearId(id) ;

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
int UCartMesh::getCellLinearId( int const &i, int const &j ,int const &k ){

    return( m_CellsInIJPlane*k + m_nc[0]*j +i );

};

/* -------------------------------------------------------------------------- */
/*! Transformation from cartesian indices to linear indices.
 * No check on bounds is performed
 *  \param[in]  id      cartesian indices
 *  \return             linear cell index
 */
int UCartMesh::getCellLinearId( std::array<int,3> const &id ){

    return( getCellLinearId(id[0], id[1], id[2] ) ); 
};

/* -------------------------------------------------------------------------- */
/*!  Get cell m_center coordinates through cartesian indices
 *   \param[in] i first index
 *   \param[in] j second index
 *   \param[in] k third index (0 for 2D mesh)
 *   \return cell m_center 
 */
std::array<double,3> UCartMesh::getCellCenter( int i, int j, int k ){

    std::array<double,3>    P ;

    P[0]= m_center[0][i] ;
    P[1]= m_center[1][j] ;
    P[2]= m_center[2][k] ;

    return P;

} ;

/* -------------------------------------------------------------------------- */
/*!  Get cell m_center coordinates through cartesian indices
 *   \param[in] id array with cartesian indices
 *   \return cell m_center 
 */
std::array<double,3> UCartMesh::getCellCenter( std::array<int,3> id ){

    std::array<double,3>    P ;

    P[0]= m_center[0][id[0]] ;
    P[1]= m_center[1][id[1]] ;
    P[2]= m_center[2][id[2]] ;

    return P;

} ;

/* -------------------------------------------------------------------------- */
/*!  Get cell m_center coordinates through linear index
 *   \param[in] J linear index
 *   \return cell m_center 
 */
std::array<double,3> UCartMesh::getCellCenter( int J ){

    return  getCellCenter( getCellCartesianId(J) );

};

/* -------------------------------------------------------------------------- */
/*!  Get min and max node coordinatess of cell through cartesian indices
 *   \param[in] i first cartesian index
 *   \param[in] j secon cartesian index
 *   \param[out] C0 min node coordinates
 *   \param[out] C1 max node coordinates
 */
void UCartMesh::getCellBoundingBox( int const &i, int const &j, std::array<double,3> &C0, std::array<double,3> &C1 ){

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
void UCartMesh::getCellBoundingBox( int const &i, int const &j, int const &k, std::array<double,3> &C0, std::array<double,3> &C1 ){

    std::array<int,3>    id({ i, j, k}) ;

    getCellBoundingBox( id, C0, C1) ;

};

/* -------------------------------------------------------------------------- */
/*!  Get min and max node coordinates of cell through cartesian indices
 *   \param[in] id array of cartesian indices
 *   \param[out] C0 min node coordinates
 *   \param[out] C1 max node coordinates
 */
void UCartMesh::getCellBoundingBox( std::array<int,3> const &id0, std::array<double,3> &C0, std::array<double,3> &C1 ){


    int         d;
    std::array<int,3>    id1(id0) ;

    for( d=0; d<m_dim; ++d){
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
void UCartMesh::getCellBoundingBox( int const &J, std::array<double,3> &C0, std::array<double,3> &C1 ){

    int         d ;
    std::array<int,3>    id0, id1 ;

    id0 = getCellCartesianId(J) ;

    for( d=0; d<m_dim; ++d){
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
std::array<double,3> UCartMesh::getNodeCoordinates( int i, int j, int k ){

    std::array<double,3>    P ;

    P[0]= m_edge[0][i] ;
    P[1]= m_edge[1][j] ;
    P[2]= m_edge[2][k] ;

    return P;

} ;

/* -------------------------------------------------------------------------- */
/*!  Get node coordinates through cartesian indices
 *   \param[in] id array with cartesian indices
 *   \return node coordinates
 */
std::array<double,3> UCartMesh::getNodeCoordinates( std::array<int,3> id ){

    std::array<double,3>    P ;

    P[0]= m_edge[0][id[0]] ;
    P[1]= m_edge[1][id[1]] ;
    P[2]= m_edge[2][id[2]] ;

    return P;

} ;

/* -------------------------------------------------------------------------- */
/*!  Get node coordinates through linear index
 *   \param[in] J linear index
 *   \return node coordinates
 */
std::array<double,3> UCartMesh::getNodeCoordinates( int J ){

    return  getNodeCoordinates( getNodeCartesianId(J) );

};

/* -------------------------------------------------------------------------- */
/*! Get linear index of closest point neighbour in given direction.
 *  If neighbour does not exist boundary point is returned
 *  \param[in]  I       reference point index
 *  \param[in]  dir     direction [0/1/2/3/4/5] 
 *  \return neighbour index
 */
int      UCartMesh::getNodeNeighbour( int const &I, int const &dir){

    int     d = m_whichDirection[dir], step = m_whichStep[dir];

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

    std::array<int,3>    i;

    i       = getNodeCartesianId(I) ;
    i[d]    += step ; 
    i[d]    = std::max( std::min( i[d], m_nc[d] ), 0 ) ;

    return  getNodeLinearId(i) ;
};

/* -------------------------------------------------------------------------- */
/*! Get linear index of closest cell neighbour in given direction.
 *  If neighbour does not exist boundary cell is returned
 *  \param[in]  I       reference point index
 *  \param[in]  dir     direction [0/1/2/3/4/5] 
 *  \return neighbour index
 */
int      UCartMesh::getCellNeighbour( int const &I, int const &dir){

    int         d = m_whichDirection[dir], step = m_whichStep[dir];

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

    std::array<int,3>    i;

    i       = getCellCartesianId(I) ;
    i[d]    += step ; 
    i[d]    = std::max( std::min( i[d], m_nc[d]-1 ), 0 ) ;

    return  getCellLinearId(i) ;
};

/* -------------------------------------------------------------------------- */
/*! Compute the cartesian indices of closest node to given point.
 *  \param[in]  P       point coordinates
 *  \return             cartesian indices
 */
std::array<int,3> UCartMesh::getNodeCartesianId( std::array<double,3> const &P ){

    int         d ;
    std::array<int,3>    id;

    id.fill(0) ;

    for( d=0; d<m_dim; ++d){
        id[d] = std::min( m_np[d]-1, std::max(0, (int) round( (P[d] - m_B0[d])/m_h[d] )) );
    };

    return id; 
};

/* -------------------------------------------------------------------------- */
/*! Compute the cartesian indices of closest node to given point in a 2D Mesh
 *  \param[in]  P       point coordinates
 *  \param[out]  i      first cartesian index
 *  \param[out]  j      second cartesian index
 */
void UCartMesh::getNodeCartesianId( std::array<double,3> const &P, int &i, int &j ){


    i = std::min( m_np[0]-1, std::max(0, (int) round( (P[0] - m_B0[0])/m_h[0] )) );
    j = std::min( m_np[1]-1, std::max(0, (int) round( (P[1] - m_B0[1])/m_h[1] )) );

    return ;

};

/* -------------------------------------------------------------------------- */
/*! Compute the cartesian indices of closest node to given point in a 3D Mesh
 *  \param[in]  P       point coordinates
 *  \param[out]  i      first cartesian index
 *  \param[out]  j      second cartesian index
 *  \param[out]  k      third cartesian index
 */
void UCartMesh::getNodeCartesianId( std::array<double,3> const &P, int &i, int &j, int &k ){


    i = std::min( m_np[0]-1, std::max(0, (int) round( (P[0] - m_B0[0])/m_h[0] )) );
    j = std::min( m_np[1]-1, std::max(0, (int) round( (P[1] - m_B0[1])/m_h[1] )) );
    k = std::min( m_np[2]-1, std::max(0, (int) round( (P[2] - m_B0[2])/m_h[2] )) );

    return ;

};

/* -------------------------------------------------------------------------- */
/*! Transformation from linear index to cartesian indices.
 * No check on bounds is performed
 *  \param[in]  J       node linear index
 *  \return             node cartesian indices
 */
std::array<int,3> UCartMesh::getNodeCartesianId( int const &J ){

    // Local variables
    std::array<int,3>    id ;

    id[0] = J % m_np[0] ;
    id[2] = J / m_NodesInIJPlane;
    id[1] =  (J - id[2] *m_NodesInIJPlane ) /m_np[0]  ;

    return id; 

};

/* -------------------------------------------------------------------------- */
/*! Transformation from linear index to cartesian indices in a 3D mesh
 * No check on bounds is performed
 *  \param[in]  J       node linear index
 *  \param[out]  i      first cartesian  index
 *  \param[out]  j      second cartesian  index
 */
void UCartMesh::getNodeCartesianId( int const &J, int &i, int &j ){

    i = J %m_np[0] ;
    j = J /m_np[0]  ;

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
void UCartMesh::getNodeCartesianId( int const &J, int &i, int &j, int &k ){

    i = J % m_np[0] ;
    k = J / m_NodesInIJPlane;
    j =  (J - k *m_NodesInIJPlane ) /m_np[0]  ;

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
int UCartMesh::getNodeLinearId( int const &i, int const &j, int const &k ){

    return (m_NodesInIJPlane *k + m_np[0]*j + i); 

};

/* -------------------------------------------------------------------------- */
/*! Transformation from cartesian indices to linear indices.
 * No check on bounds is performed
 *  \param[in]  id      node cartesian indices
 *  \return             node linear index
 */
int UCartMesh::getNodeLinearId( std::array<double,3> const &P ){

    // Counters
    std::array<int,3>    id( getNodeCartesianId(P) );

    return getNodeLinearId(id); 

};

/* -------------------------------------------------------------------------- */
/*! Transformation from cartesian indices to linear indices.
 * No check on bounds is performed
 *  \param[in]  id      node cartesian indices
 *  \return             node linear index
 */
int UCartMesh::getNodeLinearId( std::array<int,3> const &id){

    return ( getNodeLinearId(id[0], id[1], id[2] ) ); 

};

/* -------------------------------------------------------------------------- */
/*! Calculate cell subset indices form cartesian indices
 *  \param[in]  i0      min cartesian indices
 *  \param[in]  i1      max cartesian indices
 *  \return             cell linear indices of subset mesh
 */
std::vector<int> UCartMesh::extractCellSubSet( std::array<int,3> const &i0, std::array<int,3> const &i1 ){

    int                     i, j, k; 
    std::vector<int>               ids;
    std::vector<int>::iterator     it;

    i  =  i1[0]-i0[0]+1  ;
    j  =  i1[1]-i0[1]+1  ;
    k  =  i1[2]-i0[2]+1  ;

    i  =  i *j *k ;
    ids.resize(i) ;

    it = ids.begin() ;

    for( k=i0[2]; k<=i1[2]; ++k){
        for( j=i0[1]; j<=i1[1]; ++j){
            for( i=i0[0]; i<=i1[0]; ++i){

                *it = getCellLinearId( i, j, k) ;            
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
std::vector<int> UCartMesh::extractCellSubSet( int const &I0, int const &I1 ){

    return extractCellSubSet( getCellCartesianId(I0), getCellCartesianId(I1) ); 

};

/* -------------------------------------------------------------------------- */
/*! Calculate cell subset indices form min and max point.
 *  The cell conataining the points (or closest to them) are maintained
 *  \param[in]  P0      min point
 *  \param[in]  P1      max point
 *  \return             cell linear indices of subset mesh
 */
std::vector<int> UCartMesh::extractCellSubSet( std::array<double,3> const &P0, std::array<double,3> const &P1 ){

    return extractCellSubSet( getCellCartesianId(P0), getCellCartesianId(P1) ); 

};

/* -------------------------------------------------------------------------- */
/*! Calculate subset indices form cartesian indices
 *  \param[in]  i0      min cartesian indices
 *  \param[in]  i1      max cartesian indices
 *  \return             node linear indices of subset mesh
 */
std::vector<int> UCartMesh::extractNodeSubSet( std::array<int,3> const &i0, std::array<int,3> const &i1 ){

    int                     i, j, k; 
    std::vector<int>               ids;
    std::vector<int>::iterator     it;

    i  =  i1[0]-i0[0]+1  ;
    j  =  i1[1]-i0[1]+1  ;
    k  =  i1[2]-i0[2]+1  ;

    i  =  i *j *k ;
    ids.resize(i) ;

    it = ids.begin() ;

    for( k=i0[2]; k<=i1[2]; ++k){
        for( j=i0[1]; j<=i1[1]; ++j){
            for( i=i0[0]; i<=i1[0]; ++i){

                *it = getNodeLinearId( i, j, k) ;            
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
std::vector<int> UCartMesh::extractNodeSubSet( int const &I0, int const &I1 ){

    return extractNodeSubSet( getNodeCartesianId(I0), getNodeCartesianId(I1) ); 

};

/* -------------------------------------------------------------------------- */
/*! Calculate node subset indices form min and max point.
 *  The nodes closest to the points are used as limites 
 *  \param[in]  P0      min point
 *  \param[in]  P1      max point
 *  \return             cell linear indices of subset mesh
 */
std::vector<int> UCartMesh::extractNodeSubSet( std::array<double,3> const &P0, std::array<double,3> const &P1 ){

    return extractNodeSubSet( getNodeCartesianId(P0), getNodeCartesianId(P1) ); 

};

/* -------------------------------------------------------------------------- */
/*! Check if point lies within the mesh
 *  \param[in]  P      
 *  \return             true if point lies within grid
 */
bool UCartMesh::isPointInGrid(  std::array<double,3> const &P ){

    int     d;

    for( d=0; d<m_dim; ++d){

        if( P[d]< m_B0[d] || P[d] > m_B1[d] ){
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
bool UCartMesh::isPointInGrid( std::array<double,3> const &P, std::array<int,3> &I){

    int     d;

    for( d=0; d<m_dim; ++d){

        if( P[d]< m_B0[d] || P[d] > m_B1[d] ){
            return false;
        };
    };

    I = getCellCartesianId(P) ;

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
bool UCartMesh::isPointInGrid( std::array<double,3> const &P, int &i, int &j, int &k ){

    bool        inGrid ;
    std::array<int,3>    I;

    if(isPointInGrid(P,I) ){

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
bool UCartMesh::isPointInGrid( std::array<double,3> const &P, int &I ){

    int     d;

    for( d=0; d<m_dim; ++d){

        if( P[d]< m_B0[d] || P[d] > m_B1[d] ){
            return false;
        };
    };

    I = getCellLinearId(P) ;

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
void UCartMesh::convertToUnstructured( int &nV, int &nS, std::vector<std::array<double,3>> &V, std::vector<std::vector<int>> &S, std::vector<std::vector<std::vector<int>>> &A ){

    // Local variables
    int         nv, ns;

    // Counters
    int         i, j, k, J;

    // Number of new vertices/simplicies
    nv = m_nNodes ;
    ns = m_nCells ;

    // Resize vertex list
    V.resize(nV + nv);

    // Resize simplex list
    S.resize(nS + ns, std::vector<int>(pow(2,m_dim), -1));

    // Resize adjacency
    A.resize(nS + ns, std::vector<std::vector<int>>(2*m_dim, std::vector<int>(1, -1)));

    // ========================================================================== //
    // GENERATE VERTEX LIST                                                       //
    // ========================================================================== //
    for (k = 0; k < m_np[2]; k++) {
        for (j = 0; j < m_np[1]; j++) {
            for (i = 0; i < m_np[0]; i++) {
                J = getNodeLinearId(i,j,k);
                V[J][0] = m_edge[0][i] ;
                V[J][1] = m_edge[1][j] ;
                V[J][2] = m_edge[2][k] ; 
                nV++;
            };
        }
    } 

    // ========================================================================== //
    // SIMPLEX-VERTEX CONNECTIVITY                                                //
    // ========================================================================== //
    for (k = 0; k < m_nc[2]; k++) {
        for (j = 0; j < m_nc[1]; j++) {
            for (i = 0; i < m_nc[0]; i++) {
                J = getCellLinearId(i,j,k);

                S[J][0] = getNodeLinearId(i,j,k);
                S[J][1] = getNodeLinearId(i+1,j,k);
                S[J][2] = getNodeLinearId(i,j+1,k);
                S[J][3] = getNodeLinearId(i+1,j+1,k);

                if(m_dim==3){
                    S[J][4] = getNodeLinearId(i,j,k+1);
                    S[J][5] = getNodeLinearId(i+1,j,k+1);
                    S[J][6] = getNodeLinearId(i,j+1,k+1);
                    S[J][7] = getNodeLinearId(i+1,j+1,k+1);
                }

            }
        }
    } 

    // ========================================================================== //
    // SIMPLEX-VERTEX ADJACENCY                                                   //
    // ========================================================================== //
    for (k = 0; k < m_nc[2]; k++) {
        for (j = 0; j < m_nc[1]; j++) {
            for (i = 0; i < m_nc[0]; i++) {
                J = getCellLinearId(i,j);

                if (i != 0)     { A[J][0][0] = getCellLinearId(i-1,j,k); }
                if (i != m_nc[0]) { A[J][1][0] = getCellLinearId(i+1,j,k); }

                if (j != 0)     { A[J][2][0] = getCellLinearId(i,j-1,k); }
                if (j != m_nc[1]) { A[J][3][0] = getCellLinearId(i,j+1,k); }

                if (k != 0)     { A[J][4][0] = getCellLinearId(i,j,k-1); }
                if (k != m_nc[2]) { A[J][5][0] = getCellLinearId(i,j,k+1); }

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
void UCartMesh::exportVTR( std::string dir, std::string filename) {

    // Local variables

    class VTKOut : public bitpit::VTKRectilinearGrid{
        private:
            std::vector<std::vector<double>>   *ptredge ;

        public:
            VTKOut(std::string dir, std::string filename, std::vector<std::vector<double>> &edges) { 
                setNames( dir, filename) ;
                setCodex( bitpit::VTKFormat::APPENDED ) ;
                ptredge = &edges ;
                setDimensions( 0, (*ptredge)[0].size()-1, 0, (*ptredge)[1].size()-1, 0, (*ptredge)[2].size()-1 ) ;

            } ;

            void flushData(  std::fstream &str, bitpit::VTKFormat codex_, std::string name  ){

                if( name == "x_Coord"){
                    bitpit::genericIO::flushBINARY( str, (*ptredge)[0] ) ;
                }

                else if( name == "y_Coord"){
                    bitpit::genericIO::flushBINARY( str, (*ptredge)[1] ) ;
                }

                else if( name == "z_Coord"){
                    bitpit::genericIO::flushBINARY( str, (*ptredge)[2] ) ;
                };

                return ;
            };

    };

    VTKOut      myVTK( dir, filename, m_edge) ;
    myVTK.write() ;

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
void UCartMesh::exportVTR(std::string dir, std::string filename, std::string dataname, bitpit::VTKLocation location, std::vector<double> &data ) {


    class VTKOut : public bitpit::VTKRectilinearGrid{
        private:
            std::vector<std::vector<double>>   *ptredge ;
            std::vector<double>   *ptrdata ;

        public:
            VTKOut(std::string dir, std::string filename, std::vector<std::vector<double>> &edges, std::string dataname, bitpit::VTKLocation location, std::vector<double> &data){ 
                setNames( dir, filename) ;
                setCodex( bitpit::VTKFormat::APPENDED ) ;
                ptredge = &edges ;
                setDimensions( 0, (*ptredge)[0].size()-1, 0, (*ptredge)[1].size()-1, 0, (*ptredge)[2].size()-1 ) ;

                ptrdata = &data ;
                addData( dataname, bitpit::VTKFieldType::SCALAR, location, bitpit::VTKDataType::Float64 ) ;

            } ;

            void flushData(  std::fstream &str, bitpit::VTKFormat codex_, std::string name  ){

                if( name == "x_Coord"){
                    bitpit::genericIO::flushBINARY( str, (*ptredge)[0] ) ;
                }

                else if( name == "y_Coord"){
                    bitpit::genericIO::flushBINARY( str, (*ptredge)[1] ) ;
                }

                else if( name == "z_Coord"){
                    bitpit::genericIO::flushBINARY( str, (*ptredge)[2] ) ;
                }

                else{
                    bitpit::genericIO::flushBINARY( str, (*ptrdata) ) ;
                };

                return ;
            };

    };

    VTKOut      myVTK( dir, filename, m_edge, dataname, location, data) ;
    myVTK.write() ;

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
void UCartMesh::exportVTR(std::string dir, std::string filename, std::string dataname, bitpit::VTKLocation location, std::vector<int> &data ) {


    class VTKOut : public bitpit::VTKRectilinearGrid{
        private:
            std::vector<std::vector<double>>   *ptredge ;
            std::vector<int>   *ptrdata ;

        public:
            VTKOut(std::string dir, std::string filename, std::vector<std::vector<double>> &edges, std::string dataname, bitpit::VTKLocation location, std::vector<int> &data){ 
                setNames( dir, filename) ;
                setCodex( bitpit::VTKFormat::APPENDED ) ;
                ptredge = &edges ;
                setDimensions( 0, (*ptredge)[0].size()-1, 0, (*ptredge)[1].size()-1, 0, (*ptredge)[2].size()-1 ) ;

                ptrdata = &data ;
                addData( dataname, bitpit::VTKFieldType::SCALAR, location, bitpit::VTKDataType::Int32 ) ;

            } ;

            void flushData(  std::fstream &str, bitpit::VTKFormat codex_, std::string name  ){

                if( name == "x_Coord"){
                    bitpit::genericIO::flushBINARY( str, (*ptredge)[0] ) ;
                }

                else if( name == "y_Coord"){
                    bitpit::genericIO::flushBINARY( str, (*ptredge)[1] ) ;
                }

                else if( name == "z_Coord"){
                    bitpit::genericIO::flushBINARY( str, (*ptredge)[2] ) ;
                }

                else{
                    bitpit::genericIO::flushBINARY( str, (*ptrdata) ) ;
                };

                return ;
            };

    };

    VTKOut      myVTK( dir, filename, m_edge, dataname, location, data) ;
    myVTK.write() ;

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
void UCartMesh::exportVTR(std::string dir, std::string filename, std::string dataname, bitpit::VTKLocation location, std::vector<std::array<double,3>> &data ) {


    class VTKOut : public bitpit::VTKRectilinearGrid{
        private:
            std::vector<std::vector<double>>    *ptredge ;
            std::vector<std::array<double,3>>   *ptrdata ;

        public:
            VTKOut(std::string dir, std::string filename, std::vector<std::vector<double>> &edges, std::string dataname, bitpit::VTKLocation location, std::vector<std::array<double,3>> &data){ 
                setNames( dir, filename) ;
                setCodex( bitpit::VTKFormat::APPENDED ) ;
                ptredge = &edges ;
                setDimensions( 0, (*ptredge)[0].size()-1, 0, (*ptredge)[1].size()-1, 0, (*ptredge)[2].size()-1 ) ;

                ptrdata = &data ;
                addData( dataname, bitpit::VTKFieldType::VECTOR, location, bitpit::VTKDataType::Float64 ) ;

            } ;

            void flushData(  std::fstream &str, bitpit::VTKFormat codex_, std::string name ){

                if( name == "x_Coord"){
                    bitpit::genericIO::flushBINARY( str, (*ptredge)[0] ) ;
                }

                else if( name == "y_Coord"){
                    bitpit::genericIO::flushBINARY( str, (*ptredge)[1] ) ;
                }

                else if( name == "z_Coord"){
                    bitpit::genericIO::flushBINARY( str, (*ptredge)[2] ) ;
                }

                else{
                    bitpit::genericIO::flushBINARY( str, (*ptrdata) ) ;
                };

                return ;
            };

    };

    VTKOut      myVTK( dir, filename, m_edge, dataname, location, data) ;
    myVTK.write() ;

    return; 

};

/* -------------------------------------------------------------------------- */
/*!  Transform cell data to point data by calculating the mean of incident cells in each vertex
 *  \param[in]     CellData  Data on cells
 *  \param[out]    NodeData  Data on nodes
 */
void UCartMesh::convertCellDataToNodeData( std::vector<double> &CellData, std::vector<double> &NodeData ){
    // ========================================================================== //

    // Local variables
    int            K, J;
    int            ip, jp, kp;
    std::vector<int>     NodeIter ;

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
    for (k = 0; k < m_nc[2]; k++) {
        for (j = 0; j < m_nc[1]; j++) {
            for (i = 0; i < m_nc[0]; i++) {

                // Cell index
                K = getCellLinearId(i, j, k);

                for (n = 0; n < m_dim-1; n++) {
                    for (m = 0; m < 2; m++) {
                        for (l = 0; l < 2; l++) {

                            // Point index
                            ip = i + l;
                            jp = j + m;
                            kp = j + n;
                            J = getNodeLinearId(ip,jp,kp);

                            NodeData[J] = NodeData[J] + CellData[K]; 
                            NodeIter[J]++ ;
                        } //next n
                    } //next m
                } //next l

            } //next k
        } //next j
    } //next i


    for( J=0; J<m_nNodes; ++J){
        NodeData[J] = NodeData[J] / ((float) NodeIter[J]) ;
    };

    return; 

};

/* -------------------------------------------------------------------------- */
/*!  Transform node data to cell data by calculating the mean of incident vertices of each cell
 *  \param[in]    PointData  Data on nodes
 *  \param[out]   CellData  Data on cells
 */
void UCartMesh::convertNodeDataToCellData( std::vector<double> &NodeData, std::vector<double> &CellData ){


    // Local variables
    int     K, J;
    int     ip, jp, kp;
    double  factor ;

    // Counters
    int     i, j, k, l, m, n;

    factor  =   pow(0.5,m_dim) ;

    // ========================================================================== //
    // RESIZE OUTPUT VARIABLES                                                    //
    // ========================================================================== //
    CellData.resize(getNCells(), 0.0);

    // ========================================================================== //
    // CONVERT POINT DATA TO CELL DATA                                            //
    // ========================================================================== //
    for (k = 0; k < m_nc[2]; k++) {
        for (j = 0; j < m_nc[1]; j++) {
            for (i = 0; i < m_nc[0]; i++) {

                K = getCellLinearId(i,j,k);
                for (n = 0; n < 2; n++) {
                    for (m = 0; m < 2; m++) {
                        for (l = 0; l < 2; l++) {
                            ip = i + l;
                            jp = j + m;
                            kp = k + n;
                            J = getNodeLinearId(ip,jp,kp);
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
int UCartMesh::linearCellInterpolation( std::array<double,3> &P, std::vector<int> &stencil, std::vector<double> &weights ){

    int                         nStencil(0) ;
    int                         d, i, j, k;
    std::array<int,3>                    i0, i1 ;

    std::array< std::array<int,2>, 3>     cStencil ;
    std::array< std::array<double,2>, 3>  cWeights ;

    std::array<int,3>                nS({1,1,1}) ;

    std::vector<int>::iterator       itrStencil ;
    std::vector<double>::iterator    itrWeights ;

    if( isPointInGrid(P, i0) ){

        nStencil = 1 ;

        for( d=0; d<m_dim; ++d){

            // Find cell index
            if( P[d] < m_center[d][i0[d] ] ){
                i0[d] = i0[d]-1 ;
            };

            i1[d] = i0[d] +1 ;

            if( i0[d] < 0 ){
                cStencil[d][0] = 0 ; 
                cWeights[d][0] = 1. ; 

            } 

            else if( i1[d] > m_nc[d]-1 ){
                cStencil[d][0] = m_nc[d]-1 ; 
                cWeights[d][0] = 1. ; 

            } 

            else{
                nS[d] = 2 ;

                cStencil[d][0] = i0[d] ;
                cStencil[d][1] = i1[d] ;

                cWeights[d][1] = (P[d] - m_center[d][i0[d]]) /m_h[d]   ;
                cWeights[d][0] = 1.0 - cWeights[d][1] ;  
            }


        };

        for( d=m_dim; d<3; ++d){
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

                    *itrStencil =  getCellLinearId(is,js,ks) ;
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
int UCartMesh::linearNodeInterpolation( std::array<double,3> &P, std::vector<int> &stencil, std::vector<double> &weights ){


    int                         nStencil(0) ;
    int                         d, i, j, k;
    std::array<int,3>                    i0, i1 ;

    std::array< std::array<int,2>, 3>     cStencil ;
    std::array< std::array<double,2>, 3>  cWeights ;

    std::vector<int>::iterator       itrStencil ;
    std::vector<double>::iterator    itrWeights ;

    if( isPointInGrid(P,i0) ){

        nStencil = pow(2,m_dim) ;

        for(d=0; d<m_dim; ++d){
            i1[d] = i0[d] +1 ;

            cStencil[d][0] = i0[d] ;
            cStencil[d][1] = i1[d] ;

            cWeights[d][1] = ( P[d] - m_edge[d][i0[d]]) /m_h[d]  ;
            cWeights[d][0] = 1.0 - cWeights[d][1] ;  
        };

        for( d=m_dim; d<3; ++d){
            cStencil[d][0] = 0 ;
            cWeights[d][0] = 1. ;
        };

        stencil.resize(nStencil) ;
        weights.resize(nStencil) ;

        itrStencil = stencil.begin() ;
        itrWeights = weights.begin() ;

        for( k=0; k<m_dim-1; ++k){
            for( j=0; j<2; ++j){
                for( i=0; i<2; ++i){

                    int &is = cStencil[0][i] ;
                    int &js = cStencil[1][j] ;
                    int &ks = cStencil[2][k] ;

                    double &iw = cWeights[0][i] ;
                    double &jw = cWeights[1][j] ;
                    double &kw = cWeights[2][k] ;

                    *itrStencil =  getNodeLinearId(is,js,ks) ;
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
