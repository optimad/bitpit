/*!
 *	\date			10/jul/2014
 *	\authors		Alessandro Alaia
 *	\authors		Haysam Telib
 *	\authors		Edoardo Lombardi
 *	\version		0.1
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of Class_LevelSet_Stl is released under the LGPL License.
 *
 *	\brief Level Set Manager Class - 3D PABLO Octree specialization
 *
 *	Level Set Stl is a user interface class. One user should (read can...) work only
 *	with this Class and its methods to maintain a signed distance function (Sdf)
 *	computed from a piece-wise linear approximation of a d manifold in a 3D Euclidean
 *	space. Sdf is computed in a narrow band of at least 2 mesh cell centers
 *	around the geometry.
 *	Parallel implementation developed by using the features of PABLO library.
 *
 */

# include "SortAlgorithms.hpp"
# include "Operators.hpp"

# include "LevelSet.hpp"

/*!
 * Default constructor
 */
bitpit::LevelSet::LevelSet() {

    m_mesh = NULL ;

    RSearch     = 0;

    signedDF    = true ;
    propagateS  = false;
    propagateV  = false;
};

/*!
 * Constructor
 * @param[in] patch underlying mesh
 */
bitpit::LevelSet::LevelSet( bitpit::VolumeKernel *patch): LevelSet() {
    m_mesh = patch ;
};

/*!
 * Destructor of LevelSet_Stl.
*/
bitpit::LevelSet::~LevelSet(){
    m_mesh = NULL ;
};

/*!
 * Get the Sdf value of the i-th local element of the octree mesh.
 * @param[in] i Local index of target octant.
 * @return Value of the i-th local element of the octree mesh.
 */
double bitpit::LevelSet::getLS( const long &i){

    if( !info.exists(i) ){
        return 1.e18;
    } else {
        return ( info[i].value );
    };

};

/*!
 * Get the Sdf gradient vector of the i-th local element of the octree mesh.
 * @param[in] i Local index of target octant.
 * @return Array with components of the Sdf gradient of the i-th local element of the octree mesh.
 */
std::array<double,3> bitpit::LevelSet::getGradient(const long &i){

    if( !info.exists(i) ){
        return {{0.,0.,0.}};
    } else {
        return ( info[i].gradient );
    };

};

/*!
 * Get if the Sdf value of the i-th local element is exactly computed or not.
 * @param[in] i Local index of target octant.
 * @return True/false if the Sdf value is exactly computed (true) or not (false).
 */
bool bitpit::LevelSet::isInNarrowBand(const long &i){

    if( !info.exists(i) ){
        return false;
    } else {
        return ( std::abs(info[i].value) <= RSearch );
    };

};

/*!
 * Get the current size of the narrow band.
 * @return Physical size of the current narrow band to guarantee at least one element inside it.
 */
double bitpit::LevelSet::getSizeNarrowBand(){
    return RSearch;
};

/*!
 * Set if the signed or unsigned LevelSet should be computed.
 * @param[in] flag true/false for signed /unsigned Level-Set function .
 */
void bitpit::LevelSet::setSign(bool flag){
    signedDF = flag;

};

/*!
 * Set if the levelset sign has to be propagated from the narrow band to the whole domain.
 * @param[in] flag True/false to active/disable the propagation .
 */
void bitpit::LevelSet::setPropagateSign(bool flag){
    propagateS = flag;
};

/*!
 * Set if the levelset value has to be propagated from the narrow band to the whole domain.
 * @param[in] flag True/false to active/disable the propagation.
 */
void bitpit::LevelSet::setPropagateValue(bool flag){
    propagateV = flag;
};

/*!
 * Manually set the physical size of the narrow band.
 * @param[in] r Size of the narrow band.
 */
void bitpit::LevelSet::setSizeNarrowBand(double r){
    RSearch = r;
};

/*!
 * Computes the Level Set Gradient on on cell by first order finite-volume upwind stencil
 * @param[in] I index of cell 
 * @return gradient
 */
std::array<double,3> bitpit::LevelSet::computeGradientUpwind( const long &I ){

    int                     i, border ;
    std::array<double,3>    normal, gradient = {{0.,0.,0.}};
    long                    owner;

    bitpit::Cell            &cell = m_mesh->getCell(I) ;
    int                     nI = cell.getInterfaceCount() ;
    const long*             interfaces = cell.getInterfaces() ;
    const long*             neighbours = cell.getAdjacencies() ;

    long                    F, N ;
    double                  value, area;


    for( i=0; i<nI; ++i){
        F = interfaces[i] ;
        N = neighbours[i] ;

        bitpit::Interface       &interface = m_mesh->getInterface(F) ;

        owner  = interface.getOwner() ;
        border = interface.isBorder() ;

        area   = m_mesh->evalInterfaceArea(F) ;
        normal = m_mesh->evalInterfaceNormal(F) ;

        if( owner != I)
            normal = -1. *normal ;

        value = border* info[I].value + (1-border) *std::min( info[I].value , info[N].value ) ;

        gradient += area *value *normal ;

    };

    gradient = gradient /m_mesh->evalCellVolume(I) ;


    return gradient;

};

/*!
 * Computes the Level Set Gradient on on cell by second order finite-volume central stencil
 * @param[in] I index of cell 
 * @return gradient
 */
std::array<double,3> bitpit::LevelSet::computeGradientCentral( const long &I ){

    int                     i, border ;
    std::array<double,3>    normal, gradient = {{0.,0.,0.}};
    long                    owner;

    bitpit::Cell            &cell = m_mesh->getCell(I) ;
    int                     nI = cell.getInterfaceCount() ;
    const long*             interfaces = cell.getInterfaces() ;
    const long*             neighbours = cell.getAdjacencies() ;

    long                    F, N;
    double                  value, area;

    for( i=0; i<nI; ++i){
        F = interfaces[i] ;
        N = neighbours[i] ;

        bitpit::Interface       &interface = m_mesh->getInterface(F) ;

        owner  = interface.getOwner() ;
        border = interface.isBorder() ;

        area   = m_mesh->evalInterfaceArea(F) ;
        normal = m_mesh->evalInterfaceNormal(F) ;

        if( owner != I)
            normal = -1. *normal ;

        value = border* info[I].value + (1-border) *(info[I].value +info[N].value)/2.  ;

        gradient += area *value *normal ;

    };

    gradient = gradient /m_mesh->evalCellVolume(I) ;


    return gradient;

};

/*!
 * Propagate the sign of the signed distance function from narrow band to entire domain
 */
//TODO std::blablavector<bool> bitpit::LevelSet::propagateSign( bitpit::LSObject &visitor ) {
//TODO 
//TODO 
//TODO     // Local variables
//TODO     long                        seed;
//TODO     double                      s;
//TODO 
//TODO     std::vector<bool>          	flag(m_mesh->getCellCount(), true);
//TODO     bitpit::LIFOStack<int>     	stack(sqrt(m_mesh->getCellCount()));
//TODO 
//TODO     std::vector<long>			neighs;
//TODO     std::vector<long>::iterator it, itend ;
//TODO 
//TODO     visitor.seedSign( *this, seed, s) ;
//TODO 
//TODO     stack.push(seed);
//TODO 
//TODO     // PROPAGATE SIGN                                                               
//TODO     while (stack.TOPSTK > 0) {
//TODO 
//TODO         // Pop item from stack
//TODO         seed = stack.pop();
//TODO 
//TODO         // Retrieve info
//TODO         flag[seed] = false;
//TODO 
//TODO         // Loop over neighbors
//TODO         neighs  =   m_mesh->findCellFaceNeighs( seed ) ;
//TODO 
//TODO         it    = neighs.begin() ;
//TODO         itend = neighs.end() ;
//TODO 
//TODO         for ( it != itend; ++it) {
//TODO 
//TODO             LSInfo &lsInfo = info[*it] ;
//TODO 
//TODO             if ( flag[*it] ) {
//TODO 
//TODO                 LSInfo &lsInfo = info[*it] ;
//TODO 
//TODO                 if (lsInfo.value == 1.0e+18) {
//TODO                     lsInfo.value = s * lsInfo.value;
//TODO                 }
//TODO 
//TODO                 stack.push(*it);
//TODO 
//TODO             } //endif
//TODO 
//TODO         } //iterator
//TODO 
//TODO     } //stack
//TODO 
//TODO     return;
//TODO 
//TODO };

/*!
 * Driver routine for propagtaing levelset value by a fast marching method
 */
void bitpit::LevelSet::propagateValue( bitpit::LSObject *visitor ){

    BITPIT_UNUSED(visitor) ;

    // Propagate outwards ------------------------------------------------------- //
    solveEikonal(1.0, 1.0);

    // Propagate inwards -------------------------------------------------------- //
    solveEikonal(1.0, -1.0);


    return ;

};

/*! 
 * Solve the 3D Eikonal equation |grad(u)| = g, using  a fast marching method
 * (Function in the unknown region must be set to the value 1.0e+18)
 * @param[in] g Propagation speed.
 * @param[in] s Velocity sign (+1 --> propagate outwards, -1 --> propagate inwards).
 */
void bitpit::LevelSet::solveEikonal( double g, double s ){

    long                            N( m_mesh->getCellCount()) ;
    bitpit::PiercedVector<short>    active ;

    { //FLAG ALIVE, DEAD AND FAR AWAY VERTEXES
        long    m(0), myId ;
        bool    check;

        std::vector<long>               neighs ;
        std::vector<long>::iterator     it, itbeg, itend ;

        active.reserve(N) ;
        for ( const auto &cell : m_mesh->getCells() ){ 
            myId           = cell.getId() ;
            active.reclaim(myId) ;
        }

        for ( const auto &cell : m_mesh->getCells() ){ 

            myId           = cell.getId() ;

            if ( isInNarrowBand(myId) ) // dead vertex
                active[myId] = 0;

            else { // alive or far away 

                // Loop over neighbors
                check   = false;
                neighs = m_mesh->findCellFaceNeighs(myId) ;

                it    = neighs.begin() ;
                itend = neighs.end() ;

                while( !check &&  it != itend  ){

                    check = (s*info[*it].value >= 0.0) && (abs(info[*it].value) < 1.0e+18) ;
                    ++it ;

                };

                active[myId] = 2 - (int) check  ;
                m += (int) check  ;

            }
        } //next i
    }


    { // Construct min heap data structure 
        long                            m(0), I(0), myId, J ;
        double                          value ;

        std::vector<long>               neighs ;
        std::vector<long>::iterator     it, itbeg, itend ;

        std::unordered_map<long,long>   contiguos ;
        std::vector<std::array<int,2>>  map(N), *mapPtr = &map;
        bitpit::MinPQueue<double, long> heap(m, true, mapPtr);

        for ( const auto &cell : m_mesh->getCells() ){ 

            myId = cell.getId() ;
            contiguos.insert( {{myId,I}} ) ;

            if (active[myId] == 1) {

                // Solve the quadratic form
                value = updateEikonal(s, g, myId);

                // Store value into heap
                map[m][0] = I;
                map[I][1] = m; 


                heap.keys[m] = value;
                heap.labels[m] = myId;

                // Update counter
                m++;
            }

            ++I;

        } //next i

        // Build min-heap 
        heap.heap_size = m;
        heap.buildHeap();

        // FAST MARCHING
        while (heap.heap_size > 0) {


            // Extract root
            heap.extract(value, myId);

            LSInfo &lsInfo = info[myId] ;

            // Update level set value
            lsInfo.value = s*updateEikonal(s, g, myId);

            // Update m_activeNode to dead node
            active[myId] = 0;

            // Upadate far-away neighboors
            neighs = m_mesh->findCellFaceNeighs(myId) ;

            itbeg = neighs.begin() ;
            itend = neighs.end() ;

            for( it=itbeg; it != itend; ++it){

                J = *it;

                // Update the local value
                value = updateEikonal(s, g, J);

                if (active[J] == 1) { // Update value in min-heap
                    I = contiguos[J] ;
                    heap.modify( map[I][1], value, J );
                }

                else if (active[J] == 2) {

                    // Update m_activeNode for neighbor
                    active[J] = 1;

                    I = contiguos[J] ;

                    // Insert neighbor into the min heap
                    map[heap.heap_size][0] = I ;
                    map[I][1] = heap.heap_size;

                    heap.insert(value, J);

                };
            };

        } //next item
    }


    return; 
};

/*! 
 * Update scalar field value at mesh vertex on a by locally solving the 3D Eikonal equation.
 * @param[in] s Flag for inwards/outwards propagation (s = -+1).
 * @param[in] g Propagation speed for the 3D Eikonal equation.
 * @param[in] I index of the cartesian cell to be updated.
 * @return Updated value at mesh vertex
 */
double bitpit::LevelSet::updateEikonal( double s, double g, const long &I ){

    BITPIT_UNUSED(s) ;
    BITPIT_UNUSED(g) ;
    BITPIT_UNUSED(I) ;

    return 1.e18;
};

/*! 
 * Update scalar field value at mesh vertex on a by locally solving the 3D Eikonal equation.
 * @param[in] s Flag for inwards/outwards propagation (s = -+1).
 * @param[in] g Propagation speed for the 3D Eikonal equation.
 * @param[in] I index of the cartesian cell to be updated.
 * @return Updated value at mesh vertex
 */
void bitpit::LevelSet::clearAfterRefinement( std::vector<bitpit::Adaption::Info> &mapper ){

    long id ;
    for ( auto & map : mapper ){
        if( map.entity == bitpit::Adaption::Entity::ENTITY_CELL ){

            for ( auto & parent : map.previous){ //save old data and delete element
                id = (long) parent ;
                if( info.exists(id) )
                    info.erase(id,true) ;
            }
        }
    }

    info.flush() ;

    return ;
};
