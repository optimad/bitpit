/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

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

# include "bitpit_SA.hpp"
# include "bitpit_operators.hpp"

# include "levelSet.hpp"

namespace bitpit {

/*!
	@ingroup levelset
	@class  LevelSet
	@brief  Level Set on Stl Manager Class

	LevelSet is a user interface class. One user should (read can...) work only
	with this class and its methods to maintain a signed distance function (Sdf)
	computed from a piece-wise linear approximation of a d manifold in a 3D Euclidean
	space. Sdf is computed in a narrow band of at least 2 mesh cell centers
	around the geometry.
*/

/*!
 * Default constructor
 */
LevelSet::LevelSet() {

    m_mesh = NULL ;

    RSearch     = 0;
    m_userRSearch = false ;

    signedDF    = true ;
    propagateS  = false;
    propagateV  = false;
};

/*!
 * Constructor
 * @param[in] patch underlying mesh
 */
LevelSet::LevelSet( VolumeKernel *patch): LevelSet() {
    m_mesh = patch ;
};

/*!
 * Destructor of LevelSet_Stl.
*/
LevelSet::~LevelSet(){
    m_mesh = NULL ;
};

/*!
 * Get the Sdf value of the i-th local element of the octree mesh.
 * @param[in] i Local index of target octant.
 * @return Value of the i-th local element of the octree mesh.
 */
double LevelSet::getLS( const long &i)const {

    if( !info.exists(i) ){
        return levelSetDefaults::VALUE;
    } else {
        return ( info[i].value );
    };

};

/*!
 * Get the Sdf gradient vector of the i-th local element of the octree mesh.
 * @param[in] i Local index of target octant.
 * @return Array with components of the Sdf gradient of the i-th local element of the octree mesh.
 */
std::array<double,3> LevelSet::getGradient(const long &i) const {

    if( !info.exists(i) ){
        return levelSetDefaults::GRADIENT;
    } else {
        return ( info[i].gradient );
    };

};

/*!
 * Get if the Sdf value of the i-th local element is exactly computed or not.
 * @param[in] i Local index of target octant.
 * @return True/false if the Sdf value is exactly computed (true) or not (false).
 */
bool LevelSet::isInNarrowBand(const long &i)const{

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
double LevelSet::getSizeNarrowBand()const{
    return RSearch;
};

/*!
 * Set if the signed or unsigned LevelSet should be computed.
 * @param[in] flag true/false for signed /unsigned Level-Set function .
 */
void LevelSet::setSign(bool flag){
    signedDF = flag;

};

/*!
 * Set if the levelset sign has to be propagated from the narrow band to the whole domain.
 * @param[in] flag True/false to active/disable the propagation .
 */
void LevelSet::setPropagateSign(bool flag){
    propagateS = flag;
};

/*!
 * Set if the levelset value has to be propagated from the narrow band to the whole domain.
 * @param[in] flag True/false to active/disable the propagation.
 */
void LevelSet::setPropagateValue(bool flag){
    propagateV = flag;
};

/*!
 * Manually set the physical size of the narrow band.
 * @param[in] r Size of the narrow band.
 */
void LevelSet::setSizeNarrowBand(double r){
    RSearch = r;
    m_userRSearch = true ;
};

/*!
 * Computes the Level Set Gradient on on cell by first order finite-volume upwind stencil
 * @param[in] I index of cell 
 * @return gradient
 */
std::array<double,3> LevelSet::computeGradientUpwind( const long &I ){

    int                     i, border ;
    std::array<double,3>    normal, gradient = {{0.,0.,0.}};
    long                    owner;

    Cell            &cell = m_mesh->getCell(I) ;
    int                     nI = cell.getInterfaceCount() ;
    const long*             interfaces = cell.getInterfaces() ;
    const long*             neighbours = cell.getAdjacencies() ;

    long                    F, N ;
    double                  value, area;


    for( i=0; i<nI; ++i){
        F = interfaces[i] ;
        N = neighbours[i] ;

        Interface       &interface = m_mesh->getInterface(F) ;

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
std::array<double,3> LevelSet::computeGradientCentral( const long &I ){

    int                     i, border ;
    std::array<double,3>    normal, gradient = {{0.,0.,0.}};
    long                    owner;

    Cell            &cell = m_mesh->getCell(I) ;
    int                     nI = cell.getInterfaceCount() ;
    const long*             interfaces = cell.getInterfaces() ;
    const long*             neighbours = cell.getAdjacencies() ;

    long                    F, N;
    double                  value, area;

    for( i=0; i<nI; ++i){
        F = interfaces[i] ;
        N = neighbours[i] ;

        Interface       &interface = m_mesh->getInterface(F) ;

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
//TODO std::blablavector<bool> LevelSet::propagateSign( LSObject &visitor ) {
//TODO 
//TODO 
//TODO     // Local variables
//TODO     long                        seed;
//TODO     double                      s;
//TODO 
//TODO     std::vector<bool>          	flag(m_mesh->getCellCount(), true);
//TODO     LIFOStack<int>     	stack(sqrt(m_mesh->getCellCount()));
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
void LevelSet::propagateValue( LSObject *visitor ){

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
void LevelSet::solveEikonal( double g, double s ){

    long                            N( m_mesh->getCellCount()) ;
    PiercedVector<short>    active ;

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

                    check = (s*info[*it].value >= 0.0) && (abs(info[*it].value) < levelSetDefaults::VALUE) ;
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
        MinPQueue<double, long> heap(m, true, mapPtr);

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
double LevelSet::updateEikonal( double s, double g, const long &I ){

    BITPIT_UNUSED(s) ;
    BITPIT_UNUSED(g) ;
    BITPIT_UNUSED(I) ;

    return levelSetDefaults::VALUE;
};

/*! 
 * Update scalar field value at mesh vertex on a by locally solving the 3D Eikonal equation.
 * @param[in] s Flag for inwards/outwards propagation (s = -+1).
 * @param[in] g Propagation speed for the 3D Eikonal equation.
 * @param[in] I index of the cartesian cell to be updated.
 * @return Updated value at mesh vertex
 */
void LevelSet::clearAfterRefinement( std::vector<Adaption::Info> &mapper ){

    long id ;
    for ( auto & map : mapper ){
        if( map.entity == Adaption::Entity::ENTITY_CELL ){

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

}
