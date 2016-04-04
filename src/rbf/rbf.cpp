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

#include <lapacke.h>

#include <cmath>
#include "Operators.hpp"
#include "SortAlgorithms.hpp"
#include "rbf.hpp"

namespace bitpit{

/*!
 * @ingroup RBF
 * @{
 *
 * @class RBF
 * @brief Handling of Radial Basis Function with a large set of nodes
 *
 */

/*!
 * Destructor
 */
RBF::~RBF(){
    m_fPtr = NULL ;
} ;

/*! 
 * Default constructor
 */
RBF::RBF( RBFBasisFunction bfunc ) {
    m_supportRadius = 1. ;
    m_nodes         = 0 ;
    m_fields        = 0 ;

    m_node.clear() ;
    m_value.clear() ;
    m_weight.clear() ;
    m_active.clear() ;

    setFunction( bfunc ) ;
};

/*! 
 * Sets the rbf function to be used
 * @param[in] bfunc basis function to be used
 */
void RBF::setFunction( const RBFBasisFunction &bfunc ){

    switch(bfunc){

        case( RBFBasisFunction::WENDLANDC2):
            setFunction( rbf::wendlandc2);
            break;

        default:
            setFunction( rbf::wendlandc2);
            break;
    }

    return;
};

/*! 
 * Sets the rbf function to a user specified function
 * @param[in] bfunc basis function to be used
 */
void RBF::setFunction( double (&bfunc)(const double &) ){
    m_fPtr = bfunc ;
    return ;
};

/*! 
 * Set the support radius
 * @param[in] radius support radius
 */
void RBF::setSupportRadius( const double & radius ){
    m_supportRadius = radius ;
    return ;
};

/*! 
 * Get the number of fields to be interpolated
 * @return  number of fields
 */
int RBF::getFieldCount(  ){
    return m_fields ;
};

/*! 
 * Get the number of active nodes
 * @return  number of active nodes
 */
int RBF::getActiveCount(  ){

    int nActive(0);

    for( auto && active : m_active)
        nActive += (int) active ;

    return nActive ;
};

/*! 
 * Get the indices of the active nodes
 * @return  indices of active nodes
 */
std::vector<int> RBF::getActiveSet(  ){

    int                 i(0);
    std::vector<int>    activeSet ;
   
    activeSet.reserve( getActiveCount() ) ; 

    for( auto && active : m_active){
        if( active )
            activeSet.push_back(i) ;
        ++i ;
    }

    return activeSet ;
};

/*! 
 * Checks if a node is active
 * @param[in] n index of node to be checked
 * @return true if active
 */
bool RBF::isActive( const int &n ){
    return m_active[n] ;
};


/*! 
 * Adds a RBF node and sets it to active
 * @param[in] node  coordinates of node to be added
 * @return id of node within class
 */
int RBF::addNode( const std::array<double,3> &node ){
    m_node.push_back(node) ;
    m_active.push_back(true) ;

    m_nodes++ ;
    return m_nodes-1 ;
};

/*! 
 * Adds a list of RBF nodes and sets them to active
 * @param[in] node  coordinates of nodes to be added
 * @return id of node within class
 */
std::vector<int> RBF::addNode( const std::vector<std::array<double,3>> &node ){

    int                 i( m_nodes ) ;
    std::vector<int>    ids;

    ids.resize( node.size() ) ;

    for( auto & id:ids ){
        id = i ;
        ++i;
    }

    m_node.insert( m_node.end(), node.begin(), node.end() ) ;
    m_nodes += node.size() ;
    
    m_active.resize( m_nodes, true );

    return ids;
};

/*! 
 * Adds a field to be interpolated
 * @return id of field within th class
 */
int RBF::addField( ){
    m_fields++ ;
    m_value.resize(m_fields) ;
    return m_fields-1 ;
};

/*! 
 * Adds a field to be interpolated
 * @param[in] field values of field
 * @return id of field within th class
 */
int RBF::addField( const std::vector<double> & field ){
    m_value.push_back(field) ;
    m_fields++ ;
    return m_fields-1 ;
};

/*! 
 * Sets all the field values at one node
 * @param[in] id id of node
 * @param[in] value  function values to be interpolated
 */
void RBF::setNodeValue( const int &id, const std::vector<double> &value ){

    int i ;

    for( i=0; i<m_fields; ++i ){
        m_value[i][id] = value[i] ;
    }
    return ;
};

/*! 
 * Sets the values of one field at all nodes
 * @param[in] id id of field
 * @param[in] value  function values to be interpolated
 */
void RBF::setFieldValue( const int &id, const std::vector<double> &value ){

    m_value[id] = value ;

    return ;
};

/*! 
 * Evaluates the basis function
 * @param[in] dist distance
 * @return value of basis function
 */
double RBF::evalBasis( const double &dist ){
    return (*m_fPtr)(dist) ;
};

/*! 
 * Evaluates the RBF
 * @param[in] point point where to evaluate the basis
 * @return vector containing interpolated values
 */
std::vector<double> RBF::evalRBF( const std::array<double,3> &point){

    std::vector<double> values(m_fields, 0.) ;
    int                 i, j ;
    double              dist, basis ;
    
    for( i=0; i<m_nodes; ++i ){

        if( m_active[i] ) {
            dist = norm2( point - m_node[i] ) / m_supportRadius ;
            basis = evalBasis( dist ) ;

            for( j=0; j<m_fields; ++j){
                values[j] += basis * m_weight[j][i] ;       
            }
        }

    }

    return values ;
}

/*! 
 * Calculates the RBF weights using active nodes
 */
void RBF::solve(){

    int i, j, k ;
    double dist;

    int nS      = getActiveCount() ;
    int nrhs    = getFieldCount() ;

    int lda     = nS;
    int ldb     = nS;
    int info ;
    int ipiv[nS];

    std::vector<int> activeSet( getActiveSet() ) ;

    double *a = new double [lda * nS];
    double *b = new double [ldb * nrhs];


    k=0 ;
    for( j=0; j<nrhs; ++j){

        for( const auto & i : activeSet ){
            b[k] = m_value[j][i];
            ++k;
        }

    }

    k=0;
    for( const auto &i : activeSet ){
        for( const auto &j : activeSet ){

            dist = norm2(m_node[j] - m_node[i]) / m_supportRadius ;
            a[k] = evalBasis( dist ) ; 
            k++;
        }
    }


    info = LAPACKE_dgesv( LAPACK_COL_MAJOR, nS, nrhs, a, lda, ipiv, b, ldb );

    if( info > 0 ) {
        printf( "The diagonal element of the triangular factor of a,\n" );
        printf( "U(%i,%i) is zero, so that a is singular;\n", info, info );
        printf( "the solution could not be computed.\n" );
        exit( 1 );
    }


    m_weight.resize(nrhs) ;

    k=0 ;
    for( j=0; j<nrhs; ++j){
        m_weight[j].resize(m_nodes,0);

        for( const auto &i : activeSet ){
            m_weight[j][i] = b[k];
            ++k;
        }
    }


    delete[] a;
    delete[] b;

};

/*! 
 * Calculates the RBF weights using active nodes
 */
void RBF::solveLSQ(){

    int i, j, k ;
    double dist;

    int nR      = getActiveCount() ;
    int nP      = m_nodes ;
    int nrhs    = getFieldCount() ;


    std::vector<int> activeSet( getActiveSet() ) ;

    int     n ,m, lda, ldb, info, rank ;
    double  rcond = -1.0 ;
    //double  rcond = 1.e-4 ;

    m = nP;
    n = nR;

    lda = m ;
    ldb = std::max(n,m) ;

    double  *a = new double [lda * n];
    double  *b = new double [ldb * nrhs];
    double  *s = new double [m];


    for( j=0; j<nrhs; ++j){
        for( i=0; i<nP; ++i){
            k = j*ldb + i ;
            b[k] = m_value[j][i];
        }
    }

    k=0;
    for( const auto &j : activeSet ){
        for( i=0; i<nP; ++i){
            dist = norm2(m_node[j] - m_node[i]) / m_supportRadius ;
            a[k] = evalBasis( dist ) ; 
            k++;
        }
    }

    info = LAPACKE_dgelsd( LAPACK_COL_MAJOR, nP, nR, nrhs, a, lda, b, ldb, s, rcond, &rank ) ;

    if( info > 0 ) {
        exit( 1 );
    }


    m_weight.resize(nrhs) ;

    for( j=0; j<nrhs; ++j){
        m_weight[j].resize(m_nodes,0);

        k=0 ;
        for( const auto &i : activeSet ){
            m_weight[j][i] = b[j*ldb+k];
            ++k;
        }
    }


    delete[] a;
    delete[] b;
    delete[] s;

};

/*! 
 * Determines set of nodes to be used using greedy algorithm
 * @param[in] tolerance error tolerance for adding nodes
 * @return true if tolerane has been met, false if not enough nodes available
 */
bool RBF::greedy( const double &tolerance){

    int                     i, j ;
    double                  error(1.e18) ;
    std::vector<double>     local(m_fields) ;

    m_error.resize(m_nodes) ;

    for( auto && active : m_active )
        active = false ;

    for( i=0; i<m_nodes; ++i){

        for( j=0; j<m_fields; ++j)
            local[j] = m_value[j][i] ;

        m_error[i] = norm2(local) ;

    }

    while( error > tolerance){
        i = addGreedyPoint() ;

        std::cout << "added node " << i ;

        if( i != -1) {
            m_active[i] = true ;

            //solve() ;
            solveLSQ() ;

            error = evalError() ;

            std::cout << std::scientific ;
            std::cout << " error now " << error << " active nodes" << getActiveCount() << " / " << m_nodes << std::endl ;
            //std::cout << " error 621 " << m_error[621] << " error 564 " << m_error[564] << std::endl ;
        } else {
            return false ;

        };
    };

    return true;

};

/*! 
 * Determines which node to ba added to active set
 * @return index with max error; if no index available -1 is returned
 */
int RBF::addGreedyPoint( ){

    int     i(0), index(-1), nA( getActiveCount() ) ;
    double  maxError(0.), penal ;
    std::array<double,3>    myCoord ;

    std::vector<int>     active( getActiveSet() );

    for( auto error : m_error ){

        if(!m_active[i] ){

            myCoord = m_node[i] ;

            if( nA != 0){ 
                penal = 1.e18 ;
                for( auto j : active ){
                    penal = std::min(penal, norm2( myCoord - m_node[j] ))  ;
                }; 
            }

            //error += 0.0001 *penal ;

            if( error > maxError ){
                maxError = error ;
                index = i;
            };
        }

        ++i;

    };


    return index;

};

/*! 
 * Calculates the relative error between rbf interpolation and exact values at nodes.
 * @return max error
 */
double RBF::evalError( ){

    int                     i(0), j(0), index ;
    double                  maxError(0), relError, realValue, norm;
    std::vector<double>     reconValues ;

    for( auto &point : m_node ){

        reconValues = evalRBF( point ) ;

        j=0;
        norm = 0. ;
        relError = 0. ;
        for( auto &val : reconValues ){
            realValue = m_value[j][i] ;
            relError += std::pow( (val - realValue), 2  ) ; 
            norm += std::pow( realValue, 2  ) ; 

            ++j;
        };

        relError = sqrt(relError); // / std::max( sqrt(norm), 1) ;
        m_error[i] = relError ;

        if( relError > maxError ){
            maxError = relError ;
            index = i ;
        }


        ++i;
    };

    //    std::cout << std::scientific ;
    //    std::cout << std::setprecision(8) ;
    //    std::cout << " max error in node: " << index << "real values:" << m_value[0][index] << " " << m_value[1][index] << " reco values:" << evalRBF( m_node[index] ) ;


    return maxError;

};

/*! 
 * Wendland C2 function
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::wendlandc2( const double &dist ){

    if( dist > 1){
        return 0.;
    } else{
        return( pow(1.-dist,4)*(4.*dist+1.) ) ;

    }

}; 

/*! 
 * @} 
 */
}
