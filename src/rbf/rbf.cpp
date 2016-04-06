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
#include <set>
#include "Operators.hpp"
#include "SortAlgorithms.hpp"
#include "rbf.hpp"

namespace bitpit{

/*!
 * @ingroup RBF
 * @{
 *
 * @class RBF
 * @brief Handling of Radial Basis Function with a large set of nodes. 
 *
 */

/*!
 * Destructor
 */
RBF::~RBF(){
    m_fPtr = NULL ;
} ;

/*! 
 * Default constructor. Requires optionally statements of type of RBFBasisFunction 
 * which must be used and RBFType of behavior(see RBF::setType method for further information)  
 */
RBF::RBF( RBFBasisFunction bfunc, RBFType rtype ) {
    
	m_supportRadius = 1. ;
    m_nodes         = 0 ;
    m_fields        = 0 ;

    m_node.clear() ;
    m_value.clear() ;
    m_weight.clear() ;
    m_active.clear() ;

    setFunction( bfunc ) ;
	setType(rtype);
};

/*! 
 * Gets the type of RBF behaviour actually set(see RBF::setType method for further information)
 * @return RBFtype enum, type of behaviour actually set 
 */
RBFType RBF::whichType(){
	return m_rbfType;
}

/*! 
 * Sets the type of RBF behaviour of the class. Two behaviours are possible:
 * RBFType::INTERP allows the class to interpolate and automatically choose the most appropriate 
 * set the RBF weights,given an external set of data fields on prescribed RBF nodes.
 *
 * RBFType::PARAM allows the class to freely set RBF weights on given RBF nodes, without any further
 * processing.
 *  
 * According to the behaviour currently chosen, some methods of the class will be active/inactive. Such exceptions
 * will be specified in each method docs, eventually.
 *  
 * @param[in] RBFtype enum, type of chosen behaviour 
 */
void RBF::setType(RBFType type){
	m_rbfType = type;
}

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
 * Gets the number of data set attached to RBF nodes. 
 * In RBFType::PARAM mode data are meant as RBF weights
 * In RBFType::INTERP mode data are meant as fields to be interpolated
 * 
 * If no data are attached, returns -1.
 * @return  number of fields
 */
int RBF::getDataCount(  ){
	if(m_rbfType != RBFType::INTERP) return -1;
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
 * Gets the total number of nodes, active or not
 * @return  number of RBF nodes
 */
int RBF::getTotalNodesCount(  ){
	return m_nodes ;
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
 * Activate a node in your RBF node list. 
 * Method is active only in RBFType::PARAM mode.
 * @param[in] n index of node to be activated
 */
void	RBF::activateNode(const int & n){
	if(m_rbfType != RBFType::PARAM) return;
	m_active[n] = true;
};

/*! 
 * Activate a node ensamble in your RBF node list. 
 * Method is active only in RBFType::PARAM mode.
 * @param[in] list list of node indices to be activated
 */
void	RBF::activateNode(const std::vector<int> & list){
	if(m_rbfType != RBFType::PARAM) return;
	for(auto && index : list){
		 m_active[index] = true;
	}
};

/*! 
 * Activate all nodes actually available in your RBF node list. 
 * Method is active only in RBFType::PARAM mode.
 */
void	RBF::activateAllNodes(){
	if(m_rbfType != RBFType::PARAM) return;
	for(auto && active : m_active){
		active = true;
	}
};

/*! 
 * Deactivate a node in your RBF node list.
 * Method is active only in RBFType::PARAM mode.
 * @param[in] n index of node to be dactivated
 */
void	RBF::deactivateNode(const int & n ){
	if(m_rbfType != RBFType::PARAM) return;
	m_active[n] = false;
};

/*! 
 * Deactivate a node ensamble in your RBF node list.
 * Method is active only in RBFType::PARAM mode.
 * @param[in] list list of node indices to be deactivated
 */
void	RBF::deactivateNode(const std::vector<int> & list){
	if(m_rbfType != RBFType::PARAM) return;
	for(auto && index : list){
		m_active[index] = false;
	}
};
/*! 
 * Deactivate all nodes actually available in your RBF node list.
 * Method is active only in RBFType::PARAM mode.
 */
void	RBF::deactivateAllNodes(){
	if(m_rbfType != RBFType::PARAM) return;
	for(auto && active : m_active){
		active = false;
	}
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
 * Return currently set support radius
 * 
 */
double RBF::getSupportRadius(){
	return m_supportRadius ;
};


/*! 
 * Sets all the type of available data at one node. 
 * In RBFType::PARAM mode data are meant as RBF weights
 * In RBFType::INTERP mode data are meant as fields to be interpolated
 *
 * @param[in] id id of node
 * @param[in] value  data values to be set as RBF parameters for the given node
 */
void RBF::setDataToNode( const int &id, const std::vector<double> &value ){

	if(id<0 || id >= m_fields) return;

	if(value.size() != m_fields){
		std::cout<<"Mismatch dimension between value vector size and number of data attached to rbf.";
		std::cout<<"This may lead to nasty errors. Check it with getDataCount()!"<<std::endl;
		std::cout<<"Data could not be set"<<std::endl;
		return;
	}
	
	int i ;
	for( i=0; i<m_fields; ++i ){
		if(m_rbfType == RBFType::PARAM) m_weight[i][id] = value[i] ;
		if(m_rbfType == RBFType::INTERP) m_value[i][id] = value[i] ;
	}
	return ;
};

/*! 
 * Sets the values of a data set to all currently available nodes. 
 * In RBFType::PARAM mode data are meant as RBF weights
 * In RBFType::INTERP mode data are meant as fields to be interpolated
 * 
 * @param[in] id id of data
 * @param[in] value  data values 
 */
void RBF::setDataToAllNodes( const int &id, const std::vector<double> &value ){

	if(id<0 || id >= m_fields) return;
	
	int size = (m_rbfType == RBFType::PARAM)*m_weight.size() + (m_rbfType == RBFType::INTERP)*m_value[id].size();
	
	if(value.size() != size){
		std::cout<<"Mismatch dimension between data vector and current data container. One or both does not match RBF nodes count.";
		std::cout<<"This may lead to nasty errors. Use fitDataToNodes to reshape container or fit your data vector first!"<<std::endl;
		std::cout<<"Data could not be set"<<std::endl;
		return;
	}
	
	if(m_rbfType == RBFType::PARAM)		m_weight[id] = value ;
	if(m_rbfType == RBFType::INTERP)	m_value[id] = value ;
	
	return ;
};

/*! 
 * Adds a RBF node and sets it to active. Does not manage duplicated nodes 
 * @param[in] node  coordinates of node to be added
 * @return id of node within class
 */
int RBF::addNode( const std::array<double,3> &node ){
	m_node.push_back(node) ;
    m_active.push_back(true) ;

    m_nodes++ ;
    return m_nodes ;
};

/*! 
 * Adds a list of RBF nodes and sets them to active. Does not manage duplicated nodes
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

/*! Remove pre-existent node. RBF Node list is resized and renumbered after extraction. 
 * \param[in] id id of node 
 * \return boolean, true if successfully extracted, false otherwise
 */
bool RBF::removeNode(int id){
	
	if(id < 0 || id >=m_nodes) return false;
	
	m_nodes--;
	m_node.erase(m_node.begin()+id);
	m_active.erase(m_active.begin()+id);
	return(true);
}

/*! Remove pre-existent set of nodes. RBF nodal list is resized and renumbered after extraction. 
 * \param[in] list id list of candidates to extraction 
 * \return boolean, true if all nodes are successfully extracted, false if any of them or none are extracted
 */
bool RBF::removeNode(std::vector<int> & list){
	
	std::set<int> setList;
	for(auto && id : list) setList.insert(id);
	
	int extracted = 0;
	for(auto && id : setList){
		if(id>=0 && id <m_nodes){;
			m_nodes--;
			int index = id-extracted;
			m_node.erase(m_node.begin() + index);
			m_active.erase(m_active.begin() + index);
			extracted++;
			}	
	}
	return(extracted == list.size());
};

/*!
 * Remove all nodes in RBF nodal list
 */
void RBF::removeAllNodes(){
	m_nodes = 0;
	m_node.clear();
	m_active.clear();
};


/*! 
 * Increment container size for RBF control data.The RBF::fitDataToNodes() method 
 * is implicitly called, to ensure dimension consistency between data dimension 
 * and number of RBF nodes. Use RBF::setDataToAllNodes to fill them. 
 * In RBFType::PARAM mode data are meant as RBF weights
 * In RBFType::INTERP mode data are meant as fields to be interpolated
 * @return id of virtual data within the class
 */
int RBF::addData( ){
	m_fields++ ;
	fitDataToNodes(m_fields-1);
	return m_fields;
};

/*! 
 * Adds data attached to RBF nodes to current set.
 * In RBFType::PARAM mode data are meant as RBF weights
 * In RBFType::INTERP mode data are meant as fields to be interpolated
 * @param[in] data values of weight/fields for each RBF node
 * @return id of data within the class
 */
int RBF::addData( const std::vector<double> & data ){
		
	if(data.size() != m_nodes){
		std::cout<<"Mismatch dimension between data vector and actual RBF nodes count. This may lead to nasty errors."<<std::endl;
		std::cout<<"Data could not be added"<<std::endl;
		return(-1);
	}
	
	if(m_rbfType == RBFType::PARAM)	m_weight.push_back(data) ;
	if(m_rbfType == RBFType::INTERP)	m_value.push_back(data) ;
	m_fields++ ;
	return m_fields ;
};




/*! Remove pre-existent data set. Data list is resized and renumbered after extraction. 
 * In RBFType::PARAM mode data are meant as RBF weights
 * In RBFType::INTERP mode data are meant as fields to be interpolated
 * 
 * \param[in] id id of node 
 * \return boolean, true if successfully extracted, false otherwise
 */
bool RBF::removeData(int id){
	
	if(id<0 || id >=m_fields) return false;
	
	m_fields--;
	m_weight.erase(m_weight.begin() + id);
	if(m_rbfType == RBFType::INTERP) m_value.erase(m_value.begin()+id);
	return(true);
}

/*! Remove pre-existent set of data. RBF Data list is resized and renumbered after extraction. 
 * In RBFType::PARAM mode data are meant as RBF weights
 * In RBFType::INTERP mode data are meant as fields to be interpolated
 *  
 * \param[in] list id list of candidates to extraction 
 * \return boolean, true if all data set are successfully extracted, false if any of them or none are extracted
 */
bool RBF::removeData(std::vector<int> & list){
	
	std::set<int> setList;
	for(auto && id : list) setList.insert(id);
	
	int extracted = 0;
	for(auto && id : setList){
		if(id>=0 && id <m_fields){;
			m_fields--;
			int index = id-extracted; 
			m_weight.erase(m_weight.begin()+index);
			if(m_rbfType == RBFType::INTERP) m_value.erase(m_value.begin()+index);
			extracted++;
			}	
	}
	return(extracted == list.size());
};

/*!
 * Remove all data set in RBF nodal list
 */
void RBF::removeAllData(){
	m_fields = 0;
	m_weight.clear();
	m_value.clear();
};


// /*! 
//  * Sets all the field values at one node
//  * @param[in] id id of node
//  * @param[in] value  function values to be interpolated
//  */
// void RBF::setNodeValue( const int &id, const std::vector<double> &value ){
// 
//     int i ;
// 
//     for( i=0; i<m_fields; ++i ){
//         m_value[i][id] = value[i] ;
//     }
//     return ;
// };
// 
// /*! 
//  * Sets the values of one field at all nodes
//  * @param[in] id id of field
//  * @param[in] value  function values to be interpolated
//  */
// void RBF::setFieldValue( const int &id, const std::vector<double> &value ){
// 
//     m_value[id] = value ;
// 
//     return ;
// };


/*! 
 * Evaluates the RBF
 * @param[in] point point where to evaluate the basis
 * @return vector containing interpolated values. Its size matches the number of fields/weights of RBF
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
 * Calculates the RBF weights using all active nodes and just given target fields. 
 * Regular LU solver for linear system A*X=B is employed (LAPACKE dgesv). The method
 * is active only in RBFType::INTERP mode.
 */
void RBF::solve(){

	if(m_rbfType != RBFType::INTERP) return;
	
    int i, j, k ;
    double dist;

    int nS      = getActiveCount() ;
    int nrhs    = getDataCount() ;

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
 * Determines effective set of nodes to be used using greedy algorithm and calculate weights on them.
 * The method is active only in RBFType::INTERP mode.
 * @param[in] tolerance error tolerance for adding nodes
 * @return true if tolerance has been met, false if not enough nodes available
 */
bool RBF::greedy( const double &tolerance){
	
	if(m_rbfType != RBFType::INTERP) return false;
	
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
 * Check dimensions of already available data and resize them to current 
 * RBF node list dimension. The method resizes all data structures to current RBF 
 * node list dimension and does not destroy any previous stored data within such 
 * dimension. Anyway, mismatches definitions could occur. Please 
 * use RBF::setDataToAllNodesto load your data again. 
 * In RBFType::PARAM mode data are meant as RBF weights
 * In RBFType::INTERP mode data are meant as fields to be interpolated
 */
void RBF::fitDataToNodes(){

	for (int i=0;i<m_fields; ++i){
		fitDataToNodes(i);
	}
}

/*!
 * Check dimensions id-th data and resize it to current 
 * RBF node list dimension. The method resizes all data structures to current RBF 
 * node list dimension and does not destroy any previous stored data within such 
 * dimension. Anyway, mismatches definitions could occur. Please 
 * use RBF::setDataToAllNodes to load your data again. 
 * In RBFType::PARAM mode data are meant as RBF weights
 * In RBFType::INTERP mode data are meant as fields to be interpolated
 * \param[in] id id of data
 */
void RBF::fitDataToNodes(int id){
	
	m_weight[id].resize(m_nodes);
	if(m_rbfType == RBFType::INTERP){
		m_value[id].resize(m_nodes);
	}
}


//PROTECTED RBF CLASS METHODS IMPLEMENTATION

/*! 
 * Evaluates the basis function
 * @param[in] dist distance
 * @return value of basis function
 */
double RBF::evalBasis( const double &dist ){
	return (*m_fPtr)(dist) ;
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
 * Calculates the RBF weights using all active nodes and just given target fields. 
 * Compute weights as solution of a linear least squares problem (LAPACKE dglsd).
 */
void RBF::solveLSQ(){
	
	int i, j, k ;
	double dist;
	
	int nR      = getActiveCount() ;
	int nP      = m_nodes ;
	int nrhs    = getDataCount() ;
	
	
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

//RBF NAMESPACE UTILITIES 

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
