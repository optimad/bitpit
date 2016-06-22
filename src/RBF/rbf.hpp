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


#ifndef __BITPIT_RBF_HPP__
#define __BITPIT_RBF_HPP__

#include <vector>
#include <array>


namespace bitpit{

/*!
 * @ingroup RBF
 * @{
 * 
 * @enum RBFBasisFunction  
 * @brief Enum class defining types of RBF kernel functions that could be used in bitpit::RBF class
 */
enum class RBFBasisFunction{
    WENDLANDC2 = 1,/**< Compact support Wendland C2 function */
};

/*!
 * @enum RBFMode
 * @brief Enum class defining behaviour of the bitpit::RBF class
 */
enum class RBFMode{
	INTERP = 1, /**< RBF class interpolate external field data */
	PARAM =2    /**< RBF class used as pure parameterizator*/
};

/*!
 * @}
 */

class RBF{

    private:
    int     m_fields;								/**<Number of data fields defined on RBF nodes.*/
    int     m_nodes ;								/**<Number of RBF nodes.*/
	RBFMode	m_mode;									/**<Behaviour of RBF class (interpolation or parametrization).*/
    double  m_supportRadius ;						/**<Support radius of function used as Radiabl Basis Function.*/

    double  (*m_fPtr)( const double &);

    
    std::vector<bool>                   m_active ;	/**<Vector of active/inactive node (m_active[i] = true/false -> the i-th node is used/not used during RBF evaluation).*/
    std::vector<double>                 m_error ;	/**<Interpolation error of a field evaluated on each RBF node (auxiliary memeber used in Greedy algorithm).*/

	protected:
    std::vector<std::vector<double>>    m_value ;   /**< displ value to be interpolated on RBF nodes */ 
    std::vector<std::vector<double>>    m_weight ;	/**< weight of your RBF interpolation */
	int m_maxFields; 								/**< fix the maximum number of fields that can be added to your class*/
	std::vector<std::array<double,3>>   m_node ;    /**< list of RBF nodes */
	
    public:
    ~RBF();
    RBF( RBFBasisFunction = RBFBasisFunction::WENDLANDC2 ) ;
	RBF(const RBF & other);
	RBF & operator=(const RBF & other);
	
    void                    setFunction( const RBFBasisFunction & ) ;
    void                    setFunction( double (&funct)(const double &) ) ;

    int                     getDataCount();
    int                     getActiveCount();
	int 					getTotalNodesCount();
    std::vector<int>        getActiveSet() ;

    bool                    isActive( const int &) ;

	bool					activateNode(const int &);
	bool					activateNode(const std::vector<int> &);
	void					activateAllNodes();
	bool					deactivateNode(const int &);
	bool					deactivateNode(const std::vector<int> &);
	void					deactivateAllNodes();
	
    void                    setSupportRadius( const double & ) ;
	double 					getSupportRadius();
	
	void					setMode(RBFMode mode);
	RBFMode					getMode();
	
	void                    setDataToNode ( const int &, const std::vector<double> & ) ;
	void                    setDataToAllNodes( const int &, const std::vector<double> & ) ; 
	
    int                     addNode( const std::array<double,3> & ) ;
    std::vector<int>        addNode( const std::vector<std::array<double,3>> & ) ;
	bool					removeNode(int);
	bool					removeNode(std::vector<int> &);
	void					removeAllNodes();
	
	int                     addData( ) ; 
	int                     addData( const std::vector<double> & ) ;
	bool                    removeData( int) ; 
	bool					removeData(std::vector<int> &); 
	void 					removeAllData(); 
	
	void					fitDataToNodes();
	void 					fitDataToNodes(int);
	
    std::vector<double>     evalRBF( const std::array<double,3> &) ;
	double                  evalBasis( const double &) ;

	int                    solve() ; 
	int                    greedy( const double &) ;

protected:
	
	double                  evalError() ;
	int                     addGreedyPoint() ;
	int                     solveLSQ() ;	
	
};



/*!
 * @ingroup  RadialBasisFunction
 * @brief Utility fuctions for RBF
 */
namespace rbf{
double                          wendlandc2( const double &) ;
}

}


#endif
