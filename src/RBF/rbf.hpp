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
 * @ingroup VTKEnums
 * Enum class defining types of fields which may be written through class VTK
 */
enum class RBFBasisFunction{
    WENDLANDC2 = 1,
};


class RBF{

    private:
    int     m_fields ;
    int     m_nodes ;

    double  m_supportRadius ;

    double  (*m_fPtr)( const double &);

    
    std::vector<bool>                   m_active ;
    std::vector<double>                 m_error ;

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

    void                    solve() ;
	bool                    greedy( const double &) ;

protected:
	
	double                  evalError() ;
	double                  initGreedy( const int &) ;
	int                     addGreedyPoint() ;
	void                    solveLSQ() ;	
	
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
