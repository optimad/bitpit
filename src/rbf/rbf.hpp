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


/*!
 * Enum Class defining type of behaviors handled by RBF object. 
 * INTERP flag activate usage of the class as a pure interpolator of given data fields
 * PARAM flag activate usage of RBF functions for spatial 1D/2D/3D parameterization purposes
 */
enum class RBFType{
	INTERP = 1,
	PARAM  = 2
};


class RBF{

    private:
    int     m_fields ;
    int     m_nodes ;

    double  m_supportRadius ;

    double  (*m_fPtr)( const double &);

    std::vector<std::array<double,3>>   m_node ;
    std::vector<std::vector<double>>    m_value ;
    std::vector<std::vector<double>>    m_weight ;

    std::vector<bool>                   m_active ;
    std::vector<double>                 m_error ;

	RBFType m_rbfType;
	
    public:
    ~RBF();

    RBF( RBFBasisFunction = RBFBasisFunction::WENDLANDC2, RBFType = RBFType::INTERP ) ;

	RBFType					whichType();
	void					setType(RBFType);
	
    void                    setFunction( const RBFBasisFunction & ) ;
    void                    setFunction( double (&funct)(const double &) ) ;

    int                     getFieldCount();
	int                     getWeightCount();
    int                     getActiveCount() ;
    std::vector<int>        getActiveSet() ;

    bool                    isActive( const int &) ;

	void					activateNode(const int &);
	void					activateNode(const std::vector<int> &);
	void					activateAllNodes();
	void					deactivateNode(const int &);
	void					deactivateNode(const std::vector<int> &);
	void					deactivateAllNodes();
	
    void                    setSupportRadius( const double & ) ;
	
    void                    setFieldsToNode ( const int &, const std::vector<double> & ) ;
    void                    setFieldToAllNodes( const int &, const std::vector<double> & ) ;
	void 					setWeightsToNode(const int &, const std::vector<double> & );
	void 					setWeightToAllNodes(const int &, const std::vector<double> & );
	
    int                     addNode( const std::array<double,3> & ) ;
    std::vector<int>        addNode( const std::vector<std::array<double,3>> & ) ;

    int                     addField( ) ;
    int                     addField( const std::vector<double> & ) ;

	int						addWeight();
	int 					addWeight(const std::vector<double> & );
	
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
