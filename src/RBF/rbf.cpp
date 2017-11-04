/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
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

#include <cmath>
#include <set>

#include "bitpit_private_lapacke.hpp"

#include "bitpit_operators.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_SA.hpp"

#include "rbf.hpp"

namespace bitpit {

/*!
 * @class RBFKernel
 * @ingroup RBF
 * @brief Base class to handle Radial Basis Function with a large set of nodes.
 *
 * The class can be used in two different ways:
 *  - as interpolator, given a set of external fields
 *  - as parameterizator, directly modifying the weights of rbf kernels.
 *
 * The User can switch between modes, according to its needs.
 * Some internal methods of the class can change their behaviour according
 * to the class mode selected. Please check documentation of each single method to
 * appreciate the differences. Default mode of the class is the "INTERP" RBFMode::
 * The actual abstract class does not implement the type of node you want to use; it
 * can be a 3D point or an unstructured mesh. Anyway what you have to do in deriving your own
 * RBF class is to specify a type of node you want to use and implement an evaluator
 * of distances node/node or anywhere point in space/node.
 *
 */

/*!
 * Destructor
 */
RBFKernel::~RBFKernel()
{
    m_fPtr = NULL;
}

/*!
 * Default constructor. RBFBasisFunction is WENDLANDC2 by default. RBFMode is
 * INTERP, by default. Use setFunction and setMode for changing it.
 */
RBFKernel::RBFKernel()
{
    m_supportRadius = 1.;
    m_nodes         = 0;
    m_fields        = 0;

    m_mode = RBFMode::INTERP;

    m_maxFields = -1;
    m_value.clear();
    m_weight.clear();
    m_activeNodes.clear();

    setFunction( RBFBasisFunction::WENDLANDC2);
}

/*!
 * Copy Constructor
 */
RBFKernel::RBFKernel(const RBFKernel & other)
    : m_fields(other.m_fields), m_mode(other.m_mode),
      m_supportRadius(other.m_supportRadius), m_typef(other.m_typef),
      m_fPtr(other.m_fPtr), m_error(other.m_error), m_value(other.m_value),
      m_weight(other.m_weight), m_activeNodes(other.m_activeNodes),
      m_maxFields(other.m_maxFields), m_nodes(other.m_nodes)
{
}

/*!
 * Swap method. Exchange contents of each class member with those 
 * corresponding in the argument object.
 * \param[in] x object to be swapped
 */
void RBFKernel::swap(RBFKernel &x) noexcept
{
   std::swap(m_fields, x.m_fields);
   std::swap(m_mode, x.m_mode);
   std::swap(m_supportRadius, x.m_supportRadius);
   std::swap(m_typef, x.m_typef);
   std::swap(m_fPtr, x.m_fPtr);
   std::swap(m_error, x.m_error);
   std::swap(m_value, x.m_value);
   std::swap(m_weight, x.m_weight);
   std::swap(m_activeNodes, x.m_activeNodes);
   std::swap(m_maxFields, x.m_maxFields);
   std::swap(m_nodes, x.m_nodes);
}

/*!
 * Sets the rbf function to be used. Supported in both modes.
 * @param[in] bfunc basis function to be used
 */
void RBFKernel::setFunction( const RBFBasisFunction &bfunc )
{
    switch(bfunc){

    case( RBFBasisFunction::WENDLANDC2):
        setFunction( rbf::wendlandc2);
        m_typef = RBFBasisFunction::WENDLANDC2;
        break;

    case( RBFBasisFunction::LINEAR):
        setFunction( rbf::linear);
        m_typef = RBFBasisFunction::LINEAR;
        break;

    case( RBFBasisFunction::GAUSS90):
        setFunction( rbf::gauss90);
        m_typef = RBFBasisFunction::GAUSS90;
        break;

    case( RBFBasisFunction::GAUSS95):
        setFunction( rbf::gauss95);
        m_typef = RBFBasisFunction::GAUSS95;
        break;

    case( RBFBasisFunction::GAUSS99):
        setFunction( rbf::gauss99);
        m_typef = RBFBasisFunction::GAUSS99;
        break;

    case( RBFBasisFunction::C1C0):
        setFunction( rbf::c1c0);
        m_typef = RBFBasisFunction::C1C0;
        break;

    case( RBFBasisFunction::C2C0):
        setFunction( rbf::c2c0);
        m_typef = RBFBasisFunction::C2C0;
        break;

    case( RBFBasisFunction::C0C1):
        setFunction( rbf::c0c1);
        m_typef = RBFBasisFunction::C0C1;
        break;

    case( RBFBasisFunction::C1C1):
        setFunction( rbf::c1c1);
        m_typef = RBFBasisFunction::C1C1;
        break;

    case( RBFBasisFunction::C2C1):
        setFunction( rbf::c2c1);
        m_typef = RBFBasisFunction::C2C1;
        break;

    case( RBFBasisFunction::C0C2):
        setFunction( rbf::c0c2);
        m_typef = RBFBasisFunction::C0C2;
        break;

    case( RBFBasisFunction::C1C2):
        setFunction( rbf::c1c2);
        m_typef = RBFBasisFunction::C1C2;
        break;

    case( RBFBasisFunction::C2C2):
        setFunction( rbf::c2c2);
        m_typef = RBFBasisFunction::C2C2;
        break;

    default:
        setFunction( rbf::wendlandc2);
        m_typef = RBFBasisFunction::WENDLANDC2;
        break;
    }
}

/*!
 * Sets the rbf function to a user specified function. Supported in both modes.
 * @param[in] bfunc basis function to be used
 */
void RBFKernel::setFunction( double (&bfunc)(const double &) )
{
    m_fPtr = bfunc;
    m_typef = RBFBasisFunction::CUSTOM;
}

/*!
 * Gets the type of RBFBasisFunction linked to the class.
 * @return  type of RBFBasisFunction linked
 */
RBFBasisFunction RBFKernel::getFunctionType(  )
{
    return m_typef;
}

/*!
 * Gets the number of data set attached to RBFKernel nodes.
 * In INTERP mode, it is the number of different field that need to be interpolated;
 * In PARAM  mode, it identifies the dimension of the weights array for each node.
 * @return  number of data
 */
int RBFKernel::getDataCount(  )
{
    return m_fields;
}

/*!
 * Get the number of active nodes. Supported in both nodes
 * @return  number of active nodes
 */
int RBFKernel::getActiveCount(  )
{
    int nActive(0);

    for( auto && active : m_activeNodes)
        nActive += (int) active;

    return nActive;
}

/*!
 * Get the indices of the active nodes. Supported in both modes.
 * @return  indices of active nodes
 */
std::vector<int> RBFKernel::getActiveSet(  )
{
    int                 i(0);
    std::vector<int>    activeSet;

    activeSet.reserve( getActiveCount() );

    for( auto && active : m_activeNodes) {
        if( active )
            activeSet.push_back(i);
        ++i;
    }

    return activeSet;
}

/*!
 * Checks if a node is active. Supported in both modes.
 * @param[in] n index of node to be checked
 * @return true if active
 */
bool RBFKernel::isActive( const int &n )
{
    if(n<0 || n >= int(m_activeNodes.size())) {
        return false;
    }

    return m_activeNodes[n];
}

/*!
 * Activate a node in your RBFKernel node list. Supported in both modes.
 * @param[in] n index of node to be activated
 * @return    boolean, true if node is activated successfully, false if not
 */
bool RBFKernel::activateNode(const int & n)
{
    bool check = false;
    if(n>=0 && n<m_nodes) {
        m_activeNodes[n] = true;
        check = true;
    }
    return check;
}

/*!
 * Activate a node ensamble in your RBFKernel node list. Supported in both modes.
 * @param[in] list list of node indices to be activated
 * @return    boolean, true if all nodes are activated successfully, false if at least one of them is not
 */
bool RBFKernel::activateNode(const std::vector<int> & list)
{
    if(list.empty()) {
        return false;
    }

    bool check = true;
    for(auto && index : list) {
        check = check && activateNode(index);
    }
    return check;
}

/*!
 * Activate all nodes actually available in your RBFKernel node list. Supported in both modes.
 */
void RBFKernel::activateAllNodes()
{
    for(auto && active : m_activeNodes) {
        active = true;
    }
}

/*!
 * Deactivate a node in your RBFKernel node list.  Supported in both modes.
 * @param[in] n index of node to be dactivated
 * @return    boolean, true if node is activated successfully, false if not
 */
bool RBFKernel::deactivateNode(const int & n )
{
    bool check = false;
    if(n>=0 && n<m_nodes) {
        m_activeNodes[n] = false;
        check=  true;
    }
    return check;
}

/*!
 * Deactivate a node ensamble in your RBFKernel node list.  Supported in both modes.
 * @param[in] list list of node indices to be deactivated
 * @return    boolean, true if all nodes are deactivated successfully, false if at least one of them is not
 */
bool RBFKernel::deactivateNode(const std::vector<int> & list)
{
    if(list.empty()) {
        return false;
    }

    bool check = true;
    for(auto && index : list) {
        check = check && deactivateNode(index);
    }
    return check;
}

/*!
 * Deactivate all nodes actually available in your RBFKernel node list.  Supported in both modes.
 */
void RBFKernel::deactivateAllNodes()
{
    for(auto && active : m_activeNodes) {
        active = false;
    }
}

/*!
 * Set the support radius of all RBFKernel kernel functions. Supported in both modes.
 * @param[in] radius support radius
 */
void RBFKernel::setSupportRadius( const double & radius )
{
    m_supportRadius = radius;
    return;
}

/*!
 * Return currently set support radius used by all RBFKernel kernel functions.
 * Supported in both modes.
 * @return support radius
 */
double RBFKernel::getSupportRadius()
{
    return m_supportRadius;
}

/*!
 * Return currently usage mode of your class.
 * @return class mode
 */
RBFMode RBFKernel::getMode()
{
    return m_mode;
}

/*!
 * Set usage mode of your class.
 * @param[in] mode class mode. Ref to RBFMode enum
 */
void RBFKernel::setMode(RBFMode mode)
{
    m_mode = mode;
}

/*!
 * Sets all the type of available data at one node.
 * In INTERP mode, set each field value at the target node
 * In PARAM mode, set the weight coefficient array at node
 * @param[in] id id of node
 * @param[in] value  data values to be set as RBFKernel parameters for the given node
 */
void RBFKernel::setDataToNode( const int &id, const std::vector<double> &value )
{
    if(id<0 || id >= m_fields) {
        return;
    }

    if((int)(value.size()) != m_fields) {
        log::cout() << "Mismatch dimension between value vector size and number of data attached to rbf.";
        log::cout() << "This may lead to nasty errors. Check it with getDataCount()!" << std::endl;
        log::cout() << "Data could not be set" << std::endl;
        return;
    }

    int i;
    if(m_mode != RBFMode::PARAM) {
        for( i=0; i<m_fields; ++i ) {
                m_value[i][id] = value[i];
        }
    }else{
        for( i=0; i<m_fields; ++i ) {
            m_weight[i][id] = value[i];
        }
    }
}

/*!
 * Sets the values of a data set to all currently available nodes.
 * In INTERP mode, set one field to all nodes
 * In PARAM mode, set one weight array component  for all nodes.
 * @param[in] id id of data
 * @param[in] value  data values
 */
void RBFKernel::setDataToAllNodes( const int &id, const std::vector<double> &value )
{
    if(id<0 || id >= m_fields) {
        return;
    }

    int size = m_value[id].size();

    if((int)(value.size()) != size) {
        log::cout() << "Mismatch dimension between data vector and current data container. One or both does not match RBFKernel nodes count.";
        log::cout() << "This may lead to nasty errors. Use fitDataToNodes to reshape container or fit your data vector first!" << std::endl;
        log::cout() << "Data could not be set" << std::endl;
        return;
    }
    if(m_mode != RBFMode::PARAM) {
        m_value[id] = value;
    }else{
        m_weight[id] = value;
    }
}

/*!
 * Increment container size for RBFKernel control data.The RBFKernel::fitDataToNodes() method
 * is implicitly called, to ensure dimension consistency between data dimension
 * and number of RBFKernel nodes. Use RBFKernel::setDataToAllNodes to fill them.
 * In INTERP mode, increments and fits fields container, in PARAM mode, the weights one.
 * @return id of virtual data within the class
 */
int RBFKernel::addData( )
{
    if(m_fields == m_maxFields) {
        log::cout() << "Maximum number of data set reached" << std::endl;
        return -1;
    }
    m_fields++;
    fitDataToNodes(m_fields-1);
    return m_fields;
}

/*!
 * Adds data attached to RBFKernel nodes to current set, a field to be interpolated in INTERP mode,
 * a RBFKernel weight component in PARAM mode.
 * Note: data vector is added even if its size is different from actual number of RBFKernel nodes.
 * To ensure consistency use fitDataToNodes() method.
 *
 * @param[in] data values of weight/fields for each RBFKernel node
 * @return id of data within the class
 *
 */
int RBFKernel::addData( const std::vector<double> & data )
{
    if(m_fields == m_maxFields) {
        log::cout() << "Maximum number of data set reached" << std::endl;
        return -1;
    }
    if(m_mode == RBFMode::INTERP)    m_value.push_back(data);
    else                            m_weight.push_back(data);
    m_fields++;

    return m_fields;
}

/*! Remove pre-existent data set. Data list is resized and renumbered after extraction.
 *  Remove fields to be interpolated in INTERP mode, weights component in PARAM mode.
 * \param[in] id id of node
 * \return boolean, true if successfully extracted, false otherwise
 */
bool RBFKernel::removeData(int id)
{
    if(id<0 || id >=m_fields) {
        return false;
    }

    m_fields--;
    if(m_mode == RBFMode::INTERP)    m_value.erase(m_value.begin()+id);
    else                            m_weight.erase(m_value.begin()+id);

    return(true);
}

/*! Remove pre-existent set of data. RBF Data list is resized and renumbered after extraction.
 * In PARAM mode data are meant as RBF weights
 * In INTERP mode data are meant as fields to be interpolated
 *
 * \param[in] list id list of candidates to extraction
 * \return boolean, true if all data set are successfully extracted, false if any of them are not extracted
 */
bool RBFKernel::removeData(std::vector<int> & list)
{
    std::set<int> setList;
    for(auto && id : list) setList.insert(id);

    int extracted = 0;
    for(auto && id : setList) {
        if(id>=0 && id <m_fields){;
            m_fields--;
            int index = id-extracted;
            if (m_mode == RBFMode::INTERP) {
                auto valueItr = m_value.begin() + index;
                m_value.erase(valueItr);
            } else {
                auto weightItr = m_weight.begin() + index;
                m_weight.erase(weightItr);
            }

            extracted++;
        }
    }

    return(extracted == (int)(list.size()));
}

/*!
 * Remove all data set in RBF list. All fields for interpolation in INTERP mode,
 * or all RBF weights for PARAM mode
 */
void RBFKernel::removeAllData()
{
    m_fields = 0;
    if(m_mode == RBFMode::INTERP)    m_value.clear();
    else                            m_weight.clear();
}

/*!
 * Evaluates the RBF. Supported in both modes.
 * Its size matches the number of fields/weights attached to RBF.
 *
 * @param[in] point point where to evaluate the basis
 * @return vector containing interpolated/parameterized values.
 *
 */
std::vector<double> RBFKernel::evalRBF( const std::array<double,3> &point)
{
    std::vector<double> values(m_fields, 0.);
    int                 i, j;
    double              dist, basis;

    for( i=0; i<m_nodes; ++i ){
        if( m_activeNodes[i] ) {
            dist = calcDist(point, i) / m_supportRadius;
            basis = evalBasis( dist );

            for( j=0; j<m_fields; ++j) {
                values[j] += basis * m_weight[j][i];
            }
        }
    }

    return values;
}

/*!
 * Evaluates the RBF on a target RBF node. Supported in both modes.
 * Its size matches the number of fields/weights attached to RBF.
 *
 * @param[in] jnode index of the RBF node in the node list
 * @return vector containing interpolated/parameterized values.
 *
 */
std::vector<double> RBFKernel::evalRBF(int jnode)
{
    std::vector<double> values(m_fields, 0.);
    int                 i, j;
    double              dist, basis;

    if(jnode<0 || jnode>= m_nodes ) {
        return values;
    }

    for( i=0; i<m_nodes; ++i ) {
        if( m_activeNodes[i] ) {
            dist = calcDist(jnode, i) / m_supportRadius;
            basis = evalBasis( dist );

            for( j=0; j<m_fields; ++j) {
                values[j] += basis * m_weight[j][i];
            }
        }
    }

    return values;
}

/*!
 * Calculates the RBF weights using all currently active nodes and just given target fields.
 * Regular LU solver for linear system A*X=B is employed (LAPACKE dgesv).
 * Supported ONLY in INTERP mode.
 *
 * @return integer error flag . If 0-successfull computation, if 1-errors occurred, if -1 dummy method call
 */
int RBFKernel::solve()
{
    if(m_mode == RBFMode::PARAM) {
        return -1;
    }

    int  j, k;
    double dist;

    int nS      = getActiveCount();
    int nrhs    = getDataCount();

    int lda     = nS;
    int ldb     = nS;
    int info;
    int ipiv[nS];

    std::vector<int> activeSet( getActiveSet() );

    double *a = new double [lda * nS];
    double *b = new double [ldb * nrhs];

    k=0;
    for( j=0; j<nrhs; ++j) {
        for( const auto & i : activeSet ) {
            b[k] = m_value[j][i];
            ++k;
        }
    }

    k=0;
    for( const auto &i : activeSet ) {
        for( const auto &j : activeSet ){

            dist = calcDist(j,i) / m_supportRadius;
            a[k] = evalBasis( dist );
            k++;
        }
    }

    info = LAPACKE_dgesv( LAPACK_COL_MAJOR, nS, nrhs, a, lda, ipiv, b, ldb );

    if( info > 0 ) {
        printf( "The diagonal element of the triangular factor of a,\n" );
        printf( "U(%i,%i) is zero, so that a is singular;\n", info, info );
        printf( "the solution could not be computed.\n" );
        return 1;
    }

    m_weight.resize(nrhs);

    k=0;
    for( j=0; j<nrhs; ++j) {
        m_weight[j].resize(m_nodes,0);

        for( const auto &i : activeSet ) {
            m_weight[j][i] = b[k];
            ++k;
        }
    }

    delete[] a;
    delete[] b;

    return 0;
}

/*!
 * Determines effective set of nodes to be used using greedy algorithm and calculate weights on them.
 * Automatically choose which set of RBFKernel nodes is active or not, according to the given tolerance.
 * Supported ONLY in INTERP mode.
 * @param[in] tolerance error tolerance for adding nodes
 * @return integer error flag . If 0-successfull computation and tolerance met, if 1-errors occurred, not enough nodes, if -1 dummy method call
 */
int RBFKernel::greedy( const double &tolerance)
{
    if(m_mode == RBFMode::PARAM)    return -1;

    int                     i, j;
    double                  error(1.e18);
    std::vector<double>     local(m_fields);

    m_error.resize(m_nodes);

    for( auto && active : m_activeNodes )
        active = false;

    for( i=0; i<m_nodes; ++i){

        for( j=0; j<m_fields; ++j) {
            local[j] = m_value[j][i];
        }

        m_error[i] = norm2(local);

    }

    ios::fmtflags streamFlags(log::cout().flags());

    int errorFlag = 0;
    while( error > tolerance) {
        i = addGreedyPoint();

        if( i != -1) {
            m_activeNodes[i] = true;

            //solve();
            solveLSQ();

            error = evalError();

            log::cout() << std::scientific;
            log::cout() << " error now " << error << " active nodes" << getActiveCount() << " / " << m_nodes << std::endl;
        } else {
            errorFlag = 1;
            break;

        }
    }

    log::cout().flags(streamFlags);

    return errorFlag;
}

/*!
 * Check dimensions of already available data and resize them to current
 * RBFKernel node list dimension. The method resizes all data structures to current RBFKernel
 * node list dimension and does not destroy any previous stored data within such
 * dimension. Anyway, mismatches definitions could occur. Please
 * use RBFKernel::setDataToAllNodesto load your data again.
 * In RBFMode::PARAM mode data are meant as RBFKernel weights
 * In RBFMode::INTERP mode data are meant as fields to be interpolated
 */
void RBFKernel::fitDataToNodes()
{
    for (int i=0;i<m_fields; ++i) {
        fitDataToNodes(i);
    }
}

/*!
 * Check dimensions id-th data and resize it to current
 * RBFKernel node list dimension. The method resizes id-th data structure to current RBFKernel
 * node list dimension and does not destroy any previous stored data within such
 * dimension. Anyway, mismatches definitions could occur. Please
 * use RBFKernel::setDataToAllNodes to load your data again.
 * In RBFMode::PARAM mode data are meant as RBFKernel weights
 * In RBFMode::INTERP mode data are meant as fields to be interpolated
 * @param[in] id id of data
 */
void RBFKernel::fitDataToNodes(int id)
{
    if(m_mode != RBFMode::PARAM)    m_value[id].resize(m_nodes, 0.0);
    else                            m_weight[id].resize(m_nodes,0.0);
}


//PROTECTED RBFKernel CLASS METHODS IMPLEMENTATION

/*!
 * Evaluates the basis function. Supported in both modes
 * @param[in] dist distance
 * @return value of basis function
 */
double RBFKernel::evalBasis( const double &dist )
{
    return (*m_fPtr)(dist);
}

/*!
 * Determines which node has to be added to active set. Supported only in INTERP mode.
 * @return index with max error; if no index available, or dummy call -1 is returned
 */
int RBFKernel::addGreedyPoint( )
{
    if(m_mode == RBFMode::PARAM) {
        return -1;
    }

    int     i(0), index(-1);
    double  maxError(0.); 

    std::vector<int>     active( getActiveSet() );

    for( auto error : m_error ){

        if(!m_activeNodes[i] ){

            if( error > maxError ) {
                maxError = error;
                index = i;
            }
        }

        ++i;

    }

    return index;
}

/*!
 * Calculates the relative error between rbf interpolation and exact values at nodes.
 * Supported only in INTERP mode.
 * @return max error, if lesser then 0, dummy call triggered.
 */
double RBFKernel::evalError( )
{
    if(m_mode == RBFMode::PARAM) {
        return -1.0;
    }

    int                     i(0), j(0);
    //int                     index;
    double                  maxError(0), relError, realValue, norm;
    std::vector<double>     reconValues;

    for( int iX=0; iX< m_nodes; ++iX ) {
        reconValues = evalRBF( iX );

        j=0;
        norm = 0.;
        relError = 0.;
        for( auto &val : reconValues ) {
            realValue = m_value[j][i];
            relError += std::pow( (val - realValue), 2  );
            norm += std::pow( realValue, 2  );

            ++j;
        }

        relError = sqrt(relError);
        m_error[i] = relError;

        if( relError > maxError ) {
            maxError = relError;
        }

        ++i;
    }

    return maxError;
}

/*!
 * Calculates the RBFKernel weights using all active nodes and just given target fields.
 * Compute weights as solution of a linear least squares problem (LAPACKE dglsd).
 * Supported ONLY in INTERP mode.
 * \return integer error flag . If 0-successfull computation, if 1-errors occurred , -1 dummy call
 */
int RBFKernel::solveLSQ()
{
    if(m_mode == RBFMode::PARAM) {
        return -1;
    }

    int i, j, k;
    double dist;

    int nR      = getActiveCount();
    int nP      = m_nodes;
    int nrhs    = getDataCount();

    std::vector<int> activeSet( getActiveSet() );

    int     n ,m, lda, ldb, info, rank;
    double  rcond = -1.0;
    //double  rcond = 1.e-4;

    m = nP;
    n = nR;

    lda = m;
    ldb = std::max(n,m);

    double  *a = new double [lda * n];
    double  *b = new double [ldb * nrhs];
    double  *s = new double [m];

    for( j=0; j<nrhs; ++j) {
        for( i=0; i<nP; ++i) {
            k = j*ldb + i;
            b[k] = m_value[j][i];
        }
    }

    k=0;
    for( const auto &j : activeSet ) {
        for( i=0; i<nP; ++i) {
            dist = calcDist(j,i) / m_supportRadius;
            a[k] = evalBasis( dist );
            k++;
        }
    }

    info = LAPACKE_dgelsd( LAPACK_COL_MAJOR, nP, nR, nrhs, a, lda, b, ldb, s, rcond, &rank );

    if( info > 0 ) {
        return 1;
    }

    m_weight.resize(nrhs);

    for( j=0; j<nrhs; ++j) {
        m_weight[j].clear();
        m_weight[j].resize(m_nodes,0);

        k=0;
        for( const auto &i : activeSet ) {
            m_weight[j][i] = b[j*ldb+k];
            ++k;
        }
    }

    delete[] a;
    delete[] b;
    delete[] s;

    return(0);
}

// RBF derived class implementation

/*!
 * @class RBF
 * @ingroup RBF
 * @brief Class to handle Radial Basis Function with a large set of 3D points as nodes.
 *
 * The class specializes RBFKernel and employs a 3D point cloud as set of RBF nodes.
 *
 */

/*!
 * Destructor
 */
RBF::~RBF()
{
}

// /*!
//  * Default constructor. RBFBasisFunction::WENDLANDC2 is default. RBFMode is
//  * INTERP, by default. Use setMode for changing it.
//  */
// RBF::RBF(): RBFKernel()
// {
//     m_node.clear();
// }

/*!
 * Default constructor. Requires optionally statements of type of RBFBasisFunction
 * which must be used. RBFMode is INTERP, by default. Use setMode for changing it.
 */
RBF::RBF( RBFBasisFunction bfunc)
{
    m_node.clear();
    setFunction( bfunc );
}

/*!
 * Copy Constructor
 */
RBF::RBF(const RBF & other)
    : RBFKernel(other),
      m_node(other.m_node)
{
}

/*!
 * Copy Operator
 */
RBF & RBF::operator=(RBF other)
{
    swap(other);

    return *this;
}

/*!
 * Swap method. Exchange contents of each class member with those 
 * corresponding in the argument object.
 * \param[in] x object to be swapped
 */
void RBF::swap(RBF &x) noexcept
{
    RBFKernel::swap(x);

    std::swap(m_node, x.m_node);
}

/*!
 * Gets the total number of nodes, active or not. Supported in both modes.
 * @return  number of available RBF nodes
 */
int RBF::getTotalNodesCount(  )
{
    return m_nodes;
}

/*!
 * Adds a RBF node and sets it to active. Does not manage duplicated nodes.
 * Supported in both modes.
 * @param[in] node  coordinates of node to be added
 * @return id of node within class
 */
int RBF::addNode( const std::array<double,3> &node )
{
    m_node.push_back(node);
    m_activeNodes.push_back(true);
    m_nodes++;

    return m_nodes;
}

/*!
 * Adds a list of RBF nodes and sets them to active. Does not manage duplicated nodes.
 * Supported in both modes.
 * @param[in] node  coordinates of nodes to be added
 * @return id of node within class
 */
std::vector<int> RBF::addNode( const std::vector<std::array<double,3>> &node )
{
    int                 i( m_nodes );
    std::vector<int>    ids;

    ids.resize( node.size() );

    for( auto & id:ids ) {
        id = i;
        ++i;
    }

    m_node.insert( m_node.end(), node.begin(), node.end() );
    m_nodes += node.size();

    m_activeNodes.resize( m_nodes, true );

    return ids;
}

/*! Remove pre-existent node. RBF Node list is resized and renumbered after extraction.
 * Supported in both modes.
 * @param[in] id id of node
 * @return boolean, true if successfully extracted, false otherwise
 */
bool RBF::removeNode(int id)
{
    if(id < 0 || id >=m_nodes) {
        return false;
    }

    m_nodes--;
    m_node.erase(m_node.begin()+id);
    m_activeNodes.erase(m_activeNodes.begin()+id);

    return true;
}

/*! Remove pre-existent set of nodes. RBF nodal list is resized and renumbered after extraction.
 *  Supported in both modes.
 * @param[in] list id list of candidates to extraction
 * @return boolean, true if all nodes are successfully extracted, false if any of them or none are extracted
 */
bool RBF::removeNode(std::vector<int> & list)
{
    std::set<int> setList;
    for(auto && id : list) {
        setList.insert(id);
    }

    int extracted = 0;
    for(auto && id : setList) {
        if(id>=0 && id <m_nodes){;
            m_nodes--;
            int index = id-extracted;
            auto nodeItr = m_node.begin() + index;
            m_node.erase(nodeItr);
            auto activeNodeItr = m_activeNodes.begin() + index;
            m_activeNodes.erase(activeNodeItr);
            extracted++;
        }
    }

    return(extracted == (int)(list.size()));
}

/*!
 * Remove all nodes in RBF nodal list. Supported in both modes.
 */
void RBF::removeAllNodes()
{
    m_nodes = 0;
    m_node.clear();
    m_activeNodes.clear();
}

/*!
 * Calculate distances between two RBF nodes in the list of class availables
 * @param[in] i i-th node in the list
 * @param[in] j j-th node in the list
 */
double RBF::calcDist(int i, int j)
{
    return norm2(m_node[i]-m_node[j]);
}

/*!
 * Calculate distances between a 3D point in space and a RBF node in the list of class availables
 * @param[in] point std::array<double,3> coordinates of the point
 * @param[in] j j-th RBF node in the list
 */
double RBF::calcDist(const std::array<double,3>& point, int j)
{
    return norm2(point-m_node[j]);
}

// RBF NAMESPACE UTILITIES

/*!
 * Wendland C2 function
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::wendlandc2( const double &dist )
{
    if( dist > 1) {
        return 0.;
    } else{
        return std::pow(1.-dist,4)*(4.*dist+1.);
    }
}

/*!
 * Linear function
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::linear( const double &dist )
{
    if( dist > 1) {
        return 0.;
    } else{
        return (1-dist);
    }
}

/*!
 * Non compact gaussian function with 0.1 value at dist equal to 1
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::gauss90( const double &dist )
{
    double eps = std::pow(-1.0*std::log(0.1),0.5);

    return std::exp(-1.0*std::pow(dist*eps,2));
}

/*!
 * Non compact gaussian function with 0.05 value at dist equal to 1
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::gauss95( const double &dist )
{
    double eps = std::pow(-1.0*std::log(0.05),0.5);

    return std::exp(-1.0*std::pow(dist*eps,2));
}

/*!
 * Non compact gaussian function with 0.01 value at dist equal to 1
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::gauss99( const double &dist )
{
    double eps = std::pow(-1.0*std::log(0.01),0.5);

    return std::exp(-1.0*std::pow(dist*eps,2));
}

/*!
 * Polynomial function defined between 0,1. Preserve C1 continuity at dist=0, C0 continuity at dist=1.
 * At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::c1c0( const double &dist )
{
    if( dist > 1) {
        return 0.;
    } else{
        return (1.0-std::pow(dist,2));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C2 continuity at dist=0, C0 continuity at dist=1.
 * At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::c2c0( const double &dist )
{
    if( dist > 1) {
        return 0.;
    } else{
        return (1.0-std::pow(dist,3));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C0 continuity at dist=0, C1 continuity at dist=1.
 * At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::c0c1( const double &dist )
{
    if( dist > 1) {
        return 0.;
    } else{
        return (1.0- 2.0*dist + std::pow(dist,2));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C1 continuity at dist=0, C1 continuity at dist=1.
 * At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::c1c1( const double &dist )
{
    if( dist > 1) {
        return 0.;
    } else{
        return (1.0-3.0*std::pow(dist,2)+2.0*std::pow(dist,3));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C2 continuity at dist=0, C1 continuity at dist=1.
 * At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::c2c1( const double &dist )
{
    if( dist > 1) {
        return 0.;
    } else{
        return (1.0- 4.0*std::pow(dist,3) + 3.0*std::pow(dist,4));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C0 continuity at dist=0, C2 continuity at dist=1.
 * At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::c0c2( const double &dist )
{
    if( dist > 1) {
        return 0.;
    } else{
        return (1.0 -3.0*dist +3.0*std::pow(dist,2) - std::pow(dist,3));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C1 continuity at dist=0, C2 continuity at dist=1.
 * At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::c1c2( const double &dist )
{
    if( dist > 1) {
        return 0.;
    } else{
        return (1.0 -6.0*std::pow(dist,2) + 8.0*std::pow(dist,3) - 3.0*std::pow(dist,4));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C2 continuity at dist=0, C2 continuity at dist=1.
 * At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
double rbf::c2c2( const double &dist )
{
    if( dist > 1) {
        return 0.;
    } else{
        return (1.0 -10.0*std::pow(dist,3) +15.0*std::pow(dist,4) -6.0*std::pow(dist,5));
    }
}

}
