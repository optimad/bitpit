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
 * \ingroup    UCartMesh
 * @{
 * \class      UCartMesh
 * \brief Uniform Cartesian Mesh with variable spacing in each direction
 *      
 * UCartMesh provides the methods to access both Node and Cell indices either through cartesian indices or through a linear index.
 * Cartesian to Linear and vice versa.
 * SubSet extraction
 * Bi/Tri-linear interpolation methods.
 *
 */




// INCLUDES                                                                   //
# ifndef __UCARTMESH_HPP__
# define __UCARTMESH_HPP__

# include <array>
# include <vector>
# include <string>

# include <VTK.hpp>

class UCartMesh{

    // Members ============================================================== //
    private:
        int                                 m_status ;                            /**< indentifier for mesh status; is incresed each time mesh is moified */

        int                                 m_dim;                                /**< number of space dimensions*/
        int                                 m_nCells ;                            /**< number of cells in grid*/
        int                                 m_nNodes ;                            /**< number of nodes in grid*/
        int                                 m_CellsInIJPlane;                     /**< number of cells in the IJ plane*/
        int                                 m_NodesInIJPlane;                     /**< number of nodes in the IJ plane*/

        std::array<double,3>                m_B0;                                 /**< min point of axis aligned boundig box*/
        std::array<double,3>                m_B1;                                 /**< max point of axis aligned boundig box*/

        std::array<int,3>                   m_nc;                                 /**< number of cells in each direction; if 2D nc[2]=1 anyway */
        std::array<int,3>                   m_np;                                 /**< number of nodes in each direction; if 2D np[2]=1 anyway */
        std::array<double,3>                m_h;                                  /**< grid spacing in each direction*/

        std::vector<std::vector<double>>    m_center;                             /**< center coordinates in each direction */
        std::vector<std::vector<double>>    m_edge;                               /**< node coordinates in each direction */

        std::array<int,6>                   m_whichDirection ;                    /**< maps face indices [0...5] to space direction [0...2] */
        std::array<int,6>                   m_whichStep ;                         /**< maps face indices [0...5] to positive or negative steps */

        // Constructors, Destructor, assignment================================== //
    public:

        UCartMesh( ) ;
        UCartMesh( std::array<double,3> const &, std::array<double,3> const &, std::array<int,3> const &, int dim=3 ) ;
        UCartMesh( std::array<double,3> const &, std::array<double,3> const &, int const &, int const & ) ; 
        UCartMesh( std::array<double,3> const &, std::array<double,3> const &, int const &, int const &, int const & ) ; 

        ~UCartMesh( ) ;

        UCartMesh&              operator=( const UCartMesh & ) ;
        void                    setMesh( std::array<double,3> const &, std::array<double,3> const &, std::array<int,3> const &, int const &dim=3 ) ;
        void                    clearMesh( ) ;

        int                     getNCells();
        int                     getNCells( int );

        int                     getNNodes();
        int                     getNNodes( int d);

        std::array<double,3>    getSpacing( ) ;
        double                  getSpacing( int ) ;

        int                     getDimension( ) ;
        void                    getBoundingBox( std::array<double,3> &, std::array<double,3> & ) ; 

        int                     getStatus() ;

        // Transformations ------------------------------------------------------ //
        void                    translate( std::array<double,3> const & ) ;

        void                    scale( std::array<double,3> const & ) ; 
        void                    scale( std::array<double,3> const &, std::array<double,3> const & ) ; 


        // Cell information ----------------------------------------------------- //
        std::array<int,3>       getCellCartesianId( std::array<double,3> const & ) ;
        void                    getCellCartesianId( std::array<double,3> const &, int &, int & ) ;
        void                    getCellCartesianId( std::array<double,3> const &, int &, int &, int & ) ;

        std::array<int,3>       getCellCartesianId( int const & ) ; 
        void                    getCellCartesianId( int const &, int &, int & ) ; 
        void                    getCellCartesianId( int const &, int &, int &, int & ) ; 

        int                     getCellLinearId( std::array<double,3> const & ) ;
        int                     getCellLinearId( std::array<int,3> const & ) ;
        int                     getCellLinearId( int const &, int const &, int const & k=0 ) ;

        std::array<double,3>    getCellCenter( int ) ;
        std::array<double,3>    getCellCenter( std::array<int,3> ) ;
        std::array<double,3>    getCellCenter( int, int, int k=0 ) ;

        void                    getCellBoundingBox( int const &, std::array<double,3> &, std::array<double,3> &) ;
        void                    getCellBoundingBox( std::array<int,3> const &, std::array<double,3> &, std::array<double,3> &) ;
        void                    getCellBoundingBox( int const &, int const &, std::array<double,3> &, std::array<double,3> &  ) ;
        void                    getCellBoundingBox( int const &, int const &, int const &, std::array<double,3> &, std::array<double,3> &  ) ;

        int                     getCellNeighbour( int const &, int const & ) ;
        int                     getCellNeighbour( int const &, int const &, int const & ) ;


        // Node information ----------------------------------------------------- //
        std::array<int,3>       getNodeCartesianId( std::array<double,3> const & ) ;
        void                    getNodeCartesianId( std::array<double,3> const &, int &, int & ) ;
        void                    getNodeCartesianId( std::array<double,3> const &, int &, int &, int & ) ;

        std::array<int,3>       getNodeCartesianId( int const & ) ;
        void                    getNodeCartesianId( int const &, int &, int & ) ;
        void                    getNodeCartesianId( int const &, int &, int &, int & ) ;

        int                     getNodeLinearId( std::array<double,3> const & ) ;
        int                     getNodeLinearId( std::array<int,3> const & ) ;
        int                     getNodeLinearId( int const &, int const &, int const & k=0 ) ;

        std::array<double,3>    getNodeCoordinates( int ) ;
        std::array<double,3>    getNodeCoordinates( std::array<int,3> ) ;
        std::array<double,3>    getNodeCoordinates( int, int, int k=0 ) ;

        int                     getNodeNeighbour( int const &, int const & ) ;
        int                     getNodeNeighbour( int const &, int const &, int const & ) ;


        // Format conversion ---------------------------------------------------- //
        void                    convertToUnstructured( int &, int &, std::vector<std::array<double,3>> &, std::vector<std::vector<int>> &, std::vector<std::vector<std::vector<int>>> & ) ;


        // subsets  ------------------------------------------------------------- //
        std::vector<int>        extractCellSubSet( int const &, int const & ) ;
        std::vector<int>        extractCellSubSet( std::array<int,3> const &, std::array<int,3> const & ) ;
        std::vector<int>        extractCellSubSet( std::array<double,3> const &, std::array<double,3> const & ) ;

        std::vector<int>        extractNodeSubSet( int const &, int const & ) ;
        std::vector<int>        extractNodeSubSet( std::array<int,3> const &, std::array<int,3> const & ) ;
        std::vector<int>        extractNodeSubSet( std::array<double,3> const &, std::array<double,3> const & ) ;


        // Point in Grid -------------------------------------------------------- //
        bool                    isPointInGrid( std::array<double,3> const & ) ;
        bool                    isPointInGrid( std::array<double,3> const &, int &) ;
        bool                    isPointInGrid( std::array<double,3> const &, std::array<int,3> &) ;
        bool                    isPointInGrid( std::array<double,3> const &, int &, int &, int &) ;

        // Interpolation -------------------------------------------------------- //
        int                     linearCellInterpolation( std::array<double,3> &, std::vector<int> &, std::vector<double> & ) ;
        int                     linearNodeInterpolation( std::array<double,3> &, std::vector<int> &, std::vector<double> & ) ;

        void                    convertCellDataToNodeData( std::vector<double> &, std::vector<double> & ) ;
        void                    convertNodeDataToCellData( std::vector<double> &, std::vector<double> & ) ;

        // I/O methods --------------------------------------------------------- //                   
        void                    exportVTR(std::string, std::string);
        void                    exportVTR(std::string, std::string, std::string, bitpit::VTKLocation, std::vector<double> &);
        void                    exportVTR(std::string, std::string, std::string, bitpit::VTKLocation, std::vector<int> &);
        void                    exportVTR(std::string, std::string, std::string, bitpit::VTKLocation, std::vector<std::array<double,3>> &);

};


#endif

/* @} */
