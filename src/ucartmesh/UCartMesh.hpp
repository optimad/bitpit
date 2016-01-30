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

# include <bitpit_operators.hpp>
# include <bitpit_IO.hpp>


class UCartMesh{

    // Members ============================================================== //
    private:
        int                                 status ;                            /**< indentifier for mesh status; is incresed each time mesh is moified */

        int                                 dim;                                /**< number of space dimensions*/
        int                                 nCells ;                            /**< number of cells in grid*/
        int                                 nNodes ;                            /**< number of nodes in grid*/
        int                                 CellsInIJPlane;                     /**< number of cells in the IJ plane*/
        int                                 NodesInIJPlane;                     /**< number of nodes in the IJ plane*/

        std::array<double,3>                B0;                                 /**< min point of axis aligned boundig box*/
        std::array<double,3>                B1;                                 /**< max point of axis aligned boundig box*/

        std::array<int,3>                   nc;                                 /**< number of cells in each direction; if 2D nc[2]=1 anyway */
        std::array<int,3>                   np;                                 /**< number of nodes in each direction; if 2D np[2]=1 anyway */
        std::array<double,3>                h;                                  /**< grid spacing in each direction*/

        std::vector<std::vector<double>>    center;                             /**< center coordinates in each direction */
        std::vector<std::vector<double>>    edge;                               /**< node coordinates in each direction */

        std::array<int,6>                   whichDirection ;                    /**< maps face indices [0...5] to space direction [0...2] */
        std::array<int,6>                   whichStep ;                         /**< maps face indices [0...5] to positive or negative steps */

        // Constructors, Destructor, assignment================================== //
    public:

        UCartMesh( ) ;
        UCartMesh( std::array<double,3> const &, std::array<double,3> const &, std::array<int,3> const &, int dim=3 ) ;
        UCartMesh( std::array<double,3> const &, std::array<double,3> const &, int const &, int const & ) ; 
        UCartMesh( std::array<double,3> const &, std::array<double,3> const &, int const &, int const &, int const & ) ; 

        ~UCartMesh( ) ;

        UCartMesh&              operator=( const UCartMesh & ) ;
        void                    setMesh( std::array<double,3> const &, std::array<double,3> const &, std::array<int,3> const &, int const &dim=3 ) ;
        void                    ClearMesh( ) ;

    private:
        void                    ResizeMesh( ) ;

    public:
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
        void                    Translate( std::array<double,3> const & ) ;

        void                    Scale( std::array<double,3> const & ) ; 
        void                    Scale( std::array<double,3> const &, std::array<double,3> const & ) ; 


        // Cell information ----------------------------------------------------- //
        std::array<int,3>       CellCartesianId( std::array<double,3> const & ) ;
        void                    CellCartesianId( std::array<double,3> const &, int &, int & ) ;
        void                    CellCartesianId( std::array<double,3> const &, int &, int &, int & ) ;

        std::array<int,3>       CellCartesianId( int const & ) ; 
        void                    CellCartesianId( int const &, int &, int & ) ; 
        void                    CellCartesianId( int const &, int &, int &, int & ) ; 

        int                     CellLinearId( std::array<double,3> const & ) ;
        int                     CellLinearId( std::array<int,3> const & ) ;
        int                     CellLinearId( int const &, int const &, int const & k=0 ) ;

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
        std::array<int,3>       NodeCartesianId( std::array<double,3> const & ) ;
        void                    NodeCartesianId( std::array<double,3> const &, int &, int & ) ;
        void                    NodeCartesianId( std::array<double,3> const &, int &, int &, int & ) ;

        std::array<int,3>       NodeCartesianId( int const & ) ;
        void                    NodeCartesianId( int const &, int &, int & ) ;
        void                    NodeCartesianId( int const &, int &, int &, int & ) ;

        int                     NodeLinearId( std::array<double,3> const & ) ;
        int                     NodeLinearId( std::array<int,3> const & ) ;
        int                     NodeLinearId( int const &, int const &, int const & k=0 ) ;

        std::array<double,3>    getNodeCoordinates( int ) ;
        std::array<double,3>    getNodeCoordinates( std::array<int,3> ) ;
        std::array<double,3>    getNodeCoordinates( int, int, int k=0 ) ;

        int                     getNodeNeighbour( int const &, int const & ) ;
        int                     getNodeNeighbour( int const &, int const &, int const & ) ;


        // Format conversion ---------------------------------------------------- //
        void                    Cart2Unstr( int &, int &, std::vector<std::array<double,3>> &, std::vector<std::vector<int>> &, std::vector<std::vector<std::vector<int>>> & ) ;


        // subsets  ------------------------------------------------------------- //
        std::vector<int>        CellSubSet( int const &, int const & ) ;
        std::vector<int>        CellSubSet( std::array<int,3> const &, std::array<int,3> const & ) ;
        std::vector<int>        CellSubSet( std::array<double,3> const &, std::array<double,3> const & ) ;

        std::vector<int>        NodeSubSet( int const &, int const & ) ;
        std::vector<int>        NodeSubSet( std::array<int,3> const &, std::array<int,3> const & ) ;
        std::vector<int>        NodeSubSet( std::array<double,3> const &, std::array<double,3> const & ) ;


        // Point in Grid -------------------------------------------------------- //
        bool                    PointInGrid( std::array<double,3> const & ) ;
        bool                    PointInGrid( std::array<double,3> const &, int &) ;
        bool                    PointInGrid( std::array<double,3> const &, std::array<int,3> &) ;
        bool                    PointInGrid( std::array<double,3> const &, int &, int &, int &) ;

        // Interpolation -------------------------------------------------------- //
        int                     linearCellInterpolation( std::array<double,3> &, std::vector<int> &, std::vector<double> & ) ;
        int                     linearNodeInterpolation( std::array<double,3> &, std::vector<int> &, std::vector<double> & ) ;

        void                    CellData2NodeData( std::vector<double> &, std::vector<double> & ) ;
        void                    NodeData2CellData( std::vector<double> &, std::vector<double> & ) ;

        // I/O methods --------------------------------------------------------- //                   
        void                    ExportVtr(std::string, std::string);
        void                    ExportVtr(std::string, std::string, std::string, bitpit::VTKLocation, std::vector<double> &);
        void                    ExportVtr(std::string, std::string, std::string, bitpit::VTKLocation, std::vector<int> &);
        void                    ExportVtr(std::string, std::string, std::string, bitpit::VTKLocation, std::vector<std::array<double,3>> &);

};


#endif

/* @} */
