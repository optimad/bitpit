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


// NAMESPACES                                                                 //

// TYPES DEFINITIONS                                                          //

// boolean vectors
typedef std::vector< bool >                 bvector1D;
typedef std::vector< bvector1D >            bvector2D;
typedef std::vector< bvector2D >            bvector3D;
typedef std::vector< bvector3D >            bvector4D;

// characters vectors
typedef std::vector< char >                 cvector1D;
typedef std::vector< cvector1D >            cvector2D;
typedef std::vector< cvector2D >            cvector3D;
typedef std::vector< cvector3D >            cvector4D;

// integer vectors
typedef std::vector< int >                  ivector1D;
typedef std::vector< ivector1D >            ivector2D;
typedef std::vector< ivector2D >            ivector3D;
typedef std::vector< ivector3D >            ivector4D;

typedef std::array<int,3>                   iarray3E;

// double vectors
typedef std::vector< double >               dvector1D;
typedef std::vector< dvector1D >            dvector2D;
typedef std::vector< dvector2D >            dvector3D;
typedef std::vector< dvector3D >            dvector4D;

typedef std::array<double,3>                darray3E;

// string vectors
typedef std::vector< std::string >          svector1D;
typedef std::vector< svector1D >            svector2D;
typedef std::vector< svector2D >            svector3D;
typedef std::vector< svector3D >            svector4D;

class UCartMesh{

    // Members ============================================================== //
    private:

        int                 dim;                                /**< number of space dimensions*/
        int                 nCells ;                            /**< number of cells in grid*/
        int                 nNodes ;                            /**< number of nodes in grid*/
        int                 CellsInIJPlane;                     /**< number of cells in the IJ plane*/
        int                 NodesInIJPlane;                     /**< number of nodes in the IJ plane*/

        darray3E            B0;                                 /**< min point of axis aligned boundig box*/
        darray3E            B1;                                 /**< max point of axis aligned boundig box*/

        iarray3E            nc;                                 /**< number of cells in each direction; if 2D nc[2]=1 anyway */
        iarray3E            np;                                 /**< number of nodes in each direction; if 2D np[2]=1 anyway */
        darray3E            h;                                  /**< grid spacing in each direction*/

        dvector2D           center;                             /**< center coordinates in each direction */
        dvector2D           edge;                               /**< node coordinates in each direction */

        std::array<int,6>   whichDirection ;                    /**< maps face indices [0...5] to space direction [0...2] */
        std::array<int,6>   whichStep ;                         /**< maps face indices [0...5] to positive or negative steps */

        int                 status ;                            /**< indentifier for mesh status; is incresed each time mesh is moified */

        // Constructors, Destructor, assignment================================== //
    public:

        UCartMesh( ) ;
        UCartMesh( darray3E const &, darray3E const &, iarray3E const &, int dim=3 ) ;
        UCartMesh( darray3E const &, darray3E const &, int const &, int const & ) ; 
        UCartMesh( darray3E const &, darray3E const &, int const &, int const &, int const & ) ; 

        ~UCartMesh( ) ;

        UCartMesh&  operator=( const UCartMesh & ) ;
        void        setMesh( darray3E const &, darray3E const &, iarray3E const &, int const &dim=3 ) ;
        void        ClearMesh( ) ;

    private:
        void        ResizeMesh( ) ;

    public:
        int         getNCells();
        int         getNCells( int );

        int         getNNodes();
        int         getNNodes( int d);

        darray3E    getSpacing( ) ;
        double      getSpacing( int ) ;

        int         getDimension( ) ;
        void        getBoundingBox( darray3E &, darray3E & ) ; 

        int         getStatus() ;

        // Transformations ------------------------------------------------------ //
        void        Translate( darray3E const & ) ;

        void        Scale( darray3E const & ) ; 
        void        Scale( darray3E const &, darray3E const & ) ; 


        // Cell information ----------------------------------------------------- //
        iarray3E    CellCartesianId( darray3E const & ) ;
        void        CellCartesianId( darray3E const &, int &, int & ) ;
        void        CellCartesianId( darray3E const &, int &, int &, int & ) ;

        iarray3E    CellCartesianId( int const & ) ; 
        void        CellCartesianId( int const &, int &, int & ) ; 
        void        CellCartesianId( int const &, int &, int &, int & ) ; 

        int         CellLinearId( darray3E const & ) ;
        int         CellLinearId( iarray3E const & ) ;
        int         CellLinearId( int const &, int const &, int const & k=0 ) ;

        darray3E    getCellCenter( int ) ;
        darray3E    getCellCenter( iarray3E ) ;
        darray3E    getCellCenter( int, int, int k=0 ) ;

        void        getCellBoundingBox( int const &, darray3E &, darray3E &) ;
        void        getCellBoundingBox( iarray3E const &, darray3E &, darray3E &) ;
        void        getCellBoundingBox( int const &, int const &, darray3E &, darray3E &  ) ;
        void        getCellBoundingBox( int const &, int const &, int const &, darray3E &, darray3E &  ) ;

        int         getCellNeighbour( int const &, int const & ) ;
        int         getCellNeighbour( int const &, int const &, int const & ) ;


        // Node information ----------------------------------------------------- //
        iarray3E    NodeCartesianId( darray3E const & ) ;
        void        NodeCartesianId( darray3E const &, int &, int & ) ;
        void        NodeCartesianId( darray3E const &, int &, int &, int & ) ;

        iarray3E    NodeCartesianId( int const & ) ;
        void        NodeCartesianId( int const &, int &, int & ) ;
        void        NodeCartesianId( int const &, int &, int &, int & ) ;

        int         NodeLinearId( darray3E const & ) ;
        int         NodeLinearId( iarray3E const & ) ;
        int         NodeLinearId( int const &, int const &, int const & k=0 ) ;

        darray3E    getNodeCoordinates( int ) ;
        darray3E    getNodeCoordinates( iarray3E ) ;
        darray3E    getNodeCoordinates( int, int, int k=0 ) ;

        int         getNodeNeighbour( int const &, int const & ) ;
        int         getNodeNeighbour( int const &, int const &, int const & ) ;


        // Format conversion ---------------------------------------------------- //
        void        Cart2Unstr( int &, int &, std::vector<darray3E> &, ivector2D &, ivector3D & ) ;


        // subsets  ------------------------------------------------------------- //
        ivector1D   CellSubSet( int const &, int const & ) ;
        ivector1D   CellSubSet( iarray3E const &, iarray3E const & ) ;
        ivector1D   CellSubSet( darray3E const &, darray3E const & ) ;

        ivector1D   NodeSubSet( int const &, int const & ) ;
        ivector1D   NodeSubSet( iarray3E const &, iarray3E const & ) ;
        ivector1D   NodeSubSet( darray3E const &, darray3E const & ) ;


        // Point in Grid -------------------------------------------------------- //
        bool        PointInGrid( darray3E const & ) ;
        bool        PointInGrid( darray3E const &, int &) ;
        bool        PointInGrid( darray3E const &, iarray3E &) ;
        bool        PointInGrid( darray3E const &, int &, int &, int &) ;

        // Interpolation -------------------------------------------------------- //
        int         linearCellInterpolation( darray3E &, ivector1D &, dvector1D & ) ;
        int         linearNodeInterpolation( darray3E &, ivector1D &, dvector1D & ) ;

        void        CellData2NodeData( dvector1D &, dvector1D & ) ;
        void        NodeData2CellData( dvector1D &, dvector1D & ) ;

        // I/O methods --------------------------------------------------------- //                   
        void        ExportVtr(std::string, std::string);
        void        ExportVtr(std::string, std::string, std::string, std::string, std::vector<double> &);
        void        ExportVtr(std::string, std::string, std::string, std::string, std::vector<int> &);
        void        ExportVtr(std::string, std::string, std::string, std::string, std::vector<std::array<double,3>> &);

};


#endif

/* @} */
