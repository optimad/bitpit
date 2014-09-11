// ========================================================================== //
//                           - STL IO FUNCTIONS -                             //
//                                                                            //
// Input/Output function for STL format.                                      //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author      :   Alessandro Alaia                                           //
// Version     :   v3.0                                                       //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //
# ifndef __STL_IOFunct__HPP__
# define __STL_IOFunct__HPP__

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <cmath>
# include <array>
# include <vector>
# include <string>
# include <sstream>
# include <fstream>
# include <iostream>

// CC_Lib
# include "Operators.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;

// ========================================================================== //
// TYPES DEFINITIONS                                                          //
// ========================================================================== //

// boolean vectors
typedef vector< bool >                 bvector1D;
typedef vector< bvector1D >            bvector2D;
typedef vector< bvector2D >            bvector3D;
typedef vector< bvector3D >            bvector4D;

// characters vectors
typedef vector< char >                 cvector1D;
typedef vector< cvector1D >            cvector2D;
typedef vector< cvector2D >            cvector3D;
typedef vector< cvector3D >            cvector4D;

// integer vectors
typedef vector< int >                  ivector1D;
typedef vector< ivector1D >            ivector2D;
typedef vector< ivector2D >            ivector3D;
typedef vector< ivector3D >            ivector4D;

// double vectors
typedef vector< double >               dvector1D;
typedef vector< dvector1D >            dvector2D;
typedef vector< dvector2D >            dvector3D;
typedef vector< dvector3D >            dvector4D;

// string vectors
typedef vector< string >               svector1D;
typedef vector< svector1D >            svector2D;
typedef vector< svector2D >            svector3D;
typedef vector< svector3D >            svector4D;

// ========================================================================== //
// CLASSES                                                                    //
// ========================================================================== //
struct STL_data {
    int             n_solids;                                                 // number of stl solids
    svector1D       solid_names;                                              // solids names
    ivector1D       solid_facets;                                             // number of facet for each stl solid
};

class STL_obj {

    // Public members ======================================================= //
    public:

    // General info
    string          stl_name;                                                 // stl file name
    bool            stl_type;                                                 // flag for binary/ASCII stl file

    // Error flags
    unsigned int    err;                                                      // general error

    // stl content
    bvector2D       stl_errors;                                               // error flags for each stl solid
    STL_data        data;                                                     // stl data

    // Private members ====================================================== //
    private:
    ifstream        ifile_handle;                                             // input stream to stl file
    ofstream        ofile_handle;                                             // output stream to stl file

    // Constructors ========================================================= //
    public:
    STL_obj(                                                                  // Default constructor
        void                                                                  // (input) none
    );
    STL_obj(                                                                  // Custom constructor #1
        string          ,                                                     // stl file name
        bool                                                                  // flag for binary/ASCII stl file
    );

    // Destructors ========================================================== //
    // none

    // Public methods ======================================================= //
    public:
    void open(                                                                // open stream to stl file
        string                                                                // (input) stream mode (input/output/append)
    );
    void close(                                                               // close stream to from file
        string           a = ""                                               // (input) stream to be closed (input/output/append)
    );
    void display(                                                             // Display info about stl file
        ostream         &                                                     // (input/output) output stream
    );
    void clear(                                                               // Clear info and error flags
        void                                                                  // (input) none
    );
    void scan(                                                                // Scan stl file for infos
        void                                                                  // (input) none
    );
    void check(                                                               // Check data structure in stl file
        void                                                                  // (input) none
    );
    void load(                                                                // Load stl triangulation from stl file
        int                 &,                                                // (input/output) number of stl vertices
        int                 &,                                                // (input/output) number of stl facets
        dvector2D           &,                                                // (input/output) vertex coordinate list
        dvector2D           &,                                                // (input/output) unit normal to each triangle
        ivector2D           &                                                 // (input/output) triangle-vertex connectivity
    );
    template <typename ... T2>
    void load(                                                                // Load stl triangulation from stl file
        string               ,                                                // (input) solid name
        int                 &,                                                // (input/output) number of stl vertices
        int                 &,                                                // (input/output) number of stl facets
        dvector2D           &,                                                // (input/output) vertex coordinate list
        dvector2D           &,                                                // (input/output) unit normal to each triangle
        ivector2D           &,                                                // (input/output) triangle-vertex connectivity
        T2              & ...                                                 // (input/optional) other stl solids
    );
    template <typename ... T2>
    void save(                                                                // Export stl solid to stl file
        string               ,                                                // (input) solid name
        int                 &,                                                // (input/output) number of stl vertices
        int                 &,                                                // (input/output) number of stl facets
        dvector2D           &,                                                // (input) vertex coordinate list
        dvector2D           &,                                                // (input) unit normal to each triangle
        ivector2D           &,                                                // (input) triangle-vertex connectivity
        T2              & ...                                                 // (input/optional) other stl solids
    );
    template <typename ... T2>
    void append(                                                              // Append stl solid to stl file
        string               ,                                                // (input) solid name
        int                 &,                                                // (input/output) number of stl vertices
        int                 &,                                                // (input/output) number of stl facets
        dvector2D           &,                                                // (input) vertex coordinate list
        dvector2D           &,                                                // (input) unit normal to each triangle
        ivector2D           &,                                                // (input) triangle-vertex connectivity
        T2              & ...                                                 // (input/optional) other stl solids
    );

    // Private methods ====================================================== //
    private:
    void save(                                                                // Dummy function for self-recursive variadic template
        void                                                                  // (input) none
    );
    void load(                                                                // Dummy function for self recursive variadic template
        void                                                                  // (input) none
    );
};

// ========================================================================== //
// FUNCTION PROTOTYPES                                                        //
// ========================================================================== //

// Scanning functions ------------------------------------------------------- //
unsigned int Scan_STL_ASCII(                                                  // Scan ASCII stl file and returns infos
    ifstream                &,                                                // (input) input stream to stl file
    svector1D               &,                                                // (input/output) names of each stl solid
    ivector1D               &                                                 // (input/output) number of facet for each stl solid
);
unsigned int Scan_STL_bin(                                                    // Scan binary stl file and returns infos
    ifstream                &,                                                // (input) input stream to stl file
    svector1D               &,                                                // (input/output) names of each stl solid
    ivector1D               &                                                 // (input/output) number of facet for each stl solid
);
unsigned int Scan_STLsolid_ASCII(                                             // Scan stl solid from binary stl file
    ifstream                &,                                                // (input) input stream to stl file
    int                     &                                                 // (input/output) number of facets
);

// Check routines ----------------------------------------------------------- //
unsigned int Check_STL_ASCII(                                                 // Check data coeherency in ASCII stl file
    ifstream                &,                                                // (input) input stream to stl file
    bvector2D               &                                                 // (input/output) error map
);
unsigned int Check_STLsolid_ASCII(                                            // Check coherency of stl solid data in ASCII stl files
    ifstream                &,                                                // (input) input stream to stl file
    bvector1D               &                                                 // (input/output) error map
);
unsigned int Check_STLfacet_ASCII(                                            // Check coherency of stl facet data in ASCII stl files
    ifstream                &,                                                // (input) input stream to stl file
    bvector1D               &                                                 // (input/output) error map
);
unsigned int Check_STL_bin(                                                   // Check data coherency in binary stl file
    ifstream                &,                                                // (input) input stream to stl file
    bvector2D               &                                                 // (input/output) error map
);

// Input routines ----------------------------------------------------------- //
unsigned int Read_STLsolid_ASCII(                                             // Read stl solid from ASCII stl file
    ifstream                &,                                                // (input/output) input stream to stl file
    int                     &,                                                // (input/output) number of vertices
    int                     &,                                                // (input/output) number of triangles
    dvector2D               &,                                                // (input/output) vertex coordinate list
    dvector2D               &,                                                // (input/output) triangles unit normal
    ivector2D               &,                                                // (input/output) triangle-vertex connectivity
    string                   a = ""                                           // (input/optional) stl solid name
);
unsigned int Read_STLfacet_ASCII(                                             // Read stl facet from ASCII stl file
    ifstream                &,                                                // (input/output) input stream to stl file
    int                     &,                                                // (input/output) number of vertices
    int                     &,                                                // (input/output) number of triangles
    dvector2D               &,                                                // (input/output) vertex coordinate list
    dvector2D               &,                                                // (input/output) triangles unit normal
    ivector2D               &                                                 // (input/output) triangle-vertex connectivity
);
unsigned int Read_STL_ASCII(                                                  // Load stl trianuglation data from ASCII stl file
    ifstream                &,                                                // (input) input stream to stl file
    int                     &,                                                // (input/output) number of stl veritces
    int                     &,                                                // (input/output) number of stl facets
    dvector2D               &,                                                // (input/output) vertex coordinate list
    dvector2D               &,                                                // (input/output) triangles unit normals
    ivector2D               &                                                 // (input/output) triangle-vertex connectivity
);
unsigned int Read_STL_bin(                                                    // Read stl triangulation data from binary stl file
    ifstream                &,                                                // (input) input stream to stl file
    int                     &,                                                // (input/output) number of stl veritces
    int                     &,                                                // (input/output) number of stl facets
    dvector2D               &,                                                // (input/output) vertex coordinate list
    dvector2D               &,                                                // (input/output) triangles unit normals
    ivector2D               &                                                 // (input/output) triangle-vertex connectivity
);

// Output routines ---------------------------------------------------------- //
unsigned int Write_STLsolid_ASCII(                                            // Write stl triangulation in a ASCII stl file
    ofstream                &,                                                // (input) output stream
    int                     &,                                                // (input) number of triangulation vertices
    int                     &,                                                // (input) number of triangles
    dvector2D               &,                                                // (input) vertex coordinate list
    dvector2D               &,                                                // (input) triangles' unit normals
    ivector2D               &,                                                // (input) triangle-vertex connectivity
    string                   a = ""                                           // (input/optional) stl solid name
);
unsigned int Write_STLsolid_bin(                                              // Write stl triangulation in a binary stl file
    ofstream                &,                                                // (input) output stream
    int                     &,                                                // (input) number of triangulation vertices
    int                     &,                                                // (input) number of triangles
    dvector2D               &,                                                // (input) vertex coordinate list
    dvector2D               &,                                                // (input) triangles' unit normals
    ivector2D               &,                                                // (input) triangle-vertex connectivity
    string                   a = ""                                           // (input/optional) stl solid name
);

// =================================================================================== //
// TEMPLATES                                                                           //
// =================================================================================== //
# include "STL_IOFunct.tpp"

# endif