// ========================================================================== //
//                              DGF IO FUNCTIONS                              //
//                                                                            //
// I/O functions for .dgf file formats                                        //
// ========================================================================== //
// INTFO                                                                      //
// ========================================================================== //
// Author  : Alessandro Alaia                                                 //
// Version : v3.0                                                             //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //
# ifndef __DGF_IOFUNCT_HH_
# define __DGF_IOFUNCT_HH_

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //

// Standard Template Library
# include <vector>
# include <string>
# include <fstream>
# include <sstream>
# include <cstdarg>
# include <iostream>

// Classes
// none

// CC_lib
# include "Operators.hpp"

// ========================================================================== //
// NAMESPACE                                                                  //
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

typedef array<double,3>                darray3E;
typedef vector<darray3E>               dvecarr3E;

// string vectors
typedef vector< string >               svector1D;
typedef vector< svector1D >            svector2D;
typedef vector< svector2D >            svector3D;
typedef vector< svector3D >            svector4D;

// ========================================================================== //
// DATA STRUCTURES AND CLASSES                                                //
// ========================================================================== //
struct DGF_data {
    int             nV;                                                       // number of mesh vertices
    int             nS;                                                       // number of mesh simplex
    ivector1D       nV_data;                                                  // number of vertex data for each vertex data block
    svector1D       sV_data;                                                  // vertex data set names
    ivector1D       nS_data;                                                  // number of simplex data for each simplex data block
    svector1D       sS_data;                                                  // simplex data set names
};

class DGF_obj {

    // Members ============================================================== //

    // Public members ------------------------------------------------------- //
    public:
    string          dgf_name;                                                 // dgf file name
    unsigned int    err;                                                      // general error flag
    DGF_data        data;                                                     // dgf file content

    // Private members ------------------------------------------------------ //
    private:
    ifstream        ifile_handle;                                             // Input stream to dgf file
    ofstream        ofile_handle;                                             // Outuput stream to dgf file
    ivector2D       dgf_error;                                                // Error map

    // Constructor ========================================================== //
    public:
    DGF_obj(                                                                  // Standard constructor for DGF object
        void                                                                  // (input) none
    );
    DGF_obj(                                                                  // Custom constructor #1 for DGF object
        string                                                                // (input) dgf file name
    );

    // Destructor =========================================================== //
    public:
    // none

    // Methods ============================================================== //

    // Public methods ------------------------------------------------------- //
    public:
    void open(                                                                // Open input/output stream to dgf file
        string                                                                // (input) stream mode
    );
    void close(                                                               // Close input/output stream to dgf file
        string      a = "inout"                                               // (input) stream mode
    );
    void clear(                                                               // Reser members to default value
        void                                                                  // (input) none
    );
    void display(                                                             // Display dgf info
        ostream    &                                                          // (input) output stream
    );
    void scan(                                                                // Scan dgf file for infos
        void                                                                  // (input) none
    );
    void check(                                                               // Check data structure in .dgf file
        void                                                                  // (input) none
    );
    void load(                                                                // Load mesh data from .dgf file
        int        &,                                                         // (input/output) number of mesh vertices
        int        &,                                                         // (input/output) number of mesh facets
        dvector2D  &,                                                         // (input/output) vertex coordinate list
        ivector2D  &                                                          // (input/output) simplex-vertex connectivity
    );
    void load(                                                                // Load mesh data from .dgf file
        int        &,                                                         // (input/output) number of mesh vertices
        int        &,                                                         // (input/output) number of mesh facets
        dvecarr3E  &,                                                         // (input/output) vertex coordinate list
        ivector2D  &                                                          // (input/output) simplex-vertex connectivity
    );
    void save(                                                                // Save mesh data into a .dgf file
        int        &,                                                         // (input) number of mesh vertices
        int        &,                                                         // (input) number of mesh facets
        dvector2D  &,                                                         // (input) vertex coordinate list
        ivector2D  &                                                          // (input) simplex-vertex connectivity
    );
    void save(                                                                // Save mesh data into a .dgf file
        int        &,                                                         // (input) number of mesh vertices
        int        &,                                                         // (input) number of mesh facets
        dvecarr3E  &,                                                         // (input) vertex coordinate list
        ivector2D  &                                                          // (input) simplex-vertex connectivity
    );

    template< typename T, typename ... T2 >
    void load_vdata(                                                          // Load vertex data sets from dgf file
        string      data_name,                                                // (input) dataset name
        int        &n,                                                        // (input/output) number of data in the dataset
        vector< T >&data,                                                     // (input/output) loaded dataset
        T2     &... others                                                    // (input/optional) others datasets to be loaded
    );
    template< typename T, typename ... T2 >
    void load_sdata(                                                          // Load simplex data sets from dgf file
        string      data_name,                                                // (input) dataset name
        int        &n,                                                        // (input/output) number of data in the dataset
        vector< T >&data,                                                     // (input/output) loaded dataset
        T2     &... others                                                    // (input/optional) others datasets to be loaded
    );
    template < typename T, typename ... T2 >
    void append_vdata(                                                        // Append vertex data set to dgf file
        string      data_name,                                                // (input) dataset name
        int        &n,                                                        // (input) number of data in the dataset
        vector< T >&data,                                                     // (input) dataset to be exported
        T2     &... others                                                    // (input/optional) others datasets to be exported
    );
    template < typename T, typename ... T2 >
    void append_sdata(                                                        // Append simplex data set to dgf file
        string      data_name,                                                // (input) dataset name
        int        &n,                                                        // (input) number of data in the dataset
        vector< T >&data,                                                     // (input) dataset to be exported
        T2     &... others                                                    // (input/optional) others datasets to be exported
    );

    // Private methods ------------------------------------------------------ //
    private:
    void load_vdata(                                                          // Dummy function for recursive variadic-template "load_vdata"
        void                                                                  // (input) none
    );
    void load_sdata(                                                          // Dummy function for recursive variadic-template "load_sdata"
        void                                                                  // (input) none
    );
    void append_vdata(                                                        // Dummy function for recursive variadic template "append_vdata"
        void                                                                  // (input) none
    );
    void append_sdata(                                                        // Dummy function for recursive variadic template "append_sdata"
        void                                                                  // (input) none
    );

};

// ========================================================================== //
// FUNCTION PROTOTYPES                                                        //
// ========================================================================== //

// Scanning routines -------------------------------------------------------- //
unsigned int Scan_DGF_data(                                                   // Scan DGF data and returns info
    ifstream       &,                                                         // (input) input stream to dgf file
    int            &                                                          // (input/output) number of data in the dataset
);
unsigned int Scan_DGF(                                                        // Scan DGF file and returns infos
    ifstream       &,                                                         // (input) input stream to dgf file
    int            &,                                                         // (input/output) number of mesh vertices
    int            &,                                                         // (input/output) number of simplicies
    svector1D      &,                                                         // (input/output) names of vertex datasets
    svector1D      &,                                                         // (input/output) names of simplex datasets
    ivector1D      &,                                                         // (input/output) number of data in each vertex dataset
    ivector1D      &                                                          // (input/outptu) number of data in each simplex dataset
);

// Check routines ----------------------------------------------------------- //
unsigned int Check_DGF_data(                                                  // Check DGF data structure
    ifstream       &,                                                         // (input) input stream to dgf file
    int            &                                                          // (input/output) error code
);
unsigned int Check_DGF(                                                       // Check DGF
    ifstream       &,                                                         // (input) input stream to dgf file
    ivector2D      &                                                          // (input/output) error code for each dataset
);

// Input routines ----------------------------------------------------------- //
template< typename T >
unsigned int Read_DGF_data(                                                   // Read DGF data
    ifstream       &,                                                         // (input) input stream to dgf file
    int            &,                                                         // (input/output) number of data loaded
    vector< T >    &                                                          // (input/output) loaded data
);
unsigned int Read_DGF_mesh(                                                   // Read mesh data from file.
    ifstream       &,                                                         // (input) input stream to dgf file
    int            &,                                                         // (input/output) number of mesh vertices
    int            &,                                                         // (input/output) number of simplicies
    dvector2D      &,                                                         // (input/output) vertex coordinate list
    ivector2D      &                                                          // (input/output) simplex-vertex connectivity
);
unsigned int Read_DGF_mesh(                                                   // Read mesh data from file.
    ifstream       &,                                                         // (input) input stream to dgf file
    int            &,                                                         // (input/output) number of mesh vertices
    int            &,                                                         // (input/output) number of simplicies
    dvecarr3E      &,                                                         // (input/output) vertex coordinate list
    ivector2D      &                                                          // (input/output) simplex-vertex connectivity
);
template <typename T>
unsigned int Read_DGF_VERTEXDATA(                                             // Load dgf vertex data
    ifstream       &,                                                         // (input) input stream to dgf file
    int            &,                                                         // (input/output) number of loaded data
    vector< T >    &,                                                         // (input/output) loaded data
    string          a = ""                                                    // (input/optional) dataset name
);
template< typename T >
unsigned int Read_DGF_SIMPLEXDATA(                                            // Load dgf simplex data
    ifstream       &,                                                         // (input) input stream to dgf file
    int            &,                                                         // (input/output) number of loaded data
    vector< T >    &,                                                         // (input/output) loaded data
    string          a = ""                                                    // (input/optional) dataset name
);

// Output routines ---------------------------------------------------------- //
template < typename T >
unsigned int Write_DGF_data(                                                  // Export data set to dgf file
    ofstream       &,                                                         // (input) output stream to dgf file
    int            &,                                                         // (input) number of data in the dataset
    vector< T >    &                                                          // (input) dataset to be exported
);
unsigned int Write_DGF_mesh(                                                  // Export mesh data into dgf file
    ofstream       &,                                                         // (input) output stream to dgf file
    int            &,                                                         // (input) number of vertices
    int            &,                                                         // (input) number of simplicies
    dvector2D      &,                                                         // (input) vertex coordinate list
    ivector2D      &                                                          // (input) simplex-vertex connectivity
);
unsigned int Write_DGF_mesh(                                                  // Export mesh data into dgf file
    ofstream       &,                                                         // (input) output stream to dgf file
    int            &,                                                         // (input) number of vertices
    int            &,                                                         // (input) number of simplicies
    dvecarr3E      &,                                                         // (input) vertex coordinate list
    ivector2D      &                                                          // (input) simplex-vertex connectivity
);

template < typename T >
unsigned int Write_DGF_VERTEXDATA(                                            // Export vertex data to dgf file
    ofstream       &,                                                         // (input) output stream to dgf file
    int            &,                                                         // (input) number of data in the dataset
    vector< T >    &,                                                         // (input) vertex data set
    string          a = ""                                                    // (input/optional) data set name
);
template < typename T >
unsigned int Write_DGF_SIMPLEXDATA(                                           // Export simplex data to dgf file
    ofstream       &,                                                         // (input) output stream to dgf file
    int            &,                                                         // (input) number of data in the dataset
    vector< T >    &,                                                         // (input) simplex data set
    string          a = ""                                                    // (input/optional) data set name
);

// ========================================================================== //
// TEMPLATES                                                                  //
// ========================================================================== //
# include "DGF_IOFunct.tpp"

# endif
