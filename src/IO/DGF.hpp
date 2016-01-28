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

namespace bitpit{

// ========================================================================== //
// DATA STRUCTURES AND CLASSES                                                //
// ========================================================================== //
struct DGFData {
    int                                 nV;                                   // number of mesh vertices
    int                                 nS;                                   // number of mesh simplex
    std::vector<int>                    nV_data;                              // number of vertex data for each vertex data block
    std::vector<std::string>            sV_data;                              // vertex data set names
    std::vector<int>                    nS_data;                              // number of simplex data for each simplex data block
    std::vector<std::string>            sS_data;                              // simplex data set names
};

class DGFObj {

    // Members ============================================================== //

    // Public members ------------------------------------------------------- //
    public:
    std::string                         dgf_name;                             // dgf file name
    unsigned int                        err;                                  // general error flag
    DGFData                             data;                                 // dgf file content

    // Private members ------------------------------------------------------ //
    private:
    std::ifstream                       ifile_handle;                         // Input stream to dgf file
    std::ofstream                       ofile_handle;                         // Outuput stream to dgf file
    std::vector<std::vector<int> >      dgf_error;                            // Error map

    // Constructor ========================================================== //
    public:
    DGFObj(                                                                  // Standard constructor for DGF object
        void                                                                  // (input) none
    );
    DGFObj(                                                                  // Custom constructor #1 for DGF object
        std::string                                                           // (input) dgf file name
    );

    // Destructor =========================================================== //
    public:
    // none

    // Methods ============================================================== //

    // Public methods ------------------------------------------------------- //
    public:
    void open(                                                                // Open input/output stream to dgf file
        std::string                                                           // (input) stream mode
    );
    void close(                                                               // Close input/output stream to dgf file
        std::string                              a = "inout"                  // (input) stream mode
    );
    void clear(                                                               // Reser members to default value
        void                                                                  // (input) none
    );
    void display(                                                             // Display dgf info
        std::ostream                            &                             // (input) output stream
    );
    void scan(                                                                // Scan dgf file for infos
        void                                                                  // (input) none
    );
    void check(                                                               // Check data structure in .dgf file
        void                                                                  // (input) none
    );
    void load(                                                                // Load mesh data from .dgf file
        int                                     &,                            // (input/output) number of mesh vertices
        int                                     &,                            // (input/output) number of mesh facets
        std::vector<std::vector<double> >       &,                            // (input/output) vertex coordinate list
        std::vector<std::vector<int> >          &                             // (input/output) simplex-vertex connectivity
    );
    void load(                                                                // Load mesh data from .dgf file
        int                                     &,                            // (input/output) number of mesh vertices
        int                                     &,                            // (input/output) number of mesh facets
        std::vector<std::array<double,3> >      &,                            // (input/output) vertex coordinate list
        std::vector<std::vector<int> >          &                             // (input/output) simplex-vertex connectivity
    );
    void save(                                                                // Save mesh data into a .dgf file
        int                                     &,                            // (input) number of mesh vertices
        int                                     &,                            // (input) number of mesh facets
        std::vector<std::vector<double> >       &,                            // (input) vertex coordinate list
        std::vector<std::vector<int> >          &                             // (input) simplex-vertex connectivity
    );
    void save(                                                                // Save mesh data into a .dgf file
        int                                     &,                            // (input) number of mesh vertices
        int                                     &,                            // (input) number of mesh facets
        std::vector<std::array<double,3> >      &,                            // (input) vertex coordinate list
        std::vector<std::vector<int> >          &                             // (input) simplex-vertex connectivity
    );

    template< typename T, typename ... T2 >
    void loadVData(                                                          // Load vertex data sets from dgf file
        std::string                              data_name,                   // (input) dataset name
        int                                     &n,                           // (input/output) number of data in the dataset
        std::vector< T >                        &data,                        // (input/output) loaded dataset
        T2                                      &... others                   // (input/optional) others datasets to be loaded
    );
    template< typename T, typename ... T2 >
    void loadSData(                                                          // Load simplex data sets from dgf file
        std::string                              data_name,                   // (input) dataset name
        int                                     &n,                           // (input/output) number of data in the dataset
        std::vector< T >                        &data,                        // (input/output) loaded dataset
        T2                                      &... others                   // (input/optional) others datasets to be loaded
    );
    template < typename T, typename ... T2 >
    void appendVData(                                                        // Append vertex data set to dgf file
        std::string                              data_name,                   // (input) dataset name
        int                                     &n,                           // (input) number of data in the dataset
        std::vector< T >                        &data,                        // (input) dataset to be exported
        T2                                      &... others                   // (input/optional) others datasets to be exported
    );
    template < typename T, typename ... T2 >
    void appendSData(                                                        // Append simplex data set to dgf file
        std::string                              data_name,                   // (input) dataset name
        int                                     &n,                           // (input) number of data in the dataset
        std::vector< T >                        &data,                        // (input) dataset to be exported
        T2                                      &... others                   // (input/optional) others datasets to be exported
    );

    // Private methods ------------------------------------------------------ //
    private:
    void loadVData(                                                          // Dummy function for recursive variadic-template "load_vdata"
        void                                                                  // (input) none
    );
    void loadSData(                                                          // Dummy function for recursive variadic-template "load_sdata"
        void                                                                  // (input) none
    );
    void appendVData(                                                        // Dummy function for recursive variadic template "append_vdata"
        void                                                                  // (input) none
    );
    void appendSData(                                                        // Dummy function for recursive variadic template "append_sdata"
        void                                                                  // (input) none
    );

};

/*!
 * @ingroup  DuneGridFormat
 * @brief Utility fuctions for DGF IO
 */
namespace DGFUtils{

// Scanning routines -------------------------------------------------------- //
unsigned int scanData(                                                   // Scan DGF data and returns info
    std::ifstream                               &,                            // (input) input stream to dgf file
    int                                         &                             // (input/output) number of data in the dataset
);
unsigned int scan(                                                        // Scan DGF file and returns infos
    std::ifstream                               &,                            // (input) input stream to dgf file
    int                                         &,                            // (input/output) number of mesh vertices
    int                                         &,                            // (input/output) number of simplicies
    std::vector<std::string>                    &,                            // (input/output) names of vertex datasets
    std::vector<std::string>                    &,                            // (input/output) names of simplex datasets
    std::vector<int>                            &,                            // (input/output) number of data in each vertex dataset
    std::vector<int>                            &                             // (input/outptu) number of data in each simplex dataset
);

// Check routines ----------------------------------------------------------- //
unsigned int checkData(                                                  // Check DGF data structure
    std::ifstream                               &,                            // (input) input stream to dgf file
    int                                         &                             // (input/output) error code
);
unsigned int check(                                                       // Check DGF
    std::ifstream                               &,                            // (input) input stream to dgf file
    std::vector<std::vector<int> >              &                             // (input/output) error code for each dataset
);

// Input routines ----------------------------------------------------------- //
template< typename T >
unsigned int readData(                                                   // Read DGF data
    std::ifstream                               &,                            // (input) input stream to dgf file
    int                                         &,                            // (input/output) number of data loaded
    std::vector< T >                            &                             // (input/output) loaded data
);
unsigned int readMesh(                                                   // Read mesh data from file.
    std::ifstream                               &,                            // (input) input stream to dgf file
    int                                         &,                            // (input/output) number of mesh vertices
    int                                         &,                            // (input/output) number of simplicies
    std::vector<std::vector<double> >           &,                            // (input/output) vertex coordinate list
    std::vector<std::vector<int> >              &                             // (input/output) simplex-vertex connectivity
);
unsigned int readMesh(                                                   // Read mesh data from file.
    std::ifstream                               &,                            // (input) input stream to dgf file
    int                                         &,                            // (input/output) number of mesh vertices
    int                                         &,                            // (input/output) number of simplicies
    std::vector<std::array<double,3> >          &,                            // (input/output) vertex coordinate list
    std::vector<std::vector<int> >              &                             // (input/output) simplex-vertex connectivity
);
template <typename T>
unsigned int readVertexData(                                             // Load dgf vertex data
    std::ifstream                               &,                            // (input) input stream to dgf file
    int                                         &,                            // (input/output) number of loaded data
    std::vector< T >                            &,                            // (input/output) loaded data
    std::string                                  a = ""                       // (input/optional) dataset name
);
template< typename T >
unsigned int readSimplexData(                                            // Load dgf simplex data
    std::ifstream                               &,                            // (input) input stream to dgf file
    int                                         &,                            // (input/output) number of loaded data
    std::vector< T >                            &,                            // (input/output) loaded data
    std::string                                  a = ""                       // (input/optional) dataset name
);

// Output routines ---------------------------------------------------------- //
template < typename T >
unsigned int writeData(                                                  // Export data set to dgf file
    std::ofstream                               &,                            // (input) output stream to dgf file
    int                                         &,                            // (input) number of data in the dataset
    std::vector< T >                            &                             // (input) dataset to be exported
);
unsigned int writeMesh(                                                  // Export mesh data into dgf file
    std::ofstream                               &,                            // (input) output stream to dgf file
    int                                         &,                            // (input) number of vertices
    int                                         &,                            // (input) number of simplicies
    std::vector<std::vector<double> >           &,                            // (input) vertex coordinate list
    std::vector<std::vector<int> >              &                             // (input) simplex-vertex connectivity
);
unsigned int writeMesh(                                                  // Export mesh data into dgf file
    std::ofstream                               &,                            // (input) output stream to dgf file
    int                                         &,                            // (input) number of vertices
    int                                         &,                            // (input) number of simplicies
    std::vector<std::array<double,3> >          &,                            // (input) vertex coordinate list
    std::vector<std::vector<int> >              &                             // (input) simplex-vertex connectivity
);

template < typename T >
unsigned int writeVertexData(                                            // Export vertex data to dgf file
    std::ofstream                               &,                            // (input) output stream to dgf file
    int                                         &,                            // (input) number of data in the dataset
    std::vector< T >                            &,                            // (input) vertex data set
    std::string                                  a = ""                       // (input/optional) data set name
);
template < typename T >
unsigned int writeSimplexData(                                           // Export simplex data to dgf file
    std::ofstream                               &,                            // (input) output stream to dgf file
    int                                         &,                            // (input) number of data in the dataset
    std::vector< T >                            &,                            // (input) simplex data set
    std::string                                  a = ""                       // (input/optional) data set name
);

}

// ========================================================================== //
// TEMPLATES                                                                  //
// ========================================================================== //
# include "DGF.tpp"
}


# endif
