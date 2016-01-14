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
// CLASSES                                                                    //
// ========================================================================== //
struct STL_data {
    int                                 n_solids;                             // number of stl solids
    std::vector<std::string>            solid_names;                          // solids names
    std::vector<int>                    solid_facets;                         // number of facet for each stl solid
};

class STL_obj {

    // Public members ======================================================= //
    public:

    // General info
    std::string                         stl_name;                             // stl file name
    bool                                stl_type;                             // flag for binary/ASCII stl file

    // Error flags
    unsigned int                        err;                                  // general error

    // stl content
    std::vector<std::vector<bool> >     stl_errors;                           // error flags for each stl solid
    STL_data                            data;                                 // stl data

    // Private members ====================================================== //
    private:
    std::ifstream                       ifile_handle;                         // input stream to stl file
    std::ofstream                       ofile_handle;                         // output stream to stl file

    // Constructors ========================================================= //
    public:
    STL_obj(                                                                  // Default constructor
        void                                                                  // (input) none
    );
    STL_obj(                                                                  // Custom constructor #1
        std::string                      ,                                    // stl file name
        bool                                                                  // flag for binary/ASCII stl file
    );

    // Destructors ========================================================== //
    // none

    // Public methods ======================================================= //
    public:
    void open(                                                                // open stream to stl file
        std::string                                                           // (input) stream mode (input/output/append)
    );
    void close(                                                               // close stream to from file
        std::string                              a = ""                       // (input) stream to be closed (input/output/append)
    );
    void display(                                                             // Display info about stl file
        std::ostream                            &                             // (input/output) output stream
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
        int                                     &,                            // (input/output) number of stl vertices
        int                                     &,                            // (input/output) number of stl facets
        std::vector<std::vector<double> >       &,                            // (input/output) vertex coordinate list
        std::vector<std::vector<double> >       &,                            // (input/output) unit normal to each triangle
        std::vector<std::vector<int> >          &                             // (input/output) triangle-vertex connectivity
    );
    void load(                                                                // Load stl triangulation from stl file
        int                                     &,                            // (input/output) number of stl vertices
        int                                     &,                            // (input/output) number of stl facets
        std::vector<std::array<double,3> >      &,                            // (input/output) vertex coordinate list
        std::vector<std::array<double,3> >      &,                            // (input/output) unit normal to each triangle
        std::vector<std::vector<int> >          &                             // (input/output) triangle-vertex connectivity
    );

    template <typename ... T2>
    void load(                                                                // Load stl triangulation from stl file
        std::string                              ,                            // (input) solid name
        int                                     &,                            // (input/output) number of stl vertices
        int                                     &,                            // (input/output) number of stl facets
        std::vector<std::vector<double> >       &,                            // (input/output) vertex coordinate list
        std::vector<std::vector<double> >       &,                            // (input/output) unit normal to each triangle
        std::vector<std::vector<int> >          &,                            // (input/output) triangle-vertex connectivity
        T2                                      & ...                         // (input/optional) other stl solids
    );
    template <typename ... T2>
    void load(                                                                // Load stl triangulation from stl file
        std::string                              ,                            // (input) solid name
        int                                     &,                            // (input/output) number of stl vertices
        int                                     &,                            // (input/output) number of stl facets
        std::vector<std::array<double,3> >      &,                            // (input/output) vertex coordinate list
        std::vector<std::array<double,3> >      &,                            // (input/output) unit normal to each triangle
        std::vector<std::vector<int> >          &,                            // (input/output) triangle-vertex connectivity
        T2                                      & ...                         // (input/optional) other stl solids
    );
    template <typename ... T2>
    void save(                                                                // Export stl solid to stl file
        std::string                              ,                            // (input) solid name
        int                                     &,                            // (input/output) number of stl vertices
        int                                     &,                            // (input/output) number of stl facets
        std::vector<std::vector<double> >       &,                            // (input) vertex coordinate list
        std::vector<std::vector<double> >       &,                            // (input) unit normal to each triangle
        std::vector<std::vector<int> >          &,                            // (input) triangle-vertex connectivity
        T2                                      & ...                         // (input/optional) other stl solids
    );
    template <typename ... T2>
    void save(                                                                // Export stl solid to stl file
        std::string                              ,                            // (input) solid name
        int                                     &,                            // (input/output) number of stl vertices
        int                                     &,                            // (input/output) number of stl facets
        std::vector<std::array<double,3> >      &,                            // (input) vertex coordinate list
        std::vector<std::array<double,3> >      &,                            // (input) unit normal to each triangle
        std::vector<std::vector<int> >          &,                            // (input) triangle-vertex connectivity
        T2                                      & ...                         // (input/optional) other stl solids
    );
    template <typename ... T2>
    void append(                                                              // Append stl solid to stl file
        std::string                              ,                            // (input) solid name
        int                                     &,                            // (input/output) number of stl vertices
        int                                     &,                            // (input/output) number of stl facets
        std::vector<std::vector<double> >       &,                            // (input) vertex coordinate list
        std::vector<std::vector<double> >       &,                            // (input) unit normal to each triangle
        std::vector<std::vector<int> >          &,                            // (input) triangle-vertex connectivity
        T2                                      & ...                         // (input/optional) other stl solids
    );
    template <typename ... T2>
    void append(                                                              // Append stl solid to stl file
        std::string                              ,                            // (input) solid name
        int                                     &,                            // (input/output) number of stl vertices
        int                                     &,                            // (input/output) number of stl facets
        std::vector<std::array<double,3> >      &,                            // (input) vertex coordinate list
        std::vector<std::array<double,3> >      &,                            // (input) unit normal to each triangle
        std::vector<std::vector<int> >          &,                            // (input) triangle-vertex connectivity
        T2                                      & ...                         // (input/optional) other stl solids
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
    std::ifstream                               &,                            // (input) input stream to stl file
    std::vector<std::string>                    &,                            // (input/output) names of each stl solid
    std::vector<int>                            &                             // (input/output) number of facet for each stl solid
);
unsigned int Scan_STL_bin(                                                    // Scan binary stl file and returns infos
    std::ifstream                               &,                            // (input) input stream to stl file
    std::vector<std::string>                    &,                            // (input/output) names of each stl solid
    std::vector<int>                            &                             // (input/output) number of facet for each stl solid
);
unsigned int Scan_STLsolid_ASCII(                                             // Scan stl solid from binary stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    int                                         &                             // (input/output) number of facets
);

// Check routines ----------------------------------------------------------- //
unsigned int Check_STL_ASCII(                                                 // Check data coeherency in ASCII stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    std::vector<std::vector<bool> >             &                             // (input/output) error map
);
unsigned int Check_STLsolid_ASCII(                                            // Check coherency of stl solid data in ASCII stl files
    std::ifstream                               &,                            // (input) input stream to stl file
    std::vector<bool>                           &                             // (input/output) error map
);
unsigned int Check_STLfacet_ASCII(                                            // Check coherency of stl facet data in ASCII stl files
    std::ifstream                               &,                            // (input) input stream to stl file
    std::vector<bool>                           &                             // (input/output) error map
);
unsigned int Check_STL_bin(                                                   // Check data coherency in binary stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    std::vector<std::vector<bool> >             &                             // (input/output) error map
);

// Input routines ----------------------------------------------------------- //
unsigned int Read_STLsolid_ASCII(                                             // Read stl solid from ASCII stl file
    std::ifstream                               &,                            // (input/output) input stream to stl file
    int                                         &,                            // (input/output) number of vertices
    int                                         &,                            // (input/output) number of triangles
    std::vector<std::vector<double> >           &,                            // (input/output) vertex coordinate list
    std::vector<std::vector<double> >           &,                            // (input/output) triangles unit normal
    std::vector<std::vector<int> >              &,                            // (input/output) triangle-vertex connectivity
    std::string                                  a = ""                       // (input/optional) stl solid name
);
unsigned int Read_STLsolid_ASCII(                                             // Read stl solid from ASCII stl file
    std::ifstream                               &,                            // (input/output) input stream to stl file
    int                                         &,                            // (input/output) number of vertices
    int                                         &,                            // (input/output) number of triangles
    std::vector<std::array<double,3> >          &,                            // (input/output) vertex coordinate list
    std::vector<std::array<double,3> >          &,                            // (input/output) triangles unit normal
    std::vector<std::vector<int> >              &,                            // (input/output) triangle-vertex connectivity
    std::string                                  a = ""                       // (input/optional) stl solid name
);
unsigned int Read_STLfacet_ASCII(                                             // Read stl facet from ASCII stl file
    std::ifstream                               &,                            // (input/output) input stream to stl file
    int                                         &,                            // (input/output) number of vertices
    int                                         &,                            // (input/output) number of triangles
    std::vector<std::vector<double> >           &,                            // (input/output) vertex coordinate list
    std::vector<std::vector<double> >           &,                            // (input/output) triangles unit normal
    std::vector<std::vector<int> >              &                             // (input/output) triangle-vertex connectivity
);
unsigned int Read_STLfacet_ASCII(                                             // Read stl facet from ASCII stl file
    std::ifstream                               &,                            // (input/output) input stream to stl file
    int                                         &,                            // (input/output) number of vertices
    int                                         &,                            // (input/output) number of triangles
    std::vector<std::array<double,3> >          &,                            // (input/output) vertex coordinate list
    std::vector<std::array<double,3> >          &,                            // (input/output) triangles unit normal
    std::vector<std::vector<int> >              &                             // (input/output) triangle-vertex connectivity
);
unsigned int Read_STL_ASCII(                                                  // Load stl trianuglation data from ASCII stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    int                                         &,                            // (input/output) number of stl veritces
    int                                         &,                            // (input/output) number of stl facets
    std::vector<std::vector<double> >           &,                            // (input/output) vertex coordinate list
    std::vector<std::vector<double> >           &,                            // (input/output) triangles unit normals
    std::vector<std::vector<int> >              &                             // (input/output) triangle-vertex connectivity
);
unsigned int Read_STL_ASCII(                                                  // Load stl trianuglation data from ASCII stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    int                                         &,                            // (input/output) number of stl veritces
    int                                         &,                            // (input/output) number of stl facets
    std::vector<std::array<double,3> >          &,                            // (input/output) vertex coordinate list
    std::vector<std::array<double,3> >          &,                            // (input/output) triangles unit normals
    std::vector<std::vector<int> >              &                             // (input/output) triangle-vertex connectivity
);
unsigned int Read_STL_bin(                                                    // Read stl triangulation data from binary stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    int                                         &,                            // (input/output) number of stl veritces
    int                                         &,                            // (input/output) number of stl facets
    std::vector<std::vector<double> >           &,                            // (input/output) vertex coordinate list
    std::vector<std::vector<double> >           &,                            // (input/output) triangles unit normals
    std::vector<std::vector<int> >              &                             // (input/output) triangle-vertex connectivity
);
unsigned int Read_STL_bin(                                                    // Read stl triangulation data from binary stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    int                                         &,                            // (input/output) number of stl veritces
    int                                         &,                            // (input/output) number of stl facets
    std::vector<std::array<double,3> >          &,                            // (input/output) vertex coordinate list
    std::vector<std::array<double,3> >          &,                            // (input/output) triangles unit normals
    std::vector<std::vector<int> >              &                             // (input/output) triangle-vertex connectivity
);

// Output routines ---------------------------------------------------------- //
unsigned int Write_STLsolid_ASCII(                                            // Write stl triangulation in a ASCII stl file
    std::ofstream                               &,                            // (input) output stream
    int                                         &,                            // (input) number of triangulation vertices
    int                                         &,                            // (input) number of triangles
    std::vector<std::vector<double> >           &,                            // (input) vertex coordinate list
    std::vector<std::vector<double> >           &,                            // (input) triangles' unit normals
    std::vector<std::vector<int> >              &,                            // (input) triangle-vertex connectivity
    std::string                                  a = ""                       // (input/optional) stl solid name
);
unsigned int Write_STLsolid_ASCII(                                            // Write stl triangulation in a ASCII stl file
    std::ofstream                               &,                            // (input) output stream
    int                                         &,                            // (input) number of triangulation vertices
    int                                         &,                            // (input) number of triangles
    std::vector<std::array<double,3> >          &,                            // (input) vertex coordinate list
    std::vector<std::array<double,3> >          &,                            // (input) triangles' unit normals
    std::vector<std::vector<int> >              &,                            // (input) triangle-vertex connectivity
    std::string                                  a = ""                       // (input/optional) stl solid name
);
unsigned int Write_STLsolid_bin(                                              // Write stl triangulation in a binary stl file
    std::ofstream                               &,                            // (input) output stream
    int                                         &,                            // (input) number of triangulation vertices
    int                                         &,                            // (input) number of triangles
    std::vector<std::vector<double> >           &,                            // (input) vertex coordinate list
    std::vector<std::vector<double> >           &,                            // (input) triangles' unit normals
    std::vector<std::vector<int> >              &,                            // (input) triangle-vertex connectivity
    std::string                                  a = ""                       // (input/optional) stl solid name
);
unsigned int Write_STLsolid_bin(                                              // Write stl triangulation in a binary stl file
    std::ofstream                               &,                            // (input) output stream
    int                                         &,                            // (input) number of triangulation vertices
    int                                         &,                            // (input) number of triangles
    std::vector<std::array<double,3> >          &,                            // (input) vertex coordinate list
    std::vector<std::array<double,3> >          &,                            // (input) triangles' unit normals
    std::vector<std::vector<int> >              &,                            // (input) triangle-vertex connectivity
    std::string                                  a = ""                       // (input/optional) stl solid name
);

// =================================================================================== //
// TEMPLATES                                                                           //
// =================================================================================== //
# include "STL_IOFunct.tpp"

# endif
