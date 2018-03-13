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

/*---------------------------------------------------------------------------*\
*
* Input/Output function for STL format.
* Author:   Alessandro Alaia
* 
\*---------------------------------------------------------------------------*/

# ifndef __BITPIT_STL_HPP__
# define __BITPIT_STL_HPP__

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

namespace bitpit{
// ========================================================================== //
// CLASSES                                                                    //
// ========================================================================== //
/*!
 * @struct STLData
 * @ingroup STereoLithography
 * @brief structure holding meta-information
 */
struct STLData {
    int                                 n_solids;                             /**< number of stl solids */
    std::vector<std::string>            solid_names;                          /**< solids names */
    std::vector<int>                    solid_facets;                         /**< number of facet for each stl solid */
};

class STLObj {

    // Public enumerations ================================================== //
    public:

    enum FileFormat {
        FormatInvalid = -1,
        FormatASCII,
        FormatBinary
    };

    // Public members ======================================================= //
    public:

    // General info
    std::string                         stl_name;                             /**< stl file name */
    bool                                stl_type;                             /**< flag for binary/ASCII stl file */

    // Error flags
    unsigned int                        err;                                  /**< general error */

    // stl content
    std::vector<std::vector<bool> >     stl_errors;                           /**< error flags for each stl solid */
    STLData                             data;                                 /**< stl data */

    // Private members ====================================================== //
    private:
    std::ifstream                       ifile_handle;                         /**< input stream to stl file */
    std::ofstream                       ofile_handle;                         /**< output stream to stl file */

    // Constructors ========================================================= //
    public:
    STLObj(                                                                  // Default constructor
        void                                                                  // (input) none
    );
    STLObj(                                                                  // Custom constructor #1
        std::string                      ,                                    // stl file name
        bool                                                                  // flag for binary/ASCII stl file
    );
    STLObj(                                                                  // Custom constructor #1
        std::string                                                           // stl file name
    );

    // Destructors ========================================================== //
    // none

    // Public methods ======================================================= //
    public:
    FileFormat detectFileFormat(                                              // detect file format
        std::string                                                           // (input) stream mode (input/output/append)
    );
    void open(                                                                // open stream to stl file
        std::string                                                           // (input) stream mode (input/output/append)
    );
    void close(                                                               // close stream to from file
        std::string                             mode = ""                     // (input) stream to be closed (input/output/append)
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
        std::vector<std::array<int,3> >         &                             // (input/output) triangle-vertex connectivity
    );
    void loadSolid(                                                            // Load next solid stl triangulation from stl file
        int                                     &,                            // (input/output) number of stl vertices
        int                                     &,                            // (input/output) number of stl facets
        std::vector<std::vector<double> >       &,                            // (input/output) vertex coordinate list
        std::vector<std::vector<double> >       &,                            // (input/output) unit normal to each triangle
        std::vector<std::vector<int> >          &,
        std::string                             &name                         // (input/optional) stl solid name
    );
    void loadSolid(                                                            // Load next stl triangulation from stl file
        int                                     &,                            // (input/output) number of stl vertices
        int                                     &,                            // (input/output) number of stl facets
        std::vector<std::array<double,3> >      &,                            // (input/output) vertex coordinate list
        std::vector<std::array<double,3> >      &,                            // (input/output) unit normal to each triangle
        std::vector<std::array<int,3> >         &,
        std::string                             &name                         // (input/optional) stl solid name        
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
        std::vector<std::array<int,3> >         &,                            // (input/output) triangle-vertex connectivity
        T2                                      & ...                         // (input/optional) other stl solids
    );
    void saveSolid(                                                           // Export next stl solid to stl file (appending)
        std::string                              ,                            // (input) solid name
        int                                     &,                            // (input/output) number of stl vertices
        int                                     &,                            // (input/output) number of stl facets
        std::vector<std::vector<double> >       &,                            // (input) vertex coordinate list
        std::vector<std::vector<double> >       &,                            // (input) unit normal to each triangle
        std::vector<std::vector<int> >          &                             // (input) triangle-vertex connectivity
    );
    void saveSolid(                                                           // Export next stl solid to stl file (appending)
        std::string                              ,                            // (input) solid name
        int                                     &,                            // (input/output) number of stl vertices
        int                                     &,                            // (input/output) number of stl facets
        std::vector<std::array<double,3> >      &,                            // (input) vertex coordinate list
        std::vector<std::array<double,3> >      &,                            // (input) unit normal to each triangle
        std::vector<std::array<int,3> >         &                             // (input) triangle-vertex connectivity
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
        std::vector<std::array<int,3> >         &,                            // (input) triangle-vertex connectivity
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
        std::vector<std::array<int,3> >         &,                            // (input) triangle-vertex connectivity
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


/*!
 * @ingroup  STereoLithography
 * @brief Utility fuctions for STL IO
 */
namespace stl{

// Scanning functions ------------------------------------------------------- //
unsigned int scanASCII(                                                  // Scan ASCII stl file and returns infos
    std::ifstream                               &,                            // (input) input stream to stl file
    std::vector<std::string>                    &,                            // (input/output) names of each stl solid
    std::vector<int>                            &                             // (input/output) number of facet for each stl solid
);
unsigned int scanBINARY(                                                    // Scan binary stl file and returns infos
    std::ifstream                               &,                            // (input) input stream to stl file
    std::vector<std::string>                    &,                            // (input/output) names of each stl solid
    std::vector<int>                            &                             // (input/output) number of facet for each stl solid
);
unsigned int scanSolidASCII(                                             // Scan stl solid from binary stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    int                                         &                             // (input/output) number of facets
);

// Check routines ----------------------------------------------------------- //
unsigned int checkASCII(                                                 // Check data coeherency in ASCII stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    std::vector<std::vector<bool> >             &                             // (input/output) error map
);
unsigned int checkSolidASCII(                                            // Check coherency of stl solid data in ASCII stl files
    std::ifstream                               &,                            // (input) input stream to stl file
    std::vector<bool>                           &                             // (input/output) error map
);
unsigned int checkFacetASCII(                                            // Check coherency of stl facet data in ASCII stl files
    std::ifstream                               &,                            // (input) input stream to stl file
    std::vector<bool>                           &                             // (input/output) error map
);
unsigned int checkBINARY(                                                   // Check data coherency in binary stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    std::vector<std::vector<bool> >             &                             // (input/output) error map
);

// Input routines ----------------------------------------------------------- //
unsigned int readSolidASCII(                                             // Read stl solid from ASCII stl file
    std::ifstream                               &,                            // (input/output) input stream to stl file
    bool                                         ,                            // (input) wrap around when searching for the solid
    int                                         &,                            // (input/output) number of vertices
    int                                         &,                            // (input/output) number of triangles
    std::vector<std::vector<double> >           &,                            // (input/output) vertex coordinate list
    std::vector<std::vector<double> >           &,                            // (input/output) triangles unit normal
    std::vector<std::vector<int> >              &,                            // (input/output) triangle-vertex connectivity
    std::string                                 &name                         // (input/optional) stl solid name
);
unsigned int readSolidASCII(                                             // Read stl solid from ASCII stl file
    std::ifstream                               &,                            // (input/output) input stream to stl file
    bool                                         ,                            // (input) wrap around when searching for the solid
    int                                         &,                            // (input/output) number of vertices
    int                                         &,                            // (input/output) number of triangles
    std::vector<std::array<double,3> >          &,                            // (input/output) vertex coordinate list
    std::vector<std::array<double,3> >          &,                            // (input/output) triangles unit normal
    std::vector<std::array<int, 3> >            &,                            // (input/output) triangle-vertex connectivity
    std::string                                 &name                         // (input/optional) stl solid name
);
unsigned int readFacetASCII(                                             // Read stl facet from ASCII stl file
    std::ifstream                               &,                            // (input/output) input stream to stl file
    int                                         &,                            // (input/output) number of vertices
    int                                         &,                            // (input/output) number of triangles
    std::vector<std::vector<double> >           &,                            // (input/output) vertex coordinate list
    std::vector<std::vector<double> >           &,                            // (input/output) triangles unit normal
    std::vector<std::vector<int> >              &                             // (input/output) triangle-vertex connectivity
);
unsigned int readFacetASCII(                                             // Read stl facet from ASCII stl file
    std::ifstream                               &,                            // (input/output) input stream to stl file
    int                                         &,                            // (input/output) number of vertices
    int                                         &,                            // (input/output) number of triangles
    std::vector<std::array<double,3> >          &,                            // (input/output) vertex coordinate list
    std::vector<std::array<double,3> >          &,                            // (input/output) triangles unit normal
    std::vector<std::array<int,3> >             &                             // (input/output) triangle-vertex connectivity
);
unsigned int readASCII(                                                  // Load stl trianuglation data from ASCII stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    int                                         &,                            // (input/output) number of stl veritces
    int                                         &,                            // (input/output) number of stl facets
    std::vector<std::vector<double> >           &,                            // (input/output) vertex coordinate list
    std::vector<std::vector<double> >           &,                            // (input/output) triangles unit normals
    std::vector<std::vector<int> >              &                             // (input/output) triangle-vertex connectivity
);
unsigned int readASCII(                                                  // Load stl trianuglation data from ASCII stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    int                                         &,                            // (input/output) number of stl veritces
    int                                         &,                            // (input/output) number of stl facets
    std::vector<std::array<double,3> >          &,                            // (input/output) vertex coordinate list
    std::vector<std::array<double,3> >          &,                            // (input/output) triangles unit normals
    std::vector<std::array<int,3> >             &                             // (input/output) triangle-vertex connectivity
);
unsigned int readBINARY(                                                    // Read stl triangulation data from binary stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    int                                         &,                            // (input/output) number of stl veritces
    int                                         &,                            // (input/output) number of stl facets
    std::vector<std::vector<double> >           &,                            // (input/output) vertex coordinate list
    std::vector<std::vector<double> >           &,                            // (input/output) triangles unit normals
    std::vector<std::vector<int> >              &                             // (input/output) triangle-vertex connectivity
);
unsigned int readBINARY(                                                    // Read stl triangulation data from binary stl file
    std::ifstream                               &,                            // (input) input stream to stl file
    int                                         &,                            // (input/output) number of stl veritces
    int                                         &,                            // (input/output) number of stl facets
    std::vector<std::array<double,3> >          &,                            // (input/output) vertex coordinate list
    std::vector<std::array<double,3> >          &,                            // (input/output) triangles unit normals
    std::vector<std::array<int,3> >             &                             // (input/output) triangle-vertex connectivity
);

// Output routines ---------------------------------------------------------- //
unsigned int writeSolidASCII(                                            // Write stl triangulation in a ASCII stl file
    std::ofstream                               &,                            // (input) output stream
    int                                         &,                            // (input) number of triangulation vertices
    int                                         &,                            // (input) number of triangles
    std::vector<std::vector<double> >           &,                            // (input) vertex coordinate list
    std::vector<std::vector<double> >           &,                            // (input) triangles' unit normals
    std::vector<std::vector<int> >              &,                            // (input) triangle-vertex connectivity
    std::string                                  solid_name = ""              // (input/optional) stl solid name
);
unsigned int writeSolidASCII(                                            // Write stl triangulation in a ASCII stl file
    std::ofstream                               &,                            // (input) output stream
    int                                         &,                            // (input) number of triangulation vertices
    int                                         &,                            // (input) number of triangles
    std::vector<std::array<double,3> >          &,                            // (input) vertex coordinate list
    std::vector<std::array<double,3> >          &,                            // (input) triangles' unit normals
    std::vector<std::array<int,3> >             &,                            // (input) triangle-vertex connectivity
    std::string                                  solid_name = ""              // (input/optional) stl solid name
);
unsigned int writeSolidBINARY(                                              // Write stl triangulation in a binary stl file
    std::ofstream                               &,                            // (input) output stream
    int                                         &,                            // (input) number of triangulation vertices
    int                                         &,                            // (input) number of triangles
    std::vector<std::vector<double> >           &,                            // (input) vertex coordinate list
    std::vector<std::vector<double> >           &,                            // (input) triangles' unit normals
    std::vector<std::vector<int> >              &,                            // (input) triangle-vertex connectivity
    std::string                                  solid_name = ""              // (input/optional) stl solid name
);
unsigned int writeSolidBINARY(                                              // Write stl triangulation in a binary stl file
    std::ofstream                               &,                            // (input) output stream
    int                                         &,                            // (input) number of triangulation vertices
    int                                         &,                            // (input) number of triangles
    std::vector<std::array<double,3> >          &,                            // (input) vertex coordinate list
    std::vector<std::array<double,3> >          &,                            // (input) triangles' unit normals
    std::vector<std::array<int,3> >             &,                            // (input) triangle-vertex connectivity
    std::string                                 solid_name = ""               // (input/optional) stl solid name
);

}

// =================================================================================== //
// TEMPLATES                                                                           //
// =================================================================================== //
# include "STL.tpp"

}

# endif
