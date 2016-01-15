using namespace std;

// ========================================================================== //
//                           - STL IO FUNCTIONS -                             //
//                                                                            //
// Templated functions for STL format.                                        //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author      :   Alessandro Alaia                                           //
// Version     :   v3.0                                                       //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
// none

// ========================================================================== //
// TEMPLATE IMPLEMENTATIONS                                                   //
// ========================================================================== //

/*!
 * @ingroup STereoLithography
 * @{
 */

// -------------------------------------------------------------------------- //
/*!
    Load data for solid with specified label from stl ascii file.
    
    \param[in] sname label associated to the solid. If empty label is specified
    (i.e. sname = ""), data of the first solid encountered in the stl file is returned.
    If no solid is found with the specified label, no data are loaded from the stl file.
    \param[in,out] nV on input stores the current number of vertices hosted in V.
    On output stores the input values incremented by the number
    of vertices acquired from the stl file.
    \param[in,out] nT on input stores the number of facet->vertex connectivity
    entries stores in T. On output stores the input value incremented by the
    number of facets acquired from the stl file.
    \param[in,out] V vertex coordinates list. On output stores the coordinates of
    vertices vertices acquired from the stl file. New vertices are appended
    at the end of V.
    \param[in,out] N facet normals. On output stores the normal unit vector to
    each facet acquired from the stl file. New normals are appended
    at the end of N.
    \param[in,out] T facet->vertex connectivity. On output stores the facet->vertex
    connectivity entries for the facets acquired from the stl file. New connectivity entries
    are appended at the end of T.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STL_obj::load(
    std::string                                  sname,
    int                                         &nV,
    int                                         &nT,
    std::vector<std::vector<double> >           &V,
    std::vector<std::vector<double> >           &N,
    std::vector<std::vector<int> >              &T,
    T2                                          & ... others
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN INPUT STREAM                                                          //
// ========================================================================== //
open("in");
if (err == 1) {
    return;
}

// ========================================================================== //
// LOAD STL SOLID                                                             //
// ========================================================================== //

// Load solid data
if (stl_type) { Read_STL_bin(ifile_handle, nV, nT, V, N, T); }
else          { Read_STLsolid_ASCII(ifile_handle, nV, nT, V, N, T, sname); }

// ========================================================================== //
// ITERATIVELY READ STL SOLIDS                                                //
// ========================================================================== //
load(others ...);

return; }

// -------------------------------------------------------------------------- //
/*!
    Load data for solid with specified label from stl ascii file. Overloading of
    member function STL_obj::load() for container vector<array<double, 3> >
    
    \param[in] sname label associated to the solid. If empty label is specified
    (i.e. sname = ""), data of the first solid encountered in the stl file is returned.
    If no solid is found with the specified label, no data are loaded from the stl file.
    \param[in,out] nV on input stores the current number of vertices hosted in V.
    On output stores the input values incremented by the number
    of vertices acquired from the stl file.
    \param[in,out] nT on input stores the number of facet->vertex connectivity
    entries stores in T. On output stores the input value incremented by the
    number of facets acquired from the stl file.
    \param[in,out] V vertex coordinates list. On output stores the coordinates of
    vertices vertices acquired from the stl file. New vertices are appended
    at the end of V.
    \param[in,out] N facet normals. On output stores the normal unit vector to
    each facet acquired from the stl file. New normals are appended
    at the end of N.
    \param[in,out] T facet->vertex connectivity. On output stores the facet->vertex
    connectivity entries for the facets acquired from the stl file. New connectivity entries
    are appended at the end of T.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STL_obj::load(
    std::string                                  sname,
    int                                         &nV,
    int                                         &nT,
    std::vector<std::array<double,3> >          &V,
    std::vector<std::array<double,3> >          &N,
    std::vector<std::vector<int> >              &T,
    T2                                          & ... others
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN INPUT STREAM                                                          //
// ========================================================================== //
open("in");
if (err == 1) {
    return;
}

// ========================================================================== //
// LOAD STL SOLID                                                             //
// ========================================================================== //

// Load solid data
if (stl_type) { Read_STL_bin(ifile_handle, nV, nT, V, N, T); }
else          { Read_STLsolid_ASCII(ifile_handle, nV, nT, V, N, T, sname); }

// ========================================================================== //
// ITERATIVELY READ STL SOLIDS                                                //
// ========================================================================== //
load(others ...);

return; }

// -------------------------------------------------------------------------- //
/*!
    Save solid data to stl ascii file.
    
    \param[in] sname label associated to the solid.
    \param[in,out] nV number of solid vertices.
    \param[in,out] nT number of solid facets.
    \param[in,out] V vertex coordinates list.
    \param[in,out] N facet normals. 
    \param[in,out] T facet->vertex connectivity.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STL_obj::save(
    std::string                                  sname,
    int                                         &nV,
    int                                         &nT,
    std::vector<std::vector<double> >           &V,
    std::vector<std::vector<double> >           &N,
    std::vector<std::vector<int> >              &T,
    T2                                          & ... others
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN OUTPUT STREAM                                                         //
// ========================================================================== //
open("out");
if (err == 1) { return; }

// ========================================================================== //
// EXPORT STL SOLID                                                           //
// ========================================================================== //

// Export solid
if (stl_type) { Write_STLsolid_bin(ofile_handle, nV, nT, V, N, T, sname); }
else          { Write_STLsolid_ASCII(ofile_handle, nV, nT, V, N, T, sname); }

// Recursively call variadic template
save(others ...);

return; }

// -------------------------------------------------------------------------- //
/*!
    Save solid data to stl ascii file. Overloading of member function STL_obj::save()
    for container vector<array<double, 3> >
    
    \param[in] sname label associated to the solid.
    \param[in,out] nV number of solid vertices.
    \param[in,out] nT number of solid facets.
    \param[in,out] V vertex coordinates list.
    \param[in,out] N facet normals. 
    \param[in,out] T facet->vertex connectivity.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STL_obj::save(
    std::string                                  sname,
    int                                         &nV,
    int                                         &nT,
    std::vector<std::array<double,3> >          &V,
    std::vector<std::array<double,3> >          &N,
    std::vector<std::vector<int> >              &T,
    T2                                          & ... others
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN OUTPUT STREAM                                                         //
// ========================================================================== //
open("out");
if (err == 1) { return; }

// ========================================================================== //
// EXPORT STL SOLID                                                           //
// ========================================================================== //

// Export solid
if (stl_type) { Write_STLsolid_bin(ofile_handle, nV, nT, V, N, T, sname); }
else          { Write_STLsolid_ASCII(ofile_handle, nV, nT, V, N, T, sname); }

// Recursively call variadic template
save(others ...);

return; }

// -------------------------------------------------------------------------- //
/*!
    Append solid data at the end of an existing stl ascii file.
    
    \param[in] sname label associated to the solid.
    \param[in,out] nV number of solid vertices.
    \param[in,out] nT number of solid facets.
    \param[in,out] V vertex coordinates list.
    \param[in,out] N facet normals. 
    \param[in,out] T facet->vertex connectivity.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STL_obj::append(
    std::string                                  sname,
    int                                         &nV,
    int                                         &nT,
    std::vector<std::vector<double> >           &V,
    std::vector<std::vector<double> >           &N,
    std::vector<std::vector<int> >              &T,
    T2                                          & ... others
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN OUTPUT STREAM                                                         //
// ========================================================================== //
open("app");
if (err == 1) { return; }

// ========================================================================== //
// EXPORT STL SOLID                                                           //
// ========================================================================== //

// Export solid
if (stl_type) { Write_STLsolid_bin(ofile_handle, nV, nT, V, N, T, sname); }
else          { Write_STLsolid_ASCII(ofile_handle, nV, nT, V, N, T, sname); }

// Recursively call variadic template
save(others ...);

return; }

// -------------------------------------------------------------------------- //
/*!
    Append solid data at the end of an existing stl ascii file. Overloading of
    member function STL_obj::append() for container vector<array<double, 3> >.
    
    \param[in] sname label associated to the solid.
    \param[in,out] nV number of solid vertices.
    \param[in,out] nT number of solid facets.
    \param[in,out] V vertex coordinates list.
    \param[in,out] N facet normals. 
    \param[in,out] T facet->vertex connectivity.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STL_obj::append(
    std::string                                  sname,
    int                                         &nV,
    int                                         &nT,
    std::vector<std::array<double,3> >          &V,
    std::vector<std::array<double,3> >          &N,
    std::vector<std::vector<int> >              &T,
    T2                                          & ... others
) {

// ========================================================================== //
// VARIABLES DECLARATION                                                      //
// ========================================================================== //

// Local variables
// none

// Counters
// none

// ========================================================================== //
// OPEN OUTPUT STREAM                                                         //
// ========================================================================== //
open("app");
if (err == 1) { return; }

// ========================================================================== //
// EXPORT STL SOLID                                                           //
// ========================================================================== //

// Export solid
if (stl_type) { Write_STLsolid_bin(ofile_handle, nV, nT, V, N, T, sname); }
else          { Write_STLsolid_ASCII(ofile_handle, nV, nT, V, N, T, sname); }

// Recursively call variadic template
save(others ...);

return; }


/*!
 * @}
 */
