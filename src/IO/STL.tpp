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
void STLObj::load(
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
if (stl_type) { stl::readBINARY(ifile_handle, nV, nT, V, N, T); }
else          { stl::readSolidASCII(ifile_handle, true, nV, nT, V, N, T, sname); }

// ========================================================================== //
// ITERATIVELY READ STL SOLIDS                                                //
// ========================================================================== //
load(others ...);

return; }

// -------------------------------------------------------------------------- //
/*!
    Load data for solid with specified label from stl ascii file. Overloading of
    member function STLObj::load() for container vector<array<double, 3> >
    
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
void STLObj::load(
    std::string                                  sname,
    int                                         &nV,
    int                                         &nT,
    std::vector<std::array<double,3> >          &V,
    std::vector<std::array<double,3> >          &N,
    std::vector<std::array<int,3> >             &T,
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
if (stl_type) { stl::readBINARY(ifile_handle, nV, nT, V, N, T); }
else          { stl::readSolidASCII(ifile_handle, true, nV, nT, V, N, T, sname); }

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
void STLObj::save(
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
if (stl_type) { stl::writeSolidBINARY(ofile_handle, nV, nT, V, N, T, sname); }
else          { stl::writeSolidASCII(ofile_handle, nV, nT, V, N, T, sname); }

// Recursively call variadic template
save(others ...);

return; }

// -------------------------------------------------------------------------- //
/*!
    Save solid data to stl ascii file. Overloading of member function STLObj::save()
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
void STLObj::save(
    std::string                                  sname,
    int                                         &nV,
    int                                         &nT,
    std::vector<std::array<double,3> >          &V,
    std::vector<std::array<double,3> >          &N,
    std::vector<std::array<int,3> >             &T,
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
if (stl_type) { stl::writeSolidBINARY(ofile_handle, nV, nT, V, N, T, sname); }
else          { stl::writeSolidASCII(ofile_handle, nV, nT, V, N, T, sname); }

// Recursively call variadic template
save(others ...);

return; }

// -------------------------------------------------------------------------- //
/*!
    Append solid data at the end of an existing stl ascii file.
   
    \tparam T2 vriadic template 
    \param[in] sname label associated to the solid.
    \param[in,out] nV number of solid vertices.
    \param[in,out] nT number of solid facets.
    \param[in,out] V vertex coordinates list.
    \param[in,out] N facet normals. 
    \param[in,out] T facet->vertex connectivity.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STLObj::append(
    std::string                                 sname,
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
if (stl_type) { stl::writeSolidBINARY(ofile_handle, nV, nT, V, N, T, sname); }
else          { stl::writeSolidASCII(ofile_handle, nV, nT, V, N, T, sname); }

// Recursively call variadic template
save(others ...);

return; }

// -------------------------------------------------------------------------- //
/*!
    Append solid data at the end of an existing stl ascii file. Overloading of
    member function STLObj::append() for container vector<array<double, 3> >.
    
    \param[in] sname label associated to the solid.
    \param[in,out] nV number of solid vertices.
    \param[in,out] nT number of solid facets.
    \param[in,out] V vertex coordinates list.
    \param[in,out] N facet normals. 
    \param[in,out] T facet->vertex connectivity.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STLObj::append(
    std::string                                  sname,
    int                                         &nV,
    int                                         &nT,
    std::vector<std::array<double,3> >          &V,
    std::vector<std::array<double,3> >          &N,
    std::vector<std::array<int,3> >             &T,
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
if (stl_type) { stl::writeSolidBINARY(ofile_handle, nV, nT, V, N, T, sname); }
else          { stl::writeSolidASCII(ofile_handle, nV, nT, V, N, T, sname); }

// Recursively call variadic template
save(others ...);

return; }

