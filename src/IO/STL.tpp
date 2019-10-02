/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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

/*!
    Load data for solid with specified label from ASCII STL file.
    
    \param[in] name label associated to the solid. If empty label is specified
    (i.e. name = ""), data of the first solid encountered in the STL file is returned.
    If no solid is found with the specified label, no data are loaded from the STL file.
    \param[in,out] nV on input stores the current number of vertices hosted in V.
    On output stores the input values incremented by the number
    of vertices acquired from the STL file.
    \param[in,out] nT on input stores the number of facet->vertex connectivity
    entries stores in T. On output stores the input value incremented by the
    number of facets acquired from the STL file.
    \param[in,out] V vertex coordinates list. On output stores the coordinates of
    vertices vertices acquired from the STL file. New vertices are appended
    at the end of V.
    \param[in,out] N facet normals. On output stores the normal unit vector to
    each facet acquired from the STL file. New normals are appended
    at the end of N.
    \param[in,out] T facet->vertex connectivity. On output stores the facet->vertex
    connectivity entries for the facets acquired from the STL file. New connectivity entries
    are appended at the end of T.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STLObj::load(std::string name, std::size_t &nV, std::size_t &nT, std::vector<std::vector<double>> &V,
                  std::vector<std::vector<double>> &N, std::vector<std::vector<std::size_t>> &T,
                  T2 & ... others)
{
    // Open stream
    open("in");
    if (err == 1) {
        return;
    }

    // Load solid
    if (stl_type) {
        readBINARY(m_ifile_handle, nV, nT, V, N, T);
    } else {
        readSolidASCII(m_ifile_handle, true, nV, nT, V, N, T, name);
    }

    // Recursively read other solids
    load(others ...);
}

/*!
    Load data for solid with specified label from an ASCII STL file. Overloading of
    member function STLObj::load() for container vector<array<double, 3>>
    
    \param[in] name label associated to the solid. If empty label is specified
    (i.e. name = ""), data of the first solid encountered in the STL file is returned.
    If no solid is found with the specified label, no data are loaded from the STL file.
    \param[in,out] nV on input stores the current number of vertices hosted in V.
    On output stores the input values incremented by the number
    of vertices acquired from the STL file.
    \param[in,out] nT on input stores the number of facet->vertex connectivity
    entries stores in T. On output stores the input value incremented by the
    number of facets acquired from the STL file.
    \param[in,out] V vertex coordinates list. On output stores the coordinates of
    vertices vertices acquired from the STL file. New vertices are appended
    at the end of V.
    \param[in,out] N facet normals. On output stores the normal unit vector to
    each facet acquired from the STL file. New normals are appended
    at the end of N.
    \param[in,out] T facet->vertex connectivity. On output stores the facet->vertex
    connectivity entries for the facets acquired from the STL file. New connectivity entries
    are appended at the end of T.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STLObj::load(std::string name, std::size_t &nV, std::size_t &nT, std::vector<std::array<double,3>> &V,
                  std::vector<std::array<double,3>> &N, std::vector<std::array<std::size_t,3>> &T,
                  T2 & ... others)
{
    // Open stream
    open("in");
    if (err == 1) {
        return;
    }

    // Load solid
    if (stl_type) {
        readBINARY(m_ifile_handle, nV, nT, V, N, T);
    } else {
        readSolidASCII(m_ifile_handle, true, nV, nT, V, N, T, name);
    }

    // Recursively read other solids
    load(others ...);
}

/*!
    Save solid data to ASCII STL file.
    
    \param[in] name label associated to the solid.
    \param[in,out] nV number of solid vertices.
    \param[in,out] nT number of solid facets.
    \param[in,out] V vertex coordinates list.
    \param[in,out] N facet normals. 
    \param[in,out] T facet->vertex connectivity.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STLObj::save(const std::string &name, std::size_t &nV, std::size_t &nT, std::vector<std::vector<double>> &V,
                  std::vector<std::vector<double>> &N, std::vector<std::vector<std::size_t>> &T,
                  T2 & ... others)
{
    // Open stream
    open("out");
    if (err == 1) {
        return;
    }

    // Export solid
    if (stl_type) {
        writeSolidBINARY(m_ofile_handle, nV, nT, V, N, T, name);
    } else {
        writeSolidASCII(m_ofile_handle, nV, nT, V, N, T, name);
    }

    // Recursively save other solids
    save(others ...);
}

/*!
    Save solid data to ASCII STL file. Overloading of member function STLObj::save()
    for container vector<array<double, 3>>
    
    \param[in] name label associated to the solid.
    \param[in,out] nV number of solid vertices.
    \param[in,out] nT number of solid facets.
    \param[in,out] V vertex coordinates list.
    \param[in,out] N facet normals. 
    \param[in,out] T facet->vertex connectivity.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STLObj::save(const std::string &name, std::size_t &nV, std::size_t &nT, std::vector<std::array<double,3>> &V,
                  std::vector<std::array<double,3>> &N, std::vector<std::array<std::size_t,3>> &T,
                  T2 & ... others)
{
    // Open stream
    open("out");
    if (err == 1) {
        return;
    }

    // Export solid
    if (stl_type) {
        writeSolidBINARY(m_ofile_handle, nV, nT, V, N, T, name);
    } else {
        writeSolidASCII(m_ofile_handle, nV, nT, V, N, T, name);
    }

    // Recursively save other solids
    save(others ...);
}

/*!
    Append solid data at the end of an existing ASCII STL file.
   
    \tparam T2 vriadic template 
    \param[in] name label associated to the solid.
    \param[in,out] nV number of solid vertices.
    \param[in,out] nT number of solid facets.
    \param[in,out] V vertex coordinates list.
    \param[in,out] N facet normals. 
    \param[in,out] T facet->vertex connectivity.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STLObj::append(const std::string &name, std::size_t &nV, std::size_t &nT, std::vector<std::vector<double>> &V,
                    std::vector<std::vector<double>> &N, std::vector<std::vector<std::size_t>> &T,
                    T2 & ... others)
{
    // Open stream
    open("app");
    if (err == 1) {
        return;
    }

    // Append solid
    if (stl_type) {
        throw std::runtime_error("Appending data to an existing file is only supported for ASCII files.");
    } else {
        writeSolidASCII(m_ofile_handle, nV, nT, V, N, T, name);
    }

    // Recursively append other solids
    append(others ...);
}

/*!
    Append solid data at the end of an existing ASCII STL file. Overloading of
    member function STLObj::append() for container vector<array<double, 3>>.
    
    \param[in] name label associated to the solid.
    \param[in,out] nV number of solid vertices.
    \param[in,out] nT number of solid facets.
    \param[in,out] V vertex coordinates list.
    \param[in,out] N facet normals. 
    \param[in,out] T facet->vertex connectivity.
    \param[in] others parameter packs
*/
template <typename ... T2>
void STLObj::append(const std::string &name, std::size_t &nV, std::size_t &nT, std::vector<std::array<double,3>> &V,
                    std::vector<std::array<double,3>> &N, std::vector<std::array<std::size_t,3>> &T,
                    T2 & ... others)
{
    // Open stream
    open("app");
    if (err == 1) {
        return;
    }

    // Append solid
    if (stl_type) {
        throw std::runtime_error("Appending data to an existing file is only supported for ASCII files.");
    } else {
        writeSolidASCII(m_ofile_handle, nV, nT, V, N, T, name);
    }

    // Recursively append other solids
    append(others ...);
}
