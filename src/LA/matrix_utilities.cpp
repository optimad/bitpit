/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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

#include "bitpit_common.hpp"

#include "matrix_utilities.hpp"

namespace bitpit {

namespace linearalgebra {

/*!
 * @ingroup laauxiliary
 * @{
 */

/*!
    Computes the linear index of an element of a matrix stored in column-major
    ordering.

    \param[in] row is the row
    \param[in] col is the column
    \param[in] nRows is the number of rows
    \param[in] nCols is the number of columns
    \return The linear index of the specified element.
*/
int linearIndexColMajor(int row, int col, int nRows, int nCols)
{
    BITPIT_UNUSED(nCols);

    return row + col * nRows;
}

/*!
    Computes the linear index of an element of a matrix stored in row-major
    ordering.

    \param[in] row is the row
    \param[in] col is the column
    \param[in] nRows is the number of rows
    \param[in] nCols is the number of columns
    \return The linear index of the specified element.
*/
int linearIndexRowMajor(int row, int col, int nRows, int nCols)
{
    BITPIT_UNUSED(nRows);

    return row * nCols + col;
}

/*!
    Computes the linear index of an element of a symmetric matrix stored in
    column-major ordering when only either the upper or lower triangle is stored.
    E.g. if the indices (row,col) correspond to element in the lower
    triangle, but uplo indicates that upper triangle is stored, the symmetric
    index within the upper triangle is returned.

    \param[in] row is the row
    \param[in] col is the column
    \param[in] nRows is the number of rows
    \param[in] nCols is the number of columns
    \param[in] uplo defines if the matrix is an upper triangle ('U') or a
    lower triangle ('L')
    \return The linear index of the specified element.
*/
int linearIndexColMajorSymmetric(int row, int col, int nRows, int nCols, char uplo)
{
    assert(uplo == 'L' || uplo == 'U');

    if ((uplo == 'U' && col < row) || (uplo == 'L' && col > row)) {
        return linearIndexColMajor(col, row, nRows, nCols);
    } else {
        return linearIndexColMajor(row, col, nRows, nCols);
    }
}

/*!
    Computes the linear index of an element of a symmetric matrix stored in
    row-major ordering when only either the upper or lower triangle is stored.
    E.g. if the indices (row,col) correspond to element in the lower
    triangle, but uplo indicates that upper triangle is stored, the symmetric
    index within the upper triangle is returned.

    \param[in] row is the row
    \param[in] col is the column
    \param[in] nRows is the number of rows
    \param[in] nCols is the number of columns
    \param[in] uplo defines if the matrix is an upper triangle ('U') or a
    lower triangle ('L')
    \return The linear index of the specified element.
*/
int linearIndexRowMajorSymmetric(int row, int col, int nRows, int nCols, char uplo)
{
    assert(uplo == 'L' || uplo == 'U');

    if ((uplo == 'U' && col < row) || (uplo == 'L' && col > row)) {
        return linearIndexRowMajor(col, row, nRows, nCols);
    } else {
        return linearIndexRowMajor(row, col, nRows, nCols);
    }
}

/*!
 * @}
 */

}

}
