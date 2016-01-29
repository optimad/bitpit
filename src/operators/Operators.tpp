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

// Operator "sign" ------------------------------------------------------------------- //
/*!
    \ingroup MathFunctions
    Sign function.
    Given a a variable of integral type, val, returns:
    1 if val > 0
    -1, othersize.

    Template parameters can be any integral type such that operator< is defined.

    \result returns the sign of the input value.
*/
template <class T>
T sign(
        const  T                & val
      ){
    return (T(0) < val) - (val < T(0)) ;
};



