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

/*!
    \ingroup MathFunctions
    Power function with unsigned integer exponent.
    \param[in] base Input argument
    \param[in] exponent Power index

    Template parameter can be any type such that *operator is defined.

    \result returns the p-th power of the input argument.
*/
template <class T>
T uipow(const T & base, unsigned int exponent)
{
    switch (exponent) {

    case 0:
    {
        return 1;
    }

    case 1:
    {
        return base;
    }

    case 2:
    {
        return base * base;
    }

    case 3:
    {
        return base * base * base;
    }

    default:
    {
        // The integer power is evluated using the exponentiating by squaring
        // algorithm. For an explanation of the algorith see the following
        // link:
        //
        //    https://en.wikipedia.org/wiki/Exponentiation_by_squaring
        double result = 1.;
        double base_work = base;
        while (exponent) {
            if (exponent % 2) {
                result *= base_work;
            }

            exponent /= 2;
            base_work *= base_work;
        }

        return result;
    }

    }
}
