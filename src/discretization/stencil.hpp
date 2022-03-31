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

#ifndef __BTPIT_STENCIL_HPP__
#define __BTPIT_STENCIL_HPP__

#include <array>
#include <ostream>
#include <unordered_map>

#include "bitpit_common.hpp"
#include "bitpit_containers.hpp"

// Stream operators for the stencil class
namespace bitpit {

template<typename weight_t>
class DiscreteStencil;

}

template<typename weight_t>
bitpit::OBinaryStream & operator<<(bitpit::OBinaryStream &buffer, const bitpit::DiscreteStencil<weight_t> &stencil);

template<typename weight_t>
bitpit::IBinaryStream & operator>>(bitpit::IBinaryStream &buffer, bitpit::DiscreteStencil<weight_t> &stencil);

// Declaration of the stencil class
namespace bitpit {

/**
* \ingroup discretization
*
* \brief Metafunction for generating a discretization stencil.
*
* \tparam weight_t is the type of the weights stored in the stencil
*/
template <typename weight_t>
class DiscreteStencil {

template<typename U>
friend bitpit::OBinaryStream & (::operator<<) (bitpit::OBinaryStream &buffer, const DiscreteStencil<U> &stencil);
template<typename U>
friend bitpit::IBinaryStream & (::operator>>) (bitpit::IBinaryStream &buffer, DiscreteStencil<U> &stencil);

public:
    long NULL_ID = - std::numeric_limits<long>::max();

    typedef weight_t weight_type;

    /**
    * Defines an item of the stencil
    */
    struct item_type {
        long id;
        weight_type weight;
    };

    DiscreteStencil(const weight_t &zero = weight_t());
    DiscreteStencil(std::size_t nItems, const weight_t &zero = weight_t());
    DiscreteStencil(std::size_t size, const long *pattern, const weight_t &zero = weight_t());
    DiscreteStencil(std::size_t size, const long *pattern, const weight_t *weights, const weight_t &zero = weight_t());

    void initialize(std::size_t nItems, const weight_t &zero = weight_t());
    void initialize(std::size_t size, const long *pattern, const weight_t &zero = weight_t());
    void initialize(std::size_t size, const long *pattern, const weight_t *weights, const weight_t &zero = weight_t());
    void initialize(const DiscreteStencil<weight_t> &other);

    void clear(bool release = false);

    std::size_t size() const;

    void resize(std::size_t nItems);
    void reserve(std::size_t nItems);

    long & getPattern(std::size_t pos);
    const long & getPattern(std::size_t pos) const;
    long * patternData();
    const long * patternData() const;
    void setPattern(std::size_t pos, long id);

    weight_t & getWeight(std::size_t pos);
    const weight_t & getWeight(std::size_t pos) const;
    weight_t * weightData();
    const weight_t * weightData() const;
    void setWeight(std::size_t pos, const weight_t &weight);
    void setWeight(std::size_t pos, weight_t &&weight);
    void sumWeight(std::size_t pos, const weight_t &value, double factor = 1.);
    void zeroWeight(std::size_t pos);

    void setItem(std::size_t pos, long id, const weight_t &weight);
    void setItem(std::size_t pos, long id, weight_t &&weight);
    void sumItem(long id, const weight_t &value, double factor = 1.);
    void appendItem(long id, const weight_t &weight);
    void appendItem(long id, weight_t &&weight);

    weight_t & getConstant();
    const weight_t & getConstant() const;
    void setConstant(const weight_t &value);
    void setConstant(weight_t &&value);
    void sumConstant(const weight_t &value, double factor = 1.);
    void zeroConstant();

    void sum(const DiscreteStencil<weight_t> &other, double factor);

    void optimize(double tolerance = 1.e-12);
    void renumber(const std::unordered_map<long, long> &map);
    void addComplementToZero(long id);
    void zero();

    void display(std::ostream &out, double factor = 1.) const;

    size_t getBinarySize() const;

    DiscreteStencil<weight_t> & operator*=(double factor);
    DiscreteStencil<weight_t> & operator/=(double factor);
    DiscreteStencil<weight_t> & operator+=(const DiscreteStencil<weight_t> &other);
    DiscreteStencil<weight_t> & operator-=(const DiscreteStencil<weight_t> &other);

private:
    static void rawSumValue(const weight_t &value, double factor, weight_t *target);
    static void rawCopyValue(const weight_t &source, weight_t *destination);
    static void rawMoveValue(weight_t &&source, weight_t *destination);

    weight_t m_zero;
    std::vector<long> m_pattern;
    std::vector<weight_t> m_weights;
    weight_t m_constant;

    weight_t * findWeight(long id);
    const weight_t * findWeight(long id) const;

    bool optimizeWeight(std::size_t pos, double tolerance = 1.e-12);

};

// Methods for the spcializations
template<>
void DiscreteStencil<std::array<double, 3>>::rawSumValue(const std::array<double, 3> &value, double factor, std::array<double, 3> *target);

template<>
void DiscreteStencil<std::vector<double>>::rawSumValue(const std::vector<double> &value, double factor, std::vector<double> *target);

template<>
void DiscreteStencil<std::array<double, 3>>::rawCopyValue(const std::array<double, 3> &source, std::array<double, 3> *target);

template<>
void DiscreteStencil<std::vector<double>>::rawCopyValue(const std::vector<double> &source, std::vector<double> *target);

template<>
bool DiscreteStencil<std::array<double, 3>>::optimizeWeight(std::size_t pos, double tolerance);

template<>
bool DiscreteStencil<std::vector<double>>::optimizeWeight(std::size_t pos, double tolerance);

}

// Operators for the stencil class
template <typename weight_t>
bitpit::DiscreteStencil<weight_t> operator*(const bitpit::DiscreteStencil<weight_t> &stencil, double factor);

template <typename weight_t>
bitpit::DiscreteStencil<weight_t> operator*(double factor, const bitpit::DiscreteStencil<weight_t> &stencil);

template <typename weight_t>
bitpit::DiscreteStencil<weight_t> operator/(const bitpit::DiscreteStencil<weight_t> &stencil, double factor);

template <typename weight_t>
bitpit::DiscreteStencil<weight_t> operator+(const bitpit::DiscreteStencil<weight_t> &stencil_A, const bitpit::DiscreteStencil<weight_t> &stencil_B);

template <typename weight_t>
bitpit::DiscreteStencil<weight_t> operator-(const bitpit::DiscreteStencil<weight_t> &stencil_A, const bitpit::DiscreteStencil<weight_t> &stencil_B);

// Template implementation
#include "stencil.tpp"

// Declaration of the typdefs
namespace bitpit {

typedef DiscreteStencil<double> StencilScalar;
typedef DiscreteStencil<std::array<double, 3>> StencilVector;
typedef DiscreteStencil<std::vector<double>> StencilBlock;

}

// Explicit instantization
#ifndef __BTPIT_STENCIL_SRC__
namespace bitpit {

extern template class DiscreteStencil<double>;
extern template class DiscreteStencil<std::array<double, 3>>;
extern template class DiscreteStencil<std::vector<double>>;

}
#endif

// Operators for the specializations
bitpit::StencilVector operator*(const bitpit::StencilScalar &stencil, const std::array<double,3> &vector);
bitpit::StencilVector operator*(const std::array<double,3> &vector, const bitpit::StencilScalar &stencil);

// Functions for the specializations
bitpit::StencilScalar dotProduct(const bitpit::StencilVector &stencil, const bitpit::StencilVector::weight_type &vector);
void dotProduct(const bitpit::StencilVector &stencil, const bitpit::StencilVector::weight_type &vector, bitpit::StencilScalar *stencil_dotProduct);

bitpit::StencilVector project(const bitpit::StencilVector &stencil, const std::array<double, 3> &direction);
void project(const bitpit::StencilVector &stencil, const std::array<double, 3> &direction, bitpit::StencilVector *stencil_projection);

#endif
