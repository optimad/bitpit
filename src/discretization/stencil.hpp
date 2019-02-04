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
    DiscreteStencil(int nBuckets, const weight_t &zero = weight_t());
    DiscreteStencil(int nBuckets, int nBucketItems, const weight_t &zero = weight_t());
    DiscreteStencil(const std::vector<int> &bucketSizes, const weight_t &zero = weight_t());

    void initialize(const weight_t &zero = weight_t());
    void initialize(int nBuckets, const weight_t &zero = weight_t());
    void initialize(int nBuckets, int nBucketItems, const weight_t &zero = weight_t());
    void initialize(const std::vector<int> &bucketSizes, const weight_t &zero = weight_t());
    void initialize(const DiscreteStencil<weight_t> &other);

    void clear(bool release = false);

    std::size_t size() const;
    std::size_t size(int bucket) const;

    void reserve(std::size_t nItems);
    void reserve(int nBuckets, std::size_t nItems);

    int getBucketCount() const;

    long & getPattern(std::size_t pos);
    const long & getPattern(std::size_t pos) const;
    long & getPattern(int bucket, std::size_t pos);
    const long & getPattern(int bucket, std::size_t pos) const;
    long * patternData();
    const long * patternData() const;
    const FlatVector2D<long> & getPattern() const;
    void setPattern(std::size_t pos, long id);
    void setPattern(int bucket, std::size_t pos, long id);

    long & rawGetPattern(std::size_t pos);
    const long & rawGetPattern(std::size_t pos) const;
    void rawSetPattern(std::size_t pos, long id);

    weight_t & getWeight(std::size_t pos);
    const weight_t & getWeight(std::size_t pos) const;
    weight_t & getWeight(int bucket, std::size_t pos);
    const weight_t & getWeight(int bucket, std::size_t pos) const;
    weight_t * weightData();
    const weight_t * weightData() const;
    const FlatVector2D<weight_t> & getWeights() const;
    void setWeight(std::size_t pos, const weight_t &weight);
    void setWeight(int bucket, std::size_t pos, const weight_t &weight);
    void sumWeight(std::size_t pos, const weight_t &value);
    void sumWeight(int bucket, std::size_t pos, const weight_t &value);

    weight_t & rawGetWeight(std::size_t pos);
    const weight_t & rawGetWeight(std::size_t pos) const;
    void rawSetWeight(std::size_t pos, const weight_t &weight);

    void setItem(std::size_t pos, long id, const weight_t &weight);
    void setItem(int bucket, std::size_t pos, long id, const weight_t &weight);
    void sumItem(long id, const weight_t &value);
    void sumItem(int bucket, long id, const weight_t &value);
    void appendItem(long id, const weight_t &weight);
    void appendItem(int bucket, long id, const weight_t &weight);

    weight_t & getConstant();
    const weight_t & getConstant() const;
    void setConstant(const weight_t &value);
    void sumConstant(const weight_t &value);

    void flatten();
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
    weight_t m_zero;
    FlatVector2D<long> m_pattern;
    FlatVector2D<weight_t> m_weights;
    weight_t m_constant;

    weight_t * findWeight(int bucket, long id);
    const weight_t * findWeight(int bucket, long id) const;

    template<typename U = weight_t, typename std::enable_if<std::is_fundamental<U>::value>::type * = nullptr>
    bool isWeightNeglibile(int bucket, std::size_t pos, double tolerance = 1.e-12);

    template<typename U = weight_t, typename std::enable_if<!std::is_fundamental<U>::value>::type * = nullptr>
    bool isWeightNeglibile(int bucket, std::size_t pos, double tolerance = 1.e-12);

};

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

}

// Explicit instantization
#ifndef __BTPIT_STENCIL_SRC__
namespace bitpit {

extern template class DiscreteStencil<double>;
extern template class DiscreteStencil<std::array<double, 3>>;

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
