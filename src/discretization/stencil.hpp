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
#include <vector>

#include "bitpit_common.hpp"
#include "bitpit_containers.hpp"
#include "bitpit_operators.hpp"

namespace bitpit {

// Stencil weight pool
template<typename weight_t>
class DiscreteStencilWeightPool
{

public:
    DiscreteStencilWeightPool(std::size_t capacity = 128);

    std::size_t size() const;
    std::size_t capacity() const;
    void clear(bool release);

    weight_t retrieve();
    void store(weight_t &&weight);
    void store(std::vector<weight_t> *weights);

private:
    std::size_t m_capacity;
    std::vector<weight_t> m_storage;

};

// Stream operators for the stencil class
template<typename weight_t>
class DiscreteStencil;

template<typename weight_t>
OBinaryStream & operator<<(OBinaryStream &buffer, const DiscreteStencil<weight_t> &stencil);

template<typename weight_t>
IBinaryStream & operator>>(IBinaryStream &buffer, DiscreteStencil<weight_t> &stencil);

/**
* \ingroup discretization
*
* \brief Metafunction for generating a discretization stencil.
*
* \tparam weight_t is the type of the weights stored in the stencil
*/
template <typename weight_t>
class DiscreteStencil {

template<typename W>
friend OBinaryStream & (operator<<) (OBinaryStream &buffer, const DiscreteStencil<W> &stencil);
template<typename W>
friend IBinaryStream & (operator>>) (IBinaryStream &buffer, DiscreteStencil<W> &stencil);

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

    virtual ~DiscreteStencil() = default;

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
    template<typename Mapper>
    void renumber(const Mapper &mapper);
    void addComplementToZero(long id);
    void zero();

    void display(std::ostream &out, double factor = 1.) const;

    size_t getBinarySize() const;

    weight_t & at(long id);
    const weight_t & at(long id) const;

    weight_t & operator[](long id);

    DiscreteStencil<weight_t> & operator*=(double factor);
    DiscreteStencil<weight_t> & operator/=(double factor);
    DiscreteStencil<weight_t> & operator+=(const DiscreteStencil<weight_t> &other);
    DiscreteStencil<weight_t> & operator-=(const DiscreteStencil<weight_t> &other);

protected:
    weight_t m_zero;
    std::vector<long> m_pattern;
    std::vector<weight_t> m_weights;
    weight_t m_constant;

    virtual void appendWeight(weight_t &&weight);
    virtual void appendWeight(const weight_t &weight);

    virtual void clearWeights(bool release);

private:
    template<typename W>
    static void rawSumValue(const W &value, double factor, W *target);
    template<typename W = weight_t, typename V = typename W::value_type, long unsigned int D = W::size_type, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type * = nullptr>
    static void rawSumValue(const std::array<V, D> &value, double factor, std::array<V, D> *target);
    template<typename W = weight_t, typename V = typename W::value_type, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type * = nullptr>
    static void rawSumValue(const std::vector<V> &value, double factor, std::vector<V> *target);

    template<typename W>
    static void rawCopyValue(const W &source, W *target);
    template<typename W = weight_t, typename V = typename W::value_type, long unsigned int D = W::size_type, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type * = nullptr>
    static void rawCopyValue(const std::array<V, D> &source, std::array<V, D> *target);
    template<typename W = weight_t, typename V = typename W::value_type, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type * = nullptr>
    static void rawCopyValue(const std::vector<V> &source, std::vector<V> *target);

    template<typename W>
    static void rawMoveValue(W &&source, W *target);

    weight_t * findWeight(long id);
    const weight_t * findWeight(long id) const;

    template<typename W>
    bool isWeightNegligible(const W &weight, double tolerance = 1.e-12) const;
    template<typename W = weight_t, typename V = typename W::value_type, long unsigned int D = W::size_type, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type * = nullptr>
    bool isWeightNegligible(const std::array<V, D> &weight, double tolerance = 1.e-12) const;
    template<typename W = weight_t, typename V = typename W::value_type, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type * = nullptr>
    bool isWeightNegligible(const std::vector<V> &weight, double tolerance = 1.e-12) const;

};

/**
* \ingroup discretization
*
* \brief Metafunction for generating a discretization stencil with a
* memory pool (MP).
*
* \tparam weight_t is the type of the weights stored in the stencil
*/
template<typename weight_t>
class MPDiscreteStencil : public DiscreteStencil<weight_t>
{

public:
    typedef DiscreteStencilWeightPool<weight_t> weight_pool_type;

    MPDiscreteStencil(const weight_t &zero = weight_t());
    MPDiscreteStencil(std::size_t nItems, const weight_t &zero = weight_t());
    MPDiscreteStencil(std::size_t size, const long *pattern, const weight_t &zero = weight_t());
    MPDiscreteStencil(std::size_t size, const long *pattern, const weight_t *weights, const weight_t &zero = weight_t());

    void setWeightPool(weight_pool_type *pool);

protected:
    weight_pool_type *m_weightPool;

    void appendWeight(const weight_t &weight) override;

    void clearWeights(bool release) override;

};

// Declaration of the typdefs
typedef DiscreteStencil<double> StencilScalar;
typedef DiscreteStencil<std::array<double, 3>> StencilVector;
typedef DiscreteStencil<std::vector<double>> StencilBlock;

typedef MPDiscreteStencil<double> MPStencilScalar;
typedef MPDiscreteStencil<std::array<double, 3>> MPStencilVector;
typedef MPDiscreteStencil<std::vector<double>> MPStencilBlock;

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

// Operators for the specializations
template <typename V>
typename bitpit::DiscreteStencil<std::array<V, 3>> operator*(const typename bitpit::DiscreteStencil<V> &stencil, const std::array<V, 3> &vector);
template <typename V>
typename bitpit::DiscreteStencil<std::array<V, 3>> operator*(const std::array<V, 3> &vector, const typename bitpit::DiscreteStencil<V> &stencil);

// Functions for the specializations
template <typename V>
typename bitpit::DiscreteStencil<V> dotProduct(const typename bitpit::DiscreteStencil<std::array<V, 3>> &stencil, const typename bitpit::DiscreteStencil<std::array<V, 3>>::weight_type &vector);
template <typename V>
void dotProduct(const typename bitpit::DiscreteStencil<std::array<V, 3>> &stencil, const typename bitpit::DiscreteStencil<std::array<V, 3>>::weight_type &vector, typename bitpit::DiscreteStencil<V> *stencil_dotProduct);

template <typename V>
typename bitpit::DiscreteStencil<std::array<V, 3>> project(const typename bitpit::DiscreteStencil<std::array<V, 3>> &stencil, const std::array<V, 3> &direction);
template <typename V>
void project(const typename bitpit::DiscreteStencil<std::array<V, 3>> &stencil, const std::array<V, 3> &direction, typename bitpit::DiscreteStencil<std::array<V, 3>> *stencil_projection);

// Template implementation
#include "stencil.tpp"

// Explicit instantization
#ifndef __BTPIT_STENCIL_SRC__
namespace bitpit {

extern template class DiscreteStencilWeightPool<double>;
extern template class DiscreteStencilWeightPool<std::array<double, 3>>;
extern template class DiscreteStencilWeightPool<std::vector<double>>;

extern template class DiscreteStencil<double>;
extern template class DiscreteStencil<std::array<double, 3>>;
extern template class DiscreteStencil<std::vector<double>>;

extern template class MPDiscreteStencil<double>;
extern template class MPDiscreteStencil<std::array<double, 3>>;
extern template class MPDiscreteStencil<std::vector<double>>;

}
#endif

#endif
