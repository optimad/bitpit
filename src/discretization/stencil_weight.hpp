/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2024 OPTIMAD engineering Srl
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

#ifndef __BTPIT_STENCIL_WEIGHT_HPP__
#define __BTPIT_STENCIL_WEIGHT_HPP__

#include "bitpit_common.hpp"

#include <array>
#include <vector>

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

/**
* \ingroup discretization
*
* Helper class to get information about the values associated with a weight.
*/
template <typename weight_t, typename = void>
class DiscreteStencilWeightValueInfo
{
public:
    using type = weight_t;
};

template <typename weight_t>
class DiscreteStencilWeightValueInfo<weight_t, std::void_t<typename weight_t::value_type>>
{
public:
    using type = typename weight_t::value_type;
};

/**
* \ingroup discretization
*
* Helper class to handle weights.
*/
template <typename weight_t, typename value_t>
class DiscreteStencilWeightManager
{
public:
    using weight_type = weight_t;
    using value_type  = value_t;

    template<typename W>
    bool isNegligible(const W &weight, const weight_t &zero, double tolerance = 1e-12) const;
    template<typename W = weight_t, typename V = value_t, std::size_t D = std::tuple_size<W>::value, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type * = nullptr>
    bool isNegligible(const std::array<V, D> &weight, const weight_t &zero, double tolerance = 1e-12) const;
    template<typename W = weight_t, typename V = value_t, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type * = nullptr>
    bool isNegligible(const std::vector<V> &weight, const weight_t &zero, double tolerance = 1e-12) const;

    template<typename W>
    void sum(const W &weight, double factor, W *target) const;
    template<typename W = weight_t, typename V = value_t, std::size_t D = std::tuple_size<W>::value, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type * = nullptr>
    void sum(const std::array<V, D> &weight, double factor, std::array<V, D> *target) const;
    template<typename W = weight_t, typename V = value_t, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type * = nullptr>
    void sum(const std::vector<V> &weight, double factor, std::vector<V> *target) const;

    template<typename W>
    void copy(const W &weight, W *target) const;
    template<typename W = weight_t, typename V = value_t, std::size_t D = std::tuple_size<W>::value, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type * = nullptr>
    void copy(const std::array<V, D> &weight, std::array<V, D> *target) const;
    template<typename W = weight_t, typename V = value_t, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type * = nullptr>
    void copy(const std::vector<V> &weight, std::vector<V> *target) const;

    template<typename W>
    void move(W &&weight, W *target) const;

    template<typename W>
    value_t & at(const W &weight, std::size_t index);
    template<typename W>
    const value_t & at(const W &weight, std::size_t index) const;
    template<typename W = weight_t, typename V, std::size_t D = std::tuple_size<W>::value, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type * = nullptr>
    value_t & at(const std::array<V, D> &weight, std::size_t index);
    template<typename W = weight_t, typename V, std::size_t D = std::tuple_size<W>::value, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type * = nullptr>
    const value_t & at(const std::array<V, D> &weight, std::size_t index) const;
    template<typename W = weight_t, typename V, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type * = nullptr>
    value_t & at(const std::vector<V> &weight, std::size_t index);
    template<typename W = weight_t, typename V, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type * = nullptr>
    const value_t & at(const std::vector<V> &weight, int index) const;

};

// Declaration of the typdefs
typedef DiscreteStencilWeightManager<double, double> DiscreteStencilScalarWeightManager;
typedef DiscreteStencilWeightManager<std::array<double, 3>, double> DiscreteStencilVectorWeightManager;
typedef DiscreteStencilWeightManager<std::vector<double>, double> DiscreteStencilBlockWeightManager;

}

// Template implementation
#include "stencil_weight.tpp"

// Explicit instantization
#ifndef __BITPIT_STENCIL_WEIGHT_SRC__
namespace bitpit {

extern template class DiscreteStencilWeightManager<double, double>;
extern template class DiscreteStencilWeightManager<std::array<double, 3>, double>;
extern template class DiscreteStencilWeightManager<std::vector<double>, double>;

}
#endif

#endif
