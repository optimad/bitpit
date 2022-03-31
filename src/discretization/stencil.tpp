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

#ifndef __BITPIT_STENCIL_TPP__
#define __BITPIT_STENCIL_TPP__

/*!
* Output stream operator from class DiscreteStencil to communication buffer.
*
* \param[in] buffer is the output memory stream
* \param[in] stencil is the stencil to be streamed
* \result Returns the same output stream received in input.
*/
template<typename weight_t>
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const bitpit::DiscreteStencil<weight_t> &stencil)
{
    buffer << stencil.m_zero;

    std::size_t nItems = stencil.size();
    buffer << nItems;

    const long *patternData = stencil.patternData();
    const weight_t *weightData = stencil.weightData();
    for (std::size_t n = 0; n < nItems; ++n) {
        buffer << patternData[n];
        buffer << weightData[n];
    }

    buffer << stencil.m_constant;

    return buffer;
}

/*!
* Input stream operator from class DiscreteStencil to communication buffer.
*
* \param[in] buffer is the input memory stream
* \param[in] stencil is the stencil to be streamed
* \result Returns the same input stream received in input.
*/
template<typename weight_t>
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, bitpit::DiscreteStencil<weight_t> &stencil)
{
    buffer >> stencil.m_zero;

    std::size_t nItems;
    buffer >> nItems;

    stencil.resize(nItems);
    long *patternData = stencil.patternData();
    weight_t *weightData = stencil.weightData();
    for (std::size_t n = 0; n < nItems; ++n) {
        buffer >> patternData[n];
        buffer >> weightData[n];
    }

    buffer >> stencil.m_constant;

    return buffer;
}

namespace bitpit {

/*!
* Constructor
*
* \param zero is the value to be used as zero
*/
template<typename weight_t>
DiscreteStencil<weight_t>::DiscreteStencil(const weight_t &zero)
    : DiscreteStencil(0, zero)
{
}

/*!
* Initialize the stencil
*
* \param size is the stencil size, expressed in number of elements
* \param zero is the value to be used as zero
*/
template<typename weight_t>
DiscreteStencil<weight_t>::DiscreteStencil(std::size_t size, const weight_t &zero)
    : m_zero(zero),
      m_pattern(size, -1), m_weights(size, m_zero),
      m_constant(m_zero)
{
}

/*!
* Constructor
*
* \param size is the stencil size, expressed in number of elements
* \param pattern is the patterns of the stencil
* \param zero is the value to be used as zero
*/
template<typename weight_t>
DiscreteStencil<weight_t>::DiscreteStencil(std::size_t size, const long *pattern, const weight_t &zero)
    : m_zero(zero),
      m_pattern(pattern, pattern + size), m_weights(size, m_zero),
      m_constant(m_zero)
{
}

/*!
* Constructor
*
* \param size is the stencil size, expressed in number of elements
* \param pattern is the patterns of the stencil
* \param weights are the weights of the stencil
* \param zero is the value to be used as zero
*/
template<typename weight_t>
DiscreteStencil<weight_t>::DiscreteStencil(std::size_t size, const long *pattern, const weight_t *weights, const weight_t &zero)
    : m_zero(zero),
      m_pattern(pattern, pattern + size), m_weights(weights, weights + size),
      m_constant(m_zero)
{
}

/*!
* Initialize the stencil.
*
* \param size is the stencil size, expressed in number of elements
* \param zero is the value to be used as zero
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::initialize(std::size_t size, const weight_t &zero)
{
    rawCopyValue(zero, &m_zero);

    resize(size);
    for (std::size_t n = 0; n < size; ++n) {
        m_pattern[n] = -1;
        rawCopyValue(m_zero, m_weights.data() + n);
    }

    zeroConstant();
}

/*!
* Initialize the stencil.
*
* \param size is the stencil size, expressed in number of elements
* \param pattern is the patterns of the stencil
* \param zero is the value to be used as zero
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::initialize(std::size_t size, const long *pattern, const weight_t &zero)
{
    rawCopyValue(zero, &m_zero);

    resize(size);
    for (std::size_t n = 0; n < size; ++n) {
        m_pattern[n] = pattern[n];
        rawCopyValue(m_zero, m_weights.data() + n);
    }

    zeroConstant();
}

/*!
* Initialize the stencil.
*
* \param size is the stencil size, expressed in number of elements
* \param pattern is the patterns of the stencil
* \param weights are the weights of the stencil
* \param zero is the value to be used as zero
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::initialize(std::size_t size, const long *pattern, const weight_t *weights, const weight_t &zero)
{
    rawCopyValue(zero, &m_zero);

    resize(size);
    for (std::size_t n = 0; n < size; ++n) {
        m_pattern[n] = pattern[n];
        rawCopyValue(weights[n], m_weights.data() + n);
    }

    zeroConstant();
}

/*!
* Initialize the stencil.
*
* \param other is another stencil of the same time, whose items will be used
* to initialize this stencil
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::initialize(const DiscreteStencil<weight_t> &other)
{
    rawCopyValue(other.m_zero, &m_zero);

    std::size_t nItems = other.size();
    resize(nItems);
    for (std::size_t n = 0; n < nItems; ++n) {
        m_pattern[n] = other.m_pattern[n];
        rawCopyValue(other.m_weights[n], m_weights.data() + n);
    }

    setConstant(other.m_constant);
}

/*!
* Get the total size of the stencil, expressed in number of items.
*
* \result The total size of the stencil, expressed in number of items.
*/
template<typename weight_t>
std::size_t DiscreteStencil<weight_t>::size() const
{
    return m_pattern.size();
}

/*!
* Resizes the container so that it contains the specified number of items.
*
* \param size is the new stencil size, expressed in number of elements
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::resize(std::size_t size)
{
    m_pattern.resize(size, -1);
    m_weights.resize(size, m_zero);
}

/*!
* Requests a change in capacity.
*
* Requests that the stencil capacity be at least the specified value.
*
* \param capacity is the minimum number of items that the stencil should
* be able to contain
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::reserve(std::size_t capacity)
{
    m_pattern.reserve(capacity);
    m_weights.reserve(capacity);
}

/*!
* Get a reference to the specified element of the pattern.
*
* \param pos is the position of the pattern element
* \result A reference to the specified element of the pattern.
*/
template<typename weight_t>
long & DiscreteStencil<weight_t>::getPattern(std::size_t pos)
{
    return m_pattern[pos];
}

/*!
* Get a constant reference to the specified element of the pattern.
*
* \param pos is the position of the pattern element
* \result A constant reference to the specified element of the pattern.
*/
template<typename weight_t>
const long & DiscreteStencil<weight_t>::getPattern(std::size_t pos) const
{
    return m_pattern[pos];
}

/*!
* Get a pointer to the underlying array serving as pattern storage.
*
* \result A pointer to the underlying array serving as pattern storage.
*/
template<typename weight_t>
long * DiscreteStencil<weight_t>::patternData()
{
    return m_pattern.data();
}

/*!
* Get a constant pointer to the underlying array serving as pattern storage.
*
* \result A constant pointer to the underlying array serving as pattern
* storage.
*/
template<typename weight_t>
const long * DiscreteStencil<weight_t>::patternData() const
{
    return m_pattern.data();
}

/*!
* Set the index of the specified element of the pattern.
*
* \param pos is the position of the pattern element
* \param id is the index that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::setPattern(std::size_t pos, long id)
{
    m_pattern[pos] = id;
}

/*!
* Get a reference to the specified weight of the stencil.
*
* \param pos is the position of the weight
* \result A reference to the specified weight of the stencil.
*/
template<typename weight_t>
weight_t & DiscreteStencil<weight_t>::getWeight(std::size_t pos)
{
    return m_weights[pos];
}

/*!
* Get a constant reference to the specified weight of the stencil.
*
* \param pos is the position of the weight
* \result A constant reference to the specified weight of the stencil.
*/
template<typename weight_t>
const weight_t & DiscreteStencil<weight_t>::getWeight(std::size_t pos) const
{
    return m_weights[pos];
}

/*!
* Get a pointer to the underlying array serving as weight storage.
*
* \result A pointer to the underlying array serving as weight storage.
*/
template<typename weight_t>
weight_t * DiscreteStencil<weight_t>::weightData()
{
    return m_weights.data();
}

/*!
* Get a constant pointer to the underlying array serving as weight storage.
*
* \result A constant pointer to the underlying array serving as weight storage.
*/
template<typename weight_t>
const weight_t * DiscreteStencil<weight_t>::weightData() const
{
    return m_weights.data();
}

/*!
* Set the value of the specified weight of the stencil.
*
* \param pos is the position of the weight
* \param weight is the value that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::setWeight(std::size_t pos, const weight_t &weight)
{
    m_weights[pos] = weight;
}

/*!
* Set the value of the specified weight of the stencil.
*
* \param pos is the position of the weight
* \param weight is the value that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::setWeight(std::size_t pos, weight_t &&weight)
{
    m_weights[pos] = std::move(weight);
}

/*!
* Sum the given value to the weight at the specified position of the stencil.
*
* \param pos is the position of the weight
* \param value is the value that will be summed
* \param factor is the factor by which the value will be multiplied
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::sumWeight(std::size_t pos, const weight_t &value, double factor)
{
    rawSumValue(value, factor, m_weights.data() + pos);
}

/*!
* Zeros the weight at the specified position of the stencil.
*
* \param pos is the position of the weight
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::zeroWeight(std::size_t pos)
{
    rawCopyValue(m_zero, m_weights.data() + pos);
}

/*!
* Set the specified item of the stencil.
*
* \param pos is the position of the weight
* \param id is the index that will be set
* \param weight is the weight that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::setItem(std::size_t pos, long id, const weight_t &weight)
{
    setPattern(pos, id);
    setWeight(pos, weight);
}

/*!
* Set the specified item of the stencil.
*
* \param pos is the position of the weight
* \param id is the index that will be set
* \param weight is the weight that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::setItem(std::size_t pos, long id, weight_t &&weight)
{
    setPattern(pos, id);
    setWeight(pos, std::move(weight));
}

/*!
* Sum the given value to the item of the stencil with the specified index.
*
* If an item with the same id already exists, the given value will be summed
* to the weight of the existing item. Otherwise, a new item will be appended.
*
* \param id is the index of the item
* \param value is the value that will be summed
* \param factor is the factor by which the value will be multiplied
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::sumItem(long id, const weight_t &value, double factor)
{
    weight_t *weight = findWeight(id);
    if (weight) {
        rawSumValue(value, factor, weight);
    } else {
        appendItem(id, value);
        m_weights.back() *= factor;
    }
}

/*!
* Append an item the stencil.
*
* The item will be appended to the stencil also if the stencil already
* contains an item with the same id.
*
* \param id is the index that will be set
* \param weight is the weight that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::appendItem(long id, const weight_t &weight)
{
    m_pattern.push_back(id);
    m_weights.push_back(weight);
}

/*!
* Append an item the stencil.
*
* The item will be appended to the stencil also if the stencil already
* contains an item with the same id.
*
* \param id is the index that will be set
* \param weight is the weight that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::appendItem(long id, weight_t &&weight)
{
    m_pattern.push_back(id);
    m_weights.push_back(std::move(weight));
}

/*!
* Get a constant reference to the constant associated to the stencil.
*
* \result A constant reference to the constant associated to the stencil.
*/
template<typename weight_t>
const weight_t & DiscreteStencil<weight_t>::getConstant() const
{
    return m_constant;
}

/*!
* Get a reference to the constant associated to the stencil.
*
* \result A reference to the constant associated to the stencil.
*/
template<typename weight_t>
weight_t & DiscreteStencil<weight_t>::getConstant()
{
    return m_constant;
}

/*!
* Set the value of the constant associated to the stencil.
*
* \param value is the value that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::setConstant(const weight_t &value)
{
    rawCopyValue(value, &m_constant);
}

/*!
* Set the value of the constant associated to the stencil.
*
* \param value is the value that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::setConstant(weight_t &&value)
{
    rawMoveValue(std::move(value), &m_constant);
}

/*!
* Sum the specified value to the constant associated to the stencil.
*
* \param value is the value that will be summed
* \param factor is the factor by which the value will be multiplied
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::sumConstant(const weight_t &value, double factor)
{
    rawSumValue(value, factor, &m_constant);
}

/*!
* Zero the constant associated to the stencil.
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::zeroConstant()
{
    rawCopyValue(m_zero, &m_constant);
}

/*!
* Clears the items of the stencil.
*
* Removes all items from the stencil (which are destroyed), leaving it
* with a size of 0.
*
* \param release if it's true the memory hold by the stencil will be
* released, otherwise the stencil will be cleared but its memory will
* not be relased
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::clear(bool release)
{
    m_pattern.clear();
    m_weights.clear();

    if (release) {
        m_pattern.shrink_to_fit();
        m_weights.shrink_to_fit();
    }

    zeroConstant();
}

/*!
* Sum the specified stencil.
*
* \param other is the stencil that will be summed
* \param factor is the factor by which the other stencil will be multiplied
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::sum(const DiscreteStencil<weight_t> &other, double factor)
{
    const std::size_t other_nItems = other.size();
    for (std::size_t n = 0; n < other_nItems; ++n) {
        sumItem(other.m_pattern[n], other.m_weights[n], factor);
    }

    sumConstant(other.m_constant, factor);
}

/*!
* Optimize the stencil.
*
* The negligible elements will be removed from the stencil.
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::optimize(double tolerance)
{
    std::size_t nItems = size();
    for (std::size_t n = 0; n < nItems; ++n) {
        bool isWeightNeglibile = optimizeWeight(n, tolerance);
        if (isWeightNeglibile) {
            m_pattern.erase(m_pattern.begin() + n);
            m_weights.erase(m_weights.begin() + n);
            --n;
            --nItems;
        }
    }
}

/*!
* Renumber the indexes of the stencil according to the specified map.
*
* \param map is the renumbering map
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::renumber(const std::unordered_map<long, long> &map)
{
    const std::size_t nItems = size();
    for (std::size_t n = 0; n < nItems; ++n) {
        long &id = m_pattern[n];
        id = map.at(id);
    }
}

/*!
* Add a new item that complements the stencil to zero.
*
* \param id is the index associated to the new item
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::addComplementToZero(long id)
{
    const std::size_t nItems = size();
    if (nItems == 0) {
        return;
    }

    weight_t complement = m_zero;
    for (std::size_t n = 0; n < nItems; ++n) {
        rawSumValue(m_weights[n], -1., &complement);
    }

    appendItem(id, complement);
}

/*!
* Set weights and constant to zero.
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::zero()
{
    const std::size_t nItems = size();
    for (std::size_t n = 0; n < nItems; ++n) {
        weight_t &weight = m_weights[n];
        rawCopyValue(m_zero, &weight);
    }

    setConstant(m_zero);
}
/*!
* Optimize the specified weight.
*
* \param pos is the position of the weight to check
* \param tolerance is the tolerance that will be used for the check
* \result Returns true if the whole weight is neglibile
*/
template<typename weight_t>
bool DiscreteStencil<weight_t>::optimizeWeight(std::size_t pos, double tolerance)
{
    return (std::abs(m_weights[pos] - m_zero) <= tolerance);
}

/*!
* Find the weight associated to the specified id.
*
* \param id is the index associated to the weight
* \result A pointer to the weight associated to the specified id. If the weight
* is not found, a null pointer is returned.
*/
template<typename weight_t>
weight_t * DiscreteStencil<weight_t>::findWeight(long id)
{
    const std::size_t nItems = size();
    for (std::size_t n = 0; n < nItems; ++n) {
        if (m_pattern[n] == id) {
            return (m_weights.data() + n);
        }
    }

    return nullptr;
}

/*!
* Find the weight associated to the specified id.
*
* \param id is the index associated to the weight
* \result A constant pointer to the weight associated to the specified id.
* If the weight is not found, a null pointer is returned.
*/
template<typename weight_t>
const weight_t * DiscreteStencil<weight_t>::findWeight(long id) const
{
    const std::size_t nItems = size();
    for (std::size_t n = 0; n < nItems; ++n) {
        if (m_pattern[n] == id) {
            return (m_weights.data() + n);
        }
    }

    return nullptr;
}

/*!
* Sum the specified value to the target.
*
* \param value is the value that will be summed
* \param factor is the factor the value will be multiplied with
* \param target on output will contain the original weight plus the value multiplied by the factor
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::rawSumValue(const weight_t &value, double factor, weight_t *target)
{
    *target += factor * value;
}

/**
 * Copy the source value into the target.
 *
 * \param source is the value that will be copied
 * \param[out] target on output will contain the source value
 */
template<typename weight_t>
void DiscreteStencil<weight_t>::rawCopyValue(const weight_t &source, weight_t *target)
{
    *target = source;
}

/**
 * Move the source value into the target.
 *
 * \param[in,out] source is the value that will be moved
 * \param[out] target on output will contain the source value
 */
template<typename weight_t>
void DiscreteStencil<weight_t>::rawMoveValue(weight_t &&source, weight_t *target)
{
    *target = std::move(source);
}

/*!
* Display the stencil.
*
* \param out is the stream that will be used for the output
* \param factor is an optional factor the weights will be muliplied for
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::display(std::ostream &out, double factor) const
{
    const std::size_t nItems = size();

    weight_t sum = m_zero;
    for (std::size_t n = 0; n < nItems; ++n) {
        long id = m_pattern[n];
        weight_t value = factor * m_weights[n];
        out << "   id: " << id << " weight: " << value << std::endl;
        rawSumValue(value, 1., &sum);
    }

    out << " constant : " << (factor * m_constant) << std::endl;
    out << " sum      : " << sum << std::endl;
}

/*!
* Returns the buffer size (in bytes) required to store the stencil.
*
* \result The buffer size (in bytes) required to store the stencil.
*/
template<typename weight_t>
size_t DiscreteStencil<weight_t>::getBinarySize() const
{
    std::size_t nItems = size();

    return (sizeof(m_zero) + nItems * (sizeof(long) + sizeof(weight_t)) + sizeof(m_constant));
}

/*!
* The binary multiplication assignment operator.
*
* \param factor is the factor of the multiplication
* \result A reference to the stencil
*/
template<typename weight_t>
DiscreteStencil<weight_t> & DiscreteStencil<weight_t>::operator*=(double factor)
{
    const std::size_t nItems = size();
    for (std::size_t n = 0; n < nItems; ++n) {
        m_weights[n] *= factor;
    }
    m_constant *= factor;

    return *this;
}

/*!
* The binary division assignment operator.
*
* \param factor is the factor of the division
* \result A reference to the stencil
*/
template<typename weight_t>
DiscreteStencil<weight_t> & DiscreteStencil<weight_t>::operator/=(double factor)
{
    const std::size_t nItems = size();
    for (std::size_t n = 0; n < nItems; ++n) {
        m_weights[n] /= factor;
    }
    m_constant /= factor;

    return *this;
}

/*!
* The binary sum assignment operator.
*
* \param other is the stencil that will be summed
* \result A reference to the stencil
*/
template<typename weight_t>
DiscreteStencil<weight_t> & DiscreteStencil<weight_t>::operator+=(const DiscreteStencil<weight_t> &other)
{
    sum(other, 1.);

    return *this;
}

/*!
* The binary subtraction assignment operator.
*
* \param other is the stencil that will be subtracted
* \result A reference to the stencil
*/
template<typename weight_t>
DiscreteStencil<weight_t> & DiscreteStencil<weight_t>::operator-=(const DiscreteStencil<weight_t> &other)
{
    sum(other, - 1.);

    return *this;
}

}

/*!
* The multiplication operator between a stencil and a scalar value.
*
* \param stencil is the stencil
* \param factor is the factor of the multiplication
* \result The result fo the multiplication.
*/
template<typename weight_t>
bitpit::DiscreteStencil<weight_t> operator*(const bitpit::DiscreteStencil<weight_t> &stencil, double factor)
{
    return (factor * stencil);
}

/*!
* The multiplication operator between scalar value and a stencil.
*
* \param factor is the factor of the multiplication
* \param stencil is the stencil
* \result The result fo the multiplication.
*/
template<typename weight_t>
bitpit::DiscreteStencil<weight_t> operator*(double factor, const bitpit::DiscreteStencil<weight_t> &stencil)
{
    bitpit::DiscreteStencil<weight_t> stencil_result(stencil);
    stencil_result *= factor;

    return stencil_result;
}

/*!
* The division operator between a stencil and a scalar value.
*
* \param stencil is the stencil
* \param factor is the factor of the division
* \result The result fo the division.
*/
template<typename weight_t>
bitpit::DiscreteStencil<weight_t> operator/(const bitpit::DiscreteStencil<weight_t> &stencil, double factor)
{
    bitpit::DiscreteStencil<weight_t> stencil_result(stencil);
    stencil_result /= factor;

    return stencil_result;
}

/*!
* The sum operator between two stencils.
*
* \param stencil_A is the first stencil
* \param stencil_B is the second stencil
* \result The result fo the sum.
*/
template<typename weight_t>
bitpit::DiscreteStencil<weight_t> operator+(const bitpit::DiscreteStencil<weight_t> &stencil_A, const bitpit::DiscreteStencil<weight_t> &stencil_B)
{
    bitpit::DiscreteStencil<weight_t> stencil_result(stencil_A);
    stencil_result += stencil_B;

    return stencil_result;
}

/*!
* The subtraction operator between two stencils.
*
* \param stencil_A is the first stencil
* \param stencil_B is the second stencil
* \result The result fo the subtraction.
*/
template<typename weight_t>
bitpit::DiscreteStencil<weight_t> operator-(const bitpit::DiscreteStencil<weight_t> &stencil_A, const bitpit::DiscreteStencil<weight_t> &stencil_B)
{
    bitpit::DiscreteStencil<weight_t> stencil_result(stencil_A);
    stencil_result -= stencil_B;

    return stencil_result;
}

#endif
