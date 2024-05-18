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

#ifndef __BITPIT_STENCIL_WEIGHT_TPP__
#define __BITPIT_STENCIL_WEIGHT_TPP__

namespace bitpit {

/*!
* Constructor.
*
* \param capacity is the maximum number of weights that can be stored
* in the pool
*/
template<typename weight_t>
DiscreteStencilWeightPool<weight_t>::DiscreteStencilWeightPool(std::size_t capacity)
    : m_capacity(capacity)
{
}

/*!
* Get the size of the pool.
*
* The size represents the number of weights currently stored in the pool.
*
* \param size The size of the pool.
*/
template<typename weight_t>
std::size_t DiscreteStencilWeightPool<weight_t>::size() const
{
    return m_storage.size();
}

/*!
* Get the capacity of the pool.
*
* The capacity represents the maximum number of weights that can be stored
* in the pool.
*
* \param size The capacity of the pool.
*/
template<typename weight_t>
std::size_t DiscreteStencilWeightPool<weight_t>::capacity() const
{
    return m_capacity;
}

/*!
* Clear the pool.
*
* Removes all weights from the pool (which are destroyed), leaving it
* with a size of 0.
*
* \param release if it's true the memory hold by the pool will be
* released, otherwise the pool will be cleared but its memory will
* not be relased
*/
template<typename weight_t>
void DiscreteStencilWeightPool<weight_t>::clear(bool release)
{
    m_storage.clear();
    if (release) {
        m_storage.shrink_to_fit();
    }
}

/*!
* Retrieve a weight from the pool.
*
* If the pool is empty, an exception is thrown.
*
* \result The weight retrieved from the pool.
*/
template<typename weight_t>
weight_t DiscreteStencilWeightPool<weight_t>::retrieve()
{
    if (size() == 0) {
        throw std::runtime_error("Unable to retrieve a weight from the pool: the pool is empty.");
    }

    weight_t weight = std::move(m_storage.back());
    m_storage.pop_back();

    return weight;
}

/*!
* Store the given weight in the pool.
*
* \param weight is the weight that will be stored in the pool
*/
template<typename weight_t>
void DiscreteStencilWeightPool<weight_t>::store(weight_t &&weight)
{
    if (m_capacity == size()) {
        return;
    }

    m_storage.emplace_back(std::move(weight));
}

/*!
* Store the given weights in the pool.
*
* \param weights are the weight that will be stored in the pool
*/
template<typename weight_t>
void DiscreteStencilWeightPool<weight_t>::store(std::vector<weight_t> *weights)
{
    std::size_t nStorableWeights = std::min(m_capacity - size(), weights->size());
    if (nStorableWeights == 0) {
        return;
    }

    m_storage.insert(m_storage.end(),
                     std::make_move_iterator(weights->begin()),
                     std::make_move_iterator(weights->begin() + nStorableWeights));
}

/*!
 * Check if the specified weight is negligible.
 *
 * \param weight is the weight to be checked
 * \param zero is the weight to be used as zero
 * \param tolerance is the tolerance that will be used for the check
 * \result Returns true if the whole weight is negligible
 */
template<typename weight_t, typename value_t>
template<typename W>
bool DiscreteStencilWeightManager<weight_t, value_t>::isNegligible(const W &weight, const weight_t &zero, double tolerance) const
{
    return (std::abs(weight - zero) <= tolerance);
}

/*!
 * Check if the specified weight is negligible.
 *
 * \param weight is the weight to be checked
 * \param zero is the weight to be used as zero
 * \param tolerance is the tolerance that will be used for the check
 * \result Returns true if the whole weight is negligible
 */
template<typename weight_t, typename value_t>
template<typename W, typename V, std::size_t D, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type *>
bool DiscreteStencilWeightManager<weight_t, value_t>::isNegligible(const std::array<V, D> &weight, const weight_t &zero, double tolerance) const
{
    for (std::size_t k = 0; k < D; ++k) {
        if (std::abs(weight[k] - zero[k]) > tolerance) {
            return false;
        }
    }

    return true;
}

/*!
 * Check if the specified weight is negligible.
 *
 * \param weight is the weight to be checked
 * \param zero is the weight to be used as zero
 * \param tolerance is the tolerance that will be used for the check
 * \result Returns true if the whole weight is negligible
 */
template<typename weight_t, typename value_t>
template<typename W, typename V, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type *>
bool DiscreteStencilWeightManager<weight_t, value_t>::isNegligible(const std::vector<V> &weight, const weight_t &zero, double tolerance) const
{
    const int nItems = weight.size();
    for (int k = 0; k < nItems; ++k) {
        if (std::abs(weight[k] - zero[k]) > tolerance) {
            return false;
        }
    }

    return true;
}


/*!
 * Sum the specified weight to the target.
 *
 * \param weight is the weight that will be summed
 * \param factor is the factor the weight will be multiplied with
 * \param target on output will contain the original target plus the weight multiplied
 * by the specified factor
 */
template<typename weight_t, typename value_t>
template<typename W>
void DiscreteStencilWeightManager<weight_t, value_t>::sum(const W &weight, double factor, W *target) const
{
    *target += factor * weight;
}

/*!
 * Sum the specified weight to the target.
 *
 * \param weight is the weight that will be summed
 * \param factor is the factor the weight will be multiplied with
 * \param target on output will contain the original target plus the weight multiplied
 * by the specified factor
 */
template<typename weight_t, typename value_t>
template<typename W, typename V, std::size_t D, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type *>
void DiscreteStencilWeightManager<weight_t, value_t>::sum(const std::array<V, D> &weight, double factor, std::array<V, D> *target) const
{
    for (std::size_t i = 0; i < D; ++i) {
        (*target)[i] += factor * weight[i];
    }
}

/*!
 * Sum the specified weight to the target.
 *
 * The target will be resized to match the size of the weight to be summed. If the weight size is
 * greater that the target size, missing target elements will be initialized to zero before
 * summing the specified weight.
 *
 * \param weight is the weight that will be summed
 * \param factor is the factor the weight will be multiplied with
 * \param target on output will contain the original target plus the weight multiplied
 * by the specified factor
 */
template<typename weight_t, typename value_t>
template<typename W, typename V, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type *>
void DiscreteStencilWeightManager<weight_t, value_t>::sum(const std::vector<V> &weight, double factor, std::vector<V> *target) const
{
    std::size_t weightSize = weight.size();
    std::size_t targetSize = target->size();
    std::size_t commonSize = std::min(weightSize, targetSize);

    for (std::size_t i = 0; i < commonSize; ++i) {
        (*target)[i] += factor * weight[i];
    }

    if (weightSize > targetSize) {
        target->insert(target->end(), weight.cbegin() + commonSize, weight.cend());

        if (factor != 1.) {
            auto targetBegin = target->begin();
            auto targetEnd   = target->end();
            for (auto itr = targetBegin + commonSize; itr != targetEnd; ++itr) {
                *itr *= factor;
            }
        }
    }
}

/**
 * Copy the specified weight into the target.
 *
 * \param weight is the weight that will be copied
 * \param[out] target on output will contain a copy of the weight
 */
template<typename weight_t, typename value_t>
template<typename W>
void DiscreteStencilWeightManager<weight_t, value_t>::copy(const W &weight, W *target) const
{
    *target = weight;
}

/**
 * Copy the specified weight into the target.
 *
 * \param weight is the weight that will be copied
 * \param[out] target on output will contain a copy of the weight
 */
template<typename weight_t, typename value_t>
template<typename W, typename V, std::size_t D, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type *>
void DiscreteStencilWeightManager<weight_t, value_t>::copy(const std::array<V, D> &weight, std::array<V, D> *target) const
{
    std::copy_n(weight.data(), D, target->data());
}

/**
 * Copy the specified weight into the target.
 *
 * The target will be resized to match the size of the specified weight.
 *
 * \param weight is the weight that will be copied
 * \param[out] target on output will contain a copy of the weight
 */
template<typename weight_t, typename value_t>
template<typename W, typename V, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type *>
void DiscreteStencilWeightManager<weight_t, value_t>::copy(const std::vector<V> &weight, std::vector<V> *target) const
{
    std::size_t weightSize = weight.size();
    std::size_t targetSize = target->size();
    std::size_t commonSize = std::min(weightSize, targetSize);

    std::copy_n(weight.data(), commonSize, target->data());

    if (weightSize < targetSize) {
        target->resize(weightSize);
    } else if (weightSize > targetSize) {
        target->insert(target->end(), weight.begin() + commonSize, weight.end());
    }
}

/**
 * Move the specified weight into the target.
 *
 * \param weight is the weight that will be moved
 * \param[out] target on output will contain the moved weight
 */
template<typename weight_t, typename value_t>
template<typename W>
void DiscreteStencilWeightManager<weight_t, value_t>::move(W &&weight, W *target) const
{
    *target = std::move(weight);
}

/*!
 * Get the specified value.
 *
 * \param index is the position of the requested value
 * \result A reference to the specified value.
 */
template<typename weight_t, typename value_t>
template<typename W>
value_t & DiscreteStencilWeightManager<weight_t, value_t>::at(const W &weight, std::size_t index)
{
    return const_cast<value_t &>(const_cast<const DiscreteStencilWeightManager<weight_t, value_t> *>(this)->at(weight, index));
}

/*!
 * Get the specified value.
 *
 * \param index is the position of the requested value
 * \result A reference to the specified value.
 */
template<typename weight_t, typename value_t>
template<typename W>
const value_t & DiscreteStencilWeightManager<weight_t, value_t>::at(const W &weight, std::size_t index) const
{
    BITPIT_UNUSED(index);

    return weight;
}

/*!
 * Get the specified value.
 *
 * \param index is the position of the requested value
 * \result A reference to the specified value.
 */
template<typename weight_t, typename value_t>
template<typename W, typename V, std::size_t D, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type *>
value_t & DiscreteStencilWeightManager<weight_t, value_t>::at(const std::array<V, D> &weight, std::size_t index)
{
    return const_cast<value_t &>(const_cast<const DiscreteStencilWeightManager<weight_t, value_t> *>(this)->at(weight, index));
}

/*!
 * Get the specified value.
 *
 * \param index is the position of the requested value
 * \result A reference to the specified value.
 */
template<typename weight_t, typename value_t>
template<typename W, typename V, std::size_t D, typename std::enable_if<std::is_same<std::array<V, D>, W>::value>::type *>
const value_t & DiscreteStencilWeightManager<weight_t, value_t>::at(const std::array<V, D> &weight, std::size_t index) const
{
    return weight[index];
}

/*!
 * Get the specified value.
 *
 * \param index is the position of the requested value
 * \result A reference to the specified value.
 */
template<typename weight_t, typename value_t>
template<typename W, typename V, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type *>
value_t & DiscreteStencilWeightManager<weight_t, value_t>::at(const std::vector<V> &weight, std::size_t index)
{
    return const_cast<value_t &>(const_cast<const DiscreteStencilWeightManager<weight_t, value_t> *>(this)->at(weight, index));
}

/*!
 * Get the specified value.
 *
 * \param index is the position of the requested value
 * \result A reference to the specified value.
 */
template<typename weight_t, typename value_t>
template<typename W, typename V, typename std::enable_if<std::is_same<std::vector<V>, W>::value>::type *>
const value_t & DiscreteStencilWeightManager<weight_t, value_t>::at(const std::vector<V> &weight, int index) const
{
    return weight[index];
}

}

#endif
