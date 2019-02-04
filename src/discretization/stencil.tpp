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
    buffer << stencil.m_pattern;
    buffer << stencil.m_weights;
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
    buffer >> stencil.m_pattern;
    buffer >> stencil.m_weights;
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
    : DiscreteStencil(1, zero)
{
}

/*!
* Constructor
*
* \param nBuckets is the number of buckets in the stencil
* \param zero is the value to be used as zero
*/
template<typename weight_t>
DiscreteStencil<weight_t>::DiscreteStencil(int nBuckets, const weight_t &zero)
    : DiscreteStencil(nBuckets, 0, zero)
{
}

/*!
* Constructor
*
* \param nBuckets is the number of buckets in the stencil
* \param nBucketItems is the number of items contained in each bucket
* \param zero is the value to be used as zero
*/
template<typename weight_t>
DiscreteStencil<weight_t>::DiscreteStencil(int nBuckets, int nBucketItems, const weight_t &zero)
    : m_zero(zero),
      m_pattern(nBuckets, nBucketItems, NULL_ID), m_weights(nBuckets, nBucketItems, zero),
      m_constant(m_zero)
{
}

/*!
* Constructor
*
* \param bucketSizes are the sizes of the buckets
* \param zero is the value to be used as zero
*/
template<typename weight_t>
DiscreteStencil<weight_t>::DiscreteStencil(const std::vector<int> &bucketSizes, const weight_t &zero)
    : m_zero(zero),
      m_pattern(bucketSizes, NULL_ID), m_weights(bucketSizes, zero),
      m_constant(m_zero)
{
}

/*!
* Initialize the stencil
*
* \param zero is the value to be used as zero
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::initialize(const weight_t &zero)
{
    initialize(1, zero);
}

/*!
* Initialize the stencil
*
* \param nBuckets is the number of buckets in the stencil
* \param zero is the value to be used as zero
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::initialize(int nBuckets, const weight_t &zero)
{
    initialize(nBuckets, 0, zero);
}

/*!
* Initialize the stencil
*
* \param nBuckets is the number of buckets in the stencil
* \param nBucketItems is the number of items contained in each bucket
* \param zero is the value to be used as zero
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::initialize(int nBuckets, int nBucketItems, const weight_t &zero)
{
    m_zero = zero;
    m_pattern.initialize(nBuckets, nBucketItems, NULL_ID);
    m_weights.initialize(nBuckets, nBucketItems, zero);
    m_constant = m_zero;
}

/*!
* Initialize the stencil
*
* \param bucketSizes are the sizes of the buckets in the stencil
* \param zero is the value to be used as zero
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::initialize(const std::vector<int> &bucketSizes, const weight_t &zero)
{
    m_zero = zero;
    m_pattern.initialize(bucketSizes, NULL_ID);
    m_weights.initialize(bucketSizes, zero);
    m_constant = m_zero;
}

/*!
* Initialize the stencil
*
* \param other is another stencil of the same time, whose items will be used
* to initialize this stencil
* \param zero is the value to be used as zero
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::initialize(const DiscreteStencil<weight_t> &other)
{
    m_zero = other.m_zero;
    m_pattern.initialize(other.m_pattern);
    m_weights.initialize(other.m_weights);
    m_constant = other.m_zero;
}

/*!
* Get the total size of the stencil, expressed in number of items.
*
* \result The total size of the stencil, expressed in number of items.
*/
template<typename weight_t>
std::size_t DiscreteStencil<weight_t>::size() const
{
    return m_pattern.getItemCount();
}

/*!
* Get the size of the specified bucket, expressed in number of items.
*
* \param bucket is the bucket
* \result The size of the specified bucket, expressed in number of items.
*/
template<typename weight_t>
std::size_t DiscreteStencil<weight_t>::size(int bucket) const
{
    return m_pattern.getItemCount(bucket);
}

/*!
* Requests a change in capacity.
*
* Requests that the stencil capacity be at least nItems items.
*
* \param nItems is the minimum number of items that the stencil should
* be able to contain
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::reserve(std::size_t nItems)
{
    return reserve(1, nItems);
}

/*!
* Requests a change in capacity.
*
* Requests that the stencil capacity be at least enough to contain nBuckets
* buckets and nItems items.
*
* \param nBuckets is the minimum number of buckets that the stencil
* should be able to contain
* \param nItems is the minimum number of items that the stencil should
* be able to contain
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::reserve(int nBuckets, std::size_t nItems)
{
    return m_pattern.reserve(nBuckets, nItems);
}

/*!
* Get the number of buckets of the stencil.
*
* \result The number of buckets of the stencil.
*/
template<typename weight_t>
int DiscreteStencil<weight_t>::getBucketCount() const
{
    return m_pattern.size();
}

/*!
* Get a reference to the specified element of the pattern.
*
* If the stencil has more than one bucket, the first bucket will be considered.
*
* \param pos is the position of the pattern element
* \result A reference to the specified element of the pattern.
*/
template<typename weight_t>
long & DiscreteStencil<weight_t>::getPattern(std::size_t pos)
{
    return getPattern(0, pos);
}

/*!
* Get a constant reference to the specified element of the pattern.
*
* If the stencil has more than one bucket, the first bucket will be considered.
*
* \param pos is the position of the pattern element
* \result A constant reference to the specified element of the pattern.
*/
template<typename weight_t>
const long & DiscreteStencil<weight_t>::getPattern(std::size_t pos) const
{
    return getPattern(0, pos);
}

/*!
* Get a reference to the specified element of the pattern from the requested
* bucket.
*
* \param bucket is the bucket of the pattern element
* \param pos is the position of the pattern element
* \result A reference to the specified element of the pattern from the
* requested bucket.
*/
template<typename weight_t>
long & DiscreteStencil<weight_t>::getPattern(int bucket, std::size_t pos)
{
    return m_pattern.getItem(bucket, pos);
}

/*!
* Get a constant reference to the specified element of the pattern from the
* requested bucket.
*
* \param bucket is the bucket of the pattern element
* \param pos is the position of the pattern element
* \result A constant reference to the specified element of the pattern from
* the requested bucket.
*/
template<typename weight_t>
const long & DiscreteStencil<weight_t>::getPattern(int bucket, std::size_t pos) const
{
    return m_pattern.getItem(bucket, pos);
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
* Get a constant reference to the container serving as pattern storage.
*
* \result A constant reference to the container serving as pattern storage.
*/
template<typename weight_t>
const FlatVector2D<long> & DiscreteStencil<weight_t>::getPattern() const
{
    return m_pattern;
}

/*!
* Set the index of the specified element of the pattern.
*
* If the stencil has more than one bucket, the first bucket will be considered.
*
* \param pos is the position of the pattern element
* \param id is the index that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::setPattern(std::size_t pos, long id)
{
    setPattern(0, pos, id);
}

/*!
* Set the index of the specified element of the pattern.
*
* \param bucket is the bucket that will updated
* \param pos is the position of the pattern element
* \param id is the index that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::setPattern(int bucket, std::size_t pos, long id)
{
    m_pattern.setItem(bucket, pos, id);
}

/*!
* Get a constant reference to the specified element of the pattern.
*
* If the stencil has more than one bucket, the first bucket will be considered.
*
* \param pos is the raw position of the pattern element
* \result A constant reference to the specified element of the pattern.
*/
template<typename weight_t>
long & DiscreteStencil<weight_t>::rawGetPattern(std::size_t pos)
{
    return m_pattern.rawGetItem(pos);
}

/*!
* Get a constant reference to the specified element of the pattern.
*
* If the stencil has more than one bucket, the first bucket will be considered.
*
* \param pos is the raw position of the pattern element
* \result A constant reference to the specified element of the pattern.
*/
template<typename weight_t>
const long & DiscreteStencil<weight_t>::rawGetPattern(std::size_t pos) const
{
    return m_pattern.rawGetItem(pos);
}

/*!
* Set the index of the specified element of the pattern.
*
* \param pos is the raw position of the pattern element
* \param id is the index that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::rawSetPattern(std::size_t pos, long id)
{
    m_pattern.rawSetItem(pos, id);
}

/*!
* Get a reference to the specified weight of the stencil.
*
* If the stencil has more than one bucket, the first bucket will be considered.
*
* \param pos is the position of the weight
* \result A reference to the specified weight of the stencil.
*/
template<typename weight_t>
weight_t & DiscreteStencil<weight_t>::getWeight(std::size_t pos)
{
    return getWeight(0, pos);
}

/*!
* Get a constant reference to the specified weight of the stencil.
*
* If the stencil has more than one bucket, the first bucket will be considered.
*
* \param pos is the position of the weight
* \result A constant reference to the specified weight of the stencil.
*/
template<typename weight_t>
const weight_t & DiscreteStencil<weight_t>::getWeight(std::size_t pos) const
{
    return getWeight(0, pos);
}

/*!
* Get a reference to the specified weight of the stencil.
*
* \param bucket is the bucket of the weight
* \param pos is the position of the weight
* \result A reference to the specified weight of the stencil.
*/
template<typename weight_t>
weight_t & DiscreteStencil<weight_t>::getWeight(int bucket, std::size_t pos)
{
    return m_weights.getItem(bucket, pos);
}

/*!
* Get a constant reference to the specified weight of the stencil.
*
* \param bucket is the bucket of the weight
* \param pos is the position of the weight
* \result A constant reference to the specified weight of the stencil.
*/
template<typename weight_t>
const weight_t & DiscreteStencil<weight_t>::getWeight(int bucket, std::size_t pos) const
{
    return m_weights.getItem(bucket, pos);
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
* Get a constant reference to the container serving as weight storage.
*
* \result A constant reference to the container serving as weight storage.
*/
template<typename weight_t>
const FlatVector2D<weight_t> & DiscreteStencil<weight_t>::getWeights() const
{
    return m_weights;
}

/*!
* Set the value of the specified weight of the stencil.
*
* If the stencil has more than one bucket, the first bucket will be considered.
*
* \param pos is the position of the weight
* \param weight is the value that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::setWeight(std::size_t pos, const weight_t &weight)
{
    setWeight(0, pos, weight);
}

/*!
* Set the value of the specified weight of the stencil.
*
* \param bucket is the bucket that will updated
* \param pos is the position of the weight
* \param weight is the value that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::setWeight(int bucket, std::size_t pos, const weight_t &weight)
{
    m_weights.setItem(bucket, pos, weight);
}

/*!
* Sum the given value to the weight at the specified position of the stencil.
*
* If the stencil has more than one bucket, the first bucket will be considered.
*
* \param pos is the position of the weight
* \param value is the value that will be summed
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::sumWeight(std::size_t pos, const weight_t &value)
{
    sumWeight(0, pos, value);
}

/*!
* Sum the given value to the weight at the specified position of the stencil.
*
* \param bucket is the bucket that will updated
* \param pos is the position of the weight
* \param value is the value that will be summed
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::sumWeight(int bucket, std::size_t pos, const weight_t &value)
{
    weight_type &weight = m_weights.getItem(bucket, pos);
    weight += value;
}


/*!
* Get a reference to the specified weight of the stencil.
*
* If the stencil has more than one bucket, the first bucket will be considered.
*
* \param pos is the raw position of the weight
* \result A reference to the specified weight of the stencil.
*/
template<typename weight_t>
weight_t & DiscreteStencil<weight_t>::rawGetWeight(std::size_t pos)
{
    return m_weights.rawGetItem(pos);
}

/*!
* Get a constant reference to the specified weight of the stencil.
*
* If the stencil has more than one bucket, the first bucket will be considered.
*
* \param pos is the raw position of the weight
* \result A constant reference to the specified weight of the stencil.
*/
template<typename weight_t>
const weight_t & DiscreteStencil<weight_t>::rawGetWeight(std::size_t pos) const
{
    return m_weights.rawGetItem(pos);
}

/*!
* Set the value of the specified weight of the stencil.
*
* \param pos is the raw position of the weight
* \param weight is the value that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::rawSetWeight(std::size_t pos, const weight_t &weight)
{
    m_weights.rawSetItem(pos, weight);
}

/*!
* Set the specified item of the stencil.
*
* If the stencil has more than one bucket, the first bucket will be considered.
*
* \param pos is the position of the weight
* \param id is the index that will be set
* \param weight is the weight that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::setItem(std::size_t pos, long id, const weight_t &weight)
{
    setItem(0, pos, id, weight);
}

/*!
* Set the specified item of the stencil.
*
* \param bucket is the bucket that will updated
* \param pos is the position of the weight
* \param id is the index that will be set
* \param weight is the weight that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::setItem(int bucket, std::size_t pos, long id, const weight_t &weight)
{
    setPattern(bucket, pos, id);
    setWeight(bucket, pos, weight);
}

/*!
* Sum the given value to the item of the stencil with the specified index.
*
* If an item with the same id already exists, the given value will be summed
* to the weight of the existing item. Otherwise, a new item will be appended.
*
* If the stencil has more than one bucket, the item will be appended to the
* first bucket.
*
* \param id is the index of the item
* \param value is the value that will be summed
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::sumItem(long id, const weight_t &value)
{
    sumItem(0, id, value);
}

/*!
* Sum the given value to the item of the stencil with the specified index.
*
* If an item with the same id already exists, the given value will be summed
* to the weight of the existing item. Otherwise, a new item will be appended.
*
* If the stencil has more than one bucket, the item will be appended to the
* first bucket.
*
* \param bucket is the bucket that will updated
* \param id is the index of the item
* \param value is the value that will be summed
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::sumItem(int bucket, long id, const weight_t &value)
{
    weight_t *weight = findWeight(bucket, id);
    if (weight) {
        *weight += value;
    } else {
        appendItem(bucket, id, value);
    }
}


/*!
* Append an item the stencil.
*
* If the stencil has more than one bucket, the item will be appended to the
* first bucket.
*
* \param id is the index that will be set
* \param weight is the weight that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::appendItem(long id, const weight_t &weight)
{
    appendItem(0, id, weight);
}

/*!
* Append an item the stencil.
*
* \param bucket is the bucket that will updated
* \param id is the index that will be set
* \param weight is the weight that will be set
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::appendItem(int bucket, long id, const weight_t &weight)
{
    m_pattern.pushBackItem(bucket, id);
    m_weights.pushBackItem(bucket, weight);
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
    m_constant = value;
}

/*!
* Sum the specified value to the constant associated to the stencil.
*
* \param value is the value that will be summed
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::sumConstant(const weight_t &value)
{
    m_constant += value;
}

/*!
* Clears the items of the stencil.
*
* Removes all items from the stencil (which are destroyed), leaving each
* bucket with a size of 0.
*
* \param release if it's true the memory hold by the stencil will be
* released, otherwise the stencil will be cleared but its memory will
* not be relased
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::clear(bool release)
{
    m_pattern.clearItems(release);
    m_weights.clearItems(release);

    m_constant = m_zero;
}

/*!
* Merge the bucket together.
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::flatten()
{
    m_pattern.merge();
    m_weights.merge();
}

/*!
* Optimize the stencil.
*
* The negligible elements will be removed from the stencil.
*/
template<typename weight_t>
void DiscreteStencil<weight_t>::optimize(double tolerance)
{
    int nBuckets = getBucketCount();
    for (int i = 0; i < nBuckets; ++i) {
        std::size_t nBucketItems = size(i);
        for (std::size_t j = 0; j < nBucketItems; ++j) {
            if (isWeightNeglibile(i, j, tolerance)) {
                m_pattern.eraseItem(i, j);
                m_weights.eraseItem(i, j);
                --j;
                --nBucketItems;
            }
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
    for (std::size_t k = 0; k < nItems; ++k) {
        long &id = m_pattern.rawGetItem(k);
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

    weight_t complement = -1. * m_weights.rawGetItem(0);
    for (std::size_t n = 1; n < nItems; ++n) {
        complement -= m_weights.rawGetItem(n);
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
        weight_t &weight = m_weights.rawGetItem(n);
        weight = m_zero;
    }

    setConstant(m_zero);
}
/*!
* Check if the specified weight is neglibile accordingly the specified tolerance.
*
* \param bucket is the bucket of the weight to check
* \param pos is the position of the weight to check
* \param tolerance is the tolerance that will be used for the check
*/
template<typename weight_t>
template<typename U, typename std::enable_if<std::is_fundamental<U>::value>::type *>
bool DiscreteStencil<weight_t>::isWeightNeglibile(int bucket, std::size_t pos, double tolerance)
{
    return (std::abs(m_weights.getItem(bucket, pos)) <= tolerance);
}

/*!
* Check if the specified weight is neglibile accordingly the specified tolerance.
*
* \param bucket is the bucket of the weight to check
* \param pos is the position of the weight to check
* \param tolerance is the tolerance that will be used for the check
*/
template<typename weight_t>
template<typename U, typename std::enable_if<!std::is_fundamental<U>::value>::type *>
bool DiscreteStencil<weight_t>::isWeightNeglibile(int bucket, std::size_t pos, double tolerance)
{
    return (norm2(m_weights.getItem(bucket, pos)) <= tolerance);
}

/*!
* Find the weight associated to the specified id.
*
* \param bucket is the bucket of the weight
* \param id is the index associated to the weight
* \result A pointer to the weight associated to the specified id. If the weight
* is not found, a null pointer is returned.
*/
template<typename weight_t>
weight_t * DiscreteStencil<weight_t>::findWeight(int bucket, long id)
{
    std::size_t nBucketItems = size(bucket);
    for (std::size_t j = 0; j < nBucketItems; ++j) {
        long guessId = *(m_pattern.get(bucket) + j);
        if (guessId == id) {
            weight_t *weight = m_weights.get(bucket) + j;
            return weight;
        }
    }

    return nullptr;
}

/*!
* Find the weight associated to the specified id.
*
* \param bucket is the bucket of the weight
* \param id is the index associated to the weight
* \result A constant pointer to the weight associated to the specified id.
* If the weight is not found, a null pointer is returned.
*/
template<typename weight_t>
const weight_t * DiscreteStencil<weight_t>::findWeight(int bucket, long id) const
{
    std::size_t nBucketItems = size(bucket);
    for (std::size_t j = 0; j < nBucketItems; ++j) {
        long guessId = *(m_pattern.get(bucket) + j);
        if (guessId == id) {
            const weight_t *weight = m_weights.get(bucket) + j;
            return weight;
        }
    }

    return nullptr;
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
    int nBuckets = getBucketCount();

    weight_t sum = m_zero;
    for (int i = 0; i < nBuckets; ++i) {
        std::size_t nBucketItems = size(i);
        out << " bucket : " << i << " n. bucket items : " << nBucketItems <<  std::endl;
        for (std::size_t j = 0; j < nBucketItems; ++j) {
            long id = m_pattern.getItem(i, j);
            weight_t value = factor * m_weights.getItem(i, j);
            out << "   id: " << id << " weight: " << value << std::endl;
            sum += value;
        }
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
    return (sizeof(m_zero) + m_pattern.getBinarySize() + m_weights.getBinarySize() + sizeof(m_constant));
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
    for (std::size_t k = 0; k < nItems; ++k) {
        weight_t &weight = m_weights.rawGetItem(k);
        weight *= factor;
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
    for (std::size_t k = 0; k < nItems; ++k) {
        weight_t &weight = m_weights.rawGetItem(k);
        weight /= factor;
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
    int nBuckets = getBucketCount();
    int other_nBuckets = other.getBucketCount();
    if (nBuckets != other_nBuckets) {
        throw std::runtime_error("Stencil must have the same number of buckets.");
    }

    for (int i = 0; i < other_nBuckets; ++i) {
        const int other_nBucketItems = other.size(i);
        for (int j = 0; j < other_nBucketItems; ++j) {
            long id = *(other.m_pattern.get(i) + j);
            const weight_t &other_weight = *(other.m_weights.get(i) + j);
            sumItem(i, id, other_weight);
        }
    }

    m_constant += other.m_constant;

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
    int nBuckets = getBucketCount();
    int other_nBuckets = other.getBucketCount();
    if (nBuckets != other_nBuckets) {
        throw std::runtime_error("Stencil must have the same number of buckets.");
    }

    for (int i = 0; i < other_nBuckets; ++i) {
        const int other_nBucketItems = other.size(i);
        for (int j = 0; j < other_nBucketItems; ++j) {
            long id = *(other.m_pattern.get(i) + j);
            const weight_t &other_weight = *(other.m_weights.get(i) + j);
            sumItem(i, id, -1. * other_weight);
        }
    }

    m_constant -= other.m_constant;

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
