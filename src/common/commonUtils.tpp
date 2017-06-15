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

#ifndef __BITPIT_COMMON_UTILS_TPP__
#define __BITPIT_COMMON_UTILS_TPP__

/*! \file */

namespace bitpit {

namespace utils {

/*!
* Adds an id to an ordered list of unique ids.
*
* \tparam T is the type of elements contained in the list
* \tparam Comparator is the type of the binary function used for the
* comparison of the elements
*
* \param value is the value to be added
* \param list is the ordered list
* \param comparator is a binary function that accepts two arguments (the first
* of the type pointed by ForwardIterator, and the second, always val), and
* returns a value convertible to bool. The value returned indicates whether
* the first argument is considered to go before the second. The function
* shall not modify any of its arguments. This can either be a function
* pointer or a function object.
* \result Returns true is the id was added to the list, false otherwise.
*/
template <typename T, typename Comparator>
bool addToOrderedVector(const T &value, std::vector<T> &list, Comparator comparator)
{
    if (list.empty()) {
        list.push_back(value);
        return true;
    }

    typename std::vector<T>::iterator itr = std::lower_bound(list.begin(), list.end(), value, comparator);
    if (itr == list.end() || *itr != value) {
        list.insert(itr, value);
        return true;
    } else {
        return false;
    }
};

/*!
* Search a value in an ordered list of unique ids.
*
* \tparam T is the type of elements contained in the list
* \tparam Comparator is the type of the binary function used for the
* comparison of the elements
*
* \param value is the value to be searched for
* \param list is the ordered list
* \param comparator is a binary function that accepts two arguments (the first
* of the type pointed by ForwardIterator, and the second, always val), and
* returns a value convertible to bool. The value returned indicates whether
* the first argument is considered to go before the second. The function
* shall not modify any of its arguments. This can either be a function
* pointer or a function object.
* \result Returns true is the value is in the list, false otherwise.
*/
template <typename T, typename Comparator>
typename std::vector<T>::const_iterator findInOrderedVector(const T &value, const std::vector<T> &list, Comparator comparator)
{
    typename std::vector<T>::const_iterator itr = std::lower_bound(list.begin(), list.end(), value, comparator);
    if (itr == list.end() || *itr != value) {
        return list.end();
    }

    return itr;
};

/*!
* Order a vector according to a reordering vector.
*
* \tparam order_t is the type of data that needs to be reordered
* \param order is a reference to the reording vector
* \param v is a reference to the vector that will be reordered
* \param size is the size of the vector that will be reordered
*/
template<typename T>
void reorderVector(std::vector<size_t>& order, std::vector<T>& v, const size_t &size)
{
    for (size_t i = 0; i < size; i++) {
        size_t j;
        while (i != (j = order[i])) {
            size_t k = order[j];

            T temp = std::move(v[j]);
            v[j] = std::move(v[k]);
            v[k] = std::move(temp);

            std::swap(order[i], order[j]);
        }
    }
}

/*!
* Remove a element with specified value from input std::vector.
*
* \param[in,out] vec is the inptut vector set. On output stores the input
* vector deprived by the element with specified value
* \param[in] value is the value of element to be removed from the input vector
*/
template<class T>
void eraseValue(std::vector<T> &vec, const T &value)
{
    typename std::vector<T>::iterator it = find(vec.begin(), vec.end(), value);
    if (it != vec.end()) {
        vec.erase(it);
    }
}

/*!
* Compute intersection between two std::vectors.
*
* \param[in] vec_1 is the 1st argument for intersection.
* \param[in] vec_2 is the 2nd argument for intersection.
* \result Returns a vector storing the common elements of the input vectors.
*/
template<class T>
std::vector<T> intersectionVector(const std::vector<T> &vec_1, const std::vector<T> &vec_2)
{
    std::vector<T> intersect;
    std::map<T, bool> storage;

    intersect.reserve(std::min(vec_1.size(), vec_2.size()));
    for (auto cit_ = vec_2.begin(); cit_ != vec_2.end(); ++cit_) {
        storage[*cit_] = true;
    }

    for (auto it_ = vec_1.begin(); it_ != vec_1.end(); ++it_) {
        if (storage[*it_]) {
            intersect.push_back(*it_);
        }
    }

    return intersect;
}

/*!
* Check if a type is iterable
*/
template <typename T>
auto is_iterable_impl(int) -> decltype (
    std::begin(std::declval<T&>()) != std::end(std::declval<T&>()),   // begin/end and operator !=
    ++std::declval<decltype(std::begin(std::declval<T&>()))&>(),      // operator ++
    *std::begin(std::declval<T&>()),                                  // operator*
    std::true_type{}
);

template <typename T>
std::false_type is_iterable_impl(...);

template <typename T>
using is_iterable = decltype(is_iterable_impl<T>(0));

}

}

#endif
