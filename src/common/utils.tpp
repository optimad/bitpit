//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __BITP_MESH_UTILS_TPP__
#define __BITP_MESH_UTILS_TPP__

/*! \file */

#include <algorithm>
#include <vector>

namespace utils {

/*!
	Adds an id to an ordered list of unique ids.

	\tparam T is the type of elements contained in the list
	\tparam Comparator is the type of the binary function used for the
	comparison of the elements

	\param value is the value to be added
	\param list is the ordered list
	\param comparator is a binary function that accepts two arguments (the first
	of the type pointed by ForwardIterator, and the second, always val), and
	returns a value convertible to bool. The value returned indicates whether
	the first argument is considered to go before the second. The function
	shall not modify any of its arguments. This can either be a function
	pointer or a function object.
	\result Returns true is the id was added to the list, false otherwise.
*/
template <typename T, typename Comparator>
bool add_to_ordered_vector(const T &value, std::vector<T> &list, Comparator comparator)
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

}

#endif
