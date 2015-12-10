//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_COLLAPSED_VECTOR_2D_HPP__
#define __PATCHMAN_COLLAPSED_VECTOR_2D_HPP__

#include <vector>
#include <cassert>
#include <memory>

#include<iostream>

namespace pman {

/*!
	@brief Metafunction for generation of a collapsed vector of arrays.

	@details
	Usage: Use <tt>CollapsedVector2D<Type></tt> to declare a
	collapsed vector of arrays.

	@tparam T The type of the objects stored in the vector
*/

template <class T>
class CollapsedVector2D
{

public:

	/*!
		Default constructor
	*/
	CollapsedVector2D()
	{
		clear();
	}

	/*!
		Creates a new CollapsedVector2D

		\param subArraySizes is a vector with the sizes of the sub-array
		to create
		\param value is the value that will be use to initialize the
		element of the sub-arrays
	*/
	CollapsedVector2D(std::vector<int> &subArraySizes, const T &value = T())
	{
		// Initialize the vector
		clear();

		// Reserve capacity
		int nArrays = subArraySizes.size();

		int dataCapacity = 0;
		for (int i = 0; i < nArrays; i++) {
			dataCapacity += subArraySizes[i];
		}

		reserve(nArrays, dataCapacity);

		// Fill the vector
		for (int i = 0; i < nArrays; i++) {
			push_back(subArraySizes[i], value);
		}
	}

	/*!
		Creates a new CollapsedVector2D

		\param vector2D is a 2D vector that will be used to initialize
		the newly created container
	*/
	CollapsedVector2D(const std::vector<std::vector<T> > &vector2D)
	{
		// Initialize the vector
		clear();

		// Reserve capacity
		int nArrays = vector2D.size();

		int dataCapacity = 0;
		for (int i = 0; i < nArrays; i++) {
			dataCapacity += vector2D[i].size();
		}

		reserve(nArrays, dataCapacity);

		// Fill the vector
		for (int i = 0; i < nArrays; i++) {
			push_back(vector2D[i]);
		}
	}

	/*!
		Copy constructor
	*/
	CollapsedVector2D(const CollapsedVector2D &other)
	{
		// Copy the elements
		std::vector<T> new_v(other.m_v);
		std::vector<std::size_t> new_index(other.m_index);

		// Assign the new memory to the object
		m_v.swap(new_v);
		m_index.swap(new_index);
	}

	/*!
		Destructor
	*/
	~CollapsedVector2D()
	{
		clear();
	}

	/*!
		Assigns content

		Assigns new contents to the container, replacing its current
		contents, and modifying its size accordingly.
	*/
	CollapsedVector2D & operator= (CollapsedVector2D other)
	{
		if (this != &other) {
			other.swap(*this);
		}

		return *this;
	}

	/*!
		Requests a change in capacity

		Requests that the collpased-vector capacity be at least enough
		to contain nSubArrays sub-arrays and nElemnts elements.

		\param nSubArrays is the minimum number of sub-arrays that the
		collapsed-vector should be able to contain
		\param nArrays is the minimum number of elements that the
		collapsed-vector should be able to contain
	*/
	void reserve(int nSubArrays, int nElemnts = 0)
	{
		m_index.reserve(nSubArrays + 1);
		if (nElemnts > 0) {
			m_v.reserve(nElemnts);
		}
	}

	/*!
		Swaps the contents

		\param other is another collapsed-vector container of the same type
	*/
	void swap(CollapsedVector2D &other)
	{
		m_index.swap(other.m_index);
		m_v.swap(other.m_v);
	}

	/*!
		Tests whether two collapsed-vectors are equal

		\result true if the collapsed-vectors are equal, false otherwise.
	*/
	bool operator==(const CollapsedVector2D& rhs) const
	{
		return m_index == rhs.m_index && m_v == rhs.m_v;
	}

	/*!
		Tests whether collapsed-vector is empty

		\result true if the container size is 0, false otherwise.
	*/
	bool empty() const
	{
		return size() == 0;
	}

	/*!
		Clears content

		Removes all elements from the collapsed-vector (which are
		destroyed), leaving the container with a size of 0.
	*/
	void clear()
	{
		m_v.clear();

		m_index.clear();
		m_index.emplace_back(0);
	}

	/*!
		Shrinks to fit

		Requests the container to reduce its capacity to fit its size.
	*/
	void shrink_to_fit()
	{
		m_v.shrink_to_fit();
		m_index.shrink_to_fit();
	}

	/*!
		Returns a direct pointer to the memory vector used internally
		by the container to store its elements.

		\result A pointer to the first element in the vector used
		        internally by the container.

	*/
	T * data() noexcept
	{
		return m_v.data();
	}

	/*!
		Returns a constant reference to the vector used internally by the
		container to store its elements.

		\result A constant reference to the vector used internally by the
		container.

	*/
	const std::vector<T> & vector() const
	{
		return m_v;
	}

	/*!
		Adds an empty sub-array at the end

		Adds an empty element at the end of the vector, after its current
		last element.
	*/
	void push_back()
	{
		push_back(0, NULL);
	}

	/*!
		Adds a sub-array with the specified size at the end

		Adds a sub-array with the specified size at the end of the vector,
		after its current last element. The content of value is copied
		(or moved) to the new sub-array.

		\param subArraySize is the size of the sub-array
		\param value is the value to be copied (or moved) to the new
		element
	*/
	void push_back(const int &subArraySize, const T &value)
	{
		std::size_t previousLastIndex = m_index.back();
		m_index.emplace_back();
		std::size_t &lastIndex = m_index.back();
		lastIndex = previousLastIndex + subArraySize;

		m_v.resize(m_v.size() + subArraySize, value);
	}

	/*!
		Adds the specified sub-array at the end

		Adds the specified sub-array at the end of the vector, after its
		current last element.

		\param subArray is the sub-array that will be added
	*/
	void push_back(const std::vector<T> &subArray)
	{
		int subArraySize = subArray.size();

		std::size_t previousLastIndex = m_index.back();
		m_index.emplace_back();
		std::size_t &lastIndex = m_index.back();
		lastIndex = previousLastIndex + subArraySize;

		m_v.reserve(m_v.size() + subArraySize);
		for (int j = 0; j < subArraySize; j++) {
			m_v.emplace_back();
			T &storedValue = m_v.back();
			storedValue = subArray[j];
		}
	}

	/*!
		Adds an element to the last sub-array

		Adds an element at the end of to the last sub-array.

		\param value is the value that will be added
	*/
	void push_back_in_sub_array(const T& value)
	{
		m_index.back()++;

		m_v.emplace_back();
		T &storedValue = m_v.back();
		storedValue = value;
	}

	/*!
		Adds an element to the specified sub-array

		Adds an element at the end of to the specified last sub-array.

		\param i is the index of the sub-array
		\param value is the value that will be added
	*/
	void push_back_in_sub_array(const int &i, const T& value)
	{
		assert(indexValid(i));

		m_v.insert(m_v.begin() + m_index[i+1], value);

		int nIndexes = m_index.size();
		for (int k = i + 1; k < nIndexes; ++k) {
			m_index[k]++;
		}
	}

	/*!
		Deletes last sub-array

		Removes the last sub-array in the collapsed-vector, effectively
		reducing the container size by one.
	*/
	void pop_back()
	{
		if (size() == 0) {
			return;
		}

		m_index.pop_back();
		m_v.resize(m_index.back() + 1);
	}

	/*!
		Deletes last element from last sub-array

		Removes the last element from the last sub-array in the
		collapsed-vector.
	*/
	void pop_back_in_sub_array()
	{
		if (sub_array_size(size() - 1) == 0) {
			return;
		}

		m_index.back()--;
		m_v.resize(m_index.back() + 1);
	}

	/*!
		Deletes last element from specified sub-array

		Removes the last element from the specified sub-array in the
		collapsed-vector.

		\param i is the index of the sub-array
	*/
	void pop_back_in_sub_array(const int &i)
	{
		assert(indexValid(i));

		if (sub_array_size(i) == 0) {
			return;
		}

		m_v.erase(m_v.begin() + m_index[i+1] - 1);

		int nIndexes = m_index.size();
		for (int k = i + 1; k < nIndexes; ++k) {
			m_index[k]--;
		}
	}

	/*!
		Deletes specified sub-array

		Removes from the collapsed-vector the specified sub-array,
		effectively reducing the container size by one.

		\param i is the index of the sub-array
	*/
	void erase(const int &i)
	{
		assert(indexValid(i));

		m_v.erase(m_v.begin() + m_index[i], m_v.begin() + m_index[i+1] - 1);
		m_index.erase(m_index.begin() + i + 1);
	}

	/*!
		Deletes the specified element from a sub-array

		\param i is the index of the sub-array
		\param j is the index of the element that will be removed
	*/
	void erase(const int &i, const int &j)
	{
		assert(indexValid(i, j));

		m_v.erase(m_v.begin() + m_index[i] + j);

		int nIndexes = m_index.size();
		for (int k = i + 1; k < nIndexes; ++k) {
			m_index[k]--;
		}
	}

	/*!
		Sets the value of the specified element in a sub-array

		\param i is the index of the sub-array
		\param j is the index of the element that will be removed
		\param value is the value that will be set
	*/
	void set(const int &i, const int &j, T value)
	{
		assert(indexValid(i, j));
		(*this)[i][j] = value;
	}

	/*!
		Gets a reference of the specified element in a sub-array

		\param i is the index of the sub-array
		\param j is the index of the element that will be removed
		\return A reference to the requested value.
	*/
	T & get(const int &i, const int &j)
	{
		assert(indexValid(i, j));
		return (*this)[i][j];
	}

	/*!
		Gets a constant reference of the specified element in a sub-array

		\param i is the index of the sub-array
		\param j is the index of the element that will be removed
		\return A constant reference to the requested value.
	*/
	const T & get(const int &i, const int &j) const
	{
		assert(indexValid(i, j));
		return (*this)[i][j];
	}

	/*!
		Gets a pointer to the first element of the specified sub-array

		\param i is the index of the sub-array
		\return A pointer to the first element of the sub-array.
	*/
	const T * get(const int &i) const
	{
		assert(!empty());
		assert(indexValid(i));
		return (*this)[i];
	}

	/*!
		Gets a pointer to the first element of the last sub-array

		\return A pointer to the first element of the sub-array.
	*/
	T * back()
	{
		return get(size() - 1);
	}

	/*!
		Gets a pointer to the first element of the first sub-array

		\return A pointer to the first element of the sub-array.
	*/
	T * first()
	{
		return get(0);
	}

	/*!
		Returns the number of sub-arrays in the collapsed-vector

		\return The number of sub-arrays in the collapsed-vector.
	*/
	int size() const
	{
		return (m_index.size() - 1);
	}

	/*!
		Returns the total size of all the sub-arrays

		\return The total size of all the sub-arrays.
	*/
	int sub_arrays_total_size() const
	{
		return m_index[size()];
	}

	/*!
		Returns the size of the specified sub-array

		\param i is the index of the sub-array
		\return The size of the sub-array.
	*/
	int sub_array_size(int i) const
	{
		return m_index[i + 1] - m_index[i];
	}

	/*!
		Returns the size of the storage space currently allocated for
		storing sub-arrays, expressed in terms of elements.

		\return The size of the storage space currently allocated for
		storing sub-arrays, expressed in terms of elements.
	*/
	int capacity() const
	{
		return m_index.capacity();
	}

	/*!
		Returns the size of the storage space currently allocated for
		storing sub-arrays elements, expressed in terms of elements.

		\return The size of the storage space currently allocated for
		storing sub-arrays elements, expressed in terms of elements.
	*/
	int sub_array_capacity() const
	{
		return m_v.capacity();
	}

private:
	std::vector<T> m_v;
	std::vector<std::size_t> m_index;

	/*!
		Returns a constant pointer to the first element of the specified
		sub-array.

		\param i is the index of the sub-array
		\return A constant pointer to the first element of the specified
		sub-array
	*/
	const T* operator[](const int &i) const
	{
		assert(indexValid(i));

		int index = m_index[i];
		return &m_v[index];
	}

	/*!
		Returns a pointer to the first element of the specified
		sub-array.

		\param i is the index of the sub-array
		\return A pointer to the first element of the specified
		sub-array
	*/
	T* operator[](const int &i)
	{
		assert(indexValid(i));

		int index = m_index[i];
		return &m_v[index];
	}

	/*!
		Checks if the specified index is valid

		\param i is the index of the sub-array
		\return true if the index is vaid, false otherwise.
	*/
	bool indexValid(const int &i) const
	{
		return (i >= 0 && i < size());
	}

	/*!
		Checks if the specified indexes are valid

		\param i is the index of the sub-array
		\param j is the index of the element in the sub-array
		\return true if the indexes are vaid, false otherwise.
	*/
	bool indexValid(const int &i, const int &j) const
	{
		if (!indexValid(i)) {
			return false;
		}

		return (j >= 0 && j < (int) (m_index[i+1] - m_index[i]));
	}

};

}

#endif
