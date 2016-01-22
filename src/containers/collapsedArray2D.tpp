//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __BITPIT_COLLAPSED_ARRAY_2D_HPP__
#define __BITPIT_COLLAPSED_ARRAY_2D_HPP__

#include <vector>
#include <cassert>
#include <memory>

#include<iostream>

/*!
	\ingroup containers
	@{
*/

/*!
	@brief Metafunction for generation of a collapsed array of arrays.

	@details
	Usage: Use <tt>CollapsedArray2D<Type></tt> to declare a
	collapsed array of arrays.

	@tparam T The type of the objects stored in the array
*/

template <class T>
class CollapsedArray2D
{

public:

	/*!
		Default constructor
	*/
	CollapsedArray2D()
		: m_capacity(-1)
	{
	}

	/*!
		Destructor
	*/
	~CollapsedArray2D()
	{
		clear();
	}

	/*!
		Creates a new CollapsedArray2D

		\param nArrays is the number of sub-arrays that the will be
		stored in the container
		\param dataCapacity is the total capacity, expressed in number
		of elements, of the container
	*/
	CollapsedArray2D(const int &nArrays, const int &dataCapacity)
	{
		initialize(nArrays, dataCapacity);
	}

	/*!
		Creates a new CollapsedArray2D

		\param nArrays is the number of sub-arrays that the will be
		stored in the container
		\param subArraySize is the size of each sub-array stored in
		the container
	*/
	CollapsedArray2D(const int &nArrays, const int subArraySize[])
	{
		int dataCapacity = 0;
		for (int i = 0; i < nArrays; i++) {
			dataCapacity += subArraySize[i];
		}

		initialize(nArrays, dataCapacity);

		for (int i = 0; i < nArrays; i++) {
			push_back(subArraySize[i]);
		}
	}

	/*!
		Creates a new CollapsedArray2D

		\param buildFrom is a 2D vector that will be used to initialize
		the newly created container
	*/
	CollapsedArray2D(const std::vector<std::vector<T> > &buildFrom)
	{
		// Initialize the array
		int nArrays = buildFrom.size();
		m_index = std::unique_ptr<int[]>(new int[nArrays + 1]);

		int dataCapacity = 0;
		for (int i = 0; i < nArrays; i++) {
			dataCapacity += buildFrom[i].size();
		}

		initialize(nArrays, dataCapacity);

		// Fill the array
		for (int i = 0; i < nArrays; i++) {
			push_back(buildFrom[i].size());
			for(unsigned int j = 0; j < buildFrom[i].size(); j++) {
				(*this)[i][j] = buildFrom[i][j];
			}
		}
	}

	/*!
		Copy constructor

		\param other is the CollapsedArray2D from which the data will
		be copied from
	*/
	CollapsedArray2D(const CollapsedArray2D &other)
	{
		m_capacity = other.m_capacity;
		if (!other.initialized()) {
			return;
		}

		// Allocate new memory and copy the elements
		int v_size = other.m_index[m_capacity] + 1;

		std::unique_ptr<T[]> new_v = std::unique_ptr<T[]>(new T[v_size]);
		std::copy(other.m_v.get(), other.m_v.get() + v_size - 1, new_v.get());

		std::unique_ptr<int[]> new_index = std::unique_ptr<T[]>(new int[m_capacity + 1]);
		std::copy(other.m_index.get(), other.m_index.get() + m_capacity - 1, new_index.get());

		// Assign the new memory to the object
		m_v = new_v;
		m_index = new_index;
	}

	/*!
		Copy assignment operator

		Assigns new contents to the container, replacing its current
		contents, and modifying its size accordingly.
	*/
	CollapsedArray2D & operator= (CollapsedArray2D other)
	{
		if (this != &other) {
			other.swap(*this);
		}

		return *this;
	}

	/*!
		Swaps the contents

		\param other is another collapsed-array container of the same type
	*/
	void swap(CollapsedArray2D &other)
	{
		std::swap(m_index.get(), other.m_index.get());
		std::swap(m_v.get(), other.m_v.get());
		std::swap(m_capacity, other.m_capacity);
	}

	/*!
		Initializes the container

		\param nArrays is the number of sub-arrays that the will be
		stored in the container
		\param dataCapacity is the total capacity, expressed in number
		of elements, of the container
	*/
	void initialize(int nArrays, int dataCapacity)
	{
		m_capacity = nArrays;

		m_v = std::unique_ptr<T[]>(new T[dataCapacity]);

		m_index = std::unique_ptr<int[]>(new int[m_capacity + 1]);
		m_index[0] = 0;
		m_index[m_capacity] = - dataCapacity;
		std::fill_n(m_index.get() + 1, m_capacity - 1, -1);
	}

	/*!
		Tests whether two collapsed-arrays are equal

		\result true if the collapsed-arrays are equal, false otherwise.
	*/
	bool operator==(const CollapsedArray2D& rhs) const
	{
		return m_index == rhs.m_index && m_v == rhs.m_v;
	}

	/*!
		Tests whether the collapsed-array is empty

		\result true if the container size is 0, false otherwise.
	*/
	bool empty() const
	{
		return size() == 0;
	}

	/*!
		Tests whether the collapsed-array is initialized

		\result true if the container is inizializes, false otherwise.
	*/
	bool initialized() const
	{
		return m_capacity >= 0;
	}

	/*!
		Clears content

		Removes all elements from the collapsed-array (which are
		destroyed), leaving the container with a size of 0.
	*/
	void clear()
	{
		if (!initialized()) {
			return;
		}

		m_v.reset();
		m_index.reset();
		m_capacity = -1;
	}

	/*!
		Returns a direct pointer to the memory array used internally
		by the container to store its elements.

		\result A pointer to the first element in the array used
		internally by the container.

	*/
	T * data() noexcept
	{
		return m_v.get();
	}

	/*!
		Adds a sub-array with the specified size at the end

		Adds a sub-array with the specified size at the end of the vector,
		after its current last element. The specified content is copied
		to the new sub-array.

		\param subArraySize is the size of the sub-array
		\param subArray is the content to be copied to the new sub-array
	*/
	void push_back(const int &subArraySize, T * subArray = NULL)
	{
		int currentSize = size();

		assert(data_size() + subArraySize <= data_capacity());
		m_index[currentSize + 1] = m_index[currentSize] + subArraySize;

		if (subArray != NULL) {
			for (int j = 0; j < subArraySize; j++) {
				set(currentSize, j, subArray[j]);
			}
		}
	}

	/*!
		Adds an empty sub-array at the end

		Adds an empty sub-array at the end of the vector, after its
		current last sub-array.
	*/
	void push_back()
	{
		push_back(0);
	}

	/*!
		Adds an element to the last sub-array

		Adds an element at the end of to the last sub-array.

		\param value is the value that will be added
	*/
	void push_back_in_sub_array(const T& value)
	{
		assert(!empty());

		int &lastIndex = m_index[size()];
		lastIndex++;
		m_v[lastIndex] = value;
	}

	/*!
		Deletes last element from last sub-array

		Removes the last element from the last sub-array in the
		collapsed-array.
	*/
	void pop_back_in_sub_array()
	{
		assert(!empty());
		assert(subArraySize(size() - 1) > 0);

		int &lastIndex = m_index[size()];
		delete m_v[lastIndex];
		lastIndex--;
	}

	/*!
		Gets a constant copy of the specified element in a sub-array

		\param i is the index of the sub-array
		\param j is the index of the element that will be removed
		\return A constant copy of the requested value.
	*/
	const T get(const int &i, const int &j) const
	{
		assert(indexValid(i, j));
		return (*this)[i][j];
	}

	/*!
		Gets a copy of the specified element in a sub-array

		\param i is the index of the sub-array
		\param j is the index of the element that will be removed
		\return A copy of the requested value.
	*/
	T get(const int &i, const int &j)
	{
		assert(indexValid(i, j));
		return (*this)[i][j];
	}

	/*!
		Gets a pointer to the first element of the specified sub-array

		\param i is the index of the sub-array
		\return A pointer to the first element of the sub-array.
	*/
	T * get(const int &i)
	{
		assert(indexValid(i));
		return (*this)[i];
	}

	/*!
		Sets a value in the container

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
		Sets all the values of the specified sub-array

		\param i is the index of the sub-array
		\param values is a pointer to the values that will be set
	*/
	void set(const int &i, T *values)
	{
		assert(indexValid(i));
		for (int j = 1; j < sub_array_size(i); j++) {
			(*this)[i][j] = values[j];
		}
	}

	/*!
		Gets a reference to the first element of the last sub-array

		\return A reference to the first element of the last sub-array.
	*/
	const T & back() const
	{
		assert(!empty());
		return m_v[m_index[m_capacity - 1]];
	}

	/*!
		Gets a pointer to the first element of the last sub-array

		\return A pointer to the first element of the sub-array.
	*/
	T* back()
	{
		assert(!empty());
		return &m_v[m_index[m_capacity - 1]];
	}

	/*!
		Returns the number of sub-arrays in the collapsed-array

		\return The number of sub-arrays in the collapsed-array.
	*/
	int size() const
	{
		if (!initialized()) {
			return 0;
		}

		int nSubArrays = 0;
		for (int i = 1; i < m_capacity + 1; i++) {
			if (m_index[i] < 0) {
				break;
			}

			nSubArrays++;
		}

		return nSubArrays;
	}

	/*!
		Returns the number of elements stored in the sub-array

		\return The number of elements stored in the sub-array.
	*/
	int data_size() const
	{
		int currentSize = size();

		return m_index[currentSize];
	}

	/*!
		Returns the size of the storage space currently allocated for
		storing sub-arrays, expressed in terms of elements.

		\return The size of the storage space currently allocated for
		storing sub-arrays, expressed in terms of elements.
	*/
	int capacity() const
	{
		return m_capacity;
	}

	/*!
		Returns the size of the storage space currently allocated for
		storing elements, expressed in terms of elements.

		\return The size of the storage space currently allocated for
		storing elements, expressed in terms of elements.
	*/
	int data_capacity() const
	{
		if (m_capacity < 0) {
			return 0;
		}

		return abs(m_index[m_capacity]);
	}

	/*!
		Returns the size of the specified sub-array.

		\param i is the index of the sub-array
		\return The size of the specified sub-array.
	*/
	int sub_array_size(int i) const
	{
		return m_index[i + 1] - m_index[i];
	}

private:
	std::unique_ptr<T[]> m_v;
	std::unique_ptr<int[]> m_index;
	int m_capacity;

	/*!
		Gets a constant pointer to the specified sub-array

		\param i is the index of the sub-array
		\return A constant pointer to the specified sub-array.
	*/
	const T* operator[](const int &i) const
	{
		assert(indexValid(i));

		int index = m_index[i];
		return &m_v[index];
	}

	/*!
		Gets a pointer to the specified sub-array

		\param i is the index of the sub-array
		\return A pointer to the specified sub-array.
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
	bool indexValid(const int &i)
	{
		if (empty()) {
			return false;
		} else if (i < 0 || i >= m_capacity) {
			return false;
		}

		return m_index[i] >= 0;
	}

	/*!
		Checks if the specified indexes are valid

		\param i is the index of the sub-array
		\param j is the index of the element in the sub-array
		\return true if the indexes are vaid, false otherwise.
	*/
	bool indexValid(const int &i, const int &j)
	{
		if (!indexValid(i)) {
			return false;
		}

		return j < (m_index[i+1] - m_index[i]);
	}

};

/*!
	@}
*/

#endif
