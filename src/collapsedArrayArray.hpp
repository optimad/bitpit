//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//

#include <vector>
#include <cassert>
#include <memory>

#include<iostream>

/*!
	@brief Metafunction for generation of a collapsed array of arrays.

	@details
	Usage: Use <tt>CollapsedArrayArray<Type></tt> to declare a
	collapsed array of arrays.

	@tparam T The type of the objects stored in the array
*/

template <class T>
class CollapsedArrayArray
{

public:

	CollapsedArrayArray()
		: m_capacity(-1)
	{
	}

	~CollapsedArrayArray()
	{
		clear();
	}

	CollapsedArrayArray(const int &nArrays, const int &dataCapacity)
	{
		initialize(nArrays, dataCapacity);
	}

	CollapsedArrayArray(const int &nArrays, const int subArraySize[])
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

	CollapsedArrayArray(const std::vector<std::vector<T> > &buildFrom)
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

	CollapsedArrayArray(const CollapsedArrayArray &other)
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

	CollapsedArrayArray & operator= (CollapsedArrayArray other)
	{
		if (this != &other) {
			other.swap(*this);
		}

		return *this;
	}

	void swap(CollapsedArrayArray &other)
	{
		std::swap(m_index.get(), other.m_index.get());
		std::swap(m_v.get(), other.m_v.get());
		std::swap(m_capacity, other.m_capacity);
	}

	void initialize(int nArrays, int dataCapacity)
	{
		m_capacity = nArrays;

		m_v = std::unique_ptr<T[]>(new T[dataCapacity]);

		m_index = std::unique_ptr<int[]>(new int[m_capacity + 1]);
		m_index[0] = 0;
		m_index[m_capacity] = dataCapacity;
		std::fill_n(m_index.get() + 1, m_capacity - 1, -1);
	}

	bool operator==(const CollapsedArrayArray& rhs) const
	{
		return m_index == rhs.m_index && m_v == rhs.m_v;
	}

	bool empty() const
	{
		return size() == 0;
	}

	bool initialized() const
	{
		return m_capacity >= 0;
	}

	void clear()
	{
		if (!initialized()) {
			return;
		}

		m_v.reset();
		m_index.reset();
		m_capacity = -1;
	}

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

	void push_back()
	{
		push_back(0);
	}

	void push_back_in_sub_array(const T& value)
	{
		assert(!empty());

		int &lastIndex = m_index[size()];
		lastIndex++;
		m_v[lastIndex] = value;
	}

	void pop_back_in_sub_array()
	{
		assert(!empty());
		assert(subArraySize(size() - 1) > 0);

		int &lastIndex = m_index[size()];
		delete m_v[lastIndex];
		lastIndex--;
	}

	const T get(const int &i, const int &j) const
	{
		assert(indexValid(i, j));
		return (*this)[i][j];
	}

	T get(const int &i, const int &j)
	{
		assert(indexValid(i, j));
		return (*this)[i][j];
	}

	T * get(const int &i)
	{
		assert(indexValid(i));
		return (*this)[i];
	}

	void set(const int &i, const int &j, T value)
	{
		assert(indexValid(i, j));
		(*this)[i][j] = value;
	}

	void set(const int &i, T *values)
	{
		assert(indexValid(i));
		for (int j = 1; j < sub_array_size(i); j++) {
			(*this)[i][j] = values[j];
		}
	}


	const T & back() const
	{
		assert(!empty());
		return m_v[m_index[m_capacity - 1]];
	}

	T* back()
	{
		assert(!empty());
		return &m_v[m_index[m_capacity - 1]];
	}

	int size() const
	{
		if (!initialized()) {
			return 0;
		}

		int nSubArrays = 0;
		for (int i = 1; i < m_capacity; i++) {
			if (m_index[i] == -1) {
				break;
			}

			nSubArrays++;
		}

		return nSubArrays;
	}

	int data_size() const
	{
		int currentSize = size();

		return m_index[currentSize];
	}

	int capacity() const
	{
		return m_capacity;
	}

	int data_capacity() const
	{
		if (m_capacity < 0) {
			return 0;
		}

		return (m_index[m_capacity]);
	}

	int sub_array_size(int i) const
	{
		return m_index[i + 1] - m_index[i];
	}

private:
	std::unique_ptr<T[]> m_v;
	std::unique_ptr<int[]> m_index;
	int m_capacity;

	const T* operator[](const int &i) const
	{
		assert(indexValid(i));

		int index = m_index[i];
		return &m_v[index];
	}

	T* operator[](const int &i)
	{
		assert(indexValid(i));

		int index = m_index[i];
		return &m_v[index];
	}

	bool indexValid(const int &i)
	{
		if (empty()) {
			return false;
		} else if (i < 0 || i >= m_capacity) {
			return false;
		}

		return m_index[i] >= 0;
	}

	bool indexValid(const int &i, const int &j)
	{
		if (!indexValid(i)) {
			return false;
		}

		return j < (m_index[i+1] - m_index[i]);
	}

};
