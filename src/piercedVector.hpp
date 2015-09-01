//
// Written by Andrea Iob <andrea_iob@hotmail.com>
//
#ifndef __PATCHMAN_PIERCED_ARRAY_HPP__
#define __PATCHMAN_PIERCED_ARRAY_HPP__

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <deque>
#include <iostream>
#include <iterator>
#include <limits>
#include <unordered_map>
#include <memory>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>

namespace pman {

// To check if the provided template argument implements the needed methods,
// the following Stackoverflow questions are used:
//
// http://stackoverflow.com/questions/257288/is-it-possible-to-write-a-c-template-to-check-for-a-functions-existence
// http://stackoverflow.com/questions/16976720/how-to-i-restrict-a-template-class-to-certain-types

template <typename T>
class has_get_id
{
    typedef char true_type;
    typedef long false_type;

    template <typename C> static true_type test(decltype(&C::get_id));
    template <typename C> static false_type test(...);

public:
    enum { value = (sizeof(test<T>(0)) == sizeof(true_type)) };
};

template <typename T>
class has_set_id
{
    typedef char true_type;
    typedef long false_type;

    template <typename C> static true_type test(decltype(&C::set_id));
    template <typename C> static false_type test(...);

public:
    enum { value = (sizeof(test<T>(0)) == sizeof(true_type)) };
};

/*!
	@brief Iterator for the class PiercedVector

	@details
	PiercedIterator can work only with objects that are identified by a
	unique id. The objects must implement a function, called get_id,
	that returns the id of the object. The ids have to be positive, but
	the type that defines the id must be a signed integer (negative
	id are used by PiercedVector to store special information in the
	elements).

	Usage: Use <tt>PiercedVector<T>::iterator</tt> to declare an iterator
	for a pierced vector, use <tt>PiercedVector<Type>::const_iterator</tt> to
	declare a const iterator for a pierced vector.

	@tparam T The type of the objects stored in the vector
*/

template<class T, class unqualified_T = typename std::remove_cv<T>::type>
class PiercedIterator
	: public std::iterator<std::forward_iterator_tag, unqualified_T,
				std::ptrdiff_t, T*, T&>
{
	// PiercedIterator can work only with calsses that are identified by a
	// unique id (i.e., classes that implements get_id)
	static_assert(has_get_id<T>::value, "Provided class does not implement get_id");

private:
	/*!
		Iterator for the internal vector that holds the elements in
		the pierced array.
	*/
	typedef typename std::vector<unqualified_T>::iterator BaseIterator;

	/*!
		Const iterator for the internal vector that holds the elements in
		the pierced array.
	*/
	typedef typename std::vector<unqualified_T>::const_iterator BaseConstIterator;

	/*!
		Type id_type is the type of the ids.

		It is automatically defined as the type returned by the
		get_id function .
	*/
	typedef decltype(std::declval<unqualified_T>().get_id()) id_type;

	/*!
		Special id value that identifies every dummy element past
		the end of the pierced vector.
	*/
	static const id_type SENTINEL_ID;

	/*!
		Internal iterator.
	*/
	T *m_itr;

	/*!
		Creates a new iterator and initializes it with the
		specified pointer.
	*/
	explicit PiercedIterator(T *value)
		: m_itr(value)
	{
	}

public:

	/*!
		Creates a new uninitialized iterator
	*/
	PiercedIterator()
	{
	}

	/*!
		Creates a new iterator and initializes it with the position of
		the base iterator recevied in input.
	*/
	explicit PiercedIterator(BaseIterator iterator)
		: m_itr(&(*iterator))
	{
	}

	/*!
		Creates a new iterator and initializes it with the position of
		the const base iterator recevied in input.
	*/
	explicit PiercedIterator(BaseConstIterator iterator)
		: m_itr(&(*iterator))
	{
	}

	/*!
		Exchanges the values of the current iterator and
		the iterator recevied as argument.

		\param other the iterator to exchange values with
	*/
	void swap(PiercedIterator& other) noexcept
	{
		using std::swap;
		swap(m_itr, other.m_itr);
	}

	/*!
		Pre-increment operator.
	*/
	PiercedIterator& operator++ ()
	{
		m_itr++;

		id_type id = m_itr->get_id();
		if (id != SENTINEL_ID && id < 0) {
			m_itr += - m_itr->get_id();
		}

		return *this;
	}

	/*!
		Post-increment operator.
	*/
	PiercedIterator operator++ (int)
	{
		PiercedIterator tmp(m_itr);

		++(*this);

		return tmp;
	}

	/*!
		Two-way comparison.
	*/
	template<class other_T>
	bool operator == (const PiercedIterator<other_T>& rhs) const
	{
		return m_itr == rhs.m_itr;
	}

	/*!
		Two-way comparison.
	*/
	template<class other_T>
	bool operator != (const PiercedIterator<other_T>& rhs) const
	{
		return m_itr != rhs.m_itr;
	}

	/*!
		Deference operator.

		\result A reference to the element currently pointed to by the
		        iterator.
	*/
	T& operator* () const
	{
		return *m_itr;
	}

	/*!
		Deference operator.

		\result A reference to the element currently pointed to by the
		        iterator.
	*/
	T& operator-> () const
	{
		return *m_itr;
	}

	/*!
		Assignment operator.

		\param iterator is the base type iterator that holds the
		                position to be set
		\result The updated iterator.
	*/
	PiercedIterator & operator= (BaseIterator iterator)
	{
		m_itr = &(*iterator);

		return *this;
	}

	/*!
		Converts the iterator to a const_iterator.
	*/
	operator PiercedIterator<const T>() const
	{
		return PiercedIterator<const T>(m_itr);
	}
};

// Definition of static constants of PiercedIterator
template<class T, class unqualified_T>
const typename PiercedIterator<T, unqualified_T>::id_type
	PiercedIterator<T, unqualified_T>::SENTINEL_ID = std::numeric_limits<id_type>::min();

/*!
	@brief Metafunction for generating of a pierced vector.

	@details
	Usage: Use <tt>PiercedVector<T></tt> to declare a pierced vector.

	PiercedVector can work only with objects that are identified by a
	unique id. The objects must implement a function, called get_id,
	that returns the id of the object and a function, called set_id,
	that sets the id of the objects. The ids have to be positive, but
	the type that defines the id must be a signed integer (negative
	id are used by PiercedVector to store special information in the
	elements).

	@tparam T The type of the objects stored in the vector
*/

template <class T>
class PiercedVector
{
	// PiercedVector can work only with calsses that are identified by a
	// unique id (i.e., classes that implements set_id and get_id)
	static_assert(has_get_id<T>::value, "Provided class does not implement get_id");
	static_assert(has_set_id<T>::value, "Provided class does not implement set_id");

private:
	/*!
		Type size_type is an unsigned integral type.
	*/
	typedef std::size_t size_type;

	/*!
		Member type value_type is the type of the elements in the
		container, defined as an alias of the first class template
		parameter (T).
	*/
	typedef T value_type;

	/*!
		Type id_type is the type of the ids.

		It is automatically defined as the type returned by the
		get_id function .
	*/
	typedef decltype(std::declval<T>().get_id()) id_type;

	/*!
		Special id value that identifies every dummy element past
		the end of the pierced vector.
	*/
	static const id_type SENTINEL_ID;

	/*!
		Special hole value that can be used with the functions
		that insert new element. It means to insert the element
		at the end of the vector.
	*/
	static const size_type APPEND_TO_BACK;

	/*!
		Special hole value that can be used with the functions that
		insert new element. It means to insert the element
		in the first empty position, i.e., the first hole or
		at the end of the vector if there are no hole.
	*/
	static const size_type FIRST_EMPTY_POS;

	/*!
		Number of reserved hole values.
	*/
	static const size_type RESERVED_HOLE_COUNT;

	/*!
		At the end of the piecred vector, after all stored elements,
		there should always be at least one sentinel dummy element.
	*/
	static const size_type REQUIRED_SENTINEL_COUNT;

	/*!
		Number of usable positions in the vector.

		The vector must contain at least some sentinel values at
		its end. Moreover, since some hole values have a special
		meaning, the number of usable position should be decreased
		to match the maximum number of holes that is possible to
		store.
	*/
	static const size_type USABLE_POS_COUNT;

public:

	/*!
		Iterator for the pierced array.
	*/
	typedef PiercedIterator<value_type> iterator;

	/*!
		Constant iterator for the pierced array.
	*/
	typedef PiercedIterator<const value_type> const_iterator;

	/*!
		Iterator for the pierced array raw container.
	*/
	typedef typename std::vector<T>::iterator raw_iterator;

	/*!
		Constant iterator for the pierced array raw container.
	*/
	typedef typename std::vector<T>::const_iterator raw_const_iterator;

	/*!
		Constructs an empty pierced vector with no elements.
	*/
	PiercedVector()
	{
		clear();
	}

	/*!
		Constructs a pierced vector with a capacity at least enough
		to contain n elements.

		\param n the minimum capacity requested for the vector
	*/
	PiercedVector(size_type n)
	{
		clear();

		m_v.reserve(n);
	}

	/*!
		Returns a reference to the element with the specified id. If
		there is no element with the specified id, an exception is
		thrown.

		\param id the id of the element
		\result A reference to the element with the specified id.
	*/
	value_type & at(const id_type &id)
	{
		return (*this)[id];
	}

	/*!
		Returns a constant reference to the element with the specified
		id. If there is no element with the specified id, an exception
		is thrown.

		\param id the id of the element
		\result A constant reference to the element with the specified
		        id.
	*/
	const value_type & at(const id_type &id) const
	{
		return (*this)[id];
	}

	/*!
		Returns a reference to the last element of the vector. If
		the vector is empty, an exception is thrown.

		\result A reference to the last element of the vector.
	*/
	value_type & back()
	{
		if (empty()) {
			throw std::out_of_range ("Vector is empty");
		}

		return m_v[m_last_pos];
	}

	/*!
		Returns a constant reference to the last element of the vector.
		If the vector is empty, an exception is thrown.

		\result A constant reference to the last element of the vector.
	*/
	const value_type & back() const
	{
		if (empty()) {
			throw std::out_of_range ("Vector is empty");
		}

		return m_v[m_last_pos];
	}

	/*!
		Returns an iterator pointing to the first element in the
		vector.

		\result An iterator pointing to the first element in the vector.
	*/
	iterator begin()
	{
		if (empty()) {
			return end();
		}

		return iterator(raw_begin() + m_first_pos);
	}

	/*!
		Returns the size of the storage space currently allocated
		for the vector, expressed in terms of elements.

		\result The size of the currently allocated storage capacity
		        in the vector, measured in terms of the number elements
		        it can hold.
	*/
	size_type capacity()
	{
		return m_v.capacity() - REQUIRED_SENTINEL_COUNT;
	}

	/*!
		Returns an conts_iterator pointing to the first element in the
		vector.

		\result A const_iterator pointing to the first element in
		        the vector.
	*/
	const_iterator cbegin() const noexcept
	{
		if (empty()) {
			return cend();
		}

		return const_iterator(m_v.cbegin() + m_first_pos);
	}

	/*!
		Returns an const_iterator referring to the past-the-end element
		in the vector.

		\result A const_iterator referring to the past-the-end element
		        in the vector.
	*/
	const_iterator cend() const noexcept
	{
		return const_iterator(m_v.cbegin() + m_last_pos + 1);
	}

	/*!
		Removes all elements from the vector (which are destroyed),
		leaving the container with a size of 0.
	*/
	void clear()
	{
		// Reset data strucutres
		m_v.clear();

		m_holes.clear();
		std::deque<size_type>().swap(m_holes);

		m_pos.clear();
		std::unordered_map<id_type, size_type, PiercedHasher>().swap(m_pos);

		// Reset first and last counters
		m_first_pos = 0;
		m_last_pos  = 0;

		// Add sentinel
		append_sentinels(1);
	}

	/*!
		Returns a direct pointer to the memory array used internally
		by the vector to store its owned elements.

		\result A pointer to the first element in the array used
		        internally by the vector.

	*/
	value_type * data() noexcept
	{
		return m_v.data();
	}

	/*!
		The container is extended by inserting a new element. This
		new element is constructed in place using args as the
		arguments for its construction.

		\param args the arguments forwarded to construct the new element
	*/
	template <class... Args>
	void emplace(Args&&... args)
	{
		_emplace(FIRST_EMPTY_POS, std::forward<Args>(args)...);
	}

	/*!
		Inserts a new element at the end of the vector, right after
		its current last element. This new element is constructed
		in place using args as the arguments for its construction.

		\param args the arguments forwarded to construct the new element
	*/
	template <class... Args>
	void emplace_back(Args&&... args)
	{
		_emplace(APPEND_TO_BACK, std::forward<Args>(args)...);
	}

	/*!
		Returns whether the vector is empty (i.e. whether its size
		is 0).

		\result true if the container size is 0, false otherwise.
	*/
	bool empty() const
	{
		return m_pos.empty();
	}

	/*!
		Returns an iterator referring to the past-the-end element
		in the vector.

		\result An iterator referring to the past-the-end element
		        in the vector.
	*/
	iterator end()
	{
		return iterator(raw_begin() + m_last_pos + 1);
	}

	/*!
		Returns a reference to the first element of the vector. If
		the vector is empty, an exception is thrown.

		\result A reference to the first element of the vector.
	*/
	value_type & front()
	{
		if (empty()) {
			throw std::out_of_range ("Vector is empty");
		}

		return m_v[m_first_pos];
	}

	/*!
		Returns a constant reference to the first element of the vector.
		If the vector is empty, an exception is thrown.

		\result A constant reference to the first element of the vector.
	*/
	const value_type & front() const
	{
		if (empty()) {
			throw std::out_of_range ("Vector is empty");
		}

		return m_v[m_first_pos];
	}

	/*!
		Removes from the vector the element with the specified id.
		If the id does not exists the function throws an
		exception.

		\param id the id of the element to erase
		\result An iterator pointing to the new location of the
		        element that followed the element erased by the
		        function call. This is the container end if the
		        operation erased the last element in the sequence.
	*/
	iterator erase(id_type id)
	{
		return _erase(get_pos_from_id(id));
	}

	/*!
		Checks if a given id exists in the vector.

		\param id the id to look for
		\result Returns true is the given id exists in the vector,
		        otherwise it returns false.
	*/
	bool exists(id_type id)
	{
		return (m_pos.count(id) != 0);
	}

	/*!
		Gets a vector containing the ids of the elements stored in
		the vector.

		\param ordered if is true the ids will be sorted in ascending
		\              order, otherwise the ids will be in random
		               order.
		\result A vector with the id of the elements in stored in the
		        vector.
	*/
	std::vector<id_type> get_ids(bool ordered = true)
	{
		std::vector<id_type> ids;
		ids.reserve(m_pos.size());

		if (ordered) {
			for(auto item : m_pos) {
				typename std::vector<id_type>::iterator itr = lower_bound(ids.begin(), ids.end(), item.first);
				ids.insert(itr, item.first);
			}
		} else {
			for(auto const &value : *this) {
				ids.push_back(value.get_id());
			}
		}

		return ids;
	}

	/*!
		The container is extended by inserting a new element.

		\param value is the value to be copied (or moved) to the
		            inserted elements.
	*/
	iterator insert(value_type &&value)
	{
		return _insert(FIRST_EMPTY_POS, std::move(value));
	}

	/*!
		Returns the maximum number of elements that the vector can hold.

		This is the maximum potential size the container can reach due
		to known system or library implementation limitations, but the
		container is by no means guaranteed to be able to reach that
		size: it can still fail to allocate storage at any point before
		that size is reached.
	*/
	size_type max_size() const
	{
		return USABLE_POS_COUNT;
	}

	/*!
		Removes the last element in the vector, effectively reducing
		the container size by one.

		This destroys the removed element.
	*/
	void pop_back()
	{
		if (empty()) {
			throw std::out_of_range ("Vector is empty");
		}

		_erase(m_last_pos);
	}

	/*!
		Adds a new element at the end of the vector, after its current
		last element.

		The content of value is copied (or moved) to the new element.

		\param value the value to be copied (or moved) to the new
		             element
	*/
	void push_back(value_type &&value)
	{
		_insert(APPEND_TO_BACK, std::move(value));
	}

	/*!
		Returns a reference to the element at the specified position.

		\param pos the position of the element
		\result A reference to the element in the specified position.
	*/
	value_type & raw_at(const id_type &pos)
	{
		return m_v[pos];
	}

	/*!
		Returns a constant reference to the element at the specified
		position.

		\param pos the position of the element
		\result A constant reference to the element in the specified
		        position.
	*/
	const value_type & raw_at(const id_type &pos) const
	{
		return m_v[pos];
	}

	/*!
		Returns an iterator pointing to the first element in the
		raw container.

		\result An iterator pointing to the first element in the raw
		        container.
	*/
	raw_iterator raw_begin()
	{
		return m_v.begin();
	}

	/*!
		Returns an conts_iterator pointing to the first element in the
		raw container.

		\result A const_iterator pointing to the first element in
		        the raw container.
	*/
	raw_const_iterator raw_cbegin() const noexcept
	{
		return m_v.cbegin();
	}

	/*!
		Returns an const_iterator referring to the past-the-end element
		in raw container.

		\result A const_iterator referring to the past-the-end element
		        in raw container.
	*/
	raw_const_iterator raw_cend() const noexcept
	{
		return m_v.cend();
	}

	/*!
		Returns an iterator referring to the past-the-end element
		in the raw container.

		\result An iterator referring to the past-the-end element
		        in the raw container.
	*/
	raw_iterator raw_end()
	{
		return m_v.end();
	}

	/*!
		Requests that the vector capacity be at least enough to contain
		n elements.

		If n is greater than the current vector capacity, the function
		causes the container to reallocate its storage increasing its
		capacity to n (or greater).

		In all other cases, the function call does not cause a
		reallocation and the vector capacity is not affected.

		\param n the minimum capacity requested for the vector
	*/
	void reserve(size_type n)
	{
		m_v.reserve(n + REQUIRED_SENTINEL_COUNT);
	}

	/*!
		 Resizes the container so that it contains n elements.

		If n is smaller than the current container size, the content
		is reduced to its first n elements, removing those beyond
		(and destroying them).

		If n is greater than the current container size, the content
		is expanded by inserting at the end as many dummy elements as
		needed to reach a size of n. 

		If n is also greater than the current container capacity, an
		automatic reallocation of the allocated storage space takes
		place.

		Notice that this function changes the actual content of the
		container by inserting or erasing elements from it.

		\param n is the new container size, expressed in number of
		         elements.
	*/
	void resize (size_type n)
	{
		// If the size of the vector is already the requested size
		// there is nothing to do.
		if (n == size()) {
			return;
		}

		// If the requested size is greater that the current size
		// we need to add as many sentinel elements as needed
		// to reach a size of n.
		if (n > size()) {
			append_sentinels(n - size());
			return;
		}

		// If the requested size is smaller that the current size
		// we need to perform a real resize.

		// Find the updated position of the last element
		size_type updated_last_pos = n - 1;
		if (!m_holes.empty()) {
			size_type nHoles = 0;
			std::deque<size_type>::iterator it_begin = m_holes.begin();
			std::deque<size_type>::iterator it_end   = m_holes.end();
			while (true) {
				std::deque<size_type>::iterator it_hole = upper_bound(it_begin, it_end, updated_last_pos);
				nHoles = std::distance(it_begin, it_hole);
				if (nHoles == 0) {
					break;
				}

				updated_last_pos += nHoles;

				it_begin = it_hole;
			}
		}

		// Delete all ids of the elements beyond the updated position
		// of the last element
		iterator itr(raw_begin() + updated_last_pos + 1);
		while (itr != end()) {
			m_pos.erase((*itr).get_id());
			itr++;
		}

		// Delete all holes above the updated last position
		holes_delete_after(updated_last_pos);

		// Resize the vector
		m_v.resize(updated_last_pos + 1);
		append_sentinels(1);

		// Update the position of the last element
		m_last_pos = updated_last_pos;
	}

	/*!
		Returns the number of elements in the vector.

		This is the number of actual objects held in the vector,
		which is not necessarily equal to its storage capacity.

		\result The number of elements in the container.
	*/
	size_type size() const
	{
		return m_pos.size();
	}

	/*!
		Requests the container to compact the elements and reduce
		its capacity to fit its size.

		The request is non-binding, and the function can leave
		the vector with a capacity greater than its size.

		This may cause a reallocation, but has no effect on the vector
		size and cannot alter its elements.
	*/
	void squeeze()
	{
		// Compact the vector
		if (!m_holes.empty()) {
			// Move the elements
			size_type offset = 0;
			for (size_type pos = 0; pos <= m_last_pos; pos++) {
				if (m_holes[offset] == pos) {
					offset++;
					continue;
				}

				if (offset == 0) {
					continue;
				}

				id_type id = m_v[pos].get_id();
				size_type updatedPos = pos - offset;

				m_v[updatedPos] = std::move(m_v[pos]);
				m_pos[id] = updatedPos;
			}
			m_v[m_pos.size()] = std::move(m_v[m_last_pos + 1]);

			// Reset first and last counters
			m_first_pos = 0;
			m_last_pos  = m_pos.size() - 1;

			// There are no more holes
			m_holes.clear();
			std::deque<size_type>().swap(m_holes);
		}

		// Resize
		m_v.resize(m_pos.size() + 1);
		m_v.shrink_to_fit();
	}

	/*!
		Exchanges the content of the container by the content of x,
		which is another vector object of the same type. Sizes may
		differ.

		After the call to this member function, the elements in this
		container are those which were in x before the call, and the
		elements of x are those which were in this. All iterators,
		references and pointers remain valid for the swapped objects.

		\param x Another vector container of the same type (i.e.,
		         instantiated with the same template parameters, T and
		         Alloc) whose content is swapped with that of this
		         container.
	*/
	void swap(PiercedVector& x) noexcept
	{
		std::swap(x.m_first_pos, m_first_pos);
		std::swap(x.m_last_pos, m_last_pos);
		std::swap(x.m_v, m_v);
		std::swap(x.m_holes, m_holes);
		std::swap(x.m_pos, m_pos);
	}

	/*!
		Assigns new contents to the vector, replacing its current
		contents.

		\param other vector with the new contents
	*/
	PiercedVector & operator= (PiercedVector&& other)
	{
		if (this != &other) {
			other.swap(*this);
		}

		return *this;
	}


	/*!
		Returns a constant reference to the element with the
		specified id. If there is no element with the specified id,
		an exception is thrown.

		\param id the id of the element
		\result A constant reference to the element with the specified
		        id.
	*/
	const value_type & operator[](const id_type &id) const
	{
		size_type pos = get_pos_from_id(id);

		return m_v[pos];
	}

	/*!
		Returns a reference to the element with the
		specified id. If there is no element with the specified id,
		an exception is thrown.

		\param id the id of the element
		\result A reference to the element with the specified id.
	*/
	value_type & operator[](const id_type &id)
	{
		size_type pos = get_pos_from_id(id);

		return m_v[pos];
	}


private:
	/*!
		Hasher for the id map.

		Since the id are uniques, the hasher can be a function that
		takes the id and cast it to a size_t.

		The hasher is defined as a struct, because a struct can be
		passed as an object into metafunctions (meaning that the type
		deduction for the template paramenters can take place, and
		also meaning that inlining is easier for the compiler). A bare
		function would have to be passed as a function pointer.
		To transform a function template into a function pointer,
		the template would have to be manually instantiated (with a
		perhaps unknown type argument).

	*/
	struct PiercedHasher {
		/*!
			Function call operator that casts the specified
			value to a size_t.

			\tparam U type of the value
			\param value is the value to be casted
			\result Returns the value casted to a size_t.
		*/
		template<typename U>
		constexpr std::size_t operator()(U&& value) const noexcept
		{
			return static_cast<std::size_t>(std::forward<U>(value));
		}
	};

	/*!
		Vector that will hold the elements.
	*/
	std::vector<value_type>m_v;

	/*!
		Container that will hold a list of the holes present in
		the piecrecd vector.
	*/
	std::deque<size_type> m_holes;

	/*!
		Map that links the id of the elements and their position
		inside the internal vector.
	*/
	std::unordered_map<id_type, size_type, PiercedHasher> m_pos;

	/*!
		Position of the first element in the internal vector.
	*/
	size_type m_first_pos;

	/*!
		Position of the last element in the internal vector.
	*/
	size_type m_last_pos;

	/*!
		The container is extended by inserting a new element in
		the specified hole. This new element is constructed in
		place using args as the arguments for its construction.

		\param hole is the hole where the new element should be
		            inserted in. There are two special values
		            for this parameter: FIRST_EMPTY_POS and
		            APPEND_TO_BACK. The first meas that the
		            element should be inserted in the first empty
		            position (which means in the first hole or at the
		            end of the vectpr if there are no holes). The
		            latter means that the element should be inserted
		            at the end of the vector.
		\param args the arguments forwarded to construct the new element
	*/
	template <class... Args>
	size_type _emplace(size_type hole, Args&&... args)
	{
		// Position of the element
		size_type pos = get_pos_to_fill(hole);

		// Insert the element
		if (pos == (m_v.size() - 1)) {
			m_v.emplace(raw_begin() + pos, std::forward<Args>(args)...);
		} else {
			m_v[pos] = T(std::forward<Args>(args)...);
		}

		// Elements must have unique id
		id_type id = m_v[pos].get_id();
		if (exists(id)) {
			throw std::out_of_range ("Duplicate id");
		}

		// The position is now occupied by the element
		fill_pos(pos, id);

		// Return the position of the element
		return pos;
	}

	/*!
		Removes from the vector the element at the specified position.

		\param pos the position of the element to erase
		\result An iterator pointing to the new location of the
		        element that followed the element erased by the
		        function call. This is the container end if the
		        operation erased the last element in the sequence.
	*/
	iterator _erase(size_type pos)
	{
		// Id of the element
		id_type id = m_v[pos].get_id();

		// Delete the element
		m_v[pos] = T();

		// Pierce the position
		pierce_pos(pos, id);

		// Return the iterator to the element following the one erased
		iterator itr;
		if (empty() || pos > m_last_pos) {
			itr = end();
		} else {
			itr = raw_begin() + next_used_pos(pos);
		}

		return itr;
	}

	/*!
		The container is extended by inserting a new element in
		the specified hole.

		\param hole is the hole where the new element should be
		            inserted in. There are two special values
		            for this parameter: FIRST_EMPTY_POS and
		            APPEND_TO_BACK. The first meas that the
		            element should be inserted in the first empty
		            position (which means in the first hole or at the
		            end of the vectpr if there are no holes). The
		            latter means that the element should be inserted
		            at the end of the vector.
		\param value is the value to be copied (or moved) to the
		            inserted elements.
		\result An iterator that points to the the newly inserted
		        element.

	*/
	iterator _insert(size_type hole, value_type &&value)
	{
		// Elements muse have unique id
		id_type id = value.get_id();
		if (m_pos.count(id) != 0) {
			throw std::out_of_range ("Duplicate id");
		}

		// Position of the element
		size_type pos = get_pos_to_fill(hole);

		// Insert the element
		if (pos == (m_v.size() - 1)) {
			m_v.insert(m_v.begin() + pos, std::move(value));
		} else {
			m_v[pos] = std::move(value);
		}

		// The position is now occupied by the element
		fill_pos(pos, m_v[pos].get_id());

		// Return the iterator that points to the element
		iterator itr;
		itr = raw_begin() + pos;

		return itr;
	}

	/*!
		Append to the end of the vector the specified number of
		sentinel elements.

		A sentinel element is a dummy element with the special id
		SENTINEL_ID. All elements after the last non-empty position
		are sentinel elements.

		The sentinel element allows the iterator to correctly
		identify the last non-empty position.

		\param n is the number of sentinel elements to append
	*/
	void append_sentinels(size_t n)
	{
		size_type previousSize = m_v.size();
		size_type updatedSize  = previousSize + n;

		m_v.resize(updatedSize);

		for (size_type pos = previousSize; pos < updatedSize; pos++) {
			m_v[pos].set_id(SENTINEL_ID);
		}
	}


	/*!
		Mark the specified position as filled by the element with
		the given id

		\param pos the position to be marked as filled
		\param id the id of the element that occupies the position
	*/
	void fill_pos(size_type pos, id_type id)
	{
		// Add the id to the map
		m_pos[id] = pos;

		// Update first and last counters
		if (m_last_pos < pos) {
			m_last_pos = pos;
		}

		if (m_first_pos > pos) {
			m_first_pos = pos;
		}

		// If previos element is a hole, its id need to be udated
		if (pos > 0) {
			update_pos_id(pos - 1);
		}
	}

	/*!
		Gets the id of the element at the specified position. The
		function returns only the ids of real elements not dummy
		elements. If there is no element with the specified id,
		an exception is thrown.

		\param id the id to look for
		\result The position in the vector of the element with the
		        specified id.
	*/
	size_type get_pos_from_id(id_type id) const
	{
		return m_pos.at(id);
	}

	/*!
		Given a hole, the function returns an empty position that
		can be used to store an element.

		\param hole the hole that should be used to get the empty
		            position. There are two special values
		            for this parameter: FIRST_EMPTY_POS and
		            APPEND_TO_BACK. The first means that the
		            the function should return the first empty
		            position (which means the position of the first
		            hole or the position following the end of the
		            vector). The latter means that the function should
		            return the position following the end of the
		            vector.
	*/
	size_type get_pos_to_fill(size_type hole)
	{
		// If there are no hole the first avilable position is the
		// end of the vector.
		if (hole == FIRST_EMPTY_POS && m_holes.empty()) {
			hole = APPEND_TO_BACK;
		}

		// Find the position
		size_type pos;
		if (hole == APPEND_TO_BACK) {
			if (empty()) {
				pos = 0;
			} else {
				pos = m_last_pos + 1;
			}

			assert(pos < m_v.size());
		} else {
			assert(!m_v.empty());

			if (hole == FIRST_EMPTY_POS) {
				pos = holes_pop();
			} else {
				pos = holes_pop(hole);
			}

			assert(pos < m_v.size() - 1);
		}

		// Some positions are reserved
		if (pos >= USABLE_POS_COUNT) {
			std::stringstream messageStream;
			messageStream << "Positions above " << (USABLE_POS_COUNT - 1) << " are reserved";
			std::string message = messageStream.str();

			throw std::out_of_range(message);
		}

		return pos;
	}

	/*!
		Add a position to the holes.

		The container that hold the hole list has to be always
		kept in ascending order.

		\param pos the position to be added to the holes
	*/
	void holes_add(size_type pos)
	{
		std::deque<size_type>::iterator itr = lower_bound(m_holes.begin(), m_holes.end(), pos);
		m_holes.insert(itr, pos);
	}

	/*!
		Deletes all holes after a specified position.

		\param pos the position after wich all holes have to be
		           deleted
	*/
	void holes_delete_after(size_type pos)
	{
		if (m_holes.empty()) {
			return;
		} else if (m_holes.back() <= pos) {
			return;
		}

		std::deque<size_type>::iterator itr = upper_bound(m_holes.begin(), m_holes.end(), pos);
		m_holes.erase(itr, m_holes.end());
	}

	/*!
		Gets the position in the vector associated with the
		specified hole.

		\param hole the hole to pop
		\result The position of the requested hole.
	*/
	size_type holes_pop(size_type hole = 0)
	{
		size_type pos;
		if (hole == 0) {
			pos = m_holes.front();
			m_holes.pop_front();
		} else {
			pos = m_holes[hole];
			m_holes.erase(m_holes.begin() + hole);
		}

		return pos;
	}

	/*!
		Returns if the specified position is empty.

		A position is considered empty if the element in that
		position has an id less than 0.

		\param pos the position to check
		\result true is the position is empty, false otherwise.
	*/
	bool is_pos_empty(size_type pos)
	{
		return (m_v[pos].get_id() < 0);
	}

	/*!
		Returns the first non-empty position after the specified
		starting position.

		If the starting position is the last posistion, an
		exception is thrown.

		\param pos starting position
		\result The firt non-empty position after the starting
		        position.
	*/
	size_type next_used_pos(size_type pos)
	{
		if (pos == m_last_pos) {
			throw std::out_of_range ("Already in the last position");
		}

		size_type next_pos = pos + 1;
		id_type next_id    = m_v[next_pos].get_id();
		if (next_id >= 0) {
			return next_pos;
		} else {
			return next_pos - next_id;
		}
	}

	/*!
		Mark the position previously occupied the the specified id,
		as empty.

		\param pos the position to be marked as empty
		\param id the id of the element that was occuping the position
	*/
	void pierce_pos(size_type pos, id_type id)
	{
		// Delete id from map
		m_pos.erase(id);

		// Update first and last counters
		if (empty()) {
			m_last_pos  = 0;
			m_first_pos = 0;
		} else {
			if (m_last_pos == pos) {
				m_last_pos = prev_used_pos(pos);
			}

			if (m_first_pos == pos) {
				m_first_pos = next_used_pos(pos);
			}
		}

		// Update id of the dummy element
		update_pos_id(pos);

		// Hole
		//
		// An empty position is considered a hole, only if it's before
		// the last used position. All holes after the last used
		// position need to be removed.
		if (pos < m_last_pos) {
			holes_add(pos);
		} else {
			holes_delete_after(m_last_pos);
		}
	}

	/*!
		Returns the first non-empty position before the specified
		starting position.

		If the starting position is the first posistion, an
		exception is thrown.

		\param pos starting position
		\result The firt non-empty position before the starting
		        position.
	*/
	size_type prev_used_pos(size_type pos)
	{
		if (pos == m_first_pos) {
			throw std::out_of_range ("Already in the firts position");
		}

		size_type prev_pos = pos - 1;
		id_type prev_id = m_v[prev_pos].get_id();

		if (prev_id >= 0) {
			return prev_pos;
		} else {
			return prev_used_pos(prev_pos);
		}
	}

	/*!
		If the specified position is empty, updates the id of the
		dummy element in that position and, if needed, also of the
		dummy element before the specified position.

		The id of the element in an empty position contains the
		distance, measured in number of elements, between the
		dummy element and the next non-empty element (the distance
		is negative). The id of a dummy element past the last
		non-empty position is set to the special value SENTINEL_ID.

		\param pos the specified position
	*/
	void update_pos_id(size_type pos)
	{
		// Only the ids of the holes need to be updated
		if (!is_pos_empty(pos)) {
			return;
		}

		// Id associated to the position
		id_type id;
		if (empty() || pos > m_last_pos) {
			id = SENTINEL_ID;
		} else if (is_pos_empty(pos + 1)) {
			id = m_v[pos + 1].get_id() - 1;
		} else {
			id = - 1;
		}

		m_v[pos].set_id(id);

		// Update the id of the previous position
		if (pos > 0) {
			update_pos_id(pos - 1);
		}
	}

};

// Definition of static constants of PiercedVector
template<class T>
const typename PiercedVector<T>::id_type
	PiercedVector<T>::SENTINEL_ID = std::numeric_limits<id_type>::min();

template<class T>
const typename PiercedVector<T>::size_type
	PiercedVector<T>::APPEND_TO_BACK = std::numeric_limits<size_type>::max();

template<class T>
const typename PiercedVector<T>::size_type
	PiercedVector<T>::FIRST_EMPTY_POS = std::numeric_limits<size_type>::max() - 1;

template<class T>
const typename PiercedVector<T>::size_type
	PiercedVector<T>::RESERVED_HOLE_COUNT = 2;

template<class T>
const typename PiercedVector<T>::size_type
	PiercedVector<T>::REQUIRED_SENTINEL_COUNT = 1;

template<class T>
const typename PiercedVector<T>::size_type
	PiercedVector<T>::USABLE_POS_COUNT = std::numeric_limits<size_type>::max() - std::max(RESERVED_HOLE_COUNT, REQUIRED_SENTINEL_COUNT);

}

#endif
