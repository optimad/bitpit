// Communications Class for Load Balancing

#ifndef CLASS_DATA_LB_INTERFACE_HPP_
#define CLASS_DATA_LB_INTERFACE_HPP_

#include <stdint.h>

/*!
 *	\date			09/sep/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\version		0.1
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Base class for data communications
 *
 *	This class is the base class used to implement the user interface to data communications for load balance.
 *	The Curiously Recurrent Template Pattern is exploited to achieve the interface.
 *	By this way the interface is based on static polymorphism with no extra cost at runtime.
 *	The user has to implement his communication classes for load balance by deriving them from this class.
 *	The mechanism implies that the derived class derives from a template base class and that the template parameter is the derived class itself, as follow
 *	class Derived : public Base<Derived>{...}
 *	The user has to implement all the methods of the base class in his derived class.
 *	These user's methods will really do the job.
 *	An example of user derived class can be found in test folder for the case of a double datum for each grid element.
 *	Easily speaking, only the user knows his data and through the interface specialization he states the size of element data, how to write/read and move them in a communication buffer.
 */
template <class Impl>
class Class_Data_LB_Interface {
public:
	/*! Its user specification computes the specific size of data for an element.
	 * \param[in] element local index.
	 * \return the size of the data for the e element
	 */
	size_t size(const uint32_t e) const;
	/*! Its user specification computes the same data size for every element in the grid.
	 * \return the size of the data for every element
	 */
	size_t fixedSize() const;
	/*! Its user specification moves the data from the "from" element to the "to" element
	 * \param[in] element local indices from and to.
	 */
	void move(const uint32_t from, const uint32_t to);

	/*!  Its user specification writes the e element data to be communicated during the load balance in the buffer.
	 * The user has not to care about the buffer but a char buffer is available in PABLO, Class_Comm_Buffer. This class has an important method, write.
	 * This method has to be used to allocate any single element datum in the communication buffer, as follow
	 * buff.write(userdatum)
	 * where userdatum can be any MPI compatible POD variable associated to the e element.
	 * In case of a vector of double, called userdata, to store data, buff.write(userdata[e])
	 * \param[in] buff, template communication buffer
	 * \param[in] e, the element local index
	 */
	template<class Buffer>
	void gather(Buffer & buff,const uint32_t e);

	/*!  Its user specification reads the e element data from the communication buffer and store them in the user data container.
	 * The user has not to care about the buffer but a char buffer is available in PABLO, Class_Comm_Buffer. This class has an important method, read.
	 * This method has to be used to read any single element datum from the communication buffer, as follow
	 * buff.read(userdatum)
	 * where userdatum can be any MPI compatible POD variable associated to the e element.
	 * In case of a vector of double, called userdata, to store data, buff.read(userdata[e])
	 * \param[in] buff, template communication buffer
	 * \param[in] e, the element local index
	 */
	template<class Buffer>
	void scatter(Buffer & buff,const uint32_t e);


	/*!  Its user specification extracts contiguous element from a data container.
	 * 	This method is used during the very first load balance where the grid from being serial becomes parallel.
	 * 	From the initial serial container, which is the same on every process, this method takes the data relatives to the process sub-domain, storing them in the same resized parallel local container.
	 * \param[in] stride is the location where the process has to start reading from
	 * \param[in] length is the length of the data container the process has to read
	 */
	void assign(uint32_t stride, uint32_t length);

	/*!  Its user specification resizes the element data container to the newSize value.
	 * \param[in] newSize is the new size of the user data container
	 */
	void resize(uint32_t newSize);
	/*!  Its user specification reduces the capacity of the container to its size for those containers which has different values for capacity and size.
	 * 	If the user container has no capacity concept, nothing has to be done.
	 */
	void shrink();

protected:
	Class_Data_LB_Interface();

private:
	//BartonHackman trick
	Impl& getImpl();
	const Impl& getImpl() const;
};

#include "Class_Data_LB_Interface.tpp"

#endif /* CLASS_DATA_LB_INTERFACE_HPP_ */
