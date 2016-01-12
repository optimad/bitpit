// Communications Class for Ghosts Update

#ifndef CLASSDATACOMMINTERFACE_HPP_
#define CLASSDATACOMMINTERFACE_HPP_

/*!
 *  \ingroup PABLO
 *  @{
 *
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
 *	This class is the base class used to implement the user interface to data communications.
 *
 *	The Curiously Recurrent Template Pattern is exploited to achieve the interface.
 *	By this way the interface is based on static polymorphism with no extra cost at runtime.
 *
 *	The user has to implement his communication classes by deriving them from this class.
 *	The mechanism implies that the derived class derives from a template base class and that the template parameter is the derived class itself, as follow
 *	class Derived : public Base<Derived>{...}
 *
 *	The user has to implement all the methods of the base class in his derived class.
 *	These user's methods will really do the job.
 *
 *	An example of user derived class can be found, [here](https://github.com/optimad/PABLO/blob/master/test/User_Data_Comm.hpp),in test folder for the case of a double datum for each grid element.
 *
 *	Easily speaking, only the user knows his data and through the interface specialization he states the size of element data, how to write/read them in a communication buffer.
 *	Any MPI compatible POD datum can be written and read in the communication buffer.
 */

template <class Impl>
class ClassDataCommInterface {
public:
	/*! Its user specification computes the specific size of data for an element.
	 * \param[in] e Element local index.
	 * \return the size of the data for the e element
	 */
	size_t size(const uint32_t e) const;
	/*! Its user specification computes the same data size for every element in the grid.
	 * \return the size of the data for every element
	 */
	size_t fixedSize() const;


	/*!  Its user specification writes the e element data to be communicated in the buffer.
	 *
	 * The user has not to care about the buffer but a char buffer is available in PABLO, Class_Comm_Buffer. This class has an important method, Class_Comm_Buffer#write.
	 *
	 * This method has to be used to allocate any single element datum in the communication buffer, as follow
	 * ~~~~~~~~~~~~~~~~~~~{.c}
	 * buff.write(userdatum)
	 * ~~~~~~~~~~~~~~~~~~~
	 * where userdatum can be any MPI compatible POD variable associated to the e element.
	 *
	 * In case of a vector of double, called userdata, to store data,
	 * ~~~~~~~~~~~~~~~~~~~{.c}
	 * buff.write(userdata[e])
	 * ~~~~~~~~~~~~~~~~~~~
	 * \param[in] buff Template communication buffer
	 * \param[in] e The element local index
	 */
	template<class Buffer>
	void gather(Buffer & buff,const uint32_t e);

	/*!  Its user specification reads the e element data from the communication buffer and store them in the ghost user data container.
	 *
	 * The user has not to care about the buffer but a char buffer is available in PABLO, Class_Comm_Buffer. This class has an important method, Class_Comm_Buffer#read.
	 *
	 * This method has to be used to read any single element datum from the communication buffer, as follow
	 * ~~~~~~~~~~~~~~~~~~~{.c}
	 * buff.read(userdatum)
	 * ~~~~~~~~~~~~~~~~~~~
	 * where userdatum can be any MPI compatible POD variable associated to the e element.
	 *
	 * In case of a vector of double, called userdata, to store data,
	 * ~~~~~~~~~~~~~~~~~~~{.c}
	 * buff.read(userdata[e])
	 * ~~~~~~~~~~~~~~~~~~~
	 * \param[in] buff Template communication buffer
	 * \param[in] e The element local index
	 */
	template<class Buffer>
	void scatter(Buffer & buff,const uint32_t e);

protected:
	ClassDataCommInterface();

private:
	//BartonHackman trick
	Impl& getImpl();
	const Impl& getImpl() const;
};

/*  @{ */

#include "ClassDataCommInterface.tpp"

#endif /* CLASSDATACOMMINTERFACE_HPP_ */
