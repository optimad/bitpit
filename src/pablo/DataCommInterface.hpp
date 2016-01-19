// Communications Class for Ghosts Update

#ifndef DATACOMMINTERFACE_HPP_
#define DATACOMMINTERFACE_HPP_

/*!
 *  \ingroup PABLO
 *  @{
 *
 *	\date			09/sep/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
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
 *	The mechanism implies that the derived class derives from a template
 *	base class and that the template parameter is the derived class itself, as follow
 *	class Derived : public Base<Derived>{...}
 *
 *	The user has to implement all the methods of the base class in his derived class.
 *	These user's methods will really do the job.
 *
 *	Easily speaking, only the user knows his data and through the interface
 *	specialization he states the size of element data, how to write/read them
 *	 in a communication buffer.
 *	Any MPI compatible POD datum can be written and read in the communication buffer.
 */

template <class Impl>
class DataCommInterface {
public:
	size_t size(const uint32_t e) const;
	size_t fixedSize() const;

	template<class Buffer>
	void gather(Buffer & buff,const uint32_t e);

	template<class Buffer>
	void scatter(Buffer & buff,const uint32_t e);

protected:
	DataCommInterface();

private:
	//BartonHackman trick
	Impl& getImpl();
	const Impl& getImpl() const;
};

/*  @{ */

#include "DataCommInterface.tpp"

#endif /* DATACOMMINTERFACE_HPP_ */
