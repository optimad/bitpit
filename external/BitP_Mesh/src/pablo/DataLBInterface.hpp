// Communications Class for Load Balancing

#ifndef DATALBINTERFACE_HPP_
#define DATALBINTERFACE_HPP_

/*!
 *  \ingroup        PABLO
 *  @{
 *
 *	\date			09/sep/2015
 *	\authors		Edoardo Lombardi
 *	\authors		Marco Cisternino
 *	\copyright		Copyright 2014 Optimad engineering srl. All rights reserved.
 *	\par			License
 *	This version of PABLO is released under the LGPL License.
 *
 *	\brief Base class for data communications
 *
 *	This class is the base class used to implement the user interface
 *	to data communications for load balance.
 *
 *	The Curiously Recurrent Template Pattern is exploited to achieve the interface.
 *	By this way the interface is based on static polymorphism with no extra cost at runtime.
 *
 *	The user has to implement his communication classes for load balance
 *	by deriving them from this class.
 *	The mechanism implies that the derived class derives from a template base class
 *	and that the template parameter is the derived class itself, as follow
 *	~~~~~~~~~~~~~~~~~~~{.c}
 *	class Derived : public Base<Derived>{...}
 *	~~~~~~~~~~~~~~~~~~~
 *	The user has to implement all the methods of the base class in his derived class.
 *	These user's methods will really do the job.
 *
 *	An example of user derived class can be found,[here](https://github.com/optimad/PABLO/blob/master/test/User_Data_LB.hpp),in test folder for the case of a double datum for each grid element.
 *
 *	Easily speaking, only the user knows his data and through the interface specialization
 *	he states the size of element data, how to write/read and move them
 *	in a communication buffer.
 */
template <class Impl>
class DataLBInterface {
public:
	size_t size(const uint32_t e) const;
	size_t fixedSize() const;
	void move(const uint32_t from, const uint32_t to);

	template<class Buffer>
	void gather(Buffer & buff,const uint32_t e);

	template<class Buffer>
	void scatter(Buffer & buff,const uint32_t e);

	void assign(uint32_t stride, uint32_t length);
	void resize(uint32_t newSize);
	void resizeGhost(uint32_t newSize);
	void shrink();

protected:
	DataLBInterface();

private:
	//BartonHackman trick
	Impl& getImpl();
	const Impl& getImpl() const;

};

/*  @{ */
#include "DataLBInterface.tpp"

#endif /* DATA_LB_INTERFACE_HPP_ */
