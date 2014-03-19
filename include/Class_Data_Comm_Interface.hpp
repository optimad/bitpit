/*
 * Class_Data_Comm_Interface.hpp
 *
 *  Created on: 18/mar/2014
 *      Author: Marco Cisternino
 */

#ifndef CLASSDATACOMMINTERFACE_HPP_
#define CLASSDATACOMMINTERFACE_HPP_

#include <stdint.h>

template <class Impl>
class Class_Data_Comm_Interface {
public:
	size_t size(const uint32_t & e) const;
	size_t fixedSize() const;

	template<class Buffer>
	void gather(Buffer & buff,const uint32_t & e);

	template<class Buffer>
	void scatter(Buffer & buff,const uint32_t & e);

protected:
	Class_Data_Comm_Interface();

private:
	//BartonHackman trick
	Impl& getImpl();
	const Impl& getImpl() const;
};

#include "Class_Data_Comm_Interface.tpp"

#endif /* CLASSDATACOMMINTERFACE_HPP_ */
