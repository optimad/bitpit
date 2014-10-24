/*
 * Userdatacomm.hpp
 *
 *  Created on: 18/mar/2014
 *      Author: Marco Cisternino
 */

#ifndef USERDATACOMM_HPP_
#define USERDATACOMM_HPP_

#include "Class_Data_Comm_Interface.hpp"

template <class D>
class User_Data_Comm : public Class_Data_Comm_Interface< User_Data_Comm<D> > {
public:

	typedef D Data;

	Data & data;
	Data & ghostData;

	size_t fixedSize() const;
	size_t size(const uint32_t e) const;

	template<class Buffer>
	void gather(Buffer & buff, const uint32_t e);

	template<class Buffer>
	void scatter(Buffer & buff, const uint32_t e);

	User_Data_Comm(Data & data_, Data & ghostData_);
	~User_Data_Comm();
};

#include "User_Data_Comm.tpp"

#endif /* USERDATACOMM_HPP_ */
