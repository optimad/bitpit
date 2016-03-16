/*
 * User_vlsdata_LB.hpp
 *
 *  Created on: 20/jul/2015
 *      Author: Edoardo Lombardi
 */

#ifndef USER_VLSDATA_LB_HPP_
#define USER_VLSDATA_LB_HPP_

#include "DataLBInterface.hpp"

template <class D>
class User_VLSData_LB : public DataLBInterface<User_VLSData_LB<D> >{
public:

	typedef D Data;
	typedef typename D::value_type vtype;

	Data& data;

	size_t fixedSize() const;
	size_t size(const uint32_t e) const;
	void move(const uint32_t from, const uint32_t to);

	template<class Buffer>
	void gather(Buffer & buff, const uint32_t e);

	template<class Buffer>
	void scatter(Buffer & buff, const uint32_t e);

	void assign(uint32_t stride, uint32_t length);
	void resize(uint32_t newSize);
	void shrink();

	User_VLSData_LB(Data& data_);
	~User_VLSData_LB();
};

#include "User_VLSData_LB.tpp"

#endif /* USER_VLSDATA_LB_HPP_ */
