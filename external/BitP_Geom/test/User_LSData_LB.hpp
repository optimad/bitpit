/*
 * User_lsdata_LB.hpp
 *
 *  Created on: 20/jul/2015
 *      Author: Edoardo Lombardi
 */

#ifndef USER_LSDATA_LB_HPP_
#define USER_LSDATA_LB_HPP_

#include "DataLBInterface.hpp"

template <class D>
class User_LSData_LB : public DataLBInterface<User_LSData_LB<D> >{
public:

	typedef D Data;

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

	User_LSData_LB(Data& data_);
	~User_LSData_LB();
};

#include "User_LSData_LB.tpp"

#endif /* USER_LSDATA_LB_HPP_ */
