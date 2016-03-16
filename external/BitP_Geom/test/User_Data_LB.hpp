/*
 * User_data_LB.hpp
 *
 *  Created on: 27/mar/2014
 *      Author: Marco Cisternino
 */

#ifndef USER_DATA_LB_HPP_
#define USER_DATA_LB_HPP_

#include "DataLBInterface.hpp"

template <class D>
class User_Data_LB : public DataLBInterface<User_Data_LB<D> >{
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

	User_Data_LB(Data& data_);
	~User_Data_LB();
};

#include "User_Data_LB.tpp"

#endif /* USER_DATA_LB_HPP_ */
