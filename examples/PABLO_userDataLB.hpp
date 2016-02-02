/*
 * User_data_LB.hpp
 *
 *  Created on: 27/mar/2014
 *      Author: Marco Cisternino
 */

#ifndef USER_DATA_LB_HPP_
#define USER_DATA_LB_HPP_

#include "DataLBInterface.hpp"

/*!  \cond  EXAMPLE_CLASSES */
template <class D>
class UserDataLB : public DataLBInterface<UserDataLB<D> >{
public:

	typedef D Data;

	Data& data;
	Data& ghostdata;

	size_t fixedSize() const;
	size_t size(const uint32_t e) const;
	void move(const uint32_t from, const uint32_t to);

	template<class Buffer>
	void gather(Buffer & buff, const uint32_t e);

	template<class Buffer>
	void scatter(Buffer & buff, const uint32_t e);

	void assign(uint32_t stride, uint32_t length);
	void resize(uint32_t newSize);
	void resizeGhost(uint32_t newSize);
	void shrink();

	UserDataLB(Data& data_, Data& ghostdata_);
	~UserDataLB();
};
/*!  \endcond  */

#include "PABLO_userDataLB.tpp"

#endif /* USER_DATA_LB_HPP_ */
