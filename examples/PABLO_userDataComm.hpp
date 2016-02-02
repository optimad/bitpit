/*
 * Userdatacomm.hpp
 *
 *  Created on: 18/mar/2014
 *      Author: Marco Cisternino
 */

#ifndef USERDATACOMM_HPP_
#define USERDATACOMM_HPP_

#include "DataCommInterface.hpp"

/*!  \cond  EXAMPLE_CLASSES */
template <class D>
class UserDataComm : public DataCommInterface< UserDataComm<D> > {
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

	UserDataComm(Data & data_, Data & ghostData_);
	~UserDataComm();
};
/*!  \endcond  */

#include "PABLO_userDataComm.tpp"

#endif /* USERDATACOMM_HPP_ */
