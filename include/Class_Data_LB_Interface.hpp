// Communications Class for Load Balancing

#ifndef CLASS_DATA_LB_INTERFACE_HPP_
#define CLASS_DATA_LB_INTERFACE_HPP_

#include <stdint.h>

template <class Impl>
class Class_Data_LB_Interface {
public:
	size_t size(const uint32_t e) const;
	size_t fixedSize() const;
	void move(const uint32_t from, const uint32_t to);

	template<class Buffer>
	void gather(Buffer & buff,const uint32_t e);

	template<class Buffer>
	void scatter(Buffer & buff,const uint32_t e);

protected:
	Class_Data_LB_Interface();

private:
	//BartonHackman trick
	Impl& getImpl();
	const Impl& getImpl() const;
};

#include "Class_Data_LB_Interface.tpp"

#endif /* CLASS_DATA_LB_INTERFACE_HPP_ */
