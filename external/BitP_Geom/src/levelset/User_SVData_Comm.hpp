#ifndef USERSVDATACOMM_HPP_
#define USERSVDATACOMM_HPP_

#include "Class_Data_Comm_Interface.hpp"

template <class D>
class User_SVData_Comm : public Class_Data_Comm_Interface< User_SVData_Comm<D> > {
public:

        typedef D Data;
        typedef typename D::value_type::value_type vtype;

        Data & data;
        Data & ghostdata;
        int dataDim;

        size_t fixedSize() const;
        size_t size(const uint32_t e) const;

        template<class Buffer>
        void gather(Buffer & buff, const uint32_t e);

        template<class Buffer>
        void scatter(Buffer & buff, const uint32_t e);

        User_SVData_Comm(Data & data_, Data & ghostdata_);
        User_SVData_Comm(Data& data_, Data & ghostdata_, int dataDim_);
        ~User_SVData_Comm();
};

#include "User_SVData_Comm.tpp"

#endif /* USERSVDATACOMM_HPP_ */

