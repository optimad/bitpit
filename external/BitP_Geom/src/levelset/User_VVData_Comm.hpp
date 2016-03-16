#ifndef USERVVDATACOMM_HPP_
#define USERVVDATACOMM_HPP_

#include "Class_Data_Comm_Interface.hpp"

template <class D>
class User_VVData_Comm : public Class_Data_Comm_Interface< User_VVData_Comm<D> > {
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

        User_VVData_Comm(Data & data_, Data & ghostdata_);
        User_VVData_Comm(Data& data_, Data & ghostdata_, int dataDim_);
        ~User_VVData_Comm();
};

#include "User_VVData_Comm.tpp"

#endif /* USERVVDATACOMM_HPP_ */

