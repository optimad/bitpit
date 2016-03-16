
template<class D>
inline User_SVData_Comm<D>::User_SVData_Comm(Data& data_, Data & ghostdata_) : data(data_), ghostdata(ghostdata_){
    dataDim = 0;
}

template<class D>
inline User_SVData_Comm<D>::User_SVData_Comm(Data& data_, Data & ghostdata_, int dataDim_) : data(data_), ghostdata(ghostdata_), dataDim(dataDim_){
}

template<class D>
inline User_SVData_Comm<D>::~User_SVData_Comm() {
}

template<class D>
inline size_t User_SVData_Comm<D>::fixedSize() const {
        return 0;
}

template<class D>
inline size_t User_SVData_Comm<D>::size(const uint32_t e) const {
	if(dataDim){
		return sizeof(vtype)*dataDim;
	}
	else{
		return sizeof(vtype)*data[e].size() + sizeof(int);
	}
}

template<class D> template<class Buffer>
inline void User_SVData_Comm<D>::gather(Buffer& buff, const uint32_t e) {
	vtype val;
	if(dataDim){
		typename Data::value_type::iterator ib = data[e].begin();
		typename Data::value_type::iterator ie = data[e].end();
		for(typename Data::value_type::iterator i = ib; i != ie; ++i)
			val = *i;
			buff.write(val);
	}
	else{
		int length = data[e].size();
		buff.write(length);
		typename Data::value_type::iterator ib = data[e].begin();
		typename Data::value_type::iterator ie = data[e].end();
		for(typename Data::value_type::iterator i = ib; i != ie; ++i){
			val = *i;
			buff.write(val);
		}
	}
}

template<class D> template<class Buffer>
inline void User_SVData_Comm<D>::scatter(Buffer& buff, const uint32_t e) {
	vtype gdata;
	if(dataDim){
		for(int i = 0; i < dataDim; ++i){
			buff.read(gdata);
			ghostdata[e].insert(gdata);
		}
	}
	else{
		int length;
		buff.read(length);
		for(int i = 0; i < length; ++i){
			buff.read(gdata);
			ghostdata[e].insert(gdata);
		}
	}

}
