
template<class D>
inline User_VVData_Comm<D>::User_VVData_Comm(Data& data_, Data & ghostdata_) : data(data_), ghostdata(ghostdata_){
    dataDim = 0;
}

template<class D>
inline User_VVData_Comm<D>::User_VVData_Comm(Data& data_, Data & ghostdata_, int dataDim_) : data(data_), ghostdata(ghostdata_), dataDim(dataDim_){
}

template<class D>
inline User_VVData_Comm<D>::~User_VVData_Comm() {
}

template<class D>
inline size_t User_VVData_Comm<D>::fixedSize() const {
        return 0;
}

template<class D>
inline size_t User_VVData_Comm<D>::size(const uint32_t e) const {
	if(dataDim){
		return sizeof(vtype)*dataDim;
	}
	else{
		return sizeof(vtype)*data[e].size() + sizeof(int);
	}
}

template<class D> template<class Buffer>
inline void User_VVData_Comm<D>::gather(Buffer& buff, const uint32_t e) {
	if(dataDim){
		for(int i = 0; i < dataDim; ++i)
			buff.write(data[e][i]);
	}
	else{
		int length = data[e].size();
		buff.write(length);
		for(int i = 0; i < data[e].size(); ++i)
			buff.write(data[e][i]);
	}
}

template<class D> template<class Buffer>
inline void User_VVData_Comm<D>::scatter(Buffer& buff, const uint32_t e) {
	vtype gdata;
	if(dataDim){
		for(int i = 0; i < dataDim; ++i){
			buff.read(ghostdata[e][i]);
		}
	}
	else{
		int length;
		buff.read(length);
//		for(int i = 0; i < length; ++i){
//			buff.read(gdata);
//			ghostdata[e].push_back(gdata);
//		}
		for(int i = 0; i < length; ++i){
			buff.read(gdata);
			ghostdata[e][i] = gdata;
		}
	}

}
