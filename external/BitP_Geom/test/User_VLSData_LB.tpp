/*
 * User_VLSData_LB.tpp
 *
 *  Created on: 20/jul/2015
 *      Author: Edoardo Lombardi
 */

template<class D>
inline size_t User_VLSData_LB<D>::fixedSize() const {
	return 0;
}

template<class D>
inline size_t User_VLSData_LB<D>::size(const uint32_t e) const {
	size_t sizedata = 0;
	int i;

	for (i=0; i<data.size(); i++){
		if (data[i].sdfGrad){
			sizedata += sizeof(int) + sizeof(double) +
					sizeof(double)*data[i].lg[0].size() +
					sizeof(uint32_t)*data[i].et[e].size() + sizeof(int) +
					sizeof(int) +
					sizeof(bool);
		}
		else{
			sizedata += sizeof(int) + sizeof(double) +
					sizeof(uint32_t)*data[i].et[e].size() + sizeof(int) +
					sizeof(int) +
					sizeof(bool);
		}
	}
	return sizedata;
}


template<class D>
inline void User_VLSData_LB<D>::move(const uint32_t from, const uint32_t to) {
	for (int i=0; i<data.size(); i++){
		data[i].ls[to] = data[i].ls[from];
		if (data[i].sdfGrad)
			data[i].lg[to] = data[i].lg[from];
		data[i].et[to] = data[i].et[from];
		data[i].supS[to] = data[i].supS[from];
		data[i].exact[to] = data[i].exact[from];
	}
}

template<class D>
template<class Buffer>
inline void User_VLSData_LB<D>::gather(Buffer& buff, const uint32_t e) {

	uint32_t val;
	int length0 = data.size();
	buff.write(length0);
	for (int j=0; j<data.size(); j++){
		buff.write(data[j].ls[e]);
		if (data[j].sdfGrad){
			for(int i = 0; i < data[j].lg[e].size(); ++i){
				buff.write(data[j].lg[e][i]);
			}
		}
		int length = data[j].et[e].size();
		buff.write(length);
		typename std::set<uint32_t>::iterator ib = data[j].et[e].begin();
		typename std::set<uint32_t>::iterator ie = data[j].et[e].end();
		for(typename std::set<uint32_t>::iterator i = ib; i != ie; ++i){
			val = *i;
			buff.write(val);
		}
		buff.write(data[j].supS[e]);
		buff.write(data[j].exact[e]);

		ib = data[j].et[e].end();
	}

}

template<class D>
template<class Buffer>
inline void User_VLSData_LB<D>::scatter(Buffer& buff, const uint32_t e) {

	uint32_t value;
	int length0;
	buff.read(length0);
	for(int j = 0; j < length0; ++j){

		buff.read(data[j].ls[e]);
		if (data[j].sdfGrad){
			for(int i = 0; i <3; ++i){
				buff.read(data[j].lg[e][i]);
			}
		}
		int length;
		buff.read(length);
		std::set<uint32_t>().swap(data[j].et[e]);
		for(int i = 0; i < length; ++i){
			buff.read(value);
			data[j].et[e].insert(value);
		}
		buff.read(data[j].supS[e]);
		buff.read(data[j].exact[e]);
	}
}


template<class D>
inline void User_VLSData_LB<D>::assign(uint32_t stride, uint32_t length) {
	for (int i=0; i<data.size(); i++){
		vtype dataCopy = data[i];
		typename dvector1D::iterator firstls = dataCopy.ls.begin() + stride;
		typename dvector1D::iterator lastls = firstls + length;
		data[i].ls.assign(firstls,lastls);
		data[i].ls.shrink_to_fit();
		firstls = dataCopy.ls.end();
		lastls = dataCopy.ls.end();

		if (data[i].sdfGrad){
			typename std::vector<std::array<double,3>>::iterator firstlg = dataCopy.lg.begin() + stride;
			typename std::vector<std::array<double,3>>::iterator lastlg = firstlg + length;
			data[i].lg.assign(firstlg,lastlg);
			data[i].lg.shrink_to_fit();
			firstlg = dataCopy.lg.end();
			lastlg = dataCopy.lg.end();
		}

		typename std::vector<std::set<uint32_t>>::iterator firstet = dataCopy.et.begin() + stride;
		typename std::vector<std::set<uint32_t>>::iterator lastet = firstet + length;
		data[i].et.assign(firstet,lastet);
		data[i].et.shrink_to_fit();
		firstet = dataCopy.et.end();
		lastet = dataCopy.et.end();

		typename std::vector<int>::iterator firstsupS = dataCopy.supS.begin() + stride;
		typename std::vector<int>::iterator lastsupS = firstsupS + length;
		data[i].supS.assign(firstsupS,lastsupS);
		data[i].supS.shrink_to_fit();
		firstsupS = dataCopy.supS.end();
		lastsupS = dataCopy.supS.end();

		typename std::deque<bool>::iterator firstex = dataCopy.exact.begin() + stride;
		typename std::deque<bool>::iterator lastex = firstex + length;
		data[i].exact.assign(firstex,lastex);
		data[i].exact.shrink_to_fit();
		firstex = dataCopy.exact.end();
		lastex = dataCopy.exact.end();
	}
};

template<class D>
inline void User_VLSData_LB<D>::resize(uint32_t newSize) {
	for (int i=0; i<data.size(); i++){
		data[i].ls.resize(newSize);
		if (data[i].sdfGrad)
			data[i].lg.resize(newSize, darray3{0.0});
		data[i].et.resize(newSize);
		data[i].supS.resize(newSize);
		data[i].exact.resize(newSize);
	}
}

template<class D>
inline void User_VLSData_LB<D>::shrink() {
	for (int i=0; i<data.size(); i++){
		data[i].ls.shrink_to_fit();
		if (data[i].sdfGrad)
			data[i].lg.shrink_to_fit();
		data[i].et.shrink_to_fit();
		data[i].supS.shrink_to_fit();
		data[i].exact.shrink_to_fit();
	}
}

template<class D>
inline User_VLSData_LB<D>::User_VLSData_LB(Data& data_) : data(data_){}

template<class D>
inline User_VLSData_LB<D>::~User_VLSData_LB() {}



