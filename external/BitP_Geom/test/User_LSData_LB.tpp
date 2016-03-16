/*
 * User_LSData_LB.tpp
 *
 *  Created on: 20/jul/2015
 *      Author: Edoardo Lombardi
 */

template<class D>
inline size_t User_LSData_LB<D>::fixedSize() const {
	return 0;
}

template<class D>
inline size_t User_LSData_LB<D>::size(const uint32_t e) const {
	if (data.sdfGrad){
	return sizeof(double) +
			sizeof(double)*data.lg[0].size() +
			sizeof(uint32_t)*data.et[e].size() + sizeof(int) +
			sizeof(int) +
			sizeof(bool);
	}
	else{
		return sizeof(double) +
				sizeof(uint32_t)*data.et[e].size() + sizeof(int) +
				sizeof(int) +
				sizeof(bool);
	}
}

template<class D>
inline void User_LSData_LB<D>::move(const uint32_t from, const uint32_t to) {
	data.ls[to] = data.ls[from];
	if (data.sdfGrad)
		data.lg[to] = data.lg[from];
	std::set<uint32_t>(data.et[from]).swap(data.et[to]);
//	data.et[to] = data.et[from];
	data.supS[to] = data.supS[from];
	data.exact[to] = data.exact[from];

}

template<class D>
template<class Buffer>
inline void User_LSData_LB<D>::gather(Buffer& buff, const uint32_t e) {

	uint32_t val;
	buff.write(data.ls[e]);
	if (data.sdfGrad){
		for(int i = 0; i < data.lg[e].size(); ++i){
			buff.write(data.lg[e][i]);
		}
	}
	int length = data.et[e].size();
	buff.write(length);
	typename std::set<uint32_t>::iterator ib = data.et[e].begin();
	typename std::set<uint32_t>::iterator ie = data.et[e].end();
	for(typename std::set<uint32_t>::iterator i = ib; i != ie; ++i){
		val = *i;
		buff.write(val);
	}
	buff.write(data.supS[e]);
	buff.write(data.exact[e]);

	ib = data.et[e].end();

}

template<class D>
template<class Buffer>
inline void User_LSData_LB<D>::scatter(Buffer& buff, const uint32_t e) {

	uint32_t value;
	buff.read(data.ls[e]);
	if (data.sdfGrad){
		for(int i = 0; i <3; ++i){
			buff.read(data.lg[e][i]);
		}
	}
	int length;
	buff.read(length);
	std::set<uint32_t>().swap(data.et[e]);
	for(int i = 0; i < length; ++i){
		buff.read(value);
		data.et[e].insert(value);
	}
	buff.read(data.supS[e]);
	buff.read(data.exact[e]);
}

template<class D>
inline void User_LSData_LB<D>::assign(uint32_t stride, uint32_t length) {
	Data dataCopy = data;
	typename dvector1D::iterator firstls = dataCopy.ls.begin() + stride;
	typename dvector1D::iterator lastls = firstls + length;
	data.ls.assign(firstls,lastls);
	data.ls.shrink_to_fit();
	firstls = dataCopy.ls.end();
	lastls = dataCopy.ls.end();

	if (data.sdfGrad){
		typename std::vector<std::array<double,3>>::iterator firstlg = dataCopy.lg.begin() + stride;
		typename std::vector<std::array<double,3>>::iterator lastlg = firstlg + length;
		data.lg.assign(firstlg,lastlg);
		data.lg.shrink_to_fit();
		firstlg = dataCopy.lg.end();
		lastlg = dataCopy.lg.end();
	}

	typename std::vector<std::set<uint32_t>>::iterator firstet = dataCopy.et.begin() + stride;
	typename std::vector<std::set<uint32_t>>::iterator lastet = firstet + length;
	data.et.assign(firstet,lastet);
	data.et.shrink_to_fit();
	firstet = dataCopy.et.end();
	lastet = dataCopy.et.end();

	typename std::vector<int>::iterator firstsupS = dataCopy.supS.begin() + stride;
	typename std::vector<int>::iterator lastsupS = firstsupS + length;
	data.supS.assign(firstsupS,lastsupS);
	data.supS.shrink_to_fit();
	firstsupS = dataCopy.supS.end();
	lastsupS = dataCopy.supS.end();

	typename std::deque<bool>::iterator firstex = dataCopy.exact.begin() + stride;
	typename std::deque<bool>::iterator lastex = firstex + length;
	data.exact.assign(firstex,lastex);
	data.exact.shrink_to_fit();
	firstex = dataCopy.exact.end();
	lastex = dataCopy.exact.end();

};

template<class D>
inline void User_LSData_LB<D>::resize(uint32_t newSize) {
	data.ls.resize(newSize);
	if (data.sdfGrad)
		data.lg.resize(newSize, darray3{0.0});
	data.et.resize(newSize);
	data.supS.resize(newSize);
	data.exact.resize(newSize);
}

template<class D>
inline void User_LSData_LB<D>::shrink() {
	data.ls.shrink_to_fit();
	if (data.sdfGrad)
		data.lg.shrink_to_fit();
	data.et.shrink_to_fit();
	data.supS.shrink_to_fit();
	data.exact.shrink_to_fit();
}

template<class D>
inline User_LSData_LB<D>::User_LSData_LB(Data& data_) : data(data_){}

template<class D>
inline User_LSData_LB<D>::~User_LSData_LB() {}



