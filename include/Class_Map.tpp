// =================================================================================== //
// CLASS SPECIALIZATION                                                                //
// =================================================================================== //

/*!Default constructor. Origin of octree in physical domain in (0,0,0)
 * and side length 1.
 */
template <int dim>
Class_Map<dim>::Class_Map(){
	X0 = Y0 = Z0 = 0.0;
	L = 1.0;
};

/*!Customized constructor with origin of octree in physical
 * domain side length provided by the user.
 * \param[in] X Coordinate X of the origin.
 * \param[in] Y Coordinate Y of the origin.
 * \param[in] Z Coordinate Z of the origin.
 * \param[in] LL Side length of domain.
 */
template <int dim>
Class_Map<dim>::Class_Map(double & X, double & Y, double & Z, double & LL){
	X0 = X;
	Y0 = Y;
	Z0 = Z;
	L = LL;
};

// ------------------------------------------------------------------------------- //
// METHODS ----------------------------------------------------------------------- //

template <int dim>
double Class_Map<dim>::mapX(uint32_t const & X){
	return (X0 + L/double(globals.max_length) * double(X));
};

template <int dim>
double Class_Map<dim>::mapY(uint32_t const & Y){
	return (Y0 + L/double(globals.max_length) * double(Y));
};

template <int dim>
double Class_Map<dim>::mapZ(uint32_t const & Z){
	return (Z0 + L/double(globals.max_length) * double(Z));
};

template <int dim>
uint32_t Class_Map<dim>::mapX(double const & X){
	return (double(globals.max_length)/L * (double(X) - X0));
};

template <int dim>
uint32_t Class_Map<dim>::mapY(double const & Y){
	return (double(globals.max_length)/L * (double(Y) - Y0));
};

template <int dim>
uint32_t Class_Map<dim>::mapZ(double const & Z){
	return (double(globals.max_length)/L * (double(Z) - Z0));
};

template <int dim>
double Class_Map<dim>::mapSize(uint32_t const & size){
	return ((L/double(globals.max_length))*double(size));
};

template <int dim>
double Class_Map<dim>::mapArea(uint64_t const & Area){
	return ((pow(L,2.0)/pow(double(globals.max_length),2.0))*double(Area));
};

template <int dim>
double Class_Map<dim>::mapVolume(uint64_t const & Volume){
	return ((pow(L,3.0)/pow(double(globals.max_length),3.0))*double(Volume));
};

template <int dim>
void Class_Map<dim>::mapCenter(double* & center,
		vector<double> & mapcenter){
	vector<double> orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	orig.push_back(Z0);
	orig.shrink_to_fit();
	mapcenter.resize(3);
	for (int i=0; i<dim; i++){
		mapcenter[i] = orig[i] + L/double(globals.max_length) * center[i];
	}
	mapcenter.shrink_to_fit();
};

template <int dim>
void Class_Map<dim>::mapCenter(vector<double> & center,
		vector<double> & mapcenter){
	vector<double> orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	orig.push_back(Z0);
	orig.shrink_to_fit();
	mapcenter.resize(3);
	for (int i=0; i<dim; i++){
		mapcenter[i] = orig[i] + L/double(globals.max_length) * center[i];
	}
	mapcenter.shrink_to_fit();
};

template <int dim>
void Class_Map<dim>::mapNodes(uint32_t (*nodes)[3],
		vector<vector<double> > & mapnodes){
	vector<double> orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	orig.push_back(Z0);
	orig.shrink_to_fit();
	mapnodes.resize(globals.nnodes);
	for (int i=0; i<globals.nnodes; i++){
		mapnodes[i].resize(3);
		for (int j=0; j<3; j++){
			mapnodes[i][j] = orig[j] + L/double(globals.max_length) * double(nodes[i][j]);
		}
		mapnodes[i].shrink_to_fit();
	}
	mapnodes.shrink_to_fit();

};

template <int dim>
void Class_Map<dim>::mapNodes(vector<vector<uint32_t> > nodes,
		vector<vector<double> > & mapnodes){
	vector<double> orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	orig.push_back(Z0);
	orig.shrink_to_fit();
	mapnodes.resize(globals.nnodes);
	for (int i=0; i<globals.nnodes; i++){
		mapnodes[i].resize(3);
		for (int j=0; j<3; j++){
			mapnodes[i][j] = orig[j] + L/double(globals.max_length) * double(nodes[i][j]);
		}
		mapnodes[i].shrink_to_fit();
	}
	mapnodes.shrink_to_fit();

};

template <int dim>
void Class_Map<dim>::mapNodesIntersection(uint32_t (*nodes)[3],
		vector<vector<double> > & mapnodes){
	vector<double> orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	orig.push_back(Z0);
	orig.shrink_to_fit();
	mapnodes.resize(globals.nnodesperface);
	for (int i=0; i<globals.nnodesperface; i++){
		mapnodes[i].resize(3);
		for (int j=0; j<3; j++){
			mapnodes[i][j] = orig[j] + L/double(globals.max_length) * double(nodes[i][j]);
		}
		mapnodes[i].shrink_to_fit();
	}
	mapnodes.shrink_to_fit();

};

template <int dim>
void Class_Map<dim>::mapNodesIntersection(vector<vector<uint32_t> > nodes,
		vector<vector<double> > & mapnodes){
	vector<double> orig;
	orig.push_back(X0);
	orig.push_back(Y0);
	orig.push_back(Z0);
	orig.shrink_to_fit();
	mapnodes.resize(globals.nnodesperface);
	for (int i=0; i<globals.nnodesperface; i++){
		mapnodes[i].resize(3);
		for (int j=0; j<3; j++){
			mapnodes[i][j] = orig[j] + L/double(globals.max_length) * double(nodes[i][j]);
		}
		mapnodes[i].shrink_to_fit();
	}
	mapnodes.shrink_to_fit();

};

template <int dim>
void Class_Map<dim>::mapNormals(vector<int8_t> normal,
		vector<double> & mapnormal){
	mapnormal = vector<double>(normal.begin(), normal.end());
	mapnormal.shrink_to_fit();
};
