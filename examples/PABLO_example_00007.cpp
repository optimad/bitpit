/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include "bitpit_PABLO.hpp"

#if BITPIT_ENABLE_MPI==1
#include "PABLO_userDataComm.hpp"
#include "PABLO_userDataLB.hpp"
#endif

using namespace std;
using namespace bitpit;

// =================================================================================== //
/*!
  \example PABLO_example_00007.cpp
  
  \brief Parallel 2D adaptive mesh refinement (AMR) with data using PABLO
  
  The example is the parallel version of PABLO_example_00004.
  
  Here the main focus is on the load-balance of both grid and data.
  The grid is refined several times together with the data and their inheritance
  follows the same rules like in example 00004.
  Until the last refinement no parallel paradigm is in action: every process owns
  the entire grid.
  
  After this refinement, the load-balance with data is introduced, giving as result
  a parallel distribution of grid and data.
  
  Even a constant vector of weight is used, in order to show that the loadBalance
  can be performed by using a weight function for each cell.
  

 
  The user data communication interfaces are based on the Couriously Recurrent Template Pattern. The user has to implement a specification of the interface by writing a derived class.
  In the files PABLO_userDataLB.hpp and PABLO_userDataLB.tpp an example of this specification is given in the case of user data stored in a POD container similar to the STL vector.
 
  The class in PABLO_userDataLB.hpp is an example of user specification of the load blance data communication interface based on the Curiously Recurrent Template Pattern.
  The user has to implement his interface class(es) in order to define how his data have to be written and read in the communication buffer.
  These classes have to be derived from the template base class bitpit::DataLBInterface using as template argument value the derived class. Like this,
  @code
  template <class D>
  class UserDataLB : public bitpit::DataLBInterface<UserDataLB<D> >
  @endcode
  
  The choice of the members of the class is completely up to the user and they have to be useful to access both internal and ghost data container. In the example user data datatype is given as template parameter in order to pass any container similar to the STL vector.
  
  In any case, the user must at least implement all the methods reported in this example:
  - <b><code>size_t fixedsize()</code> method:</b> this method is automatically called by the manager and it is intended to define the constant size of the data to be communicated per grid element. If all the element of the grid communicate the same size of data, this method must return a value different from zero and equal to the number of byte to be communicated for every element. Otherwise, it must return zero and the different data size for each element must be specified by the size method.
  Example:
  @code
  template<class D>
  inline size_t UserDataLB<D>::size(const uint32_t e) const {
    return 0;
  }
  @endcode
  \return the size in bytes of the data to be communicated for every element
  
  - <b><code>size_t size(const uint32_t e)</code> method:</b> this method is automatically called by the manager and it is intended to define the variable size of the data to be communicated of every grid element. In order to make the manager use this method, the fixedsize method has to return zero. Implementing this method, the user can pass to the manager the specific data size to be communicated for the element e.
  @code
  template<class D>
  inline size_t UserDataLB<D>::size(const uint32_t e) const {
    return sizeof(double);
  }
  @endcode
  \param[in] index of the internal element to be communicated.
  \return the size in bytes of the data tobe communicated for the element e
  
  - <b><code>	void move(const uint32_t from, const uint32_t to)</code> method:</b> this method is automatically called by the manager and it is intended to move data inside the internal data container.
  @code
  template<class D>
  inline void UserDataLB<D>::move(const uint32_t from, const uint32_t to) {
    data[to] = data[from];
  }
  @endcode
  \param[in] from index of the element whose data have to be moved to to.
  \param[in] to index of the element where the from data have to be placed in.
  
  - <b><code>void gather(Buffer& buff, const uint32_t e)</code> method: </b> this method is automatically called by the manager and it is intended to write user data to the char communication buffer. The user has to specify in its implementation the way data can be read from the char buffer. The manager provide the user with a buffer and its read method in order to simply read POD data from the buffer. The user has to define the way his data can be read from the buffer by decomposing them in POD data and by using the buffer read method to take them from the buffer and to store them in the ghost data container. In this example we suppose that data is a container of POD data having the random access operator.
  @code
  template<class D>
  template<class Buffer>
  inline void UserDataLB<D>::gather(Buffer& buff, const uint32_t e) {
    buff.write(data[e]);
  }
  @endcode
  \param[in] e index of the ghost element to be communicated.
  \param[in] buff is the char buffer used to communicate user data. The user has not to take care of the buffer, but its method write and read. These methods are intended to write/read POD data to/fromthe buffer.
  
  - <b><code>void scatter(Buffer& buff, const uint32_t e)</code> method: </b> this method is automatically called by the manager and it is intended to read user data from the char communication buffer and store them in the data container. The user has to specify in its implementation the way data can be read from the char buffer. The manager provide the user with a buffer and its read method in order to simply read POD data from the buffer. The user has to define the way his data can be read from the buffer by decomposing them in POD data and by using the buffer read method to take them from the buffer and to store them in the ghost data container. In this example we suppose that data is a container of POD data having the random access operator.
  @code
  template<class D>
  template<class Buffer>
  inline void UserDataLB<D>::scatter(Buffer& buff, const uint32_t e) {
    buff.read(data[e]);
  }
  @endcode
  \param[in] e index of the element to be communicated.
  \param[in] buff is the char buffer used to communicate user data. The user has not to take care of the buffer, but its method write and read. These methods are intended to write/read POD data to/fromthe buffer.
  
  - <b><code>void assign(uint32_t stride, uint32_t length)</code> method: </b> this method is automatically called by the manager and it intended to be used during the static load balance. At the beginning of any application of the this manager, the entire grid is on every MPI process. At the first load balance the grid is partitioned and distributed to the processes. Data are distributed using this method, by assigning the right length to the the process starting from the right stride. The user has to tell the manager how to start reading the length of data from the stride position and how to put only these data in the same container, deleting the rest. It is a restriction operation of the container to a contiguous part of it. In this examples data is supposed to be a container with the assign operator, as std::vector.
  @code
  template<class D>
  inline void UserDataLB<D>::assign(uint32_t stride, uint32_t length) {
    Data dataCopy = data;
    typename Data::iterator first = dataCopy.begin() + stride;
    typename Data::iterator last = first + length;
    data.assign(first,last);
  #if defined(__INTEL_COMPILER)
  #else
    data.shrink_to_fit();
  #endif
    first = dataCopy.end();
    last = dataCopy.end();
  };
  @endcode
  \param[in] stride the initial position of the data to be 
  
  - <b><code>void resize(uint32_t newSize)</code> method: </b> this method is automatically called by the manager and it intended to be used during the static and the dynamic load balance. During the load balance process the user data container has to change its size to adapt itself to the new partition of the domain. The user has to tell the manager how to change the size of his data containers. In this examples data is supposed to be a container with the resize operator, as std::vector.
  @code
  template<class D>
  inline void UserDataLB<D>::resize(uint32_t newSize) {
    data.resize(newSize);
  }
  @endcode
  \param[in] newSize is the new size of the container
 
  - <b><code>void resizeGhost(uint32_t newSize)</code> method: </b> this method is automatically called by the manager and it intended to be used during the static and the dynamic load balance. During the load balance process the ghost user data container has to change its size to adapt itself to the new partition of the domain. The user has to tell the manager how to change the size of his ghost data containers. In this examples data is supposed to be a container with the resize operator, as std::vector.
  @code
  template<class D>
  inline void UserDataLB<D>::resizeGhost(uint32_t newSize) {
    ghostdata.resize(newSize);
  }
  @endcode
  \param[in] newSize is the new size of the container
 
  - <b><code>void shrink()</code> method: </b> this method is automatically called by the manager and it intended to be used during the static and the dynamic load balance. During the load balance process the user data container changes its size to adapt itself to the new partition of the domain. If the container can change its size and its capacity separately, this method is intended to get the container capacity equal to its size The user has to tell the manager how to change the capacity of his data containers to its size. In this examples data is supposed to be a container with the shrink_to_fit operator, as std::vector.
  @code
  template<class D>
  inline void UserDataLB<D>::shrink() {
  #if defined(__INTEL_COMPILER)
  #else
    data.shrink_to_fit();
  #endif
  }
 @endcode
 
 - <b><code>UserDataLB(Data& data_, Data& ghostdata_)</code> the constructor method: </b> this method has to be called by the user in his application code. The user is free to implement his constructors as he wants, but he must guarantee the access to the internal and ghost data.
 @code 
 template<class D>
 inline UserDataLB<D>::UserDataLB(Data& data_, Data& ghostdata_) : data(data_), ghostdata(ghostdata_){}
 @endcode
 
 In the code of this example application, pay attention to the use of the interface
 @code
 UserDataLB<vector<double> > data_lb(weight,weightGhost);
 pablo7.loadBalance(data_lb, &weight);
 @endcode
 
 <b>To run</b>: ./PABLO_example_00007 \n
 
 <b>To see the result visit</b>: <a href="http://optimad.github.io/PABLO/">PABLO website</a> \n
 
*/
// =================================================================================== //

int main(int argc, char *argv[]) {

#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc, &argv);

	{
#endif
		int iter = 0;

		/**<Instantation and setup of a default (named bitpit) logfile.*/
		int nproc;
		int	rank;
#if BITPIT_ENABLE_MPI==1
		MPI_Comm comm = MPI_COMM_WORLD;
		MPI_Comm_size(comm,&nproc);
		MPI_Comm_rank(comm,&rank);
#else
		nproc = 1;
		rank = 0;
#endif
		log::manager().initialize(log::SEPARATE, false, nproc, rank);
		log::cout() << fileVerbosity(log::NORMAL);
		log::cout() << consoleVerbosity(log::QUIET);

		/**<Instantation of a 2D pablo uniform object.*/
		PabloUniform pablo7(2);

		/**<Set NO 2:1 balance for the octree.*/
		int idx = 0;
		pablo7.setBalance(idx,false);

		/**<Refine globally five level and write the octree.*/
		for (iter=1; iter<6; iter++){
			pablo7.adaptGlobalRefine();
		}

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo7.getNumOctants();
		uint32_t nghosts = pablo7.getNumGhosts();
		vector<double> oct_data(nocts, 0.0), ghost_data(nghosts, 0.0);

		/**<Assign a data (distance from center of a circle) to the octants with at least one node inside the circle.*/
		for (unsigned int i=0; i<nocts; i++){
			/**<Compute the nodes of the octant.*/
			vector<array<double,3> > nodes = pablo7.getNodes(i);
			/**<Compute the center of the octant.*/
			array<double,3> center = pablo7.getCenter(i);
			for (int j=0; j<4; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					oct_data[i] = (pow((center[0]-xc),2.0)+pow((center[1]-yc),2.0));
				}
			}
		}

		/**<Update the connectivity and write the octree.*/
		iter = 0;
		pablo7.updateConnectivity();
		pablo7.writeTest("pablo00007_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

		/**<Adapt two times with data injection on new octants.*/
		int start = 1;
		/**<Weight.*/
		vector<double> weight(nocts, 1.0),weightGhost;
		for (iter=start; iter<start+2; iter++){
			for (unsigned int i=0; i<nocts; i++){
				/**<Compute the nodes of the octant.*/
				vector<array<double,3> > nodes = pablo7.getNodes(i);
				/**<Compute the center of the octant.*/
				array<double,3> center = pablo7.getCenter(i);
				for (int j=0; j<4; j++){
					weight[i] = 2.0;
					double x = nodes[j][0];
					double y = nodes[j][1];
					if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
						if (center[0]<=xc){

							/**<Set to refine to the octants in the left side of the domain inside a circle.*/
							pablo7.setMarker(i,1);
							weight[i] = 1.0;
						}
						else{

							/**<Set to coarse to the octants in the right side of the domain inside a circle.*/
							pablo7.setMarker(i,-1);
							weight[i] = 1.0;
						}
					}
				}
			}

			/**<Adapt the octree and map the data in the new octants.*/
			vector<double> oct_data_new;
			vector<double> weight_new;
			vector<uint32_t> mapper;
			vector<bool> isghost;
			pablo7.adapt(true);
			nocts = pablo7.getNumOctants();
			oct_data_new.resize(nocts, 0.0);
			weight_new.resize(nocts,0.0);

			/**<Assign to the new octant the average of the old children if it is new after a coarsening;
			 * while assign to the new octant the data of the old father if it is new after a refinement.
			 */
			for (uint32_t i=0; i<nocts; i++){
				pablo7.getMapping(i, mapper, isghost);
				if (pablo7.getIsNewC(i)){
					for (int j=0; j<4; j++){
						oct_data_new[i] += oct_data[mapper[j]]/4;
						weight_new[i] += weight[mapper[j]];
					}
				}
				else if (pablo7.getIsNewR(i)){
					oct_data_new[i] += oct_data[mapper[0]];
					weight_new[i] += weight[mapper[0]];
				}
				else{
					oct_data_new[i] += oct_data[mapper[0]];
					weight_new[i] += weight[mapper[0]];
				}
			}

			/**<Update the connectivity and write the octree.*/
			pablo7.updateConnectivity();
			pablo7.writeTest("pablo00007_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data_new);

			oct_data = oct_data_new;
			weight = weight_new;
		}

#if BITPIT_ENABLE_MPI==1
		/**<(Load)Balance the octree over the processes with communicating the data.*/
		UserDataLB<vector<double> > data_lb(weight,weightGhost);
		pablo7.loadBalance(data_lb, &weight);
#endif

		double tot = 0.0;
		for (unsigned int i=0; i<weight.size(); i++){
			tot += weight[i];
		}

		/**<Update the connectivity and write the octree.*/
		pablo7.updateConnectivity();
		pablo7.writeTest("pablo00007_iter"+to_string(static_cast<unsigned long long>(iter)), weight);

#if BITPIT_ENABLE_MPI==1
	}

	MPI_Finalize();
#endif
}



