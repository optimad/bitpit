/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
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
#endif

using namespace std;
using namespace bitpit;

// =================================================================================== //
/*!
  \example PABLO_example_00006.cpp
  
  Parallel 2D smoothing data using PABLO.
  
  
  The example is the parallel version of PABLO_example_00003.
  
  In order to perform the smoothing procedure in parallel, ghost elements and data communications towards them are needed. Now ghost quadrants and ghost data are not only instantiated, as in example 00003, but actually used to perform smoothing across processor borders.
  
  The user data communication interfaces are based on the Couriously Recurrent Template Pattern. The user has to implement a specification of the interface by writing a derived class.
  In the files PABLO_userDataComm.hpp and PABLO_userDataComm.tpp an example of this specification is given in the case of user data stored in a POD container similar to the STL vector.
  
  The class in PABLO_userDataComm.hpp is an example of user specification of the data communication interface based on the Curiously Recurrent Template Pattern.
  The user has to implement his interface class(es) in order to define how his data have to be written and read in the communication buffer.
  These classes have to be derived from the template base class bitpit::DataCommInterface using as template argument value the derived class. Like this,
  @code
  template <class D>
  class UserDataComm : public bitpit::DataCommInterface< UserDataComm<D> > {
  @endcode

  The choice of the members of the class is completely up to the user and they have to be useful to access both internal and ghost data container. In the example user data datatype is given as template parameter in order to pass any container similar to the STL vector.
  
  In any case, the user must at least implement all the methods reported in this example:
  - <b><code>size_t fixedsize()</code> method:</b> this method is automatically called by the manager and it is intended to define the constant size of the data to be communicated per grid element. If all the element of the grid communicate the same size of data, this method must return a value different from zero and equal to the number of byte to be communicated for every element. 	Otherwise, it must return zero and the different data size for each element must be specified by the size method
  @code
  template<class Data>
  inline size_t UserDataComm<Data>::fixedSize() const {
  return 0;
  };
  @endcode
  \return the size in bytes of the data tobe communicated for every element

  - <b><code>size_t size(const uint32_t e)</code> method:</b> this method is automatically called by the manager and it is intended to define the variable size of the data to be communicated of every grid element. In order to make the manager use this method, the fixedsize method has to return zero. Implementing this method, the user can pass the manager the specific data size to be communicated for the element e.
  @code
  template<class Data>
  inline size_t UserDataComm<Data>::size(const uint32_t e) const {
    return sizeof(double);
  };
  @endcode
  \param[in] index of the internal element to be communicated.
  \return the size in bytes of the data tobe communicated for the element e

  - <b><code>void gather(Buffer& buff, const uint32_t e)</code> method: </b> this method is automatically called by the manager and it is intended to write user data to the char communication buffer. The user has to specify in its implementation the way data can be written in the char buffer. The manager provide the user with a buffer and its write method in order to simply write POD data in the buffer.	The user has to define the way his data can be written in the buffer by decomposing them in POD data and by using the buffer write method to store them in the buffer. In this example we suppose that data is a container of POD data having the random access operator.
  @code
  template<class Data>
  template<class Buffer>
  inline void UserDataComm<Data>::gather(Buffer& buff, const uint32_t e) {
    buff.write(data[e]);
  };
  @endcode
  \param[in] e index of the internal element to be communicated.
  \param[in] buff is the char buffer used to communicate user data. The user has not to take care of the buffer, but its method write and read. These methods are intended to write/read POD data to/fromthe buffer.

  - <b><code>void scatter(Buffer& buff, const uint32_t e)</code> method: </b> this method is automatically called by the manager and it is intended to read user data from the char communication buffer and store them in the ghost data container. The user has to specify in its implementation the way data can be read from the char buffer.	The manager provide the user with a buffer and its read method in order to simply read POD data from the buffer. The user has to define the way his data can be read from the buffer by decomposing them in POD data and by using the buffer read method to take them from the buffer and to store them in the ghost data container. In this example we suppose that ghostData is a container of POD data having the random access operator.
  @code
  template<class Data>
  template<class Buffer>
  inline void UserDataComm<Data>::scatter(Buffer& buff,	const uint32_t e) {
    buff.read(ghostData[e]);
  };
  @endcode
  \param[in] e index of the ghost element to be communicated.
  \param[in] buff is the char buffer used to communicate user data. The user has not to take care of the buffer, but its method write and read. These methods are intended to write/read POD data to/fromthe buffer.

  - <b><code>UserDataComm(Data& data_, Data& ghostdata_)</code> the constructor method: </b> this method has to be called by the user in his application code. The user is free to implement his constructors as he wants, but he must guarantee the access to the internal and ghost data.
  @code 
  template<class D>
  inline UserDataComm<D>::UserDataComm(Data& data_, Data& ghostdata_) : data(data_), ghostdata(ghostdata_){}
  @endcode
  
  In the code of this example application, pay attention to the use of the interface
  @code
  UserDataLB<vector<double> > data_lb(weight,weightGhost);
  pablo7.loadBalance(data_lb, &weight);
  @endcode
  
  <b>To run</b>: ./PABLO_example_00006 \n
  
  <b>To see the result visit</b>: <a href="http://optimad.github.io/PABLO/">PABLO website</a> \n
  
*/
// =================================================================================== //

int main(int argc, char *argv[]) {

#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc, &argv);

	{
#endif
		int iter = 0;
		int dim = 2;

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
		PabloUniform pablo6(2);

		/**<Refine globally four level and write the octree.*/
		for (iter=1; iter<5; iter++){
			pablo6.adaptGlobalRefine();
		}

#if BITPIT_ENABLE_MPI==1
		/**<PARALLEL TEST: Call loadBalance, the octree is now distributed over the processes.*/
		pablo6.loadBalance();
#endif

		/**<Define a center point and a radius.*/
		double xc, yc;
		xc = yc = 0.5;
		double radius = 0.25;

		/**<Define vectors of data.*/
		uint32_t nocts = pablo6.getNumOctants();
		uint32_t nghosts = pablo6.getNumGhosts();
		vector<double> oct_data(nocts, 0.0), ghost_data(nghosts, 0.0);

		/**<Assign a data to the octants with at least one node inside the circle.*/
		for (unsigned int i=0; i<nocts; i++){
			vector<array<double,3> > nodes = pablo6.getNodes(i);
			for (int j=0; j<4; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					oct_data[i] = 1.0;
				}
			}
		}

		/**<Assign a data to the ghost octants (PARALLEL TEST) with at least one node inside the circle.*/
		for (unsigned int i=0; i<nghosts; i++){
			/**<Compute the nodes of the octant (Use pointer for ghost).*/
			Octant *oct = pablo6.getGhostOctant(i);
			vector<array<double,3> > nodes = pablo6.getNodes(oct);
			for (int j=0; j<4; j++){
				double x = nodes[j][0];
				double y = nodes[j][1];
				if ((pow((x-xc),2.0)+pow((y-yc),2.0) <= pow(radius,2.0))){
					ghost_data[i] = 1.0;
				}
			}
		}

		/**<Update the connectivity and write the para_tree.*/
		iter = 0;
		pablo6.updateConnectivity();
		pablo6.writeTest("pablo00006_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data);

		/**<Smoothing iterations on initial data*/
		int start = iter + 1;
		for (iter=start; iter<start+25; iter++){
			vector<double> oct_data_smooth(nocts, 0.0);
			vector<uint32_t> neigh, neigh_t;
			vector<bool> isghost, isghost_t;
			uint8_t iface, nfaces, codim;
			for (unsigned int i=0; i<nocts; i++){
				neigh.clear();
				isghost.clear();

				/**<Find neighbours through faces (codim=1) and edges (codim=2) of the octants*/
				for (codim=1; codim<dim+1; codim++){
					if (codim == 1){
						nfaces = 4;
					}
					else if (codim == 2){
						nfaces = 4;
					}
					for (iface=0; iface<nfaces; iface++){
						pablo6.findNeighbours(i,iface,codim,neigh_t,isghost_t);
						neigh.insert(neigh.end(), neigh_t.begin(), neigh_t.end());
						isghost.insert(isghost.end(), isghost_t.begin(), isghost_t.end());
					}
				}

				/**<Smoothing data with the average over the one ring neighbours of octants*/
				oct_data_smooth[i] = oct_data[i]/(neigh.size()+1);
				for (unsigned int j=0; j<neigh.size(); j++){
					if (isghost[j]){
						oct_data_smooth[i] += ghost_data[neigh[j]]/(neigh.size()+1);
					}
					else{
						oct_data_smooth[i] += oct_data[neigh[j]]/(neigh.size()+1);
					}
				}
			}

			/**<Update the connectivity and write the para_tree.*/
			pablo6.updateConnectivity();
			pablo6.writeTest("pablo00006_iter"+to_string(static_cast<unsigned long long>(iter)), oct_data_smooth);

#if BITPIT_ENABLE_MPI==1
			/**<Communicate the data of the octants and the ghost octants between the processes.*/
			UserDataComm<vector<double> > data_comm(oct_data_smooth, ghost_data);
			pablo6.communicate(data_comm);

#endif
			oct_data = oct_data_smooth;

		}
#if BITPIT_ENABLE_MPI==1
	}

	MPI_Finalize();
#endif
}


