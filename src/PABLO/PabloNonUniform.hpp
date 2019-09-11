/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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

#ifndef __BITPIT_PABLO_UNIFORM_HPP__
#define __BITPIT_PABLO_UNIFORM_HPP__

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "ParaTree.hpp"

namespace bitpit {

    // =================================================================================== //
    // TYPEDEFS			               						   //
    // =================================================================================== //
    typedef std::vector<bool>				bvector;
    typedef std::bitset<72>            			octantID;
    typedef std::vector<Octant*>			ptroctvector;
    typedef ptroctvector::iterator			octantIterator;
    typedef std::vector<darray3>			darray3vector;


    // =================================================================================== //
    // CLASS DEFINITION                                                                    //
    // =================================================================================== //
    /*!
     *  Original Version: PabloUniform
     *	\ingroup		PABLO
     *	\date			25/jan/2016
     *	\authors		Edoardo Lombardi
     *	\authors		Marco Cisternino
     *
     *	\brief PABLO Uniform is an example of user class derived from ParaTree to map
     *	ParaTree in a uniform (square/cubic) domain.
     *
     *  Modified by Luca Verzeroli into Class PabloNonUniform 
     *  \brief PABLO NonUniform is an example of user class derived from ParaTree to map
     *	ParaTree in rectangular/cuboidal domain. 
     *	PabloNonUniform takes as input in constructor
     *  1- The coordinates of the origin (X,Y,Z)
     *  2- The dimension dim=2 and dim=3, for 2D and 3D respectively.
     *  3- The length of the sides LX (in x direction), LY (in Y direction), LZ (in Z direction)
     *  4- The log file name
     *  The Uniform case is recovered by prescribing only LX or using LY=LX and LZ=LX. It assumed that LX=LY=LZ
     */

    class PabloNonUniform : public ParaTree
    {
        // =================================================================================== //
        // MEMBERS									       //
        // =================================================================================== //
    private:
        darray3 	        m_origin;		/**<Coordinate X,Y,Z of the origin of the octree in the physical domain*/
        darray3 	        m_L;		        /**Length of the octree in the physical domain in direction X,Y,Z*/
        std::vector<double>	m_area;		        /**<Area, in 3D, (length in 2D) of octree faces (edge) in the physical domain*/
        double 		        m_volume;		/**Volume, in 3D, (area in 2D) of octree in the physical domain*/

        // =================================================================================== //
        // CONSTRUCTORS AND OPERATORS                                                          //
        // =================================================================================== //

        // =================================================================================== //
        // METHODS                                                                             //
        // =================================================================================== //
        void	__reset();
    public:

// Default constructor of PabloNonUniform. It sets the Origin in (0,0,0) and sides of length 1
#if BITPIT_ENABLE_MPI==1
        PabloNonUniform(std::string logfile = DEFAULT_LOG_FILE, MPI_Comm comm = MPI_COMM_WORLD);
        PabloNonUniform(uint8_t dim, std::string logfile = DEFAULT_LOG_FILE, MPI_Comm comm = MPI_COMM_WORLD);
#else
        PabloNonUniform(std::string logfile = DEFAULT_LOG_FILE);
        PabloNonUniform(uint8_t dim, std::string logfile = DEFAULT_LOG_FILE);
#endif

// Constructors equivalent to PABLO Uniform
#if BITPIT_ENABLE_MPI==1
        PabloNonUniform(double X, double Y, double Z, double L, uint8_t dim = 2, std::string logfile = DEFAULT_LOG_FILE, MPI_Comm comm = MPI_COMM_WORLD); // equivalent to PabloUniform constructor
#else
        PabloNonUniform(double X, double Y, double Z, double L, uint8_t dim = 2, std::string logfile = DEFAULT_LOG_FILE); // equivalent to PabloUniform constructor
#endif

// Constructor PABLONonUniform
#if BITPIT_ENABLE_MPI==1
        PabloNonUniform(double X, double Y, double Z, double LX, double LY, double LZ, uint8_t dim = 2, std::string logfile = DEFAULT_LOG_FILE, MPI_Comm comm = MPI_COMM_WORLD);
#else
        PabloNonUniform(double X, double Y, double Z, double LX, double LY, double LZ, uint8_t dim = 2, std::string logfile = DEFAULT_LOG_FILE);
#endif

        // =================================================================================== //
        // METHODS                                                                             //
        // =================================================================================== //
        void	reset() override;
        int	getDumpVersion() const override;
        void	dump(std::ostream &stream, bool full = true) override;
        void	restore(std::istream &stream) override;

        // =================================================================================== //
        // BASIC GET/SET METHODS	                                                       //
        // =================================================================================== //
        darray3		getOrigin() const;
        double		getX0() const;
        double		getY0() const;
        double		getZ0() const;
        darray3		getL() const;
        void		setL(darray3 L);
        void		setOrigin(darray3 origin);
        double		levelToSize( uint8_t& level, uint8_t dir = 0);

        // =================================================================================== //
        // INDEX BASED METHODS						                       //
        // =================================================================================== //
        darray3 	getCoordinates(uint32_t idx) const;
        double 		getX(uint32_t idx) const;
        double 		getY(uint32_t idx) const;
        double 		getZ(uint32_t idx) const;
        double 		getSize(uint32_t idx, uint8_t dir = 0) const;
        double 		getArea(uint32_t idx, uint8_t iface = 0) const;
        double 		getVolume(uint32_t idx) const;
        void 		getCenter(uint32_t idx, darray3& center) const;
        darray3 	getCenter(uint32_t idx) const;
        darray3 	getFaceCenter(uint32_t idx, uint8_t iface) const;
        void 		getFaceCenter(uint32_t idx, uint8_t iface, darray3& center) const;
        darray3 	getNode(uint32_t idx, uint8_t inode) const;
        void 		getNode(uint32_t idx, uint8_t inode, darray3& node) const;
        void 		getNodes(uint32_t idx, darr3vector & nodes) const;
        darr3vector     getNodes(uint32_t idx) const;
        void 		getNormal(uint32_t idx, uint8_t iface, darray3 & normal) const;
        darray3 	getNormal(uint32_t idx, uint8_t iface) const;

        // =================================================================================== //
        // POINTER BASED METHODS							       //
        // =================================================================================== //
        darray3 	getCoordinates(const Octant* oct) const;
        double 		getX(const Octant* oct) const;
        double 		getY(const Octant* oct) const;
        double 		getZ(const Octant* oct) const;
        double 		getSize(const Octant* oct, uint8_t dir = 0) const;
        double 		getArea(const Octant* oct, uint8_t dir = 0) const;
        double 		getVolume(const Octant* oct) const;
        void 		getCenter(const Octant* oct, darray3& center) const;
        darray3 	getCenter(const Octant* oct) const;
        darray3 	getFaceCenter(const Octant* oct, uint8_t iface) const;
        void 		getFaceCenter(const Octant* oct, uint8_t iface, darray3& center) const;
        darray3 	getNode(const Octant* oct, uint8_t inode) const;
        void 		getNode(const Octant* oct, uint8_t inode, darray3& node) const;
        void 		getNodes(const Octant* oct, darr3vector & nodes) const;
        darr3vector     getNodes(const Octant* oct) const;
        void 		getNormal(const Octant* oct, uint8_t iface, darray3 & normal) const;
        darray3 	getNormal(const Octant* oct, uint8_t iface) const;

        // =================================================================================== //
        // LOCAL TREE GET/SET METHODS						               //
        // =================================================================================== //
        double	 	getLocalMaxSize() const;
        double	 	getLocalMinSize() const;
        void 		getBoundingBox(darray3 & P0, darray3 & P1) const;

        // =================================================================================== //
        // INTERSECTION GET/SET METHODS							       //
        // =================================================================================== //
        double 		getSize(const Intersection* inter, uint8_t dir = 0) const;
        double 		getArea(const Intersection* inter, uint8_t iface = 0) const;
        darray3 	getCenter(const Intersection* inter) const;
        darr3vector     getNodes(const Intersection* inter) const;
        darray3 	getNormal(const Intersection* inter) const;

        // =================================================================================== //
        // OTHER OCTANT BASED METHODS							       //
        // =================================================================================== //
        Octant*        getPointOwner(darray3 point);
        uint32_t       getPointOwnerIdx(darray3 point) const;
        Octant*        getPointOwner(darray3 point, bool & isghost);
        uint32_t       getPointOwnerIdx(darray3 point, bool & isghost) const;
        int            getPointOwnerRank(darray3 point);

        // =================================================================================== //
        // OTHER PARATREE BASED METHODS							       //
        // =================================================================================== //
        darray3 	getNodeCoordinates(uint32_t inode) const;

        // =================================================================================== //
        // TESTING OUTPUT METHODS							       //
        // =================================================================================== //
        void        write(std::string filename);
        void        writeTest(std::string filename, dvector data);
    };

}

#endif /* __BITPIT_PABLO_UNIFORM_HPP__ */
