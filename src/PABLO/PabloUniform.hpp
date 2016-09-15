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

#ifndef __BITPIT_PABLO_UNIFORM_HPP__
#define __BITPIT_PABLO_UNIFORM_HPP__

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "ParaTree.hpp"

namespace bitpit {

    // =================================================================================== //
    // TYPEDEFS																			   //
    // =================================================================================== //
    typedef std::vector<bool>				bvector;
    typedef std::bitset<72>					octantID;
    typedef std::vector<Octant*>			ptroctvector;
    typedef ptroctvector::iterator			octantIterator;
    typedef std::vector<darray3>			darray3vector;


    // =================================================================================== //
    // CLASS DEFINITION                                                                    //
    // =================================================================================== //
    /*!
     *  \ingroup        PABLO
     *  @{
     *	\date			25/jan/2016
     *	\authors		Edoardo Lombardi
     *	\authors		Marco Cisternino
     *
     *	\brief PABLO Uniform is an example of user class derived from ParaTree to map
     *	ParaTree in a uniform (square/cubic) domain.
     *	Pablo Uniform takes as input in constructor the coordinates of the origin (X,Y,Z) and the length of the side L.
     *
     *	Class PabloUniform has a dimensional parameter int dim and it accepts
     *	only two values: dim=2 and dim=3, for 2D and 3D respectively.
     */
    class PabloUniform : public ParaTree
    {
        // =================================================================================== //
        // MEMBERS																			   //
        // =================================================================================== //
    private:
        darray3 	m_origin;				/**<Coordinate X,Y,Z of the origin of the octree in the physical domain*/
        double 		m_L;					/**<Side length of octree in the physical domain*/

        // =================================================================================== //
        // CONSTRUCTORS AND OPERATORS
        // =================================================================================== //

        // =================================================================================== //
        // METHODS
        // =================================================================================== //
        void	__reset();
    public:
#if BITPIT_ENABLE_MPI==1
        PabloUniform(std::string logfile = DEFAULT_LOG_FILE, MPI_Comm comm = MPI_COMM_WORLD);
        PabloUniform(uint8_t dim, std::string logfile = DEFAULT_LOG_FILE, MPI_Comm comm = MPI_COMM_WORLD);
        PabloUniform(double X, double Y, double Z, double L, uint8_t dim = 2, std::string logfile = DEFAULT_LOG_FILE, MPI_Comm comm = MPI_COMM_WORLD);
#else
        PabloUniform(std::string logfile = DEFAULT_LOG_FILE);
        PabloUniform(uint8_t dim, std::string logfile = DEFAULT_LOG_FILE);
        PabloUniform(double X, double Y, double Z, double L, uint8_t dim = 2, std::string logfile = DEFAULT_LOG_FILE);
#endif

        // =================================================================================== //
        // METHODS
        // =================================================================================== //
        void	reset();

        int		getDumpVersion() const;
        void	dump(std::ostream &stream, bool full = true);
        void	restore(std::istream &stream);

        // =================================================================================== //
        // BASIC GET/SET METHODS															   //
        // =================================================================================== //
        darray3		getOrigin() const;
        double		getX0() const;
        double		getY0() const;
        double		getZ0() const;
        double		getL() const;
        void		setL(double L);
        void		setOrigin(darray3 origin);
        double		levelToSize( uint8_t& level);

        // =================================================================================== //
        // INDEX BASED METHODS																   //
        // =================================================================================== //
        darray3 	getCoordinates(uint32_t idx);
        double 		getX(uint32_t idx);
        double 		getY(uint32_t idx);
        double 		getZ(uint32_t idx);
        double 		getSize(uint32_t idx);
        double 		getArea(uint32_t idx);
        double 		getVolume(uint32_t idx);
        void 		getCenter(uint32_t idx, darray3& center) const;
        darray3 	getCenter(uint32_t idx) const;
        darray3 	getFaceCenter(uint32_t idx, uint8_t iface) const;
        void 		getFaceCenter(uint32_t idx, uint8_t iface, darray3& center) const;
        darray3 	getNode(uint32_t idx, uint8_t inode);
        void 		getNode(uint32_t idx, uint8_t inode, darray3& node);
        void 		getNodes(uint32_t idx, darr3vector & nodes);
        darr3vector getNodes(uint32_t idx);
        void 		getNormal(uint32_t idx, uint8_t & iface, darray3 & normal);
        darray3 	getNormal(uint32_t idx, uint8_t & iface);

        // =================================================================================== //
        // POINTER BASED METHODS															   //
        // =================================================================================== //
        darray3 	getCoordinates(Octant* oct);
        double 		getX(Octant* oct);
        double 		getY(Octant* oct);
        double 		getZ(Octant* oct);
        double 		getSize(Octant* oct);
        double 		getArea(Octant* oct);
        double 		getVolume(Octant* oct);
        void 		getCenter(const Octant* oct, darray3& center) const;
        darray3 	getCenter(const Octant* oct) const;
        darray3 	getFaceCenter(const Octant* oct, uint8_t iface) const;
        void 		getFaceCenter(const Octant* oct, uint8_t iface, darray3& center) const;
        darray3 	getNode(Octant* oct, uint8_t inode);
        void 		getNode(Octant* oct, uint8_t inode, darray3& node);
        void 		getNodes(Octant* oct, darr3vector & nodes);
        darr3vector getNodes(Octant* oct);
        void 		getNormal(Octant* oct, uint8_t & iface, darray3 & normal);
        darray3 	getNormal(Octant* oct, uint8_t & iface);

        // =================================================================================== //
        // LOCAL TREE GET/SET METHODS														   //
        // =================================================================================== //
        double	 	getLocalMaxSize();
        double	 	getLocalMinSize();
        void 		getBoundingBox(darray3 & P0, darray3 & P1);

        // =================================================================================== //
        // INTERSECTION GET/SET METHODS														   //
        // =================================================================================== //
        double 		getSize(Intersection* inter);
        double 		getArea(Intersection* inter);
        darray3 	getCenter(const Intersection* inter) const;
        darr3vector getNodes(Intersection* inter);
        darray3 	getNormal(Intersection* inter);

        // =================================================================================== //
        // OTHER OCTANT BASED METHODS												    	   //
        // =================================================================================== //
        Octant* getPointOwner(darray3 point);
        uint32_t getPointOwnerIdx(darray3 point);
        Octant* getPointOwner(darray3 point, bool & isghost);
        uint32_t getPointOwnerIdx(darray3 point, bool & isghost);

        // =================================================================================== //
        // OTHER PARATREE BASED METHODS												    	   //
        // =================================================================================== //
        darray3 	getNodeCoordinates(uint32_t inode);

    };

}

#endif /* __BITPIT_PABLO_UNIFORM_HPP__ */
