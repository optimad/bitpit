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

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "bitpit_operators.hpp"

#include "PabloNonUniform.hpp"

namespace bitpit {

    // =================================================================================== //
    // NAME SPACES                                                                         //
    // =================================================================================== //
    using namespace std;

    // =================================================================================== //
    // CLASS IMPLEMENTATION                                                                //
    // =================================================================================== //

    // =================================================================================== //
    // CONSTRUCTORS AND OPERATORS                                                          //
    // =================================================================================== //
    /*! Default empty constructor of PabloNonUniform.
     * \param[in] logfile The file name for the log of this object. PABLO.log is the default value.
     */
#if BITPIT_ENABLE_MPI==1
    /*!
     * \param[in] comm The MPI communicator used by the parallel octree. MPI_COMM_WORLD is the default value.
     */
    PabloNonUniform::PabloNonUniform(std::string logfile, MPI_Comm comm):ParaTree(logfile,comm){
#else
    PabloNonUniform::PabloNonUniform(std::string logfile):ParaTree(logfile){
#endif
        __reset();
    }

    /*! Default constructor of PabloNonUniform.
     * It sets the Origin in (0,0,0) and side of length 1.
     * \param[in] dim The space dimension of the octree.
     * \param[in] logfile The file name for the log of this object. PABLO.log is the default value.
     */
#if BITPIT_ENABLE_MPI==1
    /*!
     * \param[in] comm The MPI communicator used by the parallel octree. MPI_COMM_WORLD is the default value.
     */
    PabloNonUniform::PabloNonUniform(uint8_t dim, std::string logfile, MPI_Comm comm):ParaTree(dim,logfile,comm){
#else
    PabloNonUniform::PabloNonUniform(uint8_t dim, std::string logfile):ParaTree(dim,logfile){
#endif
        __reset();
    };

    /*! Custom constructor of PabloNonUniform. Equivalent to PabloUniform
     * It sets the Origin in (X,Y,Z) and side of length L.
     * \param[in] X x-coordinate of the origin in physical domain,
     * \param[in] Y y-coordinate of the origin in physical domain,
     * \param[in] Z z-coordinate of the origin in physical domain,
     * \param[in] L Length of the sides in physical domain.
     * \param[in] dim The space dimension of the octree.
     * \param[in] logfile The file name for the log of this object. PABLO.log is the default value.
     */
#if BITPIT_ENABLE_MPI==1
    /*!
     * \param[in] comm The MPI communicator used by the parallel octree. MPI_COMM_WORLD is the default value.
     */
    PabloNonUniform::PabloNonUniform(double X, double Y, double Z, double L, uint8_t dim, std::string logfile, MPI_Comm comm):ParaTree(dim,logfile,comm){
#else
    PabloNonUniform::PabloNonUniform(double X, double Y, double Z, double L, uint8_t dim, std::string logfile):ParaTree(dim,logfile){
#endif
        __reset();

        setOrigin({{X, Y, Z}});
        if(dim == 2)
          setL({{L, L, 0.}});
        else if(dim == 3)
          setL({{L, L, L}});
    };

    /*! Custom constructor of PabloNonUniform.
     * It sets the Origin in (X,Y,Z) and sides of length (LX, LY, LZ).
     * \param[in] X x-coordinate of the origin in physical domain,
     * \param[in] Y y-coordinate of the origin in physical domain,
     * \param[in] Z z-coordinate of the origin in physical domain,
     * \param[in] LX Length of the sides in physical domain along x direction.
     * \param[in] LY Length of the sides in physical domain along y direction.
     * \param[in] LZ Length of the sides in physical domain along z direction.
     * \param[in] dim The space dimension of the octree.
     * \param[in] logfile The file name for the log of this object. PABLO.log is the default value.
     */

#if BITPIT_ENABLE_MPI==1
    /*!
     * \param[in] comm The MPI communicator used by the parallel octree. MPI_COMM_WORLD is the default value.
     */
    PabloNonUniform::PabloNonUniform(double X, double Y, double Z, double LX, double LY, double LZ, uint8_t dim, std::string logfile, MPI_Comm comm):ParaTree(dim,logfile,comm){
#else
    PabloNonUniform::PabloNonUniform(double X, double Y, double Z, double LX, double LY, double LZ, uint8_t dim, std::string logfile):ParaTree(dim,logfile){
#endif
        __reset();

        setOrigin({{X, Y, Z}});
        setL({{LX, LY, LZ}});
     
    };

    // =================================================================================== //
    // METHODS                                                                             //
    // =================================================================================== //

    /*! Reset the octree
     */
    void
    PabloNonUniform::reset(){
        ParaTree::reset();
        __reset();
    }

    /*! Internal function to reset the octree
     */
    void
    PabloNonUniform::__reset(){
        setOrigin({{0,0,0}});
        setL({{1., 1., 0.}});
    }

    /*! Get the version associated to the binary dumps.
     *
     *  \result The version associated to the binary dumps.
     */
    int
    PabloNonUniform::getDumpVersion() const
    {
        const int DUMP_VERSION = 1;

        return (DUMP_VERSION + ParaTree::getDumpVersion());
    }

    /*! Write the octree to the specified stream.
    *
    *  \param stream is the stream to write to
    *  \param full is the flag for a complete dump with mapping structureof last operation of the tree
    */
    void
    PabloNonUniform::dump(std::ostream &stream, bool full)
    {
        ParaTree::dump(stream, full);

        utils::binary::write(stream, m_origin[0]);
        utils::binary::write(stream, m_origin[1]);
        utils::binary::write(stream, m_origin[2]);
        utils::binary::write(stream, m_L[0]);
        utils::binary::write(stream, m_L[1]);
        utils::binary::write(stream, m_L[2]);
    }

    /*! Restore the octree from the specified stream.
    *
    *  \param stream is the stream to read from
    */
    void
    PabloNonUniform::restore(std::istream &stream)
    {
        ParaTree::restore(stream);

        darray3 origin;
        utils::binary::read(stream, origin[0]);
        utils::binary::read(stream, origin[1]);
        utils::binary::read(stream, origin[2]);
        setOrigin(origin);
        
        darray3 L;
        utils::binary::read(stream, L[0]);
        utils::binary::read(stream, L[1]);
        utils::binary::read(stream, L[2]);
        setL(L);
    }

    // =================================================================================== //
    // BASIC GET/SET METHODS								   //
    // =================================================================================== //
    /*! Get the coordinates of the origin of the octree.
     * \return Coordinates of the origin.
     */
    darray3
    PabloNonUniform::getOrigin() const {
        return m_origin;
    };

    /*! Get the coordinate X of the origin of the octree.
     * \return Coordinate X of the origin.
     */
    double
    PabloNonUniform::getX0() const {
        return m_origin[0];
    };

    /*! Get the coordinate Y of the origin of the octree.
     * \return Coordinate Y of the origin.
     */
    double
    PabloNonUniform::getY0() const {
        return m_origin[1];
    };

    /*! Get the coordinate Z of the origin of the octree.
     * \return Coordinate Z of the origin.
     */
    double
    PabloNonUniform::getZ0() const {
        return m_origin[2];
    };

    /*! Get the length of the domain.
     * \return Length of the octree.
     */
    darray3
    PabloNonUniform::getL() const {
        return m_L;
    };

    /*! Set the length of the domain.
     * \param[in] L Length of the octree.
     */
    void
    PabloNonUniform::setL(darray3 L){
        m_L[0] = L[0];
        m_L[1] = L[1];
        m_L[2] = L[2];
        
        if (getDim() == 2)
        {
          m_area.resize(4);
          m_area[0] = L[1];
          m_area[1] = L[1];
          m_area[2] = L[0];
          m_area[3] = L[0];
          m_volume = L[0] * L[1];
        }
        if (getDim() == 3)
        {
          m_area.resize(6);
          m_area[0] = L[1] * L[2];
          m_area[1] = L[1] * L[2];
          m_area[2] = L[0] * L[2];
          m_area[3] = L[0] * L[2];
          m_area[4] = L[0] * L[1];
          m_area[5] = L[0] * L[1];
          m_volume = L[0] * L[1] * L[2];
        }
    };

    /*! Set the origin of the domain.
     * \param[in] origin Origin of the octree.
     */
    void
    PabloNonUniform::setOrigin(darray3 origin){
        m_origin = origin;
    };

    /*! Get the size of an octant corresponding to a target level.
     * \param[in] level Input level.
     * \return Size of an octant of input level.
     */
    double
    PabloNonUniform::levelToSize(uint8_t & level, uint8_t dir) {
        double size = ParaTree::levelToSize(level);
        return m_L[dir] *size;
    }

    // =================================================================================== //
    // INDEX BASED METHODS								   //
    // =================================================================================== //
    /*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
     * \param[in] idx Local index of target octant.
     * \return Coordinates X,Y,Z of node 0.
     */
    darray3
    PabloNonUniform::getCoordinates(uint32_t idx) const {
        darray3 coords, coords_;
        coords_ = ParaTree::getCoordinates(idx);
        for (int i=0; i<3; i++){
            coords[i] = m_origin[i] + m_L[i] * coords_[i];
        }
        return coords;
    };

    /*! Get the coordinate X of an octant, i.e. the coordinates of its node 0.
     * \param[in] idx Local index of target octant.
     * \return Coordinate X of node 0.
     */
    double
    PabloNonUniform::getX(uint32_t idx) const {
        double X, X_;
        X_ = ParaTree::getX(idx);
        X = m_origin[0] + m_L[0] * X_;
        return X;
    };

    /*! Get the coordinate Y of an octant, i.e. the coordinates of its node 0.
     * \param[in] idx Local index of target octant.
     * \return Coordinate Y of node 0.
     */
    double
    PabloNonUniform::getY(uint32_t idx) const {
        double Y, Y_;
        Y_ = ParaTree::getY(idx);
        Y = m_origin[1] + m_L[1] * Y_;
        return Y;
    };

    /*! Get the coordinate Z of an octant, i.e. the coordinates of its node 0.
     * \param[in] idx Local index of target octant.
     * \return Coordinate Z of node 0.
     */
    double
    PabloNonUniform::getZ(uint32_t idx) const {
        double Z, Z_;
        Z_ = ParaTree::getZ(idx);
        Z = m_origin[2] + m_L[2] * Z_;
        return Z;
    };

    /*! Get the size of an octant, i.e. the side length.
     * \param[in] idx Local index of target octant.
     * \return Size of octant.
     */
    double
    PabloNonUniform::getSize(uint32_t idx, uint8_t dir) const {
        return m_L[dir] * ParaTree::getSize(idx);
    };

    /*! Get the area of an octant (for 2D case the same value of getSize).
     * \param[in] idx Local index of target octant.
     * \return Area of octant.
     */
    double
    PabloNonUniform::getArea(uint32_t idx, uint8_t iface) const {
        return m_area[iface] * ParaTree::getArea(idx);
    };

    /*! Get the volume of an octant.
     * \param[in] idx Local index of target octant.
     * \return Volume of octant.
     */
    double
    PabloNonUniform::getVolume(uint32_t idx) const {
        return m_volume * ParaTree::getVolume(idx);
    };

    /*! Get the coordinates of the center of an octant.
     * \param[in] idx Local index of target octant.
     * \param[out] center Coordinates of the center of octant.
     */
    void
    PabloNonUniform::getCenter(uint32_t idx, darray3& center) const {
        darray3 center_ = ParaTree::getCenter(idx);
        for (int i=0; i<3; i++){
            center[i] = m_origin[i] + m_L[i] * center_[i];
        }
    };

    /*! Get the coordinates of the center of an octant.
     * \param[in] idx Local index of target octant.
     * \return center Coordinates of the center of octant.
     */
    darray3
    PabloNonUniform::getCenter(uint32_t idx) const {
        darray3 center, center_ = ParaTree::getCenter(idx);
        for (int i=0; i<3; i++){
            center[i] = m_origin[i] + m_L[i] * center_[i];
        }
        return center;
    };

    /*! Get the coordinates of the center of a face of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] iface Index of the target face.
     * \param[out] center Coordinates of the center of the iface-th face of octant.
     */
    void
    PabloNonUniform::getFaceCenter(uint32_t idx, uint8_t iface, darray3& center) const {
        darray3 center_ = ParaTree::getFaceCenter(idx, iface);
        for (int i=0; i<3; i++){
            center[i] = m_origin[i] + m_L[i] * center_[i];
        }
    };

    /*! Get the coordinates of the center of a face of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] iface Index of the target face.
     * \return center Coordinates of the center of the iface-th face of octant.
     */
    darray3
    PabloNonUniform::getFaceCenter(uint32_t idx, uint8_t iface) const {
        darray3 center, center_ = ParaTree::getFaceCenter(idx, iface);
        for (int i=0; i<3; i++){
            center[i] = m_origin[i] + m_L[i] * center_[i];
        }
        return center;
    };

    /*! Get the coordinates of a node of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] inode Index of the target node.
     * \return Coordinates of of the inode-th node of octant.
     */
    darray3
    PabloNonUniform::getNode(uint32_t idx, uint8_t inode) const {
        darray3 node, node_ = ParaTree::getNode(idx, inode);
        for (int i=0; i<3; i++){
            node[i] = m_origin[i] + m_L[i] * node_[i];
        }
        return node;
    };

    /*! Get the coordinates of a node of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] inode Index of the target node.
     * \param[out] node Coordinates of of the inode-th node of octant.
     */
    void
    PabloNonUniform::getNode(uint32_t idx, uint8_t inode, darray3& node) const {
        darray3 node_ = ParaTree::getNode(idx, inode);
        for (int i=0; i<3; i++){
            node[i] = m_origin[i] + m_L[i] * node_[i];
        }
    };

    /*! Get the coordinates of the nodes of an octant.
     * \param[in] idx Local index of target octant.
     * \param[out] nodes Coordinates of the nodes of octant.
     */
    void
    PabloNonUniform::getNodes(uint32_t idx, darr3vector & nodes) const {
        darray3vector nodes_ = ParaTree::getNodes(idx);
        nodes.resize(ParaTree::getNnodes());
        for (int j=0; j<ParaTree::getNnodes(); j++){
            for (int i=0; i<3; i++){
                nodes[j][i] = m_origin[i] + m_L[i] * nodes_[j][i];
            }
        }
    };

    /*! Get the coordinates of the nodes of an octant.
     * \param[in] idx Local index of target octant.
     * \return nodes Coordinates of the nodes of octant.
     */
    darr3vector
    PabloNonUniform::getNodes(uint32_t idx) const {
        darray3vector nodes, nodes_ = ParaTree::getNodes(idx);
        nodes.resize(ParaTree::getNnodes());
        for (int j=0; j<ParaTree::getNnodes(); j++){
            for (int i=0; i<3; i++){
                nodes[j][i] = m_origin[i] + m_L[i] * nodes_[j][i];
            }
        }
        return nodes;
    };

    /*! Get the normal of a face of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] iface Index of the face for normal computing.
     * \param[out] normal Coordinates of the normal of face.
     */
    void
    PabloNonUniform::getNormal(uint32_t idx, uint8_t iface, darray3 & normal) const {
        ParaTree::getNormal(idx, iface, normal);
    }

    /*! Get the normal of a face of an octant.
     * \param[in] idx Local index of target octant.
     * \param[in] iface Index of the face for normal computing.
     * \return normal Coordinates of the normal of face.
     */
    darray3
    PabloNonUniform::getNormal(uint32_t idx, uint8_t iface) const {
        return ParaTree::getNormal(idx, iface);
    }

    // =================================================================================== //
    // POINTER BASED METHODS     							   //
    // =================================================================================== //
    /*! Get the coordinates of an octant, i.e. the coordinates of its node 0.
     * \param[in] oct Pointer to the target octant
     * \return Coordinates of node 0.
     */
    darray3
    PabloNonUniform::getCoordinates(const Octant* oct) const {
        darray3 coords, coords_;
        coords_ = ParaTree::getCoordinates(oct);
        for (int i=0; i<3; i++){
            coords[i] = m_origin[i] + m_L[i] * coords_[i];
        }
        return coords;
    };

    /*! Get the coordinate X of an octant, i.e. the coordinates of its node 0.
     * \param[in] oct Pointer to the target octant
     * \return Coordinate X of node 0.
     */
    double
    PabloNonUniform::getX(const Octant* oct) const {
        double X, X_;
        X_ = ParaTree::getX(oct);
        X = m_origin[0] + m_L[0] * X_;
        return X;
    };

    /*! Get the coordinate Y of an octant, i.e. the coordinates of its node 0.
     * \param[in] oct Pointer to the target octant
     * \return Coordinate Y of node 0.
     */
    double
    PabloNonUniform::getY(const Octant* oct) const {
        double Y, Y_;
        Y_ = ParaTree::getY(oct);
        Y = m_origin[1] + m_L[1] * Y_;
        return Y;
    };

    /*! Get the coordinate Z of an octant, i.e. the coordinates of its node 0.
     * \param[in] oct Pointer to the target octant
     * \return Coordinate Z of node 0.
     */
    double
    PabloNonUniform::getZ(const Octant* oct) const {
        double Z, Z_;
        Z_ = ParaTree::getZ(oct);
        Z = m_origin[2] + m_L[2] * Z_;
        return Z;
    };

    /*! Get the size of an octant, i.e. the side length.
     * \param[in] oct Pointer to the target octant
     * \return Size of octant.
     */
    double
    PabloNonUniform::getSize(const Octant* oct, uint8_t dir) const {
        return m_L[dir] * ParaTree::getSize(oct);
    };

    /*! Get the area of an octant (for 2D case the same value of getSize).
     * \param[in] oct Pointer to the target octant
     * \return Area of octant.
     */
    double
    PabloNonUniform::getArea(const Octant* oct, uint8_t iface) const {
        return m_area[iface] * ParaTree::getArea(oct);
    };

    /*! Get the volume of an octant.
     * \param[in] oct Pointer to the target octant
     * \return Volume of octant.
     */
    double
    PabloNonUniform::getVolume(const Octant* oct) const {
        return m_volume * ParaTree::getVolume(oct);
    };

    /*! Get the coordinates of the center of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[out] center Coordinates of the center of octant.
     */
    void
    PabloNonUniform::getCenter(const Octant* oct, darray3& center) const {
        darray3 center_ = ParaTree::getCenter(oct);
        for (int i=0; i<3; i++){
            center[i] = m_origin[i] + m_L[i] * center_[i];
        }
    };

    /*! Get the coordinates of the center of an octant.
     * \param[in] oct Pointer to the target octant
     * \return center Coordinates of the center of octant.
     */
    darray3
    PabloNonUniform::getCenter(const Octant* oct) const {
        darray3 center, center_ = ParaTree::getCenter(oct);
        for (int i=0; i<3; i++){
            center[i] = m_origin[i] + m_L[i] * center_[i];
        }
        return center;
    };

    /*! Get the coordinates of the center of a face of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] iface Index of the target face.
     * \param[out] center Coordinates of the center of the iface-th face af octant.
     */
    void
    PabloNonUniform::getFaceCenter(const Octant* oct, uint8_t iface, darray3& center) const {
        darray3 center_ = ParaTree::getFaceCenter(oct, iface);
        for (int i=0; i<3; i++){
            center[i] = m_origin[i] + m_L[i] * center_[i];
        }
    };

    /*! Get the coordinates of the center of a face of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] iface Index of the target face.
     * \return center Coordinates of the center of the iface-th face af octant.
     */
    darray3
    PabloNonUniform::getFaceCenter(const Octant* oct, uint8_t iface) const {
        darray3 center, center_ = ParaTree::getFaceCenter(oct, iface);
        for (int i=0; i<3; i++){
            center[i] = m_origin[i] + m_L[i] * center_[i];
        }
        return center;
    };

    /*! Get the coordinates of single node of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] inode Index of the target node.
     * \return Coordinates of the center of the inode-th of octant.
     */
    darray3
    PabloNonUniform::getNode(const Octant* oct, uint8_t inode) const {
        darray3 node, node_ = ParaTree::getNode(oct, inode);
        for (int i=0; i<3; i++){
            node[i] = m_origin[i] + m_L[i] * node_[i];
        }
        return node;
    };

    /*! Get the coordinates of the center of a face of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] inode Index of the target node.
     * \param[out] node Coordinates of the center of the inode-th of octant.
     */
    void
    PabloNonUniform::getNode(const Octant* oct, uint8_t inode, darray3& node) const {
        darray3 node_ = ParaTree::getNode(oct, inode);
        for (int i=0; i<3; i++){
            node[i] = m_origin[i] + m_L[i] * node_[i];
        }
    };

    /*! Get the coordinates of the nodes of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[out] nodes Coordinates of the nodes of octant.
     */
    void
    PabloNonUniform::getNodes(const Octant* oct, darr3vector & nodes) const {
        darray3vector nodes_ = ParaTree::getNodes(oct);
        nodes.resize(ParaTree::getNnodes());
        for (int j=0; j<ParaTree::getNnodes(); j++){
            for (int i=0; i<3; i++){
                nodes[j][i] = m_origin[i] + m_L[i] * nodes_[j][i];
            }
        }
    };

    /*! Get the coordinates of the nodes of an octant.
     * \param[in] oct Pointer to the target octant
     * \return nodes Coordinates of the nodes of octant.
     */
    darr3vector
    PabloNonUniform::getNodes(const Octant* oct) const {
        darray3vector nodes, nodes_ = ParaTree::getNodes(oct);
        nodes.resize(ParaTree::getNnodes());
        for (int j=0; j<ParaTree::getNnodes(); j++){
            for (int i=0; i<3; i++){
                nodes[j][i] = m_origin[i] + m_L[i] * nodes_[j][i];
            }
        }
        return nodes;
    };

    /*! Get the normal of a face of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] iface Index of the face for normal computing.
     * \param[out] normal Coordinates of the normal of face.
     */
    void
    PabloNonUniform::getNormal(const Octant* oct, uint8_t iface, darray3 & normal) const {
        ParaTree::getNormal(oct, iface, normal);
    }

    /*! Get the normal of a face of an octant.
     * \param[in] oct Pointer to the target octant
     * \param[in] iface Index of the face for normal computing.
     * \return normal Coordinates of the normal of face.
     */
    darray3
    PabloNonUniform::getNormal(const Octant* oct, uint8_t iface) const {
        return ParaTree::getNormal(oct, iface);
    }

    // =================================================================================== //
    // LOCAL TREE GET/SET METHODS							   //
    // =================================================================================== //
    /*! Get the local current maximum size of the octree.
     * \return Local current maximum size of the local partition of the octree.
     */
    double
    PabloNonUniform::getLocalMaxSize() const {
        return *std::max_element(m_L.begin(), m_L.end()) * ParaTree::getLocalMaxSize();
    };

    /*! Get the local current minimum size of the octree.
     * \return Local current minimum size of the local partition of the octree.
     */
    double
    PabloNonUniform::getLocalMinSize() const {
        double minL = m_L[0];
        if (getDim() == 2)
          minL = std::min(m_L[0], m_L[1]);
        else if (getDim() == 3)
         minL = *std::min_element(m_L.begin(), m_L.end()); 

        return minL * ParaTree::getLocalMinSize();
    };


    /*! Get the coordinates of the extreme points of a bounding box containing the local tree
     *  \param[out] P0 Array with coordinates of the first point (lowest coordinates);
     *  \param[out] P1 Array with coordinates of the last point (highest coordinates).
     */
    void
    PabloNonUniform::getBoundingBox(darray3 & P0, darray3 & P1) const {
        // If there are no octants the bounding box is empty
        uint32_t nocts = ParaTree::getNumOctants();
        if (nocts == 0) {
            P0 = getOrigin();
            P1 = P0;

            return;
        }

        // If the octree is serial we can evaluate the bounding box easily
        // otherwise we need to scan all the octants
        if (getSerial()) {
            P0 = getOrigin();
            P1 = P0;
            darray3 L = getL();
            for (int i=0; i<ParaTree::getDim(); i++){
                P1[i] += L[i];
            }

            return;
        }

        // If the octree is parallel we need to scan all the octants
        darray3		cnode0, cnode1;

        uint32_t	id = 0;
        uint8_t 	nnodes = ParaTree::getNnodes();

        P0 = getNode(id, 0);
        P1 = getNode(nocts-1, nnodes-1);

        for (id=0; id<nocts; id++){
            cnode0 = getNode(id, 0);
            cnode1 = getNode(id, nnodes-1);
            for (int i=0; i<ParaTree::getDim(); i++){
                P0[i] = min(P0[i], cnode0[i]);
                P1[i] = max(P1[i], cnode1[i]);
            }
        }
    };


    // =================================================================================== //
    // INTERSECTION GET/SET METHODS							   //
    // =================================================================================== //
    /*! Get the size of an intersection.
     * \param[in] inter Pointer to target intersection.
     * \return Size of intersection.
     */
    double
    PabloNonUniform::getSize(const Intersection* inter, uint8_t dir) const {
        return m_L[dir] * ParaTree::getSize(inter);
    };

    /*! Get the area of an intersection (for 2D case the same value of getSize).
     * \param[in] inter Pointer to target intersection.
     * \return Area of intersection.
     */
    double
    PabloNonUniform::getArea(const Intersection* inter, uint8_t iface) const {
        return m_area[iface] * ParaTree::getArea(inter);
    };

    /*! Get the coordinates of the center of an intersection.
     * \param[in] inter Pointer to target intersection.
     * \return Coordinates of the center of intersection.
     */
    darray3
    PabloNonUniform::getCenter(const Intersection* inter) const {
        darray3 center = ParaTree::getCenter(inter);
        for (int i=0; i<3; i++){
            center[i] = m_origin[i] + m_L[i] * center[i];
        }
        return center;
    }

    /*! Get the coordinates of the nodes of an intersection.
     * \param[in] inter Pointer to target intersection.
     * \return Coordinates of the nodes of intersection.
     */
    darr3vector
    PabloNonUniform::getNodes(const Intersection* inter) const {
        darr3vector nodes, nodes_ = ParaTree::getNodes(inter);
        nodes.resize(ParaTree::getNnodesperface());
        for (int j=0; j<ParaTree::getNnodesperface(); j++){
            for (int i=0; i<3; i++){
                nodes[j][i] = m_origin[i] + m_L[i] * nodes_[j][i];
            }
        }
        return nodes;
    }

    /*! Get the normal of an intersection.
     * \param[in] inter Pointer to target intersection.
     * \return Coordinates of the normal of intersection.
     */
    darray3
    PabloNonUniform::getNormal(const Intersection* inter) const {
        return ParaTree::getNormal(inter);
    }

    // =================================================================================== //
    // OTHER OCTANT BASED METHODS						    	   //
    // =================================================================================== //
    /** Get the octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \return Pointer to octant owner of target point
     * (=NULL if point is outside of the domain).
     */
    Octant* PabloNonUniform::getPointOwner(darray3 point){
        for (int i=0; i<3; i++){
            point[i] = (point[i] - m_origin[i])/m_L[i];
        }
        return ParaTree::getPointOwner(point);
    };

    /** Get the octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \param[out] isghost Boolean flag, true if the octant found is ghost
     * \return Index of octant owner of target point (max uint32_t representable if point outside of the ghosted domain).
     */
    Octant* PabloNonUniform::getPointOwner(darray3 point, bool & isghost){
        for (int i=0; i<3; i++){
            point[i] = (point[i] - m_origin[i])/m_L[i];
        }
        return ParaTree::getPointOwner(point,isghost);
    };


    /** Get the octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \return Index of octant owner of target point
     * (max uint32_t representable if point outside of the domain).
     */
    uint32_t
    PabloNonUniform::getPointOwnerIdx(darray3 point) const {
        for (int i=0; i<3; i++){
            point[i] = (point[i] - m_origin[i])/m_L[i];
        }
        return ParaTree::getPointOwnerIdx(point);
    };
    
    /** Get the octant owner of an input point.
     * \param[in] point Coordinates of target point.
     * \param[out] isghost Boolean flag, true if the octant found is ghost
     * \return Index of octant owner of target point (max uint32_t representable if point outside of the ghosted domain).
     */
    uint32_t
    PabloNonUniform::getPointOwnerIdx(darray3 point, bool & isghost) const {
        for (int i=0; i<3; i++){
            point[i] = (point[i] - m_origin[i])/m_L[i];
        }
        return ParaTree::getPointOwnerIdx(point,isghost);
    };
    
    /** Get the octant owner rank of an input point.
     * \param[in] point Coordinates of target point.
     * \return Owner rank of target point (negative if out of global domain).
     */
    int PabloNonUniform::getPointOwnerRank(darray3 point){
        for (int i=0; i<3; i++){
            point[i] = (point[i] - m_origin[i])/m_L[i];
        }
        return ParaTree::getPointOwnerRank(point);
    };
    
    // =================================================================================== //
    // OTHER PARATREE BASED METHODS       					    	   //
    // =================================================================================== //
    /** Get the physical coordinates of a node
     * \param[in] inode Local index of node
     * \return Vector with the coordinates of the node.
     */
    darray3
    PabloNonUniform::getNodeCoordinates(uint32_t inode) const {
        darray3 node = ParaTree::getNodeCoordinates(inode);
        for (int i=0; i<3; i++){
            node[i] = m_origin[i] + m_L[i] * node[i];
        }
        return node;
    }

    // =================================================================================== //
    // TESTING OUTPUT METHODS		         				    	   //
    // =================================================================================== //
    /** Write the physical octree mesh in .vtu format in a user-defined file.
     * If the connectivity is not stored, the method temporary computes it.
     * If the connectivity of ghost octants is already computed, the method writes the ghosts on file.
     * \param[in] filename Name of output file (PABLO will add the total number of processes p000# and the current rank s000#).
     */
    void
    PabloNonUniform::write(string filename) {

        if (getConnectivity().size() == 0) {
            computeConnectivity();
        }

        stringstream name;
        name << "s" << std::setfill('0') << std::setw(4) << getNproc() << "-p" << std::setfill('0') << std::setw(4) << getRank() << "-" << filename << ".vtu";

        ofstream out(name.str().c_str());
        if(!out.is_open()){
            stringstream ss;
            ss << filename << "*.vtu cannot be opened and it won't be written." << endl;
            getLog() << ss.str();
            return;
        }
        int nofNodes = getNumNodes();
        int nofOctants = getNumOctants();
        int nofGhosts = getNumGhosts();
        int nofAll = nofGhosts + nofOctants;
        out << "<?xml version=\"1.0\"?>" << endl
            << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
            << "  <UnstructuredGrid>" << endl
            << "    <Piece NumberOfCells=\"" << getConnectivity().size() + getGhostConnectivity().size() << "\" NumberOfPoints=\"" << getNumNodes() << "\">" << endl;
        out << "      <Points>" << endl
            << "        <DataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\""<< 3 <<"\" format=\"ascii\">" << endl
            << "          " << std::fixed;
        for(int i = 0; i < nofNodes; i++)
            {
                const std::array<double,3> & nodeCoordinates = getNodeCoordinates(i);
                for(int j = 0; j < 3; ++j){
                    out << std::setprecision(6) << nodeCoordinates[j] << " ";
                }
                if((i+1)%4==0 && i!=nofNodes-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "      </Points>" << endl
            << "      <Cells>" << endl
            << "        <DataArray type=\"UInt64\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          ";
        for(int i = 0; i < nofOctants; i++)
            {
                for(int j = 0; j < getNnodes(); j++)
                    {
                        int jj = j;
                        if (getDim()==2){
                            if (j<2){
                                jj = j;
                            }
                            else if(j==2){
                                jj = 3;
                            }
                            else if(j==3){
                                jj = 2;
                            }
                        }
                        out << getConnectivity()[i][jj] << " ";
                    }
                if((i+1)%3==0 && i!=nofOctants-1)
                    out << endl << "          ";
            }
        for(int i = 0; i < nofGhosts; i++)
            {
                for(int j = 0; j < getNnodes(); j++)
                    {
                        int jj = j;
                        if (getDim()==2){
                            if (j<2){
                                jj = j;
                            }
                            else if(j==2){
                                jj = 3;
                            }
                            else if(j==3){
                                jj = 2;
                            }
                        }
                        out << getGhostConnectivity()[i][jj] << " ";
                    }
                if((i+1)%3==0 && i!=nofGhosts-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "        <DataArray type=\"UInt64\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          ";
        for(int i = 0; i < nofAll; i++)
            {
                out << (i+1)*getNnodes() << " ";
                if((i+1)%12==0 && i!=nofAll-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          ";
        for(int i = 0; i < nofAll; i++)
            {
                int type;
                type = 5 + (getDim()*2);
                out << type << " ";
                if((i+1)%12==0 && i!=nofAll-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "      </Cells>" << endl
            << "    </Piece>" << endl
            << "  </UnstructuredGrid>" << endl
            << "</VTKFile>" << endl;


        if(getRank() == 0){
            name.str("");
            name << "s" << std::setfill('0') << std::setw(4) << getNproc() << "-" << filename << ".pvtu";
            ofstream pout(name.str().c_str());
            if(!pout.is_open()){
                stringstream ss;
                ss << filename << "*.pvtu cannot be opened and it won't be written." << endl;
                getLog() << ss.str();
                return;
            }

            pout << "<?xml version=\"1.0\"?>" << endl
                 << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
                 << "  <PUnstructuredGrid GhostLevel=\"0\">" << endl
                 << "    <PPointData>" << endl
                 << "    </PPointData>" << endl
                 << "    <PCellData Scalars=\"\">" << endl;
            pout << "    </PCellData>" << endl
                 << "    <PPoints>" << endl
                 << "      <PDataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\"/>" << endl
                 << "    </PPoints>" << endl;
            for(int i = 0; i < getNproc(); i++)
                pout << "    <Piece Source=\"s" << std::setw(4) << std::setfill('0') << getNproc() << "-p" << std::setw(4) << std::setfill('0') << i << "-" << filename << ".vtu\"/>" << endl;
            pout << "  </PUnstructuredGrid>" << endl
                 << "</VTKFile>";

            pout.close();

        }
#if BITPIT_ENABLE_MPI==1
        if (isCommSet()) {
            MPI_Barrier(getComm());
        }
#endif

    }

    /** Write the physical octree mesh in .vtu format with data for test in a user-defined file.
     * If the connectivity is not stored, the method temporary computes it.
     * The method doesn't write the ghosts on file.
     * \param[in] filename Name of output file (PABLO will add the total number of processes p000# and the current rank s000#).
     * \param[in] data Vector of double with user data.
     */
    void
    PabloNonUniform::writeTest(string filename, vector<double> data) {

        if (getConnectivity().size() == 0) {
            computeConnectivity();
        }

        stringstream name;
        name << "s" << std::setfill('0') << std::setw(4) << getNproc() << "-p" << std::setfill('0') << std::setw(4) << getRank() << "-" << filename << ".vtu";

        ofstream out(name.str().c_str());
        if(!out.is_open()){
            stringstream ss;
            ss << filename << "*.vtu cannot be opened and it won't be written.";
            getLog() << ss.str();
            return;
        }
        int nofNodes = getNumNodes();
        int nofOctants = getNumOctants();
        int nofAll = nofOctants;
        out << "<?xml version=\"1.0\"?>" << endl
            << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
            << "  <UnstructuredGrid>" << endl
            << "    <Piece NumberOfCells=\"" << getNumOctants() << "\" NumberOfPoints=\"" << getNumNodes() << "\">" << endl;
        out << "      <CellData Scalars=\"Data\">" << endl;
        out << "      <DataArray type=\"Float64\" Name=\"Data\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          " << std::fixed;
        int ndata = getNumOctants();
        for(int i = 0; i < ndata; i++)
            {
                out << std::setprecision(6) << data[i] << " ";
                if((i+1)%4==0 && i!=ndata-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "      </CellData>" << endl
            << "      <Points>" << endl
            << "        <DataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\""<< 3 <<"\" format=\"ascii\">" << endl
            << "          " << std::fixed;
        for(int i = 0; i < nofNodes; i++)
            {
            const std::array<double,3> & nodeCoordinates = getNodeCoordinates(i);
            for(int j = 0; j < 3; ++j){
                if (j==0) out << std::setprecision(6) << m_origin[0] + m_L[0] * nodeCoordinates[j] << " ";
                if (j==1) out << std::setprecision(6) << m_origin[1] + m_L[1] * nodeCoordinates[j] << " ";
                if (j==2) out << std::setprecision(6) << m_origin[2] + m_L[2] * nodeCoordinates[j] << " ";
            }
                if((i+1)%4==0 && i!=nofNodes-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "      </Points>" << endl
            << "      <Cells>" << endl
            << "        <DataArray type=\"UInt64\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          ";
        for(int i = 0; i < nofOctants; i++)
            {
                for(int j = 0; j < getNnodes(); j++)
                    {
                        int jj = j;
                        if (getDim()==2){
                            if (j<2){
                                jj = j;
                            }
                            else if(j==2){
                                jj = 3;
                            }
                            else if(j==3){
                                jj = 2;
                            }
                        }
                        out << getConnectivity()[i][jj] << " ";
                    }
                if((i+1)%3==0 && i!=nofOctants-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "        <DataArray type=\"UInt64\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          ";
        for(int i = 0; i < nofAll; i++)
            {
                out << (i+1)*getNnodes() << " ";
                if((i+1)%12==0 && i!=nofAll-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl
            << "          ";
        for(int i = 0; i < nofAll; i++)
            {
                int type;
                type = 5 + (getDim()*2);
                out << type << " ";
                if((i+1)%12==0 && i!=nofAll-1)
                    out << endl << "          ";
            }
        out << endl << "        </DataArray>" << endl
            << "      </Cells>" << endl
            << "    </Piece>" << endl
            << "  </UnstructuredGrid>" << endl
            << "</VTKFile>" << endl;


        if(getRank() == 0){
            name.str("");
            name << "s" << std::setfill('0') << std::setw(4) << getNproc() << "-" << filename << ".pvtu";
            ofstream pout(name.str().c_str());
            if(!pout.is_open()){
                stringstream ss;
                ss << filename << "*.pvtu cannot be opened and it won't be written." << endl;
                getLog() << ss.str();
                return;
            }

            pout << "<?xml version=\"1.0\"?>" << endl
                 << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << endl
                 << "  <PUnstructuredGrid GhostLevel=\"0\">" << endl
                 << "    <PPointData>" << endl
                 << "    </PPointData>" << endl
                 << "    <PCellData Scalars=\"Data\">" << endl
                 << "      <PDataArray type=\"Float64\" Name=\"Data\" NumberOfComponents=\"1\"/>" << endl
                 << "    </PCellData>" << endl
                 << "    <PPoints>" << endl
                 << "      <PDataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\"/>" << endl
                 << "    </PPoints>" << endl;
            for(int i = 0; i < getNproc(); i++)
                pout << "    <Piece Source=\"s" << std::setw(4) << std::setfill('0') << getNproc() << "-p" << std::setw(4) << std::setfill('0') << i << "-" << filename << ".vtu\"/>" << endl;
            pout << "  </PUnstructuredGrid>" << endl
                 << "</VTKFile>";

            pout.close();

        }
#if BITPIT_ENABLE_MPI==1
        MPI_Barrier(getComm());
#endif

    }

}
