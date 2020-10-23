/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
#include "bitpit_common.hpp"

#include "LocalTree.hpp"
#include <map>
#include <unordered_map>

namespace bitpit {

    // =================================================================================== //
    // NAME SPACES                                                                         //
    // =================================================================================== //
    using namespace std;

    // =================================================================================== //
    // CLASS IMPLEMENTATION                                                                    //
    // =================================================================================== //

    // =================================================================================== //
    // CONSTRUCTORS AND OPERATORS
    // =================================================================================== //

    /*!Defaut constructor.
     */
    LocalTree::LocalTree() {
        initialize();
        reset(false);
    };

    /*!Dimensional and default constructor.
     * \param[in] dim Space dimension of octree.
     */
    LocalTree::LocalTree(uint8_t dim){
        initialize(dim);
        reset(true);
    };

    // =================================================================================== //
    // METHODS
    // =================================================================================== //

    // =================================================================================== //
    // BASIC GET/SET METHODS
    // =================================================================================== //

    /*!Get the Morton number of first descentant octant of the octree.
     * \return Constant reference to the first descendant of the octree.
     */
    uint64_t
    LocalTree::getFirstDescMorton() const{
        return m_firstDescMorton;
    };

    /*!Get the Morton number of last descentant octant of the octree.
     * \return Constant reference to the last descendant of the octree.
     */
    uint64_t
    LocalTree::getLastDescMorton() const{
        return m_lastDescMorton;
    };

    /*! Get the number of the ghosts for the local partition of the tree.
     * \return Number of ghosts.
     */
    uint32_t
    LocalTree::getNumGhosts() const{
        return m_sizeGhosts;
    };

    /*! Get the number of the octants in the local tree.
     * \return Number of local octants.
     */
    uint32_t
    LocalTree::getNumOctants() const{
        return m_sizeOctants;
    };

    /*! Get max depth reached in local tree
     *  If the tree is empty a negative number is returned.
     * \return Max depth in local partition of the octree.
     */
    int8_t
    LocalTree::getLocalMaxDepth() const{
        return m_localMaxDepth;
    };

    /** Get refinement/coarsening marker for idx-th octant
     * \param[in] idx Local index of the target octant.
     * \return Marker of the octant.
     */
    int8_t
    LocalTree::getMarker(int32_t idx) const {
        return m_octants[idx].getMarker();
    };

    /** Get refinement/coarsening marker for idx-th octant
     * \param[in] idx Local index of the target octant.
     * \return Level of the octant.
     */
    uint8_t
    LocalTree::getLevel(int32_t idx) const {
        return m_octants[idx].getLevel();
    };

    /** Compute the Morton index of the idx-th octant (without level).
     * \param[in] idx Local index of the target octant.
     * \return Morton index of the octant.
     */
    uint64_t
    LocalTree::getMorton(int32_t idx) const {
        return m_octants[idx].getMorton();
    };

    /** Compute the persistent XYZ key of the specified node of an octant
     * (without level).
     * \param[in] idx Local index of the target octant.
     * \param[in] inode Index of the target node.
     * \return persistent XYZ key of the node.
     */
    uint64_t
    LocalTree::computeNodePersistentKey(int32_t idx, uint8_t inode) const {
        return m_octants[idx].computeNodePersistentKey(inode);
    };

    /** Get refinement/coarsening marker for idx-th ghost octant
     * \param[in] idx Local index of the target ghost octant.
     * \return Level of the ghost octant.
     */
    uint8_t
    LocalTree::getGhostLevel(int32_t idx) const {
        return m_ghosts[idx].getLevel();
    };

    /** Compute the Morton index of the idx-th ghost octant (without level).
     * \param[in] idx Local index of the target octant.
     * \return Morton index of the octant.
     */
    uint64_t
    LocalTree::computeGhostMorton(int32_t idx) const {
        return m_ghosts[idx].getMorton();
    };

    /** Compute the persistent XYZ key of the specified node of a ghost octant
     * (without level).
     * \param[in] idx Local index of the target octant.
     * \param[in] inode Index of the target node.
     * \return persistent XYZ key of the node.
     */
    uint64_t
    LocalTree::computeGhostNodePersistentKey(int32_t idx, uint8_t inode) const {
        return m_ghosts[idx].computeNodePersistentKey(inode);
    };

    /** Get if balancing-blocked idx-th octant
     * \param[in] idx Local index of the target octant.
     * \return Has the octant to be balanced?
     */
    bool
    LocalTree::getBalance(int32_t idx) const{
        return m_octants[idx].getBalance();
    };

    /*! Get the codimension for 2:1 balancing
     * \return Maximum codimension of the entity through which the 2:1 balance is performed.
     */
    uint8_t
    LocalTree::getBalanceCodim() const{
        return m_balanceCodim;
    };

    /** Set refinement/coarsening marker for idx-th octant
     * \param[in] idx Local index of the target octant.
     * \param[in] marker Refinement marker for the target octant.
     */
    void
    LocalTree::setMarker(int32_t idx, int8_t marker){
        m_octants[idx].setMarker(marker);
    };

    /** Set if balancing-blocked idx-th octant
     * \param[in] idx Local index of the target octant.
     * \param[in] balance Has the octant to be balanced?
     */
    void
    LocalTree::setBalance(int32_t idx, bool balance){
        m_octants[idx].setBalance(balance);
    };

    /*! Set the codimension for 2:1 balancing
     * \param[in] b21codim codimension of the entity through which the 2:1 balance is performed.
     */
    void
    LocalTree::setBalanceCodim(uint8_t b21codim){
        m_balanceCodim = b21codim;
    };

    /*!Set the Morton number of first descentant octant of the octree.
     */
    void
    LocalTree::setFirstDescMorton(){
        if(m_sizeOctants){
            octvector::const_iterator firstOctant = m_octants.begin();
            m_firstDescMorton = firstOctant->getMorton();
        } else {
            m_firstDescMorton = std::numeric_limits<uint64_t>::max();
        }
    };

    /*!Set the Morton number of last descentant octant of the octree.
     */
    void
    LocalTree::setLastDescMorton(){
        if(m_sizeOctants){
            octvector::const_iterator lastOctant = m_octants.end() - 1;
            uint32_t x,y,z,delta;
            delta = (uint32_t)(1<<((uint8_t)TreeConstants::MAX_LEVEL - lastOctant->m_level)) - 1;
            x = lastOctant->getLogicalX() + delta;
            y = lastOctant->getLogicalY() + delta;
            z = lastOctant->getLogicalZ() + (m_dim-2)*delta;
            Octant lastDesc = Octant(m_dim, TreeConstants::MAX_LEVEL,x,y,z);
            m_lastDescMorton = lastDesc.getMorton();
        } else {
            m_lastDescMorton = 0;
        }
    };

    /*! Set the periodic condition of the boundaries.
     * \param[in] periodic Vector with the periodic conditions (true/false) of each boundary.
     */
    void
    LocalTree::setPeriodic(bvector & periodic){
        m_periodic = periodic;
    };

    // =================================================================================== //
    // OTHER GET/SET METHODS
    // =================================================================================== //

    /*!Check if the face of the specified octant is on a periodic boundary.
     * \param[in] oct Pointer to the current octant
     * \param[in] iface Index of face
     * \return Returns true if the face of the specified octant is on a
     * periodic boundary, false otherwise.
     */
    bool
    LocalTree::isPeriodic(const Octant* oct, uint8_t iface) const{
        // Check if face is on a boundary
        if (!oct->getBound(iface)) {
            return false;
        }

        // Check if boundary is periodic
        return m_periodic[iface];
    };

    /*!Check if the edge of the specified octant is on a periodic boundary.
     * \param[in] oct Pointer to the current octant
     * \param[in] iedge Index of edge
     * \return Returns true if the edge of the specified octant is on a
     * periodic boundary, false otherwise.
     */
    bool
    LocalTree::isEdgePeriodic(const Octant* oct, uint8_t iedge) const{
        // The edge needs to be on a border
        if (!oct->getEdgeBound(iedge)) {
            return false;
        }

        // If one of the edge faces is on a non-periodic border the edge is
        // non-periodic
        for (int face : m_treeConstants->edgeFace[iedge]) {
            if (oct->getBound(face) && !isPeriodic(oct, face)) {
                return false;
            }
        }

        // The edge is periodic
        return true;
    };

    /*!Check if the node of the specified octant is on a periodic boundary.
     * \param[in] oct Pointer to the current octant
     * \param[in] inode Index of node
     * \return Returns true if the node of the specified octant is on a
     * periodic boundary, false otherwise.
     */
    bool
    LocalTree::isNodePeriodic(const Octant* oct, uint8_t inode) const{
        // The node needs to be on a border
        if (!oct->getNodeBound(inode)) {
            return false;
        }

        // If one of the node faces is on a non-periodic border the node is
        // non-periodic
        for (int face : m_treeConstants->nodeFace[inode]) {
            if (oct->getBound(face) && !isPeriodic(oct, face)) {
                return false;
            }
        }

        // The edge is periodic
        return true;
    };

    // =================================================================================== //
    // OTHER METHODS
    // =================================================================================== //

    /*!Initialize a dummy octree.
     */
    void
    LocalTree::initialize() {
        initialize(0);
    }

    /*!Initialize the octree.
     * \param[in] dim Space dimension of octree.
     */
    void
    LocalTree::initialize(uint8_t dim) {
        m_dim          = dim;
        m_balanceCodim = 1;

        m_periodic.resize(m_dim*2);

        if (m_dim > 0) {
            m_treeConstants = &(TreeConstants::instance(m_dim));
        } else {
            m_treeConstants = nullptr;
        }
    }

    /*!Reset the octree.
     */
    void
    LocalTree::reset(bool createRoot){
        m_octants.clear();
        m_ghosts.clear();
        m_globalIdxGhosts.clear();

        m_lastGhostBros.clear();
        m_firstGhostBros.clear();

        clearConnectivity();
        intervector().swap(m_intersections);

        std::fill(m_periodic.begin(), m_periodic.end(), false);

        if (m_dim > 0 && createRoot) {
            m_localMaxDepth = 0;

            m_octants.push_back(Octant(m_dim));
        } else {
            m_localMaxDepth = -1;
        }

        m_sizeGhosts  = m_ghosts.size();
        m_sizeOctants = m_octants.size();

        setFirstDescMorton();
        setLastDescMorton();

    };

    /*!Extract an octant of the octree.
     * \param[in] idx Local index of the target octant.
     * \return Reference to the idx-th octant of the octree.
     */
    Octant&
    LocalTree::extractOctant(uint32_t idx){
        return m_octants[idx];
    };

    /*!Extract an octant of the octree.
     * \param[in] idx Local index of the target octant.
     * \return Constant reference to the idx-th octant of the octree.
     */
    const Octant&
    LocalTree::extractOctant(uint32_t idx) const{
        return m_octants[idx];
    };

    /*!Extract a ghost octant of the octree.
     * \param[in] idx Local index of the target ghost octant.
     * \return Reference to the idx-th ghost octant of the octree.
     */
    Octant&
    LocalTree::extractGhostOctant(uint32_t idx) {
        return m_ghosts[idx];
    };

    /*!Extract a ghost octant of the octree.
     * \param[in] idx Local index of the target ghost octant.
     * \return Constant reference to the idx-th ghost octant of the octree.
     */
    const Octant&
    LocalTree::extractGhostOctant(uint32_t idx) const{
        return m_ghosts[idx];
    };

    // =================================================================================== //

    /*! Refine local tree: refine one time octants with marker >0
     * \param[out] mapidx mapidx[i] = index in old octants vector of the new i-th octant (index of father if octant is new after refinement)
     * \return	true if refinement done
     */
    bool
    LocalTree::refine(u32vector & mapidx){

        u32vector		last_child_index;
        uint32_t 		idx, ilastch;
        uint32_t 		offset = 0, blockidx;
        uint32_t		mapsize = mapidx.size();
        uint8_t 		nchm1 = m_treeConstants->nChildren-1, ich;
        bool 			dorefine = false;

        m_sizeOctants = m_octants.size();
        for (idx=0; idx<m_sizeOctants; idx++){
            if(m_octants[idx].getMarker() > 0 && m_octants[idx].getLevel() < TreeConstants::MAX_LEVEL){
                last_child_index.push_back(idx+nchm1+offset);
                offset += nchm1;
            }
            else{
                if (m_octants[idx].m_marker > 0){
                    m_octants[idx].m_marker = 0;
                    m_octants[idx].m_info[Octant::INFO_AUX] = false;
                }
            }
        }

        if (offset > 0){
            uint8_t nChildren;
            std::array<Octant, TreeConstants::MAX_CHILDREN> children;

            // Resize octant container
            //
            // We want to be sure the container capacity is equal to its size.
            m_sizeOctants += offset;
            m_octants.reserve(m_sizeOctants);
            m_octants.resize(m_sizeOctants, Octant(m_dim));
            m_octants.shrink_to_fit();

            // Create new octants
            if(mapsize > 0){
                mapidx.resize(m_sizeOctants);
            }
            blockidx = last_child_index[0]-nchm1;
            idx = m_sizeOctants;
            ilastch = last_child_index.size()-1;

            while (idx>blockidx){
                idx--;
                if(idx == last_child_index[ilastch]){
                    m_octants[idx-offset].m_info[Octant::INFO_AUX] = false;
                    // Create children
                    nChildren = m_octants[idx-offset].countChildren();
                    m_octants[idx-offset].buildChildren(children.data());
                    //Update local max depth
                    if (children[0].getLevel() > m_localMaxDepth){
                        m_localMaxDepth = children[0].getLevel();
                    }
                    if (children[0].getMarker()){
                        //More Refinement to do
                        dorefine = true;
                    }
                    // Insert children in the tree
                    for (ich=0; ich<nChildren; ich++){
                        m_octants[idx-ich] = std::move(children[nchm1-ich]);
                        if (mapsize>0) {
                            mapidx[idx-ich] = mapidx[idx-offset];
                        }
                    }
                    offset -= nchm1;
                    idx -= nchm1;
                    if (ilastch != 0){
                        ilastch--;
                    }
                }
                else {
                    m_octants[idx] = m_octants[idx-offset];
                    if(mapsize>0) mapidx[idx]  = mapidx[idx-offset];
                }
            }
        }

        return dorefine;

    };

    // =================================================================================== //
    /*! Coarse local tree: coarse one time family of octants with marker <0
     * (if at least one octant of family has marker>=0 set marker=0 for the entire family)
     * \param[out] mapidx mpaidx[i] = index in old octants vector of the new i-th octant (index of first child if octant is new after coarsening)
     * \return true coarsening can continue (impossible for now forced one refinement for call)
     */
    bool
    LocalTree::coarse(u32vector & mapidx){

        u32vector		first_child_index;
        Octant			father(m_dim);
        uint32_t 		idx, idx2;
        uint32_t 		offset;
        uint32_t 		idx1_gh;
        uint32_t 		idx2_gh;
        uint32_t 		nidx;
        uint32_t		mapsize = mapidx.size();
        int8_t 			markerfather, marker;
        uint8_t 		nbro, nend, nstart;
        uint8_t 		nchm1 = m_treeConstants->nChildren-1;
        bool 			docoarse = false;
        bool 			wstop = false;

        //------------------------------------------ //
        // Initialization

        m_sizeOctants = m_octants.size();
        m_sizeGhosts = m_ghosts.size();

        if (m_sizeOctants == 0) return false;

        nbro = nend = nstart = 0;
        nidx = offset = 0;

        idx2_gh = 0;
        idx1_gh = m_sizeGhosts - 1;

        //------------------------------------------ //

        // Set index for start and end check for ghosts
        if (m_ghosts.size()){
            bool check = true;
            while(check){
                check = idx1_gh < m_sizeGhosts;
                if (check){
                    check = m_ghosts[idx1_gh].getMorton() > m_firstDescMorton;
                }
                if (check) idx1_gh--;
            }

            check = true;
            while(check){
                check = idx2_gh < m_sizeGhosts;
                if (check){
                    check = m_ghosts[idx2_gh].getMorton() < m_lastDescMorton;
                }
                if (check) idx2_gh++;
            }
        }

        // Check and coarse internal octants
        for (idx=0; idx<m_sizeOctants; idx++){
            if(m_octants[idx].getMarker() < 0 && m_octants[idx].getLevel() > 0){
                nbro = 0;
                father = m_octants[idx].buildFather();
                // Check if family is to be refined
                for (idx2=idx; idx2<idx+m_treeConstants->nChildren; idx2++){
                    if (idx2<m_sizeOctants){
                        if(m_octants[idx2].getMarker() < 0 && m_octants[idx2].buildFather() == father){
                            nbro++;
                        }
                    }
                }
                if (nbro == m_treeConstants->nChildren){
                    nidx++;
                    first_child_index.push_back(idx);
                    idx = idx2-1;
                }
            }
        }
        uint32_t nblock = m_sizeOctants;
        uint32_t nfchild = first_child_index.size();
        if (nidx!=0){
            nblock = m_sizeOctants - nidx*nchm1;
            nidx = 0;
            for (idx=0; idx<nblock; idx++){
                if (idx+offset < m_sizeOctants){
                    if (nidx < nfchild){
                        if (idx+offset == first_child_index[nidx]){
                            markerfather = -TreeConstants::MAX_LEVEL;
                            father = m_octants[idx+offset].buildFather();
                            for (uint32_t iii=0; iii<Octant::INFO_ITEM_COUNT; iii++){
                                father.m_info[iii] = false;
                            }
                            father.setGhostLayer(-1);
                            for(idx2=0; idx2<m_treeConstants->nChildren; idx2++){
                                if (idx2 < m_sizeOctants){
                                    if (markerfather < m_octants[idx+offset+idx2].getMarker()+1){
                                        markerfather = m_octants[idx+offset+idx2].getMarker()+1;
                                    }
                                    for (uint32_t iii=0; iii<Octant::INFO_ITEM_COUNT; iii++){
                                        father.m_info[iii] = father.m_info[iii] || m_octants[idx+offset+idx2].m_info[iii];
                                    }
                                }
                            }
                            father.m_info[Octant::INFO_NEW4COARSENING] = true;
                            father.m_info[Octant::INFO_AUX] = false;
                            father.setMarker(markerfather);
                            //Impossible in this version
//                            if (markerfather < 0 && mapsize == 0){
//                                docoarse = true;
//                            }
                            m_octants[idx] = father;
                            if(mapsize > 0) mapidx[idx] = mapidx[idx+offset];
                            offset += nchm1;
                            nidx++;
                        }
                        else{
                            m_octants[idx] = m_octants[idx+offset];
                            if(mapsize > 0) mapidx[idx] = mapidx[idx+offset];
                        }
                    }
                    else{
                        m_octants[idx] = m_octants[idx+offset];
                        if(mapsize > 0) mapidx[idx] = mapidx[idx+offset];
                    }
                }
            }
        }
        m_octants.resize(nblock, Octant(m_dim));
        m_octants.shrink_to_fit();
        m_sizeOctants = m_octants.size();
        if(mapsize > 0){
            mapidx.resize(m_sizeOctants);
        }

        //Check ghosts
        if (m_ghosts.size()){
            // Start on ghosts
            if (m_sizeOctants > 0 && idx1_gh < m_sizeGhosts){
                if (m_ghosts[idx1_gh].buildFather() == m_octants[0].buildFather()){
                    father = m_ghosts[idx1_gh].buildFather();
                    nbro = 0;
                    idx = idx1_gh;
                    marker = m_ghosts[idx].getMarker();
                    bool waszero = (idx == 0);
                    while(marker < 0 && m_ghosts[idx].buildFather() == father){
                        nbro++;
                        if(waszero){
                            break;
                        }
                        idx--;
                        if(idx == 0){
                            waszero = true;
                        }
                        marker = m_ghosts[idx].getMarker();
                    }
                    nstart = 0;
                    idx = 0;
                    marker = m_octants[idx].getMarker();
                    if (idx==m_sizeOctants-1) wstop = true;
                    while(marker < 0 && m_octants[idx].buildFather() == father){
                        nbro++;
                        nstart++;
                        if (wstop){
                            break;
                        }
                        idx++;
                        if (idx==m_sizeOctants-1){
                            wstop = true;
                        }
                        marker = m_octants[idx].getMarker();
                    }
                    if (nbro == m_treeConstants->nChildren){
                        offset = nstart;
                    }
                    else{
                        nstart = 0;
                    }
                }
                if (nstart != 0){
                    for (idx=0; idx<m_sizeOctants; idx++){
                        if (idx+offset < m_sizeOctants){
                            m_octants[idx] = m_octants[idx+offset];
                            if(mapsize > 0) mapidx[idx] = mapidx[idx+offset];
                        }
                    }
                    m_octants.resize(m_sizeOctants-offset, Octant(m_dim));
                    m_octants.shrink_to_fit();
                    m_sizeOctants = m_octants.size();
                    if(mapsize > 0){
                        mapidx.resize(m_sizeOctants);
                    }
                }
            }


            //Verify family between more then two processes
            if (m_sizeOctants > 0 && idx2_gh < m_sizeGhosts){

                if (m_ghosts[idx2_gh].buildFather() == father){

                    if (m_ghosts[idx2_gh].buildFather() == m_octants[m_sizeOctants-1].buildFather()){

                        uint64_t idx22_gh = idx2_gh;
                        marker = m_ghosts[idx22_gh].getMarker();
                        while(marker < 0 && m_ghosts[idx22_gh].buildFather() == father){
                            nbro++;
                            idx22_gh++;
                            if(idx22_gh == m_sizeGhosts){
                                break;
                            }
                            marker = m_ghosts[idx22_gh].getMarker();
                        }

                        if (nbro == m_treeConstants->nChildren){
                            m_octants.clear();
                            m_sizeOctants = 0;
                            if(mapsize > 0){
                                mapidx.clear();
                            }
                        }
                    }

                }

            }


            // End on ghosts
            if (m_sizeOctants > 0 && idx2_gh < m_sizeGhosts){
                if (m_ghosts[idx2_gh].buildFather() == m_octants[m_sizeOctants-1].buildFather()){
                    father = m_ghosts[idx2_gh].buildFather();
                    for (uint32_t iii=0; iii<Octant::INFO_ITEM_COUNT; iii++){
                        father.m_info[iii] = false;
                    }
                    father.setGhostLayer(-1);
                    markerfather = m_ghosts[idx2_gh].getMarker()+1;
                    nbro = 0;
                    idx = idx2_gh;
                    marker = m_ghosts[idx].getMarker();
                    while(marker < 0 && m_ghosts[idx].buildFather() == father){
                        nbro++;
                        if (markerfather < m_ghosts[idx].getMarker()+1){
                            markerfather = m_ghosts[idx].getMarker()+1;
                        }
                        for (uint32_t iii=0; iii<m_treeConstants->nFaces; iii++){
                            father.m_info[iii] = father.m_info[iii] || m_ghosts[idx].m_info[iii];
                        }
                        father.m_info[Octant::INFO_BALANCED] = father.m_info[Octant::INFO_BALANCED] || m_ghosts[idx].m_info[Octant::INFO_BALANCED];
                        idx++;
                        if(idx == m_sizeGhosts){
                            break;
                        }
                        marker = m_ghosts[idx].getMarker();
                    }
                    nend = 0;
                    idx = m_sizeOctants-1;
                    marker = m_octants[idx].getMarker();
                    if (idx==0) wstop = true;
                    while(marker < 0 && m_octants[idx].buildFather() == father){
                        nbro++;
                        nend++;
                        if (markerfather < m_octants[idx].getMarker()+1){
                            markerfather = m_octants[idx].getMarker()+1;
                        }
                        idx--;
                        if (wstop){
                            break;
                        }
                        if (idx==0){
                            wstop = true;
                        }
                        marker = m_octants[idx].getMarker();
                    }
                    if (nbro == m_treeConstants->nChildren){
                        offset = nend;
                    }
                    else{
                        nend = 0;
                    }
                }
                if (nend != 0){
                    for (idx=0; idx < nend; idx++){
                        for (uint32_t iii=0; iii<Octant::INFO_ITEM_COUNT - 1; iii++){
                            father.m_info[iii] = father.m_info[iii] || m_octants[m_sizeOctants-idx-1].m_info[iii];
                        }
                    }
                    father.m_info[Octant::INFO_NEW4COARSENING] = true;
                    father.m_info[Octant::INFO_AUX] = false;
                    father.setGhostLayer(-1);
                    //Impossible in this version
                    //                if (markerfather < 0 && mapsize == 0){
                        //                    docoarse = true;
                        //                }
                    father.setMarker(markerfather);
                    m_octants.resize(m_sizeOctants-offset, Octant(m_dim));
                    m_octants.push_back(father);
                    m_octants.shrink_to_fit();
                    m_sizeOctants = m_octants.size();
                    if(mapsize > 0){
                        mapidx.resize(m_sizeOctants);
                    }
                }
            }

        }//end if ghosts size

        m_sizeOctants = m_octants.size();

        // Update maximum depth
        updateLocalMaxDepth();

        // Set final last desc
        setFirstDescMorton();
        setLastDescMorton();

        return docoarse;

    };

    // =================================================================================== //

    /*! Refine local tree: refine one time all the octants
     * \param[out] mapidx mpaidx[i] = index in old octants vector of the new i-th octant (index of father if octant is new after refinement)
     * \return	true if refinement done
     */
    bool
    LocalTree::globalRefine(u32vector & mapidx){

        uint32_t 	idx;
        bool 		dorefine = false;

        for (idx=0; idx<m_sizeOctants; idx++){
            m_octants[idx].setMarker(1);
        }

        dorefine = refine(mapidx);

        return dorefine;

    };

    // =================================================================================== //

    /*! Refine local tree: corse one time all the octants
     * \param[out] mapidx mpaidx[i] = index in old octants vector of the new i-th octant (index of first child if octant is new after coarsening)
     * \return	true if coarsening can continue (impossible in this version)
     */
    bool
    LocalTree::globalCoarse(u32vector & mapidx){

        uint32_t 	idx;
        bool 		dorefine = false;

        for (idx=0; idx<m_sizeOctants; idx++){
            m_octants[idx].setMarker(-1);
        }
        for (idx=0; idx<m_sizeGhosts; idx++){
            m_ghosts[idx].setMarker(-1);
        }

        dorefine = coarse(mapidx);

        return dorefine;

    };

    /*! Delete overlapping octants after coarse local tree. Check if first octants of the partition
     * have marker = -1 (after coarse only the octants to be deleted have marker =-1).
     * \param[in] partLastDesc Last descendant of process before the actual
     * \param[out] mapidx mapidx[i] = index in old octants vector of the new i-th octant (index of first child if octant is new after coarsening)
     */
    void
    LocalTree::checkCoarse(uint64_t partLastDesc, u32vector & mapidx){

        uint32_t        idx;
        uint32_t        mapsize = mapidx.size();
        uint8_t         toDelete = 0;

        if (m_sizeOctants>0){

            idx = 0;
            if (m_octants[idx].getMorton() < partLastDesc){

                Octant father0 = m_octants[idx].buildFather();
                Octant father = father0;

                while(father == father0 && idx < m_sizeOctants){
                    toDelete++;
                    idx++;
                    if (idx<m_sizeOctants) father = m_octants[idx].buildFather();
                }

                if (m_sizeOctants>toDelete){
                    for(idx=0; idx<m_sizeOctants-toDelete; idx++){
                        m_octants[idx] = m_octants[idx+toDelete];
                        if (mapsize>0) mapidx[idx] = mapidx[idx+toDelete];
                    }
                    m_octants.resize(m_sizeOctants-toDelete, Octant(m_dim));
                    m_sizeOctants = m_octants.size();
                    if (mapsize>0){
                        mapidx.resize(m_sizeOctants);
                    }
                }
                else{
                    m_octants.clear();
                    mapidx.clear();
                }
                m_sizeOctants = m_octants.size();
                setFirstDescMorton();
            }

        }
    };

    // =================================================================================== //
    /*! Update max depth reached in local tree
     */
    void
    LocalTree::updateLocalMaxDepth(){

        uint32_t noctants = getNumOctants();
        uint32_t i;

        m_localMaxDepth = 0;
        for(i = 0; i < noctants; i++){
            if(m_octants[i].getLevel() > m_localMaxDepth){
                m_localMaxDepth = m_octants[i].getLevel();
            }
        }
    };


    /** Finds local and ghost or only local neighbours of octant(both local and ghost ones) through iface face.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] oct Pointer to the current octant
     * \param[in] iface Index of face passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     * \param[in] onlyinternal A boolean flag to specify if neighbours have to be found among all the octants (false) or only among the internal ones (true).
     */
    void
    LocalTree::findNeighbours(const Octant* oct, uint8_t iface, u32vector & neighbours, bvector & isghost, bool onlyinternal) const{

        isghost.clear();
        neighbours.clear();

        // Default if iface is nface<iface<0
        if (iface >= m_treeConstants->nFaces){
            return;
        }

        // If a face is a boundary, it can have neighbours only if periodic
        bool isperiodic = false;
        if (oct->getBound(iface)) {
            isperiodic = isPeriodic(oct, iface);
            if (!isperiodic) {
                return;
            }
        }

        // Exit if is the root octant
        uint8_t level = oct->getLevel();
        if (level == 0){
            // if periodic face return itself
            if (isperiodic){
            	neighbours.push_back(0);
            	isghost.push_back(false);
            }
            return;
        }

        // Initialize search in the octants
        uint32_t candidateIdx    = 0;
        uint64_t candidateMorton = 0;

        uint64_t neighArea = 0;
        uint64_t faceArea  = oct->getLogicalArea();

        uint32_t size = oct->getLogicalSize();
        const int8_t (&cxyz)[3] = m_treeConstants->normals[iface];

        // Build Morton number of virtual neigh of same size
        Octant sameSizeVirtualNeigh;
        if (isperiodic){
            sameSizeVirtualNeigh = oct->computePeriodicOctant(iface);
        }
        else{
            u32array3 sameSizeVirtualNeighCoords = oct->getLogicalCoordinates();
            sameSizeVirtualNeighCoords[0] += cxyz[0] * size;
            sameSizeVirtualNeighCoords[1] += cxyz[1] * size;
            sameSizeVirtualNeighCoords[2] += cxyz[2] * size;

            sameSizeVirtualNeigh = Octant(m_dim, level, sameSizeVirtualNeighCoords[0], sameSizeVirtualNeighCoords[1], sameSizeVirtualNeighCoords[2]);
        }

        uint64_t sameSizeVirtualNeighMorton = sameSizeVirtualNeigh.getMorton();

        //
        // Search in the internal octants
        //

        // Identify the index of the first neighbour candidate
        computeNeighSearchBegin(sameSizeVirtualNeighMorton, m_octants, &candidateIdx, &candidateMorton);

        // Early return if a neighbour of the same size has been found
        if(candidateMorton == sameSizeVirtualNeighMorton && m_octants[candidateIdx].m_level == level){
            neighbours.push_back(candidateIdx);
            isghost.push_back(false);
            return;
        }

        // Compute the Morton number of the last candidate
        //
        // This is the Morton number of the last discendent of the same-size
        // virtual neighbour.
        uint64_t lastCandidateMorton = sameSizeVirtualNeigh.buildLastDesc().getMorton();

        // Compute coordinates
        std::array<int64_t,3> coord;
        if (isperiodic){
            coord = oct->getPeriodicCoord(iface);
        }
        else{
            coord[0] = oct->getLogicalX();
            coord[1] = oct->getLogicalY();
            coord[2] = oct->getLogicalZ();
        }

        // Search for neighbours of different sizes
        if (candidateIdx < m_sizeOctants){
            while(true){
                // Detect if the candidate is a neighbour
                u32array3 coordtry = m_octants[candidateIdx].getLogicalCoordinates();

                bool isNeighbourCandidate = true;
                for (int idim=0; idim<m_dim; idim++){
                    int32_t Dx     = int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
                    int32_t Dxstar = int32_t((cxyz[idim]-1)/2)*(m_octants[candidateIdx].getLogicalSize()) + int32_t((cxyz[idim]+1)/2)*size;
                    if (Dx != Dxstar){
                        isNeighbourCandidate = false;
                        break;
                    }
                }

                if (isNeighbourCandidate){
                    bool isNeighbour = false;
                    uint8_t leveltry = m_octants[candidateIdx].getLevel();
                    if (leveltry > level){
                        array<int64_t,3> coord1 ={{1, 1, 1}} ;
                        for (int idim=0; idim<m_dim; idim++){
                            coord1[idim] = coord[idim] + size;
                        }

                        if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
                            isNeighbour = true;
                        }
                    }
                    else if (leveltry < level){
                        u32array3 coordtry1 = {{1, 1, 1}};
                        for (int idim=0; idim<m_dim; idim++){
                            coordtry1[idim] = coordtry[idim] + m_octants[candidateIdx].getLogicalSize();
                        }

                        if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
                            isNeighbour = true;
                        }
                    }

                    if (isNeighbour){
                        neighbours.push_back(candidateIdx);
                        isghost.push_back(false);

                        // If the neighbours already cover the whole face, we have
                        // found all the neighbours and we can exit.
                        neighArea += m_octants[candidateIdx].getLogicalArea();
                        if (neighArea == faceArea){
                            return;
                        }
                    }
                }

                // Advance to the next candidate
                candidateIdx++;
                if (candidateIdx > m_sizeOctants - 1){
                    break;
                }

                candidateMorton = m_octants[candidateIdx].getMorton();
                if (candidateMorton > lastCandidateMorton){
                    break;
                }
            }
        }

        //
        // Search in the ghost octants
        //
        // If we want also the neighbours that are ghosts, we always need to
        // search in the ghosts, the only exception is for faces of internal
        // octants that are not process boundaries.
        bool ghostSearch = !onlyinternal && (m_sizeGhosts > 0);
        if (ghostSearch){
            if (!oct->getIsGhost() && !oct->getPbound(iface)){
                ghostSearch = false;
            }
        }

        // Search in ghosts
        if(ghostSearch){
            // Identify the index of the first neighbour candidate
            computeNeighSearchBegin(sameSizeVirtualNeighMorton, m_ghosts, &candidateIdx, &candidateMorton);

            // Early return if a neighbour of the same size has been found
            if(candidateMorton == sameSizeVirtualNeighMorton && m_ghosts[candidateIdx].getLevel() == level){
                neighbours.push_back(candidateIdx);
                isghost.push_back(true);
                return;
            }

            // Search for neighbours of different sizes
            if (candidateIdx < m_sizeGhosts){
                while(true){
                    // Detect if the candidate is a neighbour
                    u32array3 coordtry = m_ghosts[candidateIdx].getLogicalCoordinates();

                    bool isNeighbourCandidate = true;
                    for (int idim=0; idim<m_dim; idim++){
                        int32_t Dx     = int32_t(int32_t(abs(cxyz[idim]))*(-coord[idim] + coordtry[idim]));
                        int32_t Dxstar = int32_t((cxyz[idim]-1)/2)*(m_ghosts[candidateIdx].getLogicalSize()) + int32_t((cxyz[idim]+1)/2)*size;
                        if (Dx != Dxstar){
                            isNeighbourCandidate = false;
                            break;
                        }
                    }

                    if (isNeighbourCandidate){
                        bool isNeighbour = false;
                        uint8_t leveltry = m_ghosts[candidateIdx].getLevel();
                        if (leveltry > level){
                            array<int64_t, 3> coord1 = {{1, 1, 1}};
                            for (int idim=0; idim<m_dim; idim++){
                                coord1[idim] = coord[idim] + size;
                            }

                            if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
                                isNeighbour = true;
                            }
                        }
                        else if (leveltry < level){
                            u32array3 coordtry1 = {{1, 1, 1}};
                            for (int idim=0; idim<m_dim; idim++){
                                coordtry1[idim] = coordtry[idim] + m_ghosts[candidateIdx].getLogicalSize();
                            }

                            if((abs(cxyz[0])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[1])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2]))) + (abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1])))){
                                isNeighbour = true;
                            }
                        }

                        if (isNeighbour){
                            neighbours.push_back(candidateIdx);
                            isghost.push_back(true);

                            // If the neighbours already cover the whole face, we have
                            // found all the neighbours and we can exit.
                            neighArea += m_ghosts[candidateIdx].getLogicalArea();
                            if (neighArea == faceArea){
                                return;
                            }
                        }
                    }

                    candidateIdx++;
                    if (candidateIdx > m_sizeGhosts - 1){
                        break;
                    }

                    candidateMorton = m_ghosts[candidateIdx].getMorton();
                    if (candidateMorton > lastCandidateMorton){
                        break;
                    }
                }
            }
        }
    };

    /** Finds local and ghost or only local neighbours of octant(both local and ghost ones) through iedge edge.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] oct Pointer to the current octant
     * \param[in] iedge Index of edge passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     * \param[in] onlyinternal A boolean flag to specify if neighbours have to be found among all the octants (false) or only among the internal ones (true).
     */
    void
    LocalTree::findEdgeNeighbours(const Octant* oct, uint8_t iedge, u32vector & neighbours, bvector & isghost, bool onlyinternal) const{

        isghost.clear();
        neighbours.clear();

        // Default if iedge is nface<iedge<0
        if (iedge >= m_treeConstants->nEdges){
            return;
        }

        // If a edge is a on boundary, it can have neighbours only if periodic
        bool isperiodic = false;
        if (oct->getEdgeBound(iedge)) {
            isperiodic = isEdgePeriodic(oct, iedge);
            if (!isperiodic) {
                return;
            }
        }

        // Exit if is the root octant
        uint8_t level = oct->getLevel();
        if (level == 0){
            // if periodic face return itself
            if (isperiodic){
                neighbours.push_back(0);
                isghost.push_back(false);
            }
            return;
        }

        // Search in the octants
        uint32_t candidateIdx    = 0;
        uint64_t candidateMorton = 0;

        uint32_t neighSize = 0;
        uint32_t edgeSize  = oct->getLogicalSize();

        uint32_t size = oct->getLogicalSize();
        const int8_t (&cxyz)[3] = m_treeConstants->edgeCoeffs[iedge];

        // Build Morton number of virtual neigh of same size
        Octant sameSizeVirtualNeigh;
        if (isperiodic){
            sameSizeVirtualNeigh = oct->computeEdgePeriodicOctant(iedge);
        }
        else{
            u32array3 sameSizeVirtualNeighCoords = oct->getLogicalCoordinates();
            sameSizeVirtualNeighCoords[0] += cxyz[0] * size;
            sameSizeVirtualNeighCoords[1] += cxyz[1] * size;
            sameSizeVirtualNeighCoords[2] += cxyz[2] * size;

            sameSizeVirtualNeigh = Octant(m_dim, level, sameSizeVirtualNeighCoords[0], sameSizeVirtualNeighCoords[1], sameSizeVirtualNeighCoords[2]);
        }

        uint64_t sameSizeVirtualNeighMorton = sameSizeVirtualNeigh.getMorton();

        //
        // Search in the internal octants
        //

        // Identify the index of the first neighbour candidate
        computeNeighSearchBegin(sameSizeVirtualNeighMorton, m_octants, &candidateIdx, &candidateMorton);

        // Early return if a neighbour of the same size has been found
        if(candidateMorton == sameSizeVirtualNeighMorton && m_octants[candidateIdx].m_level == level){
            isghost.push_back(false);
            neighbours.push_back(candidateIdx);
            return;
        }

        // Compute the Morton number of the last candidate
        //
        // This is the Morton number of the last discendent of the same-size
        // virtual neighbour.
        uint64_t lastCandidateMorton = sameSizeVirtualNeigh.buildLastDesc().getMorton();

        // Compute coordinates
        std::array<int64_t,3> coord;
        if (isperiodic){
            coord = oct->getEdgePeriodicCoord(iedge);
        }
        else{
            coord[0] = oct->getLogicalX();
            coord[1] = oct->getLogicalY();
            coord[2] = oct->getLogicalZ();
        }

        // Search for neighbours of different sizes
        if (candidateIdx < m_sizeOctants) {
            while(true){
                // Detect if the candidate is a neighbour
                u32array3 coordtry = m_octants[candidateIdx].getLogicalCoordinates();

                bool isNeighbourCandidate = true;
                for (int idim=0; idim<m_dim; idim++){
                    int32_t Dx     = int32_t(int32_t(abs(cxyz[idim]))*(-int32_t(coord[idim]) + int32_t(coordtry[idim])));
                    int32_t Dxstar = int32_t((cxyz[idim]-1)/2)*(m_octants[candidateIdx].getLogicalSize()) + int32_t((cxyz[idim]+1)/2)*size;

                    if (Dx != Dxstar) {
                        isNeighbourCandidate = false;
                        break;
                    }
                }

                if (isNeighbourCandidate){
                    bool isNeighbour = false;
                    uint8_t leveltry = m_octants[candidateIdx].getLevel();
                    if (leveltry > level){
                        u32array3 coord1 = {{1, 1, 1}};
                        for (int idim=0; idim<m_dim; idim++){
                            coord1[idim] = coord[idim] + size;
                        }

                        if((abs(cxyz[0])*abs(cxyz[2])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))) + (abs(cxyz[1])*abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))) + (abs(cxyz[0])*abs(cxyz[1])*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2])))){
                            isNeighbour = true;
                        }
                    }
                    else if (leveltry < level){
                        u32array3 coordtry1 = {{1, 1, 1}};
                        for (int idim=0; idim<m_dim; idim++){
                            coordtry1[idim] = coordtry[idim] + m_octants[candidateIdx].getLogicalSize();
                        }

                        if((abs(cxyz[0])*abs(cxyz[2])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))) + (abs(cxyz[1])*abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))) + (abs(cxyz[0])*abs(cxyz[1])*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2])))){
                            isNeighbour = true;
                        }
                    }

                    if (isNeighbour) {
                        neighbours.push_back(candidateIdx);
                        isghost.push_back(false);

                        // If the neighbours already cover the whole edge, we have
                        // found all the neighbours and we can exit.
                        neighSize += m_octants[candidateIdx].getLogicalSize();
                        if (neighSize == edgeSize){
                            return;
                        }
                    }
                }

                // Advance to the next candidate
                candidateIdx++;
                if (candidateIdx > m_sizeOctants-1){
                    break;
                }

                candidateMorton = m_octants[candidateIdx].getMorton();
                if (candidateMorton > lastCandidateMorton){
                    break;
                }
            }
        }

        //
        // Search in the ghost octants
        //
        if (m_sizeGhosts > 0 && !onlyinternal){
            // Identify the index of the first neighbour candidate
            computeNeighSearchBegin(sameSizeVirtualNeighMorton, m_ghosts, &candidateIdx, &candidateMorton);

            // Early return if a neighbour of the same size has been found
            if(candidateMorton == sameSizeVirtualNeighMorton && m_ghosts[candidateIdx].m_level == level){
                isghost.push_back(true);
                neighbours.push_back(candidateIdx);
                return;
            }

            // Search for neighbours of different sizes
            if (candidateIdx < m_sizeGhosts){
                while(true){
                    // Detect if the candidate is a neighbour
                    u32array3 coordtry = m_ghosts[candidateIdx].getLogicalCoordinates();

                    bool isNeighbourCandidate = true;
                    for (int idim=0; idim<m_dim; idim++){
                        int32_t Dx     = int32_t(int32_t(abs(cxyz[idim]))*(-int32_t(coord[idim]) + int32_t(coordtry[idim])));
                        int32_t Dxstar = int32_t((cxyz[idim]-1)/2)*(m_ghosts[candidateIdx].getLogicalSize()) + int32_t((cxyz[idim]+1)/2)*size;

                        if (Dx != Dxstar) {
                            isNeighbourCandidate = false;
                            break;
                        }
                    }

                    if (isNeighbourCandidate){
                        bool isNeighbour = false;
                        uint8_t leveltry = m_ghosts[candidateIdx].getLevel();
                        if (leveltry > level){
                            u32array3 coord1 = {{1, 1, 1}};
                            for (int idim=0; idim<m_dim; idim++){
                                coord1[idim] = coord[idim] + size;
                            }

                            if((abs(cxyz[0])*abs(cxyz[2])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))) + (abs(cxyz[1])*abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))) + (abs(cxyz[0])*abs(cxyz[1])*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2])))){
                                isNeighbour = true;
                            }
                        }
                        else if (leveltry < level){
                            u32array3 coordtry1 = {{1, 1, 1}};
                            for (int idim=0; idim<m_dim; idim++){
                                coordtry1[idim] = coordtry[idim] + m_ghosts[candidateIdx].getLogicalSize();
                            }

                            if((abs(cxyz[0])*abs(cxyz[2])*((coord[1]>=coordtry[1])*(coord[1]<coordtry1[1]))) + (abs(cxyz[1])*abs(cxyz[2])*((coord[0]>=coordtry[0])*(coord[0]<coordtry1[0]))) + (abs(cxyz[0])*abs(cxyz[1])*((coord[2]>=coordtry[2])*(coord[2]<coordtry1[2])))){
                                isNeighbour = true;
                            }
                        }

                        if (isNeighbour) {
                            neighbours.push_back(candidateIdx);
                            isghost.push_back(true);

                            // If the neighbours already cover the whole edge, we have
                            // found all the neighbours and we can exit.
                            neighSize += m_ghosts[candidateIdx].getLogicalSize();
                            if (neighSize == edgeSize){
                                return;
                            }
                        }
                    }

                    // Advance to the next candidate
                    candidateIdx++;
                    if (candidateIdx > m_sizeGhosts-1){
                        break;
                    }

                    candidateMorton = m_ghosts[candidateIdx].getMorton();
                    if (candidateMorton > lastCandidateMorton){
                        break;
                    }
                }
            }
        }
    };

    /** Finds local and ghost or only local neighbours of octant(both local and ghost ones) through inode node.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] oct Pointer to the current octant
     * \param[in] inode Index of node passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     * \param[in] onlyinternal A boolean flag to specify if neighbours have to be found among all the octants (false) or only among the internal ones (true).
     */
    void
    LocalTree::findNodeNeighbours(const Octant* oct, uint8_t inode, u32vector & neighbours, bvector & isghost, bool onlyinternal) const{

        isghost.clear();
        neighbours.clear();

        // Default if inode is nnodes<inode<0
        if (inode >= m_treeConstants->nNodes){
            return;
        }

        // If a node is a on boundary, it can have neighbours only if periodic
        bool isperiodic = false;
        if (oct->getNodeBound(inode)) {
            isperiodic = isNodePeriodic(oct, inode);
            if (!isperiodic) {
                return;
            }
        }

        // Exit if is the root octant
        uint8_t level = oct->getLevel();
        if (level == 0){
            // if periodic face return itself
            if (isperiodic){
                neighbours.push_back(0);
                isghost.push_back(false);
            }
            return;
        }

        // Search in the octants
        uint32_t candidateIdx    = 0;
        uint64_t candidateMorton = 0;

        uint32_t size = oct->getLogicalSize();
        const int8_t (&cxyz)[3] = m_treeConstants->nodeCoeffs[inode];

        // Build Morton number of virtual neigh of same size
        Octant sameSizeVirtualNeigh;
        if (isperiodic){
            sameSizeVirtualNeigh = oct->computeNodePeriodicOctant(inode);
        }
        else{
            u32array3 sameSizeVirtualNeighCoords = oct->getLogicalCoordinates();
            sameSizeVirtualNeighCoords[0] += cxyz[0] * size;
            sameSizeVirtualNeighCoords[1] += cxyz[1] * size;
            sameSizeVirtualNeighCoords[2] += cxyz[2] * size;

            sameSizeVirtualNeigh = Octant(m_dim, level, sameSizeVirtualNeighCoords[0], sameSizeVirtualNeighCoords[1], sameSizeVirtualNeighCoords[2]);
        }

        uint64_t sameSizeVirtualNeighMorton = sameSizeVirtualNeigh.getMorton();

        //
        // Search in the internal octants
        //

        // Identify the index of the first neighbour candidate
        computeNeighSearchBegin(sameSizeVirtualNeighMorton, m_octants, &candidateIdx, &candidateMorton);

        // Early return if a neighbour of the same size has been found
        if(candidateMorton == sameSizeVirtualNeighMorton && m_octants[candidateIdx].m_level == oct->m_level){
            isghost.push_back(false);
            neighbours.push_back(candidateIdx);
            return;
        }

        // Compute the Morton number of the last candidate
        //
        // This is the Morton number of the last discendent of the same-size
        // virtual neighbour.
        uint64_t lastCandidateMorton = sameSizeVirtualNeigh.buildLastDesc().getMorton();

        // Compute coordinates
        std::array<int64_t,3> coord;
        if (isperiodic){
            coord = oct->getNodePeriodicCoord(inode);
        }
        else{
            coord[0] = oct->getLogicalX();
            coord[1] = oct->getLogicalY();
            coord[2] = oct->getLogicalZ();
        }

        // Search for neighbours of different sizes
        if (candidateIdx < m_sizeOctants) {
            while(true){
                // Detect if the candidate is a neighbour
                u32array3 coordtry = m_octants[candidateIdx].getLogicalCoordinates();

                bool isNeighbour = true;
                for (int idim=0; idim<m_dim; idim++){
                    int32_t Dx     = int32_t(int32_t(abs(cxyz[idim]))*(-int32_t(coord[idim]) + int32_t(coordtry[idim])));
                    int32_t Dxstar = int32_t((cxyz[idim]-1)/2)*(m_octants[candidateIdx].getLogicalSize()) + int32_t((cxyz[idim]+1)/2)*size;

                    if (Dx != Dxstar) {
                        isNeighbour = false;
                        break;
                    }
                }

                // Advance to the next candidate
                if (isNeighbour){
                    neighbours.push_back(candidateIdx);
                    isghost.push_back(false);
                    return;
                }

                candidateIdx++;
                if (candidateIdx > m_sizeOctants-1){
                    break;
                }

                candidateMorton = m_octants[candidateIdx].getMorton();
                if (candidateMorton > lastCandidateMorton){
                    break;
                }
            }
        }

        //
        // Search in the ghost octants
        //

        if (m_sizeGhosts > 0 && !onlyinternal){
            // Identify the index of the first neighbour candidate
            computeNeighSearchBegin(sameSizeVirtualNeighMorton, m_ghosts, &candidateIdx, &candidateMorton);

            // Early return if a neighbour of the same size has been found
            if(candidateMorton == sameSizeVirtualNeighMorton && m_ghosts[candidateIdx].m_level == oct->m_level){
                isghost.push_back(true);
                neighbours.push_back(candidateIdx);
                return;
            }

            // Search for neighbours of different sizes
            if (candidateIdx < m_sizeGhosts) {
                while(true){
                    // Detect if the candidate is a neighbour
                    u32array3 coordtry = m_ghosts[candidateIdx].getLogicalCoordinates();

                    bool isNeighbour = true;
                    for (int idim=0; idim<m_dim; idim++){
                        int32_t Dx     = int32_t(int32_t(abs(cxyz[idim]))*(-int32_t(coord[idim]) + int32_t(coordtry[idim])));
                        int32_t Dxstar = int32_t((cxyz[idim]-1)/2)*(m_ghosts[candidateIdx].getLogicalSize()) + int32_t((cxyz[idim]+1)/2)*size;

                        if (Dx != Dxstar) {
                            isNeighbour = false;
                            break;
                        }
                    }

                    // Advance to the next candidate
                    if (isNeighbour){
                        neighbours.push_back(candidateIdx);
                        isghost.push_back(true);
                        return;
                    }

                    candidateIdx++;
                    if (candidateIdx > m_sizeGhosts-1){
                        break;
                    }

                    candidateMorton = m_ghosts[candidateIdx].getMorton();
                    if (candidateMorton > lastCandidateMorton){
                        break;
                    }
                }
            }
        }
    };

    /** Finds local and ghost or only local neighbours of octant(both local and ghost ones) through iface face.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] idx Index of current octant
     * \param[in] amIghost Boolean flag to specify if the octant is a ghost one
     * \param[in] iface Index of face passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     * \param[in] onlyinternal A boolean flag to specify if neighbours have to be found among all the octants (false) or only among the internal ones (true).
     */
    void
    LocalTree::findNeighbours(uint32_t idx, bool amIghost,uint8_t iface, u32vector & neighbours, bvector & isghost, bool onlyinternal) const{
        const Octant* oct = amIghost ? &m_ghosts[idx] : &m_octants[idx];
        findNeighbours(oct, iface, neighbours, isghost, onlyinternal);
    };

    /** Finds local and ghost or only local neighbours of octant(both local and ghost ones) through iedge edge.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] idx Index of current octant
     * \param[in] amIghost Boolean flag to specify if the octant is a ghost one
     * \param[in] iedge Index of edge passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     * \param[in] onlyinternal A boolean flag to specify if neighbours have to be found among all the octants (false) or only among the internal ones (true).
     */
    void
    LocalTree::findEdgeNeighbours(uint32_t idx, bool amIghost,uint8_t iedge, u32vector & neighbours, bvector & isghost, bool onlyinternal) const{
        const Octant* oct = amIghost ? &m_ghosts[idx] : &m_octants[idx];
        findEdgeNeighbours(oct, iedge, neighbours, isghost, onlyinternal);
    };

    /** Finds local and ghost or only local neighbours of octant(both local and ghost ones) through inode node.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] idx Index of current octant
     * \param[in] amIghost Boolean flag to specify if the octant is a ghost one
     * \param[in] inode Index of node passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     * \param[in] onlyinternal A boolean flag to specify if neighbours have to be found among all the octants (false) or only among the internal ones (true).
     */
    void
    LocalTree::findNodeNeighbours(uint32_t idx, bool amIghost,uint8_t inode, u32vector & neighbours, bvector & isghost, bool onlyinternal) const{
        const Octant* oct = amIghost ? &m_ghosts[idx] : &m_octants[idx];
        findNodeNeighbours(oct, inode,neighbours, isghost,onlyinternal);
    };

    /** Finds local and ghost neighbours of ghost octant through iface face.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] idx Index of the current ghost octant
     * \param[in] iface Index of face passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     */
	void
    LocalTree::findGhostNeighbours(uint32_t idx, uint8_t iface, u32vector & neighbours, bvector & isghost) const{
        const Octant* oct = &m_ghosts[idx];
        findNeighbours(oct,iface,neighbours, isghost,false);
    };

    /** Finds local and ghost neighbours of ghost octant through iedge edge.
     * Returns a vector (empty if iedge is a bound edge) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] idx Index of the current octant
     * \param[in] iedge Index of edge passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     */
	void
    LocalTree::findGhostEdgeNeighbours(uint32_t idx, uint8_t iedge, u32vector & neighbours, bvector & isghost) const{
        const Octant* oct = &m_ghosts[idx];
        findEdgeNeighbours(oct, iedge, neighbours, isghost, false);
    };

    /** Finds local and ghost neighbours of ghost octant through inode node.
     * Returns a vector (empty if inode is a bound node) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] idx Index of the current octant
     * \param[in] inode Index of node passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     */
	void
    LocalTree::findGhostNodeNeighbours(uint32_t idx, uint8_t inode, u32vector & neighbours, bvector & isghost) const{
        const Octant* oct = &m_ghosts[idx];
        findNodeNeighbours(oct, inode, neighbours, isghost, false);
    };


    /** Finds local and ghost neighbours of local octant through iface face.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] idx Index of the current octant
     * \param[in] iface Index of face passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     */
    void
    LocalTree::findNeighbours(uint32_t idx, uint8_t iface, u32vector & neighbours, bvector & isghost) const{
        findNeighbours(idx, false, iface, neighbours, isghost, false);
    };

    /** Finds local and ghost neighbours of local octant through iedge edge.
     * Returns a vector (empty if iedge is a bound edge) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] idx Index of the current octant
     * \param[in] iedge Index of edge passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     */
    void
    LocalTree::findEdgeNeighbours(uint32_t idx, uint8_t iedge, u32vector & neighbours, bvector & isghost) const{
        findEdgeNeighbours(idx, false, iedge, neighbours, isghost, false);
    };

    /** Finds local and ghost neighbours of local octant through inode node.
     * Returns a vector (empty if inode is a bound node) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] idx Index of the current octant
     * \param[in] inode Index of node passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     */
    void
    LocalTree::findNodeNeighbours(uint32_t idx, uint8_t inode, u32vector & neighbours, bvector & isghost) const{
        findNodeNeighbours(idx, false, inode, neighbours, isghost, false);
    };

    /** Finds local and ghost neighbours of an octant through iface face.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] oct Pointer to the current octant
     * \param[in] iface Index of face passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     */
    void
    LocalTree::findNeighbours(const Octant* oct, uint8_t iface, u32vector & neighbours, bvector & isghost) const{
        findNeighbours(oct, iface, neighbours, isghost, false);
    };

    /** Finds local and ghost neighbours of an octant through iedge edge.
     * Returns a vector (empty if iedge is a bound edge) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] oct Pointer to the current octant
     * \param[in] iedge Index of edge passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     */
    void
    LocalTree::findEdgeNeighbours(const Octant* oct, uint8_t iedge, u32vector & neighbours, bvector & isghost) const{
        findEdgeNeighbours(oct, iedge, neighbours, isghost, false);
    };

    /** Finds local and ghost neighbours of an octant through inode node.
     * Returns a vector (empty if inode is a bound node) with the index of neighbours
     * in their structure (octants or ghosts) and sets isghost[i] = true if the
     * i-th neighbour is ghost in the local tree.
     * \param[in] oct Pointer to the current octant
     * \param[in] inode Index of node passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     */
    void
    LocalTree::findNodeNeighbours(const Octant* oct, uint8_t inode, u32vector & neighbours, bvector & isghost) const{
        findNodeNeighbours(oct, inode, neighbours, isghost, false);
    };

    /** Finds local neighbours of a ghost octant through iface face.
     * Returns a vector (empty if iface is a bound face) with the index of neighbours
     * \param[in] idx Index of the current octant
     * \param[in] iface Index of face passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     */
	void
    LocalTree::findGhostNeighbours(uint32_t idx, uint8_t iface, u32vector & neighbours) const{
        bvector isghost;
        const Octant* oct = &m_ghosts[idx];
        findNeighbours(oct, iface, neighbours, isghost, true);
    };

    /** Finds local neighbours of a ghost octant through iedge edge.
     * Returns a vector (empty if iedge is a bound edge) with the index of neighbours
     * \param[in] idx Index of the current octant
     * \param[in] iedge Index of edge passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     */
	void
    LocalTree::findGhostEdgeNeighbours(uint32_t idx, uint8_t iedge, u32vector & neighbours) const{
        bvector isghost;
        const Octant* oct = &m_ghosts[idx];
        findEdgeNeighbours(oct, iedge, neighbours, isghost, true);
    };

    /** Finds local neighbours of a ghost octant through inode node.
     * Returns a vector (empty if inode is a bound node) with the index of neighbours
     * \param[in] idx Index of the current octant
     * \param[in] inode Index of node passed through for neighbours finding
     * \param[out] neighbours Vector of neighbours indices in octants/ghosts structure
     */
	void
    LocalTree::findGhostNodeNeighbours(uint32_t idx, uint8_t inode, u32vector & neighbours) const{
        bvector isghost;
        const Octant* oct = &m_ghosts[idx];
        findNodeNeighbours(oct, inode, neighbours, isghost, true);
    };

    // =================================================================================== //
    /*! Given the Morton number of the same-size virtual neighbour and a sorted
     *  list of octans, computes the index from which a neighbour search should
     *  begin.
     *  We are looking for the index of the last octant with a Morton number
     *  lower than the Morton number of the same-size virtual neighbour. If
     *  such an index does not exists (i.e., all the octants have a Morton
     *  number greater than the Morton number of the virtual neighbour, the
     *  search should start from the beginning the octant list).
     * \param[in] sameSizeVirtualNeighMorton Morton number of the same-size
    *  virtual neighbour
     * \param[in] octants list of octants
     * \param[out] searchBeginIdx on output will contain the index from which a
     * neighbour search should begin
     * \param[out] searchBeginMorton on output will contain the Morton of the
     * octant from which a neighbour search should begin
     */
    void
    LocalTree::computeNeighSearchBegin(uint64_t sameSizeVirtualNeighMorton, const octvector &octants, uint32_t *searchBeginIdx, uint64_t *searchBeginMorton) const {

        // Early return if there are no octants
        if (octants.empty()) {
            *searchBeginIdx    = 0;
            *searchBeginMorton = PABLO::INVALID_MORTON;
            return;
        }

        // The search should start from the lower bound if it points to the
        // first octant or to an octant whose Morton number is equal to the
        // the same-size virtual neighbour Morton number. Otherwise, the
        // search should start form the octant preceding the lower bound.
        uint32_t lowerBoundIdx;
        uint64_t lowerBoundMorton;
        findMortonLowerBound(sameSizeVirtualNeighMorton, octants, &lowerBoundIdx, &lowerBoundMorton);

        if (lowerBoundMorton == sameSizeVirtualNeighMorton || lowerBoundIdx == 0) {
            *searchBeginIdx    = lowerBoundIdx;
            *searchBeginMorton = lowerBoundMorton;
        } else {
            *searchBeginIdx    = lowerBoundIdx - 1;
            *searchBeginMorton = octants[*searchBeginIdx].getMorton();
        }

    }

    // =================================================================================== //

    /*! Pre-processing for 2:1 balancing of local tree. Check if there are broken families over processes.
     * \param[in] internal Set to true if the interior octants have to be checked.
     */
    void
    LocalTree::preBalance21(bool internal){

        Octant 			father(m_dim), lastdesc(m_dim);
        uint64_t 		mortonld;
        uint32_t 		nocts;
        uint32_t 		idx, idx2, idx0, last_idx;
        uint32_t 		idx1_gh, idx2_gh;
        int8_t 			marker;
        uint8_t 		nbro;

        //------------------------------------------ //
        // Initialization

        nbro = 0;
        idx=0;
        idx2_gh = idx0 = 0;
        idx1_gh=0;

        nocts   = m_octants.size();
        m_sizeGhosts = m_ghosts.size();
        last_idx=nocts-1;

        //Clean index of ghost brothers in case of coarsening a broken family
        m_lastGhostBros.clear();
        m_firstGhostBros.clear();

        // Set index for start and end check for ghosts
        bool checkend = true;
        bool checkstart = true;
        if (m_ghosts.size()){
            while(m_ghosts[idx2_gh].getMorton() <= m_lastDescMorton){
                idx2_gh++;
                if (idx2_gh > m_sizeGhosts-1) break;
            }
            if (idx2_gh > m_sizeGhosts-1) checkend = false;

            while(m_ghosts[idx1_gh].getMorton() <= m_octants[0].getMorton()){
                idx1_gh++;
                if (idx1_gh > m_sizeGhosts-1) break;
            }
            if (idx1_gh == 0) checkstart = false;
            idx1_gh-=1;
        }

        // Start and End on ghosts
        if (m_ghosts.size() && nocts > 0){
            if (checkstart){
                if (m_ghosts[idx1_gh].buildFather()==m_octants[0].buildFather()){
                    father = m_ghosts[idx1_gh].buildFather();
                    nbro = 0;
                    idx = idx1_gh;
                    marker = m_ghosts[idx].getMarker();
                    while(marker < 0 && m_ghosts[idx].buildFather() == father){

                        //Add ghost index to structure for mapper in case of coarsening a broken family
                        m_firstGhostBros.push_back(idx);

                        nbro++;
                        if (idx==0)
                            break;
                        idx--;
                        marker = m_ghosts[idx].getMarker();
                    }
                    idx = 0;
                    while(idx<nocts && m_octants[idx].buildFather() == father){
                        if(m_octants[idx].getMarker()<0)
                            nbro++;
                        idx++;
                        if(idx==nocts)
                            break;
                    }
                    if (nbro != m_treeConstants->nChildren && idx!=nocts){
                        for(uint32_t ii=0; ii<idx; ii++){
                            if (m_octants[ii].getMarker()<0){
                                m_octants[ii].setMarker(0);
                                m_octants[ii].m_info[Octant::INFO_AUX]=true;
                            }
                        }
                        //Clean ghost index to structure for mapper in case of coarsening a broken family
                        m_firstGhostBros.clear();
                    }
                }
            }

            if (checkend){
                bool checklocal = false;
                if (m_ghosts[idx2_gh].buildFather()==m_octants[nocts-1].buildFather()){
                    father = m_ghosts[idx2_gh].buildFather();
                    if (idx!=nocts){
                        nbro = 0;
                        checklocal = true;
                    }
                    idx = idx2_gh;
                    marker = m_ghosts[idx].getMarker();
                    while(marker < 0 && m_ghosts[idx].buildFather() == father){

                        //Add ghost index to structure for mapper in case of coarsening a broken family
                        m_lastGhostBros.push_back(idx);

                        nbro++;
                        idx++;
                        if(idx == m_sizeGhosts){
                            break;
                        }
                        marker = m_ghosts[idx].getMarker();
                    }
                    idx = nocts-1;
                    if (checklocal){
                        while(m_octants[idx].buildFather() == father ){
                            if (m_octants[idx].getMarker()<0)
                                nbro++;
                            if (idx==0)
                                break;
                            idx--;
                        }
                    }
                    last_idx=idx;
                    if (nbro != m_treeConstants->nChildren && idx!=nocts-1){
                        for(uint32_t ii=idx+1; ii<nocts; ii++){
                            if (m_octants[ii].getMarker()<0){
                                m_octants[ii].setMarker(0);
                                m_octants[ii].m_info[Octant::INFO_AUX]=true;
                            }
                        }
                        //Clean ghost index to structure for mapper in case of coarsening a broken family
                        m_lastGhostBros.clear();
                    }
                }
            }
        }

        // Check first internal octants
        if(getNumOctants()){
            if (internal){
                father = m_octants[0].buildFather();
                lastdesc = father.buildLastDesc();
                mortonld = lastdesc.getMorton();
                nbro = 0;
                for (idx=0; idx<m_treeConstants->nChildren; idx++){
                    if (idx<nocts){
                        // Check if family is complete or to be checked in the internal loop (some brother refined)
                        if (m_octants[idx].getMorton() <= mortonld){
                            nbro++;
                        }
                    }
                }
                if (nbro != m_treeConstants->nChildren)
                    idx0 = nbro;

                // Check and coarse internal octants
                for (idx=idx0; idx<nocts; idx++){
                    if(m_octants[idx].getMarker() < 0 && m_octants[idx].getLevel() > 0){
                        nbro = 0;
                        father = m_octants[idx].buildFather();
                        // Check if family is to be coarsened
                        for (idx2=idx; idx2<idx+m_treeConstants->nChildren; idx2++){
                            if (idx2<nocts){
                                if(m_octants[idx2].getMarker() < 0 && m_octants[idx2].buildFather() == father){
                                    nbro++;
                                }
                            }
                        }
                        if (nbro == m_treeConstants->nChildren){
                            idx = idx2-1;
                        }
                        else{
                            if (idx<=last_idx){
                                m_octants[idx].setMarker(0);
                                m_octants[idx].m_info[Octant::INFO_AUX]=true;
                            }
                        }
                    }
                }
            }
        }
    };

    // =================================================================================== //

    /*! Pre-processing for 2:1 balancing of local tree. Check if there are broken families over processes.
     * \param[out] newmodified Vector of indices of interior octants checked and whose marker is modified.
     */
    void
    LocalTree::preBalance21(u32vector& newmodified){

        Octant 				father(m_dim), lastdesc(m_dim);
        uint64_t 			mortonld;
        uint32_t 			nocts;
        uint32_t 			idx, idx2, idx0, last_idx;
        uint32_t 			idx1_gh, idx2_gh;
        int8_t 			marker;
        uint8_t 			nbro;

        //------------------------------------------ //
        // Initialization

        nbro = 0;
        idx=0;
        idx2_gh = idx0 = 0;
        idx1_gh=0;

        nocts   = m_octants.size();
        m_sizeGhosts = m_ghosts.size();
        last_idx=nocts-1;

        //Clean index of ghost brothers in case of coarsening a broken family
        m_lastGhostBros.clear();
        m_firstGhostBros.clear();

        // Set index for start and end check for ghosts
        bool checkend = true;
        bool checkstart = true;
        if (m_ghosts.size()){
            while(m_ghosts[idx2_gh].getMorton() <= m_lastDescMorton){
                idx2_gh++;
                if (idx2_gh > m_sizeGhosts-1) break;
            }
            if (idx2_gh > m_sizeGhosts-1) checkend = false;

            while(m_ghosts[idx1_gh].getMorton() <= m_octants[0].getMorton()){
                idx1_gh++;
                if (idx1_gh > m_sizeGhosts-1) break;
            }
            if (idx1_gh == 0) checkstart = false;
            idx1_gh-=1;
        }


        // Start and End on ghosts
        if (m_ghosts.size() && nocts > 0){
            if (checkstart){
                if (m_ghosts[idx1_gh].buildFather()==m_octants[0].buildFather()){
                    father = m_ghosts[idx1_gh].buildFather();
                    nbro = 0;
                    idx = idx1_gh;
                    marker = m_ghosts[idx].getMarker();
                    while(marker < 0 && m_ghosts[idx].buildFather() == father){

                        //Add ghost index to structure for mapper in case of coarsening a broken family
                        m_firstGhostBros.push_back(idx);

                        nbro++;
                        if (idx==0)
                            break;
                        idx--;
                        marker = m_ghosts[idx].getMarker();
                    }
                    idx = 0;
                    while(idx<nocts && m_octants[idx].buildFather() == father){
                        if (m_octants[idx].getMarker()<0)
                            nbro++;
                        idx++;
                        if(idx==nocts)
                            break;
                    }
                    if (nbro != m_treeConstants->nChildren && idx!=nocts){
                        for(uint32_t ii=0; ii<idx; ii++){
                            if (m_octants[ii].getMarker()<0){
                                m_octants[ii].setMarker(0);
                                m_octants[ii].m_info[Octant::INFO_AUX]=true;
                                newmodified.push_back(ii);
                            }
                        }
                        //Clean index of ghost brothers in case of coarsening a broken family
                        m_firstGhostBros.clear();
                    }
                }
            }

            if (checkend){
                bool checklocal = false;
                if (m_ghosts[idx2_gh].buildFather()==m_octants[nocts-1].buildFather()){
                    father = m_ghosts[idx2_gh].buildFather();
                    if (idx!=nocts){
                        nbro = 0;
                        checklocal = true;
                    }
                    idx = idx2_gh;
                    marker = m_ghosts[idx].getMarker();
                    while(marker < 0 && m_ghosts[idx].buildFather() == father){

                        //Add ghost index to structure for mapper in case of coarsening a broken family
                        m_lastGhostBros.push_back(idx);

                        nbro++;
                        idx++;
                        if(idx == m_sizeGhosts){
                            break;
                        }
                        marker = m_ghosts[idx].getMarker();
                    }
                    idx = nocts-1;
                    if (checklocal){
                        while(m_octants[idx].buildFather() == father){
                            if (m_octants[idx].getMarker()<0)
                                nbro++;
                            if (idx==0)
                                break;
                            idx--;
                        }
                    }
                    last_idx=idx;
                    if (nbro != m_treeConstants->nChildren && idx!=nocts-1){
                        for(uint32_t ii=idx+1; ii<nocts; ii++){
                            if (m_octants[ii].getMarker()<0){
                                m_octants[ii].setMarker(0);
                                m_octants[ii].m_info[Octant::INFO_AUX]=true;
                                newmodified.push_back(ii);
                            }
                        }
                        //Clean index of ghost brothers in case of coarsening a broken family
                        m_lastGhostBros.clear();
                    }
                }
            }
        }

        // Check first internal octants
        father = m_octants[0].buildFather();
        lastdesc = father.buildLastDesc();
        mortonld = lastdesc.getMorton();
        nbro = 0;
        for (idx=0; idx<m_treeConstants->nChildren; idx++){
            // Check if family is complete or to be checked in the internal loop (some brother refined)
            if (idx<nocts){
                if (m_octants[idx].getMorton() <= mortonld){
                    nbro++;
                }
            }
        }
        if (nbro != m_treeConstants->nChildren)
            idx0 = nbro;

        // Check and coarse internal octants
        for (idx=idx0; idx<nocts; idx++){
            if(m_octants[idx].getMarker() < 0 && m_octants[idx].getLevel() > 0){
                nbro = 0;
                father = m_octants[idx].buildFather();
                // Check if family is to be coarsened
                for (idx2=idx; idx2<idx+m_treeConstants->nChildren; idx2++){
                    if (idx2<nocts){
                        if(m_octants[idx2].getMarker() < 0 && m_octants[idx2].buildFather() == father){
                            nbro++;
                        }
                    }
                }
                if (nbro == m_treeConstants->nChildren){
                    idx = idx2-1;
                }
                else{
                    if (idx<=last_idx){
                        m_octants[idx].setMarker(0);
                        m_octants[idx].m_info[Octant::INFO_AUX]=true;
                        newmodified.push_back(idx);
                    }
                }
            }
        }
    };

    // =================================================================================== //

    /*! 2:1 balancing on level a local tree (refinement wins!)
     * The balance is enforced on octants with the AUX bit set and, if
     * requested, also on new octants.
     * \param[in] doNew Set to true the balance is enforced also on new octants.
     * \param[in] doInterior Set to false if the interior octants are already balanced.
     * \return True if balanced done with some markers modification.
     */
    bool
    LocalTree::localBalance(bool doNew, bool doInterior){

        uint32_t			sizeneigh, modsize;
        u32vector		 	neigh;
        u32vector		 	modified, newmodified;
        uint32_t 			i, idx;
        uint8_t				iface, iedge, inode;
        int8_t				targetmarker;
        vector<bool> 		isghost;
        bool				Bdone = false;
        bool				Bedge = ((m_balanceCodim>1) && (m_dim==3));
        bool				Bnode = (m_balanceCodim==m_dim);

        octvector::iterator 	obegin, oend, it;
        u32vector::iterator 	ibegin, iend, iit;

        //If interior octants have to be balanced
        if(doInterior){
            // First loop on the octants
            obegin = m_octants.begin();
            oend = m_octants.end();
            idx = 0;
            for (it=obegin; it!=oend; ++it){
                bool balanceOctant = it->getBalance();
                if (balanceOctant) {
                    if (doNew) {
                        balanceOctant = (it->m_info[Octant::INFO_AUX] || (it->getMarker() != 0) || it->getIsNewC() || it->getIsNewR());
                    } else {
                        balanceOctant = (it->m_info[Octant::INFO_AUX] || (it->getMarker() != 0));
                    }
                }

                if (balanceOctant){
                    targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(m_octants[idx].getLevel() + m_octants[idx].getMarker()));

                    //Balance through faces
                    for (iface=0; iface<m_treeConstants->nFaces; iface++){
						findNeighbours(idx, iface, neigh, isghost);
						sizeneigh = neigh.size();
						for(i=0; i<sizeneigh; i++){
							if (!isghost[i]){
								{
									if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
										m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
										m_octants[idx].m_info[Octant::INFO_AUX] = true;
										modified.push_back(idx);
										Bdone = true;
									}
									else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
										m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
										m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
										modified.push_back(neigh[i]);
										Bdone = true;
									}
								};
							}
							else{
								{
									if((m_ghosts[neigh[i]].getLevel() + m_ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										m_octants[idx].setMarker(m_ghosts[neigh[i]].getLevel()+m_ghosts[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
										m_octants[idx].m_info[Octant::INFO_AUX] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								};
							}
						}
						targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(m_octants[idx].getLevel() + m_octants[idx].getMarker()));
                    }

                    if (Bedge){
                        //Balance through edges
                        for (iedge=0; iedge<m_treeConstants->nEdges; iedge++){
							findEdgeNeighbours(idx, iedge, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
											m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
											m_octants[idx].m_info[Octant::INFO_AUX] = true;
											modified.push_back(idx);
											Bdone = true;
										}
										else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
											m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
											modified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
								else{
									if((m_ghosts[neigh[i]].getLevel() + m_ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										m_octants[idx].setMarker(m_ghosts[neigh[i]].getLevel()+m_ghosts[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
										m_octants[idx].m_info[Octant::INFO_AUX] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								}
							}
							targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(m_octants[idx].getLevel() + m_octants[idx].getMarker()));
                        }
                    }

                    if (Bnode){
                        //Balance through nodes
                        for (inode=0; inode<m_treeConstants->nNodes; inode++){
							findNodeNeighbours(idx, inode, neigh, isghost);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if (!isghost[i]){
									{
										if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) > (targetmarker + 1) ){
											m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
											m_octants[idx].m_info[Octant::INFO_AUX] = true;
											modified.push_back(idx);
											Bdone = true;
										}
										else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
											m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
											m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
											modified.push_back(neigh[i]);
											Bdone = true;
										}
									};
								}
								else{
									if((m_ghosts[neigh[i]].getLevel() + m_ghosts[neigh[i]].getMarker()) > (targetmarker + 1) ){
										m_octants[idx].setMarker(m_ghosts[neigh[i]].getLevel()+m_ghosts[neigh[i]].getMarker()-1-m_octants[idx].getLevel());
										m_octants[idx].m_info[Octant::INFO_AUX] = true;
										modified.push_back(idx);
										Bdone = true;
									}
								}
							}
							targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(m_octants[idx].getLevel() + m_octants[idx].getMarker()));
                        }
                    }

                }
                idx++;
            }
            // Loop on ghost octants (influence over interior borders)
            obegin = m_ghosts.begin();
            oend = m_ghosts.end();
            idx = 0;
            for (it=obegin; it!=oend; ++it){
                bool balanceOctant = it->getBalance();
                if (balanceOctant) {
                    if (doNew) {
                        balanceOctant = (it->m_info[Octant::INFO_AUX] || it->getIsNewC() || it->getIsNewR());
                    } else {
                        balanceOctant = (it->m_info[Octant::INFO_AUX] || (it->getMarker() != 0));
                    }
                }

                if (balanceOctant){
                    targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(it->getLevel()+it->getMarker()));

                    //Balance through faces
                    for (iface=0; iface<m_treeConstants->nFaces; iface++){
                        if(it->getPbound(iface) == true){
                            neigh.clear();
                            findGhostNeighbours(idx, iface, neigh);
                            sizeneigh = neigh.size();
                            for(i=0; i<sizeneigh; i++){
                                if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
                                    m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
                                    m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
                                    modified.push_back(neigh[i]);
                                    Bdone = true;
                                }
                            }
                        }
                        targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(it->getLevel()+it->getMarker()));
                    }

                    if (Bedge){
                        //Balance through edges
                        for (iedge=0; iedge<m_treeConstants->nEdges; iedge++){
							neigh.clear();
							findGhostEdgeNeighbours(idx, iedge, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
									m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
							targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(it->getLevel()+it->getMarker()));
                        }
                    }

                    if (Bnode){
                        //Balance through nodes
                        for (inode=0; inode<m_treeConstants->nNodes; inode++){
							neigh.clear();
							findGhostNodeNeighbours(idx, inode, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
									m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
							targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(it->getLevel()+it->getMarker()));
                        }
                    }

                }
                idx++;
            }

            // While loop for iterative balancing
            u32vector().swap(newmodified);
            modsize = modified.size();
            while(modsize!=0){
                ibegin = modified.begin();
                iend = modified.end();
                for (iit=ibegin; iit!=iend; ++iit){
                    idx = *iit;
                    if (m_octants[idx].getBalance()){
                        targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));

                        //Balance through faces
                        for (iface=0; iface<m_treeConstants->nFaces; iface++){
                            if(!m_octants[idx].getPbound(iface)){
                                findNeighbours(idx, iface, neigh, isghost);
                                sizeneigh = neigh.size();
                                for(i=0; i<sizeneigh; i++){
                                    if (!isghost[i]){
                                        {
                                            if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
                                                m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
                                                m_octants[idx].m_info[Octant::INFO_AUX] = true;
                                                newmodified.push_back(idx);
                                                Bdone = true;
                                            }
                                            else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
                                                m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
                                                m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
                                                newmodified.push_back(neigh[i]);
                                                Bdone = true;
                                            }
                                        };
                                    }
                                }
                            }
                            targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
                        }

                        if (Bedge){
                            //Balance through edges
                            for (iedge=0; iedge<m_treeConstants->nEdges; iedge++){
								findEdgeNeighbours(idx, iedge, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
												m_octants[idx].m_info[Octant::INFO_AUX] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
												m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
								targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
                            }
                        }

                        if (Bnode){
                            //Balance through nodes
                            for (inode=0; inode<m_treeConstants->nNodes; inode++){
								findNodeNeighbours(idx, inode, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
												m_octants[idx].m_info[Octant::INFO_AUX] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
												m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
								targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
                            }
                        }

                    }
                }
                preBalance21(newmodified);
                u32vector().swap(modified);
                swap(modified,newmodified);
                modsize = modified.size();
                u32vector().swap(newmodified);
            }// end while

        }
        else{

            // Loop on ghost octants (influence over interior borders)
            obegin = m_ghosts.begin();
            oend = m_ghosts.end();
            idx = 0;
            for (it=obegin; it!=oend; ++it){
                bool balanceOctant = it->getBalance();
                if (balanceOctant) {
                    if (doNew) {
                        balanceOctant = (it->m_info[Octant::INFO_AUX] || it->getIsNewC() || it->getIsNewR());
                    } else {
                        balanceOctant = it->m_info[Octant::INFO_AUX];
                    }
                }

                if (balanceOctant){
                    targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(it->getLevel()+it->getMarker()));

                    //Balance through faces
                    for (iface=0; iface<m_treeConstants->nFaces; iface++){
                        if(it->getPbound(iface) == true){
                            neigh.clear();
                            findGhostNeighbours(idx, iface, neigh);
                            sizeneigh = neigh.size();
                            for(i=0; i<sizeneigh; i++){
                                if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
                                    m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
                                    m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
                                    modified.push_back(neigh[i]);
                                    Bdone = true;
                                }
                            }
                        }
                        targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(it->getLevel()+it->getMarker()));
                    }

                    if (Bedge){
                        //Balance through edges
                        for (iedge=0; iedge<m_treeConstants->nEdges; iedge++){
							neigh.clear();
							findGhostEdgeNeighbours(idx, iedge, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
									m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
							targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(it->getLevel()+it->getMarker()));
                        }
                    }

                    if (Bnode){
                        //Balance through nodes
                        for (inode=0; inode<m_treeConstants->nNodes; inode++){
							neigh.clear();
							findGhostNodeNeighbours(idx, inode, neigh);
							sizeneigh = neigh.size();
							for(i=0; i<sizeneigh; i++){
								if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
									m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
									m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
									modified.push_back(neigh[i]);
									Bdone = true;
								}
							}
							targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(it->getLevel()+it->getMarker()));
                        }
                    }

                }
                idx++;
            }

            // While loop for iterative balancing
            u32vector().swap(newmodified);
            modsize = modified.size();
            while(modsize!=0){
                ibegin = modified.begin();
                iend = modified.end();
                for (iit=ibegin; iit!=iend; ++iit){
                    idx = *iit;
                    if (m_octants[idx].getBalance()){
                        targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));

                        //Balance through faces
                        for (iface=0; iface<m_treeConstants->nFaces; iface++){
                            if(!m_octants[idx].getPbound(iface)){
                                findNeighbours(idx, iface, neigh, isghost);
                                sizeneigh = neigh.size();
                                for(i=0; i<sizeneigh; i++){
                                    if (!isghost[i]){
                                        {
                                            if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
                                                m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
                                                m_octants[idx].m_info[Octant::INFO_AUX] = true;
                                                newmodified.push_back(idx);
                                                Bdone = true;
                                            }
                                            else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
                                                m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
                                                m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
                                                newmodified.push_back(neigh[i]);
                                                Bdone = true;
                                            }
                                        };
                                    }
                                }
                            }
                            targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
                        }

                        if (Bedge){
                            //Balance through edges
                            for (iedge=0; iedge<m_treeConstants->nEdges; iedge++){
								findEdgeNeighbours(idx, iedge, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
												m_octants[idx].m_info[Octant::INFO_AUX] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
												m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
								targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
                            }
                        }

                        if (Bnode){
                            //Balance through nodes
                            for (inode=0; inode<m_treeConstants->nNodes; inode++){
								findNodeNeighbours(idx, inode, neigh, isghost);
								sizeneigh = neigh.size();
								for(i=0; i<sizeneigh; i++){
									if (!isghost[i]){
										{
											if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) >  (targetmarker + 1)){
												m_octants[idx].setMarker(m_octants[neigh[i]].getLevel()+m_octants[neigh[i]].getMarker()-m_octants[idx].getLevel()-1);
												m_octants[idx].m_info[Octant::INFO_AUX] = true;
												newmodified.push_back(idx);
												Bdone = true;
											}
											else if((m_octants[neigh[i]].getLevel() + m_octants[neigh[i]].getMarker()) < (targetmarker - 1)){
												m_octants[neigh[i]].setMarker(targetmarker-m_octants[neigh[i]].getLevel()-1);
												m_octants[neigh[i]].m_info[Octant::INFO_AUX] = true;
												newmodified.push_back(neigh[i]);
												Bdone = true;
											}
										};
									}
								}
								targetmarker = min(TreeConstants::MAX_LEVEL, int8_t(m_octants[idx].getLevel()+m_octants[idx].getMarker()));
                            }
                        }

                    }
                }
                preBalance21(newmodified);
                u32vector().swap(modified);
                swap(modified,newmodified);
                modsize = modified.size();
                u32vector().swap(newmodified);
            }// end while
            obegin = oend = m_octants.end();
            ibegin = iend = modified.end();
        }
        return Bdone;
        // Pay attention : info[15] may be true after local balance for some octants
    };


    // =================================================================================== //

    /*! Compute and store in m_intersections the intersections of the local tree.
     */
    void
    LocalTree::computeIntersections() {

		octvector::iterator 	it, obegin, oend;
		Intersection 			intersection;
		u32vector 				neighbours;
		vector<bool>			isghost;
		uint32_t 				counter, idx;
		uint32_t 				i, nsize;
		uint8_t 				iface, iface2;

		m_intersections.clear();
		m_intersections.reserve(2*3*m_octants.size());

		counter = idx = 0;

		// Loop on ghosts
		obegin = m_ghosts.begin();
		oend = m_ghosts.end();
		for (it = obegin; it != oend; ++it){
			for (iface = 0; iface < m_dim; iface++){
				iface2 = iface*2;
				findGhostNeighbours(idx, iface2, neighbours);
				nsize = neighbours.size();
				if (!(it->m_info[iface2])){
					//Internal intersection
					for (i = 0; i < nsize; i++){
						intersection.m_dim = m_dim;
						intersection.m_finer = getGhostLevel(idx) >= getLevel((int)neighbours[i]);
						intersection.m_out = intersection.m_finer;
						intersection.m_outisghost = intersection.m_finer;
						intersection.m_owners[0]  = neighbours[i];
						intersection.m_owners[1] = idx;
						intersection.m_iface = m_treeConstants->oppositeFace[iface2] - (getGhostLevel(idx) >= getLevel((int)neighbours[i]));
						intersection.m_isnew = false;
						intersection.m_isghost = true;
						intersection.m_bound = false;
						intersection.m_pbound = true;
						m_intersections.push_back(intersection);
						counter++;
					}
				}
				else{
					//Periodic intersection
					for (i = 0; i < nsize; i++){
						intersection.m_dim = m_dim;
						intersection.m_finer = getGhostLevel(idx) >= getLevel((int)neighbours[i]);
						intersection.m_out = intersection.m_finer;
						intersection.m_outisghost = intersection.m_finer;
						intersection.m_owners[0]  = neighbours[i];
						intersection.m_owners[1] = idx;
						intersection.m_iface = m_treeConstants->oppositeFace[iface2] - (getGhostLevel(idx) >= getLevel((int)neighbours[i]));
						intersection.m_isnew = false;
						intersection.m_isghost = true;
						intersection.m_bound = true;
						intersection.m_pbound = true;
						m_intersections.push_back(intersection);
						counter++;
					}
				}
			}
			idx++;
		}

		// Loop on octants
		idx=0;
		obegin = m_octants.begin();
		oend = m_octants.end();
		for (it = obegin; it != oend; ++it){
			for (iface = 0; iface < m_dim; iface++){
				iface2 = iface*2;
				findNeighbours(idx, iface2, neighbours, isghost);
				nsize = neighbours.size();
				if (nsize) {
					if (!(it->m_info[iface2])){
						//Internal intersection
						for (i = 0; i < nsize; i++){
							if (isghost[i]){
								intersection.m_dim = m_dim;
								intersection.m_owners[0] = idx;
								intersection.m_owners[1] = neighbours[i];
								intersection.m_finer = (nsize>1);
								intersection.m_out = (nsize>1);
								intersection.m_outisghost = (nsize>1);
								intersection.m_iface = iface2 + (nsize>1);
								intersection.m_isnew = false;
								intersection.m_isghost = true;
								intersection.m_bound = false;
								intersection.m_pbound = true;
								m_intersections.push_back(intersection);
								counter++;
							}
							else{
								intersection.m_dim = m_dim;
								intersection.m_owners[0] = idx;
								intersection.m_owners[1] = neighbours[i];
								intersection.m_finer = (nsize>1);
								intersection.m_out = (nsize>1);
								intersection.m_outisghost = false;
								intersection.m_iface = iface2 + (nsize>1);
								intersection.m_isnew = false;
								intersection.m_isghost = false;
								intersection.m_bound = false;
								intersection.m_pbound = false;
								m_intersections.push_back(intersection);
								counter++;
							}
						}
					}
					else{
						//Periodic intersection
						for (i = 0; i < nsize; i++){
							if (isghost[i]){
								intersection.m_dim = m_dim;
								intersection.m_owners[0] = idx;
								intersection.m_owners[1] = neighbours[i];
								intersection.m_finer = (nsize>1);
								intersection.m_out = intersection.m_finer;
								intersection.m_outisghost = intersection.m_finer;
								intersection.m_iface = iface2 + (nsize>1);
								intersection.m_isnew = false;
								intersection.m_isghost = true;
								intersection.m_bound = true;
								intersection.m_pbound = true;
								m_intersections.push_back(intersection);
								counter++;
							}
							else{
								intersection.m_dim = m_dim;
								intersection.m_owners[0] = idx;
								intersection.m_owners[1] = neighbours[i];
								intersection.m_finer = (nsize>1);
								intersection.m_out = intersection.m_finer;
								intersection.m_outisghost = false;
								intersection.m_iface = iface2 + (nsize>1);
								intersection.m_isnew = false;
								intersection.m_isghost = false;
								intersection.m_bound = true;
								intersection.m_pbound = false;
								m_intersections.push_back(intersection);
								counter++;
							}
						}
					}
				}
				else{
					//Boundary intersection
					intersection.m_dim = m_dim;
					intersection.m_owners[0] = idx;
					intersection.m_owners[1] = idx;
					intersection.m_finer = 0;
					intersection.m_out = 0;
					intersection.m_outisghost = false;
					intersection.m_iface = iface2;
					intersection.m_isnew = false;
					intersection.m_isghost = false;
					intersection.m_bound = true;
					intersection.m_pbound = false;
					m_intersections.push_back(intersection);
					counter++;
				}
				if (it->m_info[iface2+1]){
					if (!(m_periodic[iface2+1])){
						//Boundary intersection
						intersection.m_dim = m_dim;
						intersection.m_owners[0] = idx;
						intersection.m_owners[1] = idx;
						intersection.m_finer = 0;
						intersection.m_out = 0;
						intersection.m_outisghost = false;
						intersection.m_iface = iface2+1;
						intersection.m_isnew = false;
						intersection.m_isghost = false;
						intersection.m_bound = true;
						intersection.m_pbound = false;
						m_intersections.push_back(intersection);
						counter++;
					}
					else{
						//Periodic intersection
						findNeighbours(idx, iface2+1, neighbours, isghost);
						nsize = neighbours.size();
						for (i = 0; i < nsize; i++){
							if (isghost[i]){
								intersection.m_dim = m_dim;
								intersection.m_owners[0] = idx;
								intersection.m_owners[1] = neighbours[i];
								intersection.m_finer = (nsize>1);
								intersection.m_out = intersection.m_finer;
								intersection.m_outisghost = intersection.m_finer;
								intersection.m_iface = iface2 + (nsize>1);
								intersection.m_isnew = false;
								intersection.m_isghost = true;
								intersection.m_bound = true;
								intersection.m_pbound = true;
								m_intersections.push_back(intersection);
								counter++;
							}
							else{
								intersection.m_dim = m_dim;
								intersection.m_owners[0] = idx;
								intersection.m_owners[1] = neighbours[i];
								intersection.m_finer = (nsize>1);
								intersection.m_out = intersection.m_finer;
								intersection.m_outisghost = false;
								intersection.m_iface = iface2 + (nsize>1);
								intersection.m_isnew = false;
								intersection.m_isghost = false;
								intersection.m_bound = true;
								intersection.m_pbound = false;
								m_intersections.push_back(intersection);
								counter++;
							}
						}
					}
				}
			}
			idx++;
		}
		intervector(m_intersections).swap(m_intersections);
	}

    // =================================================================================== //
    /*! Find an input Morton in octants and return the local idx
     * \param[in] targetMorton Morton index to be found.
     * \return Local index of the target octant (=nocts if target Morton not found).
     */
    uint32_t
    LocalTree::findMorton(uint64_t targetMorton) const {
        return findMorton(targetMorton, m_octants);
    };

    // =================================================================================== //
    /*! Find an input Morton in ghosts and return the local idx
     * \param[in] targetMorton Morton index to be found.
     * \return Index of the target ghost octant (=nghosts if target Morton not found).
     */
    uint32_t
    LocalTree::findGhostMorton(uint64_t targetMorton) const {
        return findMorton(targetMorton, m_ghosts);
    };

    // =================================================================================== //
    /*! Find the index of the octant with the specified Morton in the given
     *  sorted list of octants. If the requested octant is not in the given
     *  list, the index of the end (i.e., the element after the last element
     *  is returned).
     * \param[in] targetMorton is the Morton index to be found.
     * \param[in] octants list of octants
     * \return Local index of the target octant (=nocts if target Morton not found).
     */
    uint32_t
    LocalTree::findMorton(uint64_t targetMorton, const octvector &octants) const {

        uint32_t lowerBoundIdx;
        uint64_t lowerBoundMorton;
        findMortonLowerBound(targetMorton, octants, &lowerBoundIdx, &lowerBoundMorton);

        uint32_t targetIdx;
        if (lowerBoundMorton == targetMorton) {
            targetIdx = lowerBoundIdx;
        } else {
            targetIdx = octants.size();
        }

        return targetIdx;
    };

    // =================================================================================== //
    /*! Given a target Morton number and a sorted list of octants, finds the
     *  index of the first octant whose Morton number does not compare less
     *  than the target Morton number (in other words, the index of the first
     *  octant whose Morton number is greater or equal than the target Morton
     *  number). If the target Morton is greater than the Morton number of the
     *  last element, the index of the past-the-element element is returned.
     * \param[in] targetMorton is the Morton index to be found.
     * \param[in] octants list of octants
     * \param[out] lowerBoundIdx on output will contain the index of first
     * octant whose Morton number does not compare less than the target Morton
     * number. If the target Morton numer is greater than the Morton number of
     * the last element, the index of the past-the-element element is returned
     * \param[out] lowerBoundMorton on output will contain the Morton associated
     * with the lower bound. If the target Morton is greater than the Morton
     * number of the last element, the maximum finite value representable by
     * the numeric type is returned
     */
    void
    LocalTree::findMortonLowerBound(uint64_t targetMorton, const octvector &octants, uint32_t *lowerBoundIdx, uint64_t *lowerBoundMorton) const {

        uint32_t nOctants = octants.size();

        uint32_t lowIndex  = 0;
        uint32_t highIndex = nOctants;
        uint32_t midIndex  = nOctants;
        uint64_t midMorton = PABLO::INVALID_MORTON;
        while (lowIndex < highIndex) {
            midIndex  = lowIndex + (highIndex - lowIndex) / 2;
            midMorton = octants[midIndex].getMorton();
            if (targetMorton < midMorton) {
                highIndex = midIndex;
            }
            else if (targetMorton > midMorton) {
                lowIndex = midIndex + 1;
            }
            else {
                *lowerBoundIdx  = midIndex;
                *lowerBoundMorton = midMorton;

                return;
            }
        }

        *lowerBoundIdx = lowIndex;
        if (*lowerBoundIdx == midIndex) {
            *lowerBoundMorton = midMorton;
        }
        else if (*lowerBoundIdx < nOctants) {
            *lowerBoundMorton = octants[*lowerBoundIdx].getMorton();
        }
        else {
            *lowerBoundMorton = PABLO::INVALID_MORTON;
        }

    }

    // =================================================================================== //
    /*! Given a target Morton number and a sorted list of octants, finds the
     *  index of the first octant whose Morton number is greater than the
     *  target Morton number. If the target Morton is greater than the Morton
     *  number of the last element, the index of the past-the-element element
     *  is returned.
     * \param[in] targetMorton is the Morton index to be found.
     * \param[in] octants list of octants
     * \param[out] upperBoundIdx on output will contain the index of the first
     * octant whose Morton number is greater than the target Morton number. If
     * the target Morton numer is greater than the Morton number of the last
     * element, the index of the past-the-last-element is returned
     * \param[out] upperBoundMorton on output will contain the Morton associated
     * with the upper bound. If no upper bound was found, the maximum finite
     * value representable by the numeric type is returned
     */
    void
    LocalTree::findMortonUpperBound(uint64_t targetMorton, const octvector &octants, uint32_t *upperBoundIdx, uint64_t *upperBoundMorton) const {

        uint32_t nOctants = octants.size();

        uint32_t lowIndex  = 0;
        uint32_t highIndex = nOctants;
        uint32_t midIndex  = nOctants;
        uint64_t midMorton = PABLO::INVALID_MORTON;
        while (lowIndex < highIndex) {
            midIndex  = lowIndex + (highIndex - lowIndex) / 2;
            midMorton = octants[midIndex].getMorton();
            if (targetMorton < midMorton) {
                highIndex = midIndex;
            }
            else {
                lowIndex = midIndex + 1;
            }
        }

        *upperBoundIdx = lowIndex;
        if (*upperBoundIdx == midIndex) {
            *upperBoundMorton = midMorton;
        }
        else if (*upperBoundIdx < nOctants) {
            *upperBoundMorton = octants[*upperBoundIdx].getMorton();
        }
        else {
            *upperBoundMorton = PABLO::INVALID_MORTON;
        }

    }

    // =================================================================================== //

    /** Compute the connectivity of octants and store the coordinates of nodes.
     */
    void
    LocalTree::computeConnectivity(){
        vector<uint64_t>                             mortonList;
        unordered_map<uint64_t, array<uint32_t, 3> > nodeCoords;
        unordered_map<uint64_t, vector<uint64_t> >   nodeOctants;
        uint32_t                                     noctants = getNumOctants();
        uint32_t                                     nghosts  = m_sizeGhosts;



        // Gather node information
        mortonList.reserve(noctants);
        nodeCoords.reserve(noctants);
        nodeOctants.reserve(noctants);

        for (uint64_t n = 0; n < (noctants + nghosts); n++){
            const Octant *octant;
            if (n < noctants) {
                uint32_t octantId = n;
                octant = &(m_octants[octantId]);
            } else {
                uint32_t octantId = n - noctants;
                octant = &(m_ghosts[octantId]);
            }

            for (uint8_t i = 0; i < m_treeConstants->nNodes; ++i){
                u32array3 node;
                octant->getLogicalNode(node, i);

                uint64_t morton = octant->computeNodePersistentKey(node);
                if (nodeCoords.count(morton) == 0) {
                    mortonList.push_back(morton);
                    nodeCoords.insert({{morton, std::move(node)}});
                    nodeOctants[morton].reserve(m_treeConstants->nNodes);
                }

                nodeOctants[morton].push_back(n);
            }
        }
        std::sort(mortonList.begin(), mortonList.end());

        // Build node list and connectivity
        m_nodes.reserve(mortonList.size());
        m_connectivity.resize(noctants);
        m_ghostsConnectivity.resize(nghosts);

        uint32_t nodeId = 0;
        for (uint64_t morton : mortonList) {
            m_nodes.emplace_back(std::move(nodeCoords.at(morton)));
            for (uint64_t n : nodeOctants.at(morton)) {
                std::vector<uint32_t> *octantConnect;
                if (n < noctants) {
                    uint32_t octantId = n;
                    octantConnect = &(m_connectivity[octantId]);
                } else {
                    uint32_t octantId = n - noctants;
                    octantConnect = &(m_ghostsConnectivity[octantId]);
                }

                if (octantConnect->size() == 0) {
                    octantConnect->reserve(m_treeConstants->nNodes);
                }
                octantConnect->push_back(nodeId);
            }
            nodeId++;
        }
    };

    /*! Clear nodes vector and connectivity of octants of local tree
     */
    void
    LocalTree::clearConnectivity(){
        u32arr3vector().swap(m_nodes);
        u32vector2D().swap(m_connectivity);
        u32vector2D().swap(m_ghostsConnectivity);
    };

    /*! Updates nodes vector and connectivity of octants of local tree
     */
    void
    LocalTree::updateConnectivity(){
        clearConnectivity();
        computeConnectivity();
    };

    // =================================================================================== //

}
