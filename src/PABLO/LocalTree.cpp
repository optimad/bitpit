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
#include "LocalTree.hpp"
#include "morton.hpp"

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
        return m_ghosts.size();
    };

    /*! Get the number of the octants in the local tree.
     * \return Number of local octants.
     */
    uint32_t
    LocalTree::getNumOctants() const{
        return m_octants.size();
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
        if(getNumOctants() > 0){
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
        if(getNumOctants() > 0){
            octvector::const_iterator lastOctant = m_octants.end() - 1;
            uint32_t x,y,z,delta;
            delta = (uint32_t)(1<<((uint8_t)m_treeConstants->maxLevel - lastOctant->m_level)) - 1;
            x = lastOctant->getLogicalCoordinates(0) + delta;
            y = lastOctant->getLogicalCoordinates(1) + delta;
            z = lastOctant->getLogicalCoordinates(2) + (m_dim-2)*delta;
            Octant lastDesc = Octant(m_dim, m_treeConstants->maxLevel,x,y,z);
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
     * \return	true if additional refinement is needed in order to satisfy
     * the specified markers
     */
    bool
    LocalTree::refine(u32vector & mapidx){
        // Current number of octants
        uint32_t nOctants = m_octants.size();
        if (nOctants == 0) {
            return false;
        }

        // Validate markers
        //
        // Not all octants marked for refinement can be really refined, for
        // example octants cannot be refined further than the maximum level.
        uint32_t firstRefinedIdx   = 0;
        uint32_t nValidRefinements = 0;
        for (uint32_t idx=0; idx<nOctants; idx++) {
            // Skip octants not marked for refinement
            Octant &octant = m_octants[idx];
            if(octant.getMarker()<= 0){
                continue;
            }

            // Octants cannot be refined further than the maximum level
            if(octant.getLevel()>=m_treeConstants->maxLevel){
                octant.setMarker(0);
                continue;
            }

            // The octant will be refined
            ++nValidRefinements;
            if (nValidRefinements == 1) {
                firstRefinedIdx = idx;
            }
        }

        // Early return if no octants need to be refined
        if (nValidRefinements == 0) {
            return false;
        }

        // Resize octant container
        //
        // We want to be sure the container capacity is equal to its size.
        uint32_t nFutureOctants = nOctants + (m_treeConstants->nChildren - 1) * nValidRefinements;

        m_octants.reserve(nFutureOctants);
        m_octants.resize(nFutureOctants, Octant(m_dim));
        m_octants.shrink_to_fit();

        // Initialize mapping
        if(!mapidx.empty()){
            mapidx.resize(nFutureOctants);
        }

        // Refine the octants
        //
        // We process the octants backwards (starting from the back of their
        // container):
        //  - if an octant doesn't need refinement, it will only be moved to
        //    its new position;
        //  - if an octant needs to be refined, it will be replaced with its
        //    children.
        bool refinementCompleted = true;

        uint32_t futureIdx = nFutureOctants - 1;
        for (uint32_t n=0; n<(nOctants - firstRefinedIdx); ++n) {
            uint32_t idx = nOctants - n - 1;
            Octant &octant = m_octants[idx];
            if(octant.getMarker()<=0){
                // The octant is not refiend, we need to move it to its new
                // position in the container.
                if (futureIdx != idx) {
                    m_octants[futureIdx] = std::move(octant);
                }

                // Update the mapping
                if(!mapidx.empty()){
                    mapidx[futureIdx] = idx;
                }

                // Update future octant index
                --futureIdx;
            } else {
                // Move original octant out of the container
                Octant fatherOctant = std::move(octant);

                // Create children
                uint8_t nChildren = fatherOctant.countChildren();
                uint32_t firstChildIdx = futureIdx - (nChildren - 1);
                fatherOctant.buildChildren(m_octants.data() + firstChildIdx);

                // Update the mapping
                if(!mapidx.empty()){
                    for (uint8_t i=0; i<nChildren; i++){
                        mapidx[firstChildIdx + i] = idx;
                    }
                }

                // Check if more refinement is needed to satisfy the markers
                if (refinementCompleted) {
                    refinementCompleted = (m_octants[firstChildIdx].getMarker() <= 0);
                }

                // Update local max depth
                uint8_t childrenLevel = m_octants[firstChildIdx].getLevel();
                if (childrenLevel > m_localMaxDepth){
                    m_localMaxDepth = childrenLevel;
                }

                // Update future octant index
                futureIdx -= nChildren;
            }
        }

        // Update the mapping for cells that have not been modified
        if(!mapidx.empty()){
            for (uint32_t idx=0; idx<firstRefinedIdx; ++idx) {
                mapidx[idx] = idx;
            }
        }

        return (!refinementCompleted);

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

        uint32_t nInitialOctants = getNumOctants();
        if (nInitialOctants == 0) {
            return false;
        }

        uint32_t nInitialGhosts = getNumGhosts();

        nbro = nend = nstart = 0;
        nidx = offset = 0;

        idx2_gh = 0;
        idx1_gh = nInitialGhosts - 1;

        //------------------------------------------ //

        // Set index for start and end check for ghosts
        if (m_ghosts.size()){
            bool check = true;
            while(check){
                check = idx1_gh < nInitialGhosts;
                if (check){
                    check = m_ghosts[idx1_gh].getMorton() > m_firstDescMorton;
                }
                if (check) idx1_gh--;
            }

            check = true;
            while(check){
                check = idx2_gh < nInitialGhosts;
                if (check){
                    check = m_ghosts[idx2_gh].getMorton() < m_lastDescMorton;
                }
                if (check) idx2_gh++;
            }
        }

        // Check and coarse internal octants
        for (idx=0; idx<nInitialOctants; idx++){
            if(m_octants[idx].getMarker() < 0 && m_octants[idx].getLevel() > 0){
                nbro = 0;
                father = m_octants[idx].buildFather();
                // Check if family is to be refined
                for (idx2=idx; idx2<idx+m_treeConstants->nChildren; idx2++){
                    if (idx2<nInitialOctants){
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
        uint32_t nblock = nInitialOctants;
        uint32_t nfchild = first_child_index.size();
        if (nidx!=0){
            nblock = nInitialOctants - nidx*nchm1;
            nidx = 0;
            for (idx=0; idx<nblock; idx++){
                if (idx+offset < nInitialOctants){
                    if (nidx < nfchild){
                        if (idx+offset == first_child_index[nidx]){
                            markerfather = -m_treeConstants->maxLevel;
                            father = m_octants[idx+offset].buildFather();
                            for (uint32_t iii=0; iii<Octant::INFO_ITEM_COUNT; iii++){
                                father.m_info[iii] = false;
                            }
                            father.setGhostLayer(-1);
                            for(idx2=0; idx2<m_treeConstants->nChildren; idx2++){
                                if (idx2 < nInitialOctants){
                                    if (markerfather < m_octants[idx+offset+idx2].getMarker()+1){
                                        markerfather = m_octants[idx+offset+idx2].getMarker()+1;
                                    }
                                    for (uint32_t iii=0; iii<Octant::INFO_ITEM_COUNT; iii++){
                                        father.m_info[iii] = father.m_info[iii] || m_octants[idx+offset+idx2].m_info[iii];
                                    }
                                }
                            }
                            father.m_info[Octant::INFO_NEW4COARSENING] = true;
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
        nInitialOctants = m_octants.size();
        if(mapsize > 0){
            mapidx.resize(nInitialOctants);
        }

        //Check ghosts
        if (m_ghosts.size()){
            // Start on ghosts
            if (nInitialOctants > 0 && idx1_gh < nInitialGhosts){
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
                    if (idx==nInitialOctants-1) wstop = true;
                    while(marker < 0 && m_octants[idx].buildFather() == father){
                        nbro++;
                        nstart++;
                        if (wstop){
                            break;
                        }
                        idx++;
                        if (idx==nInitialOctants-1){
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
                    for (idx=0; idx<nInitialOctants; idx++){
                        if (idx+offset < nInitialOctants){
                            m_octants[idx] = m_octants[idx+offset];
                            if(mapsize > 0) mapidx[idx] = mapidx[idx+offset];
                        }
                    }
                    m_octants.resize(nInitialOctants-offset, Octant(m_dim));
                    m_octants.shrink_to_fit();
                    nInitialOctants = m_octants.size();
                    if(mapsize > 0){
                        mapidx.resize(nInitialOctants);
                    }
                }
            }


            //Verify family between more then two processes
            if (nInitialOctants > 0 && idx2_gh < nInitialGhosts){

                if (m_ghosts[idx2_gh].buildFather() == father){

                    if (m_ghosts[idx2_gh].buildFather() == m_octants[nInitialOctants-1].buildFather()){

                        uint64_t idx22_gh = idx2_gh;
                        marker = m_ghosts[idx22_gh].getMarker();
                        while(marker < 0 && m_ghosts[idx22_gh].buildFather() == father){
                            nbro++;
                            idx22_gh++;
                            if(idx22_gh == nInitialGhosts){
                                break;
                            }
                            marker = m_ghosts[idx22_gh].getMarker();
                        }

                        if (nbro == m_treeConstants->nChildren){
                            m_octants.clear();
                            nInitialOctants = 0;
                            if(mapsize > 0){
                                mapidx.clear();
                            }
                        }
                    }

                }

            }


            // End on ghosts
            if (nInitialOctants > 0 && idx2_gh < nInitialGhosts){
                if (m_ghosts[idx2_gh].buildFather() == m_octants[nInitialOctants-1].buildFather()){
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
                        if(idx == nInitialGhosts){
                            break;
                        }
                        marker = m_ghosts[idx].getMarker();
                    }
                    nend = 0;
                    idx = nInitialOctants-1;
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
                            father.m_info[iii] = father.m_info[iii] || m_octants[nInitialOctants-idx-1].m_info[iii];
                        }
                    }
                    father.m_info[Octant::INFO_NEW4COARSENING] = true;
                    father.setGhostLayer(-1);
                    //Impossible in this version
                    //                if (markerfather < 0 && mapsize == 0){
                        //                    docoarse = true;
                        //                }
                    father.setMarker(markerfather);
                    m_octants.resize(nInitialOctants-offset, Octant(m_dim));
                    m_octants.push_back(father);
                    m_octants.shrink_to_fit();
                    nInitialOctants = m_octants.size();
                    if(mapsize > 0){
                        mapidx.resize(nInitialOctants);
                    }
                }
            }

        }//end if ghosts size

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

        for (Octant &octant : m_octants){
            octant.setMarker(1);
        }

        return refine(mapidx);

    };

    // =================================================================================== //

    /*! Refine local tree: corse one time all the octants
     * \param[out] mapidx mpaidx[i] = index in old octants vector of the new i-th octant (index of first child if octant is new after coarsening)
     * \return	true if coarsening can continue (impossible in this version)
     */
    bool
    LocalTree::globalCoarse(u32vector & mapidx){

        for (Octant &octant : m_octants){
            octant.setMarker(-1);
        }
        for (Octant &octant : m_ghosts){
            octant.setMarker(-1);
        }

        return coarse(mapidx);

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

        uint32_t nInitialOctants = getNumOctants();
        if (nInitialOctants>0){

            idx = 0;
            if (m_octants[idx].getMorton() < partLastDesc){

                Octant father0 = m_octants[idx].buildFather();
                Octant father = father0;

                while(father == father0 && idx < nInitialOctants){
                    toDelete++;
                    idx++;
                    if (idx<nInitialOctants) father = m_octants[idx].buildFather();
                }

                if (nInitialOctants>toDelete){
                    for(idx=0; idx<nInitialOctants-toDelete; idx++){
                        m_octants[idx] = m_octants[idx+toDelete];
                        if (mapsize>0) mapidx[idx] = mapidx[idx+toDelete];
                    }
                    m_octants.resize(nInitialOctants-toDelete, Octant(m_dim));
                    if (mapsize>0){
                        mapidx.resize(m_octants.size());
                    }
                }
                else{
                    m_octants.clear();
                    mapidx.clear();
                }
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
     * \param[in,out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[in,out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     * \param[in] onlyinternal A boolean flag to specify if neighbours have to be found among all the octants (false) or only among the internal ones (true).
     * \param[in] append A boolean flag to specify if neighbours will be appended to the given vector or if the given vectors will be cleared before adding the neighbours.
     */
    void
    LocalTree::findNeighbours(const Octant* oct, uint8_t iface, u32vector & neighbours, bvector & isghost, bool onlyinternal, bool append) const{

        if (!append) {
            isghost.clear();
            neighbours.clear();
        }

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

        // Initialize search
        uint32_t candidateIdx    = 0;
        uint64_t candidateMorton = 0;

        uint64_t neighArea = 0;
        uint64_t faceArea  = oct->getLogicalArea();

        uint32_t size = oct->getLogicalSize();

        std::array<uint32_t, 3> coord = oct->getLogicalCoordinates();
        if (isperiodic) {
            std::array<int64_t, 3> periodicOffset = getPeriodicOffset(*oct, iface);
            coord[0] = static_cast<uint32_t>(coord[0] + periodicOffset[0]);
            coord[1] = static_cast<uint32_t>(coord[1] + periodicOffset[1]);
            coord[2] = static_cast<uint32_t>(coord[2] + periodicOffset[2]);
        }

        const int8_t (&cxyz)[3] = m_treeConstants->normals[iface];

        // Compute same-size virtual neighbour information
        std::array<int64_t, 3> sameSizeVirtualNeighOffset = computeFirstVirtualNeighOffset(level, iface, level);
        std::array<uint32_t, 3> sameSizeVirtualNeighCoord = coord;
        sameSizeVirtualNeighCoord[0] = static_cast<uint32_t>(sameSizeVirtualNeighCoord[0] + sameSizeVirtualNeighOffset[0]);
        sameSizeVirtualNeighCoord[1] = static_cast<uint32_t>(sameSizeVirtualNeighCoord[1] + sameSizeVirtualNeighOffset[1]);
        sameSizeVirtualNeighCoord[2] = static_cast<uint32_t>(sameSizeVirtualNeighCoord[2] + sameSizeVirtualNeighOffset[2]);
        uint64_t sameSizeVirtualNeighMorton = PABLO::computeMorton(m_dim, sameSizeVirtualNeighCoord[0], sameSizeVirtualNeighCoord[1], sameSizeVirtualNeighCoord[2]);

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
        // This is the Morton number of the last descendant of the same-size
        // virtual neighbour.
        uint8_t maxNeighLevel = getMaxNeighLevel(*oct);
        std::array<int64_t, 3> lastCandidateOffset = computeLastVirtualNeighOffset(level, iface, maxNeighLevel);
        std::array<uint32_t, 3> lastCandidateCoord = coord;
        lastCandidateCoord[0] = static_cast<uint32_t>(lastCandidateCoord[0] + lastCandidateOffset[0]);
        lastCandidateCoord[1] = static_cast<uint32_t>(lastCandidateCoord[1] + lastCandidateOffset[1]);
        lastCandidateCoord[2] = static_cast<uint32_t>(lastCandidateCoord[2] + lastCandidateOffset[2]);
        uint64_t lastCandidateMorton = PABLO::computeMorton(m_dim, lastCandidateCoord[0], lastCandidateCoord[1], lastCandidateCoord[2]);

        // Search for neighbours of different sizes
        if (candidateIdx < getNumOctants()){
            while(true){
                // Detect if the candidate is a neighbour
                u32array3 coordtry = {{0, 0, 0}};
                bool isNeighbourCandidate = true;
                for (int8_t idim=0; idim<m_dim; idim++){
                    coordtry[idim] = m_octants[candidateIdx].getLogicalCoordinates(idim);

                    int32_t Dx     = int32_t(int32_t(abs(cxyz[idim]))*(coordtry[idim] - coord[idim]));
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
                        for (int8_t idim=0; idim<m_dim; idim++){
                            coord1[idim] = coord[idim] + size;
                        }

                        if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
                            isNeighbour = true;
                        }
                    }
                    else if (leveltry < level){
                        u32array3 coordtry1 = {{1, 1, 1}};
                        for (int8_t idim=0; idim<m_dim; idim++){
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
                if (candidateIdx > getNumOctants() - 1){
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
        bool ghostSearch = !onlyinternal && (getNumGhosts() > 0);
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
            if (candidateIdx < getNumGhosts()){
                while(true){
                    // Detect if the candidate is a neighbour
                    u32array3 coordtry = {{0, 0, 0}};
                    bool isNeighbourCandidate = true;
                    for (int8_t idim=0; idim<m_dim; idim++){
                        coordtry[idim] = m_ghosts[candidateIdx].getLogicalCoordinates(idim);

                        int32_t Dx     = int32_t(int32_t(abs(cxyz[idim]))*(coordtry[idim] - coord[idim]));
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
                            for (int8_t idim=0; idim<m_dim; idim++){
                                coord1[idim] = coord[idim] + size;
                            }

                            if((abs(cxyz[0])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[1])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2]))) + (abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1])))){
                                isNeighbour = true;
                            }
                        }
                        else if (leveltry < level){
                            u32array3 coordtry1 = {{1, 1, 1}};
                            for (int8_t idim=0; idim<m_dim; idim++){
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
                    if (candidateIdx > getNumGhosts() - 1){
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
     * \param[in,out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[in,out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
     * \param[in] onlyinternal A boolean flag to specify if neighbours have to be found among all the octants (false) or only among the internal ones (true).
     * \param[in] append A boolean flag to specify if neighbours will be appended to the given vector or if the given vectors will be cleared before adding the neighbours.
     */
    void
    LocalTree::findEdgeNeighbours(const Octant* oct, uint8_t iedge, u32vector & neighbours, bvector & isghost, bool onlyinternal, bool append) const{

        if (!append) {
            isghost.clear();
            neighbours.clear();
        }

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

        std::array<uint32_t, 3> coord = oct->getLogicalCoordinates();
        if (isperiodic) {
            std::array<int64_t, 3> periodicOffset = getEdgePeriodicOffset(*oct, iedge);
            coord[0] = static_cast<uint32_t>(coord[0] + periodicOffset[0]);
            coord[1] = static_cast<uint32_t>(coord[1] + periodicOffset[1]);
            coord[2] = static_cast<uint32_t>(coord[2] + periodicOffset[2]);
        }

        const int8_t (&cxyz)[3] = m_treeConstants->edgeCoeffs[iedge];

        // Compute same-size virtual neighbour information
        std::array<int64_t, 3> sameSizeVirtualNeighOffset = computeFirstVirtualEdgeNeighOffset(level, iedge, level);
        std::array<uint32_t, 3> sameSizeVirtualNeighCoord = coord;
        sameSizeVirtualNeighCoord[0] = static_cast<uint32_t>(sameSizeVirtualNeighCoord[0] + sameSizeVirtualNeighOffset[0]);
        sameSizeVirtualNeighCoord[1] = static_cast<uint32_t>(sameSizeVirtualNeighCoord[1] + sameSizeVirtualNeighOffset[1]);
        sameSizeVirtualNeighCoord[2] = static_cast<uint32_t>(sameSizeVirtualNeighCoord[2] + sameSizeVirtualNeighOffset[2]);
        uint64_t sameSizeVirtualNeighMorton = PABLO::computeMorton(m_dim, sameSizeVirtualNeighCoord[0], sameSizeVirtualNeighCoord[1], sameSizeVirtualNeighCoord[2]);

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
        // This is the Morton number of the last descendant of the same-size
        // virtual neighbour.
        uint8_t maxEdgeNeighLevel = getMaxEdgeNeighLevel(*oct);
        std::array<int64_t, 3> lastCandidateOffset = computeLastVirtualEdgeNeighOffset(level, iedge, maxEdgeNeighLevel);
        std::array<uint32_t, 3> lastCandidateCoord = coord;
        lastCandidateCoord[0] = static_cast<uint32_t>(lastCandidateCoord[0] + lastCandidateOffset[0]);
        lastCandidateCoord[1] = static_cast<uint32_t>(lastCandidateCoord[1] + lastCandidateOffset[1]);
        lastCandidateCoord[2] = static_cast<uint32_t>(lastCandidateCoord[2] + lastCandidateOffset[2]);
        uint64_t lastCandidateMorton = PABLO::computeMorton(m_dim, lastCandidateCoord[0], lastCandidateCoord[1], lastCandidateCoord[2]);

        // Search for neighbours of different sizes
        if (candidateIdx < getNumOctants()) {
            while(true){
                // Detect if the candidate is a neighbour
                u32array3 coordtry = {{0, 0, 0}};
                bool isNeighbourCandidate = true;
                for (int8_t idim=0; idim<m_dim; idim++){
                    coordtry[idim] = m_octants[candidateIdx].getLogicalCoordinates(idim);

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
                        for (int8_t idim=0; idim<m_dim; idim++){
                            coord1[idim] = coord[idim] + size;
                        }

                        if((abs(cxyz[0])*abs(cxyz[2])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))) + (abs(cxyz[1])*abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))) + (abs(cxyz[0])*abs(cxyz[1])*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2])))){
                            isNeighbour = true;
                        }
                    }
                    else if (leveltry < level){
                        u32array3 coordtry1 = {{1, 1, 1}};
                        for (int8_t idim=0; idim<m_dim; idim++){
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
                if (candidateIdx > getNumOctants()-1){
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
        if (getNumGhosts() > 0 && !onlyinternal){
            // Identify the index of the first neighbour candidate
            computeNeighSearchBegin(sameSizeVirtualNeighMorton, m_ghosts, &candidateIdx, &candidateMorton);

            // Early return if a neighbour of the same size has been found
            if(candidateMorton == sameSizeVirtualNeighMorton && m_ghosts[candidateIdx].m_level == level){
                isghost.push_back(true);
                neighbours.push_back(candidateIdx);
                return;
            }

            // Search for neighbours of different sizes
            if (candidateIdx < getNumGhosts()){
                while(true){
                    // Detect if the candidate is a neighbour
                    u32array3 coordtry = {{0, 0, 0}};
                    bool isNeighbourCandidate = true;
                    for (int8_t idim=0; idim<m_dim; idim++){
                        coordtry[idim] = m_ghosts[candidateIdx].getLogicalCoordinates(idim);

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
                            for (int8_t idim=0; idim<m_dim; idim++){
                                coord1[idim] = int32_t(coord[idim]) + size;
                            }

                            if((abs(cxyz[0])*abs(cxyz[2])*((coordtry[1]>=coord[1])*(coordtry[1]<coord1[1]))) + (abs(cxyz[1])*abs(cxyz[2])*((coordtry[0]>=coord[0])*(coordtry[0]<coord1[0]))) + (abs(cxyz[0])*abs(cxyz[1])*((coordtry[2]>=coord[2])*(coordtry[2]<coord1[2])))){
                                isNeighbour = true;
                            }
                        }
                        else if (leveltry < level){
                            u32array3 coordtry1 = {{1, 1, 1}};
                            for (int8_t idim=0; idim<m_dim; idim++){
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
                    if (candidateIdx > getNumGhosts()-1){
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
     * \param[in,out] neighbours Vector of neighbours indices in octants/ghosts structure
     * \param[in,out] isghost Vector with boolean flag; true if the respective octant in neighbours is a ghost octant. Can be ignored in serial runs
    * \param[in] append A boolean flag to specify if neighbours will be appended to the given vector or if the given vectors will be cleared before adding the neighbours.*
     * \param[in] onlyinternal A boolean flag to specify if neighbours have to be found among all the octants (false) or only among the internal ones (true).
     */
    void
    LocalTree::findNodeNeighbours(const Octant* oct, uint8_t inode, u32vector & neighbours, bvector & isghost, bool onlyinternal, bool append) const{

        if (!append) {
            isghost.clear();
            neighbours.clear();
        }

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

        // Initialize search
        uint32_t candidateIdx    = 0;
        uint64_t candidateMorton = 0;

        uint32_t size = oct->getLogicalSize();

        std::array<uint32_t, 3> coord = oct->getLogicalCoordinates();
        if (isperiodic) {
            std::array<int64_t, 3> periodicOffset = getNodePeriodicOffset(*oct, inode);
            coord[0] = static_cast<uint32_t>(coord[0] + periodicOffset[0]);
            coord[1] = static_cast<uint32_t>(coord[1] + periodicOffset[1]);
            coord[2] = static_cast<uint32_t>(coord[2] + periodicOffset[2]);
        }

        const int8_t (&cxyz)[3] = m_treeConstants->nodeCoeffs[inode];

        // Compute same-size virtual neighbour information
        std::array<int64_t, 3> sameSizeVirtualNeighOffset = computeFirstVirtualNodeNeighOffset(level, inode, level);
        std::array<uint32_t, 3> sameSizeVirtualNeighCoord = coord;
        sameSizeVirtualNeighCoord[0] = static_cast<uint32_t>(sameSizeVirtualNeighCoord[0] + sameSizeVirtualNeighOffset[0]);
        sameSizeVirtualNeighCoord[1] = static_cast<uint32_t>(sameSizeVirtualNeighCoord[1] + sameSizeVirtualNeighOffset[1]);
        sameSizeVirtualNeighCoord[2] = static_cast<uint32_t>(sameSizeVirtualNeighCoord[2] + sameSizeVirtualNeighOffset[2]);
        uint64_t sameSizeVirtualNeighMorton = PABLO::computeMorton(m_dim, sameSizeVirtualNeighCoord[0], sameSizeVirtualNeighCoord[1], sameSizeVirtualNeighCoord[2]);

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
        // This is the Morton number of the last descendant of the same-size
        // virtual neighbour.
        uint8_t maxNodeNeighLevel = getMaxNodeNeighLevel(*oct);
        std::array<int64_t, 3> lastCandidateOffset = computeLastVirtualNodeNeighOffset(level, inode, maxNodeNeighLevel);
        std::array<uint32_t, 3> lastCandidateCoord = coord;
        lastCandidateCoord[0] = static_cast<uint32_t>(lastCandidateCoord[0] + lastCandidateOffset[0]);
        lastCandidateCoord[1] = static_cast<uint32_t>(lastCandidateCoord[1] + lastCandidateOffset[1]);
        lastCandidateCoord[2] = static_cast<uint32_t>(lastCandidateCoord[2] + lastCandidateOffset[2]);
        uint64_t lastCandidateMorton = PABLO::computeMorton(m_dim, lastCandidateCoord[0], lastCandidateCoord[1], lastCandidateCoord[2]);

        // Search for neighbours of different sizes
        if (candidateIdx < getNumOctants()) {
            while(true){
                // Detect if the candidate is a neighbour
                u32array3 coordtry = {{0, 0, 0}};
                bool isNeighbour = true;
                for (int8_t idim=0; idim<m_dim; idim++){
                    coordtry[idim] = m_octants[candidateIdx].getLogicalCoordinates(idim);

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
                if (candidateIdx > getNumOctants()-1){
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

        if (getNumGhosts() > 0 && !onlyinternal){
            // Identify the index of the first neighbour candidate
            computeNeighSearchBegin(sameSizeVirtualNeighMorton, m_ghosts, &candidateIdx, &candidateMorton);

            // Early return if a neighbour of the same size has been found
            if(candidateMorton == sameSizeVirtualNeighMorton && m_ghosts[candidateIdx].m_level == oct->m_level){
                isghost.push_back(true);
                neighbours.push_back(candidateIdx);
                return;
            }

            // Search for neighbours of different sizes
            if (candidateIdx < getNumGhosts()) {
                while(true){
                    // Detect if the candidate is a neighbour
                    u32array3 coordtry = {{0, 0, 0}};
                    bool isNeighbour = true;
                    for (int8_t idim=0; idim<m_dim; idim++){
                        coordtry[idim] = m_ghosts[candidateIdx].getLogicalCoordinates(idim);

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
                    if (candidateIdx > getNumGhosts()-1){
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

    /*! Fix markers of broken families over processes.
     * \param[out] updatedOctants If a valid pointer is provided, the pointers of the updated
     * octants will be added to the specified list.
     * \param[out] updatedGhostFlags If a valid pointer is provided, the ghost flags of the
     * updated octants will be added to the specified list.
     * \return True if some markers were modified to fix broken families.
     */
    bool
    LocalTree::fixBrokenFamiliesMarkers(std::vector<Octant *> *updatedOctants, std::vector<bool> *updatedGhostFlags){

        // Initialization
        bool updated = false;

        uint32_t nocts = m_octants.size();
        if (nocts == 0) {
            return updated;
        }

        uint32_t nghosts = m_ghosts.size();

        uint32_t idx      = 0;
        uint32_t idx0     = 0;
        uint32_t idx2     = 0;
        uint32_t idx1_gh  = 0;
        uint32_t idx2_gh  = 0;
        uint32_t last_idx =nocts-1;

        uint8_t nbro = 0;
        Octant father(m_dim);

        //Clean index of ghost brothers in case of coarsening a broken family
        m_lastGhostBros.clear();
        m_firstGhostBros.clear();

        // Process ghosts
        bool checkend = true;
        bool checkstart = true;
        if (nghosts > 0){
            // Set index for start and end check for ghosts
            while(m_ghosts[idx2_gh].getMorton() <= m_lastDescMorton){
                idx2_gh++;
                if (idx2_gh > nghosts-1) break;
            }
            if (idx2_gh > nghosts-1) checkend = false;

            while(m_ghosts[idx1_gh].getMorton() <= m_octants[0].getMorton()){
                idx1_gh++;
                if (idx1_gh > nghosts-1) break;
            }
            if (idx1_gh == 0) checkstart = false;
            idx1_gh-=1;

            // Start and End on ghosts
            if (checkstart){
                if (m_ghosts[idx1_gh].buildFather()==m_octants[0].buildFather()){
                    father = m_ghosts[idx1_gh].buildFather();
                    nbro = 0;
                    idx = idx1_gh;
                    int8_t marker = m_ghosts[idx].getMarker();
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
                            Octant &octant = m_octants[ii];
                            if (octant.getMarker()<0){
                                octant.setMarker(0);
                                if (updatedOctants) {
                                    updatedOctants->push_back(&octant);
                                }
                                if (updatedGhostFlags) {
                                    updatedGhostFlags->push_back(false);
                                }
                                updated = true;
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
                    int8_t marker = m_ghosts[idx].getMarker();
                    while(marker < 0 && m_ghosts[idx].buildFather() == father){

                        //Add ghost index to structure for mapper in case of coarsening a broken family
                        m_lastGhostBros.push_back(idx);

                        nbro++;
                        idx++;
                        if(idx == nghosts){
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
                            Octant &octant = m_octants[ii];
                            if (octant.getMarker()<0){
                                octant.setMarker(0);
                                if (updatedOctants) {
                                    updatedOctants->push_back(&octant);
                                }
                                if (updatedGhostFlags) {
                                    updatedGhostFlags->push_back(false);
                                }
                                updated = true;
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
        uint64_t mortonld = father.computeLastDescMorton();
        nbro = 0;
        for (idx=0; idx<m_treeConstants->nChildren; idx++){
            // Check if family is complete or to be checked in the internal loop (some brother refined)
            if (idx<nocts){
                if (m_octants[idx].getMorton() <= mortonld){
                    nbro++;
                }
            }
        }

        if (nbro != m_treeConstants->nChildren) {
            idx0 = nbro;
        }

        // Check and coarse internal octants
        for (idx=idx0; idx<nocts; idx++){
            Octant &octant = m_octants[idx];
            if(octant.getMarker() < 0 && octant.getLevel() > 0){
                nbro = 0;
                father = octant.buildFather();
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
                        octant.setMarker(0);
                        if (updatedOctants) {
                            updatedOctants->push_back(&octant);
                        }
                        if (updatedGhostFlags) {
                            updatedGhostFlags->push_back(false);
                        }
                        updated = true;
                    }
                }
            }
        }

        return updated;
    }

    // =================================================================================== //

    /*! 2:1 balancing on level a local tree (refinement wins!)
     * \param[in] doNew Set to true the balance is enforced also on new octants.
     * \param[in] checkInterior Set to true if interior octants should be checked.
     * \param[in] checkGhost Set to true if ghost octants should be checked.
     * \return True if balanced done with some markers modification.
     */
    bool
    LocalTree::localBalance(bool doNew, bool checkInterior, bool checkGhost){

        bool balanceEdges = ((m_balanceCodim>1) && (m_dim==3));
        bool balanceNodes = (m_balanceCodim==m_dim);

        std::vector<Octant *> processOctants;
        std::vector<bool> processGhostFlags;

        // Identify internal octants that will be processed
        if(checkInterior){
            for (Octant &octant : m_octants){
                // Skip octants that doesn't need balancing
                bool balanceOctant = octant.getBalance();
                if (balanceOctant) {
                    if (doNew) {
                        balanceOctant = (octant.getMarker() != 0) || octant.getIsNewC() || octant.getIsNewR();
                    } else {
                        balanceOctant = (octant.getMarker() != 0);
                    }
                }

                if (!balanceOctant) {
                    continue;
                }

                // Add octant to the process list
                processOctants.push_back(&octant);
                processGhostFlags.push_back(false);
            }
        }

        // Identify ghost octants that will be processed
        //
        // Ghost octants will be balanced by the process that owns them, however they may
        // affect balancing of local octants. If ghost octnts are processed, we process all
        // ghost octants of the first layer, not only the ones that need balancing. That's
        // because it's faster to process all the ghost octants that may affect balacning
        // of internal octants, rather than find the ones that actually affect balacing of
        // internal octants.
        if (checkGhost) {
            for (Octant &octant : m_ghosts){
                // Only ghosts of the first layer can affect load balance
                if (octant.getGhostLayer() > 0) {
                    continue;
                }

                // Add octant to the process list
                processOctants.push_back(&octant);
                processGhostFlags.push_back(true);
            }
        }

        // Iterative balacing
        std::vector<uint32_t> neighs;
        std::vector<bool> neighGhostFlags;

        bool updated = false;
        std::vector<Octant *> updatedProcessOctants;
        std::vector<bool> updatedProcessGhostFlags;
        while (!processOctants.empty()) {
            // Fix broken families markers
            bool brokenFamilieisFixed = fixBrokenFamiliesMarkers(&processOctants, &processGhostFlags);
            if (brokenFamilieisFixed) {
                updated = true;
            }

            // Balance
            std::size_t processSize = processOctants.size();
            for (std::size_t n = 0; n < processSize; ++n) {
                Octant &octant = *(processOctants[n]);
                bool ghostFlag = processGhostFlags[n];

                // Identify neighbours that will be used for balancing
                neighs.clear();
                neighGhostFlags.clear();

                for (int iface=0; iface < m_treeConstants->nFaces; iface++) {
                    findNeighbours(&octant, iface, neighs, neighGhostFlags, ghostFlag, true);
                }

                if (balanceNodes) {
                    for (int inode=0; inode < m_treeConstants->nNodes; inode++) {
                        findNodeNeighbours(&octant, inode, neighs, neighGhostFlags, ghostFlag, true);
                    }
                }

                if (balanceEdges) {
                    for (int iedge=0; iedge < m_treeConstants->nEdges; iedge++) {
                        findEdgeNeighbours(&octant, iedge, neighs, neighGhostFlags, ghostFlag, true);
                    }
                }

                // Balance octant
                int8_t level       = octant.getLevel();
                int8_t futureLevel = std::min(m_treeConstants->maxLevel, static_cast<int8_t>(level + octant.getMarker()));

                std::size_t nNeighs = neighs.size();
                for(std::size_t i = 0; i < nNeighs; i++){
                    bool neighGhostFlag = neighGhostFlags[i];
                    Octant &neighOctant = (neighGhostFlag ? m_ghosts[neighs[i]]: m_octants[neighs[i]]);

                    int8_t neighLevel       = neighOctant.getLevel();
                    int8_t neighFutureLevel = std::min(m_treeConstants->maxLevel, int8_t(neighLevel + neighOctant.getMarker()));
                    if(!ghostFlag && futureLevel < neighFutureLevel - 1) {
                        futureLevel = neighFutureLevel - 1;
                        octant.setMarker(futureLevel - level);
                        updatedProcessOctants.push_back(&octant);
                        updatedProcessGhostFlags.push_back(ghostFlag);
                        updated = true;
                    }
                    else if(!neighGhostFlag && neighFutureLevel < futureLevel - 1) {
                        neighOctant.setMarker((futureLevel - 1) - neighLevel);
                        if (neighOctant.getBalance()) {
                            updatedProcessOctants.push_back(&neighOctant);
                            updatedProcessGhostFlags.push_back(neighGhostFlag);
                        }
                        updated = true;
                    }
                }
            }

            // Update process list
            updatedProcessOctants.swap(processOctants);
            updatedProcessGhostFlags.swap(processGhostFlags);

            updatedProcessOctants.clear();
            updatedProcessGhostFlags.clear();
        }

        return updated;
    };

	// =================================================================================== //

	/*! Compute the offset from the origin of an octant to its first virtual face neighobur.
	* The term "virtual" means that the tree may not contain any octant at the position
	* defined by the computed offset.
	* \param[in] level The level of the octant.
	* \param[in] iface Local index of the face.
	* \param[in] neighLevel The level of the virtual neighbours.
	* \return The computed offset.
	*/
	std::array<int64_t, 3> LocalTree::computeFirstVirtualNeighOffset(uint8_t level, uint8_t iface, uint8_t neighLevel) const {

		// Get normalized offsets
		const int8_t (&normalizedOffsets)[3] = m_treeConstants->normals[iface];

		// Get octant sizes
		uint32_t size      = m_treeConstants->lengths[level];
		uint32_t neighSize = m_treeConstants->lengths[neighLevel];

		// Compute the coordinates of the virtual neighbour
		std::array<int64_t, 3> neighOffsets;
		for (int i = 0; i < 3; ++i) {
			neighOffsets[i] = normalizedOffsets[i];
			if (neighOffsets[i] > 0) {
				neighOffsets[i] *= size;
			} else {
				neighOffsets[i] *= neighSize;
			}
		}

		return neighOffsets;
	}

	/*! Compute the offset from the origin of an octant to its last virtual face neighobur.
	* The term "virtual" means that the tree may not contain any octant at the position
	* defined by the computed offset.
	* \param[in] level The level of the octant.
	* \param[in] iface Local index of the face.
	* \param[in] neighLevel The level of the virtual neighbours.
	* \return The computed offset.
	*/
	std::array<int64_t, 3> LocalTree::computeLastVirtualNeighOffset(uint8_t level, uint8_t iface, uint8_t neighLevel) const {

		// Get octant sizes
		uint32_t size      = m_treeConstants->lengths[level];
		uint32_t neighSize = m_treeConstants->lengths[neighLevel];

		// Compute the coordinates of the virtual neighbour
		std::array<int64_t, 3> neighOffsets = computeFirstVirtualNeighOffset(level, iface, neighLevel);

		uint32_t delta = size - neighSize;
		switch (iface) {
			case 0 :
			case 1 :
				neighOffsets[1] += delta;
				neighOffsets[2] += delta;
				break;

			case 2 :
			case 3 :
				neighOffsets[0] += delta;
				neighOffsets[2] += delta;
				break;

			case 4 :
			case 5 :
				neighOffsets[0] += delta;
				neighOffsets[1] += delta;
				break;
		}

		return neighOffsets;
	}

	/*! Compute the offsets from the origin of an octant to its virtual node neighoburs.
	* The term "virtual" means that the tree may not contain any octant at the positions
	* defined by the computed offsets.
	* \param[in] level The level of the octant.
	* \param[in] iface Local index of the node.
	* \param[in] neighLevel The level of the virtual neighbours.
	* \param[out] neighOffsets On output will contain the offsets.
	*/
	void LocalTree::computeVirtualNeighOffsets(uint8_t level, uint8_t iface, uint8_t neighLevel, std::vector<std::array<int64_t, 3>> *neighOffsets) const {

		// Get octant sizes
		uint32_t neighSize = m_treeConstants->lengths[neighLevel];

		// Direction of the increments
		int direction1;
		int direction2;
		switch (iface) {
			case 0:
			case 1:
				direction1 = 1;
				direction2 = 2;
				break;

			case 2:
			case 3:
				direction1 = 0;
				direction2 = 2;
				break;

			default:
				direction1 = 0;
				direction2 = 1;
				break;
		}

		// Compute the coordinates of the virtual neighbour
		int nNeighsDirection1 = 1 << (neighLevel - level);
		int nNeighsDirection2 = (m_dim > 2) ? nNeighsDirection1 : 1;
		int nNeighs           = nNeighsDirection1 * nNeighsDirection2;

		neighOffsets->assign(nNeighs, computeFirstVirtualNeighOffset(level, iface, neighLevel));

		int k = 0;
		for (int j = 0; j < nNeighsDirection2; ++j) {
			for (int i = 0; i < nNeighsDirection1; ++i) {
				(*neighOffsets)[k][direction1] += i * neighSize;
				(*neighOffsets)[k][direction2] += j * neighSize;
				++k;
			}
		}
	}

	/*! Compute the offset from the origin of an octant to its first virtual node neighobur.
	* The term "virtual" means that the tree may not contain any octant at the position
	* defined by the computed offset.
	* \param[in] level The level of the octant.
	* \param[in] inode Local index of the node.
	* \param[in] neighLevel The level of the virtual neighbours.
	* \return The computed offset.
	*/
	std::array<int64_t, 3> LocalTree::computeFirstVirtualNodeNeighOffset(uint8_t level, uint8_t inode, uint8_t neighLevel) const {

		// Get normalized offsets
		const int8_t (&normalizedOffsets)[3] = m_treeConstants->nodeCoeffs[inode];

		// Get octant sizes
		uint32_t size      = m_treeConstants->lengths[level];
		uint32_t neighSize = m_treeConstants->lengths[neighLevel];

		// Compute the coordinates of the virtual neighbour
		std::array<int64_t, 3> neighOffsets;
		for (int i = 0; i < 3; ++i) {
			neighOffsets[i] = normalizedOffsets[i];
			if (neighOffsets[i] > 0) {
				neighOffsets[i] *= size;
			} else {
				neighOffsets[i] *= neighSize;
			}
		}

		return neighOffsets;
	}

	/*! Compute the offset from the origin of an octant to its last virtual node neighobur.
	* The term "virtual" means that the tree may not contain any octant at the position
	* defined by the computed offset.
	* \param[in] level The level of the octant.
	* \param[in] inode Local index of the node.
	* \param[in] neighLevel The level of the virtual neighbours.
	* \return The computed offset.
	*/
	std::array<int64_t, 3> LocalTree::computeLastVirtualNodeNeighOffset(uint8_t level, uint8_t inode, uint8_t neighLevel) const {

		return computeFirstVirtualNodeNeighOffset(level, inode, neighLevel);
	}

	/*! Compute the offsets from the origin of an octant to its virtual node neighoburs.
	* The term "virtual" means that the tree may not contain any octant at the positions
	* defined by the computed offsets.
	* \param[in] level The level of the octant.
	* \param[in] inode Local index of the node.
	* \param[in] neighLevel The level of the virtual neighbours.
	* \param[out] neighOffsets On output will contain the offsets.
	*/
	void LocalTree::computeVirtualNodeNeighOffsets(uint8_t level, uint8_t inode, uint8_t neighLevel, std::vector<std::array<int64_t, 3>> *neighOffsets) const {

		neighOffsets->assign(1, computeFirstVirtualNodeNeighOffset(level, inode, neighLevel));
	}

	/*! Compute the offset from the origin of an octant to its first virtual edge neighobur.
	* The term "virtual" means that the tree may not contain any octant at the position
	* defined by the computed offset.
	* \param[in] level The level of the octant.
	* \param[in] iedge Local index of the edge.
	* \param[in] neighLevel The level of the virtual neighbours.
	* \return The computed offset.
	*/
	std::array<int64_t, 3> LocalTree::computeFirstVirtualEdgeNeighOffset(uint8_t level, uint8_t iedge, uint8_t neighLevel) const {

		// Get normalized offsets
		const int8_t (&normalizedOffsets)[3] = m_treeConstants->edgeCoeffs[iedge];

		// Get octant sizes
		uint32_t size      = m_treeConstants->lengths[level];
		uint32_t neighSize = m_treeConstants->lengths[neighLevel];

		// Compute the coordinates of the virtual neighbour
		std::array<int64_t, 3> neighOffsets;
		for (int i = 0; i < 3; ++i) {
			neighOffsets[i] = normalizedOffsets[i];
			if (neighOffsets[i] > 0) {
				neighOffsets[i] *= size;
			} else {
				neighOffsets[i] *= neighSize;
			}
		}

		return neighOffsets;
	}

	/*! Compute the offset from the origin of an octant to its last virtual edge neighobur.
	* The term "virtual" means that the tree may not contain any octant at the position
	* defined by the computed offset.
	* \param[in] level The level of the octant.
	* \param[in] iedge Local index of the edge.
	* \param[in] neighLevel The level of the virtual neighbours.
	* \return The computed offset.
	*/
	std::array<int64_t, 3> LocalTree::computeLastVirtualEdgeNeighOffset(uint8_t level, uint8_t iedge, uint8_t neighLevel) const {

		// Get octant sizes
		uint32_t size      = m_treeConstants->lengths[level];
		uint32_t neighSize = m_treeConstants->lengths[neighLevel];

		// Compute the coordinates of the virtual neighbour
		std::array<int64_t, 3> neighOffsets = computeFirstVirtualEdgeNeighOffset(level, iedge, neighLevel);

		uint32_t offset = size - neighSize;
		switch (iedge) {
			case 0:
			case 1:
			case 8:
			case 9:
				neighOffsets[1] += offset;
				break;

			case 2:
			case 3:
			case 10:
			case 11:
				neighOffsets[0] += offset;
				break;

			case 4:
			case 5:
			case 6:
			case 7:
				neighOffsets[2] += offset;
				break;
		}

		return neighOffsets;
	}

	/*! Compute the offsets from the origin of an octant to its virtual edge neighoburs.
	* The term "virtual" means that the tree may not contain any octant at the positions
	* defined by the computed offsets.
	* \param[in] level The level of the octant.
	* \param[in] iedge Local index of the edge.
	* \param[in] neighLevel The level of the virtual neighbours.
	* \param[out] neighOffsets On output will contain the offsets.
	*/
	void LocalTree::computeVirtualEdgeNeighOffsets(uint8_t level, uint8_t iedge, uint8_t neighLevel, std::vector<std::array<int64_t, 3>> *neighOffsets) const {

		// Get octant sizes
		uint32_t neighSize = m_treeConstants->lengths[neighLevel];

		// Direction of the increments
		int direction;
		switch (iedge) {
			case 0:
			case 1:
			case 8:
			case 9:
				direction = 1;
				break;

			case 2:
			case 3:
			case 10:
			case 11:
				direction = 0;
				break;

			default:
				direction = 2;
				break;
		}

		// Compute the coordinates of the virtual neighbours
		int nNeighs = 1 << (neighLevel - level);
		neighOffsets->assign(nNeighs, computeFirstVirtualEdgeNeighOffset(level, iedge, neighLevel));
		for (int i = 1; i < nNeighs; ++i) {
			(*neighOffsets)[i][direction] += i * neighSize;
		}
	}

	// =================================================================================== //

	/*! Get the periodic offset for the specified face of the given octant.
	* The offset is the displacement that should be applied to a point on the specified
	* face to move that point to the corresponding periodic face.
	* \param[in] octant Octant for which the offset should be evaluated
	* \param[in] iface Index of face
	* \return The periodic offset for the specified face of the given octant.
	*/
	std::array<int64_t, 3> LocalTree::getPeriodicOffset(const Octant &octant, uint8_t iface) const {

		std::array<int64_t, 3> offset = {{0, 0, 0}};
		if (!isPeriodic(&octant, iface)) {
			return offset;
		}

		int64_t maxLength = m_treeConstants->lengths[0];
		const int8_t (&referenceOffset)[3] = m_treeConstants->normals[iface];

		for (int d = 0; d < m_dim; ++d) {
			offset[d] = - maxLength * static_cast<int64_t>(referenceOffset[d]);
		}

		return offset;
	}

	/*! Get the periodic offset for the specified node of the given octant.
	* The offset is the displacement that should be applied to a point on the specified
	* node to move that point to the corresponding periodic node.
	* \param[in] octant Octant for which the offset should be evaluated
	* \param[in] inode Index of node
	* \return The periodic offset for the specified node of the given octant.
	*/
	std::array<int64_t, 3> LocalTree::getNodePeriodicOffset(const Octant &octant, uint8_t inode) const {

		std::array<int64_t, 3> offset = {{0, 0, 0}};
		for (int i = 0; i < m_dim; ++i) {
			uint8_t face = m_treeConstants->nodeFace[inode][i];
			if (!isPeriodic(&octant, face)) {
				continue;
			}

			std::array<int64_t, 3> faceOffset = getPeriodicOffset(octant, face);
			for (int d = 0; d < m_dim; ++d) {
				offset[d] += faceOffset[d];
			}
		}

		return offset;
	}

	/*! Get the periodic offset for the specified edge of the given octant.
	* The offset is the displacement that should be applied to a point on the specified
	* edge to move that point to the corresponding periodic edge.
	* \param[in] octant Octant for which the offset should be evaluated
	* \param[in] iedge Index of edge
	* \return The periodic offset for the specified edge of the given octant.
	*/
	std::array<int64_t, 3> LocalTree::getEdgePeriodicOffset(const Octant &octant, uint8_t iedge) const {

		std::array<int64_t, 3> offset = {{0, 0, 0}};
		for (uint8_t face : m_treeConstants->edgeFace[iedge]) {
			if (!isPeriodic(&octant, face)) {
				continue;
			}

			std::array<int64_t, 3> faceOffset = getPeriodicOffset(octant, face);
			for (int d = 0; d < m_dim; ++d) {
				offset[d] += faceOffset[d];
			}
		}

		return offset;
	}

    // =================================================================================== //

	/*! Get the maximum level a face neighbour of the specified octant can have.
	* If 2:1 face balance is active, a neighbour can only have one more refinement level than
	* the specified octant.
	* \param[in] octant Octant for which the check has to be performed.
	* \return The maximum level a face neighbour of the specified octant can have.
	*/
	uint8_t LocalTree::getMaxNeighLevel(const Octant &octant) const {

		if (octant.getBalance()) {
			return octant.getLevel() + 1;
		} else {
			return m_treeConstants->maxLevel;
		}
	}

	/*! Get the maximum level a node neighbour of the specified octant can have.
	* If 2:1 node balance is active, a neighbour can only have one more refinement level than
	* the specified octant.
	* \param[in] octant Octant for which the check has to be performed.
	* \return The maximum level a node neighbour of the specified octant can have.
	*/
	uint8_t LocalTree::getMaxNodeNeighLevel(const Octant &octant) const {

        if (octant.getBalance() && m_balanceCodim == m_dim) {
			return octant.getLevel() + 1;
		} else {
			return m_treeConstants->maxLevel;
		}
	}

	/*! Get the maximum level an edge neighbour of the specified octant can have.
	* If 2:1 edge balance is active, a neighbour can only have one more refinement level than
	* the specified octant.
	* \param[in] octant Octant for which the check has to be performed.
	* \return The maximum level an edge neighbour of the specified octant can have.
	*/
	uint8_t LocalTree::getMaxEdgeNeighLevel(const Octant &octant) const {

		if (octant.getBalance() && m_balanceCodim > 1) {
			return octant.getLevel() + 1;
		} else {
			return m_treeConstants->maxLevel;
		}
	}

    // =================================================================================== //

    /*! Compute and store in m_intersections the intersections of the local tree.
     */
    void
    LocalTree::computeIntersections() {

		octvector::iterator 	it, obegin, oend;
		Intersection 			intersection;
		u32vector 				neighbours;
		vector<bool>			isghost;
		uint32_t 				idx;
		uint32_t 				i, nsize;
		uint8_t 				iface, iface2;

		m_intersections.clear();
		m_intersections.reserve(2*3*m_octants.size());

		idx = 0;

		// Loop on ghosts
		obegin = m_ghosts.begin();
		oend = m_ghosts.end();
		for (it = obegin; it != oend; ++it){
			for (iface = 0; iface < m_dim; iface++){
				iface2 = iface*2;
				findNeighbours(m_ghosts.data() + idx, iface2, neighbours, isghost, true, false);
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
				findNeighbours(m_octants.data() + idx, iface2, neighbours, isghost, false, false);
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
					}
					else{
						//Periodic intersection
						findNeighbours(m_octants.data() + idx, iface2+1, neighbours, isghost, false, false);
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
        vector<uint64_t>                             nodeKeys;
        unordered_map<uint64_t, array<uint32_t, 3> > nodeCoords;
        unordered_map<uint64_t, vector<uint64_t> >   nodeOctants;
        uint32_t                                     noctants = getNumOctants();
        uint32_t                                     nghosts  = getNumGhosts();



        // Gather node information
        nodeKeys.reserve(noctants);
        nodeCoords.reserve(noctants);
        nodeOctants.reserve(noctants);

        for (uint64_t n = 0; n < (noctants + nghosts); n++){
            const Octant *octant;
            if (n < noctants) {
                uint32_t octantId = static_cast<uint32_t>(n);
                octant = &(m_octants[octantId]);
            } else {
                uint32_t octantId = static_cast<uint32_t>(n - noctants);
                octant = &(m_ghosts[octantId]);
            }

            for (uint8_t i = 0; i < m_treeConstants->nNodes; ++i){
                u32array3 node;
                octant->getLogicalNode(node, i);

                uint64_t nodeKey = octant->computeNodePersistentKey(node);
                if (nodeCoords.count(nodeKey) == 0) {
                    nodeKeys.push_back(nodeKey);
                    nodeCoords.insert({{nodeKey, std::move(node)}});
                    nodeOctants[nodeKey].reserve(m_treeConstants->nNodes);
                }

                nodeOctants[nodeKey].push_back(n);
            }
        }
        std::sort(nodeKeys.begin(), nodeKeys.end());

        // Build node list and connectivity
        m_nodes.clear();
        m_connectivity.clear();
        m_ghostsConnectivity.clear();

        m_nodes.reserve(nodeKeys.size());
        m_connectivity.resize(noctants);
        m_ghostsConnectivity.resize(nghosts);

        uint32_t nodeId = 0;
        for (uint64_t nodeKey : nodeKeys) {
            m_nodes.emplace_back(std::move(nodeCoords.at(nodeKey)));
            for (uint64_t n : nodeOctants.at(nodeKey)) {
                std::vector<uint32_t> *octantConnect;
                if (n < noctants) {
                    uint32_t octantId = static_cast<uint32_t>(n);
                    octantConnect = &(m_connectivity[octantId]);
                } else {
                    uint32_t octantId = static_cast<uint32_t>(n - noctants);
                    octantConnect = &(m_ghostsConnectivity[octantId]);
                }

                if (octantConnect->size() == 0) {
                    octantConnect->reserve(m_treeConstants->nNodes);
                }
                octantConnect->push_back(nodeId);
            }
            nodeId++;
        }

        m_nodes.shrink_to_fit();
    };

    /*! Clear nodes vector and connectivity of octants of local tree
     * \param[in] release if it's true the memory hold by the connectivity will be
     * released, otherwise the connectivity will be cleared but its memory will
     * not be released
     */
    void
    LocalTree::clearConnectivity(bool release){
        if (release) {
            u32arr3vector().swap(m_nodes);
            u32vector2D().swap(m_connectivity);
            u32vector2D().swap(m_ghostsConnectivity);
        } else {
            m_nodes.clear();
            m_connectivity.clear();
            m_ghostsConnectivity.clear();
        }
    };

    /*! Updates nodes vector and connectivity of octants of local tree
     */
    void
    LocalTree::updateConnectivity(){
        clearConnectivity(false);
        computeConnectivity();
    };

    // =================================================================================== //

}
