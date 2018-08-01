/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

# ifndef __BITPIT_LEVELSET_CACHED_HPP__
# define __BITPIT_LEVELSET_CACHED_HPP__

// Standard Template Library
# include <array>
# include <vector>
# include <unordered_set>

# include "bitpit_containers.hpp"

namespace bitpit{

namespace adaption{
    struct Info;
}

class SendBuffer;
class RecvBuffer;

class LevelSetKernel;
class LevelSetObject;

class LevelSetCachedObject : public LevelSetObject{

    private:
    static const int                            PROPAGATION_STATUS_EXTERNAL;
    static const int                            PROPAGATION_STATUS_WAITING;
    static const int                            PROPAGATION_STATUS_REACHED;

    static const int                            PROPAGATION_SIGN_UNDEFINED;
    static const int                            PROPAGATION_SIGN_DUMMY;

    void                                        setSign( long id, int sign ) ;
    void                                        propagateSeedSign( const std::vector<long> &seeds,
                                                                   PiercedStorage<int, long> *status,
                                                                   long *nWaiting, int *externalSign ) ;

    protected:
    PiercedVector<LevelSetInfo>                 m_ls ;          /**< Levelset information for each cell */
    virtual void                                getBoundingBox( std::array<double,3> &, std::array<double,3> & )const =0  ;

    void                                        _clear( ) override ;
    virtual void                                __clear() ;


    void                                        _clearAfterMeshAdaption(const std::vector<adaption::Info> & ) override ;
    virtual void                                __clearAfterMeshAdaption(const std::vector<adaption::Info> & ) ;

    void                                        _dump( std::ostream &) override ;
    virtual void                                __dump(std::ostream &) ; 
    void                                        _restore( std::istream &) override ;
    virtual void                                __restore(std::istream &) ; 

# if BITPIT_ENABLE_MPI
    void                                        _writeCommunicationBuffer( const std::vector<long> &, SendBuffer & ) override ;
    virtual void                                __writeCommunicationBuffer( const std::vector<long> &, SendBuffer & ) ; 
    void                                        _readCommunicationBuffer( const std::vector<long> &, RecvBuffer & )  override;
    virtual void                                __readCommunicationBuffer( const std::vector<long> &, RecvBuffer & ) ; 
# endif 

    public:
    virtual ~LevelSetCachedObject();
    LevelSetCachedObject(int);

    LevelSetInfo                                getLevelSetInfo(const long &) const override ;
    double                                      getLS(const long &) const override ;
    std::array<double,3>                        getGradient(const long &) const override ;

    void                                        propagateSign() override ;

};



}

#endif 
