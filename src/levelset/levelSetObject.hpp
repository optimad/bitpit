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

# ifndef __BITPIT_LEVELSET_OBJECT_HPP__
# define __BITPIT_LEVELSET_OBJECT_HPP__

// Standard Template Library
# include <iostream>
# include <array>
# include <vector>
# include <unordered_map>

# include "levelSetCommon.hpp"

namespace bitpit{

namespace adaption{
    class Info;
}
class SendBuffer;
class RecvBuffer;

class LevelSetKernel;

class LevelSetObject{

    private:
    int                                         m_id;           /**< identifier of object */
    bool                                        m_primary;      /**< identifier of object */


    protected:
    double                                      m_RSearch;      /**< Size of narrow band */

    virtual void                                _clear();
    virtual void                                _clearAfterMeshAdaption(const std::vector<adaption::Info>&) ;
    virtual void                                _filterOutsideNarrowBand(double) ;
    virtual void                                _dump(std::ostream &);
    virtual void                                _restore( std::istream &);

# if BITPIT_ENABLE_MPI
    virtual void                                writeCommunicationBuffer( const std::vector<long> &, SendBuffer & ) ; 
    virtual void                                _writeCommunicationBuffer(const std::vector<long>&, SendBuffer&)  ;
    virtual void                                readCommunicationBuffer( const std::vector<long> &, RecvBuffer & ) ; 
    virtual void                                _readCommunicationBuffer(const std::vector<long>&, RecvBuffer&)  ;
# endif 

    public:
    virtual ~LevelSetObject();
    LevelSetObject(int,bool);

    virtual LevelSetObject*                     clone() const =0;
    void                                        clear();

    int                                         getId() const ;
    bool                                        isPrimary() const ;

    virtual LevelSetInfo                        getLevelSetInfo(const long &) const =0;
    virtual double                              getLS(const long &) const =0; 
    virtual std::array<double,3>                getGradient(const long &) const =0 ; 

    virtual int                                 getPart(const long &) const ;
    virtual long                                getSupport(const long &) const;
    virtual int                                 getSupportCount(const long &) const ;

    short                                       getSign(const long &) const;
    virtual void                                propagateSign(LevelSetKernel*) ; 

    bool                                        isInNarrowBand(const long &) const;
    double                                      getSizeNarrowBand() const;
    virtual void                                setSizeNarrowBand(double) ;

    virtual double                              computeSizeNarrowBand(LevelSetKernel*)=0;
    virtual void                                computeLSInNarrowBand(LevelSetKernel*, const double &, const bool &) ;

    virtual double                              updateSizeNarrowBand(LevelSetKernel*, const std::vector<adaption::Info> &);
    virtual void                                updateLSInNarrowBand(LevelSetKernel*, const std::vector<adaption::Info> &, const double &, const bool &) ;
    void                                        clearAfterMeshAdaption(const std::vector<adaption::Info>&);
    void                                        filterOutsideNarrowBand(double);  ;

    void                                        dump(std::ostream &); 
    void                                        restore(std::istream &); 

# if BITPIT_ENABLE_MPI
    bool                                        assureMPI(LevelSetKernel * ) ;
    void                                        exchangeGhosts(LevelSetKernel * ) ;
    void                                        communicate( LevelSetKernel *, std::unordered_map<int,std::vector<long>> &, std::unordered_map<int,std::vector<long>> &, std::vector<adaption::Info> const *mapper=NULL ) ;
# endif 

};

}

#endif
