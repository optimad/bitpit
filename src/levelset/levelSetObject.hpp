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

# ifndef __BITPIT_LEVELSET_OBJECT_HPP__
# define __BITPIT_LEVELSET_OBJECT_HPP__

// Standard Template Library
# include <iostream>
# include <array>
# include <vector>
# include <unordered_map>

# include "bitpit_IO.hpp"
# include "levelSetCommon.hpp"

namespace bitpit{

namespace adaption{
    struct Info;
}
class SendBuffer;
class RecvBuffer;

class LevelSet;
class LevelSetKernel;

class LevelSetObject : public VTKBaseStreamer{

    friend LevelSet;

    private:
    int                                         m_id;           /**< identifier of object */

    protected:
    LevelSetObject(int);
    LevelSetObject(const LevelSetObject &other) = default;

    void                                        setKernel(LevelSetKernel *);
    LevelSetKernel *                            getKernel();

    void                                        clear();

    void                                        setSizeNarrowBand(double) ;

    virtual void                                computeLSInNarrowBand(bool);
    virtual void                                updateLSInNarrowBand(const std::vector<adaption::Info> &, bool);
    void                                        clearAfterMeshAdaption(const std::vector<adaption::Info>&);

    virtual void                                propagateSign() ;

    void                                        dump(std::ostream &);
    void                                        restore(std::istream &);

# if BITPIT_ENABLE_MPI
    void                                        exchangeGhosts() ;
    void                                        communicate( const std::unordered_map<int,std::vector<long>> &,
                                                             const std::unordered_map<int,std::vector<long>> &,
                                                             std::vector<adaption::Info> const *mapper=NULL );
# endif


    LevelSetKernel*                             m_kernelPtr;    /**< pointer to kernel */
    double                                      m_narrowBand;   /**< Size of narrow band */

    virtual void                                _clear();
    virtual void                                _clearAfterMeshAdaption(const std::vector<adaption::Info>&) ;
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

    const LevelSetKernel *                      getKernel() const;

    virtual LevelSetObject*                     clone() const =0;

    int                                         getId() const ;
    virtual bool                                isPrimary() const ;

    virtual LevelSetInfo                        getLevelSetInfo(const long &) const =0;
    virtual double                              getLS(const long &) const =0; 
    virtual std::array<double,3>                getGradient(const long &) const =0 ; 
    std::array<double,3>                        computeProjectionPoint(const long &) const;

    virtual int                                 getPart(const long &) const ;
    virtual std::array<double,3>                getNormal(const long &) const; 

    short                                       getSign(const long &) const;

    bool                                        isInNarrowBand(const long &) const;
    double                                      getSizeNarrowBand() const;

    LevelSetIntersectionStatus                  intersectSurface(const long &, LevelSetIntersectionMode=LevelSetIntersectionMode::FAST_FUZZY) const;
    virtual double                              getSurfaceFeatureSize(const long &) const;
    virtual double                              getMinSurfaceFeatureSize() const;
    virtual double                              getMaxSurfaceFeatureSize() const;

    void                                        enableVTKOutput(LevelSetWriteField field, bool enable=true);
    void                                        flushData(std::fstream &, std::string, VTKFormat);


};

}

#endif
