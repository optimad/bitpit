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

# ifndef __BITPIT_LEVELSET_CARTESIAN_HPP__
# define __BITPIT_LEVELSET_CARTESIAN_HPP__


namespace bitpit{

class VolCartesian;
class LevelSetKernel;

class LevelSetCartesian : public LevelSetKernel{

    private:
    VolCartesian*                               m_cartesian ;       /**< Pointer to underlying cartesian mesh*/

    public:
    virtual ~LevelSetCartesian();
    LevelSetCartesian( VolCartesian & );

    VolCartesian *                              getCartesianMesh() const;
    double                                      computeRSearchFromCell( long id ) ;

};

}

#endif
