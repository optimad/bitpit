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

#include <iostream>
#include <fstream> 
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_operators.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_RBF.hpp"

using namespace bitpit;

//Description. Directly Parametrize A Surface Mesh w/ RBFs

/*! Find bounding box of a 3D point clouds.
 * @param[in] list of points belonging to clouds
 * @result	  limits of the point cloud in abs ref system
 */
std::vector<std::array<double,2> > findBoundingBox(std::vector<std::array<double,3> > & cPoints){
	
	std::vector<std::array<double,2> > lims(3, std::array<double,2>{{0.0,0.0}});
	if (cPoints.size()==0) return lims;
	
	// Compute bounding box limits ---------------------------------------------- //
	lims[0][0] = lims[0][1] = cPoints[0][0];
	lims[1][0] = lims[1][1] = cPoints[0][1];
	lims[2][0] = lims[2][1] = cPoints[0][2];
	
	for (int i = 1; i < (int)(cPoints.size()); i++) {
		for(int loc=0; loc<3; ++loc){
			lims[loc][0] = std::min(lims[loc][0], cPoints[i][loc]);
			lims[loc][1] = std::max(lims[loc][1], cPoints[i][loc]);
		}
	} //next j
	
	return(lims);
}

/*!
 * transforms a linear index to a cartesian index
 * @param[in] I_ linear index
 * @param[in] np number of points in each direction
 * @return cartesian index
 */
std::array<int,3> fromLinearToCartesian(int I_, std::array<int,3> np ){

    std::array<int,3>    id;

    int     PointsInIJPlane( np[0]*np[1] ) ;

    id[0] = I_ % np[0] ;
    id[2] = I_ / PointsInIJPlane;
    id[1] =  (I_ - id[2] *PointsInIJPlane ) /np[0]  ;

    return id;
}

/*!
 * transforms a cartesian index to a linear index
 * @param[in] i cartesian index in first direction
 * @param[in] j cartesian index in second direction
 * @param[in] k cartesian index in third direction
 * @param[in] np number of points in each direction
 * @return linear index
 */
int fromCartesianToLinear(int i, int j, int k, std::array<int,3> np){

    int     PointsInIJPlane( np[0]*np[1] ) ;

    return PointsInIJPlane *k + np[0]*j + i;

};

/*!
 * Creates a 1D stretched grid.
 * The gird will start at coordinate =0 and have constant spacing (=delta) on the body (=body).
 * After that the spacin will increase exponantially (=ratio) until the final coordinate is reached (=final_coordinate)
 * param[in] final_coordinate end coordinate 
 * param[in] body coordinate until constant spacing is used
 * param[in] delta constant spacing
 * param[in] ratio expansion ratio after constant spacing
 * @return 1D grid point distribution
 */
std::vector<double> createStretchedGrid( double final_coordinate, double body, double delta, double ratio ){

    int                 n, i;  
    double              coord ;
    std::vector<double>      v, points;

    v.reserve(100) ;
    points.reserve(100) ;

    body = body /2. ;

    i = 0 ;
    coord = 0; 

    // crea coodinate asse y positive
    while( coord < final_coordinate){
        v.push_back( coord ) ;
        coord += delta ;

        if( coord >= body)
            delta = delta * ratio ;

        i++;
    }	

    n = i ;

    for( i=0; i< n-1 ; ++i){
        points.push_back( -v[n-i-1] ) ;
    };

    for( i=0 ; i< n  ; ++i){
        points.push_back( v[i] ) ;
    };


    return points ;
};

/*!
 * Creates a 2D cartesian mesh around a rectangle with constant spacing on the rectangle.
 * param[out] points node coordinates
 * param[out] connectivity grid connectivity
 * param[out] type type of node [-1=solid, 0=fluid, 1=solid boundary, 2=farfield boundary ]
 */
void createCMesh( std::vector< std::array<double,3> > &points, std::vector< std::vector<int> > &connectivity, std::vector<int> &type) {



    //==================================================================================================
    //			Dichiarazione
    //==================================================================================================
    double          dom_oriz, dom_vert ;
    double          delta_oriz, delta_vert ;
    double          expansion ;

    double          altezza, base ;


    std::vector<double>  xpoints, ypoints ;
    std::array<int,3>    np, nc ;
    int             i, j, k, I_ ;



    dom_oriz    =   100 ;
    dom_vert    =   50 ;

    delta_oriz  =   0.5 ;
    delta_vert  =   0.5 ;

    expansion   =   1.00 ;

    altezza     =   0.5;
    base        =   1;

    // normalizzazione
    double      lref = base ;

    dom_oriz = dom_oriz / lref ;
    dom_vert = dom_vert / lref ;

    delta_oriz  =   delta_oriz /lref ;
    delta_vert  =   delta_vert /lref ;

    altezza = altezza / lref ;
    base    = base / lref ;


    //==================================================================================================
    //			CREO DISTRIBUZIONE PUNTI
    //==================================================================================================
    {

        xpoints = createStretchedGrid( dom_oriz, base, delta_oriz, expansion ) ;
        ypoints = createStretchedGrid( dom_vert, altezza, delta_vert, expansion ) ;

        np[0] = xpoints.size() ;
        np[1] = ypoints.size() ;
        np[2] = 1 ;

        nc[0] = np[0] -1  ;
        nc[1] = np[1] -1  ;
        nc[2] = 1 ;
    }


    //==================================================================================================
    //                       CREO MATRICE GRIGLIA CONTENENTI I VETTORI DELLA GRIGLIA
    //==================================================================================================
    {

        std::array<double,3>    temp ;
        k=0;
        for(j=0; j<np[1];  j++){
            for( i=0; i<np[0]; i++){
                temp[0] = xpoints[i] ;
                temp[1] = ypoints[j] ;
                temp[2] = 0. ;

                points.push_back( temp ) ;

            }
        }

    }

    //==================================================================================================================
    { //              CONNECTIVITY - GRIGLIA NON STRUTTURATA

        connectivity.resize( nc[0]*nc[1]*nc[2] );
        I_ = 0;


        k=0;
        for(j=0; j<nc[1];  j++){
            for( i=0; i<nc[0]; i++){

                connectivity[I_].resize(4) ;
                connectivity[I_].shrink_to_fit() ;

                connectivity[I_][0]  = fromCartesianToLinear( i,   j,   k, np ) ;
                connectivity[I_][1]  = fromCartesianToLinear( i+1, j,   k, np ) ;
                connectivity[I_][2]  = fromCartesianToLinear( i,   j+1, k, np ) ;
                connectivity[I_][3]  = fromCartesianToLinear( i+1, j+1, k, np ) ;

                I_++ ;
            }
        }
    }

    //==================================================================================================
    { //           CORPO E CC

        int     d ;
        bool    found;
        std::array<int,3>    index, body_beg, body_end ;
        body_beg.fill(0);
        body_end.fill(0);

        type.resize( np[0]*np[1]*np[2] ) ;

        for( i=0; i<np[0]; ++i){
            if( xpoints[i] < -base/2 ) body_beg[0] = i+1 ;
            if( xpoints[i] <  base/2 ) body_end[0] = i   ;
        };

        for( i=0; i<np[1]; ++i){
            if( ypoints[i] < -altezza/2 ) body_beg[1] = i+1 ;
            if( ypoints[i] <  altezza/2 ) body_end[1] = i   ;
        };

        for( i=0; i<(int)(points.size()); i++){

            index   = fromLinearToCartesian( i, np ) ;
            found   = false ;

            { // check if domain boundaries
                for( d=0; d<2; ++d){
                    found = ( found || index[d] == 0 || index[d] == np[d]-1 ) ;
                };

                if( found) type[i] = 2 ;
            }

            if( !found ){ // check if solid boundary
                found = found  || ( index[0] == body_beg[0] && index[1] >= body_beg[1] && index[1] <= body_end[1]) ;
                found = found  || ( index[0] == body_end[0] && index[1] >= body_beg[1] && index[1] <= body_end[1]) ;
                found = found  || ( index[1] == body_beg[1] && index[0] >= body_beg[0] && index[0] <= body_end[0]) ;
                found = found  || ( index[1] == body_end[1] && index[0] >= body_beg[0] && index[0] <= body_end[0]) ;

                if(found) type[i] = 1 ;
            }


            if( !found){ // check if solid body

                found = true ;
                for( d=0; d<2; ++d){
                    found = found && ( index[d] > body_beg[d] ) && ( index[d] < body_end[d] );
                };

                if(found) type[i] = -1 ;

            }


            if( !found){ // // default fluid domain
                type[i] = 0 ;
            };


        }


    };


    return;
};

/*!
* Subtest 001
*
* Testing mesh deformation.
*/
int subtest_001()
{
    // Variables
    std::vector< std::array<double,3> >              points ;
    std::vector< std::vector<int> >       connectivity ;
    std::vector< int >               type ;

//    int                         nP, nC ;

    std::vector< std::array<double,3> >         controlNodes;
	std::vector< double >         				zDispl ;
	
    // create grid
    createCMesh( points, connectivity, type) ;
//     nP = points.size() ;
//     nC = connectivity.size() ;
	
	int sizeCN =10;
	controlNodes.resize( sizeCN) ;
	zDispl.resize(sizeCN);
	
    { // output of undeformed grid, deformation vector and node type
        bitpit::VTKUnstructuredGrid   output( "./", "orgGrid", bitpit::VTKElementType::PIXEL );
        output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points ) ;
        output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity ) ;
        output.setDimensions(connectivity.size(),points.size());
        output.write() ;
    }

    //adding control nodes to mesh extrema and randomly fill nodes in the mesh interior 
    //Assign random values to zDispl vector;
    {
		std::vector<std::array<double,2> > lims = findBoundingBox(points);
		
		controlNodes[0] = std::array<double,3>{{lims[0][0], lims[1][0], lims[2][1]}};
		controlNodes[1] = std::array<double,3>{{lims[0][1], lims[1][0], lims[2][1]}};
		controlNodes[2] = std::array<double,3>{{lims[0][1], lims[1][1], lims[2][1]}};
		controlNodes[3] = std::array<double,3>{{lims[0][0], lims[1][1], lims[2][1]}};
		
		srand(time(NULL));
		for(int i=4; i<sizeCN; ++i){
			for(int j=0; j<2; ++j){
				double wg =  (double)(rand())/RAND_MAX;
				double a = std::min(0.99,std::max(0.01, wg));
				controlNodes[i][j] = lims[j][0] + a*(lims[j][1] - lims[j][0]); 
			}
			controlNodes[i][2] = lims[2][1];
		}
		
		double maxDim = 0.0;
		std::array<double,3> span;
		for(int i=0; i<3; ++i){
			span[i] = lims[i][1] - lims[i][0];
			maxDim = std::max(maxDim, span[i]);
		}
		
		{
		int index = 0;
		for(auto && value : zDispl){
			//value = 0.5*maxDim*((double) (rand()) / RAND_MAX - 0.5);
			value = 0.5*maxDim*(-1.0*controlNodes[index][0]/span[0] + 0.5);
			++index;
		}
		}
	}
	
	{ // output of undeformed control cloud
		std::vector<int> conn;
		int sCN = controlNodes.size();
		for(int i=0; i<sCN; ++i) conn.push_back(i);
		bitpit::VTKUnstructuredGrid   output( "./", "orgCPointCloud", bitpit::VTKElementType::VERTEX);
        output.setGeomData( bitpit::VTKUnstructuredField::POINTS, controlNodes ) ;
        output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, conn ) ;
        output.setDimensions(conn.size(),controlNodes.size());
		output.setCodex(bitpit::VTKFormat::ASCII);
		output.write() ;
	}
	
    
     {//RBF used to directly parametrize surface
     		bitpit::RBF   paraMorph;
			RBFBasisFunction funct = RBFBasisFunction::WENDLANDC2;
			
			paraMorph.setFunction(funct);
			
			double maxVal = 0.0;
			maxval(zDispl, maxVal);
            paraMorph.setSupportRadius( 10.0*maxVal ) ;
			
			std::vector<int> nIndex = paraMorph.addNode(controlNodes);
			if(nIndex.size() != controlNodes.size())	return 1;
			int nData = paraMorph.addData(zDispl);
			if(nData != 1)	return 1;
			
			int err = paraMorph.solve() ;
			if(err > 0 ) 	return 1;
			//paraMorph.greedy(0.001) ;
			std::vector<double> disp ;
            for( auto & point : points){
                disp = paraMorph.evalRBF(point) ;
                point[2] += disp[0] ;
            };
	 }
	 
	{// output of deformed grid
        bitpit::VTKUnstructuredGrid   output( "./", "defGrid", bitpit::VTKElementType::PIXEL );
        output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points ) ;
        output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity ) ;
        output.setDimensions(connectivity.size(),points.size());
        output.write() ;
	}
	
	{ // output of undeformed control cloud
		for(int i=0; i<(int)(controlNodes.size()); ++i){
			controlNodes[i][2] += zDispl[i];
		}
		std::vector<int> conn;
		for(int i=0; i<(int)(controlNodes.size()); ++i) conn.push_back(i);
		
		bitpit::VTKUnstructuredGrid   output( "./", "defCPointCloud", bitpit::VTKElementType::VERTEX);
        output.setGeomData( bitpit::VTKUnstructuredField::POINTS, controlNodes ) ;
        output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, conn ) ;
		output.setDimensions(conn.size(),controlNodes.size());
		output.setCodex(bitpit::VTKFormat::ASCII);
		output.write() ;
	}
	
	
        return 0;
}

// ========================================================================== //
// MAIN                                                                       //
// ========================================================================== //
int main(int argc, char *argv[])
{
    // ====================================================================== //
    // INITIALIZE MPI                                                         //
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc,&argv);
#else
	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
#endif

    // ====================================================================== //
    // VARIABLES DECLARATION                                                  //
    // ====================================================================== //

    // Local variabels
    int                             status = 0;

    // ====================================================================== //
    // RUN SUB-TESTS                                                          //
    // ====================================================================== //
    try {
        status = subtest_001();
        if (status != 0) {
            return (10 + status);
        }
    } catch (const std::exception &exception) {
        bitpit::log::cout() << exception.what();
        exit(1);
    }

    // ====================================================================== //
    // FINALIZE MPI                                                           //
    // ====================================================================== //
#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

    return status;
}
