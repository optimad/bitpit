#include <iostream>
#include <fstream> 
#include <vector>
#include <array>
#include <map>


#include "Operators.hpp"
#include "VTK.hpp"

#include "rbf.hpp"

/*!
 * deformation by sine modulation
 * @param[in] P point coordinates
 * @param[in] amp amplitude in three space directions
 * @param[in] wavr wave lengththree space directions
 * @return displacement vector
 */
std::array<double,3> modulation( std::array<double,3> P, std::array<double,3> amp, std::array<double,3> wave ){

    std::array<double,3>    d;

    d[0] = amp[0]*sin( wave[0] *M_PI *P[0] ) ;
    d[1] = amp[1]*sin( wave[1] *M_PI *P[0] ) ;
    d[2] = amp[2]*sin( wave[2] *M_PI *P[0] ) ;

    return d;

}

/*!
 * deformation by translation
 * @param[in] P point coordinates
 * @param[in] tras amplitude in three space directions
 * @return displacement vector
 */
std::array<double,3> traslation( std::array<double,3> P, std::array<double,3> tras){

    return tras;
}

/*!
 * deformation rotation
 * @param[in] P point coordinates
 * @param[in] origin center of rotation
 * @param[in] rot rotation vector
 * @return displacement vector
 */
std::array<double,3> rotation( std::array<double,3> P, std::array<double,3> origin, std::array<double,3> rot ){

    return crossProduct( rot, P - origin ) ;

};

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

    std::array<double,3> origine ;
    double          altezza, base ;


    std::vector<double>  xpoints, ypoints ;
    std::array<int,3>    np, nc ;
    int             i, j, k, I_ ;



    dom_oriz    =   100 ;
    dom_vert    =   50 ;

    delta_oriz  =   0.01 ;
    delta_vert  =   0.01 ;

    expansion   =   1.2 ;

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

        type.resize( np[0]*np[1]*np[2] ) ;

        for( i=0; i<np[0]; ++i){
            if( xpoints[i] < -base/2 ) body_beg[0] = i+1 ;
            if( xpoints[i] <  base/2 ) body_end[0] = i   ;
        };

        for( i=0; i<np[1]; ++i){
            if( ypoints[i] < -altezza/2 ) body_beg[1] = i+1 ;
            if( ypoints[i] <  altezza/2 ) body_end[1] = i   ;
        };

        for( i=0; i<points.size(); i++){

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



int main() {

    // Variables
    std::vector< std::array<double,3> >              points ;
    std::vector< std::vector<int> >       connectivity ;
    std::vector< int >               type ;

    int                         nP, nC ;

    std::array<double,3>                        ampl, wave ;
    std::array<double,3>                        tras ;
    std::array<double,3>                        orig, vrot ;
    std::vector< std::array<double,3> >         deformation ;

    int                         i, j, I_ ;
    std::array<double,3>                        P, tmp0, tmp1, tmp2 ;

    // create grid
    createCMesh( points, connectivity, type) ;
    nP = points.size() ;
    nC = connectivity.size() ;

    deformation.resize( nP ) ;

    { // output of undeformed grid, deformation vector and node type
        bitpit::VTKUnstructuredGrid   output( "./", "orgGrid", bitpit::VTKElementType::PIXEL, points, connectivity );

        output.addData( "type", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, type ) ;
        output.addData( "deformation", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, deformation ) ;

        output.write() ;
    }


    for( int loop=0; loop<1; ++loop){

        std::cout << loop << std::endl;

        // introduce deformations
        ampl[0] = 0.0 ;
        ampl[1] = 0.1 ;
        ampl[2] = 0.0 ;

        wave[0] = 0. ;
        wave[1] = 2. ;
        wave[2] = 0. ;

        tras[0] = 0.0 ;
        tras[1] = 0.0 ;
        tras[2] = 0.0 ;

        orig[0] = 0.0 ;
        orig[1] = 0.0 ;
        orig[2] = 0.0 ;

        vrot[0] = 0.0 ;
        vrot[1] = 0.0 ;
        vrot[2] = 0.0 ; //M_PI/4. ;


        for( i=0; i<nP; ++i){

            if( type[i] == 1){

                P   =   points[i] ;

                tmp0 =   modulation( P, ampl, wave) ;
                tmp1 =   traslation( P, tras) ;
                tmp2 =   rotation( P, orig, vrot) ;

                deformation[i] = tmp0+tmp1+tmp2;
            }

            else{
                deformation[i].fill(0.) ;

            };
        };



        {//RBF
            int                 nNodes(0) ;
            double              maxDeform(0);
            bitpit::RBF         meshMorph;
            std::vector<double> values ;
            std::vector<int>    active ;
            std::map<int,int>   nodes ;

            for( i=0; i<nP; ++i){
                if( type[i] == 1 ){
                    j = meshMorph.addNode(points[i]) ;
                    nodes.insert( std::pair<int,int>(j,i) ) ;
                    nNodes++ ;
                }
            }

            values.resize(nNodes) ;

            for(int k=0; k<2; ++k){

                j= 0 ;
                for( i=0; i<nP; ++i){

                    if( type[i] == 1 ){
                        values[j] = deformation[i][k] ;
                        maxDeform = std::max( maxDeform, std::abs(values[j]) ) ;
                        ++j ;
                    }

                    //if( type[i] == 2 ){
                    //    values[j] = 0. ;
                    //    ++j ;
                    //}
                }

                meshMorph.addData(values) ;
            }

            meshMorph.setSupportRadius( 0.5 ) ;
            //meshMorph.solve() ;
            meshMorph.greedy(0.00001) ;

            active.resize( nP );

            for( const auto &node : nodes ){
                active[node.second] = meshMorph.isActive(node.first) ;
            };



            std::vector<double> disp ;
            for( auto & point : points){
                disp = meshMorph.evalRBF(point) ;
                point[0] += disp[0] ;
                point[1] += disp[1] ;
            };


            // output of deformed grid
            bitpit::VTKUnstructuredGrid   output( "./", "defGrid", bitpit::VTKElementType::PIXEL, points, connectivity );
            output.addData( "active", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, active ) ;
            output.addData( "type", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, type ) ;
            output.addData( "deformation", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, deformation ) ;
            output.write() ;
        }
        }
        return 0;
    }
