// ========================================================================== //
//                 - GRID MANAGER FOR CARTESIAN MESHES -                      //
//                                                                            //
// Grid manager for cartesian meshes                                          //
// ========================================================================== //
// INFO                                                                       //
// ========================================================================== //
// Author      :   Alessandro Alaia                                           //
// Version     :   v2.0                                                       //
//                                                                            //
// All rights reserved.                                                       //
// ========================================================================== //

// -------------------------------------------------------------------------- //
template < class T >
void UCartMesh::CellData2PointData(
        vector< T > &CellData,
        vector< T > &PointData
        ) {

    // ========================================================================== //
    // template < class T >                                                       //
    // void UCartMesh::CellData2PointData(                                //
    //     vector< T > &CellData,                                                 //
    //     vector< T > &PointData)                                                //
    //                                                                            //
    // Convert cell data into point data                                          //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - CellData    : vector< T >, cell data. T can be of anc[1] copy-constructible //
    //                 type s.t. summation (difference) and multiplication by a   //
    //                 constant are defined.                                      //
    // - PointData   : vector< T > with point data.                               //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //        
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int            K, J;
    int            ip, jp, kp;
    vector<int>     PointIter ;

    // Counters
    int            i, j, k, l, m, n;

    // ========================================================================== //
    // RESIZE OUTPUT VARIABLES                                                    //
    // ========================================================================== //
    PointData.resize( getNPoints() );
    PointIter.resize( getNPoints() );

    // ========================================================================== //
    // CONVERT CELL DATA INTO POINT DATA                                          //
    // ========================================================================== //
    for (k = 0; k < nc[2]; k++) {
        for (j = 0; j < nc[1]; j++) {
            for (i = 0; i < nc[0]; i++) {

                // Cell index
                K = CellLinearId(i, j, k);

                for (n = 0; n < dim-1; n++) {
                    for (m = 0; m < 2; m++) {
                        for (l = 0; l < 2; l++) {

                            // Point index
                            ip = i + l;
                            jp = j + m;
                            kp = j + n;
                            J = PointLinearId(ip,jp,kp);

                            PointData[J] = PointData[J] + CellData[K]; 
                            PointIter[J]++ ;
                        } //next n
                    } //next m
                } //next l

            } //next k
        } //next j
    } //next i


    for( J=0; J<getNPoints(); ++J){
        PointData[J] = PointData[J] / ((float) PointIter[J]) ;
    };

    return; 

};

// -------------------------------------------------------------------------- //
template < class T, typename ... T2 >
void UCartMesh::CellData2PointData(
        vector< T > &CellData,
        vector< T > &PointData,
        T2      &... others
        ) {

    // ========================================================================== //
    // template < class T, typename ... T2 >                                      //
    // void UCartMesh::CellData2PointData(                                //
    //     vector< T > &CellData,                                                 //
    //     vector< T > &PointData,                                                //
    //     T2      &... others)                                                   //
    //                                                                            //
    // Convert an arbitrary set of point cell data into point data.               //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - CellData    : vector< T >, cell data. T can be of anc[1] copy-constructible //
    //                 type s.t. summation (difference) and multiplication by a   //
    //                 constant are defined.                                      //
    // - PointData   : vector< T > with point data.                               //
    // - others      : T2 (optional) other cell data set to be converted.         //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //        
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    // none

    // ========================================================================== //
    // CONVERT FIRST CELL DATA SET INTO POINT DATA                                //
    // ========================================================================== //
    CellData2PointData(CellData, PointData);

    // ========================================================================== //
    // CONVERT OTHER CELL DATA SET INTO POINT DATA                                //
    // ========================================================================== //
    CellData2PointData(others ...);

    return; 
}

// -------------------------------------------------------------------------- //
template < class T >
void UCartMesh::PointData2CellData(
        vector< T > &PointData,
        vector< T > &CellData
        ) {

    // ========================================================================== //
    // template < class T >                                                       //
    // void UCartMesh::PointData2CellData(                                //
    //     vector< T > &PointData,                                                //
    //     vector< T > &CellData)                                                 //
    //                                                                            //
    // Convert point data into cell data.                                         //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - PointData   : vector< T >, point data. T can be anc[1] copy-constructible   //
    //                 type s.t. summation (difference), and multiplication by a  //
    //                 constant are defined.                                      //
    // - CellData    : vector< T >, point data.                                   //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    int     K, J;
    int     ip, jp, kp;
    double  factor ;

    // Counters
    int     i, j, k, l, m, n;

    factor  =   pow(0.5,dim) ;

    // ========================================================================== //
    // RESIZE OUTPUT VARIABLES                                                    //
    // ========================================================================== //
    CellData.resize(getNCells(), 0.0);

    // ========================================================================== //
    // CONVERT POINT DATA TO CELL DATA                                            //
    // ========================================================================== //
    for (k = 0; k < nc[2]; k++) {
        for (j = 0; j < nc[1]; j++) {
            for (i = 0; i < nc[0]; i++) {

                K = CellLinearId(i,j,k);
                for (n = 0; n < 2; n++) {
                    for (m = 0; m < 2; m++) {
                        for (l = 0; l < 2; l++) {
                            ip = i + l;
                            jp = j + m;
                            kp = k + n;
                            J = PointLinearId(ip,jp,kp);
                            CellData[K] = CellData[K] + factor * PointData[J];
                        } //next n
                    } //next m
                } //next l
            } //next k
        } //next j
    } //next i

    return; 
};

// -------------------------------------------------------------------------- //
template < class T, typename ... T2 >
void UCartMesh::PointData2CellData(
        vector< T > &PointData,
        vector< T > &CellData,
        T2      &... others
        ) {

    // ========================================================================== //
    // template < class T, typename ... T2 >                                      //
    // void UCartMesh::PointData2CellData(                                //
    //     vector< T > &PointData,                                                //
    //     vector< T > &CellData,                                                 //
    //     T2      &... others)                                                   //
    //                                                                            //
    // Convertes and arbitrary set of point data into cell data.                  //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - PointData   : vector< T >, point data. T can be anc[1] copy-constructible   //
    //                 type s.t. summation (difference), and multiplication by a  //
    //                 constant are defined.                                      //
    // - CellData    : vector< T >, point data.                                   //
    // - others      : T2 (optional) other data set to be converted.              //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    // none

    // ========================================================================== //
    // CONVERTES FIRST POINT DATA SET INTO CELL DATA                              //
    // ========================================================================== //
    PointData2CellData(PointData, CellData);

    // ========================================================================== //
    // ITERATIVELY CONVERTES OTHER DATA SET                                       //
    // ========================================================================== //
    PointData2CellData(others ...);

    return; 
}

// -------------------------------------------------------------------------- //
template < class T >
void UCartMesh::interpolateCellData(
        darray3E     &P,
        vector< T >  &field,
        T            &value
        ) {

    // ========================================================================== //
    // template < class T >                                                       //
    // void UCartMesh::interpolateCellData(                               //
    //     dvector1D   &P,                                                        //
    //     vector< T > &field,                                                    //
    //     T           &value)                                                    //
    //                                                                            //
    // Interpolate cell data at a given point (1st order interpolation).          //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - P      : dvector1D, point coordinates                                    //
    // - field  : vector< T >, cell data. T can be anc[1] copy-constructible type    //
    //            s.t. summation and multiplication by a constant are defined.    //
    // - value  : T, interpolation result                                         //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    if( PointInGrid(P) ){


        // Local variables
        iarray3E    i0, i1 ;
        darray3E    w0, w1 ;

        i0.fill(0);   i1.fill(0) ;
        w0.fill(0.5); w1.fill(0.5) ;

        // Counters
        int     d;

        i0 = CellCartesianId(P);

        for( d=0; d<dim; ++d){

            //            i0[d] = max( min( i0[d], nc[d]-1) ,1 );

            // Find cell index
            if( P[d] < center[d][i0[d] ] ){
                i0[d] = max(0, i0[d]-1) ; 
            };

            i1[d] = min( i0[d] +1, nc[d]-1) ;

            // Interpolation weights
            if( i0[d] == i1[d] ){
                w0[d] = 1. ;
                w1[d] = 0. ;
            }

            else{
                w1[d] = max(0.0, min(1.0, ( P[d] - center[d][i0[d]]) /h[d] ) ) ;
                w0[d] = 1.0 - w1[d] ;  
            };

        };


        value = 
            w0[0] * w0[1] * w0[2] *field[CellLinearId(i0[0],i0[1],i0[2])]
            +  w1[0] * w0[1] * w0[2] *field[CellLinearId(i1[0],i0[1],i0[2])]
            +  w0[0] * w1[1] * w0[2] *field[CellLinearId(i0[0],i1[1],i0[2])]
            +  w1[0] * w1[1] * w0[2] *field[CellLinearId(i1[0],i1[1],i0[2])]
            +  w0[0] * w0[1] * w1[2] *field[CellLinearId(i0[0],i0[1],i1[2])]
            +  w1[0] * w0[1] * w1[2] *field[CellLinearId(i1[0],i0[1],i1[2])]
            +  w0[0] * w1[1] * w1[2] *field[CellLinearId(i0[0],i1[1],i1[2])]
            +  w1[0] * w1[1] * w1[2] *field[CellLinearId(i1[0],i1[1],i1[2])];
    };

    return; 

};

// -------------------------------------------------------------------------- //
template < class T, typename ... T2 >
void UCartMesh::interpolateCellData(
        darray3E     &P,
        vector< T >  &field,
        T            &value,
        T2       &... others
        ) {

    // ========================================================================== //
    // template < class T, typename ... T2 >                                      //
    // void UCartMesh::interpolateCellData(                               //
    //    dvector1D    &P,                                                        //
    //    vector< T >  &field,                                                    //
    //    T            &value,                                                    //
    //    T2       &... others)                                                   //
    //                                                                            //
    // Interpolate a set of cell data at a specified point.                       //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - P      : dvector1D, point coordinates                                    //
    // - field  : vector< T >, cell data. T can be anc[1] copy-constructible type    //
    //            s.t. summation and multiplication by a constant are defined.    //
    // - value  : T, interpolation result                                         //
    // - others : T2 (optional) other cell data sets to be used for interpolation //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    // none

    // ========================================================================== //
    // INTERPOLATION USING THE FIRST CELL DATA SET                                //
    // ========================================================================== //
    interpolateCellData(P, field, value);

    // ========================================================================== //
    // INTERPOLATION USING THE THE OTHER DATASETS                                 //
    // ========================================================================== //
    interpolateCellData(P, others...);

    return; 
}

// -------------------------------------------------------------------------- //
template < class T >
void UCartMesh::interpolatePointData(
        darray3E    &P,
        vector<T>   &field,
        T           &value
        ) {

    // ========================================================================== //
    // template < class T >                                                       //
    // void UCartMesh::interpolatePointData(                              //
    //     dvector1D   &P,                                                        //
    //     vector<T>   &field,                                                    //
    //     T           &value)                                                    //
    //                                                                            //
    // Interpolate point data at a given point.                                   //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - P     : dvector1D, point x, y, z coordinates                             //
    // - field : vector<T>, point data to be interpolated. T can be anc[1] copy-     //
    //           constructible type s.t. summation (difference) and               //
    //           multiplication by a scalar must be defined.                      //
    // - value : T, interpolation result                                          //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    if( PointInGrid(P) ){

        // Local variables
        iarray3E    i0, i1 ;
        darray3E    w0, w1 ;

        i0.fill(0);   i1.fill(0) ;
        w0.fill(0.5); w1.fill(0.5) ;

        // Counters
        int     d;

        // ========================================================================== //
        // INTERPOLATE POINT DATA                                                     //
        // ========================================================================== //

        // Closest grid point
        i0 = CellCartesianId( P) ;


        for(d=0; d<dim; ++d){
            i1[d] = i0[d] +1 ;

            w1[d] = max(0.0, min(1.0, ( P[d] - edge[d][i0[d]]) /h[d] ) ) ;
            w0[d] = 1.0 - w1[d] ;
        };





        // Interpolation
        value = 
            w0[0] * w0[1] * w0[2] *field[PointLinearId(i0[0],i0[1],i0[2])]
            +  w1[0] * w0[1] * w0[2] *field[PointLinearId(i1[0],i0[1],i0[2])]
            +  w0[0] * w1[1] * w0[2] *field[PointLinearId(i0[0],i1[1],i0[2])]
            +  w1[0] * w1[1] * w0[2] *field[PointLinearId(i1[0],i1[1],i0[2])]
            +  w0[0] * w0[1] * w1[2] *field[PointLinearId(i0[0],i0[1],i1[2])]
            +  w1[0] * w0[1] * w1[2] *field[PointLinearId(i1[0],i0[1],i1[2])]
            +  w0[0] * w1[1] * w1[2] *field[PointLinearId(i0[0],i1[1],i1[2])]
            +  w1[0] * w1[1] * w1[2] *field[PointLinearId(i1[0],i1[1],i1[2])];
    };


    return; 

};

// -------------------------------------------------------------------------- //
template < class T, typename ... T2 >
void UCartMesh::interpolatePointData(
        darray3E    &P,
        vector<T>   &field,
        T           &value,
        T2      &... others
        ) {

    // ========================================================================== //
    // template < class T, typename ... T2 >                                      //
    // void UCartMesh::interpolatePointData(                              //
    //     dvector1D   &P,                                                        //
    //     vector<T>   &field,                                                    //
    //     T           &value,                                                    //
    //     T2      &... others)                                                   //
    //                                                                            //
    // Interpolate point datasets at a given point.                               //
    // ========================================================================== //
    // INPUT                                                                      //
    // ========================================================================== //
    // - P        : dvector1D, point x, y, z coordinates                          //
    // - field    : vector< T >, point dataset to be interpolated. T can be anc[1]   //
    //              copy-constructible type, s.t. summation (difference) and      //
    //              multiplication by a scalare are defined.                      //
    // - value    : T, interpolation result.                                      //
    // - others   : T2 (optional), interpolation results                          //
    // ========================================================================== //
    // OUTPUT                                                                     //
    // ========================================================================== //
    // - none                                                                     //
    // ========================================================================== //

    // ========================================================================== //
    // VARIABLES DECLARATION                                                      //
    // ========================================================================== //

    // Local variables
    // none

    // Counters
    // none

    // ========================================================================== //
    // INTERPOLATE THE FIRST VERTEX DATASET                                       //
    // ========================================================================== //
    interpolatePointData(P, field, value);

    // ========================================================================== //
    // ITERATIVELY INTERPOLATE THE OTHER VERTEX DATASETS                          //
    // ========================================================================== //
    interpolatePointData(P, others ...);

    return; 

}

// -------------------------------------------------------------------------------------- //
template <class T>
void UCartMesh::Export_CellData_vtr(string      filename,
        string      dataname,
        vector< T > &CellData) {

    // ====================================================================================== //
    // template <class T>                                                                     //
    // void UCartMesh::Export_CellData_vtr(string       filename,                     //
    //                                             string       dataname,                     //
    //                                             vector< T > &CellData)                     //
    //                                                                                        //
    // Export cell data into a .vtr file.                                                     //
    // ====================================================================================== //
    // INPUT                                                                                  //
    // ====================================================================================== //
    // - filename    : string, .vtr file name                                                 //
    // - dataname    : string, data name                                                      //
    // - CellData    : [nc[0]*nc[1]*nc[2]-by-1] vector< T > with cell data                             //
    // ====================================================================================== //
    // OUTPUT                                                                                 //
    // ====================================================================================== //
    // - none                                                                                 //
    // ====================================================================================== //

    // ====================================================================================== //
    // VARIABLES DECLARATION                                                                  //
    // ====================================================================================== //

    // Local variables
    ofstream        file_handle;

    // Counters

    // ====================================================================================== //
    // EXPORT CELL DATA                                                                       //
    // ====================================================================================== //

    // Open file
    Open_vtr(file_handle, trim(filename));

    // Write mesh data
    Write_vtrMeshData(file_handle, nc[0], nc[1], nc[2], edge[0], edge[1], edge[2] );

    // Write cell data
    Open_Data(file_handle, "CellData");
    Write_DataArray(file_handle, trim(dataname), 1, CellData);
    Close_Data(file_handle, "CellData");

    // Close file
    Close_vtr(file_handle);

    return; 
};

// -------------------------------------------------------------------------------------- //
template <class T>
void UCartMesh::Export_PointData_vtr(string      filename,
        string      dataname,
        vector< T > &PointData) {

    // ====================================================================================== //
    // template <class T>                                                                     //
    // void UCartMesh::Export_PointData_vtr(string       filename,                    //
    //                                              string       dataname,                    //
    //                                              vector< T > &PointData)                   //
    //                                                                                        //
    // Export point data into a .vtr file.                                                    //
    // ====================================================================================== //
    // INPUT                                                                                  //
    // ====================================================================================== //
    // - filename    : string, .vtr file name                                                 //
    // - dataname    : string, data name                                                      //
    // - PointData   : [(nc[0]+1)*(nc[1]+1)*(nc[2]+1)-by-1] vector< T > with point data                //
    // ====================================================================================== //
    // OUTPUT                                                                                 //
    // ====================================================================================== //
    // - none                                                                                 //
    // ====================================================================================== //

    // ====================================================================================== //
    // VARIABLES DECLARATION                                                                  //
    // ====================================================================================== //

    // Local variables
    ofstream        file_handle;

    // Counters

    // ====================================================================================== //
    // EXPORT CELL DATA                                                                       //
    // ====================================================================================== //

    // Open file
    Open_vtr(file_handle, trim(filename));

    // Write mesh data
    Write_vtrMeshData(file_handle, nc[0], nc[1], nc[2], edge[0], edge[1], edge[2]);

    // Write cell data
    Open_Data(file_handle, "PointData");
    Write_DataArray(file_handle, trim(dataname), 1, PointData);
    Close_Data(file_handle, "PointData");

    // Close file
    Close_vtr(file_handle);

    return; 

};

