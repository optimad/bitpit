        // =================================================================================== //
        //                         Class_SurfTri templates                                     //
        //                                                                                     //
        // Definitions of Class_SurfTri template methods.                                      //
        // =================================================================================== //
        // INFO                                                                                //
        // =================================================================================== //
        // Author   : Alessandro Alaia                                                         //
        // Company  : Optimad Engineering srl                                                  //
        // Date     : Sept 12, 2013                                                            //
        // Version  : v1.0                                                                     //
        // =================================================================================== //

        // =================================================================================== //
        // INTERPOLATION TOOLS                                                                 //
        // =================================================================================== //

        // Point data ======================================================================== //

        // ----------------------------------------------------------------------------------- //
        template <class T>
        void Class_SurfTri::PointData2CellData(
            vector<T>       &PData,
            vector<T>       &CData
        ) {

        // =================================================================================== //
        // template <class T>                                                                  //
        // void Class_SurfTri::PointData2CellData(                                             //
        //     vector<T>       &PData,                                                         //
        //     vector<T>       &CData)                                                         //
        //                                                                                     //
        // Convert point data into cell data, using linear interpolation                       //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - PData      : vector<T>, with point data at tasselation vertices                   //
        // - CData      : vector<T>, with cell  data at cell centers                           //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - none                                                                              //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        int             nedge;
        darray3E        xC, n, t;
        dvector1D       coeff(3, 0.0);

        // Counters
        int             A;

        // =================================================================================== //
        // RESHAPE INPUT VARIABLES                                                             //
        // =================================================================================== //

        // Cell Data
        CData.resize(nSimplex, 0.0);

        // Normals
        if (Normal.size() < nSimplex) {
            GenerateNormals();
        }

        // =================================================================================== //
        // CONVERT POINT DATA INTO CELL DATA.                                                  //
        // =================================================================================== //
        for (A = 0; A < nSimplex; A++) {

            // Simplex baricenter ------------------------------------------------------------ //
            xC = Baricenter(A);
            nedge = Simplex[A].size();

            // Compute value a simplex baricenter -------------------------------------------- //

            // Degenerate case (point)
            if (nedge == 1) {
                CData[A] = PData[Simplex[A][0]];
            }

            // Segment
            else if (nedge == 2) {

                // Scope variables
                int             V0, V1;
                double          d, xi;

                // Vertex global index
                V0 = Simplex[A][0];
                V1 = Simplex[A][1];
                d = norm_2(Vertex[V1] - Vertex[V0]);
                xi = norm_2(xC - Vertex[V0]);
                CData[A] = (1.0 - xi) * PData[V0] + xi * PData[V1];
            }

            // Triangle
            else if (nedge == 3) {

                // Scope variables
                int              V0, V1, V2;
                double           xi1, xi2, tau2;

                // Local reference frame
                V0 = Simplex[A][0];
                V1 = Simplex[A][1];
                V2 = Simplex[A][2];
                n = Vertex[V1] - Vertex[V0];
                xi1 = norm_2(n);
                if (xi1 < 1.0e-12) {
                    cout << "error" << endl;
                    return;
                }
                n = n/xi1;
                t = Cross_Product(Normal[A], n);

                // Vertex coordinates in local ref. frame
                xi2 = Dot_Product(Vertex[V2] - Vertex[V0], n);
                tau2 = Dot_Product(Vertex[V2] - Vertex[V0], t);
                if (abs(tau2) < 1.0e-12) {
                    cout << "error" << endl;
                    return;
                }

                // Interpolation coeff.
                coeff[0] = PData[V0];
                coeff[1] = (PData[V1] - PData[V0])/xi1;
                coeff[2] = (PData[V2] - PData[V0])/tau2 - coeff[1]*xi2/tau2;

                // Interpolated value
                CData[A] = Dot_Product(xC - Vertex[V0], n) * coeff[1]
                         + Dot_Product(xC - Vertex[V0], t) * coeff[2]
                         + coeff[0];

            }

            // Generic Simplex
            else {
            }

        } //next T

        return; };

        // ----------------------------------------------------------------------------------- //
        template <class T>
        void Class_SurfTri::RemapPointData(
            ivector1D &map,
            vector<T> &PointData
        ) {

        // =================================================================================== //
        // template <class T>                                                                  //
        // void Class_SurfTri::RemapPointData(                                                 //
        //     ivector1D &map,                                                                 //
        //     vector<T> &PointData)                                                           //
        //                                                                                     //
        //                                                                                     //
        // Remap point data onto the new vertex set.                                           //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - map        : [nVertex-by-1] ivector1D with old -> new vertex map                  //
        // - PointData  : [nVertex-by-1] vector<T> with copy-constructible point data to be    //
        //                remapped.                                                            //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - none                                                                              //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        // none

        // Counters
        int      i;

        // =================================================================================== //
        // REMAP POINT DATA                                                                    //
        // =================================================================================== //
        for (i = 0; i < map.size(); i++) {
            if (map[i] >= 0) {
                PointData[map[i]] = PointData[i];
            }
        } //next i

        // =================================================================================== //
        // RESIZE OUTPUT VARIABLE                                                              //
        // =================================================================================== //
        PointData.resize(nVertex);

        return; }

        // ----------------------------------------------------------------------------------- //
        template <class T>
        void Class_SurfTri::InterpolatePointData(
            darray3E         &P,
            int              seed,
            vector<T>        &PData,
            T                &data
        ) {

        // =================================================================================== //
        // template <class T>                                                                  //
        // void Class_SurfTri::InterpolatePointData(                                           //
        //     array3E          &P,                                                            //
        //     int              seed,                                                         //
        //     vector<T>        &PData,                                                        //
        //     T                &data)                                                         //
        //                                                                                     //
        // Interpolate point data defined on surface tasselation.                              //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - P        : darray3E, point coordinates                                            //
        // - seed     : int, seed for simplex search algorithm                                 //
        // - PData    : vector<T>, with point data                                             //
        // - data     : T, with interpolated value                                             //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - none                                                                              //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        int                  S, n;

        // Counters
        // none

        // =================================================================================== //
        // FIND SIMPLEX CONTAINING POINT P                                                     //
        // =================================================================================== //
        S = ReturnTriangleID(P, seed);
        if (S < 0) { return; };

        // =================================================================================== //
        // INTERPOLATE VALUE                                                                   //
        // =================================================================================== //
        n = Simplex[S].size();

        // Point simplex --------------------------------------------------------------------- //
        if (n == 1) {

            // Interpolate data
            data = PData[S];
        }

        // Segment simplex ------------------------------------------------------------------- //
        else if (n == 2) {

            // Scope variables
            double  xi, eta;

            // Interpolation weights (local coordinate system)
            eta = Area(S);
            xi = norm_2(P - Vertex[Simplex[S][0]])/eta;

            // Interpolate point data
            data = (1.0 - xi) * PData[Simplex[S][0]] + xi * PData[Simplex[S][1]];

        }

        // General simplex ------------------------------------------------------------------- //
        else {

            // Scope variables
            int                  i;
            double               xi, eta, scale;
            array< double, 3 >   t, n, N;
//             array< double, 3 >   coeff;
            vector<T>            coeff(3);

            // Local coordinate system
//ht            for (i = 0; i < dim; i++) {
//ht                t[i] = Vertex[Simplex[S][1]][i] - Vertex[Simplex[S][0]][i];
//ht                N[i] = Normal[S][i];
//ht            } //next i

            t = Vertex[Simplex[S][1]] - Vertex[Simplex[S][0]];
            N = Normal[S];

            scale = norm_2(t);
            t = t/scale;
            n = Cross_Product(N, t);
            n = n/norm_2(n);

            // Interpolation coeff.
//ht            for (i = 0; i < dim; i++) {
//ht                N[i] = Vertex[Simplex[S][2]][i] - Vertex[Simplex[S][0]][i];
//ht            } //next i

            N = Vertex[Simplex[S][2]] - Vertex[Simplex[S][0]];
            N = N/scale;
            xi  = Dot_Product(N, t);
            eta = Dot_Product(N, n);
            coeff[0] = PData[Simplex[S][0]];
            coeff[1] = PData[Simplex[S][1]] - coeff[0];
            coeff[2] = (PData[Simplex[S][2]] - coeff[1]*xi - coeff[0])/eta;

            // Interpolation point
//ht            for (i = 0; i < dim; i++) {
//ht                N[i] = P[i] - Vertex[Simplex[S][0]][i];
//ht            } //next i

            N = P[i] - Vertex[Simplex[S][0]];
            N = N/scale;
            xi  = Dot_Product(N, t);
            eta = Dot_Product(N, n);

            // Interpolate data
            data = coeff[0] + coeff[1]*xi + coeff[2]*eta;
        }

        return; };

        // Cell data ========================================================================= //

        // ----------------------------------------------------------------------------------- //
        template <class T>
        void Class_SurfTri::CellData2PointData(
            vector<T>       &CData,
            vector<T>       &PData
        ) {

        // =================================================================================== //
        // template <class T>                                                                  //
        // void Class_SurfTri::CellData2PointData(                                             //
        //     vector<T>       &CData,                                                         //
        //     vector<T>       &PData)                                                         //
        //                                                                                     //
        // Convert cell data into point data.                                                  //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - CData        : vector<T>, vector with cell data                                   //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - PData        : vector<T>, vector with point data                                  //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        bool                flag_open;
        int                 dim = Vertex[0].size();
        double              theta;
        dvector1D           theta_tot(nVertex, 0.0);
        array<double, 3>    u, v;

        // Counters
        int                 i, j, n, R, V;

        // =================================================================================== //
        // RESIZE INPUT/OUTPUT VARIABLES                                                       //
        // =================================================================================== //
        PData.resize(nVertex, (T) 0.0);

        // =================================================================================== //
        // CONVERT CELL DATA INTO POINT DATA                                                   //
        // =================================================================================== //
        for (R = 0; R < nSimplex; R++) {
            n = Simplex[R].size();
            for (i = 0; i < n; i++) {

                // Vertex global index ------------------------------------------------------- //
                V = Simplex[R][i];

                // Point simplex ------------------------------------------------------------- //
                if (n == 1) {

                    // Remap cell data into point data
                    PData[V] = CData[R];

                }

                // Segment simplex ----------------------------------------------------------- //
                else if (n == 2) {

                    PData[V] += CData[R];
                    theta_tot[V] += 1.0;

                }

                // General simplex ----------------------------------------------------------- //
                else {

                    u = Vertex[Simplex[R][(i+1) % n]] - Vertex[Simplex[R][i]];
                    u = u/norm_2(u);
                    v = Vertex[Simplex[R][(i+2) % n]] - Vertex[Simplex[R][i]];
                    v = v/norm_2(v);
                    theta = acos(min(1.0, max(-1.0, Dot_Product(u, v))));
                    PData[V] += ((T) theta) * CData[R];
                    theta_tot[V] += theta;

                }
            } //next i
        } //next T

        // =================================================================================== //
        // RESCALE POINT DATA                                                                  //
        // =================================================================================== //
        for (R = 0; R < nVertex; R++) {
            PData[R] = PData[R]/((T) theta_tot[R]);
        } //next T

        return; };

        // ----------------------------------------------------------------------------------- //
        template <class T>
        void Class_SurfTri::RemapCellData(ivector1D &map, vector<T> &CellData) {

        // =================================================================================== //
        // template <class T>                                                                  //
        // void Class_SurfTri::RemapCellData(ivector1D &map, vector<T> CellData)               //
        //                                                                                     //
        // Remap point data onto the new vertex set.                                           //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - map        : [nSimplex-by-1] ivector1D with old -> new simplex map                //
        // - CellData   : [nSimplex-by-1] vector<T> with copy-constructible cell data to be    //
        //                remapped.                                                            //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - none                                                                              //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        // none

        // Counters
        int      I_;

        // =================================================================================== //
        // REMAP POINT DATA                                                                    //
        // =================================================================================== //
        for (I_ = 0; I_ < map.size(); I_++) {
            if (map[I_] >= 0) {
                CellData[map[I_]] = CellData[I_];
            }
        } //next I_

        // =================================================================================== //
        // RESIZE OUTPUT VARIABLE                                                              //
        // =================================================================================== //
        CellData.resize(nSimplex);

        return; }

        // =================================================================================== //
        // I_/O FUNCTIONS                                                                       //
        // =================================================================================== //

        // ----------------------------------------------------------------------------------- //
        template <class T>
        void Class_SurfTri::ExportCellData_vtu(string                 filename,
                                               svector1D             &data_title,
                                               vector< vector< T > > &data) {

        // =================================================================================== //
        // template <class T>                                                                  //
        // void Class_SurfTri::ExportCellData_vut(string                 filename,             //
        //                                        svector1D             &data_title,           //
        //                                        vector< vector< T > > &data)                 //
        //                                                                                     //
        // Export cell data in a .vtu file.                                                    //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - filename   : string, with .vtu filename                                           //
        // - data_title : svector1D, with data titles                                          //
        // - data       : dvector2D, with data.                                                //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - none                                                                              //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        int                  n, dum, connect_size;
        vector<short int>    types(nSimplex, 0);
        ivector1D            type_lib(6, 0);
        ivector1D            offsets(nSimplex, 0);
        ivector1D            connectivity;
        dvector1D            vertex(3*nVertex, 0.0);
        vector<T>            out(nSimplex, (T) 0.0);
        ofstream             file_handle;

        // Counters
        int                  i, j, k;

        // =================================================================================== //
        // PREPARE OUTPUT                                                                      //
        // =================================================================================== //

        // Element types:
        type_lib[0] = 0;
        type_lib[1] = 1;
        type_lib[2] = 3;
        type_lib[3] = 5;
        type_lib[4] = 9;
        type_lib[5] = 7;

        // Vertex coordinate list ------------------------------------------------------------ //
        k = 0;
        for (i = 0; i < nVertex; i++) {
            for (j = 0; j < 3; j++) {
                vertex[k] = Vertex[i][j];
                k++;
            } //next j
        } //next i

        // Simplex --------------------------------------------------------------------------- //

            // offsets and element types
            connect_size = 0;
            for (i = 0; i < nSimplex; i++) {
                n = Simplex[i].size();
                connect_size += n;
                dum = type_lib[min(n,5)];
                types[i] = dum;
                offsets[i] = connect_size;
            } //next i

            // Resize connectivity array
            connectivity.resize(connect_size, 0);

            // Simplex-vertex connectivity>
            k = 0;
            for (i = 0; i < nSimplex; i++) {
                n = Simplex[i].size();
                for (j = 0; j < n; j++) {
                    connectivity[k] = Simplex[i][j];
                    k++;
                } //next j
            } //next i

        // =================================================================================== //
        // EXPORT DATA IN A .VTU FILE                                                          //
        // =================================================================================== //

        // Open .vtu file
        Open_vtu(file_handle, trim(filename));

        // Export mesh data
        Write_vtuMeshData(file_handle, nVertex, nSimplex, vertex, connectivity, offsets, types);

        // Export data
        Open_Data(file_handle, "CellData");
        for (i = 0; i < data[0].size(); i++) {
            k = 0;
            for (j = 0; j < nSimplex; j++) {
                out[k] = data[j][i];
                k++;
            } //next j
            Write_DataArray(file_handle, data_title[i], 1, out);
        } //next i
        Close_Data(file_handle, "CellData");

        // Close .vtu file
        Close_vtu(file_handle);

        return; };

        // ----------------------------------------------------------------------------------- //
        template <class T>
        void Class_SurfTri::ExportPointData_vtu(string                 filename,
                                                svector1D             &data_title,
                                                vector< vector< T > > &data) {

        // =================================================================================== //
        // template <class T>                                                                  //
        // void Class_SurfTri::ExportPointData_vut(string                 filename,            //
        //                                         svector1D             &data_title,          //
        //                                         vector< vector< T > > &data)                //
        //                                                                                     //
        // Export point data in a .vtu file.                                                   //
        // =================================================================================== //
        // INPUT                                                                               //
        // =================================================================================== //
        // - filename   : string, with .vtu filename                                           //
        // - data_title : svector1D, with data titles                                          //
        // - data       : dvector2D, with data.                                                //
        // =================================================================================== //
        // OUTPUT                                                                              //
        // =================================================================================== //
        // - none                                                                              //
        // =================================================================================== //

        // =================================================================================== //
        // VARIABLES DECLARATION                                                               //
        // =================================================================================== //

        // Local variables
        int                  n, dum, connect_size;
        vector<short int>    types(nSimplex, 0);
        ivector1D            type_lib(6, 0);
        ivector1D            offsets(nSimplex, 0);
        ivector1D            connectivity;
        dvector1D            vertex(3*nVertex, 0.0);
        vector<T>            out(nVertex, (T) 0.0);
        ofstream             file_handle;

        // Counters
        int                  i, j, k;

        // =================================================================================== //
        // PREPARE OUTPUT                                                                      //
        // =================================================================================== //

        // Element types:
        type_lib[0] = 0;
        type_lib[1] = 1;
        type_lib[2] = 3;
        type_lib[3] = 5;
        type_lib[4] = 9;
        type_lib[5] = 7;

        // Vertex coordinate list ------------------------------------------------------------ //
        k = 0;
        for (i = 0; i < nVertex; i++) {
            for (j = 0; j < 3; j++) {
                vertex[k] = Vertex[i][j];
                k++;
            } //next j
        } //next i

        // Simplex --------------------------------------------------------------------------- //

            // offsets and element types
            connect_size = 0;
            for (i = 0; i < nSimplex; i++) {
                n = Simplex[i].size();
                connect_size += n;
                dum = type_lib[min(n,5)];
                types[i] = dum;
                offsets[i] = connect_size;
            } //next i

            // Resize connectivity array
            connectivity.resize(connect_size, 0);

            // Simplex-vertex connectivity>
            k = 0;
            for (i = 0; i < nSimplex; i++) {
                n = Simplex[i].size();
                for (j = 0; j < n; j++) {
                    connectivity[k] = Simplex[i][j];
                    k++;
                } //next j
            } //next i

        // =================================================================================== //
        // EXPORT DATA IN A .VTU FILE                                                          //
        // =================================================================================== //

        // Open .vtu file
        Open_vtu(file_handle, trim(filename));

        // Export mesh data
        Write_vtuMeshData(file_handle, nVertex, nSimplex, vertex, connectivity, offsets, types);

        // Export data
        Open_Data(file_handle, "PointData");
        for (i = 0; i < data[0].size(); i++) {
            k = 0;
            for (j = 0; j < nVertex; j++) {
                out[k] = data[j][i];
                k++;
            } //next j
            Write_DataArray(file_handle, data_title[i], 1, out);
        } //next i
        Close_Data(file_handle, "PointData");

        // Close .vtu file
        Close_vtu(file_handle);

        return; };

