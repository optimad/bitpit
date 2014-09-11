      // =================================================================================== //
      //                          BASIC MATH FUNCTIONS FOR VECTORS                           //
      //                                                                                     //
      // Functions definitions.                                                              //
      //                                                                                     //
      // LIST OF FUNCTIONS                                                                   //
      // - minval        : min element on vector's column                                    //
      // - maxval        : max element on vector's column                                    //
      // - sum           : vector column-wise sum                                            //
      // - norm_1        : vector 1-norm                                                     //
      // - norm_2        : vector Euclidean norm                                             //
      // - norm_inf      : vector inf-norm                                                   //
      // - Dot_Product   : vector dot product                                                //
      // - Cross_Product : vector cross product                                              //
      // - Dyadic        : dyadic product                                                    //
      // - sign          : sign function                                                     //
      // =================================================================================== //
      // INFO                                                                                //
      // =================================================================================== //
      // Author     : Alessandro Alaia                                                       //
      // Company    : Optimad srl                                                            //
      // Date       : Jul 16, 2013                                                           //
      // Version    : v1.0                                                                   //
      //                                                                                     //
      // All rights reserved.                                                                //
      // =================================================================================== //

      // =================================================================================== //
      // TEMPLATE IMPLEMENTATIONS                                                            //
      // =================================================================================== //

      // ----------------------------------------------------------------------------------- //
      template <class T>
      vector<vector<T>> Tensor_Product(const vector<T> &x, const vector<T> &y){

      //                                                                                     //
      // Compute the tensor product between x, and y.                                         //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - x,     : vector<T>, 1st argument of the tensor product                            //
      // - y,     : vector<T>, 2nd argument of the tensor product                            //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - z      : vector<vector<T>>, tensor product between x and y                        //
      // =================================================================================== //

      int  i, j;
      int  n = x.size(); 
      int  m = y.size(); 
      vector<T>      row(m,0.0);
      vector<vector<T>> z(n,row) ;

      for( i=0; i<n; i++){
        for( j=0; j<m; j++){
          z[i][j] = x[i] *y[j] ;
        };
      };

      return (z);}

      // ----------------------------------------------------------------------------------- //
      template <class T, size_t n, size_t m>
      array<array<T,m>,n> Tensor_Product(const array<T,n> &x, const array<T,m> &y){

      //                                                                                     //
      // Compute the tensor product between x, and y.                                         //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - x,     : array<n,T>, 1st argument of the tensor product                            //
      // - y,     : array<m,T>, 2nd argument of the tensor product                            //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - z      : array<n,array<m,T>>, tensor product between x and y                        //
      // =================================================================================== //

      int  i, j;
      array<array<T,m>,n> z ;

      for( i=0; i<n; i++){
        for( j=0; j<m; j++){
          z[i][j] = x[i] *y[j] ;
        };
      };

      return (z);}

      // Matrix Vector Multiplication ------------------------------------------------------ //
      template <class T>
      vector<T> Mat_Mul( const vector< vector<T>> &M, const vector<T> &x){

      // =================================================================================== //
      // Perform Matrix Vector multiplication between M and x                                //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - M,     : vector<T, vector< T> >, Matrix                                         //
      // - x,     : vector<T> Vector                                                      //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - z      : vector<T>                                                               //
      // =================================================================================== //

        int d1 = M.size();
        int d2 = M[0].size();

        vector<T>      z(d1,0.0);

        for( int i=0; i<d1; i++){
          z[i]= Dot_Product( M[i], x );
        }

        return (z);
      };

      template <class T, size_t d1, size_t d2>
      array<T, d1> Mat_Mul( const array< array<T, d2>, d1> &M, const array<T, d2> &x){

      // =================================================================================== //
      // template <class T, size_t d>                                                        //
      // array<T, d> Mat_Mul(array<T,array<T, d> d> const &M, array<T, d> const &x);         //
      //                                                                                     //
      // Perform Matrix Vector multiplication between M and x                                //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - M,     : array<T, array< T, d2> d1>, Matrix                                       //
      // - x,     : array<T, d2>, Vector                                                     //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - z      : array<T, d1>,  product between x and y                                   //
      // =================================================================================== //

        array<T, d1>      z;

        for( int i=0; i<d1; i++){
          z[i]= Dot_Product( M[i], x);
        }

        return (z);
      };

      // Matrix Matrix Multiplication ------------------------------------------------------ //
      template <class T>
      vector< vector<T> > Mat_Mul( const vector< vector<T> > &M, const vector<vector<T> > &N){
      // =================================================================================== //
      // Perform Matrix Matrix multiplication between M and N                                //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - M,     : vector<T, vector< T> >, Matrix   d1 x d3                                 //
      // - N,     : vector<T, vector< T> >, Matrix   d3 x d2                               //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - Q,     : vector< vector< T> >, Matrix                                         //
      // =================================================================================== //

        int i, j;
        int d1, d2, d3;

        d1= M.size();
        d2= N[0].size() ;
        d3= N.size() ;

        vector< vector<T> > Q(d1, vector<double> (d2,0.0) );
        vector< vector<T> > Tr;

        Tr = Transpose( N ) ;

        for( i=0; i<d1; i++){
          for( j=0; j<d2; j++){
            Q[i][j]= Dot_Product( M[i], Tr[j] );
          };
        };

        return (Q);
      };

      template <class T>
      vector< vector<T> > Dia_Mat_Mul( const vector<T> &M, const vector<vector<T> > &N){
      // =================================================================================== //
      // Perform Matrix Matrix multiplication between M and N                                //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - M,     : vector<T>,              Diagonal Matrix, d1 elements in sparse storage   //
      // - N,     : vector<T, vector< T> >, Matrix  d1 x d2                                  //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - Q,     : vector< vector< T> >, Matrix    d1 x d2                                  //
      // =================================================================================== //

        int i;
        int d1, d2;

        d1= M.size();
        d2= N[0].size() ;

        vector< vector<T> > Q( N );

        for( i=0; i<d1; i++){
          Q[i] = N[i] * M[i] ;
        };

        return (Q);
      };

      template <class T>
      vector< vector<T> > Dia_Mat_Mul( const vector< vector<T> > &M, const vector<T> &N){
      // =================================================================================== //
      // Perform Matrix Matrix multiplication between M and N                                //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - M,     : vector<T, vector< T> >, Matrix  d1 x d2                                  //
      // - N,     : vector<T>,              Diagonal Matrix, d2 elements in sparse storage   //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - Q,     : vector< vector< T> >, Matrix    d1 x d2                                  //
      // =================================================================================== //

        int i, j;
        int d1;

        vector< vector<T> > Q( M );

        d1= M.size() ;

        for( i=0; i<d1; i++ ){
          Q[i] = M[i] * N ;
        };

        return (Q);

      };

      template <class T, size_t d1, size_t d2, size_t d3>
      array< array<T, d2> , d1> Mat_Mul( const array< array<T, d3>, d1> &M, const array<array<T, d2>, d3> &N){
      // =================================================================================== //
      // Perform Matrix Matrix multiplication between M and N                                //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - M,     : array<T, array< T, d1> d3>, Matrix                                         //
      // - N,     : array<T, array< T, d3> d2>, Matrix                                         //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - Q,     : array<T, array< T, d1> d2>, Matrix                                         //
      // =================================================================================== //

        int i, j;

        array< array<T, d2> , d1> Q;
        array< array<T, d2> , d3> Tr;

        Tr = Transpose( N ) ;

        for( i=0; i<d1; i++){
          for( j=0; j<d2; j++){
            Q[i][j]= Dot_Product( M[i], Tr[j] );
          };
        };

        return (Q);
      };

      template <class T, size_t d1, size_t d2>
      array< array<T, d2> , d1> Dia_Mat_Mul( const array< T, d1> &M, const array<array<T, d2>, d1> &N){
      // =================================================================================== //
      // Perform Matrix Matrix multiplication between M and N                                //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - M,     : array<T, d1>, Diagonal Matrix in sparse storage                          //
      // - N,     : array<T, array< T, d2> d1>, Matrix                                       //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - Q,     : array<T, array< T, d2> d1>, Matrix                                       //
      // =================================================================================== //

        int i;
        array< array<T, d2> , d1> Q;

        for( i=0; i<d1; i++){
          Q[i] = M[i] *N[i] ;
        };

        return (Q);
      };

      template <class T, size_t d1, size_t d2>
      array< array<T, d2> , d1> Dia_Mat_Mul( const array<array<T, d2>, d1> &M, const array< T, d2> &N ){
      // =================================================================================== //
      // Perform Matrix Matrix multiplication between M and N                                //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - M,     : array<T, array< T, d2> d1>, Matrix                                       //
      // - N,     : array<T, d1>, Diagonal Matrix in sparse storage                          //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - Q,     : array<T, array< T, d2> d1>, Matrix                                       //
      // =================================================================================== //

        int i;
        array< array<T, d2> , d1> Q;

        for( i=0; i<d1; i++){
          Q[i] = M[i] *N ;
        };

        return (Q);
      };

      // Matrix Transposition -------------------------------------------------------------- //
      template <class T>
      vector < vector<T> > Transpose( const vector< vector<T>> &M){

      // =================================================================================== //
      // template <class T>                                                                  //
      // vector < vector<T> > Transpose( const vector< vector<T>> &M)                        //
      //                                                                                     //
      // Perform Matrix Transposition                                                        //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - M,     : vector< vector< T> >      Matrix                                         //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - z      :  vector< vector< T> >     Transposed Metrix                              //
      // =================================================================================== //

        int n = M.size();
        int m = M[0].size() ;
        vector< vector<T> >   z(m, vector<T>(n, 0.0));

        for( int i=0; i<m; i++){
          for(int j=0; j<n; j++){
            z[i][j] = M[j][i];
          };
        };

        return (z);
      };

      template <class T, size_t m, size_t n>
      array< array<T, n>, m> Transpose( const array< array<T, m>, n> &M){

      // =================================================================================== //
      // template <class T, size_t m, size_t n>                                              //
      // array < array<T, n>, m> Transpose( array< array<T, m>, n> const &M)nst &x);         //
      //                                                                                     //
      // Perform Matrix Transposition                                                        //
      // =================================================================================== //
      // INPUT                                                                               //
      // =================================================================================== //
      // - M,     : array<T, array< T, m> n>, Matrix                                         //
      // =================================================================================== //
      // OUTPUT                                                                              //
      // =================================================================================== //
      // - z      : array < array<T, n>, m>   Transposed Metrix                              //
      // =================================================================================== //
  
        array< array<T, n> ,m>      z;
  
        for( int i=0; i<m; i++){
          for(int j=0; j<n; j++){
            z[i][j] = M[j][i];
          };
        };
  
        return (z);
      };
