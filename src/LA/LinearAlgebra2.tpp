// =================================================================================== //
// TEMPLATE IMPLEMENTATIONS                                                            //
// =================================================================================== //

// ----------------------------------------------------------------------------------- //
/*!
    Tensor product.

    \param[in] x 1st argument of the tensor product
    \param[in] y 2nd argument of the tensor product

    \result tensor product between x and y
*/
template <class T>
std::vector<std::vector<T>> Tensor_Product(
    const std::vector<T>                        &x,
    const std::vector<T>                        &y
) {

    int  i, j;
    int  n = x.size(); 
    int  m = y.size(); 
    std::vector<T>      row(m,0.0);
    std::vector<std::vector<T>> z(n,row) ;

    for( i=0; i<n; i++){
        for( j=0; j<m; j++){
            z[i][j] = x[i] *y[j] ;
        };
    };

return (z);}

// ----------------------------------------------------------------------------------- //
/*!
    Tensor product. Overloading of Tensor_Product() for container array.

    \param[in] x 1st argument of the tensor product
    \param[in] y 2nd argument of the tensor product

    \result tensor product between x and y
*/
template <class T, size_t n, size_t m>
std::array<std::array<T,m>,n> Tensor_Product(
    const std::array<T,n>                       &x,
    const std::array<T,m>                       &y
) {
    int  i, j;
    std::array<std::array<T,m>,n> z ;

    for( i=0; i<n; i++){
        for( j=0; j<m; j++){
            z[i][j] = x[i] *y[j] ;
        };
    };

    return (z);
}

// Matrix Vector Multiplication ====================================================== //

// ----------------------------------------------------------------------------------- //
/*!
    Matrix-vector multiplication.

    \param[in] M input matrix
    \param[in] x input vector

    \result product between M and x
*/
template <class T>
std::vector<T> Mat_Mul(
    const std::vector< std::vector<T>>          &M,
    const std::vector<T>                        &x
) {

    int d1 = M.size();
    int d2 = M[0].size();

    std::vector<T>      z(d1,0.0);

    for( int i=0; i<d1; i++){
        z[i]= Dot_Product( M[i], x );
    }

    return (z);
};

// ----------------------------------------------------------------------------------- //
/*!
    Matrix-vector multiplication. Overloading of Mat_Mul() for container array.

    \param[in] M input matrix
    \param[in] x input vector

    \result product between M and x
*/
template <class T, size_t d1, size_t d2>
std::array<T, d1> Mat_Mul(
    const std::array< std::array<T, d2>, d1>    &M,
    const std::array<T, d2>                     &x
) {

    std::array<T, d1>      z;

    for( int i=0; i<d1; i++){
        z[i]= Dot_Product( M[i], x);
    }

    return (z);
};

// Matrix Matrix Multiplication ====================================================== //

// ----------------------------------------------------------------------------------- //
/*!
    Matrix product.

    \param[in] M 1st argument of matrix multiplication
    \param[in] N 2nd argument of matrix multiplication

    \result product between M and N.
*/
template <class T>
std::vector< std::vector<T> > Mat_Mul(
    const std::vector< std::vector<T> >         &M,
    const std::vector<std::vector<T> >          &N
) {

    int i, j;
    int d1, d2, d3;

    d1= M.size();
    d2= N[0].size() ;
    d3= N.size() ;

    std::vector< std::vector<T> > Q(d1, std::vector<double> (d2,0.0) );
    std::vector< std::vector<T> > Tr;

    Tr = Transpose( N ) ;

    for( i=0; i<d1; i++){
        for( j=0; j<d2; j++){
            Q[i][j]= Dot_Product( M[i], Tr[j] );
        };
    };

    return (Q);
};

// ----------------------------------------------------------------------------------- //
/*!
    Diadic matrix multiplicationx

    \param[in] M 1st argument
    \param[in] N 2nd argument

    \result diadic product between M and N.
*/
template <class T>
std::vector< std::vector<T> > Dia_Mat_Mul(
    const std::vector<T>                        &M,
    const std::vector<std::vector<T> >          &N
) {

    int i;
    int d1, d2;

    d1= M.size();
    d2= N[0].size() ;

    std::vector< std::vector<T> > Q( N );

    for( i=0; i<d1; i++){
        Q[i] = N[i] * M[i] ;
    };

    return (Q);
};

// ----------------------------------------------------------------------------------- //
/*!
    Diadic matrix multiplication

    \param[in] M 1st argument
    \param[in] N 2nd argument

    \result diadic product between M and N.
*/
template <class T>
std::vector< std::vector<T> > Dia_Mat_Mul(
    const std::vector< std::vector<T> >         &M,
    const std::vector<T>                        &N
) {

    int i, j;
    int d1;

    std::vector< std::vector<T> > Q( M );

    d1= M.size() ;

    for( i=0; i<d1; i++ ){
        Q[i] = M[i] * N ;
    };

    return (Q);

};

// ----------------------------------------------------------------------------------- //
/*!
    Matrix multiplication. Overloading of Mat_Mul() function for container array.
    \param[in] M 1st argument
    \param[in] N 2nd argument

    \result product of M and N.
*/
template <class T, size_t d1, size_t d2, size_t d3>
std::array< std::array<T, d2> , d1> Mat_Mul(
    const std::array< std::array<T, d3>, d1>    &M,
    const std::array<std::array<T, d2>, d3>     &N
){
    int i, j;

    std::array< std::array<T, d2> , d1> Q;
    std::array< std::array<T, d2> , d3> Tr;

    Tr = Transpose( N ) ;

    for( i=0; i<d1; i++){
        for( j=0; j<d2; j++){
            Q[i][j]= Dot_Product( M[i], Tr[j] );
        };
    };

    return (Q);
};

// ----------------------------------------------------------------------------------- //
/*!
    Diadic matrix multiplication. Overloading of Dia_Mat_Mul() function for container array.

    \param[in] M 1st argument
    \param[in] N 2nd argument

    \result diadic product between M and N.
*/
template <class T, size_t d1, size_t d2>
std::array< std::array<T, d2> , d1> Dia_Mat_Mul(
    const std::array< T, d1>                    &M,
    const std::array<std::array<T, d2>, d1>     &N
){

    int i;
    std::array< std::array<T, d2> , d1> Q;

    for( i=0; i<d1; i++){
        Q[i] = M[i] *N[i] ;
    };

    return (Q);
};

// ----------------------------------------------------------------------------------- //
/*!
    Diadic matrix multiplication. Overloading of Dia_Mat_Mul() function for container array.

    \param[in] M 1st argument
    \param[in] N 2nd argument

    \result diadic product between M and N.
*/
template <class T, size_t d1, size_t d2>
std::array< std::array<T, d2> , d1> Dia_Mat_Mul(
    const std::array<std::array<T, d2>, d1>     &M,
    const std::array< T, d2>                    &N
) {

    int i;
    std::array< std::array<T, d2> , d1> Q;

    for( i=0; i<d1; i++){
        Q[i] = M[i] *N ;
    };

    return (Q);
};

// Matrix Transposition ============================================================== //

// ----------------------------------------------------------------------------------- //
/*!
    Compute matrix transpose.

    \param[in] M input matrix
    
    \result transpose of M
*/
template <class T>
std::vector < std::vector<T> > Transpose(
    const std::vector< std::vector<T>>          &M
) {

    int n = M.size();
    int m = M[0].size() ;
    std::vector< std::vector<T> >   z(m, std::vector<T>(n, 0.0));

    for( int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            z[i][j] = M[j][i];
        };
    };

    return (z);
};

// ----------------------------------------------------------------------------------- //
/*!
    Compute matrix transpose. Overloading of Transpose() function for container
    array.

    \param[in] M input matrix
    
    \result transpose of M
*/
template <class T, size_t m, size_t n>
std::array< std::array<T, n>, m> Transpose(
    const std::array< std::array<T, m>, n>      &M
) {

    std::array< std::array<T, n> ,m>      z;

    for( int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            z[i][j] = M[j][i];
        };
    };

    return (z);
};
