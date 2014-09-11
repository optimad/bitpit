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
      // - Tensor_Product: dyadic product                                                    //
      // - sign          : sign function                                                     //
      // - Mat_Mul       : Matrix Matrix multiplication                                      //
      // - Dia_Mat_Mul   : Diagonal Matrix Matrix multiplication                             //
      // - Mat_Mul_Dia   : Matrix Diagonal Matrix multiplication                             //
      // - Transpose     : Matrix transposition                                              //
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
      #ifndef __MATHFUNCTIONSTEMP_HH__
      #define __MATHFUNCTIONSTEMP_HH__

      // =================================================================================== //
      // INCLUDES                                                                            //
      // =================================================================================== //
      # include <cmath>
      # include <vector>

      // =================================================================================== //
      // NAME SPACES                                                                         //
      // =================================================================================== //
      using namespace std;

      // =================================================================================== //
      // FUNCTION PROTOTYPES                                                                 //
      // =================================================================================== //
      // Operator Tensor_Product ------------------------------------------------------------ //
      template <class T>
      vector<vector<T>> Tensor_Product(vector<T> const       &,
                                       vector<T> const       &);

      template <class T, size_t n, size_t m>
      array<array<T,m>, n> Tensor_Product(array<T, n> const   &,
                                          array<T, m> const   &);

      // Matrix Vector Multiplication ------------------------------------------------------ //
      template <class T>
      vector<T> Mat_Mul( const vector< vector<T>>  &, 
                         const vector<T>           &);

      template <class T, size_t d1, size_t d2>
      array<T, d1> Mat_Mul( const array< array<T, d2>, d1> &, 
                            const array<T, d2>             &);


      // Matrix Matrix Multiplication ------------------------------------------------------ //
      template <class T>
      vector< vector<T> > Mat_Mul( const vector< vector<T> > &, 
                                   const vector< vector<T> > &);

      template <class T>
      vector< vector<T> > Dia_Mat_Mul( const vector<T>           &, 
                                       const vector< vector<T> > &);

      template <class T>
      vector< vector<T> > Dia_Mat_Mul( const vector< vector<T> > &, 
                                       const vector<T>           &);
      
      template <class T, size_t d1, size_t d2, size_t d3>
      array< array<T, d2> , d1> Mat_Mul( const array< array<T, d3>, d1>     &, 
                                         const array< array<T, d2>, d3>     &) ;

      template <class T, size_t d1, size_t d2>
      array< array<T, d2> , d1> Dia_Mat_Mul( const array<T, d1>             & , 
                                             const array< array<T, d2>, d1> & );

      template <class T, size_t d1, size_t d2>
      array< array<T, d2> , d1> Dia_Mat_Mul( const array< array<T, d2>, d1> & , 
                                             const array<T, d2>             & );


      // Matrix transposition -------------------------------------------------------------- //
      template <class T>
      vector < vector<T>> Transpose( const vector< vector<T> > & );

      template <class T, size_t m, size_t n>
      array < array<T, n>, m> Transpose( const array< array<T, m>, n> & );

      // =================================================================================== //
      // TEMPLATES                                                                           //
      // =================================================================================== //
      # include "MathFunctTemp.tpp"

      #endif
