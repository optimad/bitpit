// ================================================================================== //
//                         OPERATORS - EXAMPLES OF USAGE -                            //
//                                                                                    //
// Examples of Operators usage.                                                       //
// ================================================================================== //
// INFO                                                                               //
// ================================================================================== //
// Author     : Alessandro Alaia                                                      //
// Version    : v3.0                                                                  //
//                                                                                    //
// All rights reserved.                                                               //
// ================================================================================== //

// ================================================================================== //
// INCLUDES                                                                           //
// ================================================================================== //
# include "Examples.hpp"

// ================================================================================== //
// IMPLEMENTATIONS                                                                    //
// ================================================================================== //

// ---------------------------------------------------------------------------------- //
void vectorOperators_Ex(
    void
) {

// ================================================================================== //
// void vectorOperators_Ex(                                                           //
//     void)                                                                          //
//                                                                                    //
// Examples of usage of operators for vectors                                         //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - none                                                                             //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - none                                                                             //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
string          selection;

// Counters
// none

// ================================================================================== //
// OUTPUT MESSAGE                                                                     //
// ================================================================================== //
{
    // Select sub case
    cout << "DEMO:  OPERATORS FOR STL VECTORS" << endl;
    cout << "  Select operators: " << endl;
    cout << "  ""+"" operator" << endl;
    cout << "  ""-"" operator" << endl;
    cout << "  ""*"" operator" << endl;
    cout << "  ""/"" operator" << endl;
    cin  >> selection;
}

// ================================================================================== //
// OPERATOR "+"                                                                       //
// ================================================================================== //
if (selection.compare("+") == 0) {

    // Scope variables
    double              a = 1, b = 2, c = 3;
    dvector1D           x, y, z;
    dvector2D           X, Y, Z;

    // Initialize vectors
    x.resize(3, 1.0);
    x[0] = 1.0;    x[1] = 2.0;    x[2] = 3.0;
    y = x;
    X.resize(3, dvector1D(2, 0.0));
    X[0][0] = 1.0;    X[0][1] = 2.0;
    X[1][0] = 3.0;    X[1][1] = 4.0;
    X[2][0] = 5.0;    X[2][1] = 6.0;
    Y = X;

    // Variables
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "x = ";
    cout << x << endl;
    cout << "y = ";
    cout << y << endl;
    cout << "X = ";
    cout << X << endl;
    cout << "Y = ";
    cout << Y << endl;

    // Sum between vectors
    cout << "sum between 1D vectors --------------------------- " << endl;
    cout << "x + y = ";
    z = x+y;
    cout << z << endl;
    cout << "sum between 2D vectors --------------------------- " << endl;
    cout << "X + Y = ";
    Z = X + Y;
    cout << Z << endl;
    cout << "sum between 1D/2D vectors ------------------------ " << endl;
    cout << "x + Y = ";
    Z = x + Y;
    cout << Z << endl;
    cout << "X + y = ";
    Z = X + y;
    cout << Z << endl;
    cout << "sum between const and 1D vectors ----------------- " << endl;
    cout << "a + x = ";
    z = a + x;
    cout << z << endl;
    cout << "y + b = ";
    z = y + b;
    cout << z << endl;
    cout << "sum between const and 2D vectors ----------------- " << endl;
    cout << "a + X = ";
    Z = a + X;
    cout << Z << endl;
    cout << "Y + b = ";
    Z = Y + b;
    cout << Z << endl;
}

// ================================================================================== //
// OPERATOR "-"                                                                       //
// ================================================================================== //
else if (selection.compare("-") == 0) {

    // Scope variables
    double              a = 1, b = 2, c = 3;
    dvector1D           x, y, z;
    dvector2D           X, Y, Z;

    // Initialize vectors
    x.resize(3, 1.0);
    x[0] = 1.0;    x[1] = 2.0;    x[2] = 3.0;
    y = x;
    X.resize(3, dvector1D(2, 0.0));
    X[0][0] = 1.0;    X[0][1] = 2.0;
    X[1][0] = 3.0;    X[1][1] = 4.0;
    X[2][0] = 5.0;    X[2][1] = 6.0;
    Y = X;

    // Variables
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "x = ";
    cout << x << endl;
    cout << "y = ";
    cout << y << endl;
    cout << "X = ";
    cout << X << endl;
    cout << "Y = ";
    cout << Y << endl;

    // Sum between vectors
    cout << "diff between 1D vectors --------------------------- " << endl;
    cout << "x - y = ";
    z = x - y;
    cout << z << endl;
    cout << "diff between 2D vectors --------------------------- " << endl;
    cout << "X - Y = ";
    Z = X - Y;
    cout << Z << endl;
    cout << "diff between 1D/2D vectors ------------------------ " << endl;
    cout << "x - Y = ";
    Z = x - Y;
    cout << Z << endl;
    cout << "X - y = ";
    Z = X - y;
    cout << Z << endl;
    cout << "diff between const and 1D vectors ----------------- " << endl;
    cout << "a - x = ";
    z = a - x;
    cout << z << endl;
    cout << "y - b = ";
    z = y - b;
    cout << z << endl;
    cout << "diff between const and 2D vectors ----------------- " << endl;
    cout << "a - X = ";
    Z = a - X;
    cout << Z << endl;
    cout << "Y - b = ";
    Z = Y - b;
    cout << Z << endl;
}

// ================================================================================== //
// OPERATOR "*"                                                                       //
// ================================================================================== //
else if (selection.compare("*") == 0) {

    // Scope variables
    double              a = 1, b = 2, c = 3;
    dvector1D           x, y, z;
    dvector2D           X, Y, Z;

    // Initialize vectors
    x.resize(3, 1.0);
    x[0] = 1.0;    x[1] = 2.0;    x[2] = 3.0;
    y = x;
    X.resize(3, dvector1D(2, 0.0));
    X[0][0] = 1.0;    X[0][1] = 2.0;
    X[1][0] = 3.0;    X[1][1] = 4.0;
    X[2][0] = 5.0;    X[2][1] = 6.0;
    Y = X;

    // Variables
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "x = ";
    cout << x << endl;
    cout << "y = ";
    cout << y << endl;
    cout << "X = ";
    cout << X << endl;
    cout << "Y = ";
    cout << Y << endl;

    // Sum between vectors
    cout << "prod between 1D vectors --------------------------- " << endl;
    cout << "x * y = ";
    z = x * y;
    cout << z << endl;
    cout << "prod between 2D vectors --------------------------- " << endl;
    cout << "X * Y = ";
    Z = X * Y;
    cout << Z << endl;
    cout << "prod between 1D/2D vectors ------------------------ " << endl;
    cout << "x * Y = ";
    Z = x * Y;
    cout << Z << endl;
    cout << "X * y = ";
    Z = X * y;
    cout << Z << endl;
    cout << "prod between const and 1D vectors ----------------- " << endl;
    cout << "a * x = ";
    z = a * x;
    cout << z << endl;
    cout << "y * b = ";
    z = y * b;
    cout << z << endl;
    cout << "prod between const and 2D vectors ----------------- " << endl;
    cout << "a * X = ";
    Z = a * X;
    cout << Z << endl;
    cout << "Y * b = ";
    Z = Y * b;
    cout << Z << endl;
}

// ================================================================================== //
// OPERATOR "/"                                                                       //
// ================================================================================== //
else if (selection.compare("/") == 0) {

    // Scope variables
    double              a = 1, b = 2, c = 3;
    dvector1D           x, y, z;
    dvector2D           X, Y, Z;

    // Initialize vectors
    x.resize(3, 1.0);
    x[0] = 1.0;    x[1] = 2.0;    x[2] = 3.0;
    y = x;
    X.resize(3, dvector1D(2, 0.0));
    X[0][0] = 1.0;    X[0][1] = 2.0;
    X[1][0] = 3.0;    X[1][1] = 4.0;
    X[2][0] = 5.0;    X[2][1] = 6.0;
    Y = X;

    // Variables
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "x = " << endl;
    cout << x << endl;
    cout << "y = ";
    cout << y << endl;
    cout << "X = ";
    cout << X << endl;
    cout << "Y = ";
    cout << Y << endl;

    // Sum between vectors
    cout << "div between 1D vectors --------------------------- " << endl;
    cout << "x / y = ";
    z = x / y;
    cout << z << endl;
    cout << "div between 2D vectors --------------------------- " << endl;
    cout << "X / Y = ";
    Z = X / Y;
    cout << Z << endl;
    cout << "div between 1D/2D vectors ------------------------ " << endl;
    cout << "x / Y = ";
    Z = x / Y;
    cout << Z << endl;
    cout << "X / y = ";
    Z = X / y;
    cout << Z << endl;
    cout << "div between const and 1D vectors ----------------- " << endl;
    cout << "a / x = ";
    z = a / x;
    cout << z << endl;
    cout << "y / b = ";
    z = y / b;
    cout << z << endl;
    cout << "div between const and 2D vectors ----------------- " << endl;
    cout << "a / X = ";
    Z = a / X;
    cout << Z << endl;
    cout << "Y / b = ";
    Z = Y / b;
    cout << Z << endl;
}

// ================================================================================== //
// OPERATOR "not available"                                                            //
// ================================================================================== //
else {
    cout << "  Operator not available!" << endl;
}

// ================================================================================== //
// OUTPUT MESSAGE                                                                     //
// ================================================================================== //
{
    cout << "DEMO:  done!!" << endl;
}

return; };

// ---------------------------------------------------------------------------------- //
void vectorMathFunct_Ex(
    void
) {

// ================================================================================== //
// void vectorMathFunct_Ex(                                                           //
//     void)                                                                          //
//                                                                                    //
// Examples of usage of operators for vectors                                         //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - none                                                                             //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - none                                                                             //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
string         selection;

// Counters
// none

// ================================================================================== //
// OUTPUT MESSAGE                                                                     //
// ================================================================================== //
{
    // Select sub case
    cout << "DEMO:  MATH FUNCTIONS FOR STL VECTORS" << endl;
    cout << "  Select operators: " << endl;
    cout << "  'min' operator" << endl;
    cout << "  'max' operator" << endl;
    cout << "  'minval' operator" << endl;
    cout << "  'maxval' operator" << endl;
    cout << "  'sum' operator" << endl;
    cout << "  'abs' operator" << endl;
    cout << "  'pow' operator" << endl;
    cout << "  'norm' operator" << endl;
    cout << "  'Dot_Product' operator" << endl;
    cout << "  'Cross_Product' operator" << endl;
    cin  >> selection;

}

// ================================================================================== //
// OPERATOR "min"                                                                     //
// ================================================================================== //
if (selection.compare("min") == 0) {

    // Scope variables -------------------------------------------------------------- //
    double                       a = 1.5;
    dvector1D                    x(3), t(3);
    dvector2D                    y(3, dvector1D(3)), z(3, dvector1D(3));

    // Initialize scope variables --------------------------------------------------- //
    x[0] = 3.0; x[1] = 3.0; x[2] = 3.0;
    t[0] = 1.0; t[1] = 4.0; t[2] = 2.0;
    y[0][0] = 0.0; y[0][1] = 1.0; y[0][2] = 2.0;
    y[1][0] = 3.0; y[1][1] = 4.0; y[1][2] = 5.0;
    y[2][0] = 6.0; y[2][1] = 7.0; y[2][2] = 8.0;
    z[0][0] = 2.0; z[0][1] = 2.0; z[0][2] = 2.0;
    z[1][0] = 2.0; z[1][1] = 2.0; z[1][2] = 2.0;
    z[2][0] = 2.0; z[2][1] = 2.0; z[2][2] = 2.0;
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "t = " << t << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;

    // min between scalar and 1D array ---------------------------------------------- //
    cout << "min(t, a) = " << min(t,a) << endl;

    // min between 1D arrays -------------------------------------------------------- //
    cout << "min(x, t) = " << min(x, t) << endl;

    // min between scalar and 2D array ---------------------------------------------- //
    cout << "min(y, a) = " << min(y, a) << endl;

    // min between 1D array and 2D array -------------------------------------------- //
    cout << "min(y, x) = " << min(y, x) << endl;

    // min between 2D arrays -------------------------------------------------------- //
    cout << "min(y, z) = " << min(y, z) << endl;

}

// ================================================================================== //
// OPERATOR "minval"                                                                  //
// ================================================================================== //
if (selection.compare("minval") == 0) {

    // Scope variables -------------------------------------------------------------- //
    int                nn = 0;
    double             a = 3.0, mm;
    dvector1D          x(3, 0.0);
    ivector1D          i(3, 1);
    ivector2D          y(3, ivector1D(3, 0));

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 1.1; x[2] = 3.0;
    i[0] = 4;   i[1] = 1;   i[2] = 4;
    y[0] = 2*i;
    y[1] = i;
    y[2] = 2*i;

    // Output message --------------------------------------------------------------- //
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    // Minval of double ------------------------------------------------------------- //
    minval(a,mm);
    cout << "minval(a) = " << mm << endl;

    // Minval of 1D vector ---------------------------------------------------------- //
    minval(x,mm);
    cout << "minval(x) = " << mm << endl;

    // Minval of 2D vector ---------------------------------------------------------- //
    minval(y,nn);
    cout << "minval(y) = " << nn << endl;

}

// ================================================================================== //
// OPERATOR "max"                                                                     //
// ================================================================================== //
else if (selection.compare("max") == 0) {

    // Scope variables -------------------------------------------------------------- //
    double                       a = 1.5;
    dvector1D                    x(3), t(3);
    dvector2D                    y(3, dvector1D(3)), z(3, dvector1D(3));

    // Initialize scope variables --------------------------------------------------- //
    x[0] = 3.0; x[1] = 3.0; x[2] = 3.0;
    t[0] = 1.0; t[1] = 4.0; t[2] = 2.0;
    y[0][0] = 0.0; y[0][1] = 1.0; y[0][2] = 2.0;
    y[1][0] = 3.0; y[1][1] = 4.0; y[1][2] = 5.0;
    y[2][0] = 6.0; y[2][1] = 7.0; y[2][2] = 8.0;
    z[0][0] = 2.0; z[0][1] = 2.0; z[0][2] = 2.0;
    z[1][0] = 2.0; z[1][1] = 2.0; z[1][2] = 2.0;
    z[2][0] = 2.0; z[2][1] = 2.0; z[2][2] = 2.0;
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "t = " << t << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;

    // min between scalar and 1D array ---------------------------------------------- //
    cout << "max(t, a) = " << max(t,a) << endl;

    // min between 1D arrays -------------------------------------------------------- //
    cout << "max(x, t) = " << max(x, t) << endl;

    // min between scalar and 2D array ---------------------------------------------- //
    cout << "max(y, a) = " << max(y, a) << endl;

    // min between 1D array and 2D array -------------------------------------------- //
    cout << "max(y, x) = " << max(y, x) << endl;

    // min between 2D arrays -------------------------------------------------------- //
    cout << "max(y, z) = " << max(y, z) << endl;

}

// ================================================================================== //
// OPERATOR "maxval"                                                                  //
// ================================================================================== //
if (selection.compare("maxval") == 0) {

    // Scope variables -------------------------------------------------------------- //
    int                nn = 0;
    double             a = 3.0, mm;
    dvector1D          x(3, 0.0);
    ivector1D          i(3, 1);
    ivector2D          y(3, ivector1D(3, 0));

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 1.1; x[2] = 3.0;
    i[0] = 1;   i[1] = 2;   i[2] = 4;
    y[0] = i;
    y[1] = 2*i;
    y[2] = 3*i;

    // Output message --------------------------------------------------------------- //
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    // Minval of double ------------------------------------------------------------- //
    maxval(a,mm);
    cout << "maxval(a) = " << mm << endl;

    // Minval of 1D vector ---------------------------------------------------------- //
    maxval(x,mm);
    cout << "maxval(x) = " << mm << endl;

    // Minval of 2D vector ---------------------------------------------------------- //
    maxval(y,nn);
    cout << "maxval(y) = " << nn << endl;

}

// ================================================================================== //
// OPERATOR "sum"                                                                     //
// ================================================================================== //
else if (selection.compare("sum") == 0) {

    // Scope variables -------------------------------------------------------------- //
    int                nn = -5;
    double             a = 3.0, mm;
    dvector1D          x(3, 0.0);
    ivector1D          i(3, 1);
    ivector2D          y(3, ivector1D(3, 0));

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 1.1; x[2] = 3.0;
    i[0] = 1;   i[1] = 2;   i[2] = 4;
    y[0] = i;
    y[1] = 2*i;
    y[2] = 3*i;

    // Output message --------------------------------------------------------------- //
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    // Minval of double ------------------------------------------------------------- //
    sum(a,mm);
    cout << "sum(a) = " << mm << endl;

    // Minval of 1D vector ---------------------------------------------------------- //
    sum(x,mm);
    cout << "sum(x) = " << mm << endl;

    // Minval of 2D vector ---------------------------------------------------------- //
    sum(y,nn);
    cout << "sum(y) = " << nn << endl;

}

// ================================================================================== //
// OPERATOR "abs"                                                                     //
// ================================================================================== //
else if (selection.compare("abs") == 0) {

    // Scope variables -------------------------------------------------------------- //
    int                nn = -5;
    double             a = 3.0, mm;
    dvector1D          x(3, 0.0);
    ivector1D          i(3, 1);
    ivector2D          y(3, ivector1D(3, 0));

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = -1.1; x[2] = 3.0;
    i[0] = 1;   i[1] = 2;   i[2] = 4;
    y[0] = i;
    y[1] = 2*i;
    y[2] = -1*i;

    // Output message --------------------------------------------------------------- //
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    // Minval of double ------------------------------------------------------------- //
    cout << "abs(a) = " << abs(a) << endl;

    // Minval of 1D vector ---------------------------------------------------------- //
    cout << "abs(x) = " << abs(x) << endl;

    // Minval of 2D vector ---------------------------------------------------------- //
    cout << "abs(y) = " << abs(y) << endl;

}

// ================================================================================== //
// OPERATOR "pow"                                                                     //
// ================================================================================== //
else if (selection.compare("pow") == 0) {

    // Scope variables -------------------------------------------------------------- //
    int                nn = -5;
    double             a = 3.0, mm;
    dvector1D          x(3, 0.0);
    ivector1D          i(3, 1);
    ivector2D          y(3, ivector1D(3, 0));

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 1.1; x[2] = 3.0;
    i[0] = 1;   i[1] = 2;   i[2] = 4;
    y[0] = i;
    y[1] = 2*i;
    y[2] = 3*i;

    // Output message --------------------------------------------------------------- //
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    // Minval of double ------------------------------------------------------------- //
    cout << "pow(a,2) = " << pow(a,2.0) << endl;

    // Minval of 1D vector ---------------------------------------------------------- //
    cout << "pow(x,2) = " << pow(x,2.0) << endl;

    // Minval of 2D vector ---------------------------------------------------------- //
    cout << "pow(y,2) = " << pow(y,2.0) << endl;

}

// ================================================================================== //
// OPERATOR "norm"                                                                    //
// ================================================================================== //
else if (selection.compare("norm") == 0) {

    // Scope variables -------------------------------------------------------------- //
    dvector1D          x(3, 0.0);
    ivector1D          i(3, 1);

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 1.1; x[2] = 3.0;
    i[0] = 1;   i[1] = -2;   i[2] = -4;

    // Output message --------------------------------------------------------------- //
    cout << "x = " << x << endl;
    cout << "i = " << i << endl;

    // norm of vector< double > ----------------------------------------------------- //
    cout << "norm(x,2) = " << norm(x,2) << endl;
    cout << "norm(x,3) = " << norm(x,3) << endl;
    cout << "norm_inf(x) = " << norm_inf(x) << endl;

    // norm of vector< int > -------------------------------------------------------- //
    cout << "norm(i,2) = " << norm(i,2) << endl;
    cout << "norm(i,3) = " << norm(i,3) << endl;
    cout << "norm_inf(i) = " << norm_inf(i) << endl;

}

// ================================================================================== //
// OPERATOR "Dot_Product"                                                             //
// ================================================================================== //
else if (selection.compare("Dot_Product") == 0) {

    // Scope variables -------------------------------------------------------------- //
    dvector1D          x(3, 0.0), y(3, 0.0);
    ivector1D          i(3, 1), j(3, 0);

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 1.1; x[2] = 3.0;
    y = x;
    i[0] = 1;   i[1] = -2;   i[2] = 4;
    j = i;

    // Output message --------------------------------------------------------------- //
    cout << "x = " << x << endl;
    cout << "y = " << x << endl;
    cout << "i = " << i << endl;
    cout << "j = " << i << endl;

    // Dot_Product of vector< double > ---------------------------------------------- //
    cout << "Dot_Product(x,y) = " << Dot_Product(x,y) << endl;

    // Dot_Product of vector< int > ------------------------------------------------- //
    cout << "Dot_Product(i,j) = " << Dot_Product(i,j) << endl;

}

// ================================================================================== //
// OPERATOR "Cross_Product"                                                           //
// ================================================================================== //
else if (selection.compare("Cross_Product") == 0) {

    // Scope variables -------------------------------------------------------------- //
    dvector1D          x(3, 0.0), y(3, 0.0);
    ivector1D          i(3, 1), j(3, 0);

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 0.0; x[2] = 0.0;
    y[0] = 0.00; y[1] = 2.1; y[2] = 0.0;
    i[0] = 0;   i[1] = -2;   i[2] = 0;
    j[0] = 4;   j[1] =  0;   j[2] = 0;

    // Output message --------------------------------------------------------------- //
    cout << "x = " << x << endl;
    cout << "y = " << x << endl;
    cout << "i = " << i << endl;
    cout << "j = " << i << endl;

    // Dot_Product of vector< double > ---------------------------------------------- //
    cout << "Cross_Product(x,y) = " << Cross_Product(x,y) << endl;
    cout << "Cross_Product(x,x) = " << Cross_Product(x,x) << endl;

    // Dot_Product of vector< int > ------------------------------------------------- //
    cout << "Cross_Product(i,j) = " << Cross_Product(i,j) << endl;
    cout << "Cross_Product(i,i) = " << Cross_Product(i,i) << endl;

}

// ================================================================================== //
// OUTPUT MESSAGE                                                                     //
// ================================================================================== //
{
    cout << "DEMO:  done!!" << endl;
}

return; };

// ---------------------------------------------------------------------------------- //
void arrayOperators_Ex(
    void
) {

// ================================================================================== //
// void arrayOperators_Ex(                                                            //
//     void)                                                                          //
//                                                                                    //
// Examples of usage of operators for arrays                                          //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - none                                                                             //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - none                                                                             //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
string         selection;

// Counters
// none

// ================================================================================== //
// OUTPUT MESSAGE                                                                     //
// ================================================================================== //
{
    // Select sub case
    cout << "DEMO:  OPERATORS FOR STL ARRAYS" << endl;
    cout << "  Select operators: " << endl;
    cout << "  ""+"" operator" << endl;
    cout << "  ""-"" operator" << endl;
    cout << "  ""*"" operator" << endl;
    cout << "  ""/"" operator" << endl;
    cin  >> selection;
}

// ================================================================================== //
// OPERATOR "+"                                                                       //
// ================================================================================== //
if (selection.compare("+") == 0) {

    // Scope variables
    double                      a = 1, b = 2, c = 3;
    array<double, 2>            x, y, z;
    array<array<double, 2>, 2>  X, Y, Z;

    // Initialize vectors
    x[0] = 1.0;    x[1] = 2.0;
    y = x;
    X[0][0] = 1.0;    X[0][1] = 2.0;
    X[1][0] = 3.0;    X[1][1] = 4.0;
    X[2][0] = 5.0;    X[2][1] = 6.0;
    Y = X;

    // Variables
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "x = ";
    cout << x << endl;
    cout << "y = ";
    cout << y << endl;
    cout << "X = ";
    cout << X << endl;
    cout << "Y = ";
    cout << Y << endl;

    // Sum between vectors
    cout << "sum between 1D arrays ----------------------------- " << endl;
    cout << "x + y = ";
    z = x+y;
    cout << z << endl;
    cout << "sum between 2D arrays ---------------------------- " << endl;
    cout << "X + Y = ";
    Z = X + Y;
    cout << Z << endl;
    cout << "sum between 1D/2D arrays ------------------------- " << endl;
    cout << "x + Y = ";
    Z = x + Y;
    cout << Z << endl;
    cout << "X + y = ";
    Z = X + y;
    cout << Z << endl;
    cout << "sum between const and 1D arrays ------------------ " << endl;
    cout << "a + x = ";
    z = a + x;
    cout << z << endl;
    cout << "y + b = ";
    z = y + b;
    cout << z << endl;
    cout << "sum between const and 2D arrays ------------------ " << endl;
    cout << "a + X = ";
    Z = a + X;
    cout << Z << endl;
    cout << "Y + b = ";
    Z = Y + b;
    cout << Z << endl;
}

// ================================================================================== //
// OPERATOR "-"                                                                       //
// ================================================================================== //
else if (selection.compare("-") == 0) {

    // Scope variables
    double                      a = 1, b = 2, c = 3;
    array<double, 2>            x, y, z;
    array<array<double, 2>, 2>  X, Y, Z;

    // Initialize vectors
    x[0] = 1.0;    x[1] = 2.0;
    y = x;
    X[0][0] = 1.0;    X[0][1] = 2.0;
    X[1][0] = 3.0;    X[1][1] = 4.0;
    X[2][0] = 5.0;    X[2][1] = 6.0;
    Y = X;

    // Variables
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "x = ";
    cout << x << endl;
    cout << "y = ";
    cout << y << endl;
    cout << "X = ";
    cout << X << endl;
    cout << "Y = ";
    cout << Y << endl;

    // Sum between vectors
    cout << "diff between 1D arrays ---------------------------- " << endl;
    cout << "x - y = ";
    z = x - y;
    cout << z << endl;
    cout << "diff between 2D arrays ---------------------------- " << endl;
    cout << "X - Y = ";
    Z = X - Y;
    cout << Z << endl;
    cout << "diff between 1D/2D arrays ------------------------- " << endl;
    cout << "x - Y = ";
    Z = x - Y;
    cout << Z << endl;
    cout << "X - y = ";
    Z = X - y;
    cout << Z << endl;
    cout << "diff between const and 1D arrays ------------------ " << endl;
    cout << "a - x = ";
    z = a - x;
    cout << z << endl;
    cout << "y - b = ";
    z = y - b;
    cout << z << endl;
    cout << "diff between const and 2D arrays ------------------ " << endl;
    cout << "a - X = ";
    Z = a - X;
    cout << Z << endl;
    cout << "Y - b = ";
    Z = Y - b;
    cout << Z << endl;
}

// ================================================================================== //
// OPERATOR "*"                                                                       //
// ================================================================================== //
else if (selection.compare("*") == 0) {

    // Scope variables
    double                      a = 1, b = 2, c = 3;
    array<double, 2>            x, y, z;
    array<array<double, 2>, 2>  X, Y, Z;

    // Initialize vectors
    x[0] = 1.0;    x[1] = 2.0;
    y = x;
    X[0][0] = 1.0;    X[0][1] = 2.0;
    X[1][0] = 3.0;    X[1][1] = 4.0;
    X[2][0] = 5.0;    X[2][1] = 6.0;
    Y = X;

    // Variables
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "x = ";
    cout << x << endl;
    cout << "y = ";
    cout << y << endl;
    cout << "X = ";
    cout << X << endl;
    cout << "Y = ";
    cout << Y << endl;

    // Sum between vectors
    cout << "prod between 1D arrays ---------------------------- " << endl;
    cout << "x * y = ";
    z = x * y;
    cout << z << endl;
    cout << "prod between 2D arrays ---------------------------- " << endl;
    cout << "X * Y = ";
    Z = X * Y;
    cout << Z << endl;
    cout << "prod between 1D/2D arrays ------------------------- " << endl;
    cout << "x * Y = ";
    Z = x * Y;
    cout << Z << endl;
    cout << "X * y = ";
    Z = X * y;
    cout << Z << endl;
    cout << "prod between const and 1D arrays ------------------ " << endl;
    cout << "a * x = ";
    z = a * x;
    cout << z << endl;
    cout << "y * b = ";
    z = y * b;
    cout << z << endl;
    cout << "prod between const and 2D arrays ------------------ " << endl;
    cout << "a * X = ";
    Z = a * X;
    cout << Z << endl;
    cout << "Y * b = ";
    Z = Y * b;
    cout << Z << endl;
}

// ================================================================================== //
// OPERATOR "/"                                                                       //
// ================================================================================== //
else if (selection.compare("/") == 0) {

    // Scope variables
    double                      a = 1, b = 2, c = 3;
    array<double, 2>            x, y, z;
    array<array<double, 2>, 2>  X, Y, Z;

    // Initialize vectors
    x[0] = 1.0;    x[1] = 2.0;
    y = x;
    X[0][0] = 1.0;    X[0][1] = 2.0;
    X[1][0] = 3.0;    X[1][1] = 4.0;
    X[2][0] = 5.0;    X[2][1] = 6.0;
    Y = X;

    // Variables
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "x = " << endl;
    cout << x << endl;
    cout << "y = ";
    cout << y << endl;
    cout << "X = ";
    cout << X << endl;
    cout << "Y = ";
    cout << Y << endl;

    // Sum between vectors
    cout << "div between 1D arrays ---------------------------- " << endl;
    cout << "x / y = ";
    z = x / y;
    cout << z << endl;
    cout << "div between 2D arrays ---------------------------- " << endl;
    cout << "X / Y = ";
    Z = X / Y;
    cout << Z << endl;
    cout << "div between 1D/2D arrays ------------------------- " << endl;
    cout << "x / Y = ";
    Z = x / Y;
    cout << Z << endl;
    cout << "X / y = ";
    Z = X / y;
    cout << Z << endl;
    cout << "div between const and 1D arrays ------------------ " << endl;
    cout << "a / x = ";
    z = a / x;
    cout << z << endl;
    cout << "y / b = ";
    z = y / b;
    cout << z << endl;
    cout << "div between const and 2D arrays ------------------ " << endl;
    cout << "a / X = ";
    Z = a / X;
    cout << Z << endl;
    cout << "Y / b = ";
    Z = Y / b;
    cout << Z << endl;
}

// ================================================================================== //
// OPERATOR "not available"                                                            //
// ================================================================================== //
else {
    cout << "  Operator not available!" << endl;
}

// ================================================================================== //
// OUTPUT MESSAGE                                                                     //
// ================================================================================== //
{
    cout << "DEMO:  done!!" << endl;
}

return; };

// ---------------------------------------------------------------------------------- //
void arrayMathFunct_Ex(
    void
) {

// ================================================================================== //
// void arrayMathFunct_Ex(                                                            //
//     void)                                                                          //
//                                                                                    //
// Examples of usage of operators for C++ v10.0 arrays.                               //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - none                                                                             //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - none                                                                             //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
string         selection;

// Counters
// none

// ================================================================================== //
// OUTPUT MESSAGE                                                                     //
// ================================================================================== //
{
    // Select sub case
    cout << "DEMO:  MATH FUNCTIONS FOR STL ARRAYS" << endl;
    cout << "  Select operators: " << endl;
    cout << "  'min' operator" << endl;
    cout << "  'max' operator" << endl;
    cout << "  'minval' operator" << endl;
    cout << "  'maxval' operator" << endl;
    cout << "  'sum' operator" << endl;
    cout << "  'abs' operator" << endl;
    cout << "  'pow' operator" << endl;
    cout << "  'norm' operator" << endl;
    cout << "  'Dot_Product' operator" << endl;
    cout << "  'Cross_Product' operator" << endl;
    cin  >> selection;

}

// ================================================================================== //
// OPERATOR "min"                                                                     //
// ================================================================================== //
if (selection.compare("min") == 0) {

    // Scope variables -------------------------------------------------------------- //
    double                       a = 1.5;
    array<double, 3>             x, t;
    array<array<double, 3>, 3>   y, z;

    // Initialize scope variables --------------------------------------------------- //
    x[0] = 3.0; x[1] = 3.0; x[2] = 3.0;
    t[0] = 1.0; t[1] = 4.0; t[2] = 2.0;
    y[0][0] = 0.0; y[0][1] = 1.0; y[0][2] = 2.0;
    y[1][0] = 3.0; y[1][1] = 4.0; y[1][2] = 5.0;
    y[2][0] = 6.0; y[2][1] = 7.0; y[2][2] = 8.0;
    z[0][0] = 2.0; z[0][1] = 2.0; z[0][2] = 2.0;
    z[1][0] = 2.0; z[1][1] = 2.0; z[1][2] = 2.0;
    z[2][0] = 2.0; z[2][1] = 2.0; z[2][2] = 2.0;
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "t = " << t << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;

    // min between scalar and 1D array ---------------------------------------------- //
    cout << "min(t, a) = " << min(t,a) << endl;

    // min between 1D arrays -------------------------------------------------------- //
    cout << "min(x, t) = " << min(x, t) << endl;

    // min between scalar and 2D array ---------------------------------------------- //
    cout << "min(y, a) = " << min(y, a) << endl;

    // min between 1D array and 2D array -------------------------------------------- //
    cout << "min(y, x) = " << min(y, x) << endl;

    // min between 2D arrays -------------------------------------------------------- //
    cout << "min(y, z) = " << min(y, z) << endl;

}

// ================================================================================== //
// OPERATOR "minval"                                                                  //
// ================================================================================== //
if (selection.compare("minval") == 0) {

    // Scope variables -------------------------------------------------------------- //
    int                      nn = 0;
    double                   a = 3.0, mm;
    array<double, 3>         x;
    array<int, 3>            i;
    array<array<int, 3>, 3>  y;

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 1.1; x[2] = 3.0;
    i[0] = 4;   i[1] = 1;   i[2] = 4;
    y[0] = 2*i;
    y[1] = i;
    y[2] = -3*i;

    // Output message --------------------------------------------------------------- //
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    // Minval of double ------------------------------------------------------------- //
    minval(a,mm);
    cout << "minval(a) = " << mm << endl;

    // Minval of 1D vector ---------------------------------------------------------- //
    minval(x,mm);
    cout << "minval(x) = " << mm << endl;

    // Minval of 2D vector ---------------------------------------------------------- //
    minval(y,nn);
    cout << "minval(y) = " << nn << endl;

}

// ================================================================================== //
// OPERATOR "max"                                                                     //
// ================================================================================== //
else if (selection.compare("max") == 0) {

    // Scope variables -------------------------------------------------------------- //
    double                       a = 1.5;
    array<double, 3>             x, t;
    array<array<double, 3>, 3>   y, z;

    // Initialize scope variables --------------------------------------------------- //
    x[0] = 3.0; x[1] = 3.0; x[2] = 3.0;
    t[0] = 1.0; t[1] = 4.0; t[2] = 2.0;
    y[0][0] = 0.0; y[0][1] = 1.0; y[0][2] = 2.0;
    y[1][0] = 3.0; y[1][1] = 4.0; y[1][2] = 5.0;
    y[2][0] = 6.0; y[2][1] = 7.0; y[2][2] = 8.0;
    z[0][0] = 2.0; z[0][1] = 2.0; z[0][2] = 2.0;
    z[1][0] = 2.0; z[1][1] = 2.0; z[1][2] = 2.0;
    z[2][0] = 2.0; z[2][1] = 2.0; z[2][2] = 2.0;
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "t = " << t << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;

    // min between scalar and 1D array ---------------------------------------------- //
    cout << "max(t, a) = " << max(t,a) << endl;

    // min between 1D arrays -------------------------------------------------------- //
    cout << "max(x, t) = " << max(x, t) << endl;

    // min between scalar and 2D array ---------------------------------------------- //
    cout << "max(y, a) = " << max(y, a) << endl;

    // min between 1D array and 2D array -------------------------------------------- //
    cout << "max(y, x) = " << max(y, x) << endl;

    // min between 2D arrays -------------------------------------------------------- //
    cout << "max(y, z) = " << max(y, z) << endl;

}

// ================================================================================== //
// OPERATOR "maxval"                                                                  //
// ================================================================================== //
if (selection.compare("maxval") == 0) {

    // Scope variables -------------------------------------------------------------- //
    int                          nn = 0;
    double                       a = 3.0, mm;
    array<double, 3>             x;
    array<int, 3>                i;
    array<array<int, 3>, 3>      y;

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 1.1; x[2] = 3.0;
    i[0] = 1;   i[1] = 2;   i[2] = 4;
    y[0] = i;
    y[1] = 2*i;
    y[2] = -3*i;

    // Output message --------------------------------------------------------------- //
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    // Minval of double ------------------------------------------------------------- //
    maxval(a,mm);
    cout << "maxval(a) = " << mm << endl;

    // Minval of 1D vector ---------------------------------------------------------- //
    maxval(x,mm);
    cout << "maxval(x) = " << mm << endl;

    // Minval of 2D vector ---------------------------------------------------------- //
    maxval(y,nn);
    cout << "maxval(y) = " << nn << endl;

}

// ================================================================================== //
// OPERATOR "sum"                                                                     //
// ================================================================================== //
else if (selection.compare("sum") == 0) {

    // Scope variables -------------------------------------------------------------- //
    int                      nn = -5;
    double                   a = 3.0, mm;
    array<double, 3>         x;
    array<int, 3>            i;
    array<array<int, 3>, 3>  y;

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 1.1; x[2] = 3.0;
    i[0] = 1;   i[1] = 2;   i[2] = 4;
    y[0] = i;
    y[1] = 2*i;
    y[2] = 3*i;

    // Output message --------------------------------------------------------------- //
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    // Minval of double ------------------------------------------------------------- //
    sum(a,mm);
    cout << "sum(a) = " << mm << endl;

    // Minval of 1D vector ---------------------------------------------------------- //
    sum(x,mm);
    cout << "sum(x) = " << mm << endl;

    // Minval of 2D vector ---------------------------------------------------------- //
    sum(y,nn);
    cout << "sum(y) = " << nn << endl;

}

// ================================================================================== //
// OPERATOR "abs"                                                                     //
// ================================================================================== //
else if (selection.compare("abs") == 0) {

    // Scope variables -------------------------------------------------------------- //
    int                      nn = -5;
    double                   a = 3.0, mm;
    array<double, 3>         x;
    array<int, 3>            i;
    array<array<int, 3>, 3>  y;

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = -1.1; x[2] = 3.0;
    i[0] = 1;   i[1] = 2;   i[2] = 4;
    y[0] = i;
    y[1] = 2*i;
    y[2] = -1*i;

    // Output message --------------------------------------------------------------- //
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    // Minval of double ------------------------------------------------------------- //
    cout << "abs(a) = " << abs(a) << endl;

    // Minval of 1D vector ---------------------------------------------------------- //
    cout << "abs(x) = " << abs(x) << endl;

    // Minval of 2D vector ---------------------------------------------------------- //
    cout << "abs(y) = " << abs(y) << endl;

}

// ================================================================================== //
// OPERATOR "pow"                                                                     //
// ================================================================================== //
else if (selection.compare("pow") == 0) {

    // Scope variables -------------------------------------------------------------- //
    int                      nn = -5;
    double                   a = 3.0, mm;
    array<double, 3>         x;
    array<int, 3>            i;
    array<array<int, 3>, 3>  y;

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 1.1; x[2] = 3.0;
    i[0] = 1;   i[1] = 2;   i[2] = 4;
    y[0] = i;
    y[1] = 2*i;
    y[2] = 3*i;

    // Output message --------------------------------------------------------------- //
    cout << "a = " << a << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    // Minval of double ------------------------------------------------------------- //
    cout << "pow(a,2) = " << pow(a,2.0) << endl;

    // Minval of 1D vector ---------------------------------------------------------- //
    cout << "pow(x,2) = " << pow(x,2.0) << endl;

    // Minval of 2D vector ---------------------------------------------------------- //
    cout << "pow(y,2) = " << pow(y,2.0) << endl;

}

// ================================================================================== //
// OPERATOR "norm"                                                                    //
// ================================================================================== //
else if (selection.compare("norm") == 0) {

    // Scope variables -------------------------------------------------------------- //
    array<double, 3>   x;
    array<int, 3>      i;

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 1.1; x[2] = 3.0;
    i[0] = 1;   i[1] = -2;   i[2] = -4;

    // Output message --------------------------------------------------------------- //
    cout << "x = " << x << endl;
    cout << "i = " << i << endl;

    // norm of vector< double > ----------------------------------------------------- //
    cout << "norm(x,2) = " << norm(x,2) << endl;
    cout << "norm(x,3) = " << norm(x,3) << endl;
    cout << "norm_inf(x) = " << norm_inf(x) << endl;

    // norm of vector< int > -------------------------------------------------------- //
    cout << "norm(i,2) = " << norm(i,2) << endl;
    cout << "norm(i,3) = " << norm(i,3) << endl;
    cout << "norm_inf(i) = " << norm_inf(i) << endl;

}

// ================================================================================== //
// OPERATOR "Dot_Product"                                                             //
// ================================================================================== //
else if (selection.compare("Dot_Product") == 0) {

    // Scope variables -------------------------------------------------------------- //
    array<double, 3>   x, y;
    array<int, 3>      i, j;

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 1.1; x[2] = 3.0;
    y = x;
    i[0] = 1;   i[1] = -2;   i[2] = 4;
    j = i;

    // Output message --------------------------------------------------------------- //
    cout << "x = " << x << endl;
    cout << "y = " << x << endl;
    cout << "i = " << i << endl;
    cout << "j = " << i << endl;

    // Dot_Product of vector< double > ---------------------------------------------- //
    cout << "Dot_Product(x,y) = " << Dot_Product(x,y) << endl;

    // Dot_Product of vector< int > ------------------------------------------------- //
    cout << "Dot_Product(i,j) = " << Dot_Product(i,j) << endl;

}

// ================================================================================== //
// OPERATOR "Cross_Product"                                                           //
// ================================================================================== //
else if (selection.compare("Cross_Product") == 0) {

    // Scope variables -------------------------------------------------------------- //
    dvector1D          x(3, 0.0), y(3, 0.0);
    ivector1D          i(3, 1), j(3, 0);

    // Initialize variables --------------------------------------------------------- //
    x[0] = 1.05; x[1] = 0.0; x[2] = 0.0;
    y[0] = 0.00; y[1] = 2.1; y[2] = 0.0;
    i[0] = 0;   i[1] = -2;   i[2] = 0;
    j[0] = 4;   j[1] =  0;   j[2] = 0;

    // Output message --------------------------------------------------------------- //
    cout << "x = " << x << endl;
    cout << "y = " << x << endl;
    cout << "i = " << i << endl;
    cout << "j = " << i << endl;

    // Dot_Product of vector< double > ---------------------------------------------- //
    cout << "Cross_Product(x,y) = " << Cross_Product(x,y) << endl;
    cout << "Cross_Product(x,x) = " << Cross_Product(x,x) << endl;

    // Dot_Product of vector< int > ------------------------------------------------- //
    cout << "Cross_Product(i,j) = " << Cross_Product(i,j) << endl;
    cout << "Cross_Product(i,i) = " << Cross_Product(i,i) << endl;

}

// ================================================================================== //
// OUTPUT MESSAGE                                                                     //
// ================================================================================== //
{
    cout << "DEMO:  done!!" << endl;
}

return; };

// ---------------------------------------------------------------------------------- //
void stringOperators_Ex(
    void
) {

// ================================================================================== //
// void stringOperators_Ex(                                                           //
//     void)                                                                          //
//                                                                                    //
// Examples of usage of operators for strings.                                        //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - none                                                                             //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - none                                                                             //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
string         selection;

// Counters
// none

// ================================================================================== //
// OUTPUT MESSAGE                                                                     //
// ================================================================================== //
{
    // Select sub case
    cout << "DEMO:  OPERATORS FOR STL STRINGS" << endl;
    cout << "  Select operators: " << endl;
    cout << "  'ltrim' operator" << endl;
    cout << "  'rtrim' operator" << endl;
    cout << "  'trim' operator" << endl;
    cin  >> selection;

}

// ================================================================================== //
// OPERATOR "ltrim"                                                                   //
// ================================================================================== //
if (selection.compare("ltrim") == 0) {

    // Scope variables -------------------------------------------------------------- //
    string                       s;

    // Initialize scope variables --------------------------------------------------- //
    s = "  example string  ";
    cout << "s:'" << s << "'" << endl;

    // left trimming ---------------------------------------------------------------- //
    cout << "ltrim(s) = '" << ltrim(s) << "'" << endl;
}

// ================================================================================== //
// OPERATOR "rtrim"                                                                   //
// ================================================================================== //
if (selection.compare("rtrim") == 0) {

    // Scope variables -------------------------------------------------------------- //
    string                       s;

    // Initialize scope variables --------------------------------------------------- //
    s = "  example string  ";
    cout << "s:'" << s << "'" << endl;

    // left trimming ---------------------------------------------------------------- //
    cout << "rtrim(s) = '" << rtrim(s) << "'" << endl;
}

// ================================================================================== //
// OPERATOR "trim"                                                                    //
// ================================================================================== //
if (selection.compare("trim") == 0) {

    // Scope variables -------------------------------------------------------------- //
    string                       s;

    // Initialize scope variables --------------------------------------------------- //
    s = "  example string  ";
    cout << "s:'" << s << "'" << endl;

    // left trimming ---------------------------------------------------------------- //
    cout << "trim(s) = '" << trim(s) << "'" << endl;
}


return; };