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

// Standard Template Library
# include <iostream>

// CC_Lib
// none

// Others
# include "Examples.hpp"


// ================================================================================== //
// MAIN                                                                               //
// ================================================================================== //
int main(void) {

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //
int      selection;

// ================================================================================== //
// SELECT CASE                                                                        //
// ================================================================================== //
cout << "SELECT DEMO: " << endl;
cout << "0. Basic Operators for STL vectors" << endl;
cout << "1. Math Operators for STL vectors" << endl;
cout << "2. Basic Operators for STL array" << endl;
cout << "3. Math Operators for STL array" << endl;
cout << "4. Operators for STL strings" << endl;
cin >> selection;
switch (selection) {
    case 0: { vectorOperators_Ex(); break; }
    case 1: { vectorMathFunct_Ex(); break; }
    case 2: { arrayOperators_Ex(); break; }
    case 3: { arrayMathFunct_Ex(); break; }
    case 4: { stringOperators_Ex(); break; }
};

return(0); };


