#include <iostream>
#include "examples.hpp"

using namespace std;

int main(){

int selection ;

// Output message ------------------------------------------------------- //
// cout << "===================== Input/ Output DEMO ===================== " << endl;
// cout << "'1': Basic VTK demo" << endl;
// cin >> selection;

selection =1 ;


// Run demo ------------------------------------------------------------- //
switch (selection) {
    case 1: { Demo_BasicVTK(); break; }
    default: break;
}

  return 0;
};
