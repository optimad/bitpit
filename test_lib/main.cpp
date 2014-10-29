# include<iostream>
# include <vector>
# include "Class_VolTri.hpp"

using namespace std;

int main(void) {


Class_VolTri        Mesh;

vector< vector< vector<double> > >    x;
x.resize(10);
for (int i = 0; i < x.size(); ++i) {
    x[i].resize(0, vector<double>(2, -1));
}

return 0; } 
