#include <iostream>
#include <vector>
#include <math.h>
#include <assert.h>
using namespace std;

double dot(vector<double>&, vector<double>&);

vector<double> dot(vector<vector<double> >&, vector<double>&);

vector<vector<double> > dot(vector<vector<double> >&, vector<vector<double> >&);

double dot(vector<double>& input1, vector<double>& input2){
    int n1 = input1.size(), n2 = input2.size();
    assert(n1==n2);
    double output = 0;
    for (int ii = 0; ii < n1; ii++){
        output += input1[ii]*input2[ii];
    }    
    return output;
}

vector<double> dot(vector<vector<double> >& input1, vector<double>& input2){
    int n = input1.size(), m = input1[0].size();
    int n2 = input2.size();
    assert(m==n2);
    vector <double> output(n,0);

    for (int ii = 0; ii < n; ii++){
        for (int jj =0; jj < m; jj++){
            output[ii] += input1[ii][jj]*input2[jj];
        }
    }    
    return output;
}

vector<vector<double> > dot(vector<vector<double> >& input1, vector<vector<double> >& input2){
    int n1 = input1.size(), m1 = input1[0].size();
    int n2 = input2.size(), m2 = input2[0].size();
    assert(m1==n2) ;
    
    vector<vector <double> > output;
    output.resize(n1, vector<double>(m2, 0.0));

    for (int ii = 0; ii < n1; ii++){
        for (int jj = 0; jj < m2; jj++){
            for (int kk = 0; kk < m1; kk++){
                output[ii][jj] += input1[ii][kk]*input2[kk][jj];
            }
        }
    }    
    return output;
}

// int main(){
//     vector<double> a1 = {1, 2, 3};
//     vector<double> b1 = {1, 2, 3};
//     double c1 = dot(a1,b1);
//     cout << c1 << endl; // verified

//     vector<vector<double> > a2 = {{1, 2, 3}, {2,3,4}, {4,5,6}};
//     vector<double> b2 = {1, 2, 3};
//     vector<double> c2 = dot(a2,b2);
//     for (int ii = 0; ii < 3; ii++) cout << c2[ii] << endl; // verified

//     vector<vector<double> > a3 = {{1, 2, 3}, {2,3,4}, {4,5,6}};
//     vector<vector<double> > b3 = {{1, 2, 3}, {2,3,4}, {4,5,6}};
//     vector<vector<double> > c3 = dot(a3,b3);
//     for (int ii = 0; ii < 3; ii++) {
//         cout << "\n";
//         for (int jj = 0; jj < 3; jj++) {
//             cout << c3[ii][jj] << " "; 
//         }
//         cout << "\n";
//     }

//     return 0; 
// }