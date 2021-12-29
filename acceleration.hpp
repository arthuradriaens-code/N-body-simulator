 #ifndef ACCELERATION_HPP
#define ACCELERATION_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Vec.hpp"
using namespace std;

vector<Vec> acceleration_all(const vector<double>& m, const vector<Vec>& r){
    /* Compute the acceleration of all particles
     *
     * Inputs: a vector with masses and a vector with positions (Vec) of the same length
     * Output: a vector with the acceleration (Vec) for each particle
     * standard units are assumed (G=1)
     */
    vector<Vec> a_all; // this will be the final output vector

    for (int i=0; i<r.size(); ++i) {
        Vec a(0, 0, 0); //the acceleration for particle i
        for (int j=0; j<r.size(); ++j) {
            if (i!=j) {
                a -=  m[j] * (r[i] - r[j]) / (r[i] - r[j]).norm3(); // equation (4.5) with G=1
            }
        }
        a_all.push_back(a);
    }
    return a_all;
}


#endif
