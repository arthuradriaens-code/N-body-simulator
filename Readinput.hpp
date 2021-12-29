#ifndef READINPUT_HPP
#define READINPUT_HPP

#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <numeric>
#include "Vec.hpp"
using namespace std;

// Compute the potential energy of the system of N bodies with given masses and locations.
// Input:   vector<double> m, the masses of the particles
//          vector<Vec> x,    the positions of the particles at this step
// Condition:   m.size() == x.size()
double potentialEnergy(const vector<double>& m, const vector<Vec>& x) {
    double sum = 0;
    for (int i=0; i<m.size(); ++i) {
        for (int j=0; j<m.size(); ++j) {
            if (i!=j) { 
            sum += m[i]*m[j] / (x[i] - x[j]).norm(); // Second term in equation (4.23) in the lecture notes.
            }
        }
    }
    return -sum/2.; //The prefactor of -1/2 of the second term in equation (4.23) is included here.
}

// Compute the kinetic energy of the system of N bodies with given masses and velocities.
// Input:   vector<double> m, the masses of the particles
//          vector<Vec> v,    the velocities of the particles at this step
// Condition:   m.size() == v.size()
double kineticEnergy(const vector<double>& m, const vector<Vec>& v) {
    double sum = 0;
    for (int i=0; i<m.size(); ++i) {
        sum += m[i] * v[i].norm2(); // First term in equation (4.23) in the lecture notes.
    }
    return sum/2.; //The prefactor of 1/2 of the first term in equation (4.23) is included here.
}

// Compute a location Vec of the centre of mass in the same units as the given location Vec's when given the masses and locations.
// Input:   vector<double> m, the masses of the particles
//          vector<Vec> x,    the positions of the particles at this step
// Condition:   m.size() == x.size()
Vec centreOfMass(const vector<double> m, const vector<Vec>& x) {
    Vec com(0., 0., 0.);
    double M = accumulate(m.begin(), m.end(), 0.); //uses <numeric> to sum over all entries of the vector
    for (int i=0; i<m.size(); ++i) {
        com += m[i] * x[i];
    }
    return com/M;
}

// Compute a velocity Vec of the velocitiy of the centre of mass in the same units as the given velocity Vec's when given the masses and velocities.
// Input:   vector<double> m, the masses of the particles
//          vector<Vec> v,    the velocities of the particles at this step
// Condition:   m.size() == v.size()
Vec systemicVelocity(const vector<double> m, const vector<Vec>& v) {
    Vec imp(0., 0., 0.);
    double M = accumulate(m.begin(), m.end(), 0.); //uses <numeric> to sum over all entries of the vector
    for (int i=0; i<m.size(); ++i) {
        imp += m[i] * v[i]; // There is no need to divide by the total mass, as we normalised the masses previously.
    }
    return imp/M;
}

// Translate the coordinate system by a constant to afix the origin to the centre of mass at initialisation.
// Input:   vector<double> m, the masses of the particles
//          vector<Vec> x,    the positions of the particles at this step
// Condition:   m.size() == x.size()
void centreLocations(const vector<double>& m, vector<Vec>& x) {
    Vec com = centreOfMass(m, x);
    for (int i=0; i<x.size(); ++i) {
        x[i] -= com;
    }
}

// Shift the coordinate system linearly in time such that the total initial momentum is equal to zero.
// Input:   vector<double> m, the masses of the particles
//          vector<Vec> v,    the velocities of the particles at this step
// Condition:   m.size() == v.size()
void centreVelocities(const vector<double>& m, vector<Vec>& v) {
    Vec imp = systemicVelocity(m, v);
    for (int i=0; i<m.size(); ++i) {
        v[i] -= imp;
    }
}

// Translate the coordinate system such that it will always remain centred at the centre of mass
// Input:   vector<double> m, the masses of the particles
//          vector<Vec> x,    the positions of the particles at this step
//          vector<Vec> v,    the velocities of the particles at this step
// Condition:   m.size() == x.size()
// Condition:   m.size() == v.size()
void centralize(vector<double>& m, vector<Vec>& x, vector<Vec>& v) {
    centreLocations(m, x);
    centreVelocities(m, v);
}

/*  Reads the given input file. The first line of the input file is the number of particles N.
    The following N lines each posses the initial conditions for one particle, given as x y z vx vy vz m.
    The function asks the name of the inputfile and three empty vectors which are to be filled with the initial conditions.
    The function alters these three vectors: a vector of doubles with the masses of the particles, a vector of Vec with the 
    coordinates of the particles and a vector of Vec with the velocities.    
*/
void readInput(string filename, vector<Vec>& x, vector<Vec>& v, vector<double>& m){
    ifstream inputfile(filename);

    int N;    // Number of particles.
    // Position components, velocity components and masses of the i-th particle.
    double x_i; double y_i; double z_i;
    double vx_i; double vy_i; double vz_i;
    double m_i;

    if (inputfile.is_open()){
        inputfile >> N; // Read in number of particles

        for(int i=0; i<N; ++i){ // Loop over the following N lines and read in each particle individually
            inputfile >> x_i >> y_i >> z_i >> vx_i >> vy_i >> vz_i >> m_i;

            Vec coord_i(x_i, y_i, z_i);
            Vec veloc_i(vx_i, vy_i, vz_i);
            //Add particle details into the appropriate vectors
            x.push_back(coord_i);
            v.push_back(veloc_i);
            m.push_back(m_i);
        }
    }
    else {cout << "Error: The input file cannot be opened.";}
    centralize(m,x,v); // Shift the positions and velocities of the input such that the system remains centred on the centre of mass.
    inputfile.close();
}

#endif