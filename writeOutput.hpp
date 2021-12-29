 #ifndef WRITEOUTPUT_HPP
#define WRITEOUTPUT_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Vec.hpp"

// Writes a short header for the given outfile containing the number of particles in its first line and the masses 
// of the particles in the second line, seperated by a tab.
// Input:   ofstream outfile, a previously output stream in which the header will be written. Precision must already be defined
//          vector<double> m, the masses of the particles
void writeOutputHeader(ofstream& outfile, const vector<double>& m) {
    int N = m.size();
    outfile << '#' << N << endl << '#';
    for (int i=0; i<N; ++i){outfile << m[i] << '\t';}
    outfile << endl;
}

// Writes the results of a single step of the integration into the given outfile in a single line. The first element will be the
// time of the step, the second the total energy of the system and the subsequent entrees the x, y and z coordinates of each particle.
// Input:   ofstream outfile, a previously output stream in which the data will be written. Precision must already be defined
//          double t, the time of the step
//          vector<double> m, the masses of the particles
//          vector<Vec> r, the positions of the particles at this step
//          vector<Vec> v, the velocties of the particles at this step
// Condition:   m.size() == r.size()
// Condition:   m.size() == v.size()
void writeOutputStep(ofstream& outfile, const double& t, const vector<double>& m, const vector<Vec>& r, const vector<Vec>& v) {
    int N = m.size();
    double E = kineticEnergy(m, v) + potentialEnergy(m, r);
    outfile << t << '\t' << E;
    for (int i=0; i<N; ++i) {outfile << '\t' << r[i].x() << '\t' << r[i].y() << '\t' << r[i].z();}
    outfile << endl;
}

#endif