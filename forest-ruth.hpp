 #ifndef FOREST_RUTH_HPP
#define FOREST_RUTH_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Vec.hpp"
#include "acceleration.hpp"
#include "writeOutput.hpp"
using namespace std;

// The constant defined in equation (5.34) which is used in every step.
double theta = 1./(2.-pow(2., 1./3.));


// Perform all seven steps of the Forest-Ruth algorithm as described in equations (5.33a - 5.33g)
// Input:   double h > 0,     the length of the time step
//          vector<Vec> r,    the positions of all particles in before the updating
//          vector<Vec> v,    the velocity of all particles before the updating
//          vector<double> m, the masses of all particles
// Condition:   m.size() == r.size()
// Condition:   r.size() == v.size()
void updatePositionAndVelocity(const double& h, vector<Vec>& x, vector<Vec>& v, const vector<double>& m) {
    x = x + theta*h/2*v;
    vector<Vec> a = acceleration_all(m, x);
    v = v + theta*h*a;
    x = x + (1-theta)*h/2*v;
    a = acceleration_all(m, x);
    v = v + (1-2*theta)*h*a;
    x = x + (1-theta)*h/2*v;
    a = acceleration_all(m, x);
    v = v + theta*h*a;
    x = x + theta*h/2*v;
}


// Perform the full integration using the Forest-Ruth method.
// Input:   double h > 0,     the length of each time step
//          double tmax,      the time after which the integration is terminated
//          vector<double> m, the masses of the particles
//          vector<Vec> r,    the initial positions of the particles
//          vector<Vec> v,    the initial velocities of the particles
//                            default value is "none"
//          double writeStep, (optional) the integration time that must pass between writing a step of output
//                            default value is 0.
// Condition:   m.size() == r.size()
// Condition:   r.size() == v.size()
// Creates a .txt file with the given filename
// when filename is "none" (default), no file is written
// when writeStep < h (default), every integration step will be writen out if an output file was provided
unsigned int integrate_forest_ruth(double h, const double& tmax, const vector<double>& m, vector<Vec>& r, vector<Vec>& v,
                                     const string& filename="none", const double& writeStep = 0.) {
    double t = 0.;
    unsigned int N_steps = tmax / h; //implicit conversion always rounds off down
    unsigned int N_part = r.size();
    double tWrite = 0.; // Wait until at least this time to write the next line of output

     //initialize output file (has to happen out of if-scope)
     ofstream outfile(filename);

     if (filename != "none") { //when full ouput is needed
         //open output file
         outfile << setprecision(10) << endl;

         //Write the header and the initial conditions into the output file
         writeOutputHeader(outfile, m);
         writeOutputStep(outfile, t, m, r, v);
         tWrite += writeStep; // First step is writen down, so the time for the next line must be updated
     }
    // Integration over the whole time range.
    for (int i=0; i<N_steps; ++i) {
        t += h;
        updatePositionAndVelocity(h, r, v, m);
        
        //output is writen if a filename is provided and at least one writeStep has passed since the previous line of output
        if (filename != "none" && t >= tWrite) {tWrite += writeStep; writeOutputStep(outfile, t, m, r, v);}
    }
    outfile.close();

    return N_part*(N_part - 1)*3*N_steps; // 3 driver function calls per time step
}

#endif
