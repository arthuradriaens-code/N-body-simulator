 #ifndef PEFRL_HPP
#define PEFRL_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Vec.hpp"
#include "acceleration.hpp"
#include "writeOutput.hpp"

// Some constants required in the integration 
double Xi = 0.1786178958448091;
double lambda = -0.2123418310626054;
double chi = -0.06626458266981849;

// Perform all nine steps of the PEFRL algorithm as described in equations (5.35a - 5.35i)
// Input:   double h > 0,     the length of the time step
//          vector<Vec> r,    the positions of all particles before the updating
//          vector<Vec> v,    the velocity of all particles before the updating
//          vector<double> m, the masses of all particles
// Condition:   m.size() == r.size()
// Condition:   r.size() == v.size()
void updatePositionAndVelocity_PEFRL(const double& h, const vector<double>& m, vector<Vec>& r, vector<Vec>& v) {
    // Step 1 (5.35a)
    r = r + Xi*h*v;
    // Step 2 (5.35b)
    vector<Vec> a = acceleration_all(m, r);
    v = v + (1-2*lambda)/2*h*a;
    // Step 3 (5.35c)
    r = r + chi*h*v;
    // Step 4 (5.35d)
    a = acceleration_all(m, r);
    v = v + lambda*h*a;
    // Step 5 (5.35e)
    r = r + (1 - 2*(chi+Xi))*h*v;
    // Step 6 (5.35f)
    a = acceleration_all(m, r);
    v = v + lambda*h*a;
    // Step 7 (5.35g)
    r = r + chi*h*v;
    // Step 8 (5.35h)
    a = acceleration_all(m, r);
    v = v + (1-2*lambda)/2*h*a;
    // Step 9 (5.35i)
    r = r + Xi*h*v;
}

// Perform the full integration using the Position Extended Forest-Ruth Like (PEFRL) method.
// Input:   double h > 0,     the length of each time step
//          double tmax,      the time after which the integration is terminated
//          vector<double> m, the masses of the particles
//          vector<Vec> r,    the initial positions of the particles
//          vector<Vec> v,    the initial velocities of the particles
//          string filename,  (optional) the name to be given to the file in which the integration will be recorded
//                            default value is "none"
//          double writeStep, (optional) the integration time that must pass between writing a step of output
//                            default value is 0.
// Condition:   m.size() == r.size()
// Condition:   r.size() == v.size()
// Creates a .txt file with the given filename
// when filename is "none" (default), no file is written
// when writeStep < h (default), every integration step will be writen out if an output file was provided
unsigned int integrate_pefrl(double h, const double& tmax, const vector<double>& m, vector<Vec>& r, vector<Vec>& v, 
                                const string& filename="none", const double& writeStep = 0.) {

    unsigned int N_steps = tmax / h; // implicit conversion always rounds off down
    double t = 0.;
    unsigned int N_part = r.size();
    double tWrite = 0.; // Wait until at least this time to write the next line of output

    //open output file
    ofstream outfile(filename);

    if (filename != "none") {
        outfile << setprecision(10) << endl;

        //Write the header and the initial conditions into the output file
        writeOutputHeader(outfile, m);
        writeOutputStep(outfile, t, m, r, v);
        tWrite += writeStep; // First step is writen down, so the time for the next line must be updated
    }

    // Integration over the whole time range.
    for (int i=0; i<N_steps; ++i){
        t += h;
        updatePositionAndVelocity_PEFRL(h, m, r, v);

        //output is writen if a filename is provided and at least one writeStep has passed since the previous line of output
        if (filename != "none" && t >= tWrite) {tWrite += writeStep; writeOutputStep(outfile, t, m, r, v);}
    }
    outfile.close();

    return N_part*(N_part - 1)*4*N_steps;  // 4 driver function calls per time step
}



#endif
