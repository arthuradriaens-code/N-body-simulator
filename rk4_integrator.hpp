//
// Created by arnoo on 06-Nov-21.
//

#ifndef RK4_INTEGRATOR_HPP
#define RK4_INTEGRATOR_HPP


#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Vec.hpp"
#include "writeOutput.hpp"
#include "acceleration.hpp"
using namespace std;



unsigned int integrate_rk4(double h, const double& tmax, const vector<double>& m, vector<Vec>& r, vector<Vec>& v, const string& filename="none", const double& writeStep = 0.){
    /* Perform a full integration using the fourth order Runge-Kutta scheme with time step h until tmax, using initial conditions.
     * The resulting data is written out to a text file.
     *
     * Inputs: time step h, integration time tmax, vector with masses, vector with initial positions, vector with initial velocities,
     *         an optional string with the name of the resulting file, an optional least time between writing output lines
     * Effect: performs the integration, r and v change
     *         writes out a line for each time step, in the format:
     *         time    total energy  x(particle 1)   y(particle 1)   z(particle 1) ...
     *         the file begins with a header: one line with the number of particles,
     *         a second line with the masses.
     *         when filename is "none", no file is written
     *         when writeStep > h, the Mth line is only written when t > writeStep*M
     */
    
    //check whether r, v and m represent the same number of particles    
    unsigned int N_part = r.size();
    if ((v.size() != N_part) || (m.size() != N_part)) {cout << "error: r, v and m should have an equal size";}

    unsigned int N_steps = tmax / h; //implicit conversion always rounds off down
    double t = 0;
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

    // initialize the k's
    vector<Vec> k1r, k2r, k3r, k4r;
    vector<Vec> k1v, k2v, k3v, k4v;

    //perform each time step
    for (int i=0; i<N_steps; ++i){

        // calculate the k's according to equation 4.35-4.38
        k1v = h * acceleration_all(m, r);
        k1r = h * v;

        k2v = h * acceleration_all(m, r + 0.5 * k1r);
        k2r = h * (v + 0.5 * k1v);

        k3v = h * acceleration_all(m, r + 0.5 * k2r);
        k3r = h * (v + 0.5 * k2v);

        k4v = h * acceleration_all(m, r + k3r);
        k4r = h * (v + k3v);

        // calculate the next positions and velocities, according to equation 4.39
        r = r + (k1r + 2*k2r + 2*k3r + k4r) / 6;
        v = v + (k1v + 2*k2v + 2*k3v + k4v) / 6;

        t += h;

         //Write the header and the initial conditions into the output file
        if (filename != "none" && t >= tWrite) {tWrite += writeStep; writeOutputStep(outfile, t, m, r, v);}
    }
    outfile.close();

    return N_part *(N_part - 1)*4*N_steps;// 4 driver function calls per time step
}


#endif //RK4_INTEGRATOR_HPP
