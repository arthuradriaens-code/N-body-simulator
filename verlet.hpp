#ifndef VERLET_HPP
#define VERLET_HPP


#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Vec.hpp"
#include "writeOutput.hpp"
#include "acceleration.hpp"
using namespace std;

unsigned int integrate_verlet(double h, const double& tmax, const vector<double>& m, vector<Vec>& r, vector<Vec>& v, const string& filename="none", const double& writeStep = 0.){
    /* Perform a full integration using the position Verlet scheme with time step h until tmax, using initial conditions.
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
    //declare quantities that change during the integration
    vector<Vec>r_n12,a_n12;

    //perform each time step
    for (int i=0; i<N_steps; ++i){

        r_n12 = r + 0.5*h*v;
        a_n12 = acceleration_all(m, r_n12);
        v = v + h*a_n12;
        r = r_n12 + 0.5*h*v;

        t += h;

         //Write the header and the initial conditions into the output file
        if (filename != "none" && t >= tWrite) {tWrite += writeStep; writeOutputStep(outfile, t, m, r, v);}

    }
    outfile.close();

    return N_part*(N_part - 1)*N_steps; // 1 driver function call per time step, for N particles N times
}

#endif
