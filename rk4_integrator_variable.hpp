//
// Created by arnoo on 06-Nov-21.
//

#ifndef RK4_INTEGRATOR_VARIABLE_HPP
#define RK4_INTEGRATOR_VARIABLE_HPP

#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Vec.hpp"
#include "acceleration.hpp"
using namespace std;



tuple<vector<Vec>, vector<Vec>> one_step(const double& h, const vector<Vec>& r, const vector<Vec>& v, const vector<double>& m) {
    /* Function to calculate r and v for one timestep of given length with the rk4 method.
     *
     * Input: Timestep h, vector with initial positions r, vector with initial velocities v, vector with masses.
     * Return: Tuple of new r and v.
     */

    // initialize the k's
    vector<Vec> k1r, k2r, k3r, k4r;
    vector<Vec> k1v, k2v, k3v, k4v;

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
        vector<Vec> new_r = r + (k1r + 2*k2r + 2*k3r + k4r) / 6;
        vector<Vec> new_v = v + (k1v + 2*k2v + 2*k3v + k4v) / 6;

        return make_tuple(new_r, new_v);
}

tuple<vector<Vec>, vector<Vec>, double> determine_h(const double& deltamin, const double& deltamax, const vector<Vec>& r,
                                                    const vector<Vec>& v, const vector<Vec>& r_2h, const vector<Vec>& v_2h,
                                                    const vector<double>& m, double& h, unsigned int &driver_evaluations, const float& hmin) {
    /* Determine h for the step starting at r and v.
     *
     * Input:   deltamin and deltamax are the minimum and maximum values to evaluate the degree of accuracy.
     *          r is the previous vector of position vecs.
     *          v is the previous vector of velocity vecs.
     *          r_2h is the vector of position vecs at t+2h using one timestep of size 2h.
     *          v_2h is the vector of velocity vecs at t+2h using one timestep of size 2h.
     *          m is the vector of all the masses.
     *          h is the timestep. This value is not given as a constant since we want to be able to change h.
     *          driverfunctions is how many interactions between particles have been computed so far.
     *          hmin is the minimum allowed timestep. 
     * Return:  Estimates for r and v at the determined timestep, and how much time passed in this step. The return is formatted as a tuple (r, v, timestep)
     */

    // Initial calculation with h from previous step.

    // Calculate new position and speed at t=2h using two steps of size h.

    tuple<vector<Vec>, vector<Vec>> data_h = one_step(h, r, v, m);
    vector<Vec> r_h = get<0>(data_h);  vector<Vec> v_h = get<1>(data_h);
    tuple<vector<Vec>, vector<Vec>> data_hh = one_step(h, r_h, v_h, m);
    vector<Vec> r_hh = get<0>(data_hh);  vector<Vec> v_hh = get<1>(data_hh);

    unsigned int N_part = r.size();
    driver_evaluations += 2 * 4*N_part*(N_part-1);


    // Delta is the Euclidean norm of a certain particle's position and velocity.
    // We are looking for the largest delta since this is the particle for which the smallest timestep is needed.

    // Calculate the delta for the first particle to compare with later.

    double delta = sqrt((r_hh[0] - r_2h[0]).norm2() + (v_hh[0] - v_2h[0]).norm2());

    // Calculate the delta for all other particles.
    for (int n=1; n<r.size(); ++n){
        double temp = sqrt((r_hh[n] - r_2h[n]).norm2() + (v_hh[n] - v_2h[n]).norm2());
        if (temp > delta){ delta = temp; }  // Keep largest delta.
    }

    // Check if the used step gives the desired degree of accuracy or if it dropped below hmin.
    if ((delta < deltamax && delta > deltamin) || h < hmin) {

        // Return the more precise estimate for t+2h, this is the estimate with the smallest timestep
        // (since this value was calculated anyway). Return the timestep 2h as well.
        return make_tuple(r_hh, v_hh, 2*h);
    }

    else if (delta <= deltamin) {
        h *= 2;

        // Return the more precise estimate for t+2h, this is the estimate with the smallest timestep
        // (since this value was calculated anyway). Return the timestep as well.
        // Note that we already multiplied h by two for the next step, so the actual timestep to return is h.
        return make_tuple(r_hh, v_hh, h);
    }

    else {
        h /= 2;

        // Return an estimate for h/2 (recursion).
        return determine_h(deltamin, deltamax, r, v, r_h, v_h, m, h, driver_evaluations, hmin);
    }

}


unsigned int integrate_rk4_variable(const double& deltamin, const double& deltamax, const double& tmax, vector<Vec>& r,
                                     vector<Vec>& v, const vector<double>& m, const string& filename="none", const double& hmin = 1E-10, const double& writeStep = 0.){
    /* Perform a full integration using the fourth order Runge-Kutta scheme with variable timestep until tmax, using initial conditions.
     * The resulting data is written out to a text file.
     *
     * Inputs: time step h, integration time tmax, vector with masses, vector with initial positions, vector with initial velocities,
     *         an optional string with the name of the resulting file, an optional least time between writing output lines, an optional minimum timestep
     * Effect: performs the integration, r and v change
     *         writes out a line for each time step, in the format:
     *         time    total energy  x(particle 1)   y(particle 1)   z(particle 1) ...
     *         the file begins with a header: one line with the number of particles,
     *         a second line with the masses.
     *         when filename is "none", no file is written
     *         when writeStep > h, the Mth line is only written when t > writeStep*M
     */
    unsigned int driver_evaluations = 0;

    // check whether r0, v0 and m represent the same number of particles
    unsigned int N_part = r.size();
    if ((v.size() != N_part) || (m.size() != N_part)) {cout << "error: r, v and m should have an equal size";}


    // declare quantities that change during the integration
    double t = 0;
    double h = 1;  // Initial value for h, this value will be changed in each integration step.
    double timestep;
    //double hmin = 1e-10;  // Minimal possible value for h
    double tWrite = 0.;  // Wait until at least this time to write the next line of output

    // open output file
    ofstream outfile(filename);

    if (filename != "none") { //when full ouput is needed
         outfile << setprecision(10) << endl;

         //Write the header and the initial conditions into the output file
         writeOutputHeader(outfile, m);
         writeOutputStep(outfile, t, m, r, v);
         tWrite += writeStep;
     }

    // Declare some data types used in the while loop
    tuple<vector<Vec>, vector<Vec>> data_2h;
    vector<Vec> r_2h; vector<Vec> v_2h;
    tuple<vector<Vec>, vector<Vec>, double> next;

    while (t < tmax){
        // Calculate an initial estimate for the position and speed at t=2h using one step of size 2h. This is to make the recursion easier.
        data_2h = one_step(2*h, r, v, m);
        driver_evaluations += 4*N_part*(N_part-1);

        r_2h = get<0>(data_2h); v_2h = get<1>(data_2h);

        next = determine_h(deltamin, deltamax, r, v, r_2h, v_2h, m, h, driver_evaluations, hmin);
        r = get<0>(next);  v = get<1>(next);  timestep = get<2>(next);

        t += timestep;

        if (filename != "none" && t >= tWrite) {tWrite += writeStep; writeOutputStep(outfile, t, m, r, v);}

    }
    outfile.close();

    return driver_evaluations;
}

#endif //RK4_INTEGRATOR_VARIABLE_HPP
