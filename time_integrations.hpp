#ifndef TIME_INTEGRATIONS_HPP
#define TIME_INTEGRATIONS_HPP

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <chrono>
#include <random>
#include "euler.hpp"
#include "forest-ruth.hpp"
#include "rk4_integrator.hpp"
#include "verlet.hpp"
#include "Readinput.hpp"
#include "Vec.hpp"
using namespace std;


tuple<vector<Vec>, vector<Vec>, vector<double>> random_initial_conditions(unsigned int& N) {
    /* Function to make random initial conditions for N particles

    Input:  N: Number of particles
    Return: Tuple of three vectors r, v and m, where r and v are vectors of Vecs and m is a vector of doubles
    for all particles and m is a vector of doubles for the masses.
    */

    random_device rd;  // Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<> dis(1, 10); 

    vector<Vec> r;
    vector<Vec> v;
    vector<double> m;

    for (int i=1; i<=N; ++i) { // Generate random position, velocity and mass for the i'th particle
        double xpos = dis(gen); double ypos = dis(gen); double zpos = dis(gen);
        Vec coordinate(xpos, ypos, zpos);

        double xvel = dis(gen); double yvel = dis(gen); double zvel = dis(gen);
        Vec velocity(xvel, yvel, zvel);

        double mass = dis(gen);

        r.push_back(coordinate);
        v.push_back(velocity);
        m.push_back(mass);
    }

    return make_tuple(r, v, m);
}


void time_integration(double h, const double& tmax, const vector<double>& m, vector<Vec>& r, vector<Vec>& v, const string& outputfile,
                      unsigned int integrator(double, const double&, const vector<double>&, vector<Vec>&, vector<Vec>&, const string&, const double&)) {
    /* Measure how long the integration for a given input and integrator takes and append this information to the given outputfile.
     *
     * Input:   h: integration step.
     *          tmax: maximum integration time.
     *          m: vector of doubles with the masses of all the particles.
     *          r: vector of Vecs with the coordinates of all the particles.
     *          v: vector of Vecs with the velocities of all the particles.
     *          outputfile: file to write the output, output is appended to the file.
     *          integrator: pointer to the function that has to be used, with specified argument types.
     *          When calling the time_integration function, a reference to the memory place of the integrator has to be passed with &, for example
     *          time_integration(h, tmax, m, r, v, outputfile, &integrate_rk4)
     *
     * Output:  Output is written to the outputfile as:   N   seconds
     * Return:  The function has no return.
     */

    unsigned int N_part = r.size();
    if ((v.size() != N_part) || (m.size() != N_part)) {cout << "error: r, v and m should have an equal size";}

    ofstream outfile(outputfile, std::ios_base::app | std::ios_base::out);
    outfile << setprecision(10);

    // Get current time.
    auto begintime = chrono::system_clock::now();

    // get initial energy
    double E0 = kineticEnergy(m, v) + potentialEnergy(m, r);

    // Integrate using the desired integrator.
    vector<Vec> r0 = r, v0 = v; // r and v should rest the same for the next integration
    unsigned int function_calls = integrator(h, tmax, m, r0, v0, "none", 0.);

    // get final energy
    double E = kineticEnergy(m, v0) + potentialEnergy(m, r0);

    // Get current time.
    auto endtime = chrono::system_clock::now();

    // Calculate the time needed for the integrator, and dE, and write output
    chrono::duration<double> integrationtime = endtime - begintime;
    outfile << N_part << "\t" << integrationtime.count() << "\t" << function_calls << "\t" << abs(E/E0 - 1) << endl;

    outfile.close();
}


void time_variable_N(const vector<Vec>& r0, const vector<Vec>& v0, const vector<double>& m0, double h,
                     const double& tmax, const string& outputfile, const unsigned int& Nmax,
                     unsigned int integrator(double, const double&, const vector<double>&, vector<Vec>&,
                     vector<Vec>&, const string&, const double&)) {
    /* Time the given integrator for number of particles varying between 1 and Nmax, all with the same given initial conditions.
     *
     * Input:   r0: vector of Vecs with initial conditions of the particles.
     *          v0: vector of Vecs with initial velocities of the particles.
     *          m0: vector with doubles of masses of the particles.
     *          h: integration step.
     *          tmax: maximum integration time.
     *          outputfile: file to write the output, output is appended to the file.
     *          Nmax: maximum number of particles, this has to be greater or equal to the number of particles in r0.
     *          integrator: pointer to the function that has to be used, with specified argument types.
     *          When calling the time_integration function, a reference to the memory place of the integrator has to be passed with &, for example
     *          time_integration(h, tmax, m, r, v, outputfile, &integrate_rk4)
     *
     * Output:  Output is written to the outputfile as:     N   seconds
     *          One has to make sure that the outputfile does not exist yet, as it will not be overwriten, the information will just be added.
     * 
     * Condition: r0, v0 and m0 must contain at least Nmax particles
     * 
     * Return: The function has no return.
     */

     if (r0.size() < Nmax) {cout << "Not enough particles in file" << endl; return;}

     vector<Vec> r;
     vector<Vec> v;
     vector<double> m;

     for (int N=1; N<=Nmax; ++N) {
         cout << "N = " << N << endl;

         r.push_back(r0[N-1]);   // Counting in vector starts at 0.
         v.push_back(v0[N-1]);
         m.push_back(m0[N-1]);

         time_integration(h, tmax, m, r, v, outputfile, integrator);
     }
}

void time_rk4_variable(const vector<Vec>& r0, const vector<Vec>& v0, const vector<double>& m0, double deltamin, double deltamax,
                       const double& tmax, const string& outputfile, const unsigned int& Nmax, const double& hmin=1E-10) {
    /* Time the rk4 with variable timestep for number of particles varying between 1 and Nmax, at random positions.
     *
     * Input:   r0: vector of Vecs with initial conditions of the particles.
     *          v0: vector of Vecs with initial velocities of the particles.
     *          m0: vector with doubles of masses of the particles.
     *          deltamin: parameter for the integration.
     *          deltamax: parameter for the integration.
     *          tmax: maximum integration time.
     *          outputfile: file to write the output, output is appended to the file.
     *          Nmax: maximum number of particles, this has to be greater or equal to the number of particles in r0.
     *          hmin: (optional) the minimum size of the timestep
     *
     * Output:  Output is written to the outputfile as:     N   seconds
     *          One has to make sure that the outputfile does not exist yet, the information will just be added.
     *          This is done using the function time_integration
     *
     * Condition: r0, v0 and m0 must contain at least Nmax particles
     * 
     * Return: The function has no return.
     */

     if (r0.size() < Nmax) {cout << "Not enough particles in file" << endl; return;}

     vector<Vec> r;
     vector<Vec> v;
     vector<double> m;

     ofstream outfile(outputfile);
     outfile << setprecision(10);

     for (int N=1; N<=Nmax; ++N) {
         cout << "N = " << N << endl;

         r.push_back(r0[N-1]);   // Counting in vector starts at 0.
         v.push_back(v0[N-1]);
         m.push_back(m0[N-1]);


         // Get current time.
         auto begintime = chrono::system_clock::now();

         // get initial energy
         double E0 = kineticEnergy(m, v) + potentialEnergy(m, r);

         // Integrate using the desired integrator.
         vector<Vec> r_temp = r, v_temp = v; // r and v should rest the same for the next integration
         unsigned int function_calls = integrate_rk4_variable(deltamin, deltamax, tmax, r_temp, v_temp, m, "none", hmin);

         // get final energy
         double E = kineticEnergy(m, v_temp) + potentialEnergy(m, r_temp);

         // Get current time.
         auto endtime = chrono::system_clock::now();

         // Calculate the time needed for the integrator, and dE, and write output
         chrono::duration<double> integrationtime = endtime - begintime;
         outfile << N << "\t" << integrationtime.count() << "\t" << function_calls << "\t" << abs(E/E0 - 1) << endl;

     }

     outfile.close();
}

#endif
