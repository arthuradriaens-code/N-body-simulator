#ifndef AFS_PROJECT_DE_IFV_H_HPP
#define AFS_PROJECT_DE_IFV_H_HPP

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "euler.hpp"
#include "forest-ruth.hpp"
#include "rk4_integrator.hpp"
#include "verlet.hpp"
#include "Readinput.hpp"
#include "Vec.hpp"
using namespace std;


void dE_ifv_h(const double& hmax, const double& hmin, const double& h_factor, const double tmax, const string& infilename, const string& outfilename,
              unsigned int integrator(double, const double&, const vector<double>&, vector<Vec>&, vector<Vec>&, const string&, const double&)){
    /*
     * function that makes a file with three collumns;: 
     *   one with h values,
     *   one with the final relative energy error of each integrator. 
     *   one with the number of driver functions evaluations
     * The same integration is performed several times with smaller and smaller h (divided by the same given factor each time). 
     * The result of the integrations are not written to a file, this ought to run faster than mean_dE_ifv_h(), but may be a less accurate method.
     *
     * inputs: max and min value for h, factor to multiply h each step, the end time of the integration, the name of the input file and
     *          a pointer to the integrator to be examined
     * 
     * effects: produces a file (name: dE_ifv_h + integrator) with three columns; h, dE and the number of driver functions per unit time
     */

    // the integrators change r and v; let r0 and v0 stay the same initial conditions
    vector<Vec> r, r0, v, v0;
    vector<double> m;
    long double E0, E;

    //initialize begin conditions
    readInput(infilename, r0, v0, m);
    E0 = kineticEnergy(m, v0) + potentialEnergy(m, r0);

    // open output file
    ofstream outfile(outfilename);
    outfile << setprecision(10) << endl;

    // perform iteration over h's
    for (double h=hmax; h>=hmin; h *= h_factor)
    {
        r = r0; v = v0; // reset before each integration
        unsigned int driver_functions = integrator(h, tmax, m, r, v, "none", 0.);
        E = kineticEnergy(m, v) + potentialEnergy(m, r);

        //calculate mean delta E
        long double dE = E/E0 - 1;

        outfile << h << '\t' << abs(dE) << '\t' << driver_functions/tmax << endl;
        cout << "h = " << h <<'\n'; //for keeping track of the proceeding of the integrations

    }
    outfile.close();
}


#endif //AFS_PROJECT_DE_IFV_H_HPP
