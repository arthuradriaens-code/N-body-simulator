#ifndef AFS_PROJECT_DE_IFV_EVAL_HPP
#define AFS_PROJECT_DE_IFV_EVAL_HPP

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "euler.hpp"
#include "forest-ruth.hpp"
#include "rk4_integrator.hpp"
#include "rk4_integrator_variable.hpp"
#include "verlet.hpp"
#include "Readinput.hpp"
#include "Vec.hpp"
using namespace std;

void dE_ifv_func_evals_var(const double& deltamin_max, const double& deltamin_min, const double& delta_factor, 
                            const double& tmax, const string& inputfile, const string& outfilename="dE_ifv_func_evals_rk4_var.txt"){
    /*
     * Make a file with the relative energy error in function of the driver function
     * evaluations for rk4 and rk4_var.
     * It performs a few times the same integration,
     * but with smaller and smaller deltamin.
     *
     * inputs: min and max value for deltamin, factor to multiply deltamin each step (<1), tmax,
     *          a file with input conditions, and (optionally) a name for the outputfile.
     *          the default value for the name of the outputfile is "dE_ifv_func_evals_rk4_var.txt"
     * effects: produces a file (name: dE_ifv_eval_var.txt) with three columns:
     *              one column with number of function evaluations,
     *              one column with dE/E0 for the rk4_variable integrator
     *              one with dE/E0 for the normal rk4
     *
     * example:     dE_ifv_func_evals_var(1e-3, 1e-13, 0.1, 20, "Ying-Yang 2a.txt");
     */

    // the integrator changes r and v; let r0 and v0 stay the same initial conditions
    vector<Vec> r, r0, v, v0;
    vector<double> m;
    long double E0;

    //initialize begin conditions
    readInput(inputfile, r0, v0, m);
    int N_part = r0.size();
    E0 = kineticEnergy(m, v0) + potentialEnergy(m, r0);

    // open output file
    ofstream outfile(outfilename);
    outfile << setprecision(30) << endl;

    // perform iteration over h's
    for (double delta=deltamin_max; delta>=deltamin_min; delta *= delta_factor)
    {
        r = r0; v = v0; // reset before each integration
        unsigned int N_evals = integrate_rk4_variable(delta, delta * 50, tmax, r, v, m);
        double E_var_h = kineticEnergy(m, v) + potentialEnergy(m, r); // save final energy

        r = r0; v = v0;
        double h = N_part * (N_part - 1) * tmax * 4 / N_evals; // determine h to get same N_evals with normal rk4 within a rounding error

        integrate_rk4(h, tmax, m, r, v); //deleted ...= as that's not needed
        double E_const_h = kineticEnergy(m, v) + potentialEnergy(m, r); // save final energy


        //calculate delta E for variable and constant time step
        long double dE_var_h = E_var_h/E0 - 1, dE_const_h = E_const_h / E0 - 1;

        outfile << N_evals << '\t' << abs(dE_var_h) << '\t' << abs(dE_const_h) << endl;
        cout << "deltamin = " << delta << "          h: " << h << endl; //for keeping track of the proceeding of the integrations
    }
    outfile.close();
}





#endif //AFS_PROJECT_DE_IFV_EVAL_HPP
