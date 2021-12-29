#define _USE_MATH_DEFINES
#include <vector>
#include <numeric>
#include <iostream>
#include <cmath>
#include "Vec.hpp"          //watch out order matters, hpp's have to be included last
#include "Readinput.hpp"    //This contains Ekin and Epot so this comes first
#include "euler.hpp"
#include "rk4_integrator.hpp"
#include "rk4_integrator_variable.hpp"
#include "verlet.hpp"
#include "time_integrations.hpp"
#include "dE_ifv_h.hpp"
#include "dE_ifv_eval.hpp"
#include "pefrl.hpp"
using namespace std;



int main(){
    string filename;
    double h;
    double tmax;
    double deltamin; double deltamax;
    int IntegrationChoice;
    int GlobalChoice;
    vector<Vec> x;
    vector<Vec> v;
    vector<double> m;

    cout << "Enter the name of the input file (e.g. 'input'): ";
    cin >> filename;
    filename = filename + ".txt";
    readInput(filename, x, v, m);

    //print how the input looks after centralizing (this is done in readInput).
    cout << "How the initial positions look after centralizing:" << endl;
    cout << "masses:" << endl;
    for (int i = 0; i < m.size(); ++i) { cout << m[i] << endl; }
    cout << "positions:" << endl;
    for (int i = 0; i < x.size(); ++i) { print(x[i]); }
    cout << "velocities:" << endl;
    for (int i = 0; i < v.size(); ++i) { print(v[i]); }
    cout << endl;
    // Printing Epot and Ekin.
    cout << "Epot and Ekin:" << endl;
    cout << potentialEnergy(m, x) << endl;
    cout << kineticEnergy(m, v) << endl;
    cout << endl;
    //Global menu:
    cout << "What would you like to do? The options are:" << endl;
    cout << "1. integrate" << endl;
    cout << "2. calculate energy errors in function of h \n   and in function of the number of driver function evaluations for rk4 variable" << endl;
    cout << "3. time the integrations as a function of a number of random placed particles" << endl;
    cout << "\ninsert a number from 1 to 3: ";
    cin >> GlobalChoice;
    
    if (GlobalChoice==1){
        //menu for choosing integration scheme:
        cout << "Which integrator would you like to use?\nthe options are:" << endl;
        cout << "1. euler" << endl;
        cout << "2. verlet" << endl;
        cout << "3. forest-ruth" << endl;
        cout << "4. pefrl" << endl; //what integrater is this?
        cout << "5. rk4" << endl;
        cout << "6. rk4 with variable timestep" << endl;
        cout << "Insert a number from 1 to 6: ";
        cin >> IntegrationChoice;

        if (0 < IntegrationChoice && IntegrationChoice < 6) {
            cout << "Integration step size (e.g 0.001): "; // Type a number and press enter
            cin >> h; // Get user input from the keyboard
            cout << "Integration time (e.g 10): "; // Type a number and press enter
            cin >> tmax; // Get user input from the keyboard

            cout << "Driver function evaluations per unit time: ";
            if (IntegrationChoice == 1){
                cout << integrate_euler(h,tmax,m,x,v,"euler.txt")/tmax << endl;
            }
            if (IntegrationChoice == 2){
                cout << integrate_verlet(h,tmax,m,x,v,"verlet.txt")/tmax << endl;
            }
            if (IntegrationChoice == 3){
                cout << integrate_forest_ruth(h,tmax,m,x,v,"forest-ruth.txt")/tmax << endl;
            }
            if (IntegrationChoice == 4){
                cout << integrate_pefrl(h,tmax,m,x,v,"pefrl.txt")/tmax << endl;
            }
            if (IntegrationChoice == 5){
                cout << integrate_rk4(h,tmax,m,x,v,"rk4.txt")/tmax << endl;
            }
        }
        else if (IntegrationChoice == 6) {
            cout << "Deltamin used for rk4 with variable timestep (e.g 0.00000001): ";
            cin >> deltamin;
            cout << "Deltamax used for rk4 with variable timestep (e.g 0.001): ";
            cin >> deltamax;
            cout << "Integration time (e.g 10): "; // Type a number and press enter
            cin >> tmax; // Get user input from the keyboard

            cout << "Driver function evaluations per unit time: ";
            cout << integrate_rk4_variable(deltamin,deltamax,tmax,x,v,m,"rk4-variable.txt")/tmax << endl;
        }
        else {
            cout << "That is not a valid number" << endl;
        }
    }
    
    else if (GlobalChoice==2) {
        double hmax;
        double hmin;
        double hfactor;
        double deltamin_min, deltamin_max;

        cout << "Maximum value for h: ";
        cin >> hmax;
        cout << "Minimal value for h: ";
        cin >> hmin;
        cout << "Maximal value for delta_min (variable timestep): ";
        cin >> deltamin_max;
        cout << "Minimal value for delta_min (variable timestep): ";
        cin >> deltamin_min;
        cout << "Factor to change h and deltamin (this number has to be smaller than one): ";
        cin >> hfactor;
        cout << "Integration time: ";
        cin >> tmax;

        cout << "Calculating dE for variable h for rk4..." << endl;
        dE_ifv_h(hmax, hmin, hfactor, tmax, filename, "dE_ifv_h_rk4.txt", &integrate_rk4);
        cout << "Calculating dE for variable hfor verlet..." << endl;
        dE_ifv_h(hmax, hmin, hfactor, tmax, filename, "dE_ifv_h_verlet.txt", &integrate_verlet);
        cout << "Calculating dE for variable h for euler..." << endl;
        dE_ifv_h(hmax, hmin, hfactor, tmax, filename, "dE_ifv_h_euler.txt", &integrate_euler);
        cout << "Calculating dE for variable h for forest_ruth..." << endl;
        dE_ifv_h(hmax, hmin, hfactor, tmax, filename, "dE_ifv_h_forest_ruth.txt", &integrate_forest_ruth);
        cout << "Calculating dE for variable h for pefrl..." << endl;
        dE_ifv_h(hmax, hmin, hfactor, tmax, filename, "dE_ifv_h_pefrl.txt", &integrate_pefrl);
        cout << "Calculating dE and driver function evaluations in function of delta_min for rk4 with and without variable timestep" << endl;
        dE_ifv_func_evals_var(deltamin_max, deltamin_min, hfactor, tmax, filename, "dE_ifv_func_evals_rk4_var.txt");

        // Note: the last function gives another output file format than the others

    }
    else if (GlobalChoice==3) {
        unsigned int Nmax;
        cout << "Maximum number of particles: ";
        cin >> Nmax;
        cout << "Integration step for the integrators with fixed timestep: ";
        cin >> h;
        cout << "Deltamax for rk4 with variable timestep: ";
        cin >> deltamax;
        cout << "Deltamin for rk4 with variable timestep: ";
        cin >> deltamin;
        cout << "Integration time: ";
        cin >> tmax;


        tuple<vector<Vec>, vector<Vec>, vector<double>> initial_conditions = random_initial_conditions(Nmax);
        vector<Vec> r0 = get<0>(initial_conditions);
        vector<Vec> v0 = get<1>(initial_conditions);
        vector<double> m0 = get<2>(initial_conditions);

        // Write out random coordinates
        ofstream outfile("randomcoordinates.txt");
        outfile << setprecision(10); // Precision of the output. Increase for less round-off in output, decrease to decrease size of output file.
        outfile << "# tmax = " << tmax << "\t h = " << h << "\t deltamin = " << deltamin << "\t deltamax = " << deltamax << endl;
        for (int i = 0; i < r0.size(); ++i){
            outfile << r0[i].x() << "\t" << r0[i].y() << "\t" << r0[i].z() << "\t" << v0[i].x() << "\t" << v0[i].y() << "\t" << v0[i].z() << endl;
        }
        outfile.close();

        cout << "Timing rk4 for variable N... " << endl;
        time_variable_N(r0, v0, m0, h, tmax, "time_rk4.txt", Nmax, &integrate_rk4);
        cout << "Timing verlet for variable N... " << endl;
        time_variable_N(r0, v0, m0, h, tmax, "time_verlet.txt", Nmax, &integrate_verlet);
        cout << "Timing euler for variable N... " << endl;
        time_variable_N(r0, v0, m0, h, tmax, "time_euler.txt", Nmax, &integrate_euler);
        cout << "Timing forest ruth for variable N... " << endl;
        time_variable_N(r0, v0, m0, h, tmax, "time_forest_ruth.txt", Nmax, &integrate_forest_ruth);
        cout << "Timing pefrl for variable N... " << endl;
        time_variable_N(r0, v0, m0, h, tmax, "time_pefrl.txt", Nmax, &integrate_pefrl);
        cout << "Timing rk4_variable for variable N... " << endl;
        time_rk4_variable(r0, v0, m0, deltamin, deltamax, tmax, "time_rk4_variable.txt", Nmax);
    }
    else {cout << "that's not a valid number" << endl;}
}
