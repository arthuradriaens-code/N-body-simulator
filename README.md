This is a n-body simulation program which uses G=1, created by Arno Depoorter, Mathijs Vanrespaille, Silke Van der Schueren and Arthur Adriaens. There's also a plot.py file which can be used to plot the result but is not required.

How to use:

1. Setting initial conditions:
The initial conditions should be located in  a text file (e.g input.txt) and should follow the structure presented below.


number of particles
x | y | z | vx | vy | vz | mass
-- | -- | -- | -- | -- | -- | --
x_1 | y_1 | z_1 | vx_1 | vy_1 | vz_1 | mass_1
x_2 | y_2 | z_2 | vx_2 | vy_2 | vz_2 | mass_2
x_3 | y_3 | z_3 | vx_3 | vy_3 | vz_3 | mass_3
.
.
.
x_N | y_N | z_N | vx_N | vy_N | vz_N | mass_N
---

with x_1 denoting the x coördinate of particle 1 and vy_2 denoting the velocity in the y direction of particle 2 with tabs separating the values. An example of the style of this file can be found in input.txt.
Note that the integration works in its own integration unit system where the gravitational constant G=1. To convert physical units to these units, just solve for time in the units.

Some nice systems can be found [online](http://three-body.ipb.ac.rs/) for example:

https://user-images.githubusercontent.com/66306556/147673874-b1f1dae2-aca2-42a2-bb3a-5ba904f7ab68.mp4



2. Compiling:
you should compile main.cpp to use this program, a suggested way of doing this is via navigating to this folder in bash/zsh/cmd/... and typing
"g++ -std=c++14 -Wall -O3 -o main main.cpp" (assuming you have [gnu's c compiler](https://gcc.gnu.org/) installed)
After this an execution file should be created (.sh in linux/mac and .exe in windows). Which you can then drag and drop into the terminal and hit enter to run.
Alternatively, the files can also be imported into a code editor such as Visual Studio Code and run there. 

3. Using the program
Upon starting the program you're prompted to enter the text file name in which you specified the initial conditions, here you just type in the filename without the extension (e.g input). After hitting enter you're asked to select what you want to do, the possible options are:

1. integrate, 
2. calculate energy errors in function of h and in function of the number of driver function evaluations for rk4 with variable timestep
3. time the integrations as a function of a number of random placed particles

You choose by entering a number (1 to 3) and pressing enter. If integration is selected another prompt will show up to ask what integration method. If rk4-variable is selected a minimum and maximum delta is asked, this determines the number of timestep doubling that can happen if particles get close. 
If any other integration scheme is selected you'll be prompted to enter the timestep manually. 
The output of this program is of the form "integrationmethod.txt" e.g "verlet.txt" with the following structure:

(EmptyLine)

#number of particles

#mass_1 mass_2 mass_3 ... mass_N

time | totalEnergy | x_1 | y_1 | z_1 | x_2 | y_2 | z_2 | x_3 | y_3 | z_3 | ... | x_N | y_N | z_N
-- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |
time1 | totalEnergy1 | x_1 | y_1 | z_1 | x_2 | y_2 | z_2 | x_3 | y_3 | z_3 | ... | x_N | y_N | z_N
time2 | totalEnergy2 | x_1 | y_1 | z_1 | x_2 | y_2 | z_2 | x_3 | y_3 | z_3 | ... | x_N | y_N | z_N
time3 | totalEnergy3 | x_1 | y_1 | z_1 | x_2 | y_2 | z_2 | x_3 | y_3 | z_3 | ... | x_N | y_N | z_N
.
.
.
tmax | totalEnergyFinal | x_1 | y_1 | z_1 | x_2 | y_2 | z_2 | x_3 | y_3 | z_3 | ... | x_N | y_N | z_N

Where time and totalEnergy are updated at each step, but not indexed to avoid confusion with the indices of the particles. The rest of the variables are as in input.txt. The (EmptyLine) i.e a '\n' is an artefact of setprecision() and should be ignored.

The main can be edited, for which the user is directed to the comments of each function and hpp file. Notably, there are some optional parameters in the six integrator functions.
Firstly, all filenames are optional in the integrate_ function, with the default value being "none". If this value is used, no output will be written to an outfile. 
Secondly, in rk4_variable, there is an optional parameter hmin which determines the smallest allowed timestep, to avoid the integration freezing during a direct collision where h would continually decrease. 
Finally, all integrate_ functions have an optional parameter writeStep that can be used to limit the number of output lines and thus the size of the output files. It represents the minimum time that must pass between the program writing two lines of output. If it is set to zero, its default value, every step will be writen out. If, for instance, writeStep = 1 and h=0.15, output will be written at t = {0, 1.05, 2.1, 3, 4.05 ...}
