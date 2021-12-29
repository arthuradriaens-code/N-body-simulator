from mpl_toolkits import mplot3d
import sys
import numpy as np
import matplotlib.pyplot as plt
#from mayavi import mlab

from matplotlib import animation
from IPython.display import HTML


# ### Some relevant parameters

# Select which plots and animations to compute and save
make3DtrajectoryPlot = False
makeSingleEnergyPlots = False
makeCombinedEnergyPlot = False
makeCOMdeviationPlot= False
make3Danimations = True
make_E_ifv_h_plot = False
make_E_ifv_driverfunctions = False
makeSingle2DtrajectoryPlots = False    # Plot of the xy-plane
makeCombined2DtrajectoryPlot = False   # Plot of the xy-plane

equaliseAspect = False # Set to True if you want all axes to have the same scale
# Do be careful with the 3D animation though, as it changes the camera angle
# from what you input below.

# The variables below are only relevant if make3Danimations==True

# How many steps are condensed into a single point in the 3D animation.
# Decreasing it results in a better animation resolution but a slower computation
stepsPerPoint = 2
# IF you want the animation to end earlier than the largest time in the file you
# read in, set this to when you want the animation to end.
endTime = 50 # sys.float_info.max()
fps = 10
# The elevation and azimuth of the camera in the 3D plots, in degrees.
elevation = 50
azimuth = 30
# The relative margin to be used in the limits of the axes of the 3D animation
# The formula for e.g. the lower x boundary will be
#           x_lower = x_min - relativeMargin*(x_max-x_min) + additiveMargin
# The additive margin is needed in cases where, for instance, z is always zero.
relativeMargin = 0.1
additiveMargin = 0.05



# ### Preparing to read in the data

# A list that will contain strings of all the methods selected.
method_all = []

# Lists that will contain a numpy array for each method in the order of selection.
t_all = []
E_all = []
Ereldif_all = []

# Lists that will contain a 2D array for each particle's coordinate per method.
# Yes, I know it's not pretty.
x_all = []
y_all = []
z_all = []

# ### Choose inputfiles

# Choose the first file
# Don't enter 'stop' just yet. You ought to select at least one method
method = input("Which file do you want to plot: ").lower()
bestand = "{}.txt".format(method)
while method != 'stop' :

    # ### Read inputfile
    f = open(bestand, "r")
    f.readline()    #first line is whitespace due to setprecision()

    AantalObjecten = int(f.readline()[1:])
    print("Number of objects: "+str(AantalObjecten))
    
    masses_str = f.readline()[1:]    # First character is '#'
    masses = [float(mass) for mass in (masses_str.strip('\t\n')).split('\t')]
    print('Masses: ', masses)
    f.close()
    
    # Load in the time and energy data and place them and the relative energy error in their respective lists.
    # This shouldn't need to be loaded in for every method IF they all used the same timestep. When adaptive 
    # timesteps are used, it is highly unlikely that that condition would be met. 
    load = np.loadtxt(bestand)
    t_all.append(load[:,0])
    E = load[:,1]
    E0 = E[0]
    E_all.append(E)
    Ereldif = (E - E0)/E0
    Ereldif_all.append(Ereldif)
    
    # Add empty lists to the coordinate superlists to house positions arrays of each particle for this method.
    x_all.append([])
    y_all.append([])
    z_all.append([])
    # Add the position coordinates of each particle to the 2D arrays of this method.
    for n in range(AantalObjecten):
        x_all[-1].append(load[:,(n*3+2)])
        y_all[-1].append(load[:,(n*3+3)])
        z_all[-1].append(load[:,(n*3+4)])
    
    # Add the previous method to the list, as that was not yet done.
    method_all.append(method)

    # Ask for the next method or to stop. 
    print(" ")
    print("Enter 'stop' if you do not wish to add more methods.")
    method = input("Which file do you want to plot? (e.g rk4/Verlet/Forest-Ruth): ").lower()
    bestand = "{}.txt".format(method)

print('Selected methods: ' + str(method_all))


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


# ### Make a 3D trajectory plot, energy accuracy plots and a 3D animation for each method.
# j is used as an index here, as i is already defined in the animation section.
for j, method in enumerate(method_all) : 
    
    
    # ### Make 3D figure of the trajectories
    if make3DtrajectoryPlot :
        fig = plt.figure(constrained_layout=True)
        ax = plt.axes(projection ='3d')
        for n in range(AantalObjecten):
            ax.plot3D(x_all[j][n], y_all[j][n], z_all[j][n], label='Particle'+str(n))
            #mlab.plot3d(load[:,n*3+2],load[:,n*3+3],load[:,n*3+4], t)
        #mlab.show()
        #extra dan ax2.plot(t,E)
        plt.title('Trajectories - method: ' + method, fontsize=18)
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        ax.set_zlabel('z', fontsize=16)
        ax.legend()
        #plt.show()
        if equaliseAspect: 
            set_axes_equal(ax)
        fig.savefig('3Dtrajectories_'+method, bbox_inches='tight', dpi=400)



    # ### Make subplots of the relative energy error for each method. Left panel is linear, 
    # right panel is logarithmic.
    if makeSingleEnergyPlots :
        fig, ax = plt.subplots(1, 2, figsize=(16, 6))
        ax[0].plot(t_all[j], Ereldif_all[j])
        ax[0].set_xlabel('t', fontsize=16)
        ax[0].set_ylabel("$\\frac{\Delta E}{E_0}$ (/)", fontsize=16)
        ax[0].set_title('Linear accuracy - '+method, fontsize=18)
        ax[0].grid()
        
        ax[1].plot(t_all[j], np.log10(abs(Ereldif_all[j])))
        ax[1].set_xlabel('t', fontsize=16)
        ax[1].set_ylabel("$\\log_{10}( | \\frac{\Delta E}{E_0} | )$ (/)", fontsize=16)
        ax[1].set_title('Logarithmic accuracy - '+method, fontsize=18)
        ax[1].grid()
        
        #plt.show()
        fig.savefig('Energy errors_'+method, bbox_inches='tight', dpi=400)


    
    if make3Danimations :
        # ### Create an animation of the trajectories
        
        # Select a color from the 'rainbow' colormap for each particle
        color = plt.cm.rainbow(np.linspace(0, 1, AantalObjecten))
        
        # Make a figure and axis for animation.
        fig_anim = plt.figure(constrained_layout=True, figsize=(10, 10))
        #spec = fig.add_gridspec(2, 1)
        
        ax_anim = fig_anim.add_subplot(projection='3d')
        ax_anim.view_init(elev=elevation, azim=azimuth)
        
        # Animations do not choose limits automatically, so they are computed here.
        # We seek to find the largest and smallest x, y and z of any particle in this method.
        # We begin with the minima and maxima of the first particle and then check if the others are smaller or larger respectively.
        x_min = min(x_all[j][0])
        y_min = min(y_all[j][0])
        z_min = min(z_all[j][0])
        x_max = max(x_all[j][0])
        y_max = max(y_all[j][0])
        z_max = max(z_all[j][0])
        for n in range(1, AantalObjecten) :
            if min(x_all[j][n]) < x_min :
                x_min = min(x_all[j][n])
            if min(y_all[j][n]) < y_min :
                y_min = min(y_all[j][n])
            if min(z_all[j][n]) < z_min :
                z_min = min(z_all[j][n])
                
            if max(x_all[j][n]) > x_max :
                x_max = max(x_all[j][n])
            if max(y_all[j][n]) > y_max :
                y_max = max(y_all[j][n])
            if max(z_all[j][n]) > z_max :
                z_max = max(z_all[j][n])
                
        x_dif = x_max - x_min
        y_dif = y_max - y_min
        z_dif = z_max - z_min
        
        if equaliseAspect :
            dif = max([x_dif, y_dif, z_dif])
            x_mean = (x_max+x_min)/2
            y_mean = (y_max+y_min)/2
            z_mean = (z_max+z_min)/2
            x_min = x_mean-dif/2
            x_max = x_mean+dif/2
            y_min = y_mean-dif/2
            y_max = y_mean+dif/2
            z_min = z_mean-dif/2
            z_max = z_mean+dif/2
            
        # Using the minima and maxima as limits would mean the trajectories would go right to the edges of the animation, which does not look tidy.
        # As such, we add a small margin - defined above - to either side.
        x_lower = x_min - relativeMargin*x_dif - additiveMargin
        y_lower = y_min - relativeMargin*y_dif - additiveMargin
        z_lower = z_min - relativeMargin*z_dif - additiveMargin
        x_upper = x_max + relativeMargin*x_dif + additiveMargin
        y_upper = y_max + relativeMargin*y_dif + additiveMargin
        z_upper = z_max + relativeMargin*z_dif + additiveMargin
        
        
        # Set the computed limits
        ax_anim.set(xlabel='x', ylabel='y', zlabel='z', 
                    xlim=(x_lower, x_upper), ylim=(y_lower, y_upper), zlim=(z_lower, z_upper))
        
        # Make "line" elements to update every frame.
        lines = []
        points = []
        for n in range(AantalObjecten):
            lines.append(ax_anim.plot([], [], [], color=color[n])[0])
            points.append(ax_anim.plot([], [], [], color=color[n], marker=".", ms=5)[0])
        
        
        def animate(i):
            """
            animate function is called every frame i.
            """
            
            index = int(i*stepsPerPoint)
            time = t_all[j][index]
            # print("time: {} \r".format(time), end="")
            
            fig_anim.suptitle("time: {}".format(time))    # Title changes every frame to show time.
            
            # Loop over all objects and change the figure.
            for n in range(AantalObjecten):
                lines[n].set_data(x_all[j][n][:index], y_all[j][n][:index])
                lines[n].set_3d_properties(z_all[j][n][:index])
                
                points[n].set_data(x_all[j][n][index], y_all[j][n][index])
                points[n].set_3d_properties(z_all[j][n][index])
            
            return fig_anim,
        
        
        print("Calculating frames...")

        if (endTime > t_all[j][-1]) :
            endTime = t_all[j][-1]
        endIndex = 0
        while (t_all[j][endIndex] < endTime) :
            endIndex += 1
        anim = animation.FuncAnimation(fig_anim, animate, frames=[i for i in range(round(endIndex/stepsPerPoint))], interval=1000/fps, blit=True) #50 geeft 20fps
        HTML(anim.to_jshtml())
        
        name = "3Danimation_"+ str(method) + ".mp4"
        anim.save(name)
        
           


# ### Make a single energy accuracy plot containing all methods for easy comparison
if makeCombinedEnergyPlot :
    fig, ax = plt.subplots(1, 2, figsize=(16, 6))
    for j, method in enumerate(method_all) :
        ax[0].plot(t_all[j], Ereldif_all[j], label=method)
        ax[1].plot(t_all[j], np.log10(abs(Ereldif_all[j])), label=method_all[j])
    
    ax[0].set_xlabel('t', fontsize=16)
    ax[0].set_ylabel("$\\frac{\Delta E}{E_0}$ (/)", fontsize=16)
    ax[0].set_title('Linear accuracy', fontsize=18)
    ax[0].grid()
    ax[0].legend()
    
    ax[1].set_xlabel('t', fontsize=16)
    ax[1].set_ylabel("$\\log_{10}( | \\frac{\Delta E}{E_0} | )$ (/)", fontsize=16)
    ax[1].set_title('Logarithmic accuracy', fontsize=18)
    ax[1].grid()
    ax[1].legend()
    
    #plt.show()
    fig.savefig('Energy errors_all', bbox_inches='tight', dpi=400)
    
    
    
 # Plot the total deviation of the centre of mass of the system from zero.
if makeCOMdeviationPlot : 
    
    # Make a figure that all methods will be plotted in
    fig, ax = plt.subplots(2, 2, figsize=(16, 12))
    
    for j, method in enumerate(method_all) :
        # Ik maak bewust drie van die dingen indien we plots van individuele componenten willen maken.
        xdev = []
        ydev = []
        zdev = []
        
        # Compute the three coordinates of the centre of mass for each timestep.
        for k in range(len(x_all[j][0])) :
            x = 0
            y = 0
            z = 0
            for n in range(AantalObjecten) :
                x += masses[n] * x_all[j][n][k]
                y += masses[n] * y_all[j][n][k]
                z += masses[n] * z_all[j][n][k]
            xdev.append(x)
            ydev.append(y)
            zdev.append(z)
        xdev = np.array(xdev)
        ydev = np.array(ydev)
        zdev = np.array(zdev)
        COMdev = np.sqrt(xdev**2 + ydev**2 + zdev**2)
    
        # Plot the deviations of this method in the earlier defined figure
        ax[0][0].plot(t_all[j], COMdev, label=method_all[j])
        ax[0][1].plot(t_all[j], xdev, label=method_all[j])
        ax[1][0].plot(t_all[j], ydev, label=method_all[j])
        ax[1][1].plot(t_all[j], zdev, label=method_all[j])
        
    # Labelling and formatting of the subplots
    ax[0][0].set_xlabel('t', fontsize=14)
    ax[0][0].set_ylabel('$r_{CoM}$', fontsize=14)
    ax[0][0].grid()
    ax[0][0].legend(fontsize=12)
    
    ax[0][1].set_xlabel('t', fontsize=14)
    ax[0][1].set_ylabel('$x_{CoM}$', fontsize=14)
    ax[0][1].grid()
    ax[0][1].legend(fontsize=12)

    ax[1][0].set_xlabel('t', fontsize=14)
    ax[1][0].set_ylabel('$y_{CoM}$', fontsize=14)
    ax[1][0].grid()
    ax[1][0].legend(fontsize=12)
    
    ax[1][1].set_xlabel('t', fontsize=14)
    ax[1][1].set_ylabel('$z_{CoM}$', fontsize=14)
    ax[1][1].grid()
    ax[1][1].legend(fontsize=12)
    
    fig.savefig('COMdeviations_all', dpi=400, bbox_inches='tight')


if make_E_ifv_h_plot:

    #prepare figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    ax.set(title = "E error ifv h", xlabel = "h", ylabel = "dE / E0", xscale="log", yscale = "log")
    ax.invert_xaxis()

    # load and plot file for each method
    for j, method in enumerate(method_all) :
        load = np.loadtxt("dE_ifv_h_" + method + ".txt", unpack=True)
        ax.plot(load[0], load[1], label=method)

    ax.legend()
    fig.savefig('dE_ifv_h', dpi=400, bbox_inches='tight')
    
    
if make_E_ifv_driverfunctions:
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    ax.set(title = "E error ifv number of driver function evaluations", 
           xlabel = "driver function evaluations", ylabel = "dE / E0", 
           yscale = "log")

    # load and plot file for each method
    for j, method in enumerate(method_all) :
        load = np.loadtxt("dE_ifv_h_" + method + ".txt", unpack=True)
        ax.plot(load[2], load[1], label=method)

    ax.legend()
    fig.savefig('dE_ifv_driverfunctions', dpi=400, bbox_inches='tight')
    
    
if makeCombined2DtrajectoryPlot :
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.grid()
    for j, method in enumerate(method_all) :
        for n in range(AantalObjecten) :
            plt.plot(x_all[j][n], y_all[j][n], label=method+' '+str(n+1))
    if equaliseAspect :
        ax.set_aspect('equal')
    plt.xlabel('x (a.u.)', fontsize=16)
    plt.ylabel('y (a.u.)', fontsize=16)
    plt.legend(fontsize=12)
    fig.savefig('combined2Dtrajectory', bbox_inches='tight', dpi=400)
    
if makeSingle2DtrajectoryPlots :
    for j, method in enumerate(method_all) :
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        ax.grid()
        for n in range(AantalObjecten) :
            plt.plot(x_all[j][n], y_all[j][n], label='particle'+' '+str(n+1))
        if equaliseAspect :
            ax.set_aspect('equal')
        plt.xlabel('x (a.u.)', fontsize=16)
        plt.ylabel('y (a.u.)', fontsize=16)
        plt.title(method, fontsize=18)
        fig.savefig('2Dtrajectory_'+method, bbox_inches='tight', dpi=400)
    