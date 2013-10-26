import sys
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from scipy.interpolate import griddata

def parse_parameters(lines):
    params = {"VERSION":float,
              "DIMENSIONS":int,
              "TIMESTEP":float,
              "ITERATIONS":int,
              "INTERVAL BETWEEN GAUSSIANS":int,
              "GAUSSIAN SIGMA":float,
              "GAUSSIAN W":float}
    parsedParams = {}
              
    for i in xrange(len(lines)):
        for key in params.keys():
            if key in lines[i]:
                parsedParams[key] = params[key](lines[i].split()[-1].rstrip("\n"))
    return parsedParams

def parse_trajectory(lines):
    trajectory = []
    energy = []
    for i in xrange(1, len(lines)):
        trajectory.append(map(float, lines[i].split()[1:]))
    return trajectory
    
def parse_gaussians(lines):
    tstep = []
    mu = []
    for i in xrange(1, len(lines)):
        splitLine = lines[i].split()
        tstep.append(int(splitLine[0]))
        mu.append(map(float, splitLine[1:]))
    return [tstep, mu]
    
def parse_landscape(lines, dims):
    # coords, trueEnergy and gaussianEnergies all have the same length
    coords = []
    trueEnergy = []
    gaussianEnergies = []
    for i in xrange(1, len(lines)):
        splitLine = lines[i].split()
        coords.append(map(float, splitLine[:dims]))
        trueEnergy.append(float(splitLine[dims]))
        gaussianEnergies.append(map(float, splitLine[dims+1:]))
    return [coords, trueEnergy, gaussianEnergies]

def parse_energies(lines):
    energies = []
    for i in xrange(1, len(lines)):
        energies.append(float(lines[i]))
    return energies
        
def parse_output(flist):
    parsed = {}
    i=0
    # Get the starts and ends of the sections
    dividers = [i for i in xrange(len(flist)) if flist[i] == "--------------------\n"]
    print dividers
    for i in xrange(len(dividers)-1):
        section = flist[dividers[i]+1:dividers[i+1]]
        #print section
        if section[0] == "PARAMETERS:\n":
            parsed["parameters"] = parse_parameters(section)
        elif section[0] == "TRAJECTORY (first column is timestep):\n":
            parsed["trajectory"] = parse_trajectory(section) 
        elif section[0] == "ENERGIES:\n":
            parsed["energies"] = parse_energies(section)
        elif section[0] == "GAUSSIANS' LOCATIONS:\n":
            parsed["gaussians"] = parse_gaussians(section)
        elif section[0] == "ENERGY LANDSCAPE (format: (coords), (true energy), (cumulative energy from metadynamics gaussians)):\n":
            parsed["landscape"] = parse_landscape(section, parsed["parameters"]["DIMENSIONS"])
    return parsed
    
def generate_mesh(coords, height):
    #TODO: possibly extend to non-square matrices? might be useful if a coordinate remains effectively constant
    #TODO: replace w/ numpy.reshape for optimum speed?
    size = int(math.sqrt(len(coords)))
    X = [[None for i in xrange(size)] for j in xrange(size)]
    Y = [[None for i in xrange(size)] for j in xrange(size)]
    Z = [[None for i in xrange(size)] for j in xrange(size)]
    for i in xrange(len(coords)):
        X[int(i/size)][i%size] = coords[i][0]
        Y[int(i/size)][i%size] = coords[i][1]
        Z[int(i/size)][i%size] = height[i]
        #print height[i]
    return np.array(X), np.array(Y), np.array(Z)
    
            

if __name__ == "__main__":
    fname = "C:\Users\Jamie\Documents\Physics\Rare event simulation\meta\meta\output.txt"
    f = open(fname)
    flist = list(f)
    parsed = parse_output(flist)
    parameters = parsed["parameters"]
    if parameters["DIMENSIONS"]==2:
        # Get energy surface
        coords     = parsed["landscape"][0]
        energy     = parsed["landscape"][1]
        metaEnergy = parsed["landscape"][2]
        finalEnergy = [e[-1] for e in metaEnergy]
        
        # Cut off extremes on the energy landscape
        eMax = 2E2
        for i in xrange(len(energy)):
            if energy[i] > eMax:
                energy[i] = eMax

        # Get trajectory
        trajectory = parsed["trajectory"]
        xTrajectory = np.array([r[0] for r in trajectory])
        yTrajectory = np.array([r[1] for r in trajectory])  
        zTrajectory = parsed["energies"]
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # 3D energy landscape
        #X = np.array([r[0] for r in coords])
        #Y = np.array([r[1] for r in coords])
        #xi, yi = np.mgrid[min(X):max(X):100j, min(Y):max(Y):100j]
        #xi = np.linspace(min(X), max(X), int(math.sqrt(len(X))))
        #yi = np.linspace(min(Y), max(Y), int(math.sqrt(len(Y))))
        #zi = griddata(coords, energy, (xi[None,:], yi[:,None]))
        #Xm, Ym = np.meshgrid(xi, yi)
        #ax.plot_surface(Xm, Ym, zi, rstride=1, cstride=1, cmap=cm.coolwarm, antialiased=False)
        
        X, Y, Z = generate_mesh(coords, energy)
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, antialiased=True)
        metaX, metaY, metaZ = generate_mesh(coords, finalEnergy)
        metaOffset = -5*ones_like(metaX)
        ax.plot_surface(metaX, metaY, metaOffset-metaZ, rstride=1, cstride=1, cmap=cm.coolwarm, antialiased=True)
        
        # Get energy as a function of trajectory via interpolation
#        print "Calculating trajectory"
#        zTrajectory = []
#        for i in xrange(len(xTrajectory)):
#            rx = np.array([xTrajectory[i]])
#            ry = np.array([yTrajectory[i]])
#            zTrajectory.append(griddata(np.array(coords), np.array(energy), (rx[None,:], ry[:,None]), method="linear")[0][0])
#        print "Done."
        ax.plot(xTrajectory, yTrajectory, zTrajectory, c='g')
        # Place markers at start and end
        ax.plot([xTrajectory[0]], [yTrajectory[0]], [zTrajectory[0]], c='b', marker='o', markersize=10)
        ax.plot([xTrajectory[-1]], [yTrajectory[-1]], [zTrajectory[-1]], c='r', marker='o', markersize=10)
        
        # Plot Gaussian locations
        tGauss = parsed["gaussians"][0]
        for i in xrange(len(tGauss)):
            ax.plot([xTrajectory[tGauss[i]]], [yTrajectory[tGauss[i]]], [zTrajectory[tGauss[i]]], c='y', marker='+', markersize=5)
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        
        fig.savefig('harmonic', dpi=600)
        
        #fig2d = plt.figure()
        #ax2d = fig2d.add_subplot(111)
        #ax2d.contour(xi, yi, zi)
    f.close()
    
        
    
    