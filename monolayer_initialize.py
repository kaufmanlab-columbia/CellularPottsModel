import matplotlib.pyplot as plt
from random import shuffle
import numpy as np
import scipy as sp
import sys
import os
from shapely.geometry import Polygon
from shapely.geometry import Point
import math
import argparse


eps = sys.float_info.epsilon

#### Functions ####

def check_value(a,L):
    if a<0:
        a=a+L
    elif a>L:
        a=a-L
    return a

def center_placement(L,n):
    x = np.random.random_sample(size = n)*L
    y = np.random.random_sample(size = n)*L
    return np.array(list(zip(x,y)))

def reflect_centers(points,L):
    all_points = []
    all_points.append(points)
    
    x = points[:,0]
    y = points[:,1]
    
    all_points.append(list(zip(-x,y)))
    
    all_points.append(list(zip(x,-y)))
    
    x_reflect = (2*L)-x
    all_points.append(list(zip(x_reflect,y)))
    
    y_reflect = (2*L)-y
    all_points.append(list(zip(x,y_reflect)))
    
    return np.array(all_points)

def flatten(t):
    return np.array([item for sublist in t for item in sublist])

def check_range(a,L):
    if a<0:
        a=0
    elif a>L:
        a=L-1
    else:
        pass
    return a

def randomizeCells(cells):
    l = list(range(1,cells+1))
    shuffle(l)
    return l

def readData(filepath):
    data = np.loadtxt(filepath)
    return data

def modifyLattice(d,shuffled, N, L):
    placeholder = np.zeros((L,L))
    for i in range(1,N+1):
        values = np.argwhere(d==i)
        for j in values:
            placeholder[j[0],j[1]] = shuffled[i-1]
    
    return placeholder


def in_box(cell_centers, bounding_box):
    return np.logical_and(np.logical_and(bounding_box[0] <= cell_centers[:, 0],
                                         cell_centers[:, 0] <= bounding_box[1]),
                          np.logical_and(bounding_box[2] <= cell_centers[:, 1],
                                         cell_centers[:, 1] <= bounding_box[3]))


def voronoi(cell_centers, bounding_box):
    i = in_box(cell_centers, bounding_box)
    # Mirror points
    points_center = cell_centers[i, :]
    points_left = np.copy(points_center)
    points_left[:, 0] = bounding_box[0] - (points_left[:, 0] - bounding_box[0])
    points_right = np.copy(points_center)
    points_right[:, 0] = bounding_box[1] + (bounding_box[1] - points_right[:, 0])
    points_down = np.copy(points_center)
    points_down[:, 1] = bounding_box[2] - (points_down[:, 1] - bounding_box[2])
    points_up = np.copy(points_center)
    points_up[:, 1] = bounding_box[3] + (bounding_box[3] - points_up[:, 1])
    points = np.append(points_center,
                       np.append(np.append(points_left,
                                           points_right,
                                           axis=0),
                                 np.append(points_down,
                                           points_up,
                                           axis=0),
                                 axis=0),
                       axis=0)
    # Compute Voronoi
    vor = sp.spatial.Voronoi(points)
    # Filter regions
    regions = []
    for region in vor.regions:
        flag = True
        for index in region:
            if index == -1:
                flag = False
                break
            else:
                x = vor.vertices[index, 0]
                y = vor.vertices[index, 1]
                if not(bounding_box[0] - eps <= x and x <= bounding_box[1] + eps and
                       bounding_box[2] - eps <= y and y <= bounding_box[3] + eps):
                    flag = False
                    break
        if region != [] and flag:
            regions.append(region)
    vor.filtered_points = points_center
    vor.filtered_regions = regions
    return vor

def centroid_region(vertices):
    #Polygon's signed area
    A = 0
    # Centroid's x
    C_x = 0
    # Centroid's y
    C_y = 0
    for i in range(0, len(vertices) - 1):
        s = (vertices[i, 0] * vertices[i + 1, 1] - vertices[i + 1, 0] * vertices[i, 1])
        A = A + s
        C_x = C_x + (vertices[i, 0] + vertices[i + 1, 0]) * s
        C_y = C_y + (vertices[i, 1] + vertices[i + 1, 1]) * s
    A = 0.5 * A
    C_x = (1.0 / (6.0 * A)) * C_x
    C_y = (1.0 / (6.0 * A)) * C_y
    return np.array([[C_x, C_y]])


def tessellate_monolayer(N,L,show_plot=False):
	cell_centers = center_placement(L,N)/L 
	bounding_box = np.array([-0.01, 1.01, -0.01, 1.01]) # [x_min, x_max, y_min, y_max]

	vor = voronoi(cell_centers, bounding_box)

	if show_plot:
		fig = plt.figure()
		ax = fig.gca()
		# Plot initial points
		ax.plot(vor.filtered_points[:, 0], vor.filtered_points[:, 1], 'b.')
		# Plot ridges points
		for region in vor.filtered_regions:
		    vertices = vor.vertices[region, :]
		    ax.plot(vertices[:, 0], vertices[:, 1], 'go')
		# Plot ridges
		for region in vor.filtered_regions:
		    vertices = vor.vertices[region + [region[0]], :]
		    ax.plot(vertices[:, 0], vertices[:, 1], 'k-')
		# Compute and plot centroids
		centroids = []
		for region in vor.filtered_regions:
			vertices = vor.vertices[region + [region[0]], :]
			centroid = centroid_region(vertices)
			centroids.append(list(centroid[0, :]))
			ax.plot(centroid[:, 0], centroid[:, 1], 'r.')
	plt.show()


	# Get ridges for defining shapes
	ridges = []
	for region in vor.filtered_regions:
	    vertices = vor.vertices[region + [region[0]], :]
	    ridges.append(vertices)
	    
	ridges = np.array(ridges,dtype=object)

	# Convert to tuples for use with polygon shapely
	l = []
	for i in ridges:
	    l.append(tuple(i*L)) 


	# Create polygons
	p = []
	for i in l:
	    p.append(Polygon(i))


	# Initialize empty lattice
	lattice = np.zeros(L**2).reshape(L,L) #100= lattice dim

	# Get min and max for x and y to get a bounding rectangle to narrow the coordinate search
	# Check which points lie in shapes and assign numbers
	for n,i in enumerate(p):

	    c = np.array(i.exterior.coords)
	    x_max = check_range(np.max(c[:,0]),L) #100= lattice dim
	    x_min = check_range(np.min(c[:,0]),L)
	    y_max = check_range(np.max(c[:,1]),L)
	    y_min = check_range(np.min(c[:,1]),L)
	    

	    for j in range(math.floor(x_min),math.ceil(x_max)+1):
	        for k in range(math.floor(y_min),math.ceil(y_max)+1):
	            point = Point(j,k)
	            if i.contains(point)==True:
	                lattice[j][k] = n+1



	# Randomize cells spins and export
	shuffled = randomizeCells(N)
	output = modifyLattice(lattice,shuffled, N, L)
	np.savetxt("lattice_N_"+str(N)+"_L_"+str(L)+".txt",lattice,fmt = '%d')


# Set output lattice directory
os.chdir(r'/Users/alexdevanny/Downloads')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-N", type=int, help="Specify the value for N")
    parser.add_argument("-L", type=int, help="Specify the value for L")
    args = parser.parse_args()

    if args.N is not None and args.L is not None:
        # N and L are specified explicitly
        tessellate_monolayer(args.N, args.L,True) # true for plot, False for no plot
    else:
        # Read N and L from the command line
        N = int(input("Enter the value for N: "))
        L = float(input("Enter the value for L: "))
        tessellate_monolayer(N, L,True)









 

















