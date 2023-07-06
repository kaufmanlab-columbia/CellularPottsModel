# Converts the CPM lattice into an mp4 movie

import os
import time
from os import listdir
from os.path import isfile, join
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import numpy as np
from PIL import Image

trimmed=True

# The index of the first output lattice file used to generate a movie
start_step = 0

# The index of the final output lattice file used to generate a movie
end_step = 99

# The length of the edge of the square lattice
lattice_edge = 100

# The number of cells in the simulation
# (only used to scale color values so there's a good range)
num_cells = 100

# The interval between timesteps that are printed in lattice files (e.g. if printing every 10 timesteps)
step_interval = 1

# The number of types of cells in the simulation
# Currently only works with 0 or 1, but could be extended by modifying the
# lattice_to_frame() function appropriately.
num_cell_types = 1

# Set to true to save a .mp4 file of the animation
save_animation = True

# Set to true to show the animation in a matplotlib window
show_animation = False

# Set to true to print progress messages to the terminal
verbose = True

# Array to store data read in from simulation output lattice files
sim_array = []

# Number of frames to skip when generating an animation
# (e.g. skip = 10 will skip 9 frames for every 1 that's animated
skip = 1

sim_array = []

# print(os.path.dirname('.'))
# data_path = 'movie_output/'
# data_path = join("C:/", "Users/", "Microscope PC/", "Documents/", "Lucas/", "2Types (original)_468-10A/", "output/")

def main():
    read_args()
    read_lattice_files()
    # lattices_to_animation()
    lattices_to_FuncAnimation()

def read_args():
    if (len(sys.argv) > 0):
        data_path = sys.argv[0]
        if (len(sys.argv) > 1):
            end_step = sys.argv[1]
            if (len(sys.argv) > 2):
                start_step = sys.argv[2]
                end_step = sys.argv[1] + sys.argv[2]
                if (len(sys.argv) > 3):
                    print("Too many arguments")

def read_lattice_files():
    global end_step
    global trimmed

    if trimmed==True:
        full_lattice = np.load(r'output\\lattice_select_perim.npy')
    else:
        full_lattice = np.load(r'output\\lattice_condensed.npy')

    for file in range(start_step*lattice_edge,(end_step+1)*lattice_edge,skip*lattice_edge):
        curr_arr = full_lattice[file:file+lattice_edge]
        if (verbose):
            print("reading in file: {}".format(int(file/lattice_edge)))
        curr_arr_list = [list(i) for i in curr_arr]
        sim_array.append(curr_arr_list)
    if (verbose):
        print('Lattice files read in successfully.')

# Image generation

def lattice_to_frame(step_index, _ax, colors):
    # Data for plotting
    if (verbose):
        print('Generating RGB array for timestep {0}\r'.format(step_index)),

    ID = 0
    frame = np.zeros((lattice_edge, lattice_edge, 3))
    for y in range(lattice_edge):
        for x in range(lattice_edge):
            ID = sim_array[step_index][x][y]
            frame[x,y] = colors[ID]
    return frame

# After data is/are read in from individual lattice files into the sim_array that stores all of this information,
# timepoints are processed sequentially into RGB arrays and then animated using FuncAnimation.
#
# After that, depending on settings, the plot can be saved and/or displayed.
def lattices_to_FuncAnimation():
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_ylim(0, lattice_edge)
    # colors = { i: [0.83,0.83,0.83] for i in range(1, 101)}
    # colors[79] = [1,0,0]
    # colors[45] = [1,0,0]
    # colors[48] = [1,0,0]
    # colors[32] = [1,0,0]
    # colors[51] = [1,0,0]
    # colors[85] = [1,0,0]
    # colors[26] = [1,0,0]
    # colors[99] = [1,0,0]
    colors = { i: np.random.rand(3) for i in range(1,101) }
    colors[101] = [0.33,0.33,0.33]


    frames = [lattice_to_frame(i, ax,colors) for i in range(0, len(sim_array), 1)]

    # Function used by the FuncAnimation class to generate the first frame of the video
    def init():
        return plt.imshow(frames[0], animated=True)

    # Function used by the FuncAnimation class to generate each subsequent frame of the video
    def animate(nframe):
        plt.cla()
        if (verbose):
            print("Generating artist frame number ", nframe, "of ", len(frames))
        im = plt.imshow(frames[nframe], animated=True)
        # plt.title('timestep = {}'.format((nframe + start_step) * step_interval * skip))

        ax.tick_params(
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False) # labels along the bottom edge are off
        plt.yticks([])
        # plt.gca().invert_yaxis()

    if (verbose):
        print("Generating animation.")
    ani = animation.FuncAnimation(fig, animate, frames=(int((end_step-start_step)/skip)), interval=40)  # , repeat_delay=1000

    if (save_animation):
        if (verbose):
            print("Saving animation to .mp4")
        # Writer options: imagemagick, ffmpeg
        # Imagemagick mp4 can only be opened in certain programs (not ppt), use ffmpeg as writer for more reliable saving
        ani.save('{0}_simulation_{1}_steps.mp4'.format(time.strftime("%Y%m%d-%H%M%S"), (end_step - start_step)),writer="ffmpeg",fps=10)

    if (show_animation):
        if (verbose):
            print("Displaying animation.")
        plt.show()

# Run the main method.
main()
