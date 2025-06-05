from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import os
import glob

"""
Amanatullah 'Mishu' Khan

Directory crawling and other features pulled from batch-fft-the-tables.py

Code looks at a number of directories looped over using batch mumax3 scripts
used in hpc3 for faster turnaround.

It goes into each directory, reads in the fft_table.txt, and processes as the user desires, in this case it will generate color plots and flattened 3d maps.

pseudo code is as follows:
scan for appropriate directories
for each dir:
    read in the fft data
    map it to vars extracted from directory
    map flattened data to same vars as a separate tuple

for each set of tuples made from the data crawling:
    for each fixed var1 (or var2):
        generate 2d colorplot for each var2 (or var1)
        save fftmatrix in new folder for each fixed var1 (or var2)

    generate flattened data of all colormaps together (also a colormap)
    save this flattened data matrix in same directory
    generated frequency extracted version of flattened data
    save this data matrix as well

To do:
Address cplot z scale to truncate below min cutoff frequency (now only cuts off above max)
Save raw matrix data alongside png's for magicplot import (easier to adjust and analyze)

"""
############################### Directory Inputs ##############################
# Filename related parameters, these are for regex
# Use regexr.com for help in getting proper regex
searchtoken1 = 'j2_([+-]?[\d.\d]+)'
searchtoken2 = 'j1_([+-]?[\d.\d]+)'
token1_var = 'j2'
token1_var_units = 'A/m^2'
token2_var = 'j1'
token2_var_units = 'A/m^2'
directory_template = 'YIG_SHO*.out'
file_to_process = 'fft_table.txt'

column_to_process = 2                  # Column with fft data (0 indexed)

############################ Plotting Inputs ##################################
# 2D colorplots from fixing each var1/var2
y_axis_range_max = 4.0
y_axis_range_min = 0.1
number_of_cplot_levels = 100

############################
# 2D 'flattened' colorplot from extracing max amplitude from every fft
flatten_threshold = -5

title = 'AO phaseplot from max amp across J1(wire) and J2(injector)'
z_label = 'Log Amplitude (a.u.)'
z_label_from_yaxis = 'Frequency (GHz)'
# If this is false it'll do 3d plot
contour_plot = True

# Removes this many rows from beginning of data before flattening
# Useful mainly when removing dc leakage in fft (artificial maximum at 0 f)
truncate_range = 1

# If true, returns the y axis where the max is found (like f in fft)
# when false, simply extracts the maximum value across slice
flatten_to_yaxis = False

############################### End User Inputs ###############################

def cmpkey(text, token):
    return float(re.findall(token, text)[0])

searched_directory_list = glob.glob(os.path.join(os.curdir,
                                                 directory_template))
execpath = os.getcwd()


def get_column_label(path_to_file, index_of_column):
    # var is named fields, this is the "h" part of usual fh plot
    column_extract = pd.read_table(path_to_file, nrows = 3)
    # get rid of the first # and space
    column_extract.columns = column_extract.columns.str.replace('# ', '')
    col_label = column_extract.columns[index_of_column].replace(r' \(.*\)','').replace(' ()','').replace('.','') # cleanup
    col_label = re.sub(r'\(.*\)','',col_label) # cleanup 2 since above one is dumb af
    return col_label



directory = []
for dirs in searched_directory_list:
    try:
        if not os.path.exists(execpath + '/' + dirs + '/' + file_to_process):
            print(dirs, 'Does not have', file_to_process, 'skipping!')
            continue
        directory.append((dirs, cmpkey(dirs, searchtoken1), cmpkey(dirs, searchtoken2)))

		# directory.append((dirs, cmpkey(dirs.replace('p', '.'))))
        # print(list(map(float, re.findall(searchtoken,
        #                                 dirs.replace('p', '.')))))
    except:
        print('Encountered an error trying to find\t', dirs)
        pass

for l in directory:
	print(l)
input('Is the list correct?')



data_matrix = []
xyz_flattened = []
d_index = 0
# Loads in the data that we want to fft
for d in directory:
    filename = execpath + '/' + d[0] + '/' + file_to_process
    # Extract column name using pandas subroutine
    column_name = get_column_label(filename, column_to_process)
    print("FFTing:\t", filename)
    print("Extracting:\t", column_name)
    # loaded in as data[0] = t, data[1] = mx, etc..
    data = np.loadtxt(filename, unpack=True)
    frequency = data[0]
    max_from_slice = max(data[column_to_process][truncate_range:])
    # Flatten to just raw maximum
    if not flatten_to_yaxis:
        flattened_data = (d[1], d[2], max_from_slice)
    # Otherwise flatten to the frequency where max is found
    # If max amplitude is below threshold, set to 0
    else:
        if max_from_slice < flatten_threshold:
            flattened_data = (d[1], d[2], 0)
        else:
            flattened_data = (d[1], d[2], data[0][np.where(data[column_to_process] == max_from_slice)[0]])
    xyz_flattened.append(flattened_data)
    fft_tuple = (d[1], d[2], frequency, data[column_to_process])
    data_matrix.append(fft_tuple)
    if d_index == -1:
        break   # For debugging
    d_index += 1

# Iterating over tuple list to extract axes
var1_axis = np.unique([stuff[0] for stuff in data_matrix])
var2_axis = np.unique([stuff[1] for stuff in data_matrix])

pd_data = pd.DataFrame(data_matrix, columns=[token1_var, token2_var, "frequency", "fft"])
f_pd = pd.DataFrame(frequency, columns=['frequency'])
print(pd_data)

# Makes colorplots with fixed var1's
for var1 in var1_axis:
    temp_xyz = []
    cplot_data = pd_data[pd_data[token1_var] == var1]
    print(cplot_data)
    frequency_cut = len(f_pd[f_pd['frequency'] < y_axis_range_max])
    fftmatrix = np.vstack([d for d in cplot_data["fft"]]).T[:frequency_cut]
    frequency = frequency[:frequency_cut]
    # For cplot
    fft_max = fftmatrix.max()
    fft_min = fftmatrix.min()
    fft_step = round((fft_max-fft_min)/number_of_cplot_levels, 20)
    lvls = np.arange(fft_min, fft_max+fft_step, fft_step)
    # Plotting
    plt.subplot(1, 1, 1)
    print(cplot_data[token2_var].shape, frequency.shape, fftmatrix.shape)
    cplot = plt.contourf(cplot_data[token2_var].to_numpy(), frequency, fftmatrix, levels=lvls)
    plt.ylim([y_axis_range_min, y_axis_range_max])
    plt.xlabel(token2_var + '(%s)' % (token2_var_units))
    plt.ylabel('frequency (GHz)')
    plt.title('%s with %s fixed at %s %s' % (column_name, token1_var, str(var1), token1_var_units))
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Amplitude')
    cplot_image_name = 'YIG_SF*.out'.replace('*', '%s_%s' % (token1_var, str(var1)))
    plt.savefig(cplot_image_name.replace('.out', '.png'))
    print('Saved ', cplot_image_name.replace('.out', '.png'))
    # Cleanup for next plot
    plt.clf()


# NOTE MAKE SURE THE PD IS ORDERING BY VAR2 IN THIS LOOP
# Makes colorplots with fixed var2's
for var2 in var2_axis:
    temp_xyz = []
    cplot_data = pd_data[pd_data[token2_var] == var2]
    #print(cplot_data)
    frequency_cut = len(f_pd[f_pd['frequency'] < y_axis_range_max])
    fftmatrix = np.vstack([d for d in cplot_data["fft"]]).T[:frequency_cut]
    frequency = frequency[:frequency_cut]
    # For cplot
    fft_max = fftmatrix.max()
    fft_min = fftmatrix.min()
    fft_step = round((fft_max-fft_min)/number_of_cplot_levels, 20)
    lvls = np.arange(fft_min, fft_max+fft_step, fft_step)
    # Plotting
    plt.subplot(1, 1, 1)
    #print(cplot_data[token1_var].shape, frequency.shape, fftmatrix.shape)
    cplot = plt.contourf(cplot_data[token1_var].to_numpy(), frequency, fftmatrix, levels=lvls)
    plt.ylim([y_axis_range_min, y_axis_range_max])
    plt.xlabel(token1_var + '(%s)' % (token1_var_units))
    plt.ylabel('frequency (GHz)')
    plt.title('%s with %s fixed at %s %s' % (column_name, token2_var, str(var2), token2_var_units))
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Amplitude')
    cplot_image_name = 'YIG_SF*.out'.replace('*', '%s_%s' % (token2_var, str(var2)))
    plt.savefig(cplot_image_name.replace('.out', '.png'))
    print('Saved ', cplot_image_name.replace('.out', '.png'))
    # Cleanup for next plot
    plt.clf()

# Flattened colorplot plotting
x, y, z = zip(*xyz_flattened)
z = list(map(float, z))

grid_x, grid_y = np.mgrid[min(x):max(x):100j, min(y):max(y):100j]
grid_z = griddata((x, y), z, (grid_x, grid_y), method='nearest')

for i in range(len(z)):
    if z[i] >= -2.5 and z[i] <= -2.3:
        z[i] = 1
    else:
        z[i] = 0
x = list(map(float, x))
y = list(map(float, y))
for i in range(len(x)): x[i] = x[i]/10000000000
for i in range(len(y)): y[i] = y[i]/10000000000

np.savetxt("ScikitData.txt", np.transpose(np.array([x,y,z])), newline = '\n', delimiter = ',', fmt='%f')

# Makes and saves it.
fig = plt.figure()
if contour_plot:
    ax = fig.add_subplot(111)
    grid_extent = [grid_x.min(), grid_x.max(), grid_y.min(), grid_y.max()]
    cplot = ax.imshow(grid_z.T, extent=grid_extent, aspect='auto', origin='lower')
    cbar = fig.colorbar(cplot, ax=ax)
    cbar.set_label(z_label if not flatten_to_yaxis else z_label_from_yaxis)
else:
    ax = fig.gca(projection='3d')
    ax.plot_surface(grid_x, grid_y, grid_z, cmap=plt.cm.Spectral)

plt.title(title)
plt.xlabel('%s (%s)' % (token1_var, token1_var_units))
plt.ylabel('%s (%s)' % (token2_var, token2_var_units))
if not flatten_to_yaxis:
    cplot_image_name = 'YIG_SF*.out'.replace('*', '%s_%s_%s' % (token1_var, token2_var, 'FLATTENED'))
else:
    cplot_image_name = 'YIG_SF*.out'.replace('*', '%s_%s_%s' % (token1_var, token2_var, 'FLATTENED_yaxis'))
plt.savefig(cplot_image_name.replace('.out', '.png'))

plt.show()
