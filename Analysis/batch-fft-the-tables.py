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

Bases fft from 2D_colorplot_from_mumaxtable_with_slicer_v3.py
Directory crawling and other features pulled from batch-output-ao-plotter_amplitude.py

Code looks at a number of directories looped over using batch mumax3 scripts
used in hpc3 for faster turnaround.

It goes into each directory with it's table.txt, does fft of it, and spits out
an fft_table.txt in the same directory.

Another python script will iterate through these to do 2D plots and whatnot.

"""

# Filename related parameters, these are for regex
# Use regexr.com for help in getting proper regex
searchtoken1 = 'j2_([+-]?[\d.\d]+)'
searchtoken2 = 'j1_([+-]?[\d.\d]+)'
token1_var = 'j2'
token2_var = 'j1'
directory_template = 'YIG_SHO*.out' #'Easy_Plane_SHO_v9-10-2f_jc_-z*.out'
file_to_fft = 'table.txt'

# Apply Blackman filter to reduce spectral leakage, check wikipedia for more info
debug = True                        # Additional printouts
use_bm_filter = True                # If false applies no window
interp_timestep = 10e-11            # Should set same as timestep in mumax3
skip_steps = 0e-9/interp_timestep  # To cut off transients if table from t=0
column_to_fft = 1                  # Column with m(t) data (0 indexed)
column_of_components = [1, 2, 3] # Columns of x, y, z data to process for Mvect

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
        if not os.path.exists(execpath + '/' + dirs + '/' + file_to_fft):
            print(dirs, 'Does not have', file_to_fft, 'skipping!')
            continue
        directory.append((dirs, cmpkey(dirs, searchtoken1), cmpkey(dirs, searchtoken2)))

    except Exception as e:
        print('Pulling directory vars failed!\n Reason: ', e)
        pass

for l in directory:
	print(l)
input('Is the list correct?')



xyz_data = []
d_index = 0
# Loads in the data that we want to fft
for d in directory:
    filename = execpath + '/' + d[0] + '/' + file_to_fft
    print("FFTing: ", filename)
    # Extract column name using pandas subroutine
    column_name = get_column_label(filename, column_to_fft)

    # loaded in as data[0] = t, data[1] = mx, etc..
    data = np.loadtxt(filename, unpack=True)
    tslice = data[0]
    mslice = interp1d(tslice, data[column_to_fft])
    mslice_squared = interp1d(tslice, np.square(data[column_to_fft]))

    mvect = np.sqrt(np.square(data[column_of_components[0]]) + np.square(data[column_of_components[1]]) + np.square(data[column_of_components[2]]))

    t_steps = int((tslice[-1] - tslice[0])/interp_timestep)
    t = np.linspace(tslice[0], tslice[-1], t_steps)
    # Debug
    if debug:
        print(tslice[0], tslice[-1], len(tslice), t_steps, len(mslice(t)), '\t', max(mslice(t)))

    FFTlength = len(mslice(t)) - int(skip_steps) if(int(skip_steps) < len(mslice(t))) else len(mslice(t))
    FFTstart = int(skip_steps) if (int(skip_steps) < len(mslice(t))) else 0
    print('FFT length (N-k):', FFTlength)
    print('Skipping (k) steps:', FFTstart, '=', interp_timestep*FFTstart, 'seconds')

    # Applies blackman spectral filter
    window = np.blackman(FFTlength) if use_bm_filter else 1
    mslice_remove_offset = mslice(t) - np.mean(mslice(t))
    mslice_squared_rem_off = mslice_squared(t) - np.mean(mslice_squared(t))

    # Generate linear and log scale fft, applies window, takes absolute value, keeps positive frequency portion of fft output.
    fft_linear = np.abs(np.fft.fft(mslice_remove_offset[-FFTlength:]*window))*(2.0/FFTlength)
    fft_linear = fft_linear[0:int(FFTlength/2)]
    fft_log = np.log10(np.abs(np.fft.fft(mslice_remove_offset[-FFTlength:]*window))*(2.0/FFTlength))
    fft_log = fft_log[0:int(FFTlength/2)]
    fft_squared_linear = np.abs(np.fft.fft(mslice_squared_rem_off[-FFTlength:]*window))*(2.0/FFTlength)
    fft_squared_linear = fft_squared_linear[0:int(FFTlength/2)]
    fft_squared_log = np.log10(np.abs(np.fft.fft(mslice_squared_rem_off[-FFTlength:]*window))*(2.0/FFTlength))
    fft_squared_log = fft_squared_log[0:int(FFTlength/2)]

    # Generate frequency axis in units of GHz
    frequencies = np.fft.fftfreq(FFTlength, d = interp_timestep)[0:int(FFTlength/2)]/1e9
    filename = execpath + '/' + d[0] + '/' + 'fft_' + file_to_fft
    print('Saved fft to: ', filename)
    # Save fft data as a new table
    datahdr = ('Frequency(GHz)',
               column_name + '_fft_linear',
               column_name + '_fft_log',
               column_name + '_squared_fft_linear',
               column_name + '_squared_fft_log')
    np.savetxt(filename, np.transpose(np.array([frequencies,
                                               fft_linear,
                                               fft_log,
                                               fft_squared_linear,
                                               fft_squared_log])),
    newline = '\n', header = '\t'.join(map(str, datahdr)))

    datahdr = ('t(s)', 'Mvect')
    filename = filename = execpath + '/' + d[0] + '/' + 'Mvect_table.txt'
    print('Saved mvect to: ', filename)
    np.savetxt(filename, np.transpose(np.array([tslice, mvect])),
               newline = '\n', header = '\t'.join(map(str, datahdr)))

