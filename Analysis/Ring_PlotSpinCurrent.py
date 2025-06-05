import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import scienceplots

from numpy import sin, cos, sqrt, arange


filename = 'table.txt'#'spincurrent_cmap_10deg_nodmi.txt'
timestep = 1e-10

data = np.loadtxt(filename, unpack=True)

rnum=250
reg = list(range(rnum))
mRAll = []
mzAll = []
for j in reg:
    ang = (j-1)*(360/rnum)*(np.pi/180)

    #find mr for a single region at all times
    mx = data[4+3*j]
    my = data[5+3*j]
    mz = data[6+3*j]
    mRx = []
    mRy = []
    mR = []
    for i in arange(len(data[0])):
        mR.append(-mx[i]*sin(ang)+my[i]*cos(ang))
        
    mRAll.append(mR)
    mzAll.append(mz)
    


pos = [item*(360/rnum) for item in reg]
X, Y = np.meshgrid(data[0], reg)
X=X*pow(10,9)
Y=Y*(360/rnum)


with plt.style.context(['science','no-latex']):

    fig = plt.figure(figsize=(8, 6), constrained_layout=True)
    gs = gridspec.GridSpec(1, 2, width_ratios=[20, 1], figure=fig)

    ax = fig.add_subplot(gs[0, 0])
    im = ax.pcolor(X, Y, mRAll)

    ax.set_xlabel("Time (ns)", fontsize=22)
    ax.set_ylabel("Position along ring (degrees)", fontsize=22)
    ax.tick_params(axis='both', labelsize=16)

    # Add colorbar in its own subplot so it's not clipped
    cax = fig.add_subplot(gs[0, 1])
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.set_ylabel("Normalized radial component of magnetization", rotation=90, labelpad=10, fontsize=18)
    cbar.ax.tick_params(labelsize=16)

    fig.savefig("Ring phase plot_inverse.png", dpi=600, bbox_inches='tight', transparent=False, facecolor='white')

    fig, ax = plt.subplots(constrained_layout=True)
    plt.pcolor(X, Y, mRAll)
    plt.xlabel("Time (ns)")
    plt.ylabel("Position along ring (degrees)")
    cbar = plt.colorbar()
    cbar.ax.set_ylabel("Normalized radial component of magnetization")
    plt.savefig("Ring phase plot.png", dpi=600, bbox_inches='tight',transparent='False', facecolor='white')
    plt.show()
