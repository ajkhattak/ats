#!/usr/bin/env python
"""
Plot water table as a function of time.

This is currently only useful on a single, 1D column.
"""

import h5py
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm
import parse_ats, transect_data, column_data
    

def water_table_no_surface(dirname, datum=None, patm=101325.):
    # get the columnar pressure
    col_dat = column_data.column_data(['pressure',], directory=dirname)
    wt = np.zeros((col_dat.shape[1],),'d')
    if datum is None:
        datum = column_data[0,0,-1]
    
    for k in range(col_dat.shape[1]):
        if col_dat[1,k,-1] > patm:
            wt[k] = (col_dat[1,k,-1]-patm) / 1000. / 9.81 - datum  # hydrostatic above
        else:
            wt_index = np.where(col_dat[1,k,:] >= patm)[0]
            if len(wt_index) is 0:
                wt_index = 0
            else:
                wt_index = wt_index[-1]
            wt[k] = col_dat[0,k,wt_index] - datum
    return times, wt



def water_table(dirname, datum=None, v86=False, patm=101325.):
    prefix = "" if v86 else "surface-"

    # get the ponded depth
    keys,times,dats = parse_ats.readATS(dirname, "visdump_surface_data.h5")
    pd = parse_ats.get2DSurfaceData(keys, dats, prefix+"ponded_depth")

    center = len(dats[prefix+"elevation.cell.0"][keys[0]])
    print 'Center location: ', center
    
    elev_surf = dats[prefix+"elevation.cell.0"][keys[0]][center-1]
    dats.close()
    print ('Datum (surface elevation)', elev_surf)
    if datum is None:
        datum = elev_surf
    datum_offset = elev_surf - datum

    # get the columnar pressure
    wt = np.zeros((len(keys),),'d')
    wt_ss = np.zeros((len(keys),),'d')
    wt_surf = np.zeros((len(keys),),'d')
    thaw_depth = np.zeros((len(keys),),'d')
    wt_bottom = np.zeros((len(keys),),'d')
    
    col_dat = transect_data.transect_data(['pressure', 'temperature'], directory=dirname)
        
    vars = 2
    nvar, cycles, xnum, znum = col_dat.shape
    for k in range(col_dat.shape[1]):

        if pd[k, xnum-1] > 0:
            wt[k] = pd[k, xnum-1] + datum_offset
            wt_surf[k] = pd[k, xnum-1] + datum_offset
        else:
            wt_index = np.where(col_dat[vars ,k,xnum-1, :] >= patm)[0] #var 1
            
            temp_index = np.where(col_dat[vars+1 ,k,xnum-1, :] >= 273.25)[0] ##
            temp_data = col_dat[vars + 1,k,xnum-1, :] #var 2

            wt_index0 = wt_index
            
            if len(wt_index) is 0:
                wt_index = 0
                wt_index0 = 0
            else:
                wt_index = wt_index[-1]

            if (len(wt_index0) > 2):
                wt_index0 = wt_index0[-2]
                
            if (len(temp_index) is 0):
                temp_index =0
            else:
                temp_index = temp_index[0]
            
            if (temp_data[wt_index] > 273.25):
                wt[k] = col_dat[1,k,xnum-1, wt_index] - datum
            else:
                wt_ss[k] = col_dat[1,k,xnum-1, temp_index] - datum

                
        temp_index = np.where(col_dat[vars+1 ,k,xnum-1, :] >= 273.25)[0] ##

        if len(temp_index) > 0:
            thaw_depth[k] = col_dat[1,k,xnum-1, temp_index[0]] - datum
        
    return times, wt, wt_ss, wt_bottom, thaw_depth


def water_table2D(dirname, datum=None, location='center',x_pos=-1, v86=False, patm=101325.):
    prefix = "" if v86 else "surface-"

    # get the ponded depth
    keys,times,dats = parse_ats.readATS(dirname, "visdump_surface_data.h5")
    pd = parse_ats.get2DSurfaceData(keys, dats, prefix+"ponded_depth")

    if location == 'center':
        center = len(dats[prefix+"elevation.cell.0"][keys[0]])
    elif location == 'trough':
        center = 1
    if x_pos > 0:
        center = x_pos
    
    print 'Location: ', location, center
    
    elev_cell = dats[prefix+"elevation.cell.0"][keys[0]][center-1]
    dats.close()
    print ('Datum (surface elevation)', elev_cell)
    if datum is None:
        datum = elev_cell
    datum_offset = elev_cell - datum

    # get the columnar pressure
    water_table = np.zeros((len(keys),),'d')
    thaw_depth = np.zeros((len(keys),),'d')

    
    col_dat = transect_data.transect_data(['pressure', 'temperature'], directory=dirname)
        
    vars = 2
    nvar, cycles, xnum, znum = col_dat.shape
    #print xnum
    xnum = center
    #print xnum
    for k in range(col_dat.shape[1]):

        if pd[k, xnum-1] > 0:
            water_table[k] = pd[k, xnum-1] + datum_offset
            #wt_surf[k] = pd[k, xnum-1] + datum_offset
            #print 'surface:: ',k, pd[k, xnum-1]
        else:
            #print 'subsurface: ', k, pd[k, xnum-1]
            wt_index = np.where(col_dat[vars ,k,xnum-1, :] >= patm)[0] #var 1
            
            temp_index = np.where(col_dat[vars+1 ,k,xnum-1, :] >= 273.25)[0] ##
            temp_data = col_dat[vars + 1,k,xnum-1, :] #var 2

            if len(wt_index) is 0:
                wt_index = 0
            else:
                wt_index = wt_index[-1]

            if (len(temp_index) is 0):
                temp_index =0
            else:
                temp_index = temp_index[0]
            
            if (temp_data[wt_index] > 273.25):
                water_table[k] = col_dat[1,k,xnum-1, wt_index] - datum
                
        temp_index = np.where(col_dat[vars+1 ,k,xnum-1, :] >= 273.2)[0] ##
        #print col_dat[vars+1 ,k,xnum-1, :], round(times[k]*365.25,2), k
        if len(temp_index) > 0:
            #print 'ALT: ', len(temp_index), temp_index[0], col_dat[vars+1 ,k,xnum-1, temp_index[0]]

            thaw_depth[k] = col_dat[1,k,xnum-1, temp_index[0]] - datum
        
    return times, water_table, thaw_depth, elev_cell

def plot_water_table(times, wt, ax, **kwargs):
    ax.plot(times, wt, **kwargs)

def get_axs():
    f,axs = plt.subplots(1,1)
    return f,axs
        
if __name__ == "__main__":
    import sys

    import argparse
    import colors
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("INPUT_DIRECTORIES", nargs="+", type=str,
                        help="List of directories to plot.")
    parser.add_argument("--colors", "-c", type=colors.float_list_type,
                        default=None,
                        help="List of color indices to use, of the form: --colors=[0,0.1,1], where the doubles are in the range (0,1) and are mapped to a color using the colormap.")
    parser.add_argument("--colormap", "-m", type=str,
                        default="jet",
                        help="Colormap used to pick a color.")
    args = parser.parse_args()

    cm = colors.cm_mapper(0,1,args.colormap)
    fig, axs = get_axs()

    fnames = args.INPUT_DIRECTORIES
    for i,fname in enumerate(fnames):
        if args.colors is None:
            if len(fnames) > 1:
                c = cm(float(i)/(len(fnames)-1))
            else:
                c = 'b'
        else:
            if type(args.colors[i]) is float:
                c = cm(args.colors[i])
            else:
                c = args.colors[i]

        times,wt = water_table(fname)
        plot_water_table(times, wt, axs, color=c, label=fname)

    plt.tight_layout()
    axs.legend()
    plt.show()
