"""
Mark Yamane, 3/2/2021

Helper functions used in analyzing DO data from the WOA18 dataset, including
an algorithm used to identify boundaries of OMZ cores.
"""

from netCDF4 import Dataset
import numpy as np

def pullData(fpath, fname):
    ds = Dataset(fpath+fname)
    
    # pull data
    lons = ds.variables['lon'][:]
    lats = ds.variables['lat'][:]
    deps = ds.variables['depth'][:]
    dO = ds.variables['o_an'][0][:][:][:]
    
    return lons, lats, deps, dO

def centerAtPacific(lons, ys, zs):
    # centers lons and zs to the Pacific Ocean and
    # pulls data for the specified longitudinal range
    
    # center global map at Pacific
    for i in range(len(lons)):
        for j in range(len(lons[0,:])):
            if lons[i,j] < 0:
                lons[i,j] += 360
    # reorder data (-180 - 180  --> 0 - 360)
    temp = 0
    for i in range(len(lons)):
        for j in range(180):
            temp = lons[i, j+180]
            lons[i, j+180] = lons[i,j]
            lons[i, j] = temp
            temp = zs[i, j+180]
            zs[i, j+180] = zs[i,j]
            zs[i,j] = temp
    
    # mask to keep only lons of the Pacific Ocean
    left = 125
    right = 260
    pacific_mask = (lons[0] > left) & (lons[0] < right)
    new_lon = []
    new_ys = []
    new_zs = []
    for i in range(len(lons)):
        new_lon.append(lons[i][pacific_mask])
        new_ys.append(ys[i][pacific_mask])
        new_zs.append(zs[i][pacific_mask])
    return new_lon, new_ys, new_zs

def getLatSlice(fpath, fname, target_lat):
    # pull data
    lons, lats, deps, dO = pullData(fpath, fname)

    # create full grids for each dimension
    tempLat = []
    tempLon = []
    tempDep = []

    nLons = len(lons)
    for lat in lats:
        tempLat.append(np.full(nLons, lat))
        tempLon.append(lons)
    for lon in lons:
        tempDep.append(deps)

    # cast to numpy arrays for easier handling
    lats = np.array(tempLat)
    lons = np.array(tempLon)
    deps = np.array(tempDep)
    
    # find index of latitude to perform the slice on
    lat_lo = lats[0,0]
    i_lat = int((target_lat - lat_lo)/1.)
    actual_lat = lats[i_lat,0]
    
    # create full grid for a latitudinal slice w/ depth
    lonVert = []
    for dep in deps[0,:]:
        lonVert.append(lons[i_lat,:])
    
    # format data for plot
    lonVert = np.array(lonVert)
    deps = np.transpose(deps)
    dO_slice = dO[:, i_lat, :]
    
    new_lon, new_dep, new_dO = centerAtPacific(lonVert, deps, dO_slice)
    return i_lat, actual_lat, new_lon, new_dep, new_dO

def omzbda(lons, deps, dO, threshold):
    top = 1500
    bottom = 0
    left = 360
    right = 0
    omzb_t = ()
    omzb_b = ()
    omzb_l = ()
    omzb_r = ()
    for i in range(len(dO)):
        for j in range(len(dO[0])):
            if dO[i][j] <= threshold:
                # in the OMZ core, so identify if bounds are present
                if deps[i][j] < top:
                    top = deps[i][j]
                    omzb_t = (lons[i][j], deps[i][j])
                if deps[i][j] > bottom:
                    bottom = deps[i][j]
                    omzb_b = (lons[i][j], deps[i][j])
                if lons[i][j] < left:
                    left = lons[i][j]
                    omzb_l = (lons[i][j], deps[i][j])
                if lons[i][j] > right:
                    right = lons[i][j]
                    omzb_r = (lons[i][j], deps[i][j])
                
    return top, bottom, left, right, omzb_t, omzb_b, omzb_l, omzb_r

def RUN_OMZBDA(i_fpath, i_fnames, target_lats, o_fpath, o_fname, status = True):
    # runs OMZ-BDA algorithm on all files from i_fnames and outputs to o_fname
    month = 1
    if status:
        print('Algorithm start...')
    with open(o_fpath+o_fname, 'x') as f:
        f.write('month,lat,top,bottom,west,east\n')
        for i_fname in i_fnames:
            if status:
                print('working on month:', month)
            for target_lat in target_lats:
                i_lat, actual_lat, lons, deps, dO = getLatSlice(i_fpath, i_fname, target_lat)
                top, bottom, left, right, t_pt, b_pt, l_pt, r_pt = omzbda(lons, deps, dO, 20)
                f.write(','.join([str(month), str(actual_lat), str(top), str(bottom), str(left), str(right)]) + '\n')
            month += 1
    f.close()
    if status:
        print('... algorithm finished. Output is located at: ' + o_fpath + o_fname)