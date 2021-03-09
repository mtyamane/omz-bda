"""
Mark Yamane, 3/2/2021

Helper functions used in analyzing DO data from the WOA18 dataset, including
an algorithm used to identify boundaries of OMZ cores.
"""

from netCDF4 import Dataset
import numpy as np

def pullData(filename):
    fpath = 'data/'
    ds = Dataset(fpath+filename)
    
    # pull data
    lons = ds.variables['lon'][:]
    lats = ds.variables['lat'][:]
    deps = ds.variables['depth'][:]
    dO = ds.variables['o_an'][0][:][:][:]
    
    return lons, lats, deps, dO

def getLatSlice(filename, target_lat):
    # pull data
    lons, lats, deps, dO = pullData(filename)

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
    
    # center at Pacific
    for i in range(len(lonVert)):
        for j in range(len(lonVert[0,:])):
            if lonVert[i,j] < 0:
                lonVert[i,j] += 360
    # reorder data (-180 - 180  --> 0 - 360)
    temp = 0
    for i in range(len(lonVert)):
        for j in range(180):
            temp = lonVert[i, j+180]
            lonVert[i, j+180] = lonVert[i,j]
            lonVert[i, j] = temp
            temp = dO_slice[i, j+180]
            dO_slice[i, j+180] = dO_slice[i,j]
            dO_slice[i,j] = temp
    
    # mask to keep only lons of the Pacific
    left = 125
    right = 260
    pacific_mask = (lonVert[0] > left) & (lonVert[0] < right)
    new_lon = []
    new_dep = []
    new_dO = []
    for i in range(len(lonVert)):
        new_lon.append(lonVert[i][pacific_mask])
        new_dep.append(deps[i][pacific_mask])
        new_dO.append(dO_slice[i][pacific_mask])
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