#!/usr/local/bin/python
# routines to compute strain or displacement greens functions and related matrices from okada dislocations
# Eric Lindsey, 01/2014
#
# from okada85.py:
# ue,un,uz =  displacement(e,n,depth,strike,dip,L,W,rake,slip,open,nu)
# uze,uzn =           tilt(e,n,depth,strike,dip,L,W,rake,slip,open,nu)
# unn,une,uen,uee = strain(e,n,depth,strike,dip,L,W,rake,slip,open,nu)
#

import sys
import numpy as np
import okada85
import geod_transform

#======================================================================
def displacement_greens(lat,lon,lat0,lon0,depth,strike,dip,L,W,nu=0.25):
    # arguments are lists or 1D arrays of all fault patches to include in G
    # returns numpy array
    # each set of 3 rows corresponds to [E,N,U] displacements for a single obs. point
    # each set of 2 columns corresponds to unit [left-lateral strike,along-dip thrust] displacements on a fault patch

    #observation altitudes assumed zero by default
    alt=0.*np.array(lon)

    nobs,npatch = check_lengths(lat,lon,lat0,lon0,depth,strike,dip,L,W)

    #allocate G matrix: simplest is 3rows, 2cols: [[e1_strike1 e1_dip1], [n1_strike1, n1_dip1], [u1_strike1,u1_dip1]]
    G=np.zeros((3*nobs,2*npatch))
    
    #get greens functions for strike and dip cases, these appear as columns in G
    for ipatch in range(npatch):
        # convert lat,lon to local (e,n) coordinates (in meters) relative to patch center lat0[ipatch],lon0[ipatch]
        e,n,u = geod_transform.geod2enu(lat,lon,alt,lat0[ipatch],lon0[ipatch],0.)

        # okada subroutine calls -- turns out float() is important
        strcolE,strcolN,strcolU=okada85.displacement(e, n, float(depth[ipatch]), float(strike[ipatch]), float(dip[ipatch]), float(L[ipatch]), float(W[ipatch]),  0., 1., 0., nu)
        dipcolE,dipcolN,dipcolU=okada85.displacement(e, n, float(depth[ipatch]), float(strike[ipatch]), float(dip[ipatch]), float(L[ipatch]), float(W[ipatch]), 90., 1., 0., nu)
        #interleave arrays so that each 3 observations for a given site are kept together.
        gstr = np.zeros(3*nobs)
        gdip = np.zeros(3*nobs)
        gstr[::3]  = strcolE
        gstr[1::3] = strcolN
        gstr[2::3] = strcolU
        gdip[::3]  = dipcolE
        gdip[1::3] = dipcolN
        gdip[2::3] = dipcolU
        G[:,2*ipatch  ]=gstr
        G[:,2*ipatch+1]=gdip
        
    return G

#======================================================================
def strain_greens(lat,lon,lat0,lon0,depth,strike,dip,L,W,nu=0.25):
    # greens function matrix using strains at specified points instead of displacements

    #observation altitudes assumed zero by default
    alt=0.*np.array(lon)

    nobs,npatch = check_lengths(lat,lon,lat0,lon0,depth,strike,dip,L,W)

    #allocate G matrix -- now there are 4 observations
    G=np.zeros((4*nobs,2*npatch))
    
    #get greens functions for strike and dip cases, these appear as columns in G
    for ipatch in range(npatch):
        # convert lat,lon to local (e,n) coordinates (in meters) relative to fault center lat0,lon0
        e,n,u = geod_transform.geod2enu(lat,lon,alt,lat0[ipatch],lon0[ipatch],0.)

        # okada subroutine calls -- turns out float() is important
        strcolNN,strcolNE,strcolEN,strcolEE = okada85.strain(e, n, float(depth[ipatch]), float(strike[ipatch]), float(dip[ipatch]), float(L[ipatch]), float(W[ipatch]),  0., 1., 0., nu)
        dipcolNN,dipcolNE,dipcolEN,dipcolEE = okada85.strain(e, n, float(depth[ipatch]), float(strike[ipatch]), float(dip[ipatch]), float(L[ipatch]), float(W[ipatch]),  0., 1., 0., nu)
        #interleave arrays so that each 3 observations for a given site are kept together.
        gstr = np.zeros(4*nobs) 
        gdip = np.zeros(4*nobs) 
        gstr[::4]  = strcolNN
        gstr[1::4] = strcolNE 
        gstr[2::4] = strcolEN
        gstr[3::4] = strcolEE
        gdip[::4]  = dipcolNN
        gdip[1::4] = dipcolNE
        gdip[2::4] = dipcolEN
        gdip[3::4] = dipcolEE
        G[:,2*ipatch  ]=gstr
        G[:,2*ipatch+1]=gdip
    return G

#======================================================================
def resolution(G):
    # compute resolution matrix R = pseudoinverse(G) * G, to double precision
    #R = np.dot(np.dot(np.linalg.inv(np.dot(G.T, G)), G.T),G)
    R = np.dot(np.linalg.pinv(G), G)

    return R

#======================================================================
def check_lengths(lat,lon,lat0,lon0,depth,strike,dip,L,W):
    # inputs must all be arrays/lists or this will fail
    nobs = len(lat)
    npatch = len(strike)
    #ensure two groups of lists are equal length
    if (len(set([nobs,len(lon)]))>1):
        print("Error: lists lat,lon are not the same length")
        sys.exit(1)
    if (len(set([npatch,len(lat0),len(lon0),len(depth),len(dip),len(L),len(W)]))>1):
        print("Error: lists lat0,lon0,depth,strike,dip,L,W are not the same length")
        sys.exit(1)
    return nobs, npatch

