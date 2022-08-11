# okada85.py
# python routines for calculating surface displacements and strains
# due to a finite rectangular dislocation source.

            #conversion ToDo
            # R appears in every single function... reuse?
#
# Requirements: numpy (numerical python package)
#
# Included routines:
# ue,un,uz =  displacement(e,n,depth,strike,dip,L,W,rake,slip,open,nu)
# uze,uzn =           tilt(e,n,depth,strike,dip,L,W,rake,slip,open,nu)
# unn,une,uen,uee = strain(e,n,depth,strike,dip,L,W,rake,slip,open,nu)
#
# e,n             : surface coordinates relative to the fault centroid
# depth           : depth of the fault centroid (depth > 0)
# strike          : in degrees from North, 90 points East.
# dip             : in degrees from horizontal, to the right side of the trace
# rake            : in degrees, direction the hanging wall moves. If slip > 0,
#                 : a rake of 0 is left-lateral slip, +90 is reverse fault.
# L,W             : length (along-strike) and width (along-dip) of the fault patch
# slip, open      : displacements, uniform across entire patch
# nu              : poisson's ratio (optional, default 0.25)
#
# Units of e,n,depth,L,W,slip,open should be the same (eg. everything in meters)
#
# For strains, note
#    POSITIVE = COMPRESSION
# Note that vertical strain components can be obtained with following equations:
#    uNZ = -uZN;
#    uEZ = -uZE;
#    uZZ = -(uEE + uNN)*NU/(1-NU);
#
# References:
#   Aki K., and P. G. Richards, Quantitative seismology, Freemann & Co,
#       New York, 1980.
#   Okada Y., Surface deformation due to shear and tensile faults in a
#       half-space, Bull. Seismol. Soc. Am., 75:4, 1135-1154, 1985.
#
# =================================================================
# Licensing / Copyright
#
# Copyright (c) 2014, Eric Lindsey, covered by BSD License (see text below).
# All rights reserved.
# 
# Version history:
# [01/2014] Converted to python from [08/2012] matlab version by Francois Beauducel:
# available http://www.ipgp.fr/~beaudu/matlab.html
#
# Matlab version:
# Copyright (c) 1997-2012, Francois Beauducel, covered by BSD License (see text below).
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are 
# met:
# 
#    * Redistributions of source code must retain the above copyright 
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright 
#      notice, this list of conditions and the following disclaimer in 
#      the documentation and/or other materials provided with the distribution
#                            
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE.
# =================================================================

import numpy as np

# =================================================================
def setup_args(e,n,depth,strike,dip,L,W,rake,slip,open):

    # ensure array-like behavior
    e=np.array(e)
    n=np.array(n)
    depth=np.array(depth)
    L=np.array(L)
    W=np.array(W)
    slip=np.array(slip)
    open=np.array(open)

    #convert to radians and ensure array-like behavior
    strike = np.array(strike*np.pi/180)
    dip = np.array(dip*np.pi/180)
    rake = np.array(rake*np.pi/180)

    # Defines dislocation in the fault plane system
    # 's' denotes slip scaled by the pre-multiplier 2*pi
    U1s = np.cos(rake)*slip/(2*np.pi)
    U2s = np.sin(rake)*slip/(2*np.pi)
    U3s = open/(2*np.pi)
    
    # Converts fault coordinates (E,N,DEPTH) relative to centroid
    # into Okada's reference system (X,Y,D)
    d = depth + np.sin(dip)*W/2    # d is fault's top edge
    ec = e + np.cos(strike)*np.cos(dip)*W/2
    nc = n - np.sin(strike)*np.cos(dip)*W/2
    x = np.cos(strike)*nc + np.sin(strike)*ec + L/2
    y = np.sin(strike)*nc - np.cos(strike)*ec + np.cos(dip)*W
    
    # Variable substitution (independent from xi and eta)
    p = y*np.cos(dip) + d*np.sin(dip)
    q = y*np.sin(dip) - d*np.cos(dip)

    return x,p,L,W,q,strike,dip,U1s,U2s,U3s


# =================================================================
def displacement(e,n,depth,strike,dip,L,W,rake,slip,open,nu=0.25):

    #convert coordinate systems, radians, ensure array_like, etc.
    x,p,L,W,q,strike,dip,U1s,U2s,U3s = setup_args(e,n,depth,strike,dip,L,W,rake,slip,open)

    ux = ( -U1s * chinnery(ux_ss,x,p,L,W,q,dip,nu)   # strike-slip
           -U2s * chinnery(ux_ds,x,p,L,W,q,dip,nu)   # dip-slip
           +U3s * chinnery(ux_tf,x,p,L,W,q,dip,nu) ) # tensile fault

    uy = ( -U1s * chinnery(uy_ss,x,p,L,W,q,dip,nu)   # strike-slip
           -U2s * chinnery(uy_ds,x,p,L,W,q,dip,nu)   # dip-slip
           +U3s * chinnery(uy_tf,x,p,L,W,q,dip,nu) ) # tensile fault

    uz = ( -U1s * chinnery(uz_ss,x,p,L,W,q,dip,nu)   # strike-slip
           -U2s * chinnery(uz_ds,x,p,L,W,q,dip,nu)   # dip-slip
           +U3s * chinnery(uz_tf,x,p,L,W,q,dip,nu) ) # tensile fault

    # Rotation from Okada's axes to geographic
    ue = np.sin(strike)*ux - np.cos(strike)*uy
    un = np.cos(strike)*ux + np.sin(strike)*uy

    # function call to pick these up looks like dE,dN,dU = displacement(...)
    return ue,un,uz


# =================================================================
def tilt(e,n,depth,strike,dip,L,W,rake,slip,open,nu=0.25):

    #convert coordinate systems, radians, ensure array_like, etc.
    x,p,L,W,q,strike,dip,U1s,U2s,U3s = setup_args(e,n,depth,strike,dip,L,W,rake,slip,open)
    
    uzx = ( -U1s * chinnery(uzx_ss,x,p,L,W,q,dip,nu)   # strike-slip
            -U2s * chinnery(uzx_ds,x,p,L,W,q,dip,nu)   # dip-slip
            +U3s * chinnery(uzx_tf,x,p,L,W,q,dip,nu) ) # tensile fault

    uzy = ( -U1s * chinnery(uzy_ss,x,p,L,W,q,dip,nu)   # strike-slip
            -U2s * chinnery(uzy_ds,x,p,L,W,q,dip,nu)   # dip-slip
            +U3s * chinnery(uzy_tf,x,p,L,W,q,dip,nu) ) # tensile fault

    # Rotation from Okada's axes to geographic
    uze = -np.sin(strike)*uzx + np.cos(strike)*uzy
    uzn = -np.cos(strike)*uzx - np.sin(strike)*uzy

    return uze,uzn


# =================================================================
def strain(e,n,depth,strike,dip,L,W,rake,slip,open,nu=0.25):
    # positive = compression

    #convert coordinate systems, radians, ensure array_like, etc.
    x,p,L,W,q,strike,dip,U1s,U2s,U3s = setup_args(e,n,depth,strike,dip,L,W,rake,slip,open)
    
    uxx = ( -U1s * chinnery(uxx_ss,x,p,L,W,q,dip,nu)   # strike-slip
            -U2s * chinnery(uxx_ds,x,p,L,W,q,dip,nu)   # dip-slip
            +U3s * chinnery(uxx_tf,x,p,L,W,q,dip,nu) ) # tensile fault

    uxy = ( -U1s * chinnery(uxy_ss,x,p,L,W,q,dip,nu)   # strike-slip
            -U2s * chinnery(uxy_ds,x,p,L,W,q,dip,nu)   # dip-slip
            +U3s * chinnery(uxy_tf,x,p,L,W,q,dip,nu) ) # tensile fault
            
    uyx = ( -U1s * chinnery(uyx_ss,x,p,L,W,q,dip,nu)   # strike-slip
            -U2s * chinnery(uyx_ds,x,p,L,W,q,dip,nu)   # dip-slip
            +U3s * chinnery(uyx_tf,x,p,L,W,q,dip,nu) ) # tensile fault

    uyy = ( -U1s * chinnery(uyy_ss,x,p,L,W,q,dip,nu)   # strike-slip
            -U2s * chinnery(uyy_ds,x,p,L,W,q,dip,nu)   # dip-slip
            +U3s * chinnery(uyy_tf,x,p,L,W,q,dip,nu) ) # tensile fault

    # Rotation from Okada's axes to geographic
    unn =  np.square(np.cos(strike))*uxx + np.sin(2*strike)*(uxy + uyx)/2 + np.square(np.sin(strike))*uyy
    une =  np.square(np.sin(strike))*uyx + np.sin(2*strike)*(uxx - uyy)/2 - np.square(np.cos(strike))*uxy
    uen = -np.square(np.cos(strike))*uyx + np.sin(2*strike)*(uxx - uyy)/2 + np.square(np.sin(strike))*uxy
    uee =  np.square(np.sin(strike))*uxx - np.sin(2*strike)*(uyx + uxy)/2 + np.square(np.cos(strike))*uyy

    #note, une, uen should be the same?
    return unn,une,uen,uee


# =================================================================
# Chinnery's notation [equation (24) p. 1143]
#
# -----------------------------------------------------------------
def chinnery(f,x,p,L,W,q,dip,nu):
    u = f(x,p,q,dip,nu) - f(x,p-W,q,dip,nu) - f(x-L,p,q,dip,nu) + f(x-L,p-W,q,dip,nu)
    return u


# =================================================================
# Displacement subfunctions

# strike-slip displacement subfunctions [equation (25) p. 1144]
# -----------------------------------------------------------------
def ux_ss(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    u = xi*q/(R*(R + eta)) + np.arctan(xi*eta/(q*R)) + I1(xi,eta,q,dip,nu,R)*np.sin(dip)
    return u

# -----------------------------------------------------------------
def uy_ss(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    u = (eta*np.cos(dip) + q*np.sin(dip))*q/(R*(R + eta)) + q*np.cos(dip)/(R + eta) + I2(eta,q,dip,nu,R)*np.sin(dip)
    return u

# -----------------------------------------------------------------
def uz_ss(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    db = eta*np.sin(dip) - q*np.cos(dip)
    u = (eta*np.sin(dip) - q*np.cos(dip))*q/(R*(R + eta)) + q*np.sin(dip)/(R + eta) + I4(db,eta,q,dip,nu,R)*np.sin(dip)
    return u

# dip-slip displacement subfunctions [equation (26) p. 1144]
# -----------------------------------------------------------------
def ux_ds(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    u = q/R - I3(eta,q,dip,nu,R)*np.sin(dip)*np.cos(dip)
    return u

# -----------------------------------------------------------------
def uy_ds(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    u = (eta*np.cos(dip) + q*np.sin(dip))*q/(R*(R + xi)) + np.cos(dip)*np.arctan(xi*eta/(q*R)) - I1(xi,eta,q,dip,nu,R)*np.sin(dip)*np.cos(dip)
    return u

# -----------------------------------------------------------------
def uz_ds(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    db = eta*np.sin(dip) - q*np.cos(dip)
    u = db*q/(R*(R + xi)) + np.sin(dip)*np.arctan(xi*eta/(q*R)) - I5(xi,eta,q,dip,nu,R,db)*np.sin(dip)*np.cos(dip)
    return u

# tensile fault displacement subfunctions [equation (27) p. 1144]
# -----------------------------------------------------------------
def ux_tf(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    u = np.square(q) /(R*(R + eta)) - I3(eta,q,dip,nu,R)*np.square(np.sin(dip))
    return u

# -----------------------------------------------------------------
def uy_tf(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    u = -(eta*np.sin(dip) - q*np.cos(dip))*q/(R*(R + xi)) - np.sin(dip)*(xi*q/(R*(R + eta)) - np.arctan(xi*eta/(q*R))) - I1(xi,eta,q,dip,nu,R)*np.square(np.sin(dip))
    return u

# -----------------------------------------------------------------
def uz_tf(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    db = eta*np.sin(dip) - q*np.cos(dip)
    u = (eta*np.cos(dip) + q*np.sin(dip))*q/(R*(R + xi)) + np.cos(dip)*(xi*q/(R*(R + eta)) - np.arctan(xi*eta/(q*R))) - I5(xi,eta,q,dip,nu,R,db)*np.square(np.sin(dip))
    return u


# I... displacement subfunctions [equations (28) (29) p. 1144-1145]
# -----------------------------------------------------------------
def I1(xi,eta,q,dip,nu,R):
    db = eta*np.sin(dip) - q*np.cos(dip)
    I = np.where(np.cos(dip) > np.spacing(1),
        (1 - 2*nu) * (-xi/(np.cos(dip)*(R + db))) - np.sin(dip)/np.cos(dip) *I5(xi,eta,q,dip,nu,R,db),
        -(1 - 2*nu)/2 * xi*q/np.square(R + db) )
    return I

# -----------------------------------------------------------------
def I2(eta,q,dip,nu,R):
    I = (1 - 2*nu) * (-np.log(R + eta)) - I3(eta,q,dip,nu,R)
    return I

# -----------------------------------------------------------------
def I3(eta,q,dip,nu,R):
    yb = eta*np.cos(dip) + q*np.sin(dip)
    db = eta*np.sin(dip) - q*np.cos(dip)
    I = np.where(np.cos(dip) > np.spacing(1),
        (1 - 2*nu) * (yb/(np.cos(dip)*(R + db)) - np.log(R + eta)) + np.sin(dip)/np.cos(dip) * I4(db,eta,q,dip,nu,R),
        (1 - 2*nu)/2 * (eta/(R + db) + yb*q/np.square(R + db) - np.log(R + eta)) )
    return I

# -----------------------------------------------------------------
def I4(db,eta,q,dip,nu,R):
    I = np.where(np.cos(dip) > np.spacing(1),
        (1 - 2*nu) * 1/np.cos(dip) * (np.log(R + db) - np.sin(dip)*np.log(R + eta)),
        -(1 - 2*nu) * q/(R + db) )
    return I

# -----------------------------------------------------------------
def I5(xi,eta,q,dip,nu,R,db):
    X = np.sqrt(np.square(xi) + np.square(q))
    # fix a strange intermittent zero-division warning in the np.where() clause
    xi = np.where(xi==0, xi+np.spacing(1), xi)
    I = np.where(np.cos(dip) > np.spacing(1),
        (1 - 2*nu) * 2/np.cos(dip) * np.arctan((eta*(X + q*np.cos(dip)) + X*(R + X)*np.sin(dip))/(xi*(R + X)*np.cos(dip))),
        -(1 - 2*nu) * xi*np.sin(dip)/(R + db) )
    return I


# =================================================================
# Tilt subfunctions

# strike-slip tilt subfunctions [equation (37) p. 1147]
# -----------------------------------------------------------------
def uzx_ss(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    u = -xi*np.square(q)*A(eta,R)*np.cos(dip) + ((xi*q)/np.power(R,3) - K1(xi,eta,q,dip,nu,R))*np.sin(dip)
    return u

# -----------------------------------------------------------------
def uzy_ss(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    db = eta*np.sin(dip) - q*np.cos(dip)
    yb = eta*np.cos(dip) + q*np.sin(dip)
    u = (db*q/np.power(R,3))*np.cos(dip) + (np.square(xi)*q*A(eta,R)*np.cos(dip) - np.sin(dip)/R + yb*q/np.power(R,3) - K2(xi,eta,q,dip,nu,R))*np.sin(dip)
    return u

# dip-slip tilt subfunctions [equation (38) p. 1147]
# -----------------------------------------------------------------
def uzx_ds(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    db = eta*np.sin(dip) - q*np.cos(dip)
    u = db*q/np.power(R,3) + q*np.sin(dip)/(R*(R + eta)) + K3(xi,eta,q,dip,nu,R)*np.sin(dip)*np.cos(dip)
    return u

# -----------------------------------------------------------------
def uzy_ds(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    db = eta*np.sin(dip) - q*np.cos(dip)
    yb = eta*np.cos(dip) + q*np.sin(dip)
    u = yb*db*q*A(xi,R) - (2*db/(R*(R + xi)) + xi*np.sin(dip)/(R*(R + eta)))*np.sin(dip) + K1(xi,eta,q,dip,nu,R)*np.sin(dip)*np.cos(dip)
    return u

# tensile fault tilt subfunctions [equation (39) p. 1147]
# -----------------------------------------------------------------
def uzx_tf(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    u = np.square(q)/np.power(R,3)*np.sin(dip) - np.power(q,3)*A(eta,R)*np.cos(dip) + K3(xi,eta,q,dip,nu,R)*np.square(np.sin(dip))
    return u

# -----------------------------------------------------------------
def uzy_tf(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    db = eta*np.sin(dip) - q*np.cos(dip)
    yb = eta*np.cos(dip) + q*np.sin(dip)
    u = (yb*np.sin(dip) + db*np.cos(dip))*np.square(q)*A(xi,R) + xi*np.square(q)*A(eta,R)*np.sin(dip)*np.cos(dip) - (2*q/(R*(R + xi)) - K1(xi,eta,q,dip,nu,R))*np.square(np.sin(dip))
    return u

# -----------------------------------------------------------------
def A(x,R):
    a = (2*R + x)/(np.power(R,3)*np.square(R + x))
    return a

# K... tilt subfunctions [equations (40) (41) p. 1148]
# -----------------------------------------------------------------
def K1(xi,eta,q,dip,nu,R):
    db = eta*np.sin(dip) - q*np.cos(dip)
    K = np.where(np.cos(dip) > np.spacing(1),
        (1 - 2*nu) * xi/np.cos(dip) * (1/(R*(R + db)) - np.sin(dip)/(R*(R + eta))),
        (1 - 2*nu) * xi*q/np.square(R + db) )
    return K

# -----------------------------------------------------------------
def K2(xi,eta,q,dip,nu,R):
    K = (1 - 2*nu) * (-np.sin(dip)/R + q*np.cos(dip)/(R*(R + eta))) - K3(xi,eta,q,dip,nu,R)
    return K

# -----------------------------------------------------------------
def K3(xi,eta,q,dip,nu,R):
    db = eta*np.sin(dip) - q*np.cos(dip)
    yb = eta*np.cos(dip) + q*np.sin(dip)
    K = np.where(np.cos(dip) > np.spacing(1),
        (1 - 2*nu) * 1/np.cos(dip) * (q/(R*(R + eta)) - yb/(R*(R + db))),
        (1 - 2*nu) * np.sin(dip)/(R + db) * (np.square(xi)/(R*(R + db)) - 1) )
    return K


# =================================================================
# Strain subfunctions

# strike-slip strain subfunctions [equation (31) p. 1145]
# -----------------------------------------------------------------
def uxx_ss(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    u = np.square(xi)*q*A(eta,R) - J1(xi,eta,q,dip,nu,R)*np.sin(dip)
    return u

# -----------------------------------------------------------------
def uxy_ss(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    db = eta*np.sin(dip) - q*np.cos(dip)
    u = np.power(xi,3)*db/(np.power(R,3)*(np.square(eta) + np.square(q))) - (np.power(xi,3)*A(eta,R) + J2(xi,eta,q,dip,nu,R))*np.sin(dip)
    return u

# -----------------------------------------------------------------
def uyx_ss(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    u = xi*q/np.power(R,3)*np.cos(dip) + (xi*np.square(q)*A(eta,R) - J2(xi,eta,q,dip,nu,R))*np.sin(dip)
    return u

# -----------------------------------------------------------------
def uyy_ss(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    yb = eta*np.cos(dip) + q*np.sin(dip)
    u = yb*q/np.power(R,3)*np.cos(dip) + (np.power(q,3)*A(eta,R)*np.sin(dip) - 2*q*np.sin(dip)/(R*(R + eta)) - (np.square(xi) + np.square(eta))/np.power(R,3)*np.cos(dip) - J4(xi,eta,q,dip,nu,R))*np.sin(dip)
    return u
    
# dip-slip strain subfunctions [equation (32) p. 1146]
# -----------------------------------------------------------------
def uxx_ds(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    u = xi*q/np.power(R,3) + J3(xi,eta,q,dip,nu,R)*np.sin(dip)*np.cos(dip)
    return u

# -----------------------------------------------------------------
def uxy_ds(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    yb = eta*np.cos(dip) + q*np.sin(dip)
    u = yb*q/np.power(R,3) - np.sin(dip)/R + J1(xi,eta,q,dip,nu,R)*np.sin(dip)*np.cos(dip)
    return u

# -----------------------------------------------------------------
def uyx_ds(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    yb = eta*np.cos(dip) + q*np.sin(dip)
    u = yb*q/np.power(R,3) + q*np.cos(dip)/(R*(R + eta)) + J1(xi,eta,q,dip,nu,R)*np.sin(dip)*np.cos(dip)
    return u

# -----------------------------------------------------------------
def uyy_ds(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    yb = eta*np.cos(dip) + q*np.sin(dip)
    u = np.square(yb)*q*A(xi,R) - (2*yb/(R*(R + xi)) + xi*np.cos(dip)/(R*(R + eta)))*np.sin(dip) + J2(xi,eta,q,dip,nu,R)*np.sin(dip)*np.cos(dip)
    return u

# tensile fault strain subfunctions [equation (33) p. 1146]
# -----------------------------------------------------------------
def uxx_tf(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    u = xi*np.square(q)*A(eta,R) + J3(xi,eta,q,dip,nu,R)*np.square(np.sin(dip))
    return u

# -----------------------------------------------------------------
def uxy_tf(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    db = eta*np.sin(dip) - q*np.cos(dip)
    u = -db*q/np.power(R,3) - np.square(xi)*q*A(eta,R)*np.sin(dip) + J1(xi,eta,q,dip,nu,R)*np.square(np.sin(dip))
    return u

# -----------------------------------------------------------------
def uyx_tf(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    u = np.square(q)/np.power(R,3)*np.cos(dip) + np.power(q,3)*A(eta,R)*np.sin(dip) + J1(xi,eta,q,dip,nu,R)*np.square(np.sin(dip))
    return u

# -----------------------------------------------------------------
def uyy_tf(xi,eta,q,dip,nu):
    R = np.sqrt(np.square(xi) + np.square(eta) + np.square(q))
    db = eta*np.sin(dip) - q*np.cos(dip)
    yb = eta*np.cos(dip) + q*np.sin(dip)
    u = (yb*np.cos(dip) - db*np.sin(dip))*np.square(q)*A(xi,R) - q*np.sin(2*dip)/(R*(R + xi)) - (xi*np.square(q)*A(eta,R) - J2(xi,eta,q,dip,nu,R))*np.square(np.sin(dip))
    return u


# J... tensile fault subfunctions [equations (34) (35) p. 1146-1147]
# -----------------------------------------------------------------
def J1(xi,eta,q,dip,nu,R):
    db = eta*np.sin(dip) - q*np.cos(dip)
    J = np.where(np.cos(dip) > np.spacing(1),
        (1 - 2*nu) * 1/np.cos(dip) * (np.square(xi)/(R*np.square(R + db)) - 1/(R + db)) - np.sin(dip)/np.cos(dip)*K3(xi,eta,q,dip,nu,R),
        (1 - 2*nu)/2 * q/np.square(R + db) * (2*np.square(xi)/(R*(R + db)) - 1) )
    return J

# -----------------------------------------------------------------
def J2(xi,eta,q,dip,nu,R):
    db = eta*np.sin(dip) - q*np.cos(dip)
    yb = eta*np.cos(dip) + q*np.sin(dip)
    J = np.where(np.cos(dip) > np.spacing(1),
        (1 - 2*nu) * 1/np.cos(dip) * xi*yb/(R*np.square(R + db)) - np.sin(dip)/np.cos(dip)*K1(xi,eta,q,dip,nu,R),
        (1 - 2*nu)/2 * xi*np.sin(dip)/np.square(R + db) * (2*np.square(q)/(R*(R + db)) - 1) )
    return J

# -----------------------------------------------------------------
def J3(xi,eta,q,dip,nu,R):
    J = (1 - 2*nu) * -xi/(R*(R + eta)) - J2(xi,eta,q,dip,nu,R)
    return J

# -----------------------------------------------------------------
def J4(xi,eta,q,dip,nu,R):
    J = (1 - 2*nu) * (-np.cos(dip)/R - q*np.sin(dip)/(R*(R + eta))) - J1(xi,eta,q,dip,nu,R)
    return J

