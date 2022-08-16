
import numpy as np
import fault_model

# define our fault modeling function for the sampler to fit.
def fault_model_for_fitting(params,gps_locs):
    # expand the '[params]' list into its elements
    latc,lonc,depthc,strike,dip,L,W,ss,ds = params
    # for the modeling, we need the locations as x and y, not a single vector
    Ngps=int(np.size(gps_locs)/2)
    gpslon=gps_locs[:Ngps]
    gpslat=gps_locs[Ngps:]
    # create the model
    Fmod=fault_model.FaultModel()
    Fmod.create_planar_model_centered(latc=latc,lonc=lonc,depthc=depthc,strike=strike,dip=dip,L=L,W=W,nL=1,nW=1)
    # get the "G" matrix for our particular GPS site locations (lat,lon)
    G=Fmod.get_greens(gpslat,gpslon)
    # m is the slip
    m = np.array([ss,ds])
    # predict GPS displacements with d=G*m
    predicted_displacements = np.matmul(G,m)                                  
    # return the predicted values
    return predicted_displacements

# "log likelihood" function for the sampler
def lnlike_fault(params, x, y, yerr):
    # get the predicted model values
    ypred = fault_model_for_fitting(params,x)
    # compute the misfit - sum of the residuals squared, scaled by the data uncertainties
    misfit = -0.5*np.sum(((y-ypred)/yerr)**2)
    return misfit

# set our priors - we use this to set bounds on the parameters.
# return value is set to negative infinity if any parameters are outside their bounds, otherwise it is zero. 
# this is because we have taken the log() of our probability distribution. So 10^0 = 1, while 10^-np.inf = 0.
# if we wanted gaussian priors, or other types, we could also implement them here instead of bounds.
def lnprior_fault(params):
    # define bounds here
    # order of parameters is: latc,lonc,depthc,strike,dip,L,W,ss,ds
    # from a scipy.optimize fit:
    #Best fitting solution:
    # lat 35.730677, lon -117.557194, depth 4362.581680, 
    # strike -40.846756, dip 67.706525, length 33048.330581, width 10759.567750, 
    # ss -400.061320, ds -52.384303
    # note that depth has a special variable minimum of half the fault width times the sine of the dip, 
    # since it is required that the fault does not cut the surface.
    mindepth = 0.5*params[6]*np.sin(params[4]*np.pi/180.)
    minvals = np.array([35.2,-118,mindepth,-60,45,20e3,5e3,-500,-200])
    maxvals = np.array([36.2,-117,20e3,-20,135,200e3,40e3,500,200])
    if any(params-minvals<0) or any(params-maxvals>0): # check if any bounds are exceeded
        return -np.inf
    else:
        return 0.0
    
#finally, this function puts together all the above, to determine the actual log(probability) of a set of parameters.
def lnprob_fault(params, x, y, yerr):
    prior = lnprior_fault(params)
    if np.isinf(prior):
        return -np.inf
    else:
        return prior + lnlike_fault(params, x, y, yerr) #recall if lp not -inf, its 0, so this just returns likelihood