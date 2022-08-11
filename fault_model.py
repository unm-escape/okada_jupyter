#!/usr/local/bin/python
# Eric Lindsey, 01/2014

import numpy as np
import geod_transform
import okada_greens
import sys

class FaultModel:
    
    def __init__(self):
        self.latc = np.array([])
        self.lonc = np.array([])
        self.depth= np.array([])
        self.strike = np.array([])
        self.dip = np.array([])
        self.L = np.array([])
        self.W = np.array([])        
        self.strikeid = np.array([])
        self.dipid = np.array([])
        self.npatches = 0
        #future options - rake, slip? or add this to "slip_model (extends fault_model)"

    def add_patch(self,latc,lonc,depth,strike,dip,L,W,strikeid,dipid):
        #given parameters for a single fault patch (or several patches), add them to the overall fault model
        self.latc = np.append(self.latc,latc)
        self.lonc = np.append(self.lonc,lonc)
        self.depth = np.append(self.depth,depth)
        self.strike = np.append(self.strike,strike)
        self.dip = np.append(self.dip,dip)
        self.L = np.append(self.L,L)
        self.W = np.append(self.W,W)      
        self.strikeid = np.append(self.strikeid,strikeid)
        self.dipid = np.append(self.dipid,dipid)
        self.npatches = len(self.W)
            
    def load_patches_topleft(self,fname):
        # supply patches with one corner and the width/length
        # load file containing sub-dividable fault segments for a fault model:
        # columns are [lat0 lon0 lat1 lon1 depth width dip nL nW] [optional: patchstrike, patchdip]
        filedata=np.loadtxt(fname,ndmin=2)
        
        dipid=filedata[:,1]
        strikeid=filedata[:,2]
        lon=filedata[:,3]
        lat=filedata[:,4]
        depth=filedata[:,5]
        L=filedata[:,6]
        W=filedata[:,7]
        strike=filedata[:,8]
        dip=filedata[:,9]
        
        # transform into patch-centered coordinates
        # compute East and North offsets from "top left" corner to patch center
        eoffset=(L/2.)*np.sin(np.radians(strike))+(W/2.)*np.cos(np.radians(dip))*np.cos(np.radians(strike))
        noffset=(L/2.)*np.cos(np.radians(strike))+(W/2.)*np.cos(np.radians(dip))*np.sin(np.radians(strike))
        uoffset=(W/2.)*np.sin(np.radians(dip))
        latc=[]
        lonc=[]
        depthc=[]
        for i in range(len(lat)):
            latci,lonci,depthci = geod_transform.translate_flat(lat[i],lon[i],depth[i],eoffset[i],noffset[i],uoffset[i])
            latc.append(latci)
            lonc.append(lonci)
            depthc.append(depthci)
        
        # create fault model. first look for this object stored as a pickle, 
        # using hash name determined by the contents of filedata
        #filedatahash = hashlib.sha1(filedata.view(np.uint8)).hexdigest()
        #if not os.path.isdir('hash'):
        #        os.mkdir('hash')
        #picklename="hash/load_patches_topleft_" + filedatahash + ".npy"
        #if os.path.isfile(picklename):
        #    print("loading saved copy of fault object")
        #    self.load_pickle(picklename)
        #else:
        #    print("saved fault object not found, creating new fault object")
        self.add_patch(latc,lonc,depthc,strike,dip,L,W,strikeid,dipid)           
        #self.save_pickle(picklename)
        
    def load_patches_center(self,fname):
        # supply patches with one corner and the width/length
        # load file containing sub-dividable fault segments for a fault model:
        # columns are [lat0 lon0 lat1 lon1 depth width dip nL nW] [optional: patchstrike, patchdip]
        filedata=np.loadtxt(fname,ndmin=2)
        
        dipid=filedata[:,1]
        strikeid=filedata[:,2]
        lonc=filedata[:,3]
        latc=filedata[:,4]
        depth=filedata[:,5]
        L=filedata[:,6]
        W=filedata[:,7]
        strike=filedata[:,8]
        dip=filedata[:,9]

        self.add_patch(latc,lonc,depth,strike,dip,L,W,strikeid,dipid)  

    def load_patches_comsol(self,fname):
        # supply patches with one corner and the width/length
        # load file containing sub-dividable fault segments for a fault model:
        # columns are [lat0 lon0 lat1 lon1 depth width dip nL nW] [optional: patchstrike, patchdip]
        filedata=np.loadtxt(fname,ndmin=2)
       
        lonc=filedata[:,1]
        latc=filedata[:,2]
        depth=filedata[:,5]
        L=filedata[:,3]
        W=filedata[:,4]
        strike=filedata[:,6]
        dip=filedata[:,7]
        dipid=filedata[:,8]
        strikeid=filedata[:,9]

        self.add_patch(latc,lonc,depth,strike,dip,L,W,strikeid,dipid)       
        
    def create_planar_model(self,latcorner,loncorner,depthcorner,strike,dip,L,W,nL,nW):
        #compute coordinates of every patch center and add this patch
        patchL=L/nL
        patchW=W/nW  
        for i in range(nL):
            for j in range(nW):
                tot_eoffset=(i+0.5)*patchL*np.sin(np.radians(strike))+(j+0.5)*patchW*np.cos(np.radians(dip))*np.cos(np.radians(strike))
                tot_noffset=(i+0.5)*patchL*np.cos(np.radians(strike))-(j+0.5)*patchW*np.cos(np.radians(dip))*np.sin(np.radians(strike))
                tot_uoffset=(j+0.5)*patchW*np.sin(np.radians(dip))                
                
                #get next patch-centered coordinate
                latcij,loncij,depthcij = geod_transform.translate_flat(latcorner,loncorner,depthcorner,tot_eoffset,tot_noffset,tot_uoffset)
                self.add_patch(latcij,loncij,depthcij,strike,dip,patchL,patchW,i,j)
                
    def find_patch(self,patchid,strikeoffset,dipoffset):
        '''given a patch index, find the patch offset by strikeoffset,dipoffset and
        return its index, e.g. (n,0,1) returns the next patch down-dip from n. 
        Attempting to access indices beyond the edge of the fault returns the edge.'''
        strid=min(max(self.strikeid),max(0,self.strikeid[patchid]+strikeoffset))
        dipid=min(max(self.dipid),max(0,self.dipid[patchid]+dipoffset))
        try:
            newid=np.logical_and(self.strikeid==strid, self.dipid==dipid).nonzero()[0][0]
            return newid
        except IndexError:
            print('Error in find_patch: element not found')
            
    def get_greens(self,lat,lon,kind='displacement'):
        # return the greens function matrix at specified input points
        #layout of displacement G is:
        #      ...patch1....
        #    . uE_str uE_dip
        # pt1. uN_str uN_dip
        #    . uU_str uU_dip
        #layout of strain G is:
        #      ....patch1.....
        #    . eNN_str eNN_dip
        # pt1. eNE_str eNE_dip
        #    . eEN_str eEN_dip
        #    . eEE_str eEE_dip
        if kind=='displacement': #default value
            G = okada_greens.displacement_greens(lat,lon,self.latc,self.lonc,self.depth,self.strike,self.dip,self.L,self.W)
        elif kind=='strain':
            G = okada_greens.strain_greens(lat,lon,self.latc,self.lonc,self.depth,self.strike,self.dip,self.L,self.W)    
        else:
            print("didn't understand 'kind' argument for get_greens:", kind, " exiting.")
            sys.exit(1)
        return G

    def load_pickle(self,fname):
        #get pickled arrays
        savearray=np.load(fname)
        #add them to self
        self.latc = np.append(self.latc,savearray[:,0])
        self.lonc = np.append(self.lonc,savearray[:,1])
        self.depth = np.append(self.depth,savearray[:,2])
        self.strike = np.append(self.strike,savearray[:,3])
        self.dip = np.append(self.dip,savearray[:,4])
        self.L = np.append(self.L,savearray[:,5])
        self.W = np.append(self.W,savearray[:,6])        
        self.strikeid = np.append(self.strikeid,savearray[:,7])
        self.dipid = np.append(self.dipid,savearray[:,8])
        self.npatches = len(self.W)

    def save_pickle(self,fname):
        print(np.shape(self.latc),np.shape(self.lonc),np.shape(self.depth),np.shape(self.strike),np.shape(self.dip),np.shape(self.L),np.shape(self.W),np.shape(self.strikeid),np.shape(self.dipid))
        savearray=np.column_stack((self.latc,self.lonc,self.depth,self.strike,self.dip,self.L,self.W,self.strikeid,self.dipid))
        np.save(fname,savearray)

    #def get_Ghashid(self,lat,lon,kind='displacement'):
    #    #using current object data compute the hash id 
    #    temparray=np.hstack((lat,lon,np.fromstring(kind,dtype=np.uint8), self.latc,self.lonc,self.depth,self.strike,self.dip,self.L,self.W,self.patchid,self.strikeid,self.dipid))
    #    return hashlib.sha1(temparray.view(np.uint8)).hexdigest()

    def get_patch_verts_center_3d(self):
        verts3d,verts2d=self.get_patch_verts_center_both()
        return verts3d
        
    def get_patch_verts_center_2d(self):
        verts3d,verts2d=self.get_patch_verts_center_both()
        return verts2d
                
    def get_patch_verts_center_both(self):
        verts3d=[]
        verts2d=[]
        print(self.latc)
        for i in range(len(self.latc)):
            sindip=np.sin(np.radians(self.dip[i]))
            cosdip=np.cos(np.radians(self.dip[i]))
            sinstr=np.sin(np.radians(self.strike[i]))
            cosstr=np.cos(np.radians(self.strike[i]))
            ztop=1.e-3*(self.depth[i]+(self.W[i]/2.)*sindip)
            zbot=1.e-3*(self.depth[i]-(self.W[i]/2.)*sindip)
            zs=[ztop,ztop,zbot,zbot]
            y1,x1,z1=geod_transform.translate_flat(self.latc[i],self.lonc[i],0,
                     -(self.L[i]/2.)*sinstr + (self.W[i]/2.)*cosdip*cosstr,
                     -(self.L[i]/2.)*cosstr - (self.W[i]/2.)*cosdip*sinstr,0)
            y2,x2,z2=geod_transform.translate_flat(self.latc[i],self.lonc[i],0,
                     +(self.L[i]/2.)*sinstr + (self.W[i]/2.)*cosdip*cosstr,
                     +(self.L[i]/2.)*cosstr - (self.W[i]/2.)*cosdip*sinstr,0)
            y3,x3,z3=geod_transform.translate_flat(self.latc[i],self.lonc[i],0,
                     +(self.L[i]/2.)*sinstr - (self.W[i]/2.)*cosdip*cosstr,
                     +(self.L[i]/2.)*cosstr + (self.W[i]/2.)*cosdip*sinstr,0)
            y4,x4,z4=geod_transform.translate_flat(self.latc[i],self.lonc[i],0,
                     -(self.L[i]/2.)*sinstr - (self.W[i]/2.)*cosdip*cosstr,
                     -(self.L[i]/2.)*cosstr + (self.W[i]/2.)*cosdip*sinstr,0)
            xs=[x1,x2,x3,x4]
            ys=[y1,y2,y3,y4]
            verts3d.append(list(zip(xs, ys, zs)))
            verts2d.append(list(zip(xs, ys)))
        return verts3d,verts2d