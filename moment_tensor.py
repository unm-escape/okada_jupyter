# Eric Lindsey, 04/2016
import numpy as np

def get_moment_tensor(strike,dip,rake,moment):
    #   % convert strike-dip-rake to HRV mt
    # f = strike, d = dip, l = rake
    #1 Mrr =  Mzz =  Mo sin2d sinl
    #2 Mtt =  Mxx = -Mo(sind cosl sin2f +     sin2d sinl (sinf)^2 )
    #3 Mpp =  Myy =  Mo(sind cosl sin2f -     sin2d sinl (cosf)^2 )
    #4 Mrt =  Mxz = -Mo(cosd cosl cosf  +     cos2d sinl sinf )
    #5 Mrp = -Myz =  Mo(cosd cosl sinf  -     cos2d sinl cosf )
    #6 Mtp = -Mxy = -Mo(sind cosl cos2f + 0.5 sin2d sinl sin2f )

    sin2d=np.sin(np.radians(2*dip))
    sinl=np.sin(np.radians(rake))
    sind=np.sin(np.radians(dip))
    cosd=np.cos(np.radians(dip))
    cos2d=np.cos(np.radians(2*dip))
    cosl=np.cos(np.radians(rake))
    sin2f=np.sin(np.radians(2*strike))
    cosf=np.cos(np.radians(strike))
    sinf=np.sin(np.radians(strike))
    cos2f=np.cos(np.radians(2*strike))
    
    Mrr =  moment*sin2d*sinl
    Mtt = -moment*(sind*cosl*sin2f +     sin2d*sinl*sinf**2 )
    Mpp =  moment*(sind*cosl*sin2f -     sin2d*sinl*cosf**2 )
    Mrt = -moment*(cosd*cosl*cosf  +     cos2d*sinl*sinf    )
    Mrp =  moment*(cosd*cosl*sinf  -     cos2d*sinl*cosf    )
    Mtp = -moment*(sind*cosl*cos2f + 0.5*sin2d*sinl*sin2f   )
    
    return Mrr,Mtt,Mpp,Mrt,Mrp,Mtp   
    
def get_moment(L,W,slip):
    moment=1e7*30e9*L*W*slip
    return moment    
    
def get_magnitude(L,W,slip):
    moment=get_moment(L,W,slip)
    magnitude=(2/3)*np.log10(moment) - 10.7
    return magnitude
    
def print_moment_tensor(strike,dip,rake,moment):
    Mrr,Mtt,Mpp,Mrt,Mrp,Mtp = get_moment_tensor(strike,dip,rake,moment)
    print("Mrr = ",Mrr)
    print("Mtt = ",Mtt)
    print("Mpp = ",Mpp)
    print("Mrt = ",Mrt)
    print("Mrp = ",Mrp)
    print("Mtp = ",Mtp)
    
def save_moment_tensor(fname,strike,dip,rake,L,W,slip,lonc,latc,depth,year,month,day,hour,minute,sec,duration,eventid):
    moment=get_moment(L,W,slip)
    magnitude=get_magnitude(L,W,slip)
    Mrr,Mtt,Mpp,Mrt,Mrp,Mtp = get_moment_tensor(strike,dip,rake,moment)
    f1=open(fname, 'w')
    f1.write("%d %02d %02d %02d %02d %05.2f %.4f %.4f %.4f %.1f %.1f %s\n"%(year,month,day,hour,minute,sec,latc,lonc,depth,magnitude,magnitude,eventid))
    f1.write("event name:    %s\n"%eventid)
    f1.write("time shift:    0.0000\n")
    f1.write("half duration: %.4f\n"%(duration/2.))
    f1.write("latitude:      %.4f\n"%latc)
    f1.write("longitude:     %.4f\n"%lonc)
    f1.write("depth:         %.4f\n"%depth)
    f1.write("Mrr:           %.6e\n"%Mrr)
    f1.write("Mtt:           %.6e\n"%Mtt)
    f1.write("Mpp:           %.6e\n"%Mpp)
    f1.write("Mrt:           %.6e\n"%Mrt)
    f1.write("Mrp:           %.6e\n"%Mrp)
    f1.write("Mtp:           %.6e\n"%Mtp)
    f1.write("\n")
    f1.close()
