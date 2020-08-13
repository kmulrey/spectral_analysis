import numpy as np
#import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#from radiotools.atmosphere import models as atm

# see Glaser et al., JCAP 09(2016)024 for the derivation of the formulas

#average_xmax = 669.40191244545326  # 1 EeV, 50% proton, 50% iron composition
#average_zenith = np.deg2rad(45)
#atmc = atm.Atmosphere(model=1)
average_density = 6.5e-4 #in g/cm^3#atmc.get_density(average_zenith, average_xmax) * 1e-3  # in kg/m^3
mag=2.03


rho_avg=0.72  # average density, 0.72 kg/m3
# charge excess parameters
p0a=0.2
p1a=1.27
p2a=-0.08

# slant depth parameters proton
p0s_p=1.02
p1s_p=-0.49
p2s_p=0

# slant depth parameters iron
p0s_fe=0.98
p1s_fe=-0.47
p2s_fe=0


def integral(fluence,pos):

    pos_uvw_vxb=pos[0::8]
    pos_uvw_vxvxb=pos[2::8]
    neg_uvw_vxb=pos[4::8]
    neg_uvw_vxvxb=pos[6::8]
    
    fluence_pos_vxb=fluence[0::8]
    fluence_pos_vxvxb=fluence[2::8]
    fluence_neg_vxb=fluence[4::8]
    fluence_neg_vxvxb=fluence[6::8]
        
        
    
    pos_vxvxb_all=np.concatenate([neg_uvw_vxvxb.T[1],pos_uvw_vxvxb.T[1]])
    fluence_vxvxb=np.concatenate([fluence_pos_vxvxb,fluence_neg_vxvxb])
    inds = pos_vxvxb_all.argsort()
        
    sorted_pos=pos_vxvxb_all[inds]
    sorted_fluence_vxvxb=fluence_vxvxb[inds]
    f0 = interp1d(sorted_pos, sorted_fluence_vxvxb, kind='cubic')
    r= np.linspace(0, 500, num=1000, endpoint=True)
    
    
    # integrate positive vxvxB arm
    n=len(r)
    dr=r[1]-r[0]
    integral=0
    
    for i in np.arange(n-1):
        r0=r[i]
        r1=r[i+1]
        val0=r0*(f0(r0))
        val1=r1*(f0(r1))
        integral=integral+(val0+val1)*0.5*dr

    integral=integral*2*np.pi


    return integral


def StoEm(S,A=1.683,B=2.006):
    
    Em=np.power((S/(A*1e7)),1/B)*1e18

    return Em



def return_Srd(Erad,zenith,density,aplpha,type):
    a=return_a(density[i]*1e3*np.cos(zenith[i]),rho_avg,p0a,p1a,p2a)/mag**0.9
    Srd=Erad[i]/np.sin(alpha[i])**2/mag**1.8
    Srd_1=Erad[i]/(a**2+(1-a**2)*np.sin(alpha[i])**2)/mag**1.8
    if type[i]==0:
        Srd_2=Erad/(a**2+(1-a**2)*np.sin(alpha)**2)/(1-p0s_p+p0s_p*np.exp(p1s_p*(density*1e3*np.cos(zenith)-rho_avg)))**2/mag**1.8
    else:
        Srd_2=Erad/(a**2+(1-a**2)*np.sin(alpha)**2)/(1-p0s_fe+p0s_fe*np.exp(p1s_fe*(density*1e3*np.cos(zenith)-rho_avg)))**2/mag**1.8

    return Srd_2


def return_a(rho,avg,p0,p1,p2):
    a= p2+p0*np.exp(p1*(rho-avg))
    return a



