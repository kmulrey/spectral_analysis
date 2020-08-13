import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
import sys
import glob
import process_data
from scipy import signal
import re
import matplotlib as mpl
import matplotlib.cm as cm
from scipy.optimize import curve_fit
from lmfit import Model
import fit_functions as fit
import process_spec as spec
import atmosphere as atmosphere
from scipy.interpolate import interp1d
import fluence as fluence
import radiation_energy as radiation_energy
import importlib
import scipy.interpolate as intp
import process_files




def return_Srd_info(RUNNR,event_no,primary):


    data_dir='/vol/astro7/lofar/sim/pipeline/events/'+event_no+'/1/coreas/'+primary+'/'
    steering_dir='/vol/astro7/lofar/sim/pipeline/events/'+event_no+'/1/coreas/'+primary+'/steering/'
    coreas_dir='/vol/astro7/lofar/sim/pipeline/events/'+event_no+'/1/coreas/'+primary+'/SIM'+RUNNR+'_coreas/'
    atm_dir='/vol/astro7/lofar/sim/pipeline/atmosphere_files/'
    
    coreas_files=(glob.glob(coreas_dir+'*.dat'))
    inp_file=steering_dir+'RUN'+RUNNR+'.inp'
    list_file=steering_dir+'SIM'+RUNNR+'.list'
    nFiles=len(coreas_files)
    antenna_positions,ant_pos_uvw,time,efield,poldata,zenith,az_rot,energy,xmax=process_files.get_efield(data_dir,RUNNR)
    alpha= process_files.GetAlpha(zenith,az_rot,1.1837)
    em_dep,other_dep,total_dep=process_files.getEM(data_dir,RUNNR)




    atm=process_files.get_atm(event_no)
    hi=atmosphere.get_vertical_height(xmax,atm)   # height in cm
    at=atmosphere.get_atmosphere(hi,atm)
    rho=atmosphere.return_density(hi, atm)
    rho2=atmosphere.return_density(hi/np.cos(zenith), atm)
    dmax=atmosphere.get_distance_xmax_geometric(zenith, xmax, atm)
    n_xmax=atmosphere.get_n_at_xmax(atm_dir+'ATMOSPHERE_'+event_no+'.DAT',hi)
    cherenkov_angle=np.rad2deg(np.arccos(1/n_xmax))
    alpha=process_files.GetAlpha(zenith,az_rot,1.1837)
    
    e_filt_30_80,time_filt_30_80=process_data.lofar_filter(efield,time,30.0,80.0,1.0)
    e_filt_30_200,time_filt_30_200=process_data.lofar_filter(efield,time,30.0,200.0,1.0)
    e_filt_50_350,time_filt_50_350=process_data.lofar_filter(efield,time,50.0,350.0,1.0)

    fluence_30_80=fluence.calculate_energy_fluence_vector(efield_30_80, time, signal_window=100., remove_noise=True)
    fluence_30_200=fluence.calculate_energy_fluence_vector(efield_30_200, time, signal_window=100., remove_noise=True)
    fluence_50_350=fluence.calculate_energy_fluence_vector(efield_50_350, time, signal_window=100., remove_noise=True)

    pos_uvw_vxb=ant_pos_uvw[0::8]
    pos_uvw_vxvxb=ant_pos_uvw[2::8]
    neg_uvw_vxvxb=ant_pos_uvw[6::8]

    fluence_30_80_vxvxb_0=np.concatenate([fluence_30_80[2::8].T[0],fluence_30_80[6::8].T[0]])
    fluence_30_80_vxvxb_1=np.concatenate([fluence_30_80[2::8].T[1],fluence_30_80[6::8].T[1]])

    pos_vxvxb_all=np.concatenate([neg_uvw_vxvxb.T[1],pos_uvw_vxvxb.T[1]])
            
    inds = pos_vxvxb_all.argsort()
            
    sorted_pos=pos_vxvxb_all[inds]
    sorted_fluence_30_80_gm=fluence_30_80_vxvxb_0[inds]
    sorted_fluence_30_80_ce=fluence_30_80_vxvxb_1[inds]
    xnew = np.linspace(0, 400, num=1000, endpoint=True)
    f0 = interp1d(sorted_pos, sorted_fluence_30_80_gm, kind='cubic')
    f1 = interp1d(sorted_pos, sorted_fluence_30_80_ce, kind='cubic')

    Erad=radiation_energy.integrate(xnew,f0(xnew),f1(xnew))
    Erad_gm=radiation_energy.integrate_one_pol(xnew,f0(xnew))
    Erad_ce=radiation_energy.integrate_one_pol(xnew,f1(xnew))
    clip_ratio=radiation_energy.get_clipping(dmax)



    return energy,zenith,az_rot,xmax,hi,rho,rho2,dmax,n_xmax,alpha,clip_ratio,Erad,Erad_gm,Erad_ce,em_dep,total_dep
