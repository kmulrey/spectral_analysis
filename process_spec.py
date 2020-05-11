import numpy as np
import os
import pickle
import matplotlib.pyplot as plt
import sys
import glob
import process_files
import re
import matplotlib as mpl
import matplotlib.cm as cm
from scipy.optimize import curve_fit
import fit_functions as fit


# process for characterizing LOFAR bandwidth 30-80 MHz
def fit_spec_LOFAR(freq,spec):

    fmin=30
    fmax=80
    fmax2=80
    e_min=-3.3
    min_points=8
    flag=0

    n_points=len(freq[(freq>fmin)*(freq<fmax)*(spec>e_min)])

    if n_points<min_points:
        flag=1

    try:
        a,m,n,chi2,redchi2=fit.fit_routine_second_order(freq[(freq>fmin)*(freq<fmax2)*(spec>e_min)],spec[(freq>fmin)*(freq<fmax2)*(spec>e_min)])
    except:
        print("error doin LOFAR spec fit")
        flag=1
        a=0
        m=0
        n=0
        chi2=100
        redchi2=100
        n_points=0
        

    return a,m,n,chi2,redchi2,n_points,flag


def fit_slope_LOFAR(freq,spec):
    
    fmin=30
    fmax=80
    fmax2=80
    e_min=-3.3
    min_points=8
    flag=0

    n_points=len(freq[(freq>fmin)*(freq<fmax)*(spec>e_min)])

    if n_points<min_points:
        flag=1

    try:
        a,m,chi2,redchi2=fit.fit_routine_linear(freq[(freq>fmin)*(freq<fmax2)*(spec>e_min)],spec[(freq>fmin)*(freq<fmax2)*(spec>e_min)])
    except:
        flag=1

    return a,m,chi2,redchi2,n_points,flag



# process for characterizing SKA bandwidth 50-350 MHz
def fit_spec_SKA(freq,spec):


    i=14
    fmin=50
    fmax=350
    fmax2=250
    min_points=8

    f_0=50
    f_1=110

    e_min=-3.3
    logspec1=np.log10(np.abs(spec))
    n_points=len(freq[(freq>fmin)*(freq<fmax)*(spec>e_min)])

    if n_points<min_points:
        flag=1

    flag=0

    if len(freq[(freq>fmin)*(freq<fmax)*(logspec1>e_min)])<10:
        flag=1


    try:
        a_init,m_init,n_init,chi2_init,redchi2_init=fit.fit_routine_second_order(freq[(freq>f_0)*(freq<f_1)*(logspec1>e_min)], logspec1[(freq>f_0)*(freq<f_1)*(logspec1>e_min)])
    except:
        a_init=0
        m_init=0
        n_init=0
        chi2_init=0
        redchi2_init=0

    regime=-1
    if n_init<-2.5e-5:
        regime=0 # close to core
    elif n_init>=-2.5e-5 and n_init<0:
        regime=1  # double 2nd order poly
    else:
        regime=2 # 2nd order poly + flat


    a_use=0
    m_use=0
    n_use=0
    p_use=0
    q_use=0
    t_use=0
    chi2_use=100
    redchi2_use=100


    if regime==0:
        a_use,m_use,n_use,t_use,b_use,p_use,q_use,chi2_use,redchi2_use=fit.fit_routine_piecewise2(freq[(freq>fmin)*(freq<fmax2)*(logspec1>e_min)], logspec1[(freq>fmin)*(freq<fmax2)*(logspec1>e_min)])
    elif regime==1:
        a_use,m_use,n_use,t_use,b_use,p_use,q_use,chi2_use,redchi2_use=fit.fit_routine_piecewise2(freq[(freq>fmin)*(freq<fmax)*(logspec1>e_min)], logspec1[(freq>fmin)*(freq<fmax)*(logspec1>e_min)])
    elif regime==2:
        a_use,m_use,n_use,t_use,chi2_use,redchi2_use=fit.fit_routine_piecewise1(freq[(freq>fmin)*(freq<fmax)*(logspec1>e_min)], logspec1[(freq>fmin)*(freq<fmax)*(logspec1>e_min)])
    else:
        flag=1


    return a_use,m_use,n_use,chi2_use,redchi2_use,n_points,flag





# process for characterizing CODALEMA bandwidth 30-250 MHz
def fit_spec_CODALEMA(freq,spec):


    i=14
    fmin=30
    fmax=250
    fmax2=250
    min_points=8

    f_0=30
    f_1=110


    e_min=-3.3
    logspec1=np.log10(np.abs(spec))
    n_points=len(freq[(freq>fmin)*(freq<fmax)*(spec>e_min)])

    if n_points<min_points:
        flag=1

    flag=0

    if len(freq[(freq>fmin)*(freq<fmax)*(logspec1>e_min)])<10:
        flag=1


    try:
        a_init,m_init,n_init,chi2_init,redchi2_init=fit.fit_routine_second_order(freq[(freq>f_0)*(freq<f_1)*(logspec1>e_min)], logspec1[(freq>f_0)*(freq<f_1)*(logspec1>e_min)])
    except:
        a_init=0
        m_init=0
        n_init=0
        chi2_init=0
        redchi2_init=0

    regime=-1
    if n_init<-2.5e-5:
        regime=0 # close to core
    elif n_init>=-2.5e-5 and n_init<0:
        regime=1  # double 2nd order poly
    else:
        regime=2 # 2nd order poly + flat


    a_use=0
    m_use=0
    n_use=0
    p_use=0
    q_use=0
    t_use=0
    chi2_use=100
    redchi2_use=100


    if regime==0:
        a_use,m_use,n_use,t_use,b_use,p_use,q_use,chi2_use,redchi2_use=fit.fit_routine_piecewise2(freq[(freq>fmin)*(freq<fmax2)*(logspec1>e_min)], logspec1[(freq>fmin)*(freq<fmax2)*(logspec1>e_min)])
    elif regime==1:
        a_use,m_use,n_use,t_use,b_use,p_use,q_use,chi2_use,redchi2_use=fit.fit_routine_piecewise2(freq[(freq>fmin)*(freq<fmax)*(logspec1>e_min)], logspec1[(freq>fmin)*(freq<fmax)*(logspec1>e_min)])
    elif regime==2:
        a_use,m_use,n_use,t_use,chi2_use,redchi2_use=fit.fit_routine_piecewise1(freq[(freq>fmin)*(freq<fmax)*(logspec1>e_min)], logspec1[(freq>fmin)*(freq<fmax)*(logspec1>e_min)])
    else:
        flag=1


    return a_use,m_use,a_use,m_use,n_use,t_use,chi2_use,redchi2_use,n_points,flag

