#collation of fit functions and fit routines to handle characterizing specta

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
from lmfit import Model

# first order polynomial fit

def linear(x,b,m):
    y=np.zeros(len(x))
    for i in np.arange(len(x)):
        y[i]=b+x[i]*m
    return y



def first_order(x,a,m):
    y=np.zeros(len(x))
    for i in np.arange(len(x)):
        y[i]=a+x[i]*m
    return y

# second order polynomial fit
def second_order(x,a,m,n):
    y=np.zeros(len(x))
    for i in np.arange(len(x)):
        y[i]=a+x[i]*m+(x[i]-80)**2*n
    return y


# piecewise: second order polynomial with flat line, transition at t
def piecewise1(x,a,m,n,t):
    y=np.zeros(len(x))
    h=0
    for i in np.arange(len(x)):
        if x[i]<t:
            y[i]=a+x[i]*m+(x[i]-80)**2*n
            h=y[i]
        else:
            y[i]=h
    return y

# piecewise: second order polynomial  x 2, transition at t
def piecewise2(x,a,m,n,t,b,p,q):
    y=np.zeros(len(x))
    for i in np.arange(len(x)):
        if x[i]<t:
            y[i]=a+x[i]*m+(x[i]-80)**2*n
        else:
            y[i]=b+x[i]*p+(x[i]-80)**2*q
    return y



def fit_routine_first_order(freq,spec):
    mymodel = Model(first_order)
    params = mymodel.make_params(a=-2, m=-1.2e-2)
    result = mymodel.fit(spec, params, x=freq)
    chi2=result.chisqr
    red_chi2=result.redchi
    a=result.best_values['a']
    m=result.best_values['m']
    return a,m,chi2,red_chi2

def fit_routine_second_order(freq,spec):
    mymodel = Model(second_order)
    params = mymodel.make_params(a=-2, m=-1.2e-2, n=3e-5)
    result = mymodel.fit(spec, params, x=freq)
    chi2=result.chisqr
    red_chi2=result.redchi
    a=result.best_values['a']
    m=result.best_values['m']
    n=result.best_values['n']
    return a,m,n,chi2,red_chi2




def fit_routine_piecewise1(freq,spec):
    chi2s=100*np.zeros([320])
    turnover=np.zeros([320])
    mymodel = Model(piecewise1)

    for i in np.arange(320):
        turnover[i]=i+30
        params = mymodel.make_params(a=-2, m=-1.2e-2, n=3e-5,t=turnover[i])
        params['t'].vary = False
        result = mymodel.fit(spec, params, x=freq)
        chi2s[i]=result.chisqr
        
    t0=turnover[np.argmin(chi2s)]
    params = mymodel.make_params(a=-2, m=-1.2e-2, n=3e-5,t=t0)
    params['t'].vary = False
    result = mymodel.fit(spec, params, x=freq)

    a0=result.best_values['a']
    m0=result.best_values['m']
    n0=result.best_values['n']
    chi2=result.chisqr
    red_chi2=result.redchi



    return a0,m0,n0,t0,chi2,red_chi2



def fit_routine_piecewise2(freq,spec):
    chi2s=100*np.zeros([320])
    turnover=np.zeros([320])
    mymodel = Model(piecewise2)

    for i in np.arange(320):
        turnover[i]=i+30
        params = mymodel.make_params(a=-2, m=-1.2e-2, n=3e-5,b=-4, p=-1.2e-2, q=3e-5,t=turnover[i])
        params['t'].vary = False
        result = mymodel.fit(spec, params, x=freq)
        chi2s[i]=result.chisqr

    t0=turnover[np.argmin(chi2s)]
    params = mymodel.make_params(a=-2, m=-1.2e-2, n=3e-5,b=-4, p=-1.2e-2, q=3e-5,t=t0)
    params['t'].vary = False
    result = mymodel.fit(spec, params, x=freq)

    a0=result.best_values['a']
    m0=result.best_values['m']
    n0=result.best_values['n']
    b0=result.best_values['b']
    p0=result.best_values['p']
    q0=result.best_values['q']
    chi2=result.chisqr
    red_chi2=result.redchi

    return a0,m0,n0,t0,b0,p0,q0,chi2,red_chi2

