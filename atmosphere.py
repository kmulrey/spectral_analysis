import numpy as np
import re


r_e = 6.371 * 1e6  # radius of Earth
h_max = 11282900.2  # height above sea level where the mass overburden vanishes, cm

def get_atmosphere(h,atm):

    #h and layers in cm
    y=0
    layers=atm[0]
    a=atm[1]
    b=atm[2]
    c=atm[3]

    if h>=layers[0] and h<layers[1]:
        y=a[0] + b[0] * np.exp(-1 * h / c[0])
    if h>=layers[1] and h<layers[2]:
        y=a[1] + b[1] * np.exp(-1 * h / c[1])
    if h>=layers[2] and h<layers[3]:
        y=a[2] + b[2] * np.exp(-1 * h / c[2])
    if h>=layers[3] and h<layers[4]:
        y=a[3] + b[3] * np.exp(-1 * h / c[3])
    if h>=layers[4]:
        y=atmA[4]-1*atmB[4]*h/atmC[4]
    return y




def get_i_at(at,atm):
    layers=atm[0]
    a=atm[1]
    b=atm[2]
    c=atm[3]
    
    if at > get_atmosphere(layers[0], atm):
        i = 0
    elif at > get_atmosphere(layers[1], atm):
        i = 1
    elif at > get_atmosphere(layers[2], atm):
        i = 2
    elif at > get_atmosphere(layers[3], atm):
        i = 3
    else:
        i = 4
    if i == 4:
        h = -1. * c[i] * (at - a[i]) / b[i]
    else:
        h = -1. * c[i] * np.log((at - a[i]) / b[i])
    return h


    
def get_vertical_height(at,atm):
    return get_i_at(at,atm)




def return_density(h,atm):
    rho=0
 
    layers=atm[0]
    atmA=atm[1]
    atmB=atm[2]
    atmC=atm[3]

    if h>=layers[0] and h<layers[1]:
        rho=atmB[0]/atmC[0]*np.exp(-1*h/atmC[0])
     
    if h>layers[1] and h<layers[2]:
        rho=atmB[1]/atmC[1]*np.exp(-1*h/atmC[1])
    if h>layers[2] and h<layers[3]:
        rho=atmB[2]/atmC[2]*np.exp(-1*h/atmC[2])
    if h>layers[3] and h<layers[4]:
        rho=atmB[3]/atmC[3]*np.exp(-1*h/atmC[3])
    if h>=layers[4]:
        rho=-1*atmB[4]/atmC[4]

    return rho


def get_distance_xmax(zenith, xmax, observation_level=760.): # obs cm
    '''
    input:
        - xmax in g/cm^2
        - zenith in radians
        output: distance to xmax in g/cm^2
    '''
    return get_atmosphere(observation_level) - xmax

def get_distance_for_height_above_ground(h, zenith, observation_level=760.):
    """ inverse of get_height_above_ground() """
    r = r_e + observation_level
    return (h ** 2 + 2 * r * h + r ** 2 * np.cos(zenith) ** 2) ** 0.5 - r * np.cos(zenith)

    
def get_distance_xmax_geometric(zenith, xmax, atm,observation_level=760.):
    '''
    input:
        - xmax in g/cm^2
        - zenith in radians
        - atm params
        output: distance to xmax in cm
    '''

    h = get_vertical_height(xmax,atm) - observation_level
    return get_distance_for_height_above_ground(h, zenith, observation_level)


def get_vertical_height(T,atm):
    layers=atm[0]
    a=atm[1]
    b=atm[2]
    c=atm[3]
    
    if T > get_atmosphere(layers[1], atm):
        i = 0
    elif T > get_atmosphere(layers[2], atm):
        i = 1
    elif T > get_atmosphere(layers[3], atm):
        i = 2
    elif T > get_atmosphere(layers[4], atm):
        i = 3
    else:
        i = 4
    if i == 4:
        h = -1. * c[i] * (T - a[i]) / b[i]
    else:
        h = -1. * c[i] * np.log((T - a[i]) / b[i])
    return h
