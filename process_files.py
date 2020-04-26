from scipy.signal import hilbert
import numpy as np
import sys
from optparse import OptionParser
import pickle
import re
from scipy.signal import hilbert
from scipy.signal import resample
import scipy.fftpack as fftp
import os
from scipy import signal
from scipy.optimize import curve_fit

conversion_factor_integrated_signal = (1/376.73)#* 6.24150934e18  # to convert V**2/m**2 * s -> J/m**2 -> eV/m**2


def GetUVW(pos, cx, cy, cz, zen, az, Binc):
    relpos = pos-np.array([cx,cy,cz])
    B = np.array([0,np.cos(Binc),-np.sin(Binc)])
    v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
    vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
    vxB = vxB/np.linalg.norm(vxB)
    vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
    return np.array([np.inner(vxB,relpos),np.inner(vxvxB,relpos),np.inner(v,relpos)]).T

def GetAlpha(zen,az,Binc):
    B = np.array([0,np.cos(Binc),-np.sin(Binc)])
    v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
    vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
    vxB = vxB/np.linalg.norm(vxB)
    vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
    return np.arccos(np.inner(np.asarray(B), np.asarray(v)) / (np.linalg.norm(B) * np.linalg.norm(v)))

def get_efield(datadir,fileno):

    dlength=4082
    longfile = '{0}/DAT{1}.long'.format(datadir,str(fileno).zfill(6))
    steerfile = '{0}/steering/RUN{1}.inp'.format(datadir,str(fileno).zfill(6))
    listfile = open('{0}/steering/SIM{1}.list'.format(datadir,str(fileno).zfill(6)))
    lines = listfile.readlines()
    nTotalAnt=len(lines)

    antenna_positions=np.zeros([0,3])
    antenna_files=[]

    efield=np.zeros([nTotalAnt,dlength,3])
    time=np.zeros([nTotalAnt,dlength])

    for l in np.arange(nTotalAnt):
        antenna_position_hold=np.asarray([float(lines[l].split(" ")[2]),float(lines[l].split(" ")[3]),float(lines[l].split(" ")[4])])#read antenna position...
        antenna_file_hold=(lines[l].split(" ")[5].split()[0])   #... and output filename from the antenna list file
        antenna_files.append(antenna_file_hold)
        antenna_positions=np.concatenate((antenna_positions,[antenna_position_hold]))


    nantennas=len(antenna_files)

    '''
    hillas = np.genfromtxt(re.findall("PARAMETERS.*",open(longfile,'r').read()))[2:]
    zenith=(np.genfromtxt(re.findall("THETAP.*",open(steerfile,'r').read()))[1])*np.pi/180. #rad; CORSIKA coordinates
    azimuth=np.mod(np.genfromtxt(re.findall("PHIP.*",open(steerfile,'r').read()))[1],360)*np.pi/180.  #rad; CORSIKA coordinates
    az_rot=3*np.pi/2+azimuth    #conversion from CORSIKA coordinates to 0=east, pi/2=north

    energy=np.genfromtxt(re.findall("ERANGE.*",open(steerfile,'r').read()))[1] #GeV

    for j in np.arange(nantennas):

        antenna_file = lines[j].split(" ")[5]
        coreasfile = '{0}/SIM{1}_coreas/raw_{2}.dat'.format(datadir,str(fileno).zfill(6),antenna_files[j])

        data=np.genfromtxt(coreasfile)
        data[:,1:]*=2.99792458e4 # convert Ex, Ey and Ez (not time!) to Volt/meter
        dlength=data.shape[0]
    
        XYZ=np.zeros([dlength,3])
        XYZ[:,0]=-data[:,2] #conversion from CORSIKA coordinates to 0=east, pi/2=north
        XYZ[:,1]=data[:,1]
        XYZ[:,2]=data[:,3]

        XYZ[:,0]=np.roll(XYZ[:,0], 800)
        XYZ[:,1]=np.roll(XYZ[:,1], 800)
        XYZ[:,2]=np.roll(XYZ[:,2], 800)

        UVW=GetUVW(XYZ,0,0,0,zenith,az_rot,1.1837)
    
    
    
    
        efield[j]=UVW#data[:,1:]#UVW#
        time[j]=data.T[0]

        temp=np.copy(antenna_positions)
    antenna_positions[:,0], antenna_positions[:,1], antenna_positions[:,2] = -1*(temp[:,1])/100.,(temp[:,0])/100., temp[:,2]/100.
    '''
    return antenna_positions,time,efield,zenith,az_rot,energy,hillas[2]
