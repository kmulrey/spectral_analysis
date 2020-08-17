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
atm_dir=atm_dir='/vol/astro7/lofar/sim/pipeline/atmosphere_files/'

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
    polfield=np.zeros([nTotalAnt,dlength,2])

    time=np.zeros([nTotalAnt,dlength])

    for l in np.arange(nTotalAnt):
        antenna_position_hold=np.asarray([float(lines[l].split(" ")[2]),float(lines[l].split(" ")[3]),float(lines[l].split(" ")[4])])#read antenna position...
        antenna_file_hold=(lines[l].split(" ")[5].split()[0])   #... and output filename from the antenna list file
        antenna_files.append(antenna_file_hold)
        antenna_positions=np.concatenate((antenna_positions,[antenna_position_hold]))


    nantennas=len(antenna_files)

    file=open(longfile,'r')
    param_list=(re.findall("PARAMETERS.*",file.read()))[0]
    xmax=(float(param_list.split()[4]))
    file.close()
    file=open(steerfile,'r')
    az_list=re.findall("PHI.*",file.read())[0]
    azimuth=np.mod(float(az_list.split()[1]),360.0)*np.pi/180.#rad; CORSIKA coordinates
    az_rot=3*np.pi/2+azimuth
     
    file.seek(0)
    zenith_list=(re.findall("THETAP.*",file.read()))[0]
    zenith=float(zenith_list.split()[1])*np.pi/180. #rad; CORSIKA coordinates

    file.seek(0)
    energy_list=(re.findall("ERANGE.*",file.read()))[0]
    energy=float(energy_list.split()[1])#GeV
    file.close()
     
    for j in np.arange(nantennas):

        antenna_file = lines[j].split(" ")[5]
        coreasfile = '{0}/SIM{1}_coreas/raw_{2}.dat'.format(datadir,str(fileno).zfill(6),antenna_files[j])

        data=np.genfromtxt(coreasfile)
        data[:,1:]*=2.99792458e4 # convert Ex, Ey and Ez (not time!) to Volt/meter
        dlength=data.shape[0]
        poldata=np.ndarray([dlength,2])

        XYZ=np.zeros([dlength,3])
        XYZ[:,0]=-data[:,2] #conversion from CORSIKA coordinates to 0=east, pi/2=north
        XYZ[:,1]=data[:,1]
        XYZ[:,2]=data[:,3]

        XYZ[:,0]=np.roll(XYZ[:,0], 800)
        XYZ[:,1]=np.roll(XYZ[:,1], 800)
        XYZ[:,2]=np.roll(XYZ[:,2], 800)

        UVW=GetUVW(XYZ,0,0,0,zenith,az_rot,1.1837)
        poldata[:,0] = -1.0/np.sin(zenith)*XYZ[:,2] # -1/sin(theta) *z
        poldata[:,1] = -1*np.sin(az_rot)*XYZ[:,0] + np.cos(az_rot)*XYZ[:,1] # -sin(phi) *x + cos(phi)*y in coREAS 0=positive y, 1=negative x
    
    
        polfield[j]=poldata
        efield[j]=UVW#data[:,1:]#UVW#
        time[j]=data.T[0]

        temp=np.copy(antenna_positions)
    antenna_positions[:,0], antenna_positions[:,1], antenna_positions[:,2] = -1*(temp[:,1])/100.,(temp[:,0])/100., temp[:,2]/100.
    ant_pos_uvw=GetUVW(antenna_positions, 0, 0, 0, zenith, az_rot,1.1837)

    return antenna_positions,ant_pos_uvw,time,efield,polfield,zenith,az_rot,energy,xmax,XYZ
    



def ProcessSim(datadir,fileno):
 
    lSample=128*2
    lFFT=lSample/2+1
 
    lFFTkeep=20
 
    longfile = '{0}/DAT{1}.long'.format(datadir,str(fileno).zfill(6))
    steerfile = '{0}/steering/RUN{1}.inp'.format(datadir,str(fileno).zfill(6))
    listfile = open('{0}/steering/SIM{1}.list'.format(datadir,str(fileno).zfill(6)))
    lines = listfile.readlines()
    nTotalAnt=len(lines)
 
    antenna_positions=np.zeros([0,3])
    antenna_files=[]
 
    for l in np.arange(nTotalAnt):
        antenna_position_hold=np.asarray([float(lines[l].split(" ")[2]),float(lines[l].split(" ")[3]),float(lines[l].split(" ")[4])])#read antenna position...
        antenna_file_hold=(lines[l].split(" ")[5].split()[0])   #... and output filename from the antenna list file
        antenna_files.append(antenna_file_hold)
        antenna_positions=np.concatenate((antenna_positions,[antenna_position_hold]))

    lowco=30.0
    hico=80.0
    nantennas=len(antenna_files)

    onskypower=np.zeros([nantennas,2])
    #antenna_position=np.zeros([nantennas,3])
    filteredpower=np.zeros([nantennas,2])
    power=np.zeros([nantennas,2])
    power11=np.zeros([nantennas,2])
    power21=np.zeros([nantennas,2])
    power41=np.zeros([nantennas,2])
    peak_time=np.zeros([nantennas,2])
    peak_bin=np.zeros([nantennas,2])
    peak_amplitude=np.zeros([nantennas,2])
    pol_angle=np.zeros([nantennas])
    pol_angle_filt=np.zeros([nantennas])


    file=open(longfile,'r')
    param_list=(re.findall("PARAMETERS.*",file.read()))[0]
    xmax=(float(param_list.split()[4]))
    file.close()
    file=open(steerfile,'r')
    az_list=re.findall("PHI.*",file.read())[0]
    azimuth=np.mod(float(az_list.split()[1]),360.0)*np.pi/180.#rad; CORSIKA coordinates
    az_rot=3*np.pi/2+azimuth
     
    file.seek(0)
    zenith_list=(re.findall("THETAP.*",file.read()))[0]
    zenith=float(zenith_list.split()[1])*np.pi/180. #rad; CORSIKA coordinates

    file.seek(0)
    energy_list=(re.findall("ERANGE.*",file.read()))[0]
    energy=float(energy_list.split()[1])#GeV
    file.close()

    fullFFT=np.zeros([nantennas,lFFTkeep])
    fullPower=np.zeros([nantennas])
    fullFrequencies=np.zeros([lFFTkeep])
    wfX=np.zeros([nantennas,81])
    wfY=np.zeros([nantennas,81])
    wfZ=np.zeros([nantennas,81])
    time_all=np.zeros([81])
    wfX=np.zeros([nantennas,81])
    wfY=np.zeros([nantennas,81])
    wfZ=np.zeros([nantennas,81])

    for j in np.arange(nantennas):
     #for j in np.arange(1):

        antenna_file = lines[j].split(" ")[5]
        coreasfile = '{0}/SIM{1}_coreas/raw_{2}.dat'.format(datadir,str(fileno).zfill(6),antenna_files[j])

        data=np.genfromtxt(coreasfile)
        data[:,1:]*=2.99792458e4 # convert Ex, Ey and Ez (not time!) to Volt/meter
        dlength=data.shape[0]
        poldata=np.ndarray([dlength,2])
        XYZdata=np.ndarray([dlength,2])
        az_rot=3*np.pi/2+azimuth    #conversion from CORSIKA coordinates to 0=east, pi/2=north
        zen_rot=zenith
        XYZ=np.zeros([dlength,3])
        XYZ[:,0]=-data[:,2] #conversion from CORSIKA coordinates to 0=east, pi/2=north
        XYZ[:,1]=data[:,1]
        XYZ[:,2]=data[:,3]
     

        # Convert to, v, vxB, vxvxB coordinates to compute Stokes parameters and polarization angle
        UVW=GetUVW(XYZ,0,0,0,zen_rot,az_rot,1.1837)
        alpha= GetAlpha(zen_rot,az_rot,1.1837)
     
        poldata[:,0] = UVW[:,0]
        poldata[:,1] = UVW[:,1]

        spec=np.fft.rfft(poldata, axis=-2)
        #print spec.shape
        # Apply antenna model
        tstep = data[1,0]-data[0,0]
        onskypower[j]=np.array([np.sum(poldata[:,0]*poldata[:,0]),np.sum(poldata[:,1]*poldata[:,1])])*tstep


        freqhi = 0.5/tstep/1e6 # MHz
        freqstep = freqhi/(dlength/2+1) # MHz
        frequencies = np.arange(0,freqhi,freqstep)*1e6 # Hz
        frequencies = np.arange(0,dlength/2+1)*freqstep*1e6

        #Apply window and reduce maximum frequency to acquire downsampled signal
        fb = int(np.floor(lowco/freqstep))
        lb = int(np.floor(hico/freqstep)+1)
        window = np.zeros([1,dlength/2+1,1])
        window[0,fb:lb+1,0]=1
     
        pow0=np.abs(spec[:,0])*np.abs(spec[:,0])
        pow1=np.abs(spec[:,1])*np.abs(spec[:,1])
     
     
        ospow0=np.abs(spec[:,0])*np.abs(spec[:,0])
        ospow1=np.abs(spec[:,1])*np.abs(spec[:,1])
        power[j]=np.array([np.sum(pow0[fb:lb+1]),np.sum(pow1[fb:lb+1])])/(dlength/2.)*tstep
        filteredpower[j]=np.array([np.sum(ospow0[fb:lb+1]),np.sum(ospow1[fb:lb+1])])/(dlength/2.)*tstep
     
        # assume that simulated time resolution is higher than LOFAR time resolution (t_step=5 ns)
        maxfreqbin= int(np.floor(tstep/5e-9 * dlength/2.)+1)
        shortspec=np.array([spec[0:maxfreqbin,0]*window[0,0:maxfreqbin,0],spec[0:maxfreqbin,1]*window[0,0:maxfreqbin,0]])
        filt=np.fft.irfft(shortspec, axis=-1)
     
        # after downsampling, renormalize the signal!
        dlength_new=filt.shape[1]
        filt=filt*1.0*dlength_new/dlength
        # to calculate the time of arrival upsample with a factor 5
        filt_upsampled=resample(filt,5*dlength_new,axis=-1)
        # compute hilbert enevelope
        hilbenv=np.abs(hilbert(filt,axis=-1))
        hilbenv_upsampled=np.abs(hilbert(filt_upsampled,axis=-1))
     
        # peak_time is the bin where the maximum is located; NOT the actual time of the peak!
        peak_bin[j]=np.argmax(hilbenv,axis=-1)
        peak_time[j]=np.argmax(hilbenv_upsampled,axis=-1)*1e-9 #in seconds
        peak_amplitude[j]=np.max(hilbenv_upsampled,axis=-1)
        if (peak_amplitude[j,0]>peak_amplitude[j,1]):
            pt=peak_bin[j,0]
        else:
            pt=peak_bin[j,1]
     
        # for 3 different window size, the total power is calculated. The window is allowed to `wrap around', so some voodoo is needed to determine the range:
        d=filt.shape[1]
        rng=5
        a=int(np.max([0,pt-rng]))
        b=int(pt+rng+1)
        c=int(np.min([d,pt+d-rng]))
        power11[j]=(np.sum(np.square(filt[:,a:b]),axis=-1)+np.sum(np.square(filt[:,c:d]),axis=-1))*5e-9
        rng=10
        a=int(np.max([0,pt-rng]))
        b=int(pt+rng+1)
        c=int(np.min([d,pt+d-rng]))
        power21[j]=(np.sum(np.square(filt[:,a:b]),axis=-1)+np.sum(np.square(filt[:,c:d]),axis=-1))*5e-9
        rng=20
        a=int(np.max([0,pt-rng]))
        b=int(pt+rng+1)
        c=int(np.min([d,pt+d-rng]))
        power41[j]=(np.sum(np.square(filt[:,a:b]),axis=-1)+np.sum(np.square(filt[:,c:d]),axis=-1))*5e-9

    temp=np.copy(antenna_positions)
    antenna_positions[:,0], antenna_positions[:,1], antenna_positions[:,2] = -1*(temp[:,1])/100.,(temp[:,0])/100., temp[:,2]/100.

    azimuth=3*np.pi/2+azimuth #  +x = east (phi=0), +y = north (phi=90)
 
    ant_pos_uvw=GetUVW(antenna_positions, 0, 0, 0, zenith, az_rot,1.1837)

    return zenith, azimuth, alpha, energy, hillas, antenna_positions, ant_pos_uvw, power11/(377.0),power21/(377.0),power41/(377.0)


def getEM(datadir,fileno):

    longfile = '{0}/DAT{1}.long'.format(datadir,str(fileno).zfill(6))
    lookup1='LONGITUDINAL ENERGY DEPOSIT'
    lookup2='FIT OF THE HILLAS CURVE'
    start=-1
    stop=-1

    with open(longfile) as myFile:
        for num, line in enumerate(myFile, 1):
            if lookup1 in line:
                start=num+2
            if lookup2 in line:
                stop=num-1

    myFile.close()
    #myFile=open(longfile,'r')
    longinfo=np.genfromtxt(longfile,skip_header=start-1,max_rows=(stop-start+1))

    #myFile.close()


    total_dep= np.sum(longinfo.T[9])
    em_dep=np.sum(longinfo.T[1]+longinfo.T[2]+longinfo.T[3])
    other_dep=np.sum(longinfo.T[4]+longinfo.T[5]+longinfo.T[6]+longinfo.T[7]+longinfo.T[8])

    return em_dep,other_dep,total_dep


def get_atm(event):

    # read atmospheric parameters ATMLAY, A, B, C respectively for event
    filename=atm_dir+'ATMOSPHERE_'+event+'.DAT'
    #file=open(filename,'r')
    atm=np.genfromtxt(filename,skip_header=1,max_rows=4)
    #file.close()
    return atm

