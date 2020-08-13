import numpy as np
import os
import pickle
import sys
import glob
from scipy import signal
import re
import collect
from os import path
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-e", "--event", default = "132649890", help = "event)
(options, args) = parser.parse_args()
event = int(options.eventid)



base_dir='/vol/astro7/lofar/sim/pipeline/events/'

energy=[]
zenith=[]
azimuth=[]
xmax=[]
hi=[]
rho=[]
rho2=[]
dmax=[]
n_xmax=[]
alpha=[]
clip_ratio=[]
Erad_30_80=[]
Erad_gm_30_80=[]
Erad_ce_30_80=[]
Erad_30_200=[]
Erad_gm_30_200=[]
Erad_ce_30_200=[]
Erad_50_350=[]
Erad_gm_50_350=[]
Erad_ce_50_350=[]
em_dep=[]
total_dep=[]
prim=[]

dir_list_proton=[]
dir_list_iron=[]

if path.isdir(base_dir+'/'+event+'/0/coreas/proton/'):
    dir_list_proton.append(base_dir+'/'+event+'/0/coreas/proton/')
if path.isdir(base_dir+'/'+event+'/1/coreas/proton/'):
    dir_list_proton.append(base_dir+'/'+event+'/1/coreas/proton/')
if path.isdir(base_dir+'/'+event+'/2/coreas/proton/'):
    dir_list_proton.append(base_dir+'/'+event+'/1/coreas/proton/')
                      
                      
if path.isdir(base_dir+'/'+event+'/0/coreas/iron/'):
    dir_list_iron.append(base_dir+'/'+event+'/0/coreas/iron/')
if path.isdir(base_dir+'/'+event+'/1/coreas/iron/'):
    dir_list_iron.append(base_dir+'/'+event+'/1/coreas/iron/')
if path.isdir(base_dir+'/'+event+'/2/coreas/iron/'):
    dir_list_iron.append(base_dir+'/'+event+'/2/coreas/iron/')

for j in np.arange(len(dir_list_proton));
    direc=dir_list_proton[j]
    os.chdir(direc)
    runlist=[]
    primary='proton'

    
    for file in glob.glob("*.long"):
        runlist.append(file.split('.')[0].split('DAT')[1])


    for i in np.arange(len(runlist)):
        try:
            base=direc
            print (event,runlist[i])
            primary='proton'
            energy1,zenith1,azimuth1,xmax1,hi1,rho1,rho21,dmax1,n_xmax1,alpha1,clip_ratio1,Erad_30_80_1,Erad_gm_30_80_1,Erad_ce_30_80_1,Erad_30_200_1,Erad_gm_30_200_1,Erad_ce_30_200_1,Erad_50_350_1,Erad_gm_50_350_1,Erad_ce_50_350_1,em_dep1,total_dep1=collect.return_Srd_info(runlist[i],event,primary,base)
            if primary=='proton':
                prim.append(0)
            if primary=='iron':
                prim.append(1)
            energy.append(energy1)
            zenith.append(zenith1)
            azimuth.append(azimuth1)
            xmax.append(xmax1)
            hi.append(hi1)
            rho.append(rho1)
            rho2.append(rho21)
            dmax.append(dmax1)
            n_xmax.append(n_xmax1)
            alpha.append(alpha1)
            clip_ratio.append(clip_ratio1)
            Erad_30_80.append(Erad_30_80_1)
            Erad_gm_30_80.append(Erad_gm_30_80_1)
            Erad_ce_30_80.append(Erad_ce_30_80_1)
            Erad_30_200.append(Erad_30_200_1)
            Erad_gm_30_200.append(Erad_gm_30_200_1)
            Erad_ce_30_200.append(Erad_ce_30_200_1)
            Erad_50_350.append(Erad_50_350_1)
            Erad_gm_50_350.append(Erad_gm_50_350_1)
            Erad_ce_50_350.append(Erad_ce_50_350_1)
            em_dep.append(em_dep1)
            total_dep.append(total_dep1)
            
            
       except:
        print('error__________')
            
for j in np.arange(len(dir_list_iron));
    direc=dir_list_iron[j]
    os.chdir(direc)
    runlist=[]
    primary='iron'

                  
    for file in glob.glob("*.long"):
        runlist.append(file.split('.')[0].split('DAT')[1])


    for i in np.arange(len(runlist)):
        try:
            base=direc
            print (event,runlist[i])
            primary='iron'
            energy1,zenith1,azimuth1,xmax1,hi1,rho1,rho21,dmax1,n_xmax1,alpha1,clip_ratio1,Erad_30_80_1,Erad_gm_30_80_1,Erad_ce_30_80_1,Erad_30_200_1,Erad_gm_30_200_1,Erad_ce_30_200_1,Erad_50_350_1,Erad_gm_50_350_1,Erad_ce_50_350_1,em_dep1,total_dep1=collect.return_Srd_info(runlist[i],event,primary,base)
            if primary=='proton':
                prim.append(0)
            if primary=='iron':
                prim.append(1)
            energy.append(energy1)
            zenith.append(zenith1)
            azimuth.append(azimuth1)
            xmax.append(xmax1)
            hi.append(hi1)
            rho.append(rho1)
            rho2.append(rho21)
            dmax.append(dmax1)
            n_xmax.append(n_xmax1)
            alpha.append(alpha1)
            clip_ratio.append(clip_ratio1)
            Erad_30_80.append(Erad_30_80_1)
            Erad_gm_30_80.append(Erad_gm_30_80_1)
            Erad_ce_30_80.append(Erad_ce_30_80_1)
            Erad_30_200.append(Erad_30_200_1)
            Erad_gm_30_200.append(Erad_gm_30_200_1)
            Erad_ce_30_200.append(Erad_ce_30_200_1)
            Erad_50_350.append(Erad_50_350_1)
            Erad_gm_50_350.append(Erad_gm_50_350_1)
            Erad_ce_50_350.append(Erad_ce_50_350_1)
            em_dep.append(em_dep1)
            total_dep.append(total_dep1)
        
        except:
            print('error___________')
    
energy=np.asarray(energy)
zenith=np.asarray(zenith)
azimuth=np.asarray(azimuth)
xmax=np.asarray(xmax)
hi=np.asarray(hi)
rho=np.asarray(rho)
rho2=np.asarray(rho2)
dmax=np.asarray(dmax)
n_xmax=np.asarray(n_xmax)
alpha=np.asarray(alpha)
clip_ratio=np.asarray(clip_ratio)
Erad_30_80=np.asarray(Erad_30_80)
Erad_gm_30_80=np.asarray(Erad_gm_30_80)
Erad_ce_30_80=np.asarray(Erad_ce_30_80)

Erad_30_200=np.asarray(Erad_30_200)
Erad_gm_30_200=np.asarray(Erad_gm_30_200)
Erad_ce_30_200=np.asarray(Erad_ce_30_200)

Erad_50_350=np.asarray(Erad_50_350)
Erad_gm_50_350=np.asarray(Erad_gm_50_350)
Erad_ce_50_350=np.asarray(Erad_ce_50_350)

em_dep=np.asarray(em_dep)
total_dep=np.asarray(total_dep)
prim=np.asarray(prim)

info={'energy':energy,'zenith':zenith,'azimuth':azimuth,'xmax':xmax,'hi':hi,'rho':rho,'rho2':rho2,'dmax':dmax,'n_xmax':n_xmax,'alpha':alpha,'clip_ratio':clip_ratio,'Erad_30_80':Erad_30_80,'Erad_gm_30_80':Erad_gm_30_80,'Erad_ce_30_80':Erad_ce_30_80,'Erad_30_200':Erad_30_200,'Erad_gm_30_200':Erad_gm_30_200,'Erad_ce_30_200':Erad_ce_30_200,'Erad_50_350':Erad_50_350,'Erad_gm_50_350':Erad_gm_50_350,'Erad_ce_50_350':Erad_ce_50_350,'em_dep':em_dep,'total_dep':total_dep,'prim':prim}

outfilename='/vol/astro3/lofar/sim/kmulrey/spectral_analysis/Srd_Data/'+event+'.p'
outfile=open(outfilename,'wb')
pickle.dump(info,outfile)
outfile.close()