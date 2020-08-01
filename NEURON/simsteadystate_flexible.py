#simsteadystate_flexible.py: Simulate the Ca binding and protein activation time courses for steady Ca inputs with varying Ca input fluxes
#Input: the amplitude of the fluxes of other ligands: beta-adrenergic, glutamate and ACh
#Tuomo Maki-Marttunen, 2019-2020
#CC BY 4.0
import matplotlib
matplotlib.use('Agg')
import numpy as np
from pylab import *
import pickle
import mytools
import time
import random
import sys
import scipy.io
from os.path import exists
import os
import calcconds

globalcounter = 0
             
mesh_input_file = open('mesh_general.out','r')
mesh_firstline = mesh_input_file.readline()
mesh_secondline = mesh_input_file.readline()
mesh_values = mesh_secondline.split()
my_volume = float(mesh_values[-2])*1e-15 #litres
mesh_input_file.close()

conds_hom1 = [12.4, 18.9]
conds_hom2 = 2.2
conds_het = 2.5

for iargv in range(0,len(sys.argv)):
  print "argv["+str(iargv)+"] = "+sys.argv[iargv]

Duration = 4940000
tolerance = 1e-6
Ca_input_onset   = 4040000.0
Ca_input_N       = 1
Ca_input_freq    = 1.0
Ca_input_dur     = 300000.0
Ca_input_fluxes  = [2.5*i for i in range(0,101)] + [250+5*i for i in range(0,101)] + [750+10*i for i in range(0,101)] + [1750+20*i for i in range(0,101)]
Ntrains          = 1
cols = ['#0000FF','#FF0000','#00BB00']
trainT = 100000.0
nrnfactor = 6.022e23*my_volume*1e-9 #how many molecules is 1nM in the given volume (0.5 fl)
ligandFluxes = 0.005
if len(sys.argv) > 1:
  ligandFluxes = float(sys.argv[1])

MAXERR = 1e8 # Maximum error for missing data

def run_model(parameters,saveFig="",deleteFiles=True,rankID=0):
  data = []
  
  paramkeys = parameters.keys()
  L_input_flux   = parameters['Lflux'] if 'Lflux' in parameters.keys() else ligandFluxes
  Glu_input_flux = parameters['AChGluflux'] if 'AChGluflux' in parameters.keys() else ligandFluxes
  ACh_input_flux = parameters['AChGluflux'] if 'AChGluflux' in parameters.keys() else ligandFluxes
  GluR1_ratio = parameters['GluR1_ratio'] if 'GluR1_ratio' in parameters.keys() else 0.5 
  addition_IC = ''
  addition_IC_values = ''
  addition_ks = ''
  addition_ks_values = ''
  for iparam in range(0,len(paramkeys)):
    if paramkeys[iparam] == 'IC_MGluRM1GqPLC':
      addition_IC = addition_IC + '-MGluR-M1-Gqabg-PLC'
      addition_IC_values = addition_IC_values + '-'+str(parameters[paramkeys[iparam]]) + '-'+str(parameters[paramkeys[iparam]])+'-'+str(parameters[paramkeys[iparam]]) + '-'+str(parameters[paramkeys[iparam]])
    elif paramkeys[iparam] == 'IC_RGsAC1AC8':
      addition_IC = addition_IC + '-R-Gs-AC1-AC8'
      addition_IC_values = addition_IC_values + '-'+str(parameters[paramkeys[iparam]]) + '-'+str(parameters[paramkeys[iparam]]) + '-'+str(parameters[paramkeys[iparam]]) + '-'+str(parameters[paramkeys[iparam]])
    elif paramkeys[iparam][0:3] == 'IC_':
      addition_IC = addition_IC + '-'+paramkeys[iparam][3:]
      addition_IC_values = addition_IC_values + '-'+str(parameters[paramkeys[iparam]])
      if paramkeys[iparam] == 'IC_Calbin': #If Calbin is changed, change CalbinC in the same proportion
        addition_IC = addition_IC + '-CalbinC'
        addition_IC_values = addition_IC_values + '-'+str(parameters[paramkeys[iparam]])        
    elif paramkeys[iparam][0:3] == 'ks[':
      ks_ids = paramkeys[iparam][3:paramkeys[iparam].find(']')].split(',')
      for ik in range(0,len(ks_ids)):
        addition_ks = addition_ks + '-'+ks_ids[ik]
        addition_ks_values = addition_ks_values + '-'+str(parameters[paramkeys[iparam]])
    elif paramkeys[iparam].find('flux') == -1 and paramkeys[iparam].find('GluR1_ratio') == -1:
      print 'Param '+paramkeys[iparam]+' not recognized, rankID='+str(rankID)
  
  addition_filename = ''
  addition_IC = addition_IC + '-GluR1-GluR1_memb-GluR2-GluR2_memb'
  addition_IC_values = addition_IC_values + '-'+str(2.*GluR1_ratio)+'-'+str(2.*GluR1_ratio)+'-'+str(2*(1-GluR1_ratio))+'-'+str(2*(1-GluR1_ratio))
  if len(addition_IC) > 0:
    addition_IC = addition_IC[1:]
    addition_IC_values = addition_IC_values[1:] 
    addition_filename = addition_filename+'_'+addition_IC+'x'+addition_IC_values
  else:
    addition_IC = 'Ca'
    addition_IC_values = '1.0'
    addition_filename = addition_filename+'_Cax1.0'
  if len(addition_ks) > 0:
    addition_ks = addition_ks[1:]
    addition_ks_values = addition_ks_values[1:]
    addition_filename = addition_filename+'_k'+addition_ks+'x'+addition_ks_values

  timesAll = []
  timeCoursesAll = []
  maxValsAll = []
  if len(saveFig) > 0:
    close("all")
    f,axarr = subplots(len(Ca_input_fluxes),2)
  filenames = []
  timenow = time.time()
  DATA_all_all_all = []
  for iexperiment in range(0,len(Ca_input_fluxes)):
    randomString = ''.join(random.choice('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789') for _ in range(0,10))+'_'+str(rankID)
    print('python model_nrn_altered_noU_extfilename_lowmem_recall_local.py '+str(Duration)+' '+str(tolerance)+' '+str(Ca_input_onset)+' '+str(Ca_input_N)+' '+str(Ca_input_freq)+' '+str(Ca_input_dur)+' '+
           str(Ca_input_fluxes[iexperiment])+' '+str(L_input_flux)+' '+str(Glu_input_flux)+' '+str(ACh_input_flux)+' '+
           str(Ntrains)+' '+str(trainT)+' None fit'+randomString+' '+addition_IC+' '+addition_IC_values+' '+addition_ks+' '+addition_ks_values)
    os.system('python model_nrn_altered_noU_extfilename_lowmem_recall_local.py '+str(Duration)+' '+str(tolerance)+' '+str(Ca_input_onset)+' '+str(Ca_input_N)+' '+str(Ca_input_freq)+' '+str(Ca_input_dur)+' '+
           str(Ca_input_fluxes[iexperiment])+' '+str(L_input_flux)+' '+str(Glu_input_flux)+' '+str(ACh_input_flux)+' '+
           str(Ntrains)+' '+str(trainT)+' None fit'+randomString+' '+addition_IC+' '+addition_IC_values+' '+addition_ks+' '+addition_ks_values)
    print('Exp. '+str(iexperiment)+', ID='+str(rankID)+' done in '+str(time.time()-timenow)+' sec')

    filename = 'fit'+randomString
    print "thisfile = "+filename
    filenames.append(filename)
    #print 'loading filename = '+filename
    if not exists(filename+'.mat'):
      print 'Error: filename = '+filename+' does not exists, Exp. ='+str(iexperiment)
      timesAll.append([])
      timeCoursesAll.append([])
      maxValsAll.append([])
      continue
    
    conds, times = calcconds.calcconds_nrn(filename+'.mat')

    DATA_all_all = scipy.io.loadmat(filename+'.mat')
    if DATA_all_all['maxDATA'].shape[0] != 1:
      DATA_all_all['maxDATA'] = DATA_all_all['maxDATA'].T
    DATA_all = {}
    maxDATA_all = {}
    for ikey in range(0,len(DATA_all_all['headers'])):
      first_space = DATA_all_all['headers'][ikey].find(' ')
      mykey = DATA_all_all['headers'][ikey]
      if first_space > -1:
        mykey = DATA_all_all['headers'][ikey][0:first_space]
      DATA_all[mykey] = DATA_all_all['DATA'][ikey]
      maxDATA_all[mykey] = DATA_all_all['maxDATA'][0][ikey]

    timesAll.append(times[:])
    timeCoursesAll.append(conds[:])
    maxValsAll.append(maxDATA_all['Ca'])
    if len(saveFig) > 0:
      axarr[iexperiment,0].plot(times,conds,label='g')
      axarr[iexperiment,1].plot(times,DATA_all['Ca'],label='Ca')
      axarr[iexperiment,1].text(mean(times),mean(DATA_all['Ca']),'max(Ca)='+str(maxDATA_all['Ca']))
      for iax in range(0,2):
        for tick in axarr[iexperiment,iax].xaxis.get_major_ticks() + axarr[iexperiment,iax].yaxis.get_major_ticks():
          tick.label.set_fontsize(4)
    DATA_all_all_all.append(DATA_all_all.copy())
  if len(saveFig) > 0:
    f.savefig(saveFig+'_run.eps')
  if deleteFiles:
    for filename in filenames:
      os.system('rm '+filename+'.mat')
      #print('rm '+filename+'.mat')
  return [timesAll, timeCoursesAll, maxValsAll,  DATA_all_all_all]

if True:
  A = run_model({}, "", True, 0)
  picklelist = A
  file=open('steadystate_flux'+str(ligandFluxes)+'_raw.sav', 'w')
  pickle.dump(picklelist,file)
  file.close()

  timesThis = A[0]
  condsThis = A[1]
  maxCasThis = A[2]
  DATA_all_all_all = A[3]

  conds = [condsThis[i][-1] for i in range(0,len(condsThis))]
  baselines = [condsThis[i][0] for i in range(0,len(condsThis))]

