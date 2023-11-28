#cp simsteadystate_someonly.py simsteadystate_li2020.py
#simsteadystate.py: Simulate the Ca binding and protein activation time courses for steady Ca inputs with varying Ca input fluxes
#Input: the amplitude of the fluxes of other ligands: beta-adrenergic, glutamate and ACh
#Tuomo Maki-Marttunen, 2019-2020
#CC BY 4.0
import matplotlib
matplotlib.use('Agg')
import numpy as np
import emoo
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
Ca_input_onset   = 40000.0
Ca_input_N       = 1
Ca_input_freq    = 1.0
Ca_input_dur     = 300000.0

Ca_input_fluxes = [0.000025,0.00005,0.0001,0.00025,0.0005,0.001,0.0015,0.002]+[0.025*i for i in range(1,10)]+[0.25,0.3,0.35,0.4,0.45]+[0.25*i for i in range(2,10)] +[2.5,3.0,3.5,4.0,4.5]+ [2.5*i for i in range(2,10)]
Ntrains          = 1
cols = ['#0000FF','#FF0000','#00BB00']
trainT = 100000.0
nrnfactor = 6.022e23*my_volume*1e-9 #how many molecules is 1nM in the given volume (0.5 fl)

MAXERR = 1e8 # Maximum error for missing data

fb_coeff = 1.0
CaM_tot = 5.0
calbin_coeff = 1.0

if len(sys.argv) > 1:
  fb_coeff = float(sys.argv[1])
if len(sys.argv) > 2:
  CaM_tot = float(sys.argv[2])
if len(sys.argv) > 3:
  calbin_coeff = float(sys.argv[3])


def run_model(parameters,saveFig="",deleteFiles=False,rankID=0):
  data = []
  
  paramkeys = parameters.keys()
  L_input_flux   = 0.0
  Glu_input_flux = 0.0
  ACh_input_flux = 0.0


  species_alts = [['CaOut', 1.0],
                  ['Leak', 0.0],
                  ['Calbin', calbin_coeff],
                  ['CalbinC', calbin_coeff],
                  ['LOut', 0.0],
                  ['Epac1', 0.0],
                  ['PMCA', 0.0],
                  ['NCX', 0.0],
                  ['L', 0.0],
                  ['R', 0.0],
                  ['Gs', 0.0],
                  ['Gi', 0.0],
                  ['AC1', 0.0],
                  ['ATP', 0.0],
                  ['AC8', 0.0],
                  ['PDE1', 0.0],
                  ['AMP', 0.0],
                  ['Ng', 2.5],            # originally 20 uM
                  ['CaM', CaM_tot/60.0],  # originally 60 uM
                  ['PP2B', 21.739],       # originally 2.3 uM
                  ['CK', 2.1739],         # originally 23 uM
                  ['PKA', 0.0],
                  ['I1', 0.0],
                  ['PP1', 0.0],
                  ['GluR1', 0.0],
                  ['GluR1_memb', 0.0],
                  ['PDE4', 0.0],
                  ['fixedbuffer', fb_coeff],
                  ['MGluR', 0.0],
                  ['GluOut', 0.0],
                  ['Gqabg', 0.0],
                  ['PLC', 0.0],
                  ['Pip2', 0.0],
                  ['PIkinase', 0.0],
                  ['Ip3degPIk', 0.0],
                  ['PKC', 0.0],
                  ['DAG', 0.0],
                  ['DAGK', 0.0],
                  ['DGL', 0.0],
                  ['CaDGL', 0.0],
                  ['DAGCaDGL', 0.0],
                  ['Ip3degrad', 0.0],
                  ['GluR2', 0.0],
                  ['GluR2_memb', 0.0],
                  ['PP2A', 0.0],
                  ['M1R', 0.0],
                  ['PLA2', 0.0]]

  addition_IC = ''
  addition_IC_values = ''
  for ispecalt in range(0,len(species_alts)):
    if species_alts[ispecalt][1] != 1.0:
      addition_IC = addition_IC + species_alts[ispecalt][0] + '-'
      addition_IC_values = addition_IC_values + str(species_alts[ispecalt][1]) + '-'
  addition_IC = addition_IC[0:-1]
  addition_IC_values = addition_IC_values[0:-1]

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
    myString = 'nrn_li2020_'+str(fb_coeff)+'_'+str(CaM_tot)+'_'+str(calbin_coeff)+'_caflux'+str(iexperiment)

    filename = myString
    print "thisfile = "+filename
    filenames.append(filename)
    #print 'loading filename = '+filename
    if not exists(filename+'.mat'):
      print('python model_nrn_altered_noU_extfilename_lowmem_recall.py '+str(Duration)+' '+str(tolerance)+' '+str(Ca_input_onset)+' '+str(Ca_input_N)+' '+str(Ca_input_freq)+' '+str(Ca_input_dur)+' '+
             str(Ca_input_fluxes[iexperiment])+' '+str(L_input_flux)+' '+str(Glu_input_flux)+' '+str(ACh_input_flux)+' '+
             str(Ntrains)+' '+str(trainT)+' None '+myString+' '+addition_IC+' '+addition_IC_values)
      os.system('python model_nrn_altered_noU_extfilename_lowmem_recall.py '+str(Duration)+' '+str(tolerance)+' '+str(Ca_input_onset)+' '+str(Ca_input_N)+' '+str(Ca_input_freq)+' '+str(Ca_input_dur)+' '+
             str(Ca_input_fluxes[iexperiment])+' '+str(L_input_flux)+' '+str(Glu_input_flux)+' '+str(ACh_input_flux)+' '+
             str(Ntrains)+' '+str(trainT)+' None '+myString+' '+addition_IC+' '+addition_IC_values)
      print('Exp. '+str(iexperiment)+', ID='+str(rankID)+' done in '+str(time.time()-timenow)+' sec')
      if not exists(filename+'.mat'):
        print 'Error: filename = '+filename+' does not exists, Exp. ='+str(iexperiment)
        timesAll.append([])
        timeCoursesAll.append([])
        maxValsAll.append([])
        DATA_all_all_all.append([])
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

    itime = argmin(abs(times - 4339000))
    timesAll.append(times[itime])
    timeCoursesAll.append(conds[itime])
    maxValsAll.append(maxDATA_all['Ca'])
    DATA_all_all_all.append([DATA_all_all['DATA'][ikey][itime] for ikey in range(0,len(DATA_all_all['DATA']))])
  if deleteFiles:
    for filename in filenames:
      os.system('rm '+filename+'.mat')
      #print('rm '+filename+'.mat')
  return [timesAll, timeCoursesAll, maxValsAll,  DATA_all_all_all, DATA_all_all['headers']]

if True:
  A = run_model({}, "", False, 0)
  picklelist = [A, Ca_input_fluxes]
  file=open('steadystate_new_li2020_'+str(fb_coeff)+'_'+str(CaM_tot)+'_'+str(calbin_coeff)+'_processed.sav', 'w')
  pickle.dump(picklelist,file)
  file.close()

