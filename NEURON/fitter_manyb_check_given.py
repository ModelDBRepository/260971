#16.8.2019: fitter_manyb.py: Fix GluR1_ratio, allow the important molecules to vary
#14.8.2019: Use custom Ca limits, separate for each imeas. Use simple flux.
#14.7.2019: Fixed the itimepoint in calculation of the error function. Fixed the duration: earlier the second time point was not even included in the simulated time
#cp ../spineM7b/fitter_manyks.py fitter_saez-briones.py
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
import protocols_many
import scipy.io
from os.path import exists
import os
import calcconds

MAXSIMTIME = 2000

#try:
#  from pylab import *
#except:
#  print "pylab not available"

#pars = c_[[r**pow*c for c in basicparams]] + c_[[(1-r**pow)*c for c in params]]

VARIABLES = [["Caflux",50,5000], #the upper limit of Caflux will be changed according to imeas
             ["Lflux",0.001,5.0],
             ["Gluflux",2,200],
             ["GluR1_ratio",0.0,1.0],
             ["IC_fixedbuffer",0.0,2.0],
             ["IC_NCX",0.0,2.0],
             ["IC_CaM",0.0,5.0],
             ["IC_CK",0.0,5.0],
             ["IC_PKC",0.0,5.0],
             ["IC_MGluR",0.0,5.0],
             ["IC_Gqabg",0.0,5.0],
             ["IC_AC1",0.0,2.0],
             ["IC_AC8",0.0,2.0],
             ["IC_PP1",0.0,2.0],
             ["IC_PP2B",0.0,2.0],
             ["IC_PKA",0.0,2.0],
             ["IC_PDE1",0.0,2.0],
             ["IC_PDE4",0.0,2.0],
             ["IC_I1",0.0,2.0],
             ["IC_PLC",0.0,2.0],
             ["IC_PP2A",0.0,2.0],
             ]
             
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

imeas = 0
maxPerGroupCheck = 1
maxerr = 1.0
maxcaerr = 0.0
istart = 0

argv_i_imeas = 1
argv_i_myseed = 2
argv_i_Nsamp = 3
argv_i_maxPerGroup = 4
argv_i_maxerr = 5
argv_i_maxcaerr = 6
if len(sys.argv) > 4 and sys.argv[0].find("nrniv") > -1:
  argv_i_imeas = 3
  argv_i_myseed = 4
  argv_i_Nsamp = 5
  argv_i_maxPerGroup = 6
  argv_i_maxerr = 7
  argv_i_maxcaerr = 8
if len(sys.argv) > 5 and (sys.argv[1].find("mpi") > -1 or sys.argv[2].find("mpi") > -1):
  argv_i_imeas = 4
  argv_i_myseed = 5
  argv_i_Nsamp = 6
  argv_i_maxPerGroup = 7
  argv_i_maxerr = 8
  argv_i_maxcaerr = 9

Nsamp = 100
myseed = 1
print "argv_i_myseed = "+str(argv_i_myseed)+", argv_i_Nsamp = "+str(argv_i_Nsamp)
if len(sys.argv) > argv_i_imeas:
  imeas = int(float(sys.argv[argv_i_imeas]))
if len(sys.argv) > argv_i_myseed:
  myseed = int(float(sys.argv[argv_i_myseed]))
if len(sys.argv) > argv_i_Nsamp:
  Nsamp = int(float(sys.argv[argv_i_Nsamp]))
if len(sys.argv) > argv_i_maxPerGroup:
  maxPerGroupCheck = int(float(sys.argv[argv_i_maxPerGroup]))
if len(sys.argv) > argv_i_maxerr:
  maxerr = float(sys.argv[argv_i_maxerr])
if len(sys.argv) > argv_i_maxcaerr:
  maxcaerr = float(sys.argv[argv_i_maxcaerr])
if len(sys.argv) > argv_i_maxcaerr+1:
  istart = int(sys.argv[argv_i_maxcaerr+1])

#Caflux_limits = [20000, 20000, 5000, 20000, 20000, 20000, 50000, 5000, 5000, 20000, 20000]
Caflux_limits = [20000, 20000, 13000, 13000, 50000, 10000, 40000, 20000, 20000, 16000, 16000] #Planned so that [Ca flux]*T_total_input is around 2e6, but for imeas=7,8, two different protocols used - something in the middle taken
VARIABLES[0][2] = Caflux_limits[imeas]

Measurement_protocol = protocols_many.get_measurement_protocol()
MeasurementsAll =      Measurement_protocol[0]
Experiments =          Measurement_protocol[1]
protoparams_fixed =    Measurement_protocol[2]
protoparams_var =      Measurement_protocol[3]
Measured_species     = Measurement_protocol[4]
Quantification_types = Measurement_protocol[5]
Measurements_txtsAll = Measurement_protocol[6]

#These are fixed across experiments:
Duration         = protoparams_fixed['Duration']
tolerance        = protoparams_fixed['tolerance']
Ca_input_onset   = protoparams_fixed['Ca_input_onset']

#These are indexed 0-7 according to experiments:
Ca_input_NsAll    = protoparams_var['Ca_input_Ns']
Ca_input_freqsAll = protoparams_var['Ca_input_freqs']
Ca_input_dursAll  = protoparams_var['Ca_input_durs']
NtrainsAll        = protoparams_var['Ca_input_Ntrains']
trainTsAll        = protoparams_var['Ca_input_trainTs']

Measurements =     MeasurementsAll[imeas]
iExperiments =     Measurements[0]
targetTs =         Measurements[1]
targetVals =       Measurements[2]
Measurement_txts = Measurements_txtsAll[imeas]
OBJECTIVES = ['f'+str(i) for i in range(0,len(Measurements[0])+1)]

Ca_input_flux  = 2500.
L_input_flux   = 10.
Glu_input_flux = 20.
ACh_input_flux = 20.

cols = ['#0000FF','#FF0000','#00BB00']
nrnfactor = 6.022e23*my_volume*1e-9 #how many molecules is 1nM in the given volume (0.5 fl)

MAXERR = 1e8 # Maximum error for missing data

def run_model(parameters,saveFig="",deleteFiles=True,rankID=0):
  data = []
  
  paramkeys = parameters.keys()
  Ca_input_flux  = parameters['Caflux'] if 'Caflux' in parameters.keys() else 0.0
  L_input_flux   = parameters['Lflux'] if 'Lflux' in parameters.keys() else 0.0
  Glu_input_flux = parameters['Gluflux'] if 'Gluflux' in parameters.keys() else 0.0
  ACh_input_flux = parameters['AChflux'] if 'AChflux' in parameters.keys() else 0.0
  GluR1_ratio = parameters['GluR1_ratio'] if 'GluR1_ratio' in parameters.keys() else 0.5
  addition_IC = ''
  addition_IC_values = ''
  addition_ks = ''
  addition_ks_values = ''
  for iparam in range(0,len(paramkeys)):
    if paramkeys[iparam][0:3] == 'IC_':
      addition_IC = addition_IC + '-'+paramkeys[iparam][3:]
      addition_IC_values = addition_IC_values + '-'+str(parameters[paramkeys[iparam]])
    elif paramkeys[iparam][0:3] == 'ks[':
      ks_ids = paramkeys[iparam][3:paramkeys[iparam].find(']')].split(',')
      for ik in range(0,len(ks_ids)):
        addition_ks = addition_ks + '-'+ks_ids[ik]
        addition_ks_values = addition_ks_values + '-'+str(parameters[paramkeys[iparam]])
    elif paramkeys[iparam].find('flux') == -1 and paramkeys[iparam].find('GluR1_ratio') == -1:
      print 'Param '+paramkeys[iparam]+' not recognized, rankID='+str(rankID)
  
  addition_IC = addition_IC + '-GluR1-GluR1_memb-GluR2-GluR2_memb'
  addition_IC_values = addition_IC_values + '-'+str(2.*GluR1_ratio)+'-'+str(2.*GluR1_ratio)+'-'+str(2*(1-GluR1_ratio))+'-'+str(2*(1-GluR1_ratio))
  if len(addition_IC) > 0:
    addition_IC = addition_IC[1:]
    addition_IC_values = addition_IC_values[1:] 
  else:
    addition_IC = 'Ca'
    addition_IC_values = '1.0'
  if len(addition_ks) > 0:
    addition_ks = addition_ks[1:]
    addition_ks_values = addition_ks_values[1:]
  #else:
  #  addition_ks = '0'
  #  addition_ks_values = '1.0'
 
  timesAll = []
  timeCoursesAll = []
  maxValsAll = []
  if len(saveFig) > 0:
    close("all")
    f,axarr = subplots(len(iExperiments),2)
    if len(iExperiments) == 1:
      axarr = array([axarr])
  filenames = []
  timenow = time.time()
  timetmp = time.time()
  print str(parameters)
  doQuit = False
  for iiexperiment in range(0,len(iExperiments)):
    iexperiment = iExperiments[iiexperiment]
    iStimulusProtocol = Experiments[iexperiment][0]
    Ca_input_flux_coeff = Experiments[iexperiment][1]
    L_input_flux_coeff = Experiments[iexperiment][2]
    Glu_input_flux_coeff = Experiments[iexperiment][3]
    ACh_input_flux_coeff = Experiments[iexperiment][4]
    Blocked = Experiments[iexperiment][5]
    Altered = Experiments[iexperiment][6]

    if doQuit:
      print 'Skipping, time = '+str(time.time()-timenow)
      timesAll.append([])
      timeCoursesAll.append([])
      maxValsAll.append([])
      continue

    Ca_input_N    = Ca_input_NsAll[iStimulusProtocol]
    Ca_input_freq = Ca_input_freqsAll[iStimulusProtocol]
    Ca_input_dur  = Ca_input_dursAll[iStimulusProtocol]
    Ntrains       = NtrainsAll[iStimulusProtocol]
    trainT        = trainTsAll[iStimulusProtocol]

    addition_IC_this = addition_IC
    addition_IC_values_this = addition_IC_values
    addition_ks_this = addition_ks
    addition_ks_values_this = addition_ks_values
    
    if Blocked != 'None':
      xInd = Blocked.rfind('x')
      if xInd == -1:
        print "Something's wrong with the Blocked!!"
        time.sleep(5)
      addition_IC_this = addition_IC_this + '-' + Blocked[0:xInd]
      addition_IC_values_this = addition_IC_values_this + '-' + Blocked[xInd+1:]

    if len(Altered) > 0:
      for ialtered in range(0,len(Altered[0])):
        addition_ks_this = addition_ks_this + '-'*(ialtered != 0 or len(addition_ks_this) > 0) + str(Altered[0][ialtered])
        addition_ks_values_this = addition_ks_values_this + '-'*(ialtered != 0 or len(addition_ks_values_this) > 0) + str(Altered[1])

    randomString = ''.join(random.choice('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789') for _ in range(0,10))+'_'+str(rankID)
    print('python model_nrn_altered_noU_simpleflux_extfilename_lowmem.py '+str(Duration)+' '+str(tolerance)+' '+str(Ca_input_onset)+' '+str(Ca_input_N)+' '+str(Ca_input_freq)+' '+str(Ca_input_dur)+' '+
           str(Ca_input_flux*Ca_input_flux_coeff)+' '+str(Ca_input_onset-600000)+' '+str(L_input_flux*L_input_flux_coeff)+' '+str(Glu_input_flux*Glu_input_flux_coeff)+' '+str(ACh_input_flux*ACh_input_flux_coeff)+' '+
           str(Ntrains)+' '+str(trainT)+' None fit'+randomString+' '+addition_IC_this+' '+addition_IC_values_this+' '+addition_ks_this+' '+addition_ks_values_this+' > /dev/null 2>&1')
    os.system('python model_nrn_altered_noU_simpleflux_extfilename_lowmem.py '+str(Duration)+' '+str(tolerance)+' '+str(Ca_input_onset)+' '+str(Ca_input_N)+' '+str(Ca_input_freq)+' '+str(Ca_input_dur)+' '+
           str(Ca_input_flux*Ca_input_flux_coeff)+' '+str(Ca_input_onset-600000)+' '+str(L_input_flux*L_input_flux_coeff)+' '+str(Glu_input_flux*Glu_input_flux_coeff)+' '+str(ACh_input_flux*ACh_input_flux_coeff)+' '+
           str(Ntrains)+' '+str(trainT)+' None fit'+randomString+' '+addition_IC_this+' '+addition_IC_values_this+' '+addition_ks_this+' '+addition_ks_values_this+' > /dev/null 2>&1')
    print('Exp. '+str(iiexperiment)+' ('+str(iexperiment)+'), ID='+str(rankID)+' done in '+str(time.time()-timenow)+' sec')
    if time.time()-timetmp > MAXSIMTIME:
      print "MAXSIMTIME exceeded"
      doQuit = True
    timetmp = time.time()
    filename = 'fit'+randomString
    filenames.append(filename)
    #print 'loading filename = '+filename
    if not exists(filename+'.mat'):
      print 'Error: filename = '+filename+' does not exists, Exp. ='+str(iiexperiment)
      timesAll.append([])
      timeCoursesAll.append([])
      maxValsAll.append([])
      continue

    conds, times, cas = calcconds.calcconds_nrn_withcas(filename+'.mat')
    A = scipy.io.loadmat(filename+'.mat')
    timesAll.append(times[:])
    timeCoursesAll.append(conds[:])
    if A['maxDATA'].shape[0] == 1:
      maxValsAll.append(A['maxDATA'][0][1])
    else:
      maxValsAll.append(A['maxDATA'][1][0])

    if len(saveFig) > 0:
      axarr[iiexperiment,0].plot(times,cas,label='Ca')
      axarr[iiexperiment,0].plot([times[0],times[-1]],[maxValsAll[-1]]*2,'k--')
      axarr[iiexperiment,1].plot(times,conds,label='g')
      for iax in range(0,2):
        for tick in axarr[iiexperiment,iax].xaxis.get_major_ticks() + axarr[iiexperiment,iax].yaxis.get_major_ticks():
          tick.label.set_fontsize(4)


  if len(saveFig) > 0:
    f.savefig(saveFig+'_run.eps')
  if deleteFiles:
    for filename in filenames:
      os.system('rm '+filename+'.mat')
      #print('rm '+filename+'.mat')
  return [timesAll, timeCoursesAll, maxValsAll]


def func_to_optimize(parameters,saveFig="",deleteFiles=True,rankID=0):
  A = run_model(parameters, saveFig, deleteFiles,rankID)
  timesAll = A[0]
  timeCoursesAll = A[1]
  maxValsAll = A[2]

  mydict = {}
  if len(saveFig) > 0:
    close("all")
    f,axarr = subplots(len(iExperiments)+1,1)

  objs = []
  mydict['f'+str(len(iExperiments))] = 0
  for iobjective in range(0,len(iExperiments)):
    mydict['f'+str(iobjective)] = 0
    if len(timesAll[iobjective]) == 0:
      mydict['f'+str(iobjective)] = MAXERR
      mydict['f'+str(len(iExperiments))] = mydict['f'+str(len(iExperiments))] + MAXERR
      continue
      
    valsThis = []
    for itarget in range(0,len(targetTs)):
      itime = argmin(abs(timesAll[iobjective]-targetTs[itarget]))
      myval = timeCoursesAll[iobjective][itime]/timeCoursesAll[iobjective][0]
      valsThis.append(myval)
      if isnan(targetVals[iobjective][itarget]):
        continue
      mydict['f'+str(iobjective)] = mydict['f'+str(iobjective)] + abs(targetVals[iobjective][itarget] - myval)
    if len(saveFig) > 0:
      axarr[iobjective].plot(targetTs,valsThis,'b.-')
      axarr[iobjective].plot(targetTs,targetVals[iobjective],'r.-')
      axarr[iobjective].set_title(Measurement_txts[1+iobjective],fontsize=7)
      axarr[iobjective].text(mean(targetTs),mean(axarr[iobjective].get_ylim()),'f = '+str(mydict['f'+str(iobjective)]),fontsize=4)
    mydict['f'+str(len(iExperiments))] = mydict['f'+str(len(iExperiments))] + (maxValsAll[iobjective] > 0.002)*maxValsAll[iobjective]
  if len(saveFig) > 0:
    for iobjective in range(0,len(iExperiments)):
      axarr[len(iExperiments)].bar(iobjective+1,maxValsAll[iobjective])
    axarr[len(iExperiments)].set_xlim([-1,len(iExperiments)+1])
    axarr[len(iExperiments)].set_ylim([0,0.002])
    axarr[len(iExperiments)].set_title('max Ca')
    axarr[len(iExperiments)].text(1.0,mean(axarr[len(iExperiments)].get_ylim()),str(parameters),fontsize=3)
    f.suptitle(Measurement_txts[0])
    axarr[len(iExperiments)].text(1.0,axarr[len(iExperiments)].get_ylim()[1]*0.9,'f = '+str(mydict['f'+str(len(iExperiments))]),fontsize=4)
    f.text(0.01,0.01,str(mydict),fontsize=2)
    

  if len(saveFig) > 0:
    f.savefig(saveFig+'_obj.eps')
  return mydict, A


# After each generation this function is called
def checkpopulation(population, columns, gen, gensdonealready):
  picklelist = [population, columns]
  print 'Generation '+str(gen)+' done, saving to manyb'+str(imeas)+'_seed'+str(myseed)+'_N'+str(N_samples)+'_tmp'+str(gensdonealready+gen)+'.sav'
  file=open('fitfiles/manyb'+str(imeas)+'_seed'+str(myseed)+'_N'+str(N_samples)+'_tmp'+str(gensdonealready+gen)+'.sav', 'w')
  pickle.dump(picklelist,file)
  file.close() 

# Parameters:
# N: size of population
# C: size of capacity 
N_samples = Nsamp
C_samples = 2*Nsamp
N_generations = 50
# eta_m_0, eta_c_0: defines the initial strength of the mution and crossover parameter (large values mean weak effect)
# p_m: probabily of mutation of a parameter (holds for each parameter independently)

random.seed(imeas*1000000+myseed+istart)
filename = 'fitfiles/manyb'+str(imeas)+'_seed'+str(myseed)+'_N'+str(N_samples)
filenamerun = 'fitfiles/rungiven_manyb'+str(imeas)+'_seed'+str(myseed)+'_N'+str(N_samples)

myemoo = emoo.Emoo(N = N_samples, C = C_samples, variables = VARIABLES, objectives = OBJECTIVES)
myemoo.setup(eta_m_0 = 20, eta_c_0 = 20, p_m = 0.5)

A=scipy.io.loadmat('fits_goodparams_manyb.mat')
if A['params'].shape[1] == 1:
  A['params'] = A['params'].T
goodparams = A['params'][0][imeas-7]

counter = 0
for isamp in range(istart,istart+min(maxPerGroupCheck,len(goodparams))):
  paramdict = {}
  for iparam in range(0,len(VARIABLES)):
    paramdict[VARIABLES[iparam][0]] = goodparams[isamp][iparam]
  if counter%myemoo.comm.size == myemoo.comm.rank:
    print "rank = "+str(myemoo.comm.rank)+", running isamp = "+str(isamp)+", counter = "+str(counter)+", file "+'fitfiles/rungiven_'+filename+'_maxerr'+str(maxerr)+'_maxcaerr'+str(maxcaerr)+'_'+str(isamp)
    mydict,A = func_to_optimize(paramdict,filenamerun+'_maxerr'+str(maxerr)+'_maxcaerr'+str(maxcaerr)+'_'+str(isamp),True)
    picklelist = [mydict, A]
    file=open(filenamerun+'_maxerr'+str(maxerr)+'_maxcaerr'+str(maxcaerr)+'_'+str(isamp)+'.sav','w')
    pickle.dump(picklelist,file)
    file.close()       
  counter = counter + 1
