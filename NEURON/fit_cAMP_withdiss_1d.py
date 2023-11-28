#Check the time course of PKA activation using the single-step model (that follows the law of mass action)
import matplotlib
matplotlib.use('Agg')
from pylab import *
import pickle
import mytools
import time
import random
import sys
import scipy.io
from os.path import exists
import os

cAMP_flux = 0.2

for iargv in range(0,len(sys.argv)):
  print "argv["+str(iargv)+"] = "+sys.argv[iargv]

argv_i_flux = 1
if len(sys.argv) > 4 and sys.argv[0].find("nrniv") > -1:
  argv_i_flux = 3
if len(sys.argv) > 5 and (sys.argv[1].find("mpi") > -1 or sys.argv[2].find("mpi") > -1):
  argv_i_flux = 4

if len(sys.argv) > argv_i_flux:
  cAMP_flux = float(sys.argv[argv_i_flux])

A = scipy.io.loadmat('cAMP_withdiss_test_tstop22000_tol1e-08_onset800.0_n1_freq1.0_dur16000.0_flux'+str(cAMP_flux)+'.mat')
if A['times'].shape[1] == 1:
  A['times'] = A['times'].T

OBJECTIVE_ts = [1000+1000*i for i in range(0,21)]
OBJECTIVE_DATA = [mytools.interpolate(A['times'][0],A['tcDATA'][0]+A['tcDATA'][1],OBJECTIVE_ts), #Time course of cAMP
                  mytools.interpolate(A['times'][0],A['tcDATA'][5],OBJECTIVE_ts)] #Time course of PKAc
OBJECTIVE_ts_highres = [1000+100*i for i in range(0,201)]
OBJECTIVE_DATA_highres = [mytools.interpolate(A['times'][0],A['tcDATA'][0]+A['tcDATA'][1],OBJECTIVE_ts_highres), #Time course of cAMP
                          mytools.interpolate(A['times'][0],A['tcDATA'][5],OBJECTIVE_ts_highres)] #Time course of PKAc

mesh_input_file = open('mesh_general.out','r')
mesh_firstline = mesh_input_file.readline()
mesh_secondline = mesh_input_file.readline()
mesh_values = mesh_secondline.split()
my_volume = float(mesh_values[-2])*1e-15 #litres
mesh_input_file.close()

titles = ['cAMP_flux']
cols = ['#0000FF','#FF0000','#00BB00']
trainT = 1.0

def run_model(parameters,saveFig="",deleteFiles=True,rankID=0):
  data = []
  
  paramkeys = parameters.keys()
  k0 = parameters['k[0]']
  k1 = parameters['k[1]']
  timesAll = []
  timeCoursesAll = []
  filenames = []
  timenow = time.time()

  randomString = ''.join(random.choice('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789') for _ in range(0,10))+'_'+str(rankID)
  print('python model_nrn_testPKA_withdiss_williamson_varyrates.py 22000 1e-8 800 1 1 16000 '+str(cAMP_flux)+' '+str(k0)+' '+str(k1)+' fit'+randomString)
  os.system('python model_nrn_testPKA_withdiss_williamson_varyrates.py 22000 1e-8 800 1 1 16000 '+str(cAMP_flux)+' '+str(k0)+' '+str(k1)+' fit'+randomString)
  print('Exp. 0, ID='+str(rankID)+' done in '+str(time.time()-timenow)+' sec')
  filename = 'fit'+randomString
  if not exists(filename+'.mat'):
    print 'Error: filename = '+filename+' does not exists, Exp. = 0'
    timesAll.append([])
    timeCoursesAll.append([])
      
  DATA_all = scipy.io.loadmat(filename+'.mat')
  times = DATA_all['times']
  tcs = DATA_all['tcDATA']

  timesAll.append(times)
  timeCoursesAll.append(tcs[:])
  if deleteFiles:
    os.system('rm '+filename+'.mat')
  return [timesAll, timeCoursesAll]

def func_to_optimize(parameters,saveFig="",deleteFiles=True,rankID=0):
  A = run_model(parameters, saveFig, deleteFiles,rankID)
  timesAll = A[0]
  timeCoursesAll = A[1]

  mydict = {}
  if len(saveFig) > 0:
    close("all")
    f,axarr = subplots(len(OBJECTIVE_DATA),1)
  objs = []
  imeas = [0, 5] # cAMP and PKAc traces                                                                  
  for iobjective in range(0,len(OBJECTIVE_DATA)):
    objs_thisobj = []
    times = timesAll[0]
    if times.shape[0] == 1:
      times = times.T
    timeCourses = mytools.interpolate(times[0],timeCoursesAll[0][imeas[iobjective]],OBJECTIVE_ts)
    mydict['f'+str(iobjective)] = sum(abs(array(timeCourses)-array(OBJECTIVE_DATA[iobjective])))/max(array(OBJECTIVE_DATA[iobjective]))

    if len(saveFig) > 0:
     try:
      axarr[iobjective].plot(OBJECTIVE_ts_highres,OBJECTIVE_DATA_highres[iobjective],'b--')
      axarr[iobjective].plot(OBJECTIVE_ts,OBJECTIVE_DATA[iobjective],'bx')
      axarr[iobjective].plot(times,timeCoursesAll[0][imeas[iobjective]],'r--')
      axarr[iobjective].plot(OBJECTIVE_ts,timeCourses,'rx')
      #axarr[iobjective].text(times[0],Ylims[1]*0.2+Ylims[0]*0.8,str(objs[iobjective]),fontsize=5)
      #axarr[iobjective].set_ylabel(Objective_titles[iobjective],fontsize=5)
      for tick in axarr[iobjective].xaxis.get_major_ticks() + axarr[iobjective].yaxis.get_major_ticks():
        tick.label.set_fontsize(3.5)
     except:
      print "Something wrong here..."
    Ylims = axarr[iobjective].get_ylim()
    if iobjective == 0:
      axarr[0].text(OBJECTIVE_ts[0],Ylims[1]*0.8+Ylims[0]*0.2,str(parameters),color='#0000FF',fontsize=5)
    else:
      axarr[1].text(OBJECTIVE_ts[0],Ylims[1]*0.8+Ylims[0]*0.2,str(mydict),color='#0000FF',fontsize=5)

  if len(saveFig) > 0:
    picklelist = A
    file=open(saveFig+'_data.sav', 'w')
    pickle.dump(picklelist,file)
    file.close()   
    f.savefig(saveFig+'_obj.eps')
  return mydict

ks = [0.4e9, 1.0e9, 1.6e9, 2.2e9, 2.8e9, 4.0e9]
for ik in range(0,len(ks)):
  paramdict = {}
  paramdict['k[0]'] = ks[ik]
  paramdict['k[1]'] = 6e-5
  A = func_to_optimize(paramdict,"run_cAMP_withdiss_flux"+str(cAMP_flux)+"_k"+str(ks[ik]),False)



