import matplotlib
matplotlib.use('Agg')
from pylab import *
import sys
import os
import time
from os.path import exists
import analyzeoutfile
from mpi4py import MPI
import scipy.io
myrank = MPI.COMM_WORLD.rank
MPIsize = MPI.COMM_WORLD.size
if myrank==0:
  print "Starting runmodel.py, myrank="+str(myrank)
Ca_input_flux  = 600.0   #What is the unit?
Ca_input_dur   = 0.005   #ms
Ca_input_freq  = 100      #Hz
Ca_input_N     = 100     #how many successive inputs
Ca_input_onset = 800     #s

L_input_flux   = 2.0
Glu_input_flux   = 50.0
ACh_input_flux   = 20.0
Ntrains        = 1
Duration = 1000          #ms
tolerance = 0.2          
OutputInterval = 1
Nseeds = 10
check_if_exists = 0

Ntrains = 1
trainT = 3000
initfile = ''
Nseeds = 10
addition = ''
if len(sys.argv) > 1:
  Duration = int(sys.argv[1])
if len(sys.argv) > 2:
  tolerance = float(sys.argv[2])
if len(sys.argv) > 3:
  Ca_input_onset = float(sys.argv[3])
if len(sys.argv) > 4:
  Ca_input_N     = int(sys.argv[4])
if len(sys.argv) > 5:
  Ca_input_freq  = float(sys.argv[5])
if len(sys.argv) > 6:
  Ca_input_dur   = float(sys.argv[6])
if len(sys.argv) > 7:
  Ca_input_flux  = float(sys.argv[7])
if len(sys.argv) > 8:
  L_input_flux   = float(sys.argv[8])
if len(sys.argv) > 9:
  Glu_input_flux   = float(sys.argv[9])
if len(sys.argv) > 10:
  ACh_input_flux   = float(sys.argv[10])
if len(sys.argv) > 11:
  Ntrains  = int(float(sys.argv[11]))
if len(sys.argv) > 12:
  trainT  = float(sys.argv[12])
if len(sys.argv) > 13:
  initfile = sys.argv[13]
if len(sys.argv) > 14:
  OutputInterval = float(sys.argv[14])
if len(sys.argv) > 15:
  Nseeds = int(sys.argv[15])
if len(sys.argv) > 16:
  addition = sys.argv[16]


T = 1000./Ca_input_freq
DATA_all = []
for iseed in range(0,Nseeds):
  if iseed%MPIsize != myrank:
    continue
  myseed = 1+iseed
  templ_file = open('Model_template.xml')
  bodyname = 'tstop'+str(Duration)+'_tol'+str(tolerance)+addition+'_onset'+str(Ca_input_onset)+'_n'+str(Ca_input_N)+'_freq'+str(Ca_input_freq)+'_dur'+str(Ca_input_dur)+'_flux'+str(Ca_input_flux)+'_Lflux'+str(L_input_flux)+'_Gluflux'+str(Glu_input_flux)+'_AChflux'+str(ACh_input_flux)+'_Ntrains'+str(Ntrains)+'_trainT'+str(trainT)+'_seed'+str(myseed)
  print 'rank='+str(myrank)+',opened Model_template.xml, writing the model file to model_'+str(bodyname)+'.xml, stimulus file to Stim_'+str(bodyname)+'.xml'
  output_file = open('model_'+str(bodyname)+'.xml','w')
  for line in templ_file:
    newline = line
    if line.find('Stimulations_template.xml') > -1:
      newline = '    <xi:include href="Stim_'+str(bodyname)+'.xml" />\n'
    if line.find('simulationSeed') > -1:
      newline = '    <simulationSeed>    '+str(myseed)+'         </simulationSeed>\n'
    if line.find('runtime') > -1:
      newline = '    <runtime>'+str(Duration)+'</runtime>\n'
    if line.find('outputInterval') > -1:
      newline = '    <outputInterval>'+str(OutputInterval)+'</outputInterval>\n'
    if line.find('tolerance') > -1:
      newline = '    <tolerance>'+str(tolerance)+'</tolerance>\n'
    if line.find('IC_singlecompartment.xml') > -1 and len(initfile) > 0 and initfile != 'None':
      if not exists('IC_from_'+initfile+'.xml'):
        if myrank == 0:
          print 'rank='+str(myrank)+', Generating IC_from_'+initfile+'.xml...'
          #print 'rank='+str(myrank)+', initvalues = analyzeoutfile.getlastvalues(initfile), initfile='+initfile
          DATA = analyzeoutfile.getlastvalues(initfile)
          #print 'rank='+str(myrank)+', analyzeoutfile.writeICfile(DATA[0],DATA[1]), DATA[0]='+str(DATA[0])+', DATA[1]='+str(DATA[1])
          analyzeoutfile.writeICfile(DATA[0],DATA[1],'IC_from_'+initfile+'.xml')
        else:
          time.sleep(3)
      #print 'rank='+str(myrank)+', Using IC_from_'+initfile+'.xml as the initial condiction file'
      newline = '    <xi:include href="IC_from_'+initfile+'.xml" />\n'
    output_file.write(newline)
  templ_file.close()
  output_file.close()

  output_file = open('Stim_'+str(bodyname)+'.xml','w')
  output_file.write('<StimulationSet>\n')

  output_file.write('<InjectionStim specieID="Ca" injectionSite="dendrite:submembrane">\n')
  output_file.write('  <onset>'+str(Ca_input_onset)+'</onset>\n')
  output_file.write('  <duration>'+str(Ca_input_dur)+'</duration>\n')
  output_file.write('  <rate>'+str(Ca_input_flux)+'</rate>\n')
  output_file.write('  <period>'+str(T)+'</period>\n')
  output_file.write('  <end>'+str(Ca_input_onset+T*Ca_input_N+Ca_input_dur)+'</end>\n')
  output_file.write('  <interTrainInterval>'+str(trainT - ((Ca_input_N-1)*T+Ca_input_dur))+'</interTrainInterval>\n')
  output_file.write('  <numTrains>'+str(Ntrains)+'</numTrains>\n')
  output_file.write('</InjectionStim>\n')
  output_file.write('<InjectionStim specieID="L" injectionSite="dendrite:submembrane">\n')
  output_file.write('  <onset>'+str(Ca_input_onset)+'</onset>\n')
  output_file.write('  <duration>'+str(Ca_input_dur)+'</duration>\n')
  output_file.write('  <rate>'+str(L_input_flux)+'</rate>\n')
  output_file.write('  <period>'+str(T)+'</period>\n')
  output_file.write('  <end>'+str(Ca_input_onset+T*Ca_input_N+Ca_input_dur)+'</end>\n')
  output_file.write('  <interTrainInterval>'+str(trainT - ((Ca_input_N-1)*T+Ca_input_dur))+'</interTrainInterval>\n')
  output_file.write('  <numTrains>'+str(Ntrains)+'</numTrains>\n')
  output_file.write('</InjectionStim>\n')
  output_file.write('<InjectionStim specieID="Glu" injectionSite="dendrite:submembrane">\n')
  output_file.write('  <onset>'+str(Ca_input_onset)+'</onset>\n')
  output_file.write('  <duration>'+str(Ca_input_dur)+'</duration>\n')
  output_file.write('  <rate>'+str(Glu_input_flux)+'</rate>\n')
  output_file.write('  <period>'+str(T)+'</period>\n')
  output_file.write('  <end>'+str(Ca_input_onset+T*Ca_input_N+Ca_input_dur)+'</end>\n')
  output_file.write('  <interTrainInterval>'+str(trainT - ((Ca_input_N-1)*T+Ca_input_dur))+'</interTrainInterval>\n')
  output_file.write('  <numTrains>'+str(Ntrains)+'</numTrains>\n')
  output_file.write('</InjectionStim>\n')
  output_file.write('<InjectionStim specieID="ACh" injectionSite="dendrite:submembrane">\n')
  output_file.write('  <onset>'+str(Ca_input_onset)+'</onset>\n')
  output_file.write('  <duration>'+str(Ca_input_dur)+'</duration>\n')
  output_file.write('  <rate>'+str(ACh_input_flux)+'</rate>\n')
  output_file.write('  <period>'+str(T)+'</period>\n')
  output_file.write('  <end>'+str(Ca_input_onset+T*Ca_input_N+Ca_input_dur)+'</end>\n')
  output_file.write('  <interTrainInterval>'+str(trainT - ((Ca_input_N-1)*T+Ca_input_dur))+'</interTrainInterval>\n')
  output_file.write('  <numTrains>'+str(Ntrains)+'</numTrains>\n')
  output_file.write('</InjectionStim>\n')

  output_file.write('</StimulationSet>\n')
  output_file.close()


  timenow = time.time()
  if not check_if_exists or not exists('model_'+str(bodyname)+'.out'):
    if check_if_exists:
      print 'rank='+str(myrank)+', model_'+str(bodyname)+'.out does not yet exist'
    print 'rank='+str(myrank)+', java -jar neurord-3.2.3-all-deps.jar -Dneurord.writers=text model_'+str(bodyname)+'.xml'
    os.system('java -jar neurord-3.2.3-all-deps.jar -Dneurord.writers=text model_'+str(bodyname)+'.xml')
    print 'rank='+str(myrank)+', Simulation done in '+str(time.time()-timenow)+' seconds'
  else:
    print 'rank='+str(myrank)+', model_'+str(bodyname)+'.out exists'

#Synchronize
data = (myrank+1)
data = MPI.COMM_WORLD.gather(data, root=0)
if myrank == 0:
  for i in range(0,MPIsize):
    assert data[i] == (i+1)
else:
  assert data is None
time.sleep(15)

DATA_all = []
bodyname_noseed = 'tstop'+str(Duration)+'_tol'+str(tolerance)+addition+'_onset'+str(Ca_input_onset)+'_n'+str(Ca_input_N)+'_freq'+str(Ca_input_freq)+'_dur'+str(Ca_input_dur)+'_flux'+str(Ca_input_flux)+'_Lflux'+str(L_input_flux)+'_Gluflux'+str(Glu_input_flux)+'_AChflux'+str(ACh_input_flux)+'_Ntrains'+str(Ntrains)+'_trainT'+str(trainT)

if myrank == 0:
 seeds = []
 for iseed in range(0,Nseeds):
  myseed = 1 + iseed
  bodyname = 'tstop'+str(Duration)+'_tol'+str(tolerance)+addition+'_onset'+str(Ca_input_onset)+'_n'+str(Ca_input_N)+'_freq'+str(Ca_input_freq)+'_dur'+str(Ca_input_dur)+'_flux'+str(Ca_input_flux)+'_Lflux'+str(L_input_flux)+'_Gluflux'+str(Glu_input_flux)+'_AChflux'+str(ACh_input_flux)+'_Ntrains'+str(Ntrains)+'_trainT'+str(trainT)+'_seed'+str(myseed)
  if not exists('model_'+str(bodyname)+'.out'):
    continue
  seeds.append(myseed)
  input_file = open('model_'+str(bodyname)+'.out','r')
  #input_file = open(output_data_filename)
  firstline = input_file.readline()
  line = input_file.readline()
  Nvoxels = int(line[11:])
  coords = []
  for ivox in range(0,Nvoxels):
    coords.append(input_file.readline())
  headers_all = []
  values_all = []
  line = input_file.readline()
  times = []
  if len(line) == 0:
    print "iseed="+str(iseed)+": first line of model_"+str(bodyname)+".out empty!!"
  while len(line) > 0:
    header_strs = line.split()
    #print "Appending output t = "+str(header_strs[3])
    headers_all.append(header_strs[:])
    values_thist = []
    for ivox in range(0,Nvoxels):
      line = input_file.readline()
      value_strs = line.split()
      values_thist.append([int(x) for x in value_strs])
    values_all.append(values_thist[:])
    line = input_file.readline()
    times.append(float(header_strs[3]))
  input_file.close()

  mesh_input_file = open('model_'+str(bodyname)+'-mesh.txt.out','r')
  mesh_firstline = mesh_input_file.readline()
  mesh_secondline = mesh_input_file.readline()
  mesh_values = mesh_secondline.split()
  my_volume = float(mesh_values[-2])*1e-15 #litres
  mesh_input_file.close()

  DATA = array(values_all)[:,0,:] #DATA in the first (and only) compartment
  DATA_all.append(DATA[:])
 DATA = mean(array(DATA_all),axis=0)
 scipy.io.savemat(bodyname_noseed+'_'+str(Nseeds)+'seeds.mat', {'DATA': DATA, 'headers': header_strs, 'seeds': seeds})

if myrank == 0:
 for iseed in range(0,Nseeds):
  myseed = 1 + iseed
  bodyname = 'tstop'+str(Duration)+'_tol'+str(tolerance)+addition+'_onset'+str(Ca_input_onset)+'_n'+str(Ca_input_N)+'_freq'+str(Ca_input_freq)+'_dur'+str(Ca_input_dur)+'_flux'+str(Ca_input_flux)+'_Lflux'+str(L_input_flux)+'_Gluflux'+str(Glu_input_flux)+'_AChflux'+str(ACh_input_flux)+'_Ntrains'+str(Ntrains)+'_trainT'+str(trainT)+'_seed'+str(myseed)
  print 'rm model_'+str(bodyname)+'.out'
  os.system('rm model_'+str(bodyname)+'.out')
  print 'rm Stim_'+str(bodyname)+'.xml'
  os.system('rm Stim_'+str(bodyname)+'.xml')
  print 'rm model_'+str(bodyname)+'.log'
  os.system('rm model_'+str(bodyname)+'.log')
  print 'rm model_'+str(bodyname)+'-mesh.txt.out'
  os.system('rm model_'+str(bodyname)+'-mesh.txt.out')
  print 'rm model_'+str(bodyname)+'.xml'
  os.system('rm model_'+str(bodyname)+'.xml')
