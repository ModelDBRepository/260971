#cp model_nrn_testPKA_williamson_varyrates.py model_nrn_testPKA_withdiss_williamson_varyrates.py
#cp model_nrn_testPKA_williamson.py model_nrn_testPKA_williamson_varyrates.py
#cp model_nrn_testPKA.py model_nrn_testPKA_williamson.py
#cp model_nrn_altered_template.py model_nrn_testPKA.py
#Based on Model C in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2719611/
from neuron import h, rxd
from pylab import *
from matplotlib import pyplot
import scipy.io
import time
import re
import mytools

h.load_file('stdrun.hoc')

dend = h.Section(name='dend')
dend.L=1
dend.diam=0.79788
cyt = rxd.Region([dend], name='cyt', nrn_region='i')

mesh_input_file = open('mesh_general.out','r')
mesh_firstline = mesh_input_file.readline()
mesh_secondline = mesh_input_file.readline()
mesh_values = mesh_secondline.split()
my_volume = float(mesh_values[-2])*1e-15 #litres
mesh_input_file.close()

Duration = 3000
tolerance = 1e-6
cAMP_input_onset = 800
cAMP_input_N     = 100
cAMP_input_freq  = 100
cAMP_input_dur   = 0.005
cAMP_input_flux  = 600.0
Ntrains        = 1
trainT = 0
extfilename = ''

ks = [1.0]*10
ks[0] = 8.72e9   # PKA + cAMP*4 <-> PKAcAMP4 (forward) #(8.72 e-17 in Williamson et al. 2009)
ks[1] = 0.0       # PKA + cAMP*4 <-> PKAcAMP4 (backward)
ks[4] = 0.00024   # PKAcAMP4 <-> PKAr + PKAc*2 (forward)
ks[5] = 25.5      # PKAcAMP4 <-> PKAr + PKAc*2 (backward)

ks[6] = 0.001        # AMP --> ATP (forward)
ks[7] = 21.66        # PDE4 + cAMP <-> PDE4cAMP (forward)
ks[8] = 0.0034656    # PDE4 + cAMP <-> PDE4cAMP (backward)
ks[9] = 0.017233     # PDE4cAMP --> PDE4 + AMP (forward)

if len(sys.argv) > 1:
  Duration = int(sys.argv[1])
if len(sys.argv) > 2:
  tolerance = float(sys.argv[2])
if len(sys.argv) > 3:
  cAMP_input_onset = float(sys.argv[3])
if len(sys.argv) > 4:
  cAMP_input_N     = int(sys.argv[4])
if len(sys.argv) > 5:
  cAMP_input_freq  = float(sys.argv[5])
if len(sys.argv) > 6:
  cAMP_input_dur   = float(sys.argv[6])
if len(sys.argv) > 7:
  cAMP_input_flux  = float(sys.argv[7])
if len(sys.argv) > 8:
  ks[0] = float(sys.argv[8])
if len(sys.argv) > 9:
  ks[1] = float(sys.argv[9])
if len(sys.argv) > 10:
  extfilename = sys.argv[10]

initvalues = [0.0, 0.0, 0.0, 0.0064, 0.0, 0.0, 0.00098, 2.0, 0.00067, 0.0]

species = ['cAMP', 'PKAcAMP2', 'PKAcAMP4', 'PKA', 'PKAr', 'PKAc', 'AMP', 'ATP','PDE4','PDE4cAMP']
windows = [0, 0, 0, 2, 1, 1, 3, 3, 4, 4]
tolscales = [1.0 for i in range(0,len(species))]

print "my_volume = "+str(my_volume)+" l ?= "+str(dend.L*(dend.diam/2)**2*3.14159265358)+" um3"
specs = []
for ispec in range(0,len(species)):
  specs.append(rxd.Species(cyt, name='spec'+str(ispec), charge=0, initial=initvalues[ispec], atolscale=tolscales[ispec]))
cAMP_flux_rate = rxd.Parameter(cyt, initial=0)

#specs[0] cAMP
#specs[1] PKA2cAMP
#specs[2] PKA4cAMP
#specs[3] PKA
#specs[4] PKAr
#specs[5] PKAc
#specs[6] AMP
#specs[7] ATP
#specs[8] PDE4
#specs[9] PDE4cAMP
reaction1 = rxd.Reaction(specs[3] + specs[0]*4 <> specs[2], ks[0], ks[1])
reaction3 = rxd.Reaction(specs[2] <> specs[4] + specs[5]*2, ks[4]*specs[2], ks[5]*specs[4]*specs[5], custom_dynamics=True)

reaction4 = rxd.Reaction(specs[6] > specs[7], ks[6])
reaction5 = rxd.Reaction(specs[8] + specs[0] <> specs[9], ks[7], ks[8])
reaction6 = rxd.Reaction(specs[9] > specs[7] + specs[8], ks[9])


reaction_cAMP_flux = rxd.Rate(specs[0], cAMP_flux_rate) # cAMP
vec_t = h.Vector()

vecs = []
vec_t = h.Vector()
vec_t.record(h._ref_t)
for ispec in range(0,len(species)):
  vecs.append(h.Vector())
  vecs[ispec].record(specs[ispec].nodes(dend)(0.5)[0]._ref_concentration)

cvode = h.CVode()
cvode.active(1)
hmax = cvode.maxstep(100)
hmin = cvode.minstep(1e-10)
cvode.atol(tolerance)

h.finitialize(-65)
def set_param(param, val):
    param.nodes.value = val
    h.cvode.re_init()

### Set on and off the inputs to the spine
T = 1000./cAMP_input_freq
unow = 1
tnow = 0
for itrain in range(0,Ntrains):
    for istim in range(0,cAMP_input_N):
      tnew = cAMP_input_onset + istim*T + trainT*itrain
      h.cvode.event(tnew, lambda: set_param(cAMP_flux_rate, cAMP_input_flux/6.022e23/my_volume*1e3))
      h.cvode.event(tnew+cAMP_input_dur, lambda: set_param(cAMP_flux_rate, 0))
      tnow = tnew
timenow = time.time()
h.continuerun(Duration)
print "Simulation done in "+str(time.time()-timenow)+" seconds"
def isFlux(t):
  for itrain in range(0,Ntrains):
    for istim in range(0,cAMP_input_N):
      tnew = cAMP_input_onset + istim*T + trainT*itrain
      if t >= tnew and t < tnew+cAMP_input_dur:
        return 1
  return 0
tvec = array(vec_t)
minDT_nonFlux = 20.0
minDT_Flux = 1.0
lastt = -inf
itvec2 = []
for it in range(0,len(tvec)):
  if tvec[it] - lastt > minDT_nonFlux or (isFlux(tvec[it]) and tvec[it] - lastt > minDT_Flux):
    itvec2.append(it)
    lastt = tvec[it]

headers = [ 'tvec', 'cAMP', 'PKAcAMP2', 'PKAcAMP4', 'PKA', 'PKAr', 'PKAc']
if len(extfilename) > 0:
  scipy.io.savemat(extfilename+'.mat', {'tcDATA': array(vecs), 'times': array(tvec)})
else:
  scipy.io.savemat('cAMP_withdiss_testwilliamson_tstop'+str(Duration)+'_tol'+str(tolerance)+'_onset'+str(cAMP_input_onset)+'_n'+str(cAMP_input_N)+'_freq'+str(cAMP_input_freq)+'_dur'+str(cAMP_input_dur)+'_flux'+str(cAMP_input_flux)+'.mat',
    {'tcDATA': array(vecs), 'times': array(tvec)})

if False:
  f,axarr = subplots(1,max(windows)+1)
  for i in range(0,len(vecs)):
    axarr[windows[i]].plot(tvec,array(vecs[i]),label=species[i])
  for i in range(0,len(axarr)):
    axarr[i].legend(fontsize=6)
  if len(extfilename) > 0:
    f.savefig('cAMP_withdiss_testwilliamson_'+extfilename+'.eps')
  else:
    f.savefig('cAMP_withdiss_testwilliamson_tstop'+str(Duration)+'_tol'+str(tolerance)+'_onset'+str(cAMP_input_onset)+'_n'+str(cAMP_input_N)+'_freq'+str(cAMP_input_freq)+'_dur'+str(cAMP_input_dur)+'_flux'+str(cAMP_input_flux)+'.eps')
