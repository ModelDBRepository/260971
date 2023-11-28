#drawfig5.py: Draws the STDP curves
#Tuomo Maki-Marttunen, 2019-2020
#CC BY 4.0
import matplotlib
matplotlib.use('Agg')
from pylab import *
import scipy.io
import sys
import itertools
from os.path import exists
import mytools
import pickle
from matplotlib.collections import PatchCollection

close("all")

species = [ [ ['GluR1_S831', 'GluR1_S845_S831', 'GluR1_S831_PKAc', 'GluR1_S845_S831_PP1', 'GluR1_S831_PP1', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'] ],
            [ ['GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'] ],
            [ ['GluR1_S845', 'GluR1_S845_S831', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_S845_CKp', 'GluR1_S845_PKCt', 'GluR1_S845_PKCp', 'GluR1_S845_PP1', 'GluR1_S845_S831_PP1', 'GluR1_S845_PP2B', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S845', 'GluR1_memb_S845_S831', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCt', 'GluR1_memb_S845_PKCp', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'] ],
            [ ['GluR1_memb', 'GluR1_memb_S845', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_PKAc', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_CKp', 'GluR1_memb_PKCt', 'GluR1_memb_PKCp', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCt', 'GluR1_memb_S845_PKCp', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'] ],
            [ ['GluR2_S880', 'GluR2_S880_PP2A', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'] ],
            [ ['GluR2_memb', 'GluR2_memb_PKCt', 'GluR2_memb_PKCp', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'] ] ]

titles = ['GluR1 S831 phosph.', 'GluR1-S831 at memb.', 'GluR1 S845 phosph.', 'GluR1 at memb.', 'GluR2 S880 phosph.', 'GluR2 at memb.']

protocol = 'pair'
Gluflux = 0.0
Lfluxmax = 0.05
AChfluxmax = 0.05
fluxtxt = '0.0-0.05-0.05'
CaCoeff = 1.0
icell_l23pc = 1
Econ_l23pc = 0.001
wNMDA_l23pc = 3.2
location_l23pc = 'apic250-300'
Nsyn_l23pc = 10
pulseamp_l23pc = 5.0
Nskip = 1
#TRAINISISALL=[-200.0, -100.0, -50.0, -30.0, -10.0, 0.0, 10.0, 30.0, 50.0, 100.0, 200.0]
TRAINISISALL=[-200.0, -180.0, -160.0, -140.0, -120.0, -100.0, -80.0, -60.0, -50.0, -40.0, -30.0, -25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0]
cols = ['#440154', '#470f62', '#481d6f', '#472a79', '#453681', '#414387', '#3c4f8a', '#37598c', '#32648e', '#2d6f8e', '#29788e', '#26828e', '#228b8d', '#1f958b', '#1e9f88', '#22a884', '#2bb17e', '#3bbb75', '#4dc36b', '#62cb5f', '#7ad251'][::10] + ['#dddd00']

cols_nrn = ['#360043', '#390c4f', '#391759', '#392161', '#372b67', '#34366c', '#303f6f', '#2c4770', '#285071', '#245872', '#216072', '#1e6872', '#1b6f71', '#19776f', '#187f6d', '#1b866a', '#238e65', '#2f955e', '#3e9c56', '#4fa24c', '#61a841'][::10] + ['#cccc00']
dashes_3 = [(2,1), (2,3.5), (3.5,5.0), ()]

irows = [1,1,2,2]
icols = [0,1,0,1]
idashes = [0,3,0,3]
#styles = ['x-','.-','x-','.-']
styles = ['x','.','x','.']
mews = [0.4,0.8,0.4,0.8]

conds_hom1 = [12.4, 18.9]
conds_hom2 = 2.2
conds_het = 2.5

if len(sys.argv) > 1:
  fluxtxt = sys.argv[1]
if len(sys.argv) > 2:
  CaCoeff = float(sys.argv[2])
if len(sys.argv) > 3:
  icell_l23pc  = int(sys.argv[3])
if len(sys.argv) > 4:
  Econ_l23pc  = float(sys.argv[4])
if len(sys.argv) > 5:
  wNMDA_l23pc  = float(sys.argv[5])
if len(sys.argv) > 6:
  location_l23pc = sys.argv[6]
if len(sys.argv) > 7:
  Nsyn_l23pc = int(float(sys.argv[7]))
if len(sys.argv) > 8:
  pulseamp_l23pc = float(sys.argv[8])
if len(sys.argv) > 9:
  Nskip = int(sys.argv[9])

if fluxtxt.find('-') > -1:
  Gluflux = float(fluxtxt.split('-')[0])
  Lfluxmax = float(fluxtxt.split('-')[1])
  AChfluxmax = float(fluxtxt.split('-')[2])
else:
  Gluflux = float(fluxtxt)

ylims = [[0,200],[0,90],[0,200],[90,150],[0,200],[0,260]]
yticks = [[0,100,200],[0,40,80],[0,50,100,150,200],[90,110,130,150],[0,50,100,150,200],[0,50,100,150,200,250]]

mesh_input_file = open('mesh_general.out','r')
mesh_firstline = mesh_input_file.readline()
mesh_secondline = mesh_input_file.readline()
mesh_values = mesh_secondline.split()
my_volume = float(mesh_values[-2])*1e-15 #litres
mesh_input_file.close()
locator_params(nticks=4)

T = 500000.0
Duration = 500000.0
tstop = 500000/1000

f,axs = subplots(4,3)
axarr = [axs[0,0],axs[0,1],axs[0,2],axs[1,0],axs[1,1],axs[1,2],axs[2,0],axs[2,1],axs[2,2],axs[3,0],axs[3,1],axs[3,2]]
axnew = []
axnew.append(f.add_axes([0.015, 0.70, 0.31, 0.23]))
for iax in range(0,3):
  for iay in range(0,2):
    axs[3*iay,iax].set_position([0.3+iax*0.2, 0.735+0.13*iay, 0.2, 0.11])
    axs[3*iay,iax].spines['top'].set_visible(False)
    axs[3*iay,iax].spines['right'].set_visible(False)
    if iay == 1:
      axnew.append(f.add_axes([0.13+0.3+iax*0.2-0.115*(iax==2), 0.02+0.735+0.13*1, 0.055, 0.06]))
      axnew[-1].spines['top'].set_visible(False)
      axnew[-1].spines['right'].set_visible(False)

  for iay in range(0,2):
    axs[iay+1,iax].set_position([0.14+iax*0.287, 0.45-0.17*iay, 0.19, 0.16])
    axs[iay+1,iax].spines['top'].set_visible(False)
    axs[iay+1,iax].spines['right'].set_visible(False)

foundone = 0


#Pre-saved file with the data for plotting the morphologies
unpicklefile = open('../L23PC/morph_accurate_segdata_icell'+str(icell_l23pc)+'.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
segdata = unpickledlist[:]

########################## Draw the morphology of the L23PC #################################
for ipos in range(0,len(segdata[0])):
  coord1 = segdata[0][ipos][0]
  coord2 = segdata[0][ipos][1]
  diam = segdata[0][ipos][2]
  mydist = segdata[0][ipos][3]
  mytree = segdata[0][ipos][4]
  if coord1[1] < -180 or coord2[1] < -180:
    continue
  hilight = mydist >= 250 and mydist <= 300 and mytree == 0
  if diam < 4:
    axnew[0].plot([coord1[0], coord2[0]],[coord1[1], coord2[1]],'k-',color='#000000' if hilight else '#808080', linewidth=diam*(0.5+0.5*hilight))
  else:
    axnew[0].plot([coord1[0], coord2[0]],[coord1[1], coord2[1]],'k-',color='#808080', linewidth=1.0)
axnew[0].axis('equal')

if icell_l23pc == 1:
  polygon = Polygon(array([[-400,-110,-390],[240,180,260]]).T, True)
else:
  polygon = Polygon(array([[-300,10,-298],[300,270,320]]).T, True)
  if icell_l23pc != 5:
    print "polygon for icell ="+str(icell_l23pc)+" not adjusted"
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor('#000000')
p.set_edgecolor('none')
axnew[0].add_collection(p)
xs = range(-480,-280)
xs0 = [i*0.02 for i in range(0,200)]
ys0 = [(1-exp(-xs0[i]/0.2))*exp(-xs0[i]/1.5) for i in range(0,len(xs0))]
ys = [270 + ys0[i]*100 for i in range(0,200)]
axnew[0].plot([-510,-230],[270,270],'k-',color='#808080', linewidth=0.5)
axnew[0].plot(xs,ys,'k-', linewidth=0.5)

polygon = Polygon(array([[-300,0,-290],[60,0,80]]).T, True)
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor('#000000')
p.set_edgecolor('none')
axnew[0].add_collection(p)
#axnew[0].plot([-510,-460,-460,-460,-410,-410,-410,-360,-360,-360,-310,-310,-310,-230],[100,100,160,100,100,160,100,100,160,100,100,160,100,100],'k-', linewidth=0.5)
axnew[0].plot([-510,-460,-460,-450,-450,-410,-410,-400,-400,-360,-360,-350,-350,-310,-310,-300,-300,-230],[100,100,160,160,100,100,160,160,100,100,160,160,100,100,160,160,100,100],'k-', linewidth=0.5)

if icell_l23pc == 1:
  #polygon = Polygon(array([[-200,10,-182],[400,180,418]]).T, True)
  polygon = Polygon(array([[-320,-110,-302],[360,180,375]]).T, True)
  #polygon = Polygon(array([[-320,-110,-302],[360,270,375]]).T, True)
else:
  #polygon = Polygon(array([[-200,10,-182],[400,180,418]]).T, True)
  polygon = Polygon(array([[-200,10,-182],[360,270,378]]).T, True)
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor('#808080')
p.set_edgecolor('none')
axnew[0].add_collection(p)
xs = range(-480,-180)
xs0 = [i*0.02 for i in range(0,300)]
ys0 = [(1-exp(-xs0[i]/0.6))*exp(-xs0[i]/2.0) for i in range(0,len(xs0))]
ts0 = [0.4,1.4,2.4,3.4]
for it in range(0,len(ts0)):
  for iy in range(0,len(ys0)):
    ys0[iy] = ys0[iy] + 0.25*(1-exp(-(xs0[iy]-ts0[it])/0.1))*exp(-(xs0[iy]-ts0[it])/0.3)*(xs0[iy] >= ts0[it])
ys = [390 + ys0[i]*100 for i in range(0,300)]
axnew[0].plot([-510,-130],[390,390],'k-',color='#808080', linewidth=0.5)
axnew[0].plot(xs,ys,'k-',color='#808080', linewidth=0.5)
axnew[0].set_axis_off()
axnew[0].patch.set_alpha(0.0)
axnew[0].set_xlim([-570, 830])

axnew[0].plot([-400-80*(icell_l23pc==5),-200-80*(icell_l23pc==5)],[-120,-120],'k-',linewidth=3)

ISIdts = [-80.0, -30.0, 20.0]
########################## Draw the Ca inputs at -dt, 0, and dt ISIs #################################
if True: 
  Ca_input_freq = 1.0
  Ca_input_coeff = 1.0
  frac_Ca = 0.1
  for iISI in range(0,len(ISIdts)):
    ISIdt = ISIdts[iISI]
    A = scipy.io.loadmat('nrn_noninterp_tstop5000000_tol1e-06_onset4040000.0_n120_freq'+str(Ca_input_freq)+'_dur3.0_CaCoeff1.0_Lflux0.05_Gluflux0.0_AChflux0.05_contnm3560000.0_600000.0_Ntrains1_trainT100000.0_pair'+str(ISIdt)+'_icell'+str(icell_l23pc)+'_pulseamp5.0_Nsyn'+str(Nsyn_l23pc)+'_Ninputs1_Econ'+str(Econ_l23pc)+'_wNMDA'+str(wNMDA_l23pc)+'_'+location_l23pc+'.mat')
    istart = argmin(abs(A['DATA'][0]-4158800)) #Take the second to last (119th) second of stimulation
    iend = argmin(abs(A['DATA'][0]-4159800))
    Ca_concs = A['DATA'][1,istart:iend]

    MAT = scipy.io.loadmat('../L23PC/currClips'+str(icell_l23pc-1)+'_imut0_neckLen0.5_neckDiam0.1_stimfreq'+str(Ca_input_freq)+'_pulseamp5.0_Nsyn'+str(Nsyn_l23pc)+'_Ninputs1'+location_l23pc+'_Econ'+str(Econ_l23pc)+'_wNMDA'+str(wNMDA_l23pc)+'_Npulses4_ISI'+str(ISIdt)+'_withV.mat')
    print "loaded "+'../L23PC/currClips'+str(icell_l23pc-1)+'_imut0_neckLen0.5_neckDiam0.1_stimfreq'+str(Ca_input_freq)+'_pulseamp5.0_Nsyn'+str(Nsyn_l23pc)+'_Ninputs1'+location_l23pc+'_Econ'+str(Econ_l23pc)+'_wNMDA'+str(wNMDA_l23pc)+'_Npulses4_ISI'+str(ISIdt)+'_withV.mat'
    Ca_rates0 = minimum(MAT['currClips'],0)*frac_Ca*Ca_input_coeff*1e-9/(-2*1.602e-19)/1e3 #Current (nA) -> Current (nA) carried by Ca -> Current (A) carried by Ca -> Ca ions/sec -> Ca ions/msec
    Ca_Ts = range(-100,-50,4)+range(-50,-40,2)+range(-40,90)+range(90,100,2)+range(100,120,4)+range(120,160,8)+range(160,240,16)+range(240,496,32)+range(496,996,100)+[900]
    Ca_rates = []
    Ca_times = []
    Vspine = []

    for i in range(0,len(Ca_Ts)-1):
      Ca_rates.append(mean(Ca_rates0[100+Ca_Ts[i]:100+Ca_Ts[i+1]])) #+mean(Ca_rates0[2000+Ca_Ts[i]:2000+Ca_Ts[i+1]]))) #earlier it was taken as a mean across two curves, one at 1000ms and another at 2000ms
      Vspine.append(mean(MAT['VspineClips'][100+Ca_Ts[i]:100+Ca_Ts[i+1]]))
      Ca_times.append(mean(range(Ca_Ts[i],Ca_Ts[i+1])))
    #axs[0,iISI].plot(Ca_times, Ca_rates, 'k-', label = 'ISI = '+str(ISIdts[iISI]+30)+' ms',lw=1.0)
    axs[0,iISI].plot(A['DATA'][0][istart:iend]-4159000,A['DATA'][1][istart:iend]*1e6,'k-', label = 'ISI = '+str(ISIdts[iISI]+30)+' ms',lw=1.0)

    axs[3,iISI].plot(range(-100,900),MAT['VspineClips'],'k-',lw=1.0)
    axnew[1+iISI].plot(range(-100,900),1/(1+(1.0/4.1)*exp(0.063*(-MAT['VspineClips']))),'r-',lw=0.5)
    #axs[3,iISI].plot(Ca_times,Vspine,'k-',lw=1.0)
    print('ISI = '+str(ISIdts[iISI]+30)+' ms')

#f.savefig('fig6withvspine_icell'+str(icell_l23pc)+'_Econ'+str(Econ_l23pc)+'_wNMDA'+str(wNMDA_l23pc)+'_contnm_flux'+fluxtxt+'_pulseamp'+str(pulseamp_l23pc)+'_Nsyn'+str(Nsyn_l23pc)+'.eps')
f.savefig('fig5.eps')

########################## Draw the LTP/LTD windows w.r.t. pairingISI ################################
Lfluxes = [0.0, Lfluxmax, 0.0, Lfluxmax]
AChfluxes = [0.0, 0.0, AChfluxmax, AChfluxmax]
labels = ['No $\\beta$-adr. ligand','With $\\beta$-adr. ligand']
for ifile in [1,0,3,2]:
  Lflux = Lfluxes[ifile]
  AChflux = AChfluxes[ifile]
  LTPcoeffs_saved = []
  LTDcoeffs_saved = []
  gLTPLTDcoeffs_saved = []
  ISIs_saved = []
  maxLTP = 1.0
  maxLTD = 1.0
  for iISI in range(0,len(TRAINISISALL)):
    ISI = TRAINISISALL[iISI]
    MAT = {}

    filename_nrn = 'nrn_tstop5000000_tol1e-06_onset4040000.0_n120_freq1.0_dur3.0_CaCoeff'+str(CaCoeff)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(AChflux)+'_contnm3560000.0_600000.0_Ntrains1_trainT100000.0_pair'+str(ISI)+'_icell'+str(icell_l23pc)+'_pulseamp'+str(pulseamp_l23pc)+'_Nsyn'+str(Nsyn_l23pc)+'_Ninputs1_Econ'+str(Econ_l23pc)+'_wNMDA'+str(wNMDA_l23pc)+'_'+location_l23pc+'.mat'

    DATANRN_all = {}
    if not exists(filename_nrn):
      print filename_nrn+' does not exist'
      continue
    foundone = 1
    print "Found "+filename_nrn
    DATANRN_all_all = scipy.io.loadmat(filename_nrn)
    for ikey in range(0,len(DATANRN_all_all['headers'])):
      mykey = DATANRN_all_all['headers'][ikey][0:DATANRN_all_all['headers'][ikey].find(' ')]
      DATANRN_all[mykey] = DATANRN_all_all['DATA'][ikey]

    if len(MAT) > 0:
      times = [T/(DATA_all.shape[0]-1)*i for i in range(0,DATA_all.shape[0])]
    if len(DATANRN_all) > 0:
      times_nrn = DATANRN_all['tvec']

    ISIs_saved.append(ISI+30) #Measure ISI from onset of presyn activation until the onset of the _last_ post-synaptic stimulus (within the burst), not the the first one as before
    TCs_all = []
    TCsN_all = []
    for iax in range(0,len(species)):
      for ispecgroup in range(0,len(species[iax])):
        specgroup = species[iax][ispecgroup]
        if len(MAT) > 0:
          mytimecourse = zeros(DATA_all[:,0].shape[0])
        if len(DATANRN_all) > 0:
          mytimecourse_nrn = zeros(times_nrn.shape[0])
        if type(specgroup) is not list:
          specgroup = [specgroup]
        for ispec in range(0,len(specgroup)):
          specfactor = 1.0
          if len(MAT) > 0:
            mytimecourse = mytimecourse + specfactor*DATA_all[:,inddict[specgroup[ispec]]]
          if len(specgroup[ispec]) > 24 and len(DATANRN_all) > 0:
            mytimecourse_nrn = mytimecourse_nrn + DATANRN_all[specgroup[ispec][:24]]
          elif len(DATANRN_all) > 0:
            mytimecourse_nrn = mytimecourse_nrn + DATANRN_all[specgroup[ispec]]
        if iax == 3:
          LTPcoeffs_saved.append(mytimecourse_nrn[-1]/mytimecourse_nrn[0])
        elif iax == 5:
          LTDcoeffs_saved.append(mytimecourse_nrn[-1]/mytimecourse_nrn[0])

        factor = 1.0/6.022e23/my_volume*1e9
        nrnfactor = 1.0
        TCs_all.append(mytimecourse_nrn[::Nskip]*1e6*nrnfactor)
        TCsN_all.append(mytimecourse_nrn[::Nskip]*1e6*nrnfactor/factor)

    ENhom1_np = (TCsN_all[3] + TCsN_all[5])/4.0 * (TCsN_all[3]-TCsN_all[1])**4/(TCsN_all[3] + TCsN_all[5])**4                                 #Number of complexes times the probability of a complex consisting of 4 non-phos GluR1s
    ENhom1_p = (TCsN_all[3] + TCsN_all[5])/4.0 * (TCsN_all[3]**4 - (TCsN_all[3]-TCsN_all[1])**4)/(TCsN_all[3] + TCsN_all[5])**4               #The same, but use prob. of having 4 GluR1s, minus the cases where all of them are non-phos
    ENhom2 = (TCsN_all[3] + TCsN_all[5])/4.0 * (TCsN_all[5]/(TCsN_all[3] + TCsN_all[5]))**4
    ENhet = (TCsN_all[3] + TCsN_all[5])/4.0 * (1 - (TCsN_all[3]/(TCsN_all[3] + TCsN_all[5]))**4 - (TCsN_all[5]/(TCsN_all[3] + TCsN_all[5]))**4)
    Egtot = ENhom1_np*conds_hom1[0] + ENhom1_p*conds_hom1[1] + ENhom2*conds_hom2 + ENhet*conds_het

    gLTPLTDcoeffs_saved.append(Egtot[-1]/Egtot[0])

  if icols[ifile] == 0:
    for iax in range(0,3):
      axs[irows[ifile],iax].plot([min(ISIs_saved),max(ISIs_saved)],[1.0,1.0],'k--',dashes=[4,2],lw=0.5)
  #axs[irows[ifile],0].plot(ISIs_saved,LTPcoeffs_saved, 'k'+styles[ifile], markersize=2.5, mew=0.8, color=cols[icols[ifile]],label=labels[icols[ifile]],dashes=dashes_3[idashes[ifile]],lw=1.0)
  #axs[irows[ifile],1].plot(ISIs_saved,LTDcoeffs_saved, 'k'+styles[ifile], markersize=2.5, mew=0.8, color=cols[icols[ifile]],label=labels[icols[ifile]],dashes=dashes_3[idashes[ifile]],lw=1.0)
  #axs[irows[ifile],2].plot(ISIs_saved,gLTPLTDcoeffs_saved, 'k'+styles[ifile], markersize=2.5, mew=0.8, color=cols[icols[ifile]],label=labels[icols[ifile]],dashes=dashes_3[idashes[ifile]],lw=1.0)
  axs[irows[ifile],0].plot(ISIs_saved,LTPcoeffs_saved, 'k'+styles[ifile], markersize=2.5, mew=mews[ifile], color=cols[icols[ifile]],label=labels[icols[ifile]],lw=0.1)
  axs[irows[ifile],1].plot(ISIs_saved,LTDcoeffs_saved, 'k'+styles[ifile], markersize=2.5, mew=mews[ifile], color=cols[icols[ifile]],label=labels[icols[ifile]],lw=0.1)
  axs[irows[ifile],2].plot(ISIs_saved,gLTPLTDcoeffs_saved, 'k'+styles[ifile], markersize=2.5, mew=mews[ifile], color=cols[icols[ifile]],label=labels[icols[ifile]],lw=0.1)
  maxLTP = max(maxLTP,max(LTPcoeffs_saved))
  maxLTD = min(maxLTD,min(LTDcoeffs_saved))

  axs[irows[ifile],0].set_ylim([0.9,2.0])
  axs[irows[ifile],0].set_yticks([1.0,1.4,1.8])
  axs[irows[ifile],1].set_ylim([0,1.1])
  axs[irows[ifile],1].set_yticks([0.2,0.6,1.0])

axs[1,2].set_ylim([0.7,1.9])
axs[1,2].set_yticks([1.0,1.4])
axs[2,2].set_ylim([0.55,1.75])
axs[2,2].set_yticks([0.6,1.0,1.4])
print "ifile = "+str(ifile)+", max gcoeff = "+str(maxLTP)+", min gcoeff = "+str(maxLTD)

for iay in range(0,4):
  for iax in range(0,3):
    for tick in axs[iay,iax].xaxis.get_major_ticks() + axs[iay,iax].yaxis.get_major_ticks():
      tick.label.set_fontsize(6)
    for axis in ['top','bottom','left','right']:
      axs[iay,iax].spines[axis].set_linewidth(0.5)
    if iay == 0:    
      for tick in axnew[1+iax].xaxis.get_major_ticks() + axnew[1+iax].yaxis.get_major_ticks():
        tick.label.set_fontsize(5)
      for axis in ['top','bottom','left','right']:
        axnew[1+iax].spines[axis].set_linewidth(0.5)

for iax in range(0,len(axarr)):
  pos = axarr[iax].get_position()
  if iax < 3:
    f.text(pos.x0 + 0.01, pos.y1 - 0.04, chr(ord('A')+iax+1+3), fontsize=12)
    #f.text(pos.x0 + 0.18, pos.y1 - 0.02, 'ISI = '+str(ISIdts[iax]+30)+' ms', fontsize=6, ha='right')
  elif iax >= 9:
    f.text(pos.x0 + 0.01, pos.y1 - 0.02, chr(ord('A')+iax-9+1), fontsize=12)
    f.text(pos.x0 + 0.19, pos.y1 - 0.01, 'ISI = '+str(ISIdts[iax-9]+30)+' ms', fontsize=6, ha='right')
    axarr[iax].set_xticklabels([])
  else:
    f.text(pos.x0 - 0.08 - 0.03*(iax%3==0), pos.y1 - 0.04, chr(ord('A')+iax+1+3), fontsize=12)
pos = axnew[0].get_position()
f.text(pos.x0, pos.y1 + 0.02, 'A', fontsize=12)

Ca_titles = ['Ca$^{2+}$ influx, ISI = '+str(int(ISIdts[iISI]+30))+' ms' for iISI in range(0,3)]
titles = ['GluR1 at memb.','GluR2 at memb.','total syn. strength']
for iax in range(0,3):
  #axs[0,iax].set_xlim([-85,140])
  #axs[0,iax].set_xticks([-50, 0, 50, 100])
  #axs[0,iax].set_ylim([0,3500])
  axs[0,iax].set_xlabel('time (ms)', fontsize = 7)
  axs[0,iax].set_ylabel('[Ca$^{2+}$] (nM)', fontsize = 7)
  axs[0,iax].set_ylim([0,1000])
  axs[0,iax].set_yticks([0,500,1000])
  #axs[0,iax].set_xlim([-200,800])
  #axs[0,iax].set_xticks([0,400,800])
  axs[0,iax].set_xlim([-85,140])
  axs[0,iax].set_xticks([-50,0,50,100])

  axs[3,iax].set_xlim([-85,140])
  axs[3,iax].set_xticks([-50, 0, 50, 100])
  axs[3,iax].set_ylabel('$V_{\mathrm{spine}}$ (mV)', fontsize = 7)
  axs[3,iax].set_ylim([-80,-15])
  axnew[1+iax].set_xlim([-80,140])
  axnew[1+iax].set_xticks([]) #([-50, 150])
  axnew[1+iax].set_ylim([0,0.5])
  axnew[1+iax].set_yticks([0,0.25,0.5])
  axnew[1+iax].spines['bottom'].set_color('#FF0000')
  axnew[1+iax].spines['left'].set_color('#FF0000')
  axnew[1+iax].tick_params(axis='y',color='#FF0000',length=1.6)
  [t.set_color('#FF0000') for t in axnew[1+iax].yaxis.get_ticklabels()]
  if iax > 0:
    axnew[1+iax].set_yticklabels([])
    for iay in [0,3]:
      axs[iay,iax].set_ylabel('', fontsize = 8)
      axs[iay,iax].set_yticklabels([])
    
  axs[1,iax].set_title(titles[iax], fontsize = 8)
  axs[2,iax].set_xlabel('ISI (ms)', fontsize = 7)
  for iay in range(1,2):
    axs[iay,iax].set_xticklabels([])

for iax in range(0,3):
  for iay in [1,2]:
    axs[iay,iax].set_xlim([-170.0,230.0])
    axs[iay,iax].set_xticks([-150.0,0,150.0])

axs[1,2].legend(fontsize=6,frameon=False)
for ax in axarr+axnew:
  for line in ax.xaxis.get_ticklines()+ax.yaxis.get_ticklines():
    line.set_markeredgewidth(0.5)


conditions = ['(No ACh)','(With ACh)']
for iay in range(0,2):
  axs[1+iay,0].set_ylabel(conditions[iay]+'\nGluR1',fontsize=7)
  axs[1+iay,1].set_ylabel('GluR2',fontsize=7)
  axs[1+iay,2].set_ylabel('$g_{\mathrm{rel}}$',fontsize=7)

f.savefig('fig5.eps')


