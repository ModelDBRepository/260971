#drawfig11.py: Draws the figure of model calibration
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
rc('axes',linewidth=0.5)


f,axarr = subplots(12,1)
for iax in range(0,12):
  axarr[iax].spines['top'].set_visible(False)
  axarr[iax].spines['right'].set_visible(False)
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.5)
  for line in axarr[iax].xaxis.get_ticklines()+axarr[iax].yaxis.get_ticklines():
    line.set_markeredgewidth(0.5)

axnew = []
axarr[0].set_position([0.08,0.8,0.2,0.16])
axarr[1].set_position([0.36,0.8,0.2,0.16])
axarr[2].set_position([0.64,0.8,0.3,0.16]); axnew.append(f.add_axes([0.76, 0.8+0.07, 0.16, 0.09]))

axarr[3].set_position([0.08,0.56,0.16,0.16]) 
axarr[4].set_position([0.32,0.56,0.16,0.16]) 
axarr[5].set_position([0.56,0.56,0.16,0.16]) 
axarr[6].set_position([0.80,0.56,0.16,0.16]) 

axarr[7].set_position([0.08,0.32,0.4,0.16]); axnew.append(f.add_axes([0.14, 0.32+0.07, 0.085, 0.1]))
axarr[8].set_position([0.56,0.32,0.16,0.16])
axarr[9].set_position([0.8,0.32,0.16,0.16])

axarr[10].set_position([0.08,0.08,0.3,0.16])
axarr[11].set_position([0.46,0.08,0.48,0.16]);

axarr[10].set_visible(False)
axarr[11].set_visible(False)

HoffmanData = array([[0.024525077566201113, 0.010894660894660957],
 [0.05026882121429447, 0.007422147744728269],
 [0.10150641893704561, 0.014699995345156758],
 [0.2526954386274699, 0.014841968067774713],
 [0.49522800371345205, 0.02211516082483822],
 [0.970538989520764, 0.03476469766792345],
 [2.019271916170997, 0.05996834706512111],
 [2.9786945720946925, 0.08691058045896738],
 [3.8407459778154385, 0.1174160964483546],
 [4.806380863064389, 0.15508541637573892],
 [10, 0.3863822557370946],
 [14.751329666105498, 0.5943280733603316],
 [19.306977288832496, 0.7269864544058093],
 [38.98603702549072, 0.9009309686729041],
 [57.50958845380086, 0.9529628078015175],
 [98.5159372650316, 0.9888888888888889],
 [145.32410679618488, 0.9997020900246707],
 [193.06977288832496, 0.9907857375599312],
 [407.7459071436882, 0.9962784527300657]]) #Obtained with https://apps.automeris.io/wpd/

for iax in range(0,len(axnew)):
  axnew[iax].spines['top'].set_visible(False)
  axnew[iax].spines['right'].set_visible(False)
  for tick in axnew[iax].xaxis.get_major_ticks() + axnew[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  for axis in ['top','bottom','left','right']:
    axnew[iax].spines[axis].set_linewidth(0.25)
  for line in axnew[iax].xaxis.get_ticklines()+axnew[iax].yaxis.get_ticklines():
    line.set_markeredgewidth(0.25)

mesh_input_file = open('mesh_general.out','r')
mesh_firstline = mesh_input_file.readline()
mesh_secondline = mesh_input_file.readline()
mesh_values = mesh_secondline.split()
my_volume = float(mesh_values[-2])*1e-15 #litres
Nskip = 1

def mydraw(ax,species,DATAsets,mydashes=[],colors=[]):
  DATA = []
  tcs = []
  for idata in range(0,len(DATAsets)):    
    DATA.append({})
    for ikey in range(0,len(DATAsets[idata]['headers'])):
      mykey = DATAsets[idata]['headers'][ikey][0:DATAsets[idata]['headers'][ikey].find(' ')]
      DATA[idata][mykey] = DATAsets[idata]['DATA'][ikey]

    times_nrn = DATA[idata]['tvec']
    mytimecourse_nrn = zeros(times_nrn.shape[0])
    spectext = ''
    if type(species) is not list:
      species = [species]
    for ispec in range(0,len(species)):
      specfactor = 1.0
      if len(species[ispec]) > 24:
        mytimecourse_nrn = mytimecourse_nrn + DATA[idata][species[ispec][:24]]
      else:
        mytimecourse_nrn = mytimecourse_nrn + DATA[idata][species[ispec]]
      spectext = spectext+species[ispec]+'+'
      if len(spectext)>30:
        spectext = spectext[0:-1]+'+\n'+'+'
  
    spectext = spectext[0:-1]
    nrnfactor = 1.0
    mycolor = '#000000'
    if len(colors) > idata and type(colors[idata]) is str:
      mycolor = colors[idata]
    if len(colors) > idata and type(colors[idata]) is type(None):
      print "  not plotting idata = "+str(idata)
    else:
      if len(mydashes) > idata and len(mydashes[idata]) > 0:
        ax.plot([x/1000-4000 for x in times_nrn[::Nskip]],mytimecourse_nrn[::Nskip]*1e6*nrnfactor,'k--',dashes=mydashes[idata],color=mycolor,linewidth=1.0)
      else:
        ax.plot([x/1000-4000 for x in times_nrn[::Nskip]],mytimecourse_nrn[::Nskip]*1e6*nrnfactor,'k-',color=mycolor,linewidth=1.0)
    tcs.append(mytimecourse_nrn[::Nskip]*1e6*nrnfactor)
  return tcs

def mydrawmaxes(ax,species,parVals,DATAsetsets,mydashes=[],colors=[]):
  for iset in range(0,len(DATAsetsets)):
    DATAsets = DATAsetsets[iset] 
    maxVals = []
    DATA = []
    nrnfactor = 1.0
    for idata in range(0,len(DATAsets)):    
      DATA.append({})
      for ikey in range(0,len(DATAsets[idata]['headers'])):
        mykey = DATAsets[idata]['headers'][ikey][0:DATAsets[idata]['headers'][ikey].find(' ')]
        DATA[idata][mykey] = DATAsets[idata]['DATA'][ikey]
      mytimecourse_nrn = zeros(DATA[idata]['tvec'].shape[0])
      if type(species) is not list:
        species = [species]
      for ispec in range(0,len(species)):
        specfactor = 1.0
        if len(species[ispec]) > 24:
          mytimecourse_nrn = mytimecourse_nrn + DATA[idata][species[ispec][:24]]
        else:
          mytimecourse_nrn = mytimecourse_nrn + DATA[idata][species[ispec]]
      maxVals.append(max(mytimecourse_nrn)*1e6*nrnfactor)
  
    mycolor = '#000000'
    if len(colors) > iset and type(colors[iset]) is str:
      mycolor = colors[iset]
    if len(mydashes) > iset and len(mydashes[iset]) > 0:
      ax.plot(parVals,maxVals,'k--',dashes=mydashes[iset],color=mycolor,linewidth=1.0)
    else:
      ax.plot(parVals,maxVals,'k-',color=mycolor,linewidth=1.0)


filename_nrn = 'nrn_tstop5000000_tol1e-06_onset4040000.0_n100_freq100.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains4_trainT4000.0.mat'
assert exists(filename_nrn), filename_nrn+' does not exist'
DATANRN_orig_HFS = scipy.io.loadmat(filename_nrn)
print "loaded "+filename_nrn

########## A ############ (New vs. old membrane-insertion of non-S880-phosphorylated GluR2)
if True:
  filename_nrn = 'nrn_tstop5000000_tol1e-06_Cax1.0_k385-387-389x22.4-22.4-22.4_onset4040000.0_n100_freq100.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains4_trainT4000.0.mat'
  assert exists(filename_nrn), filename_nrn+' does not exist'
  DATANRN_HFS = scipy.io.loadmat(filename_nrn)
  print "loaded "+filename_nrn

  mydraw(axarr[0],['GluR2_memb', 'GluR2_memb_PKCt', 'GluR2_memb_PKCp', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'],[DATANRN_orig_HFS,DATANRN_HFS],[[],[2,2]],['#000000','#808080'])

########## B ############ (old vs. new CaM Ca-binding, CaMCa4 time courses)
if True:
  filename_nrn = 'nrn_oldCaM_tstop5000000_tol1e-06_CaMx0.55_onset4040000.0_n100_freq100.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains4_trainT4000.0.mat'
  assert exists(filename_nrn), filename_nrn+' does not exist'
  DATANRN_HFS = scipy.io.loadmat(filename_nrn)
  print "loaded "+filename_nrn

  tcs_B = mydraw(axarr[2],['AC1GsaGTPCaMCa4', 'AC1GsaGTPCaMCa4ATP', 'AC1GiaGTPCaMCa4', 'AC1GiaGTPCaMCa4ATP', 'AC1GsaGTPGiaGTPCaMCa4', 'AC1GsGiCaMCa4ATP', 'AC1CaMCa4', 'AC1CaMCa4ATP', 'AC8CaMCa4', 'AC8CaMCa4ATP', 'PDE1CaMCa4', 'PDE1CaMCa4cAMP', 'CaMCa4', 'PP2BCaMCa4', 'CKCaMCa4', 'CKpCaMCa4', 'CKpCaMCa4PP1', 'Ip35PP2BCaMCa4', 'Ip35PP1PP2BCaMCa4', 'PP1PP2BCaMCa4'],[DATANRN_orig_HFS,DATANRN_HFS],[[],[2,2]],['#000000','#808080'])
  maxCa_B = DATANRN_HFS['maxDATA'][1]
  maxCa_orig_B = DATANRN_orig_HFS['maxDATA'][1]
  mydraw(axnew[0],['CaMCa4'],[DATANRN_orig_HFS,DATANRN_HFS],[[],[2,2]],['#000000','#808080'])

########## C ############ (old vs. new CaM Ca-binding, Ca-sensitivity)
if True:
  for iset in range(0,2):
    filename_nrn = 'steadystate_new_li2020_0.0_50.0_0.0_processed.sav' if iset == 0 else 'steadystate_new_oldCaM_li2020_0.0_50.0_0.0_processed.sav'
    mycolor = '#000000' if iset == 0 else '#808080'
    assert exists(filename_nrn), filename_nrn+' does not exist'
    unpicklefile = open(filename_nrn,'r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    timesAll = unpickledlist[0][0]
    DATA_all_all = unpickledlist[0][3]  
    headers = unpickledlist[0][4]
    CaFluxes = unpickledlist[1]
    print "loaded "+filename_nrn

    species = ['AC1GsaGTPCaMCa4', 'AC1GsaGTPCaMCa4ATP', 'AC1GiaGTPCaMCa4', 'AC1GiaGTPCaMCa4ATP', 'AC1GsaGTPGiaGTPCaMCa4', 'AC1GsGiCaMCa4ATP', 'AC1CaMCa4', 'AC1CaMCa4ATP', 'AC8CaMCa4', 'AC8CaMCa4ATP', 'PDE1CaMCa4', 'PDE1CaMCa4cAMP', 'CaMCa4', 'PP2BCaMCa4', 'CKCaMCa4', 'CKpCaMCa4', 'CKpCaMCa4PP1', 'Ip35PP2BCaMCa4', 'Ip35PP1PP2BCaMCa4', 'PP1PP2BCaMCa4']
    ispecies = [[i for i in range(0,len(headers)) if species[j] in headers[i] and (len(headers[i]) == species[j] or headers[i][len(species[j])] == ' ')][0] for j in range(0,len(species))]

    axarr[1].semilogx([DATA_all_all[i][1]*1e6 for i in range(0,len(CaFluxes))], [sum([DATA_all_all[i][ispec]*1e6 for ispec in ispecies]) for i in range(0,len(CaFluxes))],'k-',color=mycolor,linewidth=1.0,mew=0.9,ms=3.0)
    axarr[1].set_xlim([3e1,3e6])
    axarr[1].set_yticks([0,20000,40000])
  axarr[1].plot(HoffmanData[:,0]*1e3, HoffmanData[:,1]*50000, 'r.',ms=1.4,mew=1.0,fillstyle='full')

  #axarr[1].plot(1e6*maxCa_B, max(tcs_B[0]), 'k.')
  #axarr[1].plot(1e6*maxCa_orig_B, max(tcs_B[1]), 'k.', color='#808080')
  #axarr[1].set_xlim([0,10000])


########## D ############ (fit the PKCp activation) 
if True:
  filenames_nrn = ['nrn_tstop5000000_tol1e-06_Cax1.0_k411x1.0_onset4040000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains1_trainT100000.0.mat', #fitted value
                   'nrn_tstop5000000_tol1e-06_Cax1.0_k411x0.1_onset4040000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains1_trainT100000.0.mat', #slower PKCp activation (close to original)
                   'nrn_tstop5000000_tol1e-06_Cax1.0_k411x0.3_onset4040000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains1_trainT100000.0.mat', #slower PKCp activation
                   'nrn_tstop5000000_tol1e-06_Cax1.0_k411x3.0_onset4040000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains1_trainT100000.0.mat', #faster PKCp activation
                   'nrn_tstop5000000_tol1e-06_Cax1.0_k411x10.0_onset4040000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains1_trainT100000.0.mat']#faster PKCp activation
  DATANRNs = []
  for ifile in range(0,len(filenames_nrn)):
    filename_nrn = filenames_nrn[ifile]
    assert exists(filename_nrn), filename_nrn+' does not exist'
    DATANRNs.append(scipy.io.loadmat(filename_nrn))
    print "loaded "+filename_nrn
  GluR2_S880_tcs = mydraw(axarr[3],['GluR2_S880','GluR2_S880_PP2A','GluR2_memb_S880','GluR2_memb_S880_PP2A'],DATANRNs,[[],[3,1],[],[],[2,2]],[None,None,None,None,None]) # (just load, don't draw)
  xs = [3,1,2,4,5]
  kfvals = [0.0005*x for x in [1.0,0.1,0.3,3.0,10.0]]
  for i in range(0,5):
    axarr[3].bar(xs[i],GluR2_S880_tcs[i][-1]/270.,edgecolor='#808080' if i != 0 else '#000000',facecolor='#ffffff') ############# D ################# (amount of GluR2 S880 after LFS (270 nM is the total concentration of GluR2))
    axarr[3].text(xs[i]-0.07,0.04+GluR2_S880_tcs[i*(i<3)][-1]/270,'k$_f$ = '+str(kfvals[i]),fontsize=6,rotation='vertical',va='bottom')
  axarr[3].set_xticks([])
  axarr[3].set_ylim([0,1.6])

########## E -- G ############
  filenames_nrn = ['nrn_tstop5000000_tol1e-06_Cax1.0_k411x1.0_onset4040000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains1_trainT100000.0.mat', #fitted value
                   'nrn_tstop5000000_tol1e-06_Cax1.0_k411x0.3_onset4040000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains1_trainT100000.0.mat', #slower PKCp activation
                   'nrn_tstop5000000_tol1e-06_Cax1.0_k411x3.0_onset4040000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains1_trainT100000.0.mat'] #faster PKCp activation
  DATANRNs = []
  for ifile in range(0,len(filenames_nrn)):
    filename_nrn = filenames_nrn[ifile]
    assert exists(filename_nrn), filename_nrn+' does not exist'
    DATANRNs.append(scipy.io.loadmat(filename_nrn))
    print "loaded "+filename_nrn
  mydraw(axarr[4],['PKCt'],DATANRNs,[[],[3,1],[1,3]],['#000000','#808080','#808080']) ############# E ################# (PKCt)
  mydraw(axarr[5],['PKCp'],DATANRNs,[[],[3,1],[1,3]],['#000000','#808080','#808080']) ############# F ################# (PKCp)
  GluR2_S880_tcs = mydraw(axarr[6],['GluR2_S880','GluR2_S880_PP2A','GluR2_memb_S880','GluR2_memb_S880_PP2A'],DATANRNs,[[],[3,1],[1,3]],['#000000','#808080','#808080']) ############# G ################# (GluR2 S880 time course)

f.savefig('fig11.eps')
########## H ############
if True:
  A = scipy.io.loadmat('cAMP_withdiss_test_tstop22000_tol1e-08_onset800.0_n1_freq1.0_dur16000.0_flux0.64.mat')
  if A['times'].shape[1] == 1:
    A['times'] = A['times'].T
  OBJECTIVE_ts = [1000+1000*i for i in range(0,18)]
  OBJECTIVE_DATA = [mytools.interpolate(A['times'][0],A['tcDATA'][0]+A['tcDATA'][1],OBJECTIVE_ts), #Time course of cAMP                                                                                              
                    mytools.interpolate(A['times'][0],A['tcDATA'][5],OBJECTIVE_ts)] #Time course of PKAc                                                                                         
  ks = [0.4e9, 1.0e9, 1.6e9, 2.2e9, 2.8e9]
  fs = []

  polygon = Polygon(array([[1000,17000,17000,1000],[0,0,100,100]]).T, True)
  p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
  p.set_facecolor('#E0E0E0')
  p.set_edgecolor('none')
  axarr[7].add_collection(p)

  mydashes = [[2,2],[1,3],[1,1],[1,3,2,2],[3,1]]
  for ik in range(0,len(ks)):
    unpicklefile = open('run_cAMP_withdiss_flux0.64_k'+str(ks[ik])+'_data.sav','r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()

    timesAll = unpickledlist[0]
    timeCoursesAll = unpickledlist[1]
    mydict = {}
    objs = []
    imeas = [0, 5] # cAMP and PKAc traces
    for iobjective in range(0,len(OBJECTIVE_DATA)):
      objs_thisobj = []
      times = timesAll[0]
      if times.shape[1] == 1:
        times = times.T
      timeCourses = mytools.interpolate(times[0],timeCoursesAll[0][imeas[iobjective]],OBJECTIVE_ts)
      mydict['f'+str(iobjective)] = sum(abs(array(timeCourses)-array(OBJECTIVE_DATA[iobjective])))
    fs.append(mydict['f1'])
    if ik != 2:
      #axarr[7].plot(times[0],1e6*timeCoursesAll[0][imeas[1]],'b--',color='#808080',dashes=mydashes[ik])
      axarr[7].plot(times[0],1e6*timeCoursesAll[0][imeas[1]],'b-',color='#A0A0A0',linewidth=1.0)
      axnew[1].bar(ik,1e6*mydict['f1'],edgecolor='#A0A0A0',facecolor='#A0A0A0')
      print str(ik)+","+str(1e6*mydict['f1'])
    else:
      axarr[7].plot(times[0],1e6*timeCoursesAll[0][imeas[1]],'k-',color='#000000',linewidth=1.0)
      control_data = [times[0],1e6*timeCoursesAll[0][imeas[1]]]
      axnew[1].bar(ik,1e6*mydict['f1'],edgecolor='#000000',facecolor='#000000')
      print str(ik)+","+str(1e6*mydict['f1'])
    axnew[1].text(ik-0.3,(1e6*mydict['f1']*(ik > 0))+40,'k$_f$ = '+"{:.1E}".format(ks[ik]),fontsize=6,rotation='vertical',va='bottom')

  axarr[7].plot(OBJECTIVE_ts,[x*1e6 for x in OBJECTIVE_DATA[1]],'kx',color='#808080')
  axarr[7].plot(control_data[0],control_data[1],'k-',color='#000000',linewidth=1.0)
  axarr[7].set_xlim([-2000,18500])
  axarr[7].set_ylim([0,100])
  axnew[1].set_xticks([])
  axnew[1].set_ylim([0,800])

########## I ############ (old PKA, 4xHFS with Lflux 10.0) (Change to response curve to cAMP in the stub model?)
if True:
  filename_nrn = 'nrn_oldPKA_tstop5000000_tol1e-06_onset4040000.0_n100_freq100.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains4_trainT4000.0.mat'
  assert exists(filename_nrn), filename_nrn+' does not exist'
  DATANRN_HFS = scipy.io.loadmat(filename_nrn)
  print "loaded "+filename_nrn

  mydraw(axarr[8],['GluR1_S845', 'GluR1_S845_S831', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_S845_CKp', 'GluR1_S845_PKCt', 'GluR1_S845_PKCp', 'GluR1_S845_PP1', 'GluR1_S845_S831_PP1', 'GluR1_S845_PP2B', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S845', 'GluR1_memb_S845_S831', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCt', 'GluR1_memb_S845_PKCp', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'],[DATANRN_orig_HFS,DATANRN_HFS],[[],[2,2]],['#000000','#808080'])
  
########## J ############ (PKC does not phosphorylate S831)
if True:
  filename_nrn = 'nrn_tstop5000000_tol1e-06_Cax1.0_k166-169-181-184-218-221-233-236x0.0-0.0-0.0-0.0-0.0-0.0-0.0-0.0_onset4040000.0_n100_freq100.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains4_trainT4000.0.mat'
  assert exists(filename_nrn), filename_nrn+' does not exist'
  DATANRN_HFS = scipy.io.loadmat(filename_nrn)
  print "loaded "+filename_nrn

  mydraw(axarr[9],['GluR1_S831', 'GluR1_S845_S831', 'GluR1_S831_PKAc', 'GluR1_S845_S831_PP1', 'GluR1_S831_PP1', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'],[DATANRN_orig_HFS,DATANRN_HFS],[[],[2,2]],['#000000','#808080'])
  

for iax in range(0,10):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.06, pos.y1 - 0.02, chr(ord('A')+iax), fontsize=15)

for ax in [axarr[0], axarr[2], axnew[0], axarr[1], axarr[4], axarr[5], axarr[6], axarr[7], axnew[1], axarr[8], axarr[9], axarr[10], axarr[11]]:
  ylab = ax.set_ylabel('(nM)',fontsize = 6)
  if ax == axnew[0]:
    xlab = ax.set_xlabel('time (s)',fontsize = 6)
    ax.xaxis.set_label_coords(0.5, -0.2)
  elif ax == axarr[1]:
    xlab = ax.set_xlabel('flux (part./ms)',fontsize = 6)
  elif ax != axnew[1]:
    xlab = ax.set_xlabel('time (s)',fontsize = 6)

axarr[0].text(mean(axarr[0].get_xlim()), [0.01*x[0]+0.99*x[1] for x in [axarr[0].get_ylim()]][0],'[GluR2 at memb.]', fontsize = 6, ha = 'center')
axarr[1].text(mean(axarr[1].get_xlim()), [0.01*x[0]+0.99*x[1] for x in [axarr[1].get_ylim()]][0],'max. [CaMCa4]', fontsize = 6, ha = 'center')
axarr[2].text(mean(axarr[2].get_xlim()), [-0.04*x[0]+1.04*x[1] for x in [axarr[2].get_ylim()]][0],'[CaMCa4]', fontsize = 6, ha = 'center')
axarr[3].text(mean(axarr[3].get_xlim()), [0.01*x[0]+0.99*x[1] for x in [axarr[3].get_ylim()]][0],'[GluR2 S880p]/[tot. GluR2]', fontsize = 6, ha = 'center')
axarr[4].text(mean(axarr[4].get_xlim()), [0.01*x[0]+0.99*x[1] for x in [axarr[4].get_ylim()]][0],'[PKC tr.]', fontsize = 6, ha = 'center')
axarr[5].text(mean(axarr[5].get_xlim()), [0.01*x[0]+0.99*x[1] for x in [axarr[5].get_ylim()]][0],'[PKC pers.]', fontsize = 6, ha = 'center')
axarr[6].text(mean(axarr[6].get_xlim()), [0.01*x[0]+0.99*x[1] for x in [axarr[6].get_ylim()]][0],'[GluR2 S880p]', fontsize = 6, ha = 'center')
axarr[7].text(mean(axarr[7].get_xlim()), [0.01*x[0]+0.99*x[1] for x in [axarr[7].get_ylim()]][0],'[act. PKA]', fontsize = 6, ha = 'center')
axarr[8].text(mean(axarr[8].get_xlim()), [0.01*x[0]+0.99*x[1] for x in [axarr[8].get_ylim()]][0],'[GluR1 S845p]', fontsize = 6, ha = 'center')
axarr[9].text(mean(axarr[9].get_xlim()), [0.01*x[0]+0.99*x[1] for x in [axarr[9].get_ylim()]][0],'[GluR1 S831p]', fontsize = 6, ha = 'center')

f.savefig('fig11.eps')
