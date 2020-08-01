#drawfig3.py: Draws the figure of steady-state activation by Ca.
#Tuomo Maki-Marttnen, 2019-2020
#CC BY 4.0
import matplotlib
matplotlib.use('Agg')
from pylab import *
import pickle
from os.path import exists
from matplotlib.collections import PatchCollection
import scipy.io

#Check whether the processed data already exists. If not, load them from the huuge *_raw.sav files saved by simsteadystates.py.
for flux in [0.0, 0.005, 0.05]:
  if not exists('steadystate_flux'+str(flux)+'.sav'):
    print 'loading steadystate_flux'+str(flux)+'_raw.sav'
    unpicklefile = open('steadystate_flux'+str(flux)+'_raw.sav', 'r')
    unpickledlist = pickle.load(unpicklefile)
    unpicklefile.close()
    timesThis = unpickledlist[0]
    condsThis = unpickledlist[1]
    maxCasThis = unpickledlist[2]
    DATA_all_all_all = unpickledlist[3]
    Ca_input_fluxes  = [2.5*i for i in range(0,101)] + [250+5*i for i in range(0,101)] + [750+10*i for i in range(0,101)] + [1750+20*i for i in range(0,101)]
    try:
      conds = [condsThis[i][-1] for i in range(0,len(condsThis))]
      baselines = [condsThis[i][0] for i in range(0,len(condsThis))]
    except:
      print 'Warning: missing data in steadystate_flux'+str(flux)+'_raw.sav'
      conds = []
      baselines = []
      missings = []
      for i in range(0,len(condsThis)):
        if len(condsThis[i]) == 0:
          missings.append(i)
          Ca_input_fluxes = Ca_input_fluxes[0:i]+Ca_input_fluxes[i+1:]
        else:
          conds.append(condsThis[i][-1])
          baselines.append(condsThis[i][0])

    picklelist = [Ca_input_fluxes, conds, maxCasThis,baselines,DATA_all_all_all]
    file=open('steadystate_flux'+str(flux)+'.sav', 'w')
    pickle.dump(picklelist,file)
    file.close()
    print 'saved steadystate_flux'+str(flux)+'.sav'

flux = 0.05
print 'loading steadystate_flux'+str(flux)+'.sav'
unpicklefile = open('steadystate_flux'+str(flux)+'.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
Ca_input_fluxes_2 = unpickledlist[0]
DATA_all_2 = unpickledlist[4]

print 'loading steadystate_flux0.0.sav'
unpicklefile = open('steadystate_flux0.0.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
Ca_input_fluxes_0 = unpickledlist[0]
DATA_all_0 = unpickledlist[4]

print 'loading steadystate_flux0.005.sav'
unpicklefile = open('steadystate_flux0.005.sav', 'r')
unpickledlist = pickle.load(unpicklefile)
unpicklefile.close()
Ca_input_fluxes_1 = unpickledlist[0]
DATA_all_1 = unpickledlist[4]

Ca_input_fluxes = Ca_input_fluxes_2
DATA_all = DATA_all_2

mesh_input_file = open('mesh_general.out','r')
mesh_firstline = mesh_input_file.readline()
mesh_secondline = mesh_input_file.readline()
mesh_values = mesh_secondline.split()
my_volume = float(mesh_values[-2])*1e-15 #litres
mesh_input_file.close()

#species: 
#  0: buffers
#  1: pumps
#  2: PKC
#  3: CaM
#  4: free Ca

timecourse_iinputs = [100,80,60]
cols = ['#440154', '#470f62', '#481d6f', '#472a79', '#453681', '#414387', '#3c4f8a', '#37598c', '#32648e', '#2d6f8e', '#29788e', '#26828e', '#228b8d', '#1f958b', '#1e9f88', '#22a884', '#2bb17e', '#3bbb75', '#4dc36b', '#62cb5f', '#7ad251'][::10] + ['#dddd00']
nrncols = ['#360043', '#390c4f', '#391759', '#392161', '#372b67', '#34366c', '#303f6f', '#2c4770', '#285071', '#245872', '#216072', '#1e6872', '#1b6f71', '#19776f', '#187f6d', '#1b866a', '#238e65', '#2f955e', '#3e9c56', '#4fa24c', '#61a841'][::10] + ['#cccc00']
rc('axes',linewidth=0.5)
f,axs = subplots(4,4)
axarr = sum([axs[i].tolist() for i in range(0,len(axs))]+[[]])
for iax in range(0,len(axarr)):
  for line in axarr[iax].xaxis.get_ticklines()+axarr[iax].yaxis.get_ticklines():
    line.set_markeredgewidth(0.5)

DATA_nrd = []
headers_nrd =  []
for iflux in range(0,len(timecourse_iinputs)):
  print 'Loading ../NeuroRD/tstop500000_tol0.01_onset40000.0_n1_freq1.0_dur300000.0_flux'+str(Ca_input_fluxes[timecourse_iinputs[iflux]])+'_Lflux0.05_Gluflux0.05_AChflux0.05_Ntrains1_trainT100000.0_8seeds.mat'
  A = scipy.io.loadmat('../NeuroRD/tstop500000_tol0.01_onset40000.0_n1_freq1.0_dur300000.0_flux'+str(Ca_input_fluxes[timecourse_iinputs[iflux]])+'_Lflux0.05_Gluflux0.05_AChflux0.05_Ntrains1_trainT100000.0_8seeds.mat')
  DATA_nrd.append(A['DATA'])
  headers_nrd.append(A['headers'])

timecourse_species_titles = ['buffers', 'pumps', 'PKC pathway', 'CaM', 'free Ca']
timecourse_species = [['fixedbufferCa','CalbinC','fixedbuffer','Calbin'],
                      ['NCXCa','PMCACa','NCX','PMCA'],
                      ['PLC', 'PLCCa', 'PLCCaGqaGTP', 'PLCGqaGTP', 'PLCCaPip2', 'PLCCaGqaGTPPip2', 'PLCCaDAG', 'PLCCaGqaGTPDAG', 'PLA2', 'CaPLA2', 'CaPLA2Pip2'],
                      ['AC1GsaGTPCaMCa4', 'AC1GsaGTPCaMCa4ATP', 'AC1GiaGTPCaMCa4', 'AC1GiaGTPCaMCa4ATP', 'AC1GsaGTPGiaGTPCaMCa4', 'AC1GsGiCaMCa4ATP', 'AC1CaMCa4', 'AC1CaMCa4ATP', 'AC8CaMCa4', 'AC8CaMCa4ATP', 'PDE1CaMCa4', 'PDE1CaMCa4cAMP', 'NgCaM', 'CaM', 'CaMCa2', 'CaMCa3', 'CaMCa4', 'PP2BCaM', 'PP2BCaMCa2', 'PP2BCaMCa3', 'PP2BCaMCa4', 'CKCaMCa4', 'CKpCaMCa4', 'Complex', 'pComplex', 'CKpCaMCa4PP1', 'Ip35PP2BCaMCa4', 'Ip35PP1PP2BCaMCa4', 'PP1PP2BCaMCa4', 'GluR1_CKCam', 'GluR1_CKpCam', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam'],
                      ['Ca']
                     ]

timecourse_coeffs = [[1,1,0,0],
                     [1,1,0,0],
                     [0,1,1,0,1,1,1,1,0,1,1],
                     [4,4,4,4,4,4,4,4,4,4,4,4,0,0,2,3,4,0,2,3,4,4,4,8,8,4,4,4,4,4,4,4,4,4,4,4,4],
                     [1]
                    ]

timecourse_allCa_species = ['Ca', 'CaOutLeak', 'CalbinC', 'PMCACa', 'NCXCa', 'AC1GsaGTPCaMCa4', 'AC1GsaGTPCaMCa4ATP', 'AC1GiaGTPCaMCa4', 'AC1GiaGTPCaMCa4ATP', 'AC1GsaGTPGiaGTPCaMCa4', 'AC1GsGiCaMCa4ATP', 'AC1CaMCa4', 'AC1CaMCa4ATP', 'AC8CaMCa4', 'AC8CaMCa4ATP', 'PDE1CaMCa4', 'PDE1CaMCa4cAMP', 'CaMCa2', 'CaMCa3', 'CaMCa4', 'PP2BCaMCa2', 'PP2BCaMCa3', 'PP2BCaMCa4', 'CKCaMCa4', 'CKpCaMCa4', 'Complex', 'pComplex', 'CKpCaMCa4PP1', 'Ip35PP2BCaMCa4', 'Ip35PP1PP2BCaMCa4', 'PP1PP2BCaMCa4', 'GluR1_CKCam', 'GluR1_CKpCam', 'GluR1_PKCt', 'GluR1_PKCp', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_S845_PKCt', 'GluR1_S845_PKCp', 'GluR1_S845_PP2B', 'GluR1_S845_S831_PP2B', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_PKCt', 'GluR1_memb_PKCp', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_PKCt', 'GluR1_memb_S845_PKCp', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B', 'fixedbufferCa', 'PLCCa', 'PLCCaGqaGTP', 'PLCCaPip2', 'PLCCaGqaGTPPip2', 'PLCCaDAG', 'PLCCaGqaGTPDAG', 'PKCCa', 'PKCt', 'PKCp', 'CaDGL', 'DAGCaDGL', 'GluR2_PKCt', 'GluR2_PKCp', 'GluR2_memb_PKCt', 'GluR2_memb_PKCp', 'CaPLA2', 'CaPLA2Pip2'] #,'CaOut']
timecourse_allCa_coeffs = [1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 3, 4, 2, 3, 4, 4, 4, 8, 8, 4, 4, 4, 4, 4, 4, 1, 1, 4, 4, 1, 1, 4, 4, 4, 4, 1, 1, 4, 4, 1, 1, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

tc_scpoint = 80
tc_xlim = 120
tc_sccoeff = 4.25
conc2parts = 1e-3*my_volume*6.022e23

for iinput in range(0,len(timecourse_iinputs)):
  DATA_this = DATA_all[timecourse_iinputs[iinput]]['DATA']  
  headers = DATA_all[iinput]['headers']
  conc_arr_ref = zeros(len(DATA_this[0,:]))
  itc = argmin(abs(DATA_this[0]-(4040000+tc_scpoint*1000)))
  itend = argmin(abs(DATA_this[0]-(4340000)))
  for ispec in [-1]+range(0,len(timecourse_species)):
    conc_arr = zeros(len(DATA_this[0,:]))
    mytimecourse = timecourse_species[ispec] if ispec >= 0 else timecourse_allCa_species
    mytimecourse_coeffs = timecourse_coeffs[ispec] if ispec >= 0 else timecourse_allCa_coeffs
    conc_arr_nrd = zeros(DATA_nrd[iinput][:,0].shape[0])
    for iispecie in range(0,len(mytimecourse)):
      ispecieind = -1
      ispecieind_nrd = -1
      for iispecieDATA in range(0,len(headers)):
        mystr = headers[iispecieDATA]
        firstspace = mystr.find(' ')
        if firstspace >= 0:
          mystr = mystr[0:firstspace]
        if mystr == mytimecourse[iispecie]:
          ispecieind = iispecieDATA
      for iispecieDATA in range(0,len(headers_nrd[iinput])):
        mystr = headers_nrd[iinput][iispecieDATA]
        firstspace = mystr.find(' ')
        if firstspace >= 0:
          mystr = mystr[0:firstspace]
        if mystr == mytimecourse[iispecie]:
          ispecieind_nrd = iispecieDATA - 4
      if ispecieind == -1 or ispecieind_nrd == -1:
        print mytimecourse[iispecie]+' not found in headers or headers_nrd, ispecieind = '+str(ispecieind)+', ispecieind_nrd = '+str(ispecieind_nrd)
        continue
      conc_arr = conc_arr + DATA_this[ispecieind,:]*mytimecourse_coeffs[iispecie]
      conc_arr_nrd = conc_arr_nrd + DATA_nrd[iinput][:,ispecieind_nrd]*mytimecourse_coeffs[iispecie]
    if ispec == -1:
      conc_arr_ref = conc_arr[:]
      continue
    if len(conc_arr_nrd) < 500001:
      conc_arr_nrd = array(conc_arr_nrd.tolist() + conc_arr_nrd[200000:280000].tolist() + conc_arr_nrd[200000:280000].tolist() + conc_arr_nrd[200000:280000].tolist() + conc_arr_nrd[200000:280000].tolist())[0:500001]
    axarr[1+ispec].plot([(i-40000.)/1000 for i in range(0,500001)],conc_arr_nrd,'k-',color=cols[iinput],lw=0.375)
    axarr[1+ispec].plot([(x-4040000.0)/1000 for x in [0]+DATA_this[0,0:itc].tolist()],[conc_arr[0]*conc2parts]+(conc_arr[0:itc]*conc2parts).tolist(),'k--',color=nrncols[iinput],lw=1.0,dashes=(1,2))         #before scale change point
    axarr[1+ispec].plot([tc_scpoint + ((x-4040000.0)/1000 - tc_scpoint)*tc_sccoeff for x in DATA_this[0,itc:].tolist()],(conc_arr[itc:]*conc2parts).tolist(),'k--',color=nrncols[iinput],lw=1.0,dashes=(1,2)) #after scale change point
    ireached95 = next((i for i,x in enumerate(conc_arr.tolist()) if x >= conc_arr[itend]*0.95))
    print 'Ca_flux = '+str(Ca_input_fluxes[timecourse_iinputs[iinput]])+': 95% of '+timecourse_species_titles[ispec]+' reached at time '+str((DATA_this[0,ireached95]-4040000.0)/1000)

    #Plot the curve discontinuity markers
    yl = axarr[1+ispec].get_ylim(); yeps = (yl[1]-yl[0])/30.0; xeps = 1.2
    axarr[1+ispec].plot([tc_scpoint-2*xeps+xeps*1,tc_scpoint+xeps*1],[conc_arr[itc]*conc2parts-yeps,conc_arr[itc]*conc2parts+yeps],'k-',lw=2)
    axarr[1+ispec].plot([tc_scpoint-xeps+xeps*1-xeps*1.6,tc_scpoint-xeps+xeps*1+xeps*1.6],[conc_arr[itc]*conc2parts-yeps*1.6,conc_arr[itc]*conc2parts+yeps*1.6],'w-',lw=1)
    print "species ="+str(timecourse_species[ispec][0:min(len(timecourse_species[ispec]),3)])+", max "+str(max([conc_arr[0]]+conc_arr.tolist()))+", yeps = "+str(yeps)+", iax="+str(1+ispec)
  axarr[0].plot([-100,0,0,900],[0,0,Ca_input_fluxes[timecourse_iinputs[iinput]],Ca_input_fluxes[timecourse_iinputs[iinput]]],'k-',color=cols[iinput],lw=1.0)
  yl = axarr[0].get_ylim()
  yeps = (yl[1]-yl[0])/30.0
  axarr[0].plot([tc_scpoint-2*xeps+xeps*1,tc_scpoint+xeps*1],[Ca_input_fluxes[timecourse_iinputs[iinput]]-yeps,Ca_input_fluxes[timecourse_iinputs[iinput]]+yeps],'k-',lw=2)
  axarr[0].plot([tc_scpoint-xeps+xeps*1-xeps*1.6,tc_scpoint-xeps+xeps*1+xeps*1.6],[Ca_input_fluxes[timecourse_iinputs[iinput]]-yeps*1.6,Ca_input_fluxes[timecourse_iinputs[iinput]]+yeps*1.6],'w-',lw=1)

for i in range(0,6):
  axarr[i].set_xlim([-25,100])
  axarr[i].set_xticks([0,30,60,100])
axarr[5].set_xticklabels(['0','30','60','250'])
for i in range(0,len(axarr)):
  axarr[i].set_position([0.09,0.88-0.096*i,0.11,0.08])
  for tick in axarr[i].xaxis.get_major_ticks() + axarr[i].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  axarr[i].spines['top'].set_visible(False)
  axarr[i].spines['right'].set_visible(False)


species = [[['PLC', 'PLCCa', 'PLCCaGqaGTP', 'PLCGqaGTP', 'PLCCaPip2', 'PLCCaGqaGTPPip2', 'PLCCaDAG', 'PLCCaGqaGTPDAG'], ['PLA2', 'CaPLA2', 'CaPLA2Pip2'], ['DGL', 'CaDGL', 'DAGCaDGL'], ['AC1GsaGTPCaMCa4', 'AC1GsaGTPCaMCa4ATP', 'AC1GiaGTPCaMCa4', 'AC1GiaGTPCaMCa4ATP', 'AC1GsaGTPGiaGTPCaMCa4', 'AC1GsGiCaMCa4ATP', 'AC1CaMCa4', 'AC1CaMCa4ATP', 'AC8CaMCa4', 'AC8CaMCa4ATP', 'PDE1CaMCa4', 'PDE1CaMCa4cAMP', 'NgCaM', 'CaM', 'CaMCa2', 'CaMCa3', 'CaMCa4', 'PP2BCaM', 'PP2BCaMCa2', 'PP2BCaMCa3', 'PP2BCaMCa4', 'CKCaMCa4', 'CKpCaMCa4', 'Complex', 'pComplex', 'CKpCaMCa4PP1', 'Ip35PP2BCaMCa4', 'Ip35PP1PP2BCaMCa4', 'PP1PP2BCaMCa4', 'GluR1_CKCam', 'GluR1_CKpCam', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam']], #ca bindings
           [['PLC', 'PLCCa', 'PLCCaGqaGTP', 'PLCGqaGTP', 'PLCCaPip2', 'PLCCaGqaGTPPip2', 'PLCCaDAG', 'PLCCaGqaGTPDAG'], ['PLA2', 'CaPLA2', 'CaPLA2Pip2'], ['DGL', 'CaDGL', 'DAGCaDGL'], ['AC1GsaGTPCaMCa4', 'AC1GsaGTPCaMCa4ATP', 'AC1GiaGTPCaMCa4', 'AC1GiaGTPCaMCa4ATP', 'AC1GsaGTPGiaGTPCaMCa4', 'AC1GsGiCaMCa4ATP', 'AC1CaMCa4', 'AC1CaMCa4ATP', 'AC8CaMCa4', 'AC8CaMCa4ATP', 'PDE1CaMCa4', 'PDE1CaMCa4cAMP', 'NgCaM', 'CaM', 'CaMCa2', 'CaMCa3', 'CaMCa4', 'PP2BCaM', 'PP2BCaMCa2', 'PP2BCaMCa3', 'PP2BCaMCa4', 'CKCaMCa4', 'CKpCaMCa4', 'Complex', 'pComplex', 'CKpCaMCa4PP1', 'Ip35PP2BCaMCa4', 'Ip35PP1PP2BCaMCa4', 'PP1PP2BCaMCa4', 'GluR1_CKCam', 'GluR1_CKpCam', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam']], #ca bindings, copy, inset
           [['PKAcLR', 'PKAcpLR', 'PKAcppLR', 'PKAcpppLR', 'PKAcR', 'PKAcpR', 'PKAcppR', 'PKAcpppR', 'PKA', 'PKAcAMP4', 'PKAc', 'I1PKAc', 'GluR1_PKAc', 'GluR1_S831_PKAc', 'GluR1_memb_PKAc', 'GluR1_memb_S831_PKAc', 'PKAcPDE4', 'PKAc_PDE4_cAMP']],
           [['CK', 'CKCaMCa4', 'CKpCaMCa4', 'CKp', 'Complex', 'pComplex', 'CKpPP1', 'CKpCaMCa4PP1', 'GluR1_CKCam', 'GluR1_CKpCam', 'GluR1_CKp', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_S845_CKp', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_CKp', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp']],
           [['GluR1_PKCt', 'GluR1_PKCp', 'GluR1_S845_PKCt', 'GluR1_S845_PKCp', 'GluR1_memb_PKCt', 'GluR1_memb_PKCp', 'GluR1_memb_S845_PKCt', 'GluR1_memb_S845_PKCp', 'PKC', 'PKCCa', 'PKCt', 'PKCp', 'GluR2_PKCt', 'GluR2_PKCp', 'GluR2_memb_PKCt', 'GluR2_memb_PKCp']],
           [['Ca']]
           ]

coeffs = [[[0,1,1,0,1,1,1,1],[0,1,1],[0,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,0,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1]],       #PLC,PLA2,CaM,DGL
          [[0,1,1,0,1,1,1,1],[0,1,1],[0,1,1],[1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,0,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1]],       #PLC,PLA2,CaM,DGL (copy)
          [[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0,0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]],                 #PKA
          [[0,0,1,1,2,2,1,1,0,1,1,0,1,1,0,1,1,0,1,1]],                                             #CK
          [[1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1]], #PKC
          [[1 for i in range(0,len(species[5][0]))]],
          ]
ifluxes_to_draw = [[2],[2],[0,1,2],[2],[0,1,2],[0,1,2]]
fluxes = [0.0, 0.005, 0.05]

#titles = ['CaM', 'immob. buffer', 'PP2B', 'PP1', 'AC1', 'AC8', 'PDE1', 'PDE4', 'PKAc', 'PLC', 'DAG', 'PKC', 'PP2A', 'CK']
titles = [['PLC','PLA2','DGL','CaM'], ['','','',''], ['PKAc'], ['CaMKII'], ['PKC'], ['Ca']]

itend = argmin(abs(DATA_all[0]['DATA'][0]-(4340000)))
itpost = argmin(abs(DATA_all[0]['DATA'][0]-(4940000)))
firstAbove1mMCa = -1

#Zoomed in area:
polygon = Polygon(array([[0, 150, 150, 0],[0,0,0.05,0.05]]).T, True)
p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
p.set_facecolor('#FFD9D9')
p.set_edgecolor('none')
axarr[6].add_collection(p)

ligands = ['','','NE','','ACh+Glu','']

for ispec in range(0,len(species)):
 for iiflux in range(0,len(ifluxes_to_draw[ispec])):
  iflux = ifluxes_to_draw[ispec][iiflux]
  if iflux == 0:
   DATA_all = DATA_all_0
   Ca_input_fluxes = Ca_input_fluxes_0
  if iflux == 1:
   DATA_all = DATA_all_1
   Ca_input_fluxes = Ca_input_fluxes_1
  if iflux == 2:
   DATA_all = DATA_all_2
   Ca_input_fluxes = Ca_input_fluxes_2
  for ispecgroup in range(0,len(species[ispec])):
   concs_end = []
   concs_post = []
   concs_sum = []
   concs_arr =   []
   concs_max =   []
   concs_end_ref = []
   if len(ifluxes_to_draw[ispec]) == 1:
     mylabel = titles[ispec][ispecgroup]
     mycol = cols[ispecgroup]
   else:
     mylabel = str(fluxes[iflux])+' '+ligands[ispec]
     mycol = cols[iflux]
   for iinput in range(0,len(DATA_all)):
    DATA = DATA_all[iinput]['DATA']
    headers = DATA_all[iinput]['headers']
    conc_end = 0
    conc_post = 0
    conc_sum = 0
    conc_arr = zeros(DATA.shape[1])
    conc_end_ref = 0
    for iispecie in range(0,len(species[ispec][ispecgroup])):
      ispecieind = -1
      for iispecieDATA in range(0,len(headers)):
        mystr = headers[iispecieDATA]
        firstspace = mystr.find(' ')
        if firstspace >= 0:
          mystr = mystr[0:firstspace]
        if mystr == species[ispec][ispecgroup][iispecie]:
          ispecieind = iispecieDATA
          break
      if ispecieind == -1:
        print species[ispec][ispecgroup][iispecie]+' not found in headers'
        continue
      conc_end = conc_end + DATA[ispecieind,itend]*coeffs[ispec][ispecgroup][iispecie]
      conc_post = conc_post + DATA[ispecieind,itpost]*coeffs[ispec][ispecgroup][iispecie]
      conc_arr = conc_arr + DATA[ispecieind,:]*coeffs[ispec][ispecgroup][iispecie]
      conc_sum = conc_sum + sum(DATA[ispecieind,:])*coeffs[ispec][ispecgroup][iispecie]
      conc_end_ref = conc_end_ref + DATA[ispecieind,itend]*max(1,coeffs[ispec][ispecgroup][iispecie])
    concs_end.append(conc_end)
    concs_post.append(conc_post)
    concs_sum.append(conc_sum)
    concs_arr.append(conc_arr)
    concs_max.append(max(conc_arr))
    concs_end_ref.append(conc_end_ref)
   axarr[6+ispec].plot(Ca_input_fluxes, array(concs_end)/concs_end_ref[0],'k-',color=mycol,lw=1.0,label=mylabel,zorder=5)    #Relative to max. theoretical concentration
   if ispec == 0:
     axarr[6+ispec].legend(fontsize=6,loc=2,frameon=False)
   elif ispec == 2:
     axarr[6+ispec].legend(fontsize=6,loc="center", bbox_to_anchor=(0.84,0.22),frameon=False)
   elif ispec == 4:
     axarr[6+ispec].legend(fontsize=6,loc="center", bbox_to_anchor=(0.75,0.6),frameon=False)
   print titles[ispec][ispecgroup]+": min,max concs_end_ref = "+str(min(concs_end_ref))+", "+str(max(concs_end_ref))+", max chosen = "+str(max(concs_end))
   if species[ispec]==[['Ca']]:
     firstAbove1mMCa = [i for i,x in enumerate(concs_end_ref) if x >= 1.0][0]

axarr[6].set_position([0.27, 0.4, 0.4, 0.56])
axarr[7].set_position([0.52, 0.435, 0.14, 0.125])
axarr[8].set_position([0.74, 0.78, 0.18, 0.17])
axarr[9].set_position([0.74, 0.59, 0.18, 0.17])
axarr[10].set_position([0.74, 0.4, 0.18, 0.17])

axarr[6].set_ylabel('fraction of Ca$^{2+}$-bound in steady state',fontsize=6)
axarr[8].set_yticks([0,0.001,0.002])
axarr[9].set_yticks([0,0.5,1.0])
axarr[10].set_yticks([0,0.5,1.0])
axarr[8].set_ylabel('PKA',fontsize=6)
axarr[9].set_ylabel('CaMKII',fontsize=6)
axarr[10].set_ylabel('PKC',fontsize=6)

for iax in range(0,len(axarr)-6):
  polygon = Polygon(array([[Ca_input_fluxes[firstAbove1mMCa], 3000, 3000, Ca_input_fluxes[firstAbove1mMCa]],[0,0,1,1]]).T, True)
  p = PatchCollection([polygon], cmap=matplotlib.cm.jet)
  p.set_facecolor('#E0E0E0')
  p.set_edgecolor('none')
  axarr[6+iax].add_collection(p)

  axarr[6+iax].set_xticks([0,500,1000])
  axarr[6+iax].set_xlim([0,1200])
  axarr[6+iax].set_ylim([0,axarr[6+iax].get_ylim()[1]])
  for tick in axarr[6+iax].xaxis.get_major_ticks() + axarr[6+iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  axarr[6+iax].spines['top'].set_visible(False)
  axarr[6+iax].spines['right'].set_visible(False)

axarr[6].set_xticks([0,250,500,750,1000])
axarr[7].set_xticks([0,50,100,150])
axarr[7].set_xlim([0,150])
axarr[7].set_ylim([0,0.05])

for iax in [0,1,2,3,4]+[8,9]:
  axarr[iax].set_xticklabels([])
axarr[5].set_xlabel('time (s)',fontsize=6)
for iax in [6,10]:
  axarr[iax].set_xlabel('Ca$^{2+}$ flux (particles/ms)',fontsize=6)
axarr[0].set_ylabel('Ca$^{2+}$ ions/ms',fontsize=6)
for iax in [1,2,3,4,5]:
  axarr[iax].set_ylabel(timecourse_species_titles[iax-1],fontsize=6)

for iax in range(11,len(axarr)):
  axarr[iax].set_visible(False)

for iax in range(0,6):
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.09, pos.y1 - 0.02, chr(ord('A')+iax), fontsize=10)
myiax = 6
for iax in range(6,11):
  if iax == 7:
    continue
  pos = axarr[iax].get_position()
  f.text(pos.x0 - 0.05, pos.y1 - 0.02, chr(ord('A')+myiax), fontsize=10)
  myiax = myiax + 1
  
f.savefig("fig3.eps")

