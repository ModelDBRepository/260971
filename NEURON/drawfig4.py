#drawfig4.py: Draws the outcomes of 4xHFS and LFS protocols when different fluxes of neurotransmitters are absent.
#Tuomo Maki-Marttunen, 2019-2020
#CC BY 4.0
import matplotlib
matplotlib.use('Agg')
from pylab import *
import scipy.io
import sys
import calcconds

filename_control =     '../NeuroRD/tstop500000_tol0.01_onset40000.0_n100_freq100.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains4_trainT4000.0_8seeds.mat'
filename_control_LFS = '../NeuroRD/tstop500000_tol0.01_onset40000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains1_trainT100000.0_8seeds.mat'
filenames = ['../NeuroRD/tstop500000_tol0.01_onset40000.0_n100_freq100.0_dur3.0_flux0.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains4_trainT4000.0_8seeds.mat',
             '../NeuroRD/tstop500000_tol0.01_onset40000.0_n100_freq100.0_dur3.0_flux1900.0_Lflux0.0_Gluflux20.0_AChflux20.0_Ntrains4_trainT4000.0_8seeds.mat',
             '../NeuroRD/tstop500000_tol0.01_onset40000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux0.0_AChflux20.0_Ntrains1_trainT100000.0_8seeds.mat',
             '../NeuroRD/tstop500000_tol0.01_onset40000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux0.0_AChflux0.0_Ntrains1_trainT100000.0_8seeds.mat',
             '../NeuroRD/tstop500000_tol0.01_onset40000.0_n100_freq100.0_dur3.0_flux1900.0_Lflux10.0_Gluflux0.0_AChflux0.0_Ntrains4_trainT4000.0_8seeds.mat']
filename_nrn_control =     'nrn_tstop5000000_tol1e-06_onset4040000.0_n100_freq100.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains4_trainT4000.0.mat'
filename_nrn_control_LFS = 'nrn_tstop5000000_tol1e-06_onset4040000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains1_trainT100000.0.mat'
filenames_nrn = ['nrn_tstop5000000_tol1e-06_onset4040000.0_n100_freq100.0_dur3.0_flux0.0_Lflux10.0_Gluflux20.0_AChflux20.0_Ntrains4_trainT4000.0.mat',
                 'nrn_tstop5000000_tol1e-06_onset4040000.0_n100_freq100.0_dur3.0_flux1900.0_Lflux0.0_Gluflux20.0_AChflux20.0_Ntrains4_trainT4000.0.mat',
                 'nrn_tstop5000000_tol1e-06_onset4040000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux0.0_AChflux20.0_Ntrains1_trainT100000.0.mat',
                 'nrn_tstop5000000_tol1e-06_onset4040000.0_n900_freq5.0_dur3.0_flux1900.0_Lflux10.0_Gluflux0.0_AChflux0.0_Ntrains1_trainT100000.0.mat',
                 'nrn_tstop5000000_tol1e-06_onset4040000.0_n100_freq100.0_dur3.0_flux1900.0_Lflux10.0_Gluflux0.0_AChflux0.0_Ntrains4_trainT4000.0.mat']

species = [ [ ['GluR1_memb', 'GluR1_memb_S845', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_PKAc', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_CKp', 'GluR1_memb_PKCt', 'GluR1_memb_PKCp', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCt', 'GluR1_memb_S845_PKCp', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'] ],
            [ ['GluR2_memb', 'GluR2_memb_PKCt', 'GluR2_memb_PKCp', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'] ] ]
labels = [['membr. GluR1','cyt. GluR1'],['membr. GluR2','cyt. GluR2']]
titles = ['No Ca$^{2+}$','No $\\beta$-adrenergic ligand','No mGluR activation','No M1 or mGluR activation','No M1 or mGluR activation']
irows = [0, 0, 1, 1, 0]
icols = [1, 3, 1, 2, 2]

f,axs = subplots(2,3)
for ifile in range(0,2):
  for iax in range(0,3):
    axs[ifile,iax].set_position([0.415+0.245*iax,0.74-0.28*ifile-(0.04*ifile==4),0.165,0.20])
  axs[ifile,2].set_position([0.1,0.74-0.28*ifile-(0.04*ifile==4),0.245,0.20])


mesh_input_file = open('mesh_general.out','r')
mesh_firstline = mesh_input_file.readline()
mesh_secondline = mesh_input_file.readline()
mesh_values = mesh_secondline.split()
my_volume = float(mesh_values[-2])*1e-15 #litres
mesh_input_file.close()

locator_params(nticks=4)

#Control HFS
MAT = scipy.io.loadmat(filename_control)
DATA_all_control = MAT['DATA']
header_strs_control = MAT['headers']
for i in range(0,len(header_strs_control)):
  first_space = header_strs_control[i].find(' ')
  header_strs_control[i] = header_strs_control[i][:first_space]
inddict_control = {}
for iheader in range(4,len(header_strs_control)):
  inddict_control[header_strs_control[iheader]] = iheader-4

#Control LFS
MAT_LFS = scipy.io.loadmat(filename_control_LFS)
DATA_all_control_LFS = MAT_LFS['DATA']
header_strs_control_LFS = MAT_LFS['headers']
for i in range(0,len(header_strs_control_LFS)):
  first_space = header_strs_control_LFS[i].find(' ')
  header_strs_control[i] = header_strs_control_LFS[i][:first_space]
inddict_control_LFS = {}
for iheader in range(4,len(header_strs_control_LFS)):
  inddict_control_LFS[header_strs_control_LFS[iheader]] = iheader-4

DATANRN_all_all = scipy.io.loadmat(filename_nrn_control)
DATANRN_all_control = {}
for ikey in range(0,len(DATANRN_all_all['headers'])):
  mykey = DATANRN_all_all['headers'][ikey][0:DATANRN_all_all['headers'][ikey].find(' ')]
  DATANRN_all_control[mykey] = DATANRN_all_all['DATA'][ikey]
times_nrn_control = DATANRN_all_control['tvec']

DATANRN_all_all_LFS = scipy.io.loadmat(filename_nrn_control_LFS)
DATANRN_all_control_LFS = {}
for ikey in range(0,len(DATANRN_all_all_LFS['headers'])):
  mykey = DATANRN_all_all_LFS['headers'][ikey][0:DATANRN_all_all_LFS['headers'][ikey].find(' ')]
  DATANRN_all_control_LFS[mykey] = DATANRN_all_all_LFS['DATA'][ikey]
times_nrn_control_LFS = DATANRN_all_control_LFS['tvec']

conds_control = calcconds.calcconds(filename_control, filename_nrn_control)
conds_control_LFS = calcconds.calcconds(filename_control_LFS, filename_nrn_control_LFS)

cols_all = ['#440154', '#470f62', '#481d6f', '#472a79', '#453681', '#414387', '#3c4f8a', '#37598c', '#32648e', '#2d6f8e', '#29788e', '#26828e', '#228b8d', '#1f958b', '#1e9f88', '#22a884', '#2bb17e', '#3bbb75', '#4dc36b', '#62cb5f', '#7ad251'][::10] + ['#dddd00']
cols_nrn_all = ['#360043', '#390c4f', '#391759', '#392161', '#372b67', '#34366c', '#303f6f', '#2c4770', '#285071', '#245872', '#216072', '#1e6872', '#1b6f71', '#19776f', '#187f6d', '#1b866a', '#238e65', '#2f955e', '#3e9c56', '#4fa24c', '#61a841'][::10] + ['#cccc00']
iicols = [[0,1,3,2,4],[0,1,2,3,4]]
cols_all = [cols_all[i] for i in [0,1,3,2]]
cols_nrn_all = [cols_nrn_all[i] for i in [0,1,3,2]]


GluR1R2ratios = []
GluR1R2ratios_controls = []
conds_all = []
conds_control_all = []
for ifile in range(0,len(filenames)):
  filename = filenames[ifile]
  print "Loading "+filename
  MAT = scipy.io.loadmat(filename)

  axarr = [axs[irows[ifile],0],axs[irows[ifile],1],axs[irows[ifile],2]]
  DATA_all = MAT['DATA']
  header_strs = MAT['headers']
  for i in range(0,len(header_strs)):
    first_space = header_strs[i].find(' ')
    header_strs[i] = header_strs[i][:first_space]
  inddict = {}
  for iheader in range(4,len(header_strs)):
    inddict[header_strs[iheader]] = iheader-4

  DATANRN_all_all = scipy.io.loadmat(filenames_nrn[ifile])
  DATANRN_all = {}
  for ikey in range(0,len(DATANRN_all_all['headers'])):
    mykey = DATANRN_all_all['headers'][ikey][0:DATANRN_all_all['headers'][ikey].find(' ')]
    DATANRN_all[mykey] = DATANRN_all_all['DATA'][ikey]
  times_nrn = DATANRN_all['tvec']

  conds = calcconds.calcconds(filename, filenames_nrn[ifile])
  conds_all.append(conds[:])

  Nskip = 1
  if len(sys.argv) > 3:
    Nskip = int(sys.argv[3])

  basalmaxes = []

  DATA_this_control  = DATA_all_control 
  DATANRN_this_control  = DATANRN_all_control 
  times_nrn_control_this = times_nrn_control
  conds_control_this = conds_control
  if ifile > 1 and ifile != 4:
    DATA_this_control  = DATA_all_control_LFS 
    DATANRN_this_control  = DATANRN_all_control_LFS 
    times_nrn_control_this = times_nrn_control_LFS
    conds_control_this = conds_control_LFS

  conds_control_all.append(conds_control_this[:])

  T = float(filename[filename.find('tstop')+5:filename.find('_',filename.find('tstop'))])
  Duration = float(filename[filename.find('tstop')+5:filename.find('_',filename.find('tstop'))])
  times = [T/(DATA_all.shape[0]-1)*i for i in range(0,DATA_all.shape[0])]

  cols = [cols_all[icols[ifile]]]
  cols_nrn = [cols_nrn_all[icols[ifile]]]
  col_control = cols_all[0]
  col_control_nrn = cols_nrn_all[0]

  for iax in range(0,2):
    times = conds[4]
    times_nrn = conds[5]
    mytimecourse_control = conds_control_this[2][3+2*iax]     #[3] for GluR1, [5] for GluR2
    mytimecourse_nrn_control = conds_control_this[3][3+2*iax] #[3] for GluR1, [5] for GluR2
    mytimecourse = conds[2][3+2*iax]                          #[3] for GluR1, [5] for GluR2
    mytimecourse_nrn = conds[3][3+2*iax]                      #[3] for GluR1, [5] for GluR2

    axarr[iax].plot(array(times[::Nskip])/1000,mytimecourse_control[::Nskip]/6.022e23/my_volume*1e9,'k-',color=col_control,linewidth=0.7) #,label='control')
    axarr[iax].plot([x/1000-4000 for x in times_nrn_control_this[::Nskip]],mytimecourse_nrn_control[::Nskip]/6.022e23/my_volume*1e9,'k--',color=col_control_nrn,linewidth=1.5,dashes=[1,2],zorder=100)

    axarr[iax].plot(array(times[::Nskip])/1000,mytimecourse[::Nskip]/6.022e23/my_volume*1e9,color=cols[0],linewidth=0.7) #,label=titles[ifile])
    axarr[iax].plot([x/1000-4000 for x in times_nrn[::Nskip]],mytimecourse_nrn[::Nskip]/6.022e23/my_volume*1e9,'k--',color=cols_nrn[0],linewidth=1.5,dashes=[1+0.1*ifile,2+0.2*ifile],zorder=100)

  axarr[2].plot([x/1000 for x in times],conds[0],color=cols[0],linewidth=0.7,label=titles[ifile])
  axarr[2].plot([x/1000-4000 for x in times_nrn[::Nskip]],conds[1],'k--',color=cols_nrn[0],linewidth=1.5,dashes=[1+0.1*ifile,2+0.2*ifile],zorder=100)
  if ifile == 3 or ifile == 4:
    axarr[2].plot([x/1000 for x in times],conds_control_this[0],color=col_control,linewidth=0.7,label='control')
    axarr[2].plot([x/1000-4000 for x in times_nrn_control_this[::Nskip]],conds_control_this[1],'k--',color=col_control_nrn,linewidth=1.5,dashes=[1,2],zorder=100)
  GluR1R2ratios.append(1.0*conds[3][3]/(conds[3][3]+conds[3][5]))
  GluR1R2ratios_controls.append(1.0* conds_control_this[3][3]/( conds_control_this[3][3]+ conds_control_this[3][5]))
  print("conds baseline = "+str(conds[1][0]))
  
  for iax in range(0,3):
    axarr[iax].set_xlim([0,500])
    basalmaxes.append(max(mytimecourse[0:499])/6.022e23/my_volume*1e9)
    for tick in axarr[iax].yaxis.get_major_ticks()+axarr[iax].xaxis.get_major_ticks():
      tick.label.set_fontsize(8)
      tick.label2.set_fontsize(8)
    xticks = axarr[iax].get_xticks()
    if len(xticks) > 6:
      axarr[iax].set_xticks(xticks[::2])
    yticks = axarr[iax].get_yticks()
    if len(yticks) > 4:
      axarr[iax].set_yticks(yticks[::2])
    if iax == 2:
      if ifile < 2 or ifile == 4:
        axarr[iax].set_ylabel('$\\mathbf{(HFS)}$\npS',fontsize=8)
      else:  
        axarr[iax].set_ylabel('$\\mathbf{(LFS)}$\npS',fontsize=8)
    else:
      axarr[iax].set_ylabel('nM',fontsize=8)
    for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
      tick.label.set_fontsize(6)
    axarr[iax].spines['top'].set_visible(False)
    axarr[iax].spines['right'].set_visible(False)
    axarr[iax].tick_params(top='off', right='off')
    #axarr[iax].set_position([0.08+0.48*iax,0.82-0.18*ifile,0.4,0.16])
    #0.4+0.1+0.28*iax
    if iax == 2:
      myloc = [2,10,1,10,10][ifile]
      mybb = [nan, (0.37,0.5), nan, (0.63,0.94),(0.4,0.84)][ifile]
      try: #this gave errors in one python setup
        if myloc == 10:
          leg = axarr[iax].legend(fontsize=6,loc=myloc,bbox_to_anchor=mybb,facecolor='None',handletextpad=0.05,handlelength=1.0,labelspacing=0.2) #,loc=1+3*(iax != 1))
        else:
          leg = axarr[iax].legend(fontsize=6,loc=myloc,facecolor='None',handletextpad=0.05,handlelength=1.0,labelspacing=0.2) #,loc=1+3*(iax != 1))
        leg.get_frame().set_linewidth(0)
      except:
        print "Legend frames not made invisible"
    #axarr[iax].set_title(titles[ifile],fontsize=8)
  axarr[0].set_ylim([0,195])
  axarr[0].set_yticks([0,80,160])
  axarr[1].set_ylim([0,210])
  axarr[1].set_yticks([0,100,200])
  if ifile < 2 or ifile == 4:
    axarr[2].set_ylim([0,230])
  else:
    axarr[2].set_ylim([20,42])

axs[0,0].set_title('GluR1 at membrane',fontsize=8)
axs[0,1].set_title('GluR2 at membrane',fontsize=8)
axs[0,2].set_title('Total synaptic conductance',fontsize=8)

for iax in range(0,3):
  axs[0,iax].set_xticklabels([])
  axs[1,iax].set_xlabel('time (s)',fontsize=8)
for ifile in range(0,2):
  for iax in range(0,2):
    f.text(0.36+0.245*iax, 0.82-0.28*ifile+0.13-0.04*(ifile==4),chr(65+ifile*4+1+iax),fontsize=11)
  f.text(0.04, 0.82-0.21*ifile+0.13-0.04*(ifile==4),chr(65+ifile*4),fontsize=11)


axs[0,2].set_yticks([0,100,200])
axs[1,2].set_yticks([20,30,40])

conds_check_control = calcconds.calcconds_nrn_withTCsN(filename_nrn_control)

conds_hom1 = [12.4, 18.9]
conds_hom2 = 2.2
conds_het = 2.5

axnew = []
ifiles2 = [4,3]
for iifile in range(0,len(ifiles2)):
  ifile = ifiles2[iifile]
  axnew.append(f.add_axes([0.905,0.78-0.21*iifile,0.083,0.18]))
  axnew[iifile].bar(1,100.*GluR1R2ratios[ifile][-1],facecolor=cols_all[icols[ifile]])
  axnew[iifile].bar(2,100.*GluR1R2ratios_controls[ifile][-1],facecolor=col_control)

  conds_control_this = conds_control_LFS if ifile > 1 and ifile != 4 else conds_control
  ENhom1_np_nrn_control = (conds_control_this[2][3] + conds_control_this[2][5])/4.0 * (conds_control_this[2][3]-conds_control_this[2][1])**4/(conds_control_this[2][3] + conds_control_this[2][5])**4                       
  ENhom1_p_nrn_control = (conds_control_this[2][3] + conds_control_this[2][5])/4.0 * (conds_control_this[2][3]**4 - (conds_control_this[2][3]-conds_control_this[2][1])**4)/(conds_control_this[2][3] + conds_control_this[2][5])**4 
  ENhom2_nrn_control = (conds_control_this[2][3] + conds_control_this[2][5])/4.0 * (conds_control_this[2][5]/(conds_control_this[2][3] + conds_control_this[2][5]))**4
  ENhet_nrn_control = (conds_control_this[2][3] + conds_control_this[2][5])/4.0 * (1 - (conds_control_this[2][3]/(conds_control_this[2][3] + conds_control_this[2][5]))**4 - (conds_control_this[2][5]/(conds_control_this[2][3] + conds_control_this[2][5]))**4)
  Egtot_nrn_control = ENhom1_np_nrn_control*conds_hom1[0] + ENhom1_p_nrn_control*conds_hom1[1] + ENhom2_nrn_control*conds_hom2 + ENhet_nrn_control*conds_het

  endvals = [conds_all[ifile][3][3][-1],conds_all[ifile][3][5][-1]]
  endvals_control = [conds_control_all[ifile][3][3][-1],conds_control_all[ifile][3][5][-1]]
  phom = (endvals[0]/(endvals[0] + endvals[1]))**4
  phom_control = (endvals_control[0]/(endvals_control[0] + endvals_control[1]))**4

  ENhom1_np_nrn = (conds_all[ifile][3][3] + conds_all[ifile][3][5])/4.0 * (conds_all[ifile][3][3]-conds_all[ifile][3][1])**4/(conds_all[ifile][3][3] + conds_all[ifile][3][5])**4
  ENhom1_p_nrn = (conds_all[ifile][3][3] + conds_all[ifile][3][5])/4.0 * (conds_all[ifile][3][3]**4 - (conds_all[ifile][3][3]-conds_all[ifile][3][1])**4)/(conds_all[ifile][3][3] + conds_all[ifile][3][5])**4
  ENhom2_nrn = (conds_all[ifile][3][3] + conds_all[ifile][3][5])/4.0 * (conds_all[ifile][3][5]/(conds_all[ifile][3][3] + conds_all[ifile][3][5]))**4
  ENhet_nrn = (conds_all[ifile][3][3] + conds_all[ifile][3][5])/4.0 * (1 - (conds_all[ifile][3][3]/(conds_all[ifile][3][3] + conds_all[ifile][3][5]))**4 - (conds_all[ifile][3][5]/(conds_all[ifile][3][3] + conds_all[ifile][3][5]))**4)
  homgrel = (ENhom1_np_nrn*conds_hom1[0] + ENhom1_p_nrn*conds_hom1[1])/(ENhom1_np_nrn*conds_hom1[0] + ENhom1_p_nrn*conds_hom1[1] + ENhom2_nrn*conds_hom2 + ENhet_nrn*conds_het)
  homgrel_control = (ENhom1_np_nrn_control*conds_hom1[0] + ENhom1_p_nrn_control*conds_hom1[1])/(ENhom1_np_nrn_control*conds_hom1[0] + ENhom1_p_nrn_control*conds_hom1[1] + ENhom2_nrn_control*conds_hom2 + ENhet_nrn_control*conds_het)
  phom = (conds_all[ifile][3][3]/(conds_all[ifile][3][3] + conds_all[ifile][3][5]))**4
  phom_control = (conds_control_all[ifile][3][3]/(conds_control_all[ifile][3][3] + conds_control_all[ifile][3][5]))**4

  axnew[iifile].bar(4,100.*phom[-1],facecolor=cols_all[icols[ifile]])
  axnew[iifile].bar(5,100.*phom_control[-1],facecolor=col_control)

  axnew[iifile].bar(7,100.*homgrel[-1],facecolor=cols_all[icols[ifile]])
  axnew[iifile].bar(8,100.*homgrel_control[-1],facecolor=col_control)
  print "phom_control = "+str(phom_control)
  print "homgrel_control[-1] = "+str(homgrel_control[-1])

  for tick in axnew[iifile].xaxis.get_major_ticks() + axnew[iifile].yaxis.get_major_ticks():
    tick.label.set_fontsize(6)
  axnew[iifile].spines['top'].set_visible(False)
  axnew[iifile].spines['right'].set_visible(False)
  axnew[iifile].tick_params(top='off', right='off')
  axnew[iifile].set_xticks([1.5,4.5,7.5])
  if iifile == 1:
    axnew[iifile].set_xticklabels(['R1 fraction','$p_{\mathrm{hom-R1}}$','hom-R1 contr.'],rotation=90,ha='center',fontsize=8)
  else:
    axnew[iifile].set_xticklabels([])
  axnew[iifile].set_yticks([0,50,100])
  axnew[iifile].set_ylabel('%',rotation=0,fontsize=8)
  axnew[iifile].tick_params(axis='x',width=0.2,length=0.0)
  #f.text(0.82, 0.2,'P',fontsize=11)
  f.text(0.84, 0.82-0.24*iifile+0.13,chr(68+iifile*4),fontsize=11)

f.savefig('fig4.eps')



