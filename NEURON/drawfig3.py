#drawfig3.py: Draws the pathway activations and synaptic conductance+GluR membrane expressions for 4xHFS (fig3A.eps) and the synaptic conductance+GluR membrane expressions for LFS
#Tuomo Maki-Marttunen, 2019-2020
#CC BY 4.0
import matplotlib
matplotlib.use('Agg')
from pylab import *
import scipy.io
import sys
import itertools
from os.path import exists

##################### Fig 4, 4xHFS panels #######################

close("all")

species = [ ['Ca'], ['CaMCa2'], [['AC1GsaGTPCaMCa4','AC1GsaGTPCaMCa4ATP','AC1GiaGTPCaMCa4','AC1GiaGTPCaMCa4ATP','AC1GsaGTPGiaGTPCaMCa4','AC1GsGiCaMCa4ATP','AC1CaMCa4','AC1CaMCa4ATP','AC8CaMCa4','AC8CaMCa4ATP','PDE1CaMCa4','PDE1CaMCa4cAMP','CaMCa4','PP2BCaMCa4','CKCaMCa4','CKpCaMCa4','CKpCaMCa4PP1','Ip35PP2BCaMCa4','Ip35PP1PP2BCaMCa4','PP1PP2BCaMCa4']], [['CKp','CKpCaMCa4']], [['GluR1_S831', 'GluR1_S845_S831', 'GluR1_S831_PKAc', 'GluR1_S845_S831_PP1', 'GluR1_S831_PP1', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'],
                      ['GluR1_S845_S831', 'GluR1_S845_S831_PP1', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S845_S831', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S845_S831_PP2B'] ],
            [['L','LR','LRGs','PKAcLR','PKAcpLR','PKAcppLR','PKAcpppLR','pLR','ppLR','pppLR','ppppLR','ppppLRGi','ppppLRGibg','LRGsbg']], [['GsaGTP','AC1GsaGTP','AC1GsaGTPCaMCa4'],['GiaGTP','AC1GiaGTP','AC1GiaGTPCaMCa4']], 
             ['cAMP'], ['PKAc'], [['GluR1_S845', 'GluR1_S845_S831', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_S845_CKp', 'GluR1_S845_PKCt', 'GluR1_S845_PKCp', 'GluR1_S845_PP1', 'GluR1_S845_S831_PP1', 'GluR1_S845_PP2B', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S845', 'GluR1_memb_S845_S831', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCt', 'GluR1_memb_S845_PKCp', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'],
                                  ['GluR1_S845_S831', 'GluR1_S845_S831_PP1', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S845_S831', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S845_S831_PP2B'] ],
            [['ACh','AChM1R','AChM1RGq'], ['Glu','MGluR_Glu','MGluR_Gqabg_Glu']], [['GqaGTP','PLCGqaGTP','PLCCaGqaGTP','PLCCaGqaGTPPip2']], ['DAG'], ['PKCt','PKCp'], [['GluR2_S880', 'GluR2_S880_PP2A', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'] ],
            [['GluR1_memb', 'GluR1_memb_S845', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_PKAc', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_CKp', 'GluR1_memb_PKCt', 'GluR1_memb_PKCp', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCt', 'GluR1_memb_S845_PKCp', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'], 
             ['GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B']],
            [['GluR2_memb', 'GluR2_memb_PKCt', 'GluR2_memb_PKCp', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'] ] ]

ispecies_needed_for_cond = [[15,0], [15,1], [16,0]]
conds_hom1 = [12.4, 18.9]
conds_hom2 = 2.2
conds_het = 2.5
factorExceptions = [[15, [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]]]

inparticles = [ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2 ]

titles = ['','','','',['S831p', '2x phosph.'],
          '',['Gs','Gi'],'','',['S845p', '2x phosph.'],
          ['Glu','ACh'],'','',['transient','persistent'],'', ['total', 'S831p'], '', '']
ylabels = ['[Ca$^{2+}$]','[CaMCa2]','[CaMCa4]','[act. CaMKII]','[GluR1 phos.]',
          '[$\\beta$-adr. ligand]','[act. G prot.]','[cAMP]','[act. PKA]','[GluR1 phos]',
          '[Glu],[ACh]','[act. Gq]','[DAG]','[act. PKC]','[GluR2 S880p]', '[GluR1 at memb.]', '[GluR2 at memb.]', 'max. cond.']

yticks_all = [[], [], [0,10000], [0,3000], [],    [], [0,400,800], [], [], [0,40,80],    [], [], [], [], [],    [0,100,200], [0,100,200], [] ]

grouptitles = ['CaMKII pathway','Gs/Gi-cAMP-PKA pathway','Gq-PLC-PKC pathway']
protocol = 'HFS'
Caflux = 1900.0
Lflux = 10.0
Gluflux = 20.0
AChflux = 20.0
if len(sys.argv) > 1:
  Caflux = float(sys.argv[1])
if len(sys.argv) > 2:
  Lflux = float(sys.argv[2])
if len(sys.argv) > 3:
  Gluflux = float(sys.argv[3])
if len(sys.argv) > 4:
  AChflux = float(sys.argv[4])

ylims = [nan]*(len(species)+1)
ylims[3] = [0,4500]
ylims[4] = [0,290]
ylims[15] = [0,270]
ylims[16] = [0,270]

filename = '../NeuroRD/tstop500000_tol0.01_onset40000.0_n100_freq100.0_dur3.0_flux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(AChflux)+'_Ntrains4_trainT4000.0_8seeds.mat'
MAT = {}
if exists(filename):
  MAT = scipy.io.loadmat(filename)
  print "loaded "+filename

  DATA_all = MAT['DATA']
  header_strs = MAT['headers']
  for i in range(0,len(header_strs)):
    first_space = header_strs[i].find(' ')
    if first_space > -1:
      header_strs[i] = header_strs[i][:first_space]

  inddict = {}
  for iheader in range(4,len(header_strs)):
    inddict[header_strs[iheader]] = iheader-4
else:
  print filename+' not found. Variable MAT will remain empty and NeuroRD results will not be plotted'

Nskip = 1
if len(sys.argv) > 5:
  Nskip = int(sys.argv[5])

filename_nrn = 'nrn_noninterp_tstop5000000_tol1e-06_onset4040000.0_n100_freq100.0_dur3.0_flux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(AChflux)+'_Ntrains4_trainT4000.0.mat'
DATANRN_all = {}
if exists(filename_nrn):
  DATANRN_all_all = scipy.io.loadmat(filename_nrn)
  print "loaded "+filename_nrn
  for ikey in range(0,len(DATANRN_all_all['headers'])):
    mykey = DATANRN_all_all['headers'][ikey][0:DATANRN_all_all['headers'][ikey].find(' ')]
    DATANRN_all[mykey] = DATANRN_all_all['DATA'][ikey]
else:
  print filename_nrn+' not found. Variable DATANRN_all will remain empty and NEURON results will not be plotted'
nrncols_all = [ ['#0000AA','#00AA00','#AA0000','#770077'], ['#0000AA','#770077'] ]

mesh_input_file = open('mesh_general.out','r')
mesh_firstline = mesh_input_file.readline()
mesh_secondline = mesh_input_file.readline()
mesh_values = mesh_secondline.split()
my_volume = float(mesh_values[-2])*1e-15 #litres
mesh_input_file.close()

locator_params(nticks=4)

basalmaxes = []
 
T = float(filename[filename.find('tstop')+5:filename.find('_',filename.find('tstop'))])
Duration = float(filename[filename.find('tstop')+5:filename.find('_',filename.find('tstop'))])
times = []
times_nrn = []
if len(MAT) > 0:
  times = [T/(DATA_all.shape[0]-1)*i for i in range(0,DATA_all.shape[0])]
if len(DATANRN_all) > 0:
  times_nrn = DATANRN_all['tvec']
tstop = 500000/1000

f,axs = subplots(6,3)
axarr = [axs[i,j] for j,i in list(itertools.product([0,1,2], [0,1,2,3,4],repeat=1))]+[axs[5,j] for j in range(0,3)]

axnews = []

TCs_memb = []
TCs_memb_nrn = []

for iax in range(0,len(species)+1):
  ipathway = iax/5
  ispecinpathway = iax%5
  mylabel = ylabels[iax]
  if iax < 15:
    axarr[iax].set_position([0.11+0.33*ipathway, 0.79-0.129*(1.6+ispecinpathway), 0.215, 0.129])
    f.text(0.02+0.33*ipathway, 0.79-0.13*(1.6+ispecinpathway)+0.12,chr(68+iax))
  else:
    if iax < 17:
      axarr[iax].set_position([0.11+0.33*(iax-14), 0.82-0.13*0, 0.215, 0.15])
    else:
      axarr[iax].set_position([0.11, 0.82-0.13*0, 0.215, 0.15])
    f.text(0.02+0.33*(iax-15), 0.82-0.13*0+0.12,chr(65+iax-15))
  if iax < len(species):
    if ispecinpathway < 2 and iax < 15:
      axnews.append(f.add_axes([0.18+0.33*ipathway,0.82-0.13*(1.6+ispecinpathway)+0.04,0.11-0.015*(ipathway==1)-0.03*(ipathway==2),0.048]))
    if len(species[iax]) == 1:
      cols = ['#3333FF']
      nrncols = ['#0000AA']
    else:
      cols = ['#3333FF','#00FF00','#FF7777','#FF00FF']
      nrncols = ['#0000BB','#009900','#BB3333','#990099']
    for ispecgroup in range(0,len(species[iax])):
      specgroup = species[iax][ispecgroup]
      if len(MAT) > 0:
        mytimecourse = zeros(DATA_all[:,0].shape[0])
      if len(DATANRN_all) > 0:
        mytimecourse_nrn = zeros(times_nrn.shape[0])
      spectext = ''
      if type(specgroup) is not list:
        specgroup = [specgroup]
      for ispec in range(0,len(specgroup)):
        specfactor = 1.0
        for iexception in range(0,len(factorExceptions)):
          if iax == factorExceptions[iexception][0]:
            specfactor = factorExceptions[iexception][1][ispec]
            print "Using factor "+str(factorExceptions[iexception][1][ispec])+" for iax="+str(iax)+", ispecgroup="+str(ispecgroup)+", ispec="+str(ispec)
            break
        if len(MAT) > 0:
          mytimecourse = mytimecourse + specfactor*DATA_all[:,inddict[specgroup[ispec]]]
        if len(DATANRN_all) > 0:
          if len(specgroup[ispec]) > 24:
            mytimecourse_nrn = mytimecourse_nrn + DATANRN_all[specgroup[ispec][:24]]
          else:
            mytimecourse_nrn = mytimecourse_nrn + DATANRN_all[specgroup[ispec]]
        spectext = spectext+specgroup[ispec]+'+'
        if len(spectext)>30:
          spectext = spectext+specgroup[ispec]+'+\n'

      spectext = spectext[0:-1]
      factor = 1.0
      nrnfactor = 1.0/(1.0/6.022e23/my_volume*1e9)
      if not inparticles[iax]:
        factor = 1.0/6.022e23/my_volume*1e9
        nrnfactor = 1.0
      if ispecinpathway < 2 and iax < 15:
        if len(MAT) > 0:
          axarr[iax].plot([39.5,45.5,45.5,39.5,39.5],[0,0,max(mytimecourse)*factor*1.1,max(mytimecourse)*factor*1.1,0],'r--',linewidth=0.1,dashes=(1,1))
          axnews[ipathway*2+ispecinpathway].plot([x/1000 for x in times[::Nskip]],mytimecourse[::Nskip]*factor,color=cols[ispecgroup],linewidth=0.5,label=titles[iax][ispecgroup] if type(titles[iax]) is list else None)
        if len(DATANRN_all) > 0:
          axnews[ipathway*2+ispecinpathway].plot([x/1000-4000 for x in times_nrn],mytimecourse_nrn*1e6*nrnfactor,'k--',color=nrncols[ispecgroup],linewidth=1.0,dashes=[1,2])
        axnews[ipathway*2+ispecinpathway].set_xlim([39.5,45.5])
      if len(MAT) > 0:
        axarr[iax].plot([x/1000 for x in times[::Nskip]],mytimecourse[::Nskip]*factor,color=cols[ispecgroup],linewidth=0.5,label=titles[iax][ispecgroup] if type(titles[iax]) is list else None)
      if len(DATANRN_all) > 0:
        axarr[iax].plot([x/1000-4000 for x in times_nrn[::Nskip]],mytimecourse_nrn[::Nskip]*1e6*nrnfactor,'k--',color=nrncols[ispecgroup],linewidth=1.0,dashes=[1,2])
      print "mylabel = "+mylabel+", ipathway = "+str(ipathway)+", iax = "+str(iax)+", species="+str(species[iax])
      for iisspecial in range(0,len(ispecies_needed_for_cond)):
        if iax == ispecies_needed_for_cond[iisspecial][0] and ispecgroup == ispecies_needed_for_cond[iisspecial][1]:
          if len(MAT) > 0:
            TCs_memb.append(mytimecourse[::Nskip])
          if len(DATANRN_all) > 0:
            TCs_memb_nrn.append(mytimecourse_nrn[::Nskip]*1e6*nrnfactor/factor)

    print "iax = "+str(iax)
  else:
   if len(MAT) > 0:
    ENhom1_np = (TCs_memb[0] + TCs_memb[2])/4.0 * (TCs_memb[0]-TCs_memb[1])**4/(TCs_memb[0] + TCs_memb[2])**4                                   #Number of complexes times the probability of a complex consisting of 4 non-phos GluR1s
    ENhom1_p = (TCs_memb[0] + TCs_memb[2])/4.0 * (TCs_memb[0]**4 - (TCs_memb[0]-TCs_memb[1])**4)/(TCs_memb[0] + TCs_memb[2])**4                 #The same, but use prob. of having 4 GluR1s, minus the cases where all of them are 
                                                                                                                                                #  non-phos
    ENhom2 = (TCs_memb[0] + TCs_memb[2])/4.0 * (TCs_memb[2]/(TCs_memb[0] + TCs_memb[2]))**4                                                     #Number of complexes times the probability of a complex consisting of 4 GluR2s
    ENhet = (TCs_memb[0] + TCs_memb[2])/4.0 * (1 - (TCs_memb[0]/(TCs_memb[0] + TCs_memb[2]))**4 - (TCs_memb[2]/(TCs_memb[0] + TCs_memb[2]))**4) #Number of complexes times the probability of a complex NOT consisting 
                                                                                                                                                  #  of 4 GluR1s or of 4 GluR2s
    Egtot = ENhom1_np*conds_hom1[0] + ENhom1_p*conds_hom1[1] + ENhom2*conds_hom2 + ENhet*conds_het
    axarr[17].plot([x/1000 for x in times[::Nskip]], Egtot, color=cols[0],linewidth=0.5,label=None)

   if len(DATANRN_all) > 0:
    ENhom1_np = (TCs_memb_nrn[0] + TCs_memb_nrn[2])/4.0 * (TCs_memb_nrn[0]-TCs_memb_nrn[1])**4/(TCs_memb_nrn[0] + TCs_memb_nrn[2])**4                                          
    ENhom1_p = (TCs_memb_nrn[0] + TCs_memb_nrn[2])/4.0 * (TCs_memb_nrn[0]**4 - (TCs_memb_nrn[0]-TCs_memb_nrn[1])**4)/(TCs_memb_nrn[0] + TCs_memb_nrn[2])**4                 
                                                                                                                                               
    ENhom2 = (TCs_memb_nrn[0] + TCs_memb_nrn[2])/4.0 * (TCs_memb_nrn[2]/(TCs_memb_nrn[0] + TCs_memb_nrn[2]))**4                                
    ENhet = (TCs_memb_nrn[0] + TCs_memb_nrn[2])/4.0 * (1 - (TCs_memb_nrn[0]/(TCs_memb_nrn[0] + TCs_memb_nrn[2]))**4 - (TCs_memb_nrn[2]/(TCs_memb_nrn[0] + TCs_memb_nrn[2]))**4)
                                                                                                                                                 
    Egtot = ENhom1_np*conds_hom1[0] + ENhom1_p*conds_hom1[1] + ENhom2*conds_hom2 + ENhet*conds_het
    axarr[17].plot([x/1000-4000 for x in times_nrn], Egtot, 'k--',color=nrncols[0],linewidth=1.0,dashes=[1,2])

  thisylim = axarr[iax].get_ylim()[1]
  if type(ylims[iax]) is list:
    if len(ylims[iax]) > 1:
      axarr[iax].set_ylim(ylims[iax])
  elif not isnan(ylims[iax]):
    axarr[iax].set_ylim([0,ylims[iax]])
  else:
    axarr[iax].set_ylim([0,thisylim])

  try: #this gave errors in one python setup - plot legends without facecolor='None' and wherever you please in that case
    if iax < 15 and iax !=1 and iax != 3 and iax != 4 and iax !=5 and iax != 6 and iax != 11:
      leg = axarr[iax].legend(fontsize=6,facecolor='None',handletextpad=0.1,handlelength=1.0,labelspacing=0.2)
    elif iax == 1:
      leg = axarr[iax].legend(fontsize=6,loc="center",bbox_to_anchor=(0.9,0.8),facecolor='None',handletextpad=0.1,handlelength=1.0,labelspacing=0.2) #upper left
    elif iax == 4:
      leg = axarr[iax].legend(fontsize=6,loc="center",bbox_to_anchor=(0.3,0.75),facecolor='None',handletextpad=0.1,handlelength=1.0,labelspacing=0.2) #upper left
    elif iax == 5:
      leg = axarr[iax].legend(fontsize=6,loc="center",bbox_to_anchor=(0.8,0.95),facecolor='None',handletextpad=0.1,handlelength=1.0,labelspacing=0.2) #upper left
    elif iax == 6:
      leg = axarr[iax].legend(fontsize=6,loc="center",bbox_to_anchor=(0.9,0.8),facecolor='None',handletextpad=0.1,handlelength=1.0,labelspacing=0.2) #upper left
    elif iax == 11:
      leg = axarr[iax].legend(fontsize=6,loc="center",bbox_to_anchor=(0.8,0.85),facecolor='None',handletextpad=0.1,handlelength=1.0,labelspacing=0.2) #upper left
    elif iax == 15:
      leg = axarr[iax].legend(fontsize=6,loc="center",bbox_to_anchor=(0.8,0.2),facecolor='None',handletextpad=0.1,handlelength=1.0,labelspacing=0.2) #lower right      
    elif iax == 16:
      leg = axarr[iax].legend(fontsize=6,loc=5,facecolor='None',handletextpad=0.1,handlelength=1.0,labelspacing=0.2) #right      
    else:
      leg = axarr[iax].legend(fontsize=6,loc=4,facecolor='None',handletextpad=0.1,handlelength=1.0,labelspacing=0.2) #lower right
  except:
    print "Legend not printed where it should go, please see what's wrong with the try block"
    leg = axarr[iax].legend(fontsize=6,handletextpad=0.1,handlelength=1.0,labelspacing=0.2) #lower right
  if len(MAT) > 0:
    try: #this gave errors in one python setup - we can live without it
      leg.get_frame().set_linewidth(0)
    except:
      print "Legend frames not made invisible"
  if ispecinpathway == 0 and iax < 15:
    axarr[iax].set_title(grouptitles[ipathway], fontsize=7)
  if not inparticles[iax]:
    axarr[iax].set_ylabel(mylabel+'\n(nM)', fontsize=6) 
  elif inparticles[iax] == 2:
    axarr[iax].set_ylabel(mylabel+'\n(pS)', fontsize=6) 
  else:
    axarr[iax].set_ylabel(mylabel+'\n(part.)', fontsize=6)
  if iax % 5 < 4 and iax < 15:
    axarr[iax].set_xticklabels([])
  for tick in axarr[iax].yaxis.get_major_ticks()+axarr[iax].xaxis.get_major_ticks():
    tick.label.set_fontsize(6)
    tick.label2.set_fontsize(6)
  xticks = axarr[iax].get_xticks()
  if len(xticks) > 6:
    axarr[iax].set_xticks(xticks[::2])
  if len(yticks_all[iax]) == 0:
    yticks = axarr[iax].get_yticks()
    if len(yticks) > 4:
      axarr[iax].set_yticks(yticks[::2])
  else:
    axarr[iax].set_yticks(yticks_all[iax])
  axarr[iax].set_xlabel('time (s)', fontsize=6)
  if type(ylims[iax]) is not list and not isnan(ylims[iax]):
    axarr[iax].text(Duration*0.6,ylims[iax]*0.82,titles[iax],fontsize=7.0)
  elif type(ylims[iax]) is list:
    axarr[iax].text(Duration*0.6,ylims[iax][1]*0.82+ylims[iax][0]*0.18,titles[iax],fontsize=7.0)    
  elif not isnan(ylims[iax]):
    axarr[iax].text(Duration*0.6,ylims[iax][-1]*0.82,titles[iax],fontsize=7.0)
  if ispecinpathway < 2 and iax < 15:
    for tick in axnews[ipathway*2+ispecinpathway].xaxis.get_major_ticks() + axnews[ipathway*2+ispecinpathway].yaxis.get_major_ticks():
      tick.label.set_fontsize(5)
    for axis in ['top','bottom','left','right']:
      axnews[ipathway*2+ispecinpathway].spines[axis].set_linewidth(0.5)
    axnews[ipathway*2+ispecinpathway].spines['top'].set_visible(False)
    axnews[ipathway*2+ispecinpathway].spines['right'].set_visible(False)
    axnews[ipathway*2+ispecinpathway].tick_params(top='off', right='off')
    xticks = axnews[ipathway*2+ispecinpathway].get_xticks()
    if len(xticks) > 6:
      axnews[ipathway*2+ispecinpathway].set_xticks(xticks[::2])
    yticks = axnews[ipathway*2+ispecinpathway].get_yticks()
    if len(yticks) > 4:
      axnews[ipathway*2+ispecinpathway].set_yticks(yticks[::2])
      
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.5)
  axarr[iax].spines['top'].set_visible(False)
  axarr[iax].spines['right'].set_visible(False)
  axarr[iax].tick_params(top='off', right='off')
  axarr[iax].set_xlim([0,tstop])

for ax in axarr+axnews:
  ax.tick_params(axis='both',direction='in')


if len(MAT) > 0 or len(DATANRN_all) > 0:
  f.savefig('fig3A.eps')
else:
  print "No files found, fig3A.eps not saved"




##################### Fig 4, LFS panels #######################


close("all")

protocol = 'LFS'
ylims[15] = [-10,150]
ylims[16] = [0,150]

filename = '../NeuroRD/tstop500000_tol0.01_onset40000.0_n900_freq5.0_dur3.0_flux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(AChflux)+'_Ntrains1_trainT100000.0_8seeds.mat'
MAT = {}
if exists(filename):
  MAT = scipy.io.loadmat(filename)
  print "loaded "+filename

  DATA_all = MAT['DATA']
  header_strs = MAT['headers']
  for i in range(0,len(header_strs)):
    first_space = header_strs[i].find(' ')
    if first_space > -1:
      header_strs[i] = header_strs[i][:first_space]

  inddict = {}
  for iheader in range(4,len(header_strs)):
    inddict[header_strs[iheader]] = iheader-4
else:
  print filename+' not found. Variable MAT will remain empty and NeuroRD results will not be plotted'

filename_nrn = 'nrn_noninterp_tstop5000000_tol1e-06_onset4040000.0_n900_freq5.0_dur3.0_flux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(AChflux)+'_Ntrains1_trainT100000.0.mat'
DATANRN_all = {}
if exists(filename_nrn):
  DATANRN_all_all = scipy.io.loadmat(filename_nrn)
  print "loaded "+filename_nrn
  for ikey in range(0,len(DATANRN_all_all['headers'])):
    mykey = DATANRN_all_all['headers'][ikey][0:DATANRN_all_all['headers'][ikey].find(' ')]
    DATANRN_all[mykey] = DATANRN_all_all['DATA'][ikey]
else:
  print filename_nrn+' not found. Variable DATANRN_all will remain empty and NEURON results will not be plotted'
nrncols_all = [ ['#0000AA','#00AA00','#AA0000','#770077'], ['#0000AA','#770077'] ]

locator_params(nticks=4)

times = []
times_nrn = []
if len(MAT) > 0:
  times = [T/(DATA_all.shape[0]-1)*i for i in range(0,DATA_all.shape[0])]
if len(DATANRN_all) > 0:
  times_nrn = DATANRN_all['tvec']

f,axs = subplots(6,3)
axarr = [axs[i,j] for j,i in list(itertools.product([0,1,2], [0,1,2,3,4],repeat=1))]+[axs[5,j] for j in range(0,3)]

axnews = []

TCs_memb = []
TCs_memb_nrn = []

for iax in range(0,len(species)+1):
  ipathway = iax/5
  ispecinpathway = iax%5
  mylabel = ylabels[iax]
  if iax < 15:
    axarr[iax].set_position([0.11+0.33*ipathway, 0.79-0.129*(1.6+ispecinpathway), 0.215, 0.129])
    axarr[iax].set_visible(False)
  else:
    if iax < 17:
      axarr[iax].set_position([0.11+0.33*(iax-14), 0.82-0.13*0, 0.215, 0.15])
    else:
      axarr[iax].set_position([0.11, 0.82-0.13*0, 0.215, 0.15])
    f.text(0.02+0.33*(iax-15), 0.82-0.13*0+0.12,chr(65+18+iax-15))
  if iax < len(species):
    if ispecinpathway < 2 and iax < 15:
      axnews.append(f.add_axes([0.18+0.33*ipathway,0.82-0.13*(1.6+ispecinpathway)+0.04,0.1,0.048]))
    if len(species[iax]) == 1:
      cols = ['#3333FF']
      nrncols = ['#0000AA']
    else:
      cols = ['#3333FF','#00FF00','#FF7777','#FF00FF']
      nrncols = ['#0000BB','#009900','#BB3333','#990099']
    for ispecgroup in range(0,len(species[iax])):
      specgroup = species[iax][ispecgroup]
      if len(MAT) > 0:
        mytimecourse = zeros(DATA_all[:,0].shape[0])
      if len(DATANRN_all) > 0:
        mytimecourse_nrn = zeros(times_nrn.shape[0])
      spectext = ''
      if type(specgroup) is not list:
        specgroup = [specgroup]
      for ispec in range(0,len(specgroup)):
        specfactor = 1.0
        for iexception in range(0,len(factorExceptions)):
          if iax == factorExceptions[iexception][0]:
            specfactor = factorExceptions[iexception][1][ispec]
            print "Using factor "+str(factorExceptions[iexception][1][ispec])+" for iax="+str(iax)+", ispecgroup="+str(ispecgroup)+", ispec="+str(ispec)
            break
        if len(MAT) > 0:
          mytimecourse = mytimecourse + specfactor*DATA_all[:,inddict[specgroup[ispec]]]
        if len(DATANRN_all) > 0:
          if len(specgroup[ispec]) > 24:
            mytimecourse_nrn = mytimecourse_nrn + DATANRN_all[specgroup[ispec][:24]]
          else:
            mytimecourse_nrn = mytimecourse_nrn + DATANRN_all[specgroup[ispec]]
        spectext = spectext+specgroup[ispec]+'+'
        if len(spectext)>30:
          spectext = spectext+specgroup[ispec]+'+\n'

      spectext = spectext[0:-1]
      factor = 1.0
      nrnfactor = 1.0/(1.0/6.022e23/my_volume*1e9)
      if not inparticles[iax]:
        factor = 1.0/6.022e23/my_volume*1e9
        nrnfactor = 1.0
      if ispecinpathway < 2 and iax < 15:
        if len(MAT) > 0:
          axarr[iax].plot([39.5,45.5,45.5,39.5,39.5],[0,0,max(mytimecourse)*factor*1.1,max(mytimecourse)*factor*1.1,0],'r--',linewidth=0.1,dashes=(1,1))
          axnews[ipathway*2+ispecinpathway].plot([x/1000 for x in times[::Nskip]],mytimecourse[::Nskip]*factor,color=cols[ispecgroup],linewidth=0.5,label=titles[iax][ispecgroup] if type(titles[iax]) is list else None)
        if len(DATANRN_all) > 0:
          axnews[ipathway*2+ispecinpathway].plot([x/1000-4000 for x in times_nrn],mytimecourse_nrn*1e6*nrnfactor,'k--',color=nrncols[ispecgroup],linewidth=1.0,dashes=[1,2])
        axnews[ipathway*2+ispecinpathway].set_xlim([39.5,45.5])
      if len(MAT) > 0:
        axarr[iax].plot([x/1000 for x in times[::Nskip]],mytimecourse[::Nskip]*factor,color=cols[ispecgroup],linewidth=0.5,label=titles[iax][ispecgroup] if type(titles[iax]) is list else None)
      if len(DATANRN_all) > 0:
        axarr[iax].plot([x/1000-4000 for x in times_nrn[::Nskip]],mytimecourse_nrn[::Nskip]*1e6*nrnfactor,'k--',color=nrncols[ispecgroup],linewidth=1.0,dashes=[1,2])
      print "mylabel = "+mylabel+", ipathway = "+str(ipathway)+", iax = "+str(iax)+", species="+str(species[iax])
      for iisspecial in range(0,len(ispecies_needed_for_cond)):
        if iax == ispecies_needed_for_cond[iisspecial][0] and ispecgroup == ispecies_needed_for_cond[iisspecial][1]:
          if len(MAT) > 0:
            TCs_memb.append(mytimecourse[::Nskip])
          if len(DATANRN_all) > 0:
            TCs_memb_nrn.append(mytimecourse_nrn[::Nskip]*1e6*nrnfactor/factor)

    print "iax = "+str(iax)
  else:
   if len(MAT) > 0:
    ENhom1_np = (TCs_memb[0] + TCs_memb[2])/4.0 * (TCs_memb[0]-TCs_memb[1])**4/(TCs_memb[0] + TCs_memb[2])**4                                   #Number of complexes times the probability of a complex consisting of 4 non-phos GluR1s
    ENhom1_p = (TCs_memb[0] + TCs_memb[2])/4.0 * (TCs_memb[0]**4 - (TCs_memb[0]-TCs_memb[1])**4)/(TCs_memb[0] + TCs_memb[2])**4                 #The same, but use prob. of having 4 GluR1s, minus the cases where all of them are 
                                                                                                                                                #  non-phos
    ENhom2 = (TCs_memb[0] + TCs_memb[2])/4.0 * (TCs_memb[2]/(TCs_memb[0] + TCs_memb[2]))**4                                                     #Number of complexes times the probability of a complex consisting of 4 GluR2s
    ENhet = (TCs_memb[0] + TCs_memb[2])/4.0 * (1 - (TCs_memb[0]/(TCs_memb[0] + TCs_memb[2]))**4 - (TCs_memb[2]/(TCs_memb[0] + TCs_memb[2]))**4) #Number of complexes times the probability of a complex NOT consisting 
                                                                                                                                                  #  of 4 GluR1s or of 4 GluR2s
    Egtot = ENhom1_np*conds_hom1[0] + ENhom1_p*conds_hom1[1] + ENhom2*conds_hom2 + ENhet*conds_het
    axarr[17].plot([x/1000 for x in times[::Nskip]], Egtot, color=cols[0],linewidth=0.5,label=None)

   if len(DATANRN_all) > 0:
    ENhom1_np = (TCs_memb_nrn[0] + TCs_memb_nrn[2])/4.0 * (TCs_memb_nrn[0]-TCs_memb_nrn[1])**4/(TCs_memb_nrn[0] + TCs_memb_nrn[2])**4                                          
    ENhom1_p = (TCs_memb_nrn[0] + TCs_memb_nrn[2])/4.0 * (TCs_memb_nrn[0]**4 - (TCs_memb_nrn[0]-TCs_memb_nrn[1])**4)/(TCs_memb_nrn[0] + TCs_memb_nrn[2])**4                 
                                                                                                                                               
    ENhom2 = (TCs_memb_nrn[0] + TCs_memb_nrn[2])/4.0 * (TCs_memb_nrn[2]/(TCs_memb_nrn[0] + TCs_memb_nrn[2]))**4                                
    ENhet = (TCs_memb_nrn[0] + TCs_memb_nrn[2])/4.0 * (1 - (TCs_memb_nrn[0]/(TCs_memb_nrn[0] + TCs_memb_nrn[2]))**4 - (TCs_memb_nrn[2]/(TCs_memb_nrn[0] + TCs_memb_nrn[2]))**4)
                                                                                                                                                 
    Egtot = ENhom1_np*conds_hom1[0] + ENhom1_p*conds_hom1[1] + ENhom2*conds_hom2 + ENhet*conds_het
    axarr[17].plot([x/1000-4000 for x in times_nrn], Egtot, 'k--',color=nrncols[0],linewidth=1.0,dashes=[1,2])
    print "g_rel: "+str(Egtot[-1]/Egtot[[i for i in range(0,len(times_nrn)) if times_nrn[i] < 4040000][-1]])


  thisylim = axarr[iax].get_ylim()[1]
  if type(ylims[iax]) is list:
    if len(ylims[iax]) > 1:
      axarr[iax].set_ylim(ylims[iax])
  elif not isnan(ylims[iax]):
    axarr[iax].set_ylim([0,ylims[iax]])
  else:
    axarr[iax].set_ylim([0,thisylim])

  if iax < 15 and iax != 3 and iax != 4:
    leg = axarr[iax].legend(fontsize=6)
  elif iax == 4:
    leg = axarr[iax].legend(fontsize=6,loc=2) #upper left
  elif iax == 15:
    leg = axarr[iax].legend(fontsize=6,loc=1) #upper right      
  elif iax == 16:
    leg = axarr[iax].legend(fontsize=6,loc=1) #upper right      
  else:
    leg = axarr[iax].legend(fontsize=6,loc=1) #upper right
  if len(MAT) > 0:
    try: #this gave errors in one python setup - we can live without it
      leg.get_frame().set_linewidth(0)
    except:
      print "Legend frames not made invisible"
  if ispecinpathway == 0 and iax < 15:
    axarr[iax].set_title(grouptitles[ipathway], fontsize=7)
  if not inparticles[iax]:
    axarr[iax].set_ylabel(mylabel+'\n(nM)', fontsize=6) 
  elif inparticles[iax] == 2:
    axarr[iax].set_ylabel(mylabel+'\n(pS)', fontsize=6) 
  else:
    axarr[iax].set_ylabel(mylabel+'\n(part.)', fontsize=6)
  if iax % 5 < 4 and iax < 15:
    axarr[iax].set_xticklabels([])
  for tick in axarr[iax].yaxis.get_major_ticks()+axarr[iax].xaxis.get_major_ticks():
    tick.label.set_fontsize(6)
    tick.label2.set_fontsize(6)
  xticks = axarr[iax].get_xticks()
  if len(xticks) > 6:
    axarr[iax].set_xticks(xticks[::2])
  if len(yticks_all[iax]) == 0:
    yticks = axarr[iax].get_yticks()
    if len(yticks) > 4:
      axarr[iax].set_yticks(yticks[::2])
  else:
    axarr[iax].set_yticks(yticks_all[iax])
  axarr[iax].set_xlabel('time (s)', fontsize=6)
  if type(ylims[iax]) is not list and not isnan(ylims[iax]):
    axarr[iax].text(Duration*0.6,ylims[iax]*0.82,titles[iax],fontsize=7.0)
  elif type(ylims[iax]) is list:
    axarr[iax].text(Duration*0.6,ylims[iax][1]*0.82+ylims[iax][0]*0.18,titles[iax],fontsize=7.0)    
  elif not isnan(ylims[iax]):
    axarr[iax].text(Duration*0.6,ylims[iax][-1]*0.82,titles[iax],fontsize=7.0)
  if ispecinpathway < 2 and iax < 15:
    for tick in axnews[ipathway*2+ispecinpathway].xaxis.get_major_ticks() + axnews[ipathway*2+ispecinpathway].yaxis.get_major_ticks():
      tick.label.set_fontsize(5)
    for axis in ['top','bottom','left','right']:
      axnews[ipathway*2+ispecinpathway].spines[axis].set_linewidth(0.5)
    axnews[ipathway*2+ispecinpathway].spines['top'].set_visible(False)
    axnews[ipathway*2+ispecinpathway].spines['right'].set_visible(False)
    axnews[ipathway*2+ispecinpathway].tick_params(top='off', right='off')
    xticks = axnews[ipathway*2+ispecinpathway].get_xticks()
    if len(xticks) > 6:
      axnews[ipathway*2+ispecinpathway].set_xticks(xticks[::2])
    yticks = axnews[ipathway*2+ispecinpathway].get_yticks()
    if len(yticks) > 4:
      axnews[ipathway*2+ispecinpathway].set_yticks(yticks[::2])
      
  for tick in axarr[iax].xaxis.get_major_ticks() + axarr[iax].yaxis.get_major_ticks():
    tick.label.set_fontsize(5)
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.5)
  axarr[iax].spines['top'].set_visible(False)
  axarr[iax].spines['right'].set_visible(False)
  axarr[iax].tick_params(top='off', right='off')
  axarr[iax].set_xlim([0,tstop])

for ax in axarr+axnews:
  ax.tick_params(axis='both',direction='in')
for ax in axnews:
  ax.set_visible(False)


if len(MAT) > 0 or len(DATANRN_all) > 0:
  f.savefig('fig3B.eps')
else:
  print "No files found, fig3B.eps not saved"


