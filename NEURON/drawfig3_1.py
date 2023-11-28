#drawfig3_1.py: Draws the synaptic conductance for 4xHFS and LFS when GluR1 or GluR2 subunits are absent using the normal tetramer formation rule, or 
#               when GluR1 or GluR2 subunits are present in proportion 35:65 using the alternative dimer-of-like-dimers tetramer formation rule.
#Tuomo Maki-Marttunen, 2019-2020
#CC BY 4.0
import matplotlib
matplotlib.use('Agg')
from pylab import *
import scipy.io
import sys
import itertools
from os.path import exists
import calcconds
import calcconds_dimerdimer


close("all")

species = [ ['Ca'], ['CaMCa4'], [['Complex', 'pComplex']], [['CKp','CKpCaMCa4']], [['GluR1_S831', 'GluR1_S845_S831', 'GluR1_S831_PKAc', 'GluR1_S845_S831_PP1', 'GluR1_S831_PP1', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'],
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

titles = ['[Ca$^{2+}$]','[CaMCa4]','[CaMKII complex]','[CaMKII phos.]',['[GluR1 S831p]', '[GluR1 2x phosph.]'],
          '[$\\beta$-adr. ligand]',['[act. Gs]','[act. Gi]'],'[cAMP]','[act. PKA]',['[GluR1 S845p]', '[GluR1 2x phosph.]'],
          ['[glutamate]','[ACh]'],'[act. Gq]','[DAG]',['[PKC tr.]','[PKC pers.]'],'[GluR2 S880p]', ['[GluR1 at memb.]', '[GluR1 S831p at memb.]'], '[GluR2 at memb.]', 'max. cond.']
yticks_all = [[], [], [0,400], [0,2000], [],    [], [0,400,800], [], [], [0,40,80],    [], [], [], [], [],    [0,100,200], [0,100,200], [] ]

grouptitles = ['CaMKII pathway','Gs/Gi-cAMP-PKA pathway','Gq-PLC-PKC pathway']
protocols = ['HFS','LFS']
Caflux = 1900.0
Lflux = 10.0
Gluflux = 20.0
AChflux = 20.0

myratio = 0.35
if len(sys.argv) > 1:
  myratio = float(sys.argv[1])

ratios = [0.0, 1.0, myratio]
dimers = [False, False, True]
titles = ['3.1A','3.1B','3.2A']

mesh_input_file = open('mesh_general.out','r')
mesh_firstline = mesh_input_file.readline()
mesh_secondline = mesh_input_file.readline()
mesh_values = mesh_secondline.split()
my_volume = float(mesh_values[-2])*1e-15 #litres
mesh_input_file.close()

locator_params(nticks=4)
cols = ['#0000FF', '#BBBB22']

f,axarr = subplots(1,3)

for iprotocol in [0,1]:
  protocol = protocols[iprotocol]
  for iratio in range(0,len(ratios)):
    axarr[iratio].set_position([0.06+0.32*iratio, 0.63, 0.25, 0.3])
    ratio = ratios[iratio]
    Nskip = 1
    if iprotocol == 0:
      filename_nrn = 'nrn_tstop5000000_tol1e-06_GluR1-GluR1_memb-GluR2-GluR2_membx'+str(2*ratio)+'-'+str(2*ratio)+'-'+str(2-2*ratio)+'-'+str(2-2*ratio)+'_onset4040000.0_n100_freq100.0_dur3.0_flux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(AChflux)+'_Ntrains4_trainT4000.0.mat'
    else:
      filename_nrn = 'nrn_tstop5000000_tol1e-06_GluR1-GluR1_memb-GluR2-GluR2_membx'+str(2*ratio)+'-'+str(2*ratio)+'-'+str(2-2*ratio)+'-'+str(2-2*ratio)+'_onset4040000.0_n900_freq5.0_dur3.0_flux'+str(Caflux)+'_Lflux'+str(Lflux)+'_Gluflux'+str(Gluflux)+'_AChflux'+str(AChflux)+'_Ntrains1_trainT100000.0.mat'
    DATANRN_all = {}
    if exists(filename_nrn):
      DATANRN_all_all = scipy.io.loadmat(filename_nrn)
      print "loaded "+filename_nrn
      for ikey in range(0,len(DATANRN_all_all['headers'])):
        mykey = DATANRN_all_all['headers'][ikey][0:DATANRN_all_all['headers'][ikey].find(' ')]
        DATANRN_all[mykey] = DATANRN_all_all['DATA'][ikey]
    else:
      print filename_nrn+" not found"

    basalmaxes = []
 
    times_nrn = []
    if len(DATANRN_all) > 0:
      times_nrn = DATANRN_all['tvec']
    tstop = 500000/1000

    if dimers[iratio]:
      DATA = calcconds_dimerdimer.calcconds_nrn(filename_nrn)
    else:
      DATA = calcconds.calcconds_nrn(filename_nrn)
    Egtot = DATA[0]
    axarr[iratio].plot([x/1000-4000 for x in times_nrn], Egtot, 'k-',color=cols[iprotocol],linewidth=1.0,label=protocol)
      
    if iprotocol == 1:
      for tick in axarr[iratio].xaxis.get_major_ticks() + axarr[iratio].yaxis.get_major_ticks():
        tick.label.set_fontsize(5)
      for axis in ['top','bottom','left','right']:
        axarr[iratio].spines[axis].set_linewidth(0.5)
      axarr[iratio].spines['top'].set_visible(False)
      axarr[iratio].spines['right'].set_visible(False)
      axarr[iratio].tick_params(top='off', right='off')
      axarr[iratio].set_xlim([0,tstop])
      axarr[iratio].set_xlabel('time (s)',fontsize=8)
      axarr[iratio].set_ylabel('pS',fontsize=8)
      pos = axarr[iratio].get_position()
      f.text(pos.x0 - 0.05, pos.y1 + 0.02, titles[iratio], fontsize=14)

      axarr[iratio].legend(fontsize=7)
axarr[0].set_title('0% GluR1, 100% GluR2', fontsize=8)
axarr[1].set_title('100% GluR1, 0% GluR2', fontsize=8)
axarr[2].set_title(str(ratio*100)+'% GluR1, '+str(100-ratio*100)+'% GluR2,\ndimers of dimers', fontsize=8)
f.savefig('fig3_1.eps')



