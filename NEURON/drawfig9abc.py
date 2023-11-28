#drawfig9abc.py: Draws the model outputs for the parameters fitted to cortical LTP/LTD data.
#Tuomo Maki-Marttunen, 2019-2020
#CC BY 4.0
import matplotlib
matplotlib.use('Agg')
from pylab import *
import mytools
import pickle
import protocols_many
import protocols_many_78withoutCK
import protocols_many_78withoutCK_1withCK
from os.path import exists
import time
import scipy.stats
def plotmybox(ax,ys,x=0,w=0.5,lw=0.5,col='#000000'): #ys: vector of 5 elements: min, prc-25, median, prc-75, max
  ax.plot([x-w,x+w,x,x,x-w,x-w,x+w,x+w,x,nan,x-w,x+w,nan,x,x,x-w,x+w],[ys[0],ys[0],ys[0],ys[1],ys[1],ys[3],ys[3],ys[1],ys[1],nan,ys[2],ys[2],nan,ys[3],ys[4],ys[4],ys[4]],'k-',linewidth=lw,color=col)

VARIABLES = [["Caflux",0,5000], #the upper limit of Caflux will be changed according to imeas
             ["Lflux",0.0,5.0],
             ["Gluflux",0,200],
             ["GluR1_ratio",0.0,1.0],
             ["IC_MGluRM1GqPLC",0.0,2.0],
             ["IC_RGsAC1AC8",0.0,2.0],
             ["IC_CaMCK",0.0,2.0],
             ["IC_NCX",0.0,2.0],
             ["IC_PKC",0.0,5.0],
             ["IC_PKA",0.0,2.0],
             ["IC_PP1PP2B",0.0,2.0],
             ["IC_PDE1PDE4",0.0,2.0],
             ]

Caflux_limits = [20000, 20000, 13000, 13000, 50000, 10000, 40000, 20000, 20000, 16000, 16000] #Planned so that [Ca flux]*T_total_input is around 2e6, but for imeas=7,8, two different protocols used - something in the middle taken

imeass = [0,1,2,3,4,5,6,7,8,9,10,7,8]
captions = ['EC-1','EC-2','PFC-1','PFC-2','BC','ACC','PFC-3','VC-1','VC-2','AC-1','AC-2']
myseeds = [1,1,1,1,1,1,1,30,30,1,1,1,1]
exts = ['fewer','fewer','fewer','fewer','fewer','fewer','fewer','manyb','manyb','fewer','fewer','fewer','fewer',]
rundexts = ['fewer','fewerCK1imeas','fewer','fewer','fewer','fewer','fewer','manyb','manyb','fewer','fewer','fewer','fewer',]

isamps = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
Nsamps = [1000,1000,1000,1000,1000,1000,1000,500,500,1000,1000,1000,1000]

Measurement_protocol = protocols_many.get_measurement_protocol()
Measurement_protocol_78withoutCK = protocols_many_78withoutCK.get_measurement_protocol()
Measurement_protocol_78withoutCK_1withCK = protocols_many_78withoutCK_1withCK.get_measurement_protocol()

maxerr = 1.0
maxcaerr = 0.0

def clamp(x): 
  return max(0, min(int(256*x), 255))

def col2hexcol(rgb,brightness=1.0,dim=0.0):
  meanrgb = mean(rgb[0:3])
  r = (rgb[0]*(1-dim)+meanrgb*dim)*brightness
  g = (rgb[1]*(1-dim)+meanrgb*dim)*brightness
  b = (rgb[2]*(1-dim)+meanrgb*dim)*brightness
  return "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))  
try:
  cmap = matplotlib.cm.get_cmap('viridis')
  colors = [col2hexcol(cmap(0.31*i)) for i in range(0,4)]
except:
  colors = ['#440154', '#33628d', '#26ad81', '#d3e21b']
  

f,axarr = subplots(16,1)
for iax in range(0,16):
  for axis in ['top','bottom','left','right']:
    axarr[iax].spines[axis].set_linewidth(0.5)
  axarr[iax].spines['top'].set_visible(False)
  axarr[iax].spines['right'].set_visible(False)
  axarr[iax].tick_params(width=0.2,length=2.0,labelsize=5)
  
for iax in range(0,7):
  axarr[iax].set_position([0.05, 0.86-0.1*iax, 0.11, 0.08])
for iax in range(0,2):
  axarr[9+iax].set_position([0.05, 0.16-0.1*iax, 0.11, 0.08])
  axarr[7+iax].set_position([0.23, 0.75-0.19*iax, 0.11, 0.17])
  axarr[11+iax].set_position([0.41, 0.75-0.19*iax, 0.11, 0.17])

labels78 = ['control (HFS)','CaMKII blocked (HFS)', 'control (LFS)','CaMKII blocked (LFS)']
labels78_B = ['control (HFS)','control (LFS)']

for iimeas in range(0,13):
  imeas = imeass[iimeas]

  MeasurementsAll =      Measurement_protocol_78withoutCK[0]
  Measurements_stdsAll = Measurement_protocol_78withoutCK[7]
  if iimeas == 7 or iimeas == 8:
    MeasurementsAll =      Measurement_protocol[0]
    Measurements_stdsAll = Measurement_protocol[7]
    print "iimeas = "+str(iimeas)

  Measurements =     MeasurementsAll[imeas]
  targetTs =         Measurements[1]
  targetVals =       Measurements[2]
  Measurement_stds = Measurements_stdsAll[imeas]
  OBJECTIVES = ['f'+str(i) for i in range(0,len(Measurements[0])+1)]
  VARIABLES[0][2] = Caflux_limits[imeas]

  mylw = 0.6 
  myms = 1.1 

  finalThrsAbsolute = [maxerr*nansum(Measurement_stds[i]) for i in range(0,len(Measurement_stds))]+[maxcaerr]
  goodparams = []
  gooddata = []
  IDs = []
  coeffs = rand(len(VARIABLES),)
  
  Nall = 0
  filename = 'fitfiles/'+exts[iimeas]+str(imeas)+'_seed'+str(myseeds[iimeas])+'_N'+str(Nsamps[iimeas])
  fitnesses = []
  for gen in range(24,0,-1):
    if exists(filename+'_tmp'+str(gen)+'.sav'):
      gensdone = gen
      print 'loading '+filename+'_tmp'+str(gen)+'.sav'
      unpicklefile = open(filename+'_tmp'+str(gen)+'.sav', 'r')
      unpickledlist = pickle.load(unpicklefile)
      unpicklefile.close()
      params_all = unpickledlist[0]
      columns = unpickledlist[1]
      for iparam in range(0,params_all.shape[0]):
        Nall = Nall + 1
        isbelowMed = True
        fitness = 0
        for iobj in range(0,len(OBJECTIVES)):
          if finalThrsAbsolute[iobj] > 0:
            fitness = fitness + params_all[iparam,len(VARIABLES)+iobj]/finalThrsAbsolute[iobj]
          if params_all[iparam,len(VARIABLES)+iobj] > finalThrsAbsolute[iobj]:
            isbelowMed = False
            break
        if isbelowMed:
          myID = sum([coeffs[i]*params_all[iparam,i] for i in range(0,len(VARIABLES))])
          if myID not in IDs:
            gooddata.append(params_all[iparam,:])
            goodparams.append([(params_all[iparam,i] - VARIABLES[i][1])/(VARIABLES[i][2] - VARIABLES[i][1]) for i in range(0,len(VARIABLES))])
            IDs.append(myID)
            fitnesses.append(fitness)

  myord = [i[0] for i in sorted(enumerate(fitnesses), key=lambda x:x[1])]
  filename = 'fitfiles/rungiven_'+rundexts[iimeas]+str(imeas)+'_seed'+str(myseeds[iimeas])+'_N'+str(Nsamps[iimeas])+'_maxerr1.0_maxcaerr0.0_'+str(isamps[iimeas])+'.sav'
  print filename

  if not exists(filename):
    print filename+' does not exist'
    time.sleep(0.02)
    continue
  print 'loading '+filename
  unpicklefile = open(filename,'r')
  unpickledlist = pickle.load(unpicklefile)
  unpicklefile.close()
  mydict = unpickledlist[0]
  A = unpickledlist[1]

  timesAll = A[0]
  timeCoursesAll = A[1]
  maxValsAll = A[2]

  MeasurementsAll =      Measurement_protocol_78withoutCK_1withCK[0]
  Measurements_stdsAll = Measurement_protocol_78withoutCK_1withCK[7]
  if iimeas == 7 or iimeas == 8:
    MeasurementsAll =      Measurement_protocol[0]
    Measurements_stdsAll = Measurement_protocol[7]
  if iimeas == 1: # Do not plot the data for the neglected objective function
    MeasurementsAll =      Measurement_protocol[0]
    Measurements_stdsAll = Measurement_protocol[7]
  #All data sets have same checking protocol file

  Measurements =     MeasurementsAll[imeas]
  targetTs =         Measurements[1]
  targetVals =       Measurements[2]
  Measurement_stds = Measurements_stdsAll[imeas]
  OBJECTIVES = ['f'+str(i) for i in range(0,len(Measurements[0])+1)]

  errSum = 0
  for iobj in range(0,len(targetVals)):
    mycolor=colors[iobj] if iimeas != 11 and iimeas != 12 else colors[2*iobj]
    errThis = 0
    for itarget in range(0,len(targetTs)):
      itime = argmin(abs(timesAll[iobj]-targetTs[itarget]))
      myval = timeCoursesAll[iobj][itime]/timeCoursesAll[iobj][0]
      if isnan(targetVals[iobj][itarget]):
        continue
      errThis = errThis + abs(targetVals[iobj][itarget] - myval)
    errSum = errSum + errThis
    print "imeas = "+str(imeas)+", ext="+exts[iimeas]+",iobj = "+str(iobj)+", errSum = "+str(errSum)
    if iimeas == 7 or iimeas == 8:
      axarr[iimeas].plot([3e6]+timesAll[iobj].tolist(),[1.0]+[timeCoursesAll[iobj][k]/timeCoursesAll[iobj][0] for k in range(0,len(timesAll[iobj]))],'b-',color=mycolor,lw=mylw,label=labels78[iobj])
    elif iimeas == 11 or iimeas == 12:
      axarr[iimeas].plot([3e6]+timesAll[iobj].tolist(),[1.0]+[timeCoursesAll[iobj][k]/timeCoursesAll[iobj][0] for k in range(0,len(timesAll[iobj]))],'b-',color=mycolor,lw=mylw,label=labels78[iobj])
    else:
      axarr[iimeas].plot([3e6]+timesAll[iobj].tolist(),[1.0]+[timeCoursesAll[iobj][k]/timeCoursesAll[iobj][0] for k in range(0,len(timesAll[iobj]))],'b-',color=mycolor,lw=mylw,label='{:.3f}'.format(errThis)+', '+'{:.4f}'.format(max(A[2])))
    axarr[iimeas].plot([x-10000*(iobj-1) for x in targetTs],Measurements[2][iobj],'r.',color=mycolor,mew=myms,ms=myms)
    if imeas == 7 or imeas == 8:
      for itarget in range(0,len(targetTs)):
        axarr[iimeas].plot([targetTs[itarget]-10000*(iobj-1),targetTs[itarget]-10000*(iobj-1)],[Measurements[2][iobj][itarget]-Measurement_stds[iobj][itarget],Measurements[2][iobj][itarget]+Measurement_stds[iobj][itarget]],'r-',color=mycolor,lw=mylw)

legax = mytools.mylegend(f,[0.04,0.945,0.15,0.037],['b-','b-','b-'],['1st experiment', '2nd experiment', '3rd experiment'],1,2,0.5,0.35,colors[0:3],dashes=[],linewidths=[],myfontsize=4.5)
for axis in ['top','bottom','left','right']:
  legax.spines[axis].set_visible(False) #linewidth(0.5)
legax = mytools.mylegend(f,[0.22,0.93,0.15,0.05],['b-','b-','b-','b-'],labels78,1,2,0.5,0.35,colors,dashes=[],linewidths=[],myfontsize=4.5)
for axis in ['top','bottom','left','right']:
  legax.spines[axis].set_visible(False) #linewidth(0.5)
legax = mytools.mylegend(f,[0.4,0.93,0.1,0.027],['b-','b-'],labels78_B,1,2,0.5,0.35,[colors[i] for i in [0,2]],dashes=[],linewidths=[],myfontsize=4.5)
f.text(0.005,0.95,'A',fontsize=12)
f.text(0.175,0.95,'B',fontsize=12)
f.text(0.375,0.95,'C',fontsize=12)

for axis in ['top','bottom','left','right']:
  legax.spines[axis].set_visible(False) #.set_linewidth(0.5)

ylims = [[0.98,1.7],[0.98,1.7],[0.98,2.05],[0.98,1.74],[0.98,1.6],[0.9,1.6],[0.98,1.44],[0.7,1.65],[0.7,1.45],[0.51,2.15],[0.51,2.15],[0.4,1.4],[0.4,1.4]]
yticks = [[1.0,1.25,1.5],[1.0,1.25,1.5],[1.0,1.4,1.8],[1.0,1.3,1.6],[1.0,1.25,1.5],[1.0,1.25,1.5],[1.0,1.25],[0.75,1.0,1.25],[0.75,1.0,1.25],[0.6,1.0,1.4,1.8],[0.6,1.0,1.4,1.8]]
for iimeas in range(0,13):
  imeas = imeass[iimeas]
  axarr[iimeas].set_xlim([3.5e6,5.3e6])
  axarr[iimeas].set_xticks([3.44e6, 4.04e6, 4.64e6, 5.24e6])
  if imeas in [8,10]:
    axarr[iimeas].set_xticklabels(['-10', '0', '10', '20'],fontsize=5)
  else:
    axarr[iimeas].set_xticklabels([])

  axarr[iimeas].set_ylim(ylims[iimeas])
  axarr[iimeas].set_yticks(yticks[imeas])
  axarr[iimeas].set_yticklabels([str(x) for x in yticks[imeas]],fontsize=5)
  for tick in axarr[iimeas].xaxis.get_major_ticks() + axarr[iimeas].yaxis.get_major_ticks():
    tick.label.set_fontsize(3.5)

  axarr[iimeas].text(3.48e6,ylims[iimeas][0]*0.01+ylims[iimeas][1]*0.99,captions[imeas],fontsize=6,va='top')

axarr[10].set_xlabel('time (min)',fontsize=5)
axarr[8].set_xlabel('time (min)',fontsize=5)
axarr[12].set_xlabel('time (min)',fontsize=5)
axarr[4].set_ylabel('relative conductance',fontsize=5)
axarr[8].set_ylabel('                                       relative conductance',fontsize=5)
axarr[12].set_ylabel('                                       relative conductance',fontsize=5)

axarr[13].set_visible(False)
axarr[14].set_visible(False)
axarr[15].set_visible(False)

f.savefig("fig9abc.eps")
