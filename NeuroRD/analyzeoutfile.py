from pylab import *
import scipy.io
exceptions = [['R','GluR']]


def getlastvalues(filename):
  if filename[-4:] == '.mat':
    DATA = getlastvaluesmat(filename)
    return DATA
  input_file = open(filename+'.out','r')
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

  mesh_input_file = open(filename+'-mesh.txt.out','r')
  mesh_firstline = mesh_input_file.readline()
  mesh_secondline = mesh_input_file.readline()
  mesh_values = mesh_secondline.split()
  my_volume = float(mesh_values[-2])*1e-15 #litres                                                                                                                                                                                           
  mesh_input_file.close()

  DATA = array(values_all)[:,0,:]/6.022e23/my_volume*1e9 #DATA in the first (and only) compartment
  return [DATA[-1,:], header_strs]

def getlastvaluesmat(filename):
  MAT = scipy.io.loadmat(filename)
  DATA_all = MAT['DATA']
  header_strs = MAT['headers']
  mesh_input_file = open('mesh_general.out','r')
  mesh_firstline = mesh_input_file.readline()
  mesh_secondline = mesh_input_file.readline()
  mesh_values = mesh_secondline.split()
  my_volume = float(mesh_values[-2])*1e-15 #litres
  mesh_input_file.close()
  headers = [s + ' ' for s in header_strs]
  return [DATA_all[-1,:]/6.022e23/my_volume*1e9, [s[:s.find(' ')] for s in headers]]

def getinitvalues():
  #return [73, 1877026, 216332, 161, 145869, 14658, 0, 0, 2502689, 490, 13, 16767, 537911, 5711, 4060, 9, 54, 3820, 2623, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1508, 0, 0, 0, 0, 0, 0, 0, 4166, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1994911, 55, 5, 81, 3709, 2, 0, 11482, 297, 0, 976, 2143, 18013, 11485, 180, 4, 0, 2060, 266, 3, 21291, 7, 1059, 868, 0, 0, 7, 6, 1826, 238, 120, 0, 0, 729, 0, 6, 1121, 442, 0, 0, 11, 270, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 505, 0, 0, 260, 0, 0, 1.0e6, 0, 0, 2472, 0, 0, 0, 1019075, 1441, 0, 0, 493, 0, 0, 0, 23783, 0, 0, 0, 0, 0, 293, 402, 14994, 0, 0, 0, 695, 0, 0, 1556, 253, 0, 0, 0, 601]
  #cat IC_singlecompartment.xml | grep specie | cut -d'"' -f2 |while read F; do printf "'$F', ";done
  species = ['Ca', 'CaOut', 'CaOutLeak', 'Leak', 'Calbin', 'CalbinC', 'LOut', 'Epac1', 'Epac1cAMP', 'PMCA', 'NCX', 'PMCACa', 'NCXCa', 'L', 'R', 'Gs', 'Gi', 'LR', 'LRGs', 'PKAcLR', 'PKAcpLR', 'PKAcppLR', 'PKAcpppLR', 'pLR', 'ppLR', 'pppLR', 'ppppLR', 'ppppLRGi', 'ppppLRGibg', 'PKAcR', 'PKAcpR', 'PKAcppR', 'PKAcpppR', 'pR', 'ppR', 'pppR', 'ppppR', 'ppppRGi', 'ppppRGibg', 'GsR', 'GsaGTP', 'GsaGDP', 'GiaGTP', 'GiaGDP', 'Gibg', 'Gsbg', 'LRGsbg', 'AC1', 'AC1GsaGTP', 'AC1GsaGTPCaMCa4', 'AC1GsaGTPCaMCa4ATP', 'AC1GiaGTP', 'AC1GiaGTPCaMCa4', 'AC1GiaGTPCaMCa4ATP', 'AC1GsaGTPGiaGTP', 'AC1GsaGTPGiaGTPCaMCa4', 'AC1GsGiCaMCa4ATP', 'ATP', 'cAMP', 'AC1CaMCa4', 'AC1CaMCa4ATP', 'AC8', 'AC8CaMCa4', 'AC8CaMCa4ATP', 'PDE1', 'PDE1CaMCa4', 'PDE1CaMCa4cAMP', 'AMP', 'Ng', 'NgCaM', 'CaM', 'CaMCa2', 'CaMCa3', 'CaMCa4', 'PP2B', 'PP2BCaM', 'PP2BCaMCa2', 'PP2BCaMCa3', 'PP2BCaMCa4', 'CK', 'CKCaMCa4', 'CKpCaMCa4', 'CKp', 'Complex', 'pComplex', 'CKpPP1', 'CKpCaMCa4PP1', 'PKA', 'PKAcAMP4', 'PKAr', 'PKAc', 'I1', 'I1PKAc', 'Ip35', 'PP1', 'Ip35PP1', 'Ip35PP2BCaMCa4', 'Ip35PP1PP2BCaMCa4', 'PP1PP2BCaMCa4', 'GluR1', 'GluR1_S845', 'GluR1_S831', 'GluR1_S845_S831', 'GluR1_PKAc', 'GluR1_CKCam', 'GluR1_CKpCam', 'GluR1_CKp', 'GluR1_PKCt', 'GluR1_PKCp', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_S845_CKp', 'GluR1_S845_PKCt', 'GluR1_S845_PKCp', 'GluR1_S831_PKAc', 'GluR1_S845_PP1', 'GluR1_S845_S831_PP1', 'GluR1_S831_PP1', 'GluR1_S845_PP2B', 'GluR1_S845_S831_PP2B', 'GluR1_memb', 'GluR1_memb_S845', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_PKAc', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_CKp', 'GluR1_memb_PKCt', 'GluR1_memb_PKCp', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCt', 'GluR1_memb_S845_PKCp', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B', 'PDE4', 'PDE4cAMP', 'PKAcPDE4', 'pPDE4', 'pPDE4cAMP', 'PKAc_PDE4_cAMP', 'fixedbuffer', 'fixedbufferCa', 'Glu', 'MGluR', 'MGluR_Glu', 'MGluR_Glu_desens', 'MGluR_Gqabg_Glu', 'GluOut', 'Gqabg', 'GqaGTP', 'GqaGDP', 'PLC', 'PLCCa', 'PLCCaGqaGTP', 'PLCGqaGTP', 'Pip2', 'PLCCaPip2', 'PLCCaGqaGTPPip2', 'Ip3', 'PLCCaDAG', 'PLCCaGqaGTPDAG', 'PIkinase', 'Ip3degPIk', 'PKC', 'PKCCa', 'PKCt', 'PKCp', 'DAG', 'DAGK', 'DAGKdag', 'PA', 'DGL', 'CaDGL', 'DAGCaDGL', '2AG', '2AGdegrad', 'Ip3degrad', 'GluR2', 'GluR2_PKCt', 'GluR2_PKCp', 'GluR2_S880', 'GluR2_S880_PP2A', 'GluR2_memb', 'GluR2_memb_PKCt', 'GluR2_memb_PKCp', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A', 'PP2A', 'ACh', 'M1R', 'AChM1R', 'M1RGq', 'AChM1RGq', 'PLA2', 'CaPLA2', 'CaPLA2Pip2', 'AA']
  #species = ['Ca','CaOut','CaOutLeak','Leak','Calbin','CalbinC','CaBCa','CaB','LOut','Epac1','Epac1cAMP','pmca','ncx','pmcaCa','ncxCa','L','R','Gs','Gi','LR','LRGs','PKAcLR','PKAcpLR','PKAcppLR','PKAcpppLR','pLR',
  #'ppLR','pppLR','ppppLR','ppppLRGi','ppppLRGibg','PKAcR','PKAcpR','PKAcppR','PKAcpppR','pR','ppR','pppR','ppppR','ppppRGi','ppppRGibg','GsR','GasGTP','GasGDP','GaiGTP','GaiGDP','Gibg','Gsbg','LRGsbg','AC1','AC1GasGTP',
  #'AC1GasGTPCaMCa4','AC1GasGTPCaMCa4ATP','AC1GaiGTP','AC1GaiGTPCaMCa4','AC1GaiGTPCaMCa4ATP','AC1GasGTPGaiGTP','AC1GasGTPGaiGTPCaMCa4','AC1GsGiCaMCa4ATP','ATP','cAMP','AC1CaMCa4','AC1CaMCa4ATP','AC8','AC8CaMCa4','AC8CaMCa4ATP',
  #'PDE1','PDE1CaMCa4','PDE1CaMCa4cAMP','AMP','Ng','NgCaM','CaM','CaMCa2','CaMCa4','PP2B','PP2BCaM','PP2BCaMCa2','PP2BCaMCa4','CK','CKCaMCa4','CKpCaMCa4','CKp','Complex','pComplex','CKpPP1','CKpCaMCa4PP1','PKA','PKAcAMP2',
  #'PKAcAMP4','PKAr','PKAc','I1','I1PKAc','Ip35','PP1','Ip35PP1','Ip35PP2BCaMCa4','Ip35PP1PP2BCaMCa4','PP1PP2BCaMCa4','GluR1','GluR1_S845','GluR1_S831','GluR1_S845_S831','GluR1_PKAc','GluR1_CKCam','GluR1_CKpCam','GluR1_CKp',
  #'GluR1_S845_CKCam','GluR1_S845_CKpCam','GluR1_S845_CKp','GluR1_S831_PKAc','GluR1_S845_PP1','GluR1_S845_S831_PP1','GluR1_S831_PP1','GluR1_S845_PP2B','GluR1_S845_S831_PP2B','GluR1_memb','GluR1_memb_S845','GluR1_memb_S831',
  #'GluR1_memb_S845_S831','GluR1_memb_PKAc','GluR1_memb_CKCam','GluR1_memb_CKpCam','GluR1_memb_CKp','GluR1_memb_S845_CKCam','GluR1_memb_S845_CKpCam','GluR1_memb_S845_CKp','GluR1_memb_S831_PKAc','GluR1_memb_S845_PP1',
  #'GluR1_memb_S845_S831_PP1','GluR1_memb_S831_PP1','GluR1_memb_S845_PP2B','GluR1_memb_S845_S831_PP2B','PDE4','PDE4cAMP','PKAcPDE4','pPDE4','pPDE4cAMP','PKAc_PDE4_cAMP','fixedbuffer','fixedbufferCa','Glu','MGluR','MGluR_Glu',
  #'MGluR_Glu_desens','MGluR_Gabgq_Glu','GluOut','Gabgq','GaqGTP','GaqGDP','Plc','PlcCa','PlcCaGaqGTP','PlcGaqGTP','Pip2','PlcCaPip2','PlcCaGaqGTPPip2','Ip3','PlcCaDag','PlcCaGaqGTPDag','PIkinase','Ip3degPIk','Pkc','PkcCa',
  #'Pkct','Dag','DagK','DagKdag','PA','Dgl','CaDgl','DagCaDgl','2ag','2agDegrad','Ip3degrad','GluR2','GluR2_Pkct','GluR2_S880','GluR2_S880_PP2A','GluR2_memb','GluR2_memb_Pkct','GluR2_memb_S880','GluR2_memb_S880_PP2A','PP2A',
  #'ACh','m1R','AChm1R','m1RGq','AChm1RGq']

  input_file = open('IC_singlecompartment.xml','r')
  IC_dict = {}
  line = input_file.readline()
  while len(line) > 0:
    if line.find('NanoMolarity specieID=') > -1 and (line.find('<!--') == -1 or line.find('<!--') > line.find('NanoMolarity specieID=')):
      split_by_quotes = line.split('"')
      IC_dict[split_by_quotes[1]] = int(float(split_by_quotes[3]))
    line = input_file.readline()
  input_file.close()
  initvalues = []
  for ispec in range(0,len(species)):
    initvalues.append(IC_dict[species[ispec]])
  return initvalues


def writeICfile(initvalues,header_strs,ICfilename='IC_this.xml',blocks=[]):
  #cat IC_singlecompartment.xml | grep specie | cut -d'"' -f2 |while read F; do printf "'$F', ";done
  species = ['Ca', 'CaOut', 'CaOutLeak', 'Leak', 'Calbin', 'CalbinC', 'LOut', 'Epac1', 'Epac1cAMP', 'PMCA', 'NCX', 'PMCACa', 'NCXCa', 'L', 'R', 'Gs', 'Gi', 'LR', 'LRGs', 'PKAcLR', 'PKAcpLR', 'PKAcppLR', 'PKAcpppLR', 'pLR', 'ppLR', 'pppLR', 'ppppLR', 'ppppLRGi', 'ppppLRGibg', 'PKAcR', 'PKAcpR', 'PKAcppR', 'PKAcpppR', 'pR', 'ppR', 'pppR', 'ppppR', 'ppppRGi', 'ppppRGibg', 'GsR', 'GsaGTP', 'GsaGDP', 'GiaGTP', 'GiaGDP', 'Gibg', 'Gsbg', 'LRGsbg', 'AC1', 'AC1GsaGTP', 'AC1GsaGTPCaMCa4', 'AC1GsaGTPCaMCa4ATP', 'AC1GiaGTP', 'AC1GiaGTPCaMCa4', 'AC1GiaGTPCaMCa4ATP', 'AC1GsaGTPGiaGTP', 'AC1GsaGTPGiaGTPCaMCa4', 'AC1GsGiCaMCa4ATP', 'ATP', 'cAMP', 'AC1CaMCa4', 'AC1CaMCa4ATP', 'AC8', 'AC8CaMCa4', 'AC8CaMCa4ATP', 'PDE1', 'PDE1CaMCa4', 'PDE1CaMCa4cAMP', 'AMP', 'Ng', 'NgCaM', 'CaM', 'CaMCa2', 'CaMCa3', 'CaMCa4', 'PP2B', 'PP2BCaM', 'PP2BCaMCa2', 'PP2BCaMCa3', 'PP2BCaMCa4', 'CK', 'CKCaMCa4', 'CKpCaMCa4', 'CKp', 'Complex', 'pComplex', 'CKpPP1', 'CKpCaMCa4PP1', 'PKA', 'PKAcAMP4', 'PKAr', 'PKAc', 'I1', 'I1PKAc', 'Ip35', 'PP1', 'Ip35PP1', 'Ip35PP2BCaMCa4', 'Ip35PP1PP2BCaMCa4', 'PP1PP2BCaMCa4', 'GluR1', 'GluR1_S845', 'GluR1_S831', 'GluR1_S845_S831', 'GluR1_PKAc', 'GluR1_CKCam', 'GluR1_CKpCam', 'GluR1_CKp', 'GluR1_PKCt', 'GluR1_PKCp', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_S845_CKp', 'GluR1_S845_PKCt', 'GluR1_S845_PKCp', 'GluR1_S831_PKAc', 'GluR1_S845_PP1', 'GluR1_S845_S831_PP1', 'GluR1_S831_PP1', 'GluR1_S845_PP2B', 'GluR1_S845_S831_PP2B', 'GluR1_memb', 'GluR1_memb_S845', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_PKAc', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_CKp', 'GluR1_memb_PKCt', 'GluR1_memb_PKCp', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCt', 'GluR1_memb_S845_PKCp', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B', 'PDE4', 'PDE4cAMP', 'PKAcPDE4', 'pPDE4', 'pPDE4cAMP', 'PKAc_PDE4_cAMP', 'fixedbuffer', 'fixedbufferCa', 'Glu', 'MGluR', 'MGluR_Glu', 'MGluR_Glu_desens', 'MGluR_Gqabg_Glu', 'GluOut', 'Gqabg', 'GqaGTP', 'GqaGDP', 'PLC', 'PLCCa', 'PLCCaGqaGTP', 'PLCGqaGTP', 'Pip2', 'PLCCaPip2', 'PLCCaGqaGTPPip2', 'Ip3', 'PLCCaDAG', 'PLCCaGqaGTPDAG', 'PIkinase', 'Ip3degPIk', 'PKC', 'PKCCa', 'PKCt', 'PKCp', 'DAG', 'DAGK', 'DAGKdag', 'PA', 'DGL', 'CaDGL', 'DAGCaDGL', '2AG', '2AGdegrad', 'Ip3degrad', 'GluR2', 'GluR2_PKCt', 'GluR2_PKCp', 'GluR2_S880', 'GluR2_S880_PP2A', 'GluR2_memb', 'GluR2_memb_PKCt', 'GluR2_memb_PKCp', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A', 'PP2A', 'ACh', 'M1R', 'AChM1R', 'M1RGq', 'AChM1RGq', 'PLA2', 'CaPLA2', 'CaPLA2Pip2', 'AA']
  output_file = open(ICfilename,'w')
  output_file.write('<InitialConditions>\n')
  output_file.write('  <ConcentrationSet>\n')
  print "blocks = "+str(blocks)
  for ispecies in range(0,len(species)):
    found = 0
    for ispecies2 in range(0,len(header_strs)):
      if header_strs[ispecies2] == species[ispecies]:
        isblocked = -1
        for iblock in range(0,len(blocks)):
          if header_strs[ispecies2].find(blocks[iblock][0]) > -1:
            isException = 0
            for iexc in range(0,len(exceptions)):
              if blocks[iblock][0] == exceptions[iexc][0] and header_strs[ispecies2].find(exceptions[iexc][1]) > -1:
                isException = 1
            if isException == 0:
              isblocked = iblock
            elif blocks[iblock][0] == 'GluR' and header_strs[ispecies2].find(blocks[iblock][0]) == 0:
              print 'blocks[iblock][0] == ''GluR'' and header_strs[ispecies2].find(blocks[iblock][0]) == '+str(header_strs[ispecies2].find(blocks[iblock][0]))+' == 0'
              isblocked = iblock
        if isblocked < 0:
          output_file.write('    <NanoMolarity specieID="'+str(species[ispecies])+'" value="'+str(initvalues[ispecies2-4])+'"/>\n')
        else:
          output_file.write('    <NanoMolarity specieID="'+str(species[ispecies])+'" value="'+str(initvalues[ispecies2-4]*blocks[isblocked][1])+'"/>\n')
        found = 1
        break
    if not found:
      print 'Error: '+species[ispecies]+' not found'
    else:
      print "ispecies="+str(ispecies)+", ispecies2="+str(ispecies2)+" found"
  output_file.write('  </ConcentrationSet>\n')
  output_file.write('</InitialConditions>\n')
  output_file.close()






