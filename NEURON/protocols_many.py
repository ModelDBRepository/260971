import numpy

protoparams_fixed = { 'Duration': 5260000, 'tolerance': 1e-6, 'Ca_input_onset': 4040000}
# Stimulus protocols:
protoparams_var = {
  'Ca_input_Ns':      [100,      156,      4,        4,        50,       4,        180,      5],
  'Ca_input_freqs':   [100,      312,      100,      100,      0.1,      100,      1,        100],
  'Ca_input_Ntrains': [1,        1,        10,       50,       1,        10,       1,        25],
  'Ca_input_trainTs': [1,        1,        200,      200,      1,        200,      1,        1000],
  'Ca_input_durs':    [3,        3,        3,        3,        3,        3,        15,       3]    #used a standard 3-ms Ca flux for all cases. 900 x 5Hz x 3ms replaced by 180 x 1Hz x 15ms for speed
}

#Measure 0) Ca
#        1) S845-phos GluR1
#        2) S831-phos GluR1
#        3) double-phos GluR1
#        4) membrane-inserted GluR1
#        5) membrane-inserted S831-phos GluR1
#        6) S880-phos GluR2
#        7) membrane-inserted GluR2
#        8) synaptic maximal conductance, which is calculated from 4), 5), and 7)
Measured_species = [ ['Ca'],
                     ['GluR1_S845', 'GluR1_S845_S831', 'GluR1_S845_CKCam', 'GluR1_S845_CKpCam', 'GluR1_S845_CKp', 'GluR1_S845_PKCt', 'GluR1_S845_PKCp', 'GluR1_S845_PP1', 'GluR1_S845_S831_PP1', 'GluR1_S845_PP2B', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S845', 'GluR1_memb_S845_S831', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCt', 'GluR1_memb_S845_PKCp', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'],
                     ['GluR1_S831', 'GluR1_S845_S831', 'GluR1_S831_PKAc', 'GluR1_S845_S831_PP1', 'GluR1_S831_PP1', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'],
                     ['GluR1_S845_S831', 'GluR1_S845_S831_PP1', 'GluR1_S845_S831_PP2B', 'GluR1_memb_S845_S831', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S845_S831_PP2B'],
                     ['GluR1_memb', 'GluR1_memb_S845', 'GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_PKAc', 'GluR1_memb_CKCam', 'GluR1_memb_CKpCam', 'GluR1_memb_CKp', 'GluR1_memb_PKCt', 'GluR1_memb_PKCp', 'GluR1_memb_S845_CKCam', 'GluR1_memb_S845_CKpCam', 'GluR1_memb_S845_CKp', 'GluR1_memb_S845_PKCt', 'GluR1_memb_S845_PKCp', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_PP1', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_PP2B', 'GluR1_memb_S845_S831_PP2B'],
                     ['GluR1_memb_S831', 'GluR1_memb_S845_S831', 'GluR1_memb_S831_PKAc', 'GluR1_memb_S845_S831_PP1', 'GluR1_memb_S831_PP1', 'GluR1_memb_S845_S831_PP2B'],
                     ['GluR2_S880', 'GluR2_S880_PP2A', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'],
                     ['GluR2_memb', 'GluR2_memb_PKCt', 'GluR2_memb_PKCp', 'GluR2_memb_S880', 'GluR2_memb_S880_PP2A'],
                     'syncond' ]

Quantification_types = ['abs(target-max val)', 'abs(target-last val)', 'abs(target-(last val/baseline))', 'abs(target-(val(t)/baseline))']

# Experiments: [ [STIMULUS PROTOCOL INDEX], [CAFLUX COEFF], [LFLUX COEFF], [GLUFLUX COEFF], [ACHFLUX COEFF], [BLOCKED], [ALTERED] ]
Experiments = [ [0, 1.0,  0.0, 1.0, 0.0, 'None', []],                  #0; Ma 2008 differential
                [0, 1.0,  0.0, 1.0, 0.0, 'None', [[125,126,127],0.0]], #1; (CK phosphorylation blocked)
                [0, 0.01, 0.0, 1.0, 0.0, 'None', []],                  #2; (post-syn Ca blocked)
                [0, 1.0,  0.0, 1.0, 0.0, 'None', [[317],0]],           #3; (PKAc separation from PKAcAMP4 blocked)
                [1, 1.0,  1.0, 1.0, 0.0, 'None', []],                  #4; Saez-Briones 2015 b2-Adrenoceptor and Flores 2011 hidden
                [1, 1.0,  0.0, 1.0, 0.0, 'None', []],                  #5;
                [2, 1.0,  0.0, 1.0, 0.0, 'None', []],                  #6; Hardingham 2003 neocortical 
                [2, 1.0,  0.0, 1.0, 0.0, 'None', [[125,126,127],0.0]], #7; (CK phosphorylation blocked)
                [3, 1.0,  0.0, 1.0, 0.0, 'None', []],                  #8; Song 2017 selective 
                [3, 1.0,  0.0, 1.0, 0.0, 'None', [[154,187,206,239],0.0]], #9 (s845 phosphorylation by PKA blocked)
                [3, 1.0,  0.0, 1.0, 0.0, 'None', [[157,160,163,166,169,172,175,178,181,184,209,212,215,218,221,224,227,230,233,236],0.0]], #10 (s831 phosphorylation by PKC and CK blocked)
                [4, 1.0,  1.0, 1.0, 0.0, 'None', []],                  #11; Zhou 2013 activation
                [4, 1.0,  0.0, 1.0, 0.0, 'None', []],                  #12;
                [5, 1.0,  0.0, 1.0, 0.0, 'None', []],                  #13; Kirkwood 1997 age-dependent
                [5, 1.0,  0.0, 1.0, 0.0, 'CKx0.0', []],                #14; (CK knockout)
                [6, 1.0,  0.0, 1.0, 0.0, 'None', []],                  #15;
                [6, 1.0,  0.0, 1.0, 0.0, 'CKx0.0', []],                #16; (CK knockout)
                [7, 1.0,  0.0, 1.0, 0.0, 'None', []] ]                 #17; Kotak 2007 developmental 

#Measurement: [ [EXPERIMENT_INDEX], [TARGET_T], [TARGET_VAL] ]
Measurements = [ [ [0,1,2],       [4640000, 4940000, 5240000], [[1.3,1.4,1.3],[1.05,1.02,0.95],[1.05,1.05,1.1]] ],                              #Ma 2008 differential, horizontal
                 [ [0,3,2],       [4640000, 4940000, 5240000], [[1.6,1.6,1.6],[1.4,1.4,1.4],[1.3,1.4,1.4]] ],                                   #Ma 2008 differential, ascending
                 [ [4,5],         [4640000, 4940000, 5240000], [[2.0,1.98,1.9],[1.34,1.4,1.36]] ],                                              #Saez-Briones 2015 b2-Adrenoceptor
                 [ [4,5],         [4640000, 4940000, 5240000], [[1.7,1.6,1.64],[1.43,1.45,1.43]] ],                                             #Flores 2011 hidden
                 [ [6,7],         [4640000, 4940000, 5240000], [[1.35,1.4,1.3],[1.25,1.2,1.1]] ],                                               #Hardingham 2003 neocortical
                 [ [8,9,10],      [4640000, 4940000, 5240000], [[1.55,1.4,1.4],[1.1,1.05,1.05],[1.35,1.4,1.3]] ],                               #Song 2017 selective
                 [ [11,12],       [4640000, 4940000, 5240000], [[1.3,1.4,1.4],[1.1,1.2,1.2]] ],                                                 #Zhou 2013 activation
                 [ [13,14,15,16], [4640000, 4940000, 5240000], [[1.3,1.26,1.26],[1.02,1.02,1.02],[numpy.nan,0.95,0.95],[numpy.nan,0.88,0.93]] ],#Kirkwood 1997 age-dependent, adult neurons
                 [ [13,14,15,16], [4640000, 4940000, 5240000], [[1.2,1.18,1.18],[1.07,1.09,1.08],[numpy.nan,0.79,0.82],[numpy.nan,0.82,0.89]] ],#Kirkwood 1997 age-dependent, 4-5 week old neurons
                 [ [17],          [4640000, 4940000, 5240000], [[1.98,1.58,1.93]] ],                                                            #Kotak 2007 developmental, LTP-expressing neurons
                 [ [17],          [4640000, 4940000, 5240000], [[0.77,0.68,0.67]] ] ]                                                           #Kotak 2007 developmental, LTD-expressing neurons
Measurements_txt = [['Ma 2008 differential, horizontal', 'CONTROL', 'CK BLOCKED',  'POST-SYN CA BLOCKED'],
                    ['Ma 2008 differential, ascending', 'CONTROL', 'PKA BLOCKED',  'POST-SYN CA BLOCKED'],
                    ['Saez-Briones 2015 b2-Adrenoceptor', 'WITH L', 'WITHOUT L'],
                    ['Flores 2011 hidden', 'WITH L', 'WITHOUT L'],
                    ['Hardingham 2003 neocortical', 'CONTROL', 'CK BLOCKED'],
                    ['Song 2017 selective', 'CONTROL', 'S845 BLOCKED', 'S831 BLOCKED'],
                    ['Zhou 2013 activation', 'WITH L', 'WITHOUT L'],
                    ['Kirkwood 1997 age-dependent, adult', 'TBS, CONTROL', 'TBS, CK KNOCKOUT', 'LFS, CONTROL', 'LFS, CK KNOCKOUT'],
                    ['Kirkwood 1997 age-dependent, 4-5 week', 'TBS, CONTROL', 'TBS, CK KNOCKOUT', 'LFS, CONTROL', 'LFS, CK KNOCKOUT'],
                    ['Kotak 2007 developmental', 'LTP'],
                    ['Kotak 2007 developmental', 'LTD'] ]
Measurements_stds = [[[0.1, 0.1, 0.1], [0.07, 0.08, 0.05], [0.08, 0.1, 0.08]],
                     [[0.12, 0.1, 0.12], [0.11, 0.13, 0.15], [0.15, 0.12, 0.13]],
                     [[0.08, 0.09, 0.08], [0.08, 0.09, 0.1]],
                     [[0.13, 0.13, 0.1], [0.11, 0.1, 0.09]],
                     [[0.07, 0.1, 0.09], [0.09, 0.12, 0.07]],
                     [[0.05, 0.05, 0.05], [0.07, 0.07, 0.07], [0.1, 0.1, 0.1]], #Figs. 1J, 1E, 4E
                     [[0.15, 0.17, 0.1], [0.1, 0.1, 0.18]],                     #Figs. 2D and 2C
                     [[0.08, 0.07, 0.07], [0.015, 0.02, 0.015], [numpy.nan, 0.05, 0.04], [numpy.nan, 0.03, 0.03]], #Figs. 1A, 1B, 3A, 3C, white
                     [[0.05, 0.05, 0.05], [0.02, 0.03, 0.03], [numpy.nan, 0.03, 0.02], [numpy.nan, 0.035, 0.03]],  #Figs. 1A, 1B, 3A, 3C, black
                     [[0.25, 0.11, 0.21]],                                                                         #Figs. 4A, 4B
                     [[0.09, 0.1, 0.08]]]                                                                          #Figs. 4A, 4B

def get_measurement_protocol():
  return [Measurements, Experiments, protoparams_fixed, protoparams_var, Measured_species, Quantification_types, Measurements_txt, Measurements_stds]
