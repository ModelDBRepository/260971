#printtab1: This file prints the Table 1 based on Reactions.xml
#Tuomo Maki-Marttunen, 2019

from pylab import *
import scipy.io
import sys
import time

class bcolors:
  OKBLUE = '\033[94m'
  OKGREEN = '\033[92m'
  FAIL = '\033[91m'
  ENDC = '\033[0m'
  BOLD = '\033[1m'
  UNDERLINE = '\033[4m'

Reacs = []
Reactants = []
Reactant_powers = []
Products = []
Product_powers = []
Forwardrates = []
Backwardrates = []
reac_file = open('Reactions.xml')
firstline = reac_file.readline()
line = reac_file.readline()
nsAll = []
nsProdAll = []

print_reac_N = True
if len(sys.argv) > 1 and int(sys.argv[1]) == 0:
  print_reac_N = False

### Load reactions ###
toPrint = []
reacRatesToPrint = []
while len(line) > 0:
  if line.find('Reaction name') > -1 and line[0:4] != '<!--':
    newline = reac_file.readline()
    #print "newline = "+newline

    Reactants_this = []
    Products_this = []
    Forwardrate = 0
    Backwardrate = 0
    Powers = []
    PowersProd = []
    ns = []
    nsProd = []
    while newline.find('</Reaction>') == -1:
      if newline.find('Reactant specieID') > -1:
        Reactants_this.append(newline[newline.find('"')+1:newline.find('"')+newline[newline.find('"')+1:].find('"')+1])
        if newline.find('power') > -1:
          newline2 = newline[newline.find('power'):]
          Powers.append(int(newline2[newline2.find('"')+1:newline2.find('"')+newline2[newline2.find('"')+1:].find('"')+1]))
        else:
          Powers.append(1)
        if newline.find('n=') > -1:
          newline2 = newline[newline.find('n='):]
          ns.append(int(newline2[newline2.find('"')+1:newline2.find('"')+newline2[newline2.find('"')+1:].find('"')+1]))
        else:
          ns.append(Powers[-1])
      if newline.find('Product specieID') > -1 or newline.find('Product  specieID') > -1:
        Products_this.append(newline[newline.find('"')+1:newline.find('"')+newline[newline.find('"')+1:].find('"')+1])
        if newline.find('power') > -1:
          newline2 = newline[newline.find('power'):]
          PowersProd.append(int(newline2[newline2.find('"')+1:newline2.find('"')+newline2[newline2.find('"')+1:].find('"')+1]))
        else:
          PowersProd.append(1)
        if newline.find('n=') > -1:
          newline2 = newline[newline.find('n='):]
          nsProd.append(int(newline2[newline2.find('"')+1:newline2.find('"')+newline2[newline2.find('"')+1:].find('"')+1]))
        else:
          nsProd.append(PowersProd[-1])
      if newline.find('forwardRate') > -1:
        newline2 = newline[newline.find('>')+1:]
        Forwardrate = float(newline2[:newline2.find('<')])
      if newline.find('reverseRate') > -1:
        newline2 = newline[newline.find('>')+1:]
        Backwardrate = float(newline2[:newline2.find('<')])
      newline = reac_file.readline()
    Forwardrates.append(Forwardrate)
    Backwardrates.append(Backwardrate)
    Reactants.append(Reactants_this[:])
    Products.append(Products_this[:])
    Reactant_powers.append(Powers[:])
    Product_powers.append(PowersProd[:])
    nsAll.append(ns[:])
    nsProdAll.append(nsProd[:])
    Reac_txt = ''
    for i in range(0,len(Reactants_this)):
      if ns[i] != 1:
        Reac_txt = Reac_txt + str(ns[i])+"*"+Reactants_this[i]
      else:
        Reac_txt = Reac_txt + Reactants_this[i]
      if i < len(Reactants_this) - 1:
        Reac_txt = Reac_txt + " + "
    Reac_txt = Reac_txt + " $\\rightleftharpoons$ "
    for i in range(0,len(Products_this)):
      if nsProd[i] != 1:
        Reac_txt = Reac_txt + str(nsProd[i])+"*"+Products_this[i]
      else:
        Reac_txt = Reac_txt + Products_this[i]
      if i < len(Products_this) - 1:
        Reac_txt = Reac_txt + " + "
    addition = ""
    reacRatesToPrint.append([Forwardrate,Backwardrate])
    if ns == Powers and nsProd == PowersProd:
      toPrint.append(addition+Reac_txt+" & "+str(Forwardrate)+" & "+str(Backwardrate) + " & \\\\")
    elif ns != Powers:
      toPrint.append(addition+Reac_txt+" & "+str(Forwardrate)+" & "+str(Backwardrate) + " & $*$ \\\\")
    elif nsProd != PowersProd:
      toPrint.append(addition+Reac_txt+" & "+str(Forwardrate)+" & "+str(Backwardrate) + " & $\dagger$ \\\\")
    Reacs.append(Reac_txt)
  line = reac_file.readline()    
reac_file.close()

### Divide reactions to groups ###
toPrintHere = [[i] for i in range(0,len(toPrint))]
extraConditions = 2
for ireac in range(0,len(toPrint)):
  for ireacgone in range(0,ireac):
    if reacRatesToPrint[ireac][0] == reacRatesToPrint[ireacgone][0] and reacRatesToPrint[ireac][1] == reacRatesToPrint[ireacgone][1]:
      if extraConditions > 0:
        #Extra conditions:
        #1) There has to be same number of reactants and products
        if len(Reactants[ireac]) != len(Reactants[ireacgone]) or len(Products[ireac]) != len(Products[ireacgone]): 
          continue
        
      if extraConditions > 1:
        #2) There has to be a sequence of at least three letters that appears in the two corresponding reactants/products
        foundOne = False
        for wordpair in zip(Reactants[ireac]+Products[ireac],Reactants[ireacgone]+Products[ireacgone])+zip(Reactants[ireacgone]+Products[ireacgone],Reactants[ireac]+Products[ireac]):
          if (wordpair[0] == 'R' and 'R' in wordpair[1]) or (wordpair[0] == 'pR' and 'pR' in wordpair[1]) or (wordpair[0] == 'ppR' and 'ppR' in wordpair[1]):
            foundOne = True
          for itriple in range(0,len(wordpair[0])-3):
            if wordpair[0][itriple:itriple+3] in wordpair[1]:
              foundOne = True
              #print wordpair[0]+" has a triplet included in "+wordpair[1]
              break
        #Exceptions: don't let PKAc, if being the second product, be replaced by Gibg
        if len(Products[ireac]) == 2 and len(Products[ireacgone]) == 2 and (Products[ireac][1] == 'PKAc' and Products[ireacgone][1] == 'Gibg' or Products[ireac][1] == 'Gibg' and Products[ireacgone][1] == 'PKAc'):
          foundOne = False
        if not foundOne:
          continue
      if extraConditions > 2:
        #3) One of the reactants or products has to be included in the name of the corresponding reactant or product in the other reaction
        reacNamesOK = 0 
        for ireactant in range(0,len(Reactants[ireac])):
          if Reactants[ireac][ireactant] in Reactants[ireacgone][ireactant] or Reactants[ireacgone][ireactant] in Reactants[ireac][ireactant]:
            reacNamesOK = 1
        for iproduct in range(0,len(Products[ireac])):
          if Products[ireac][iproduct] in Products[ireacgone][iproduct] or Products[ireacgone][iproduct] in Products[ireac][iproduct]:
            reacNamesOK = 1
        if not reacNamesOK:
          continue

      if extraConditions > 3:
        #4) The length of the first (or nth) reactant (or product) has to be a fixed number of letters shorter or longer than the name of the first product (or reactant)
        lengthsOK = 0
        for iproduct in range(0,len(Products[ireac])):
          if len(Reactants[ireac][0]) - len(Reactants[ireacgone][0]) == len(Products[ireac][iproduct]) - len(Products[ireacgone][iproduct]):
            lengthsOK = 1
        for ireactant in range(0,len(Reactants[ireac])):
          if len(Products[ireac][0]) - len(Products[ireacgone][0]) == len(Reactants[ireac][ireactant]) - len(Reactants[ireacgone][ireactant]):
            lengthsOK = 1
        if not lengthsOK:
          continue

      #Only add the reaction if the reaction is printed where it stands
      if ireacgone in toPrintHere[ireacgone]:
        toPrintHere[ireacgone].append(ireac)
        toPrintHere[ireac] = []
      else: #otherwise, find where it is printed and print it there
        foundOne = False
        for ireacgone2 in range(0,ireacgone):
          if ireacgone in toPrintHere[ireacgone2]:
            foundOne = True
            break
        if foundOne:
          toPrintHere[ireacgone2].append(ireac)
          toPrintHere[ireac] = []          
      break


### Order the reactions by the groups they belong to ###
toPrintFinal = []
originalLocation = []
for ireac in range(0,len(toPrint)):
  for ireachere in range(0,len(toPrintHere[ireac])):
    toPrintFinal.append(toPrint[toPrintHere[ireac][ireachere]].replace('_','\_'))
    originalLocation.append(toPrintHere[ireac][ireachere])
originalLocationInv = []
for i in range(0,len(originalLocation)):
  originalLocationInv.append(find(array(originalLocation)==i)[0])

# Optional: print all reactions by the groups
#for ireac in range(0,len(toPrint)):
#  print str(ireac) + " (" + str(originalLocation[ireac])+") "+toPrintFinal[ireac]

print "\n"
print "\n"
print "\n"
diffCounter = 1
variables = []
variableNames = []
ReacsGrouped = []
for ireac in range(0,len(toPrintHere)):
  if len(toPrintHere[ireac]) == 0:
    continue
  if len(toPrintHere[ireac]) == 1:
    ReacsGrouped.append(toPrintFinal[originalLocationInv[ireac]])
    continue
  print "    ireac="+str(ireac)+", toPrintHere[ireac] = "+str(toPrintHere[ireac])
  thisReac = ''
  NsameLetters_all = []
  for iword in range(0,len(Reactants[toPrintHere[ireac][0]]) + len(Products[toPrintHere[ireac][0]])):
    NsameLetters = 0
    if iword < len(Reactants[toPrintHere[ireac][0]]):
      target = Reactants
      iiword = iword
    else:
      target = Products
      iiword = iword - len(Reactants[toPrintHere[ireac][0]])
    for iletter in range(0,100):
      allSame = True
      for ireac2 in range(1,len(toPrintHere[ireac])):
        if iletter >= min(len(target[toPrintHere[ireac][0]][iiword]),len(target[toPrintHere[ireac][ireac2]][iiword])) or target[toPrintHere[ireac][0]][iiword][iletter] != target[toPrintHere[ireac][ireac2]][iiword][iletter]:
          allSame = False
          break
      if allSame:
        NsameLetters = NsameLetters + 1
    #Exceptions: If the common string is short and cuts the name of the species at an awkward place, then use 0 as NsameLetters. Other exceptions: make shorter NsameLetters
    if target[toPrintHere[ireac][0]][iiword][0:NsameLetters] in ['AC','C','G','P']:
      NsameLetters = 0
    elif target[toPrintHere[ireac][0]][iiword][NsameLetters-1] == '_':
      NsameLetters = NsameLetters-1
    if target[toPrintHere[ireac][0]][iiword][0:NsameLetters] == 'PP2BCaMCa':
      NsameLetters = 4
    if target[toPrintHere[ireac][0]][iiword][0:NsameLetters] == 'AC1GsaGT':
      NsameLetters = 6
    if 'GluR' in target[toPrintHere[ireac][0]][iiword][0:NsameLetters] and target[toPrintHere[ireac][0]][iiword][NsameLetters-2:NsameLetters] == '_S':
      NsameLetters = NsameLetters - 2
    NsameLetters_all.append(NsameLetters)

    #Check if the word is same in all reactions, i.e., if NsameLetters is the same as maximum length along the targets:
    if NsameLetters == max([len(target[toPrintHere[ireac][j]][iiword]) for j in range(0,len(toPrintHere[ireac]))]):
      thisReac = thisReac + target[toPrintHere[ireac][j]][iiword] + ' + '
      if iword == len(Reactants[toPrintHere[ireac][0]]) - 1:
        thisReac = thisReac[0:-2] + '$\\rightleftharpoons$ '
      continue

    #Check if the word is the same as one of the previous words:
    iwordsame = -1
    for iword2 in range(0,iword):
      allSame = True
      if iword2 < len(Reactants[toPrintHere[ireac][0]]):
        target2 = Reactants
        iiword2 = iword2
      else:
        target2 = Products
        iiword2 = iword2 - len(Reactants[toPrintHere[ireac][0]])
      NsameLetters2 = NsameLetters_all[iword2]
      for j in range(0,len(toPrintHere[ireac])):
        if target[toPrintHere[ireac][j]][iiword][NsameLetters:] != target2[toPrintHere[ireac][j]][iiword2][NsameLetters2:]:
          allSame = False
      if allSame:
        iwordsame = iword2
        break
      
    if iwordsame > -1:
      thisReac = thisReac + target[toPrintHere[ireac][0]][iiword][0:NsameLetters] + '$\mathbf{'+chr(ord('X')+iwordsame)+'}_{'+str(diffCounter)+'}$ + '      
      if iword == len(Reactants[toPrintHere[ireac][0]]) - 1:
        thisReac = thisReac[0:-2] + '$\\rightleftharpoons$ '
      continue

    thisReac = thisReac + target[toPrintHere[ireac][0]][iiword][0:NsameLetters] + '$\mathbf{'+chr(ord('X')+iword)+'}_{'+str(diffCounter)+'}$ + '      
    if iword == len(Reactants[toPrintHere[ireac][0]]) - 1:
      thisReac = thisReac[0:-2] + '$\\rightleftharpoons$ '
    variablesThis = []
    for ireac2 in range(0,len(toPrintHere[ireac])):
      if len(target[toPrintHere[ireac][ireac2]][iiword][NsameLetters:]) == 0:
        variablesThis.append('$\\{\\}$')
      else:
        variablesThis.append(target[toPrintHere[ireac][ireac2]][iiword][NsameLetters:])
    variables.append(variablesThis[:])
    variableNames.append('$\mathbf{'+chr(ord('X')+iword)+'}_{'+str(diffCounter)+'}$')
    targetNames = ''
    for i in range(0,len(toPrintHere[ireac])):
      targetNames = targetNames + target[toPrintHere[ireac][i]][iiword] + ', '
    #print "thisReac = "+str(thisReac)+", NsameLetters = "+str(NsameLetters)+" ("+targetNames[0:-2]+")"
    #time.sleep(0.5)

  ns = nsAll[toPrintHere[ireac][0]]
  nsProd = nsProdAll[toPrintHere[ireac][0]]
  Powers = Reactant_powers[toPrintHere[ireac][0]]
  PowersProd = Product_powers[toPrintHere[ireac][0]]
             
  reacStart = thisReac[0:-2].replace('_{',':::{').replace('_','\_').replace(':::{','_{')   #Change all '_'s to '\_'s except when followed by '{'
  if ns == Powers and nsProd == PowersProd:
    thisReac = reacStart + ' & '+str(Forwardrates[toPrintHere[ireac][0]])+' & '+str(Backwardrates[toPrintHere[ireac][0]]) + ' & \\\\'
  elif ns != Powers:
    thisReac = reacStart + ' & '+str(Forwardrates[toPrintHere[ireac][0]])+' & '+str(Backwardrates[toPrintHere[ireac][0]]) + ' & $*$ \\\\'
  elif nsProd != PowersProd:
    thisReac = reacStart + ' & '+str(Forwardrates[toPrintHere[ireac][0]])+' & '+str(Backwardrates[toPrintHere[ireac][0]]) + ' & $\dagger$ \\\\'

  ReacsGrouped.append(thisReac)

  diffCounter = diffCounter + 1

print ""
print "VARIABLES PRINTED WHERE THEY BELONG (not used in the article):"
ReacLens = [len(x) for x in ReacsGrouped]
for ireac in range(0,len(ReacsGrouped)):
  print ReacsGrouped[ireac]
  icheck = 0                                              # Check the numbers of 'mathbf's in the ReacsGrouped.
  while ReacsGrouped[ireac][icheck:].find('mathbf') > -1: # This is important as we want to get an idea how long the text is in PDF so we can divide to two rows if needed.
    ReacLens[ireac] = ReacLens[ireac] - 13                # Each mathbf means an extra 13 characters that don't affect the length.
    icheck = icheck + ReacsGrouped[ireac][icheck:].find('mathbf') + 1
  for ivariable in range(0,len(variables)):
    if variableNames[ivariable] in ReacsGrouped[ireac]:
      myText = '  '+variableNames[ivariable]+' $\in$ $\\{$'
      for ivar in range(0,len(variables[ivariable])):
        myText = myText + bcolors.OKBLUE + variables[ivariable][ivar] + bcolors.ENDC + ', '
      myText = myText[0:-2]+'$\\}$'
      print myText
      
ReacsLined = []
print ""
print ""
print "VARIABLES PRINTED IN THE END (Table 1)"
for ireac in range(0,len(ReacsGrouped)):
  if ReacLens[ireac] > 80:
    ind = ReacsGrouped[ireac].find('$\\right')
    ReacsLined.append(str(ireac+1)+' & '+ReacsGrouped[ireac][0:ind] + ' & & & \\\\')
    ReacsLined.append(' &   '+ReacsGrouped[ireac][ind:])
  else:
    ReacsLined.append(str(ireac+1)+' & '+ReacsGrouped[ireac])

for ireacper2 in range(0,(len(ReacsLined)+1)/2):
  ind = ReacsLined[ireacper2].find('\\\\')
  if ireacper2 + (len(ReacsLined)+1)/2 < len(ReacsLined):
    secondLine = ReacsLined[ireacper2 + (len(ReacsLined)+1)/2]
  else:
    secondLine = ' & \\\\'
  if ind > -1:
    print ReacsLined[ireacper2][0:ind] + ' '*(115-len(ReacsLined[ireacper2][0:ind]))+'& ' + secondLine
  else:
    print ReacsLined[ireacper2] + ' '*(115-len(ReacsLined[ireacper2]))+'& ' + secondLine
  

print ""
ivariable = 0
varLines = []
while ivariable < len(variables):
  if not 'X' in variableNames[ivariable]:
    print 'error, variableNames[ivariable]='+str(variableNames[ivariable])
  NvarsSameReac = 1
  while ivariable + NvarsSameReac < len(variables) and 'X' not in variableNames[ivariable+NvarsSameReac]:
    NvarsSameReac = NvarsSameReac + 1
  if NvarsSameReac == 1:
    myText = variableNames[ivariable]+' $\in$ $\\{$'
    for ivar in range(0,len(variables[ivariable])):
      myText = myText + variables[ivariable][ivar].replace('_','\_') + ', '
    myText = myText[0:-2]+'$\\}$'
    varLines.append(myText)
  else:
    myText = '('
    for ivar in range(0,NvarsSameReac):
      myText = myText + variableNames[ivariable+ivar] + ', '
    myText = myText[0:-2] + ') $\in$ $\\{$ ('
    for ivar in range(0,len(variables[ivariable])):
      myText = myText[0:-2]+' ('
      for ivarX in range(0,NvarsSameReac):
        myText = myText + variables[ivariable+ivarX][ivar].replace('_','\_') + ', '
      myText = myText[0:-2]+'), ('
    myText = myText[0:-3]+' $\\}$'
    varLines.append(myText)
  ivariable = ivariable + NvarsSameReac

varsLined = []
VarLens = [len(x) for x in varLines]
for ivar in range(0,len(varLines)):
  varLines[ivar] = varLines[ivar]+'\\\\'

for ivar in range(0,len(varLines)):
  myLine = varLines[ivar]
  while len(myLine) > 0:
    nMathbf = 0
    toDoLater = ''
    icheck = 0                                             # Check the numbers of 'mathbf's in the varLines.
    while myLine[icheck:].find('mathbf') > -1:             # This is important as we want to get an idea how long the text is in PDF so we can divide to two rows if needed.
      nMathbf = nMathbf + 1                                # Each mathbf means an extra 13 (here 18) characters that don't affect the length.
      icheck = icheck + myLine[icheck:].find('mathbf') + 1
    while len(myLine) > 72 + 13*nMathbf:
      s = myLine[:]
      ind = len(s)-1-(s[::-1].find(','))
      toDoLater = myLine[ind:]+toDoLater
      myLine = myLine[0:ind]
      nMathbf = 0
      icheck = 0                                             # Check the numbers of 'mathbf's in the varLines.
      while myLine[icheck:].find('mathbf') > -1:             # This is important as we want to get an idea how long the text is in PDF so we can divide to two rows if needed.
        nMathbf = nMathbf + 1                                # Each mathbf means an extra 13 (here 18) characters that don't affect the length.
        icheck = icheck + myLine[icheck:].find('mathbf') + 1
    if len(toDoLater) > 0 and toDoLater[0] == ',':
      varsLined.append(myLine+',\\\\')
      myLine = toDoLater[1:]
    else:
      varsLined.append(myLine)
      myLine = toDoLater


for ivarper2 in range(0,(len(varsLined)+1)/2):
  ind = varsLined[ivarper2].find('\\\\')
  if ivarper2 + (len(varsLined)+1)/2 < len(varsLined):
    secondLine = varsLined[ivarper2 + (len(varsLined)+1)/2]
  else:
    secondLine = '\\\\'
  if ind > -1:
    print varsLined[ivarper2][0:ind] + ' '*(96-len(varsLined[ivarper2][0:ind]))+' & ' + secondLine
  else:
    print varsLined[ivarper2] + ' '*(96-len(varsLined[ivarper2]))+' & ' + secondLine


      
