from dict_CrossSections import *
from dict_acceptances import *
from dict_DijetLimits_MoriondPaper import *
from dict_TLALimits_ATLAS_CONF_2016_030 import *

# di-B results
from dict_acceptanceTimesTrigEff_DiBJet import *
from dict_taggingEff_DiBJet import *
from dict_DiBJetLimits_ATLAS_CONF_2016_060 import *

import rpcrpv_udd_mapping

#import matplotlib.pyplot as plt
import numpy as np

###############################
## Settings

# These are the lambda'' values to scan.
# My original values
#lambdaVals = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
# Values matching Javier's BR table
logVals = [-4.0,-3.0,-2.0,-1.5,-1.2,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,-0.0, 0.2]
lambdaVals = []
for logVal in logVals :
  lambdaVals.append(pow(10.0,logVal))

# This is the width to use. 0.07 is conservative
# according to the estimate in the paper.
widthVal = "0.07"

###############################
# Some helper functions for extrapolations

# Find values I'm in between
def getLowHighVals(list, val) :
  lowVal = -1
  highVal = -1
  for item in sorted(list) :
    if item < val :
      lowVal = item
    if item > val and highVal < 0 :
      highVal = item
  return lowVal, highVal

# Linear extrapolation between value at ylow and value at yhigh
def getExtrapolation(low, ylow, high, yhigh, val) :

  # y = mx + b
  slope = (yhigh - ylow)/(float(high)-float(low))
  rise = slope * (float(val)-float(low))

  yval = ylow + rise

  #print "for",val,"between",ylow,"at",low,"and",yhigh,"at",high,"found",yval
  return yval

###############################
# Calculate exclusions

# Make a dictionary in which to organise our information
limitDict = {}

# The masses we can use are restricted by those
# for which we have acceptance values. So use that to set the keys.
for mass in dict_acceptances :

  acc = dict_acceptances[mass]
  if acc < 0 : continue

  # TEMP: until we get a response from Angelo,
  # we also only have cross sections up to 3 TeV.
  if eval(mass) > 3000 : continue

  # So now any mass surviving is one we can use.
  limitDict[mass] = {}

# Things we need to store:
# Excluded cross section
for mass in limitDict.keys() :

  exclSigma_dijet = 1e10
  if mass in dict_DijetLimits.keys() :
    exclSigma_dijet = dict_DijetLimits[mass][widthVal]/1000.0 # dijet limits are fb not pb

  exclSigma_TLA = 1e10
  if mass in dict_TLALimits.keys() :
    exclSigma_TLA = dict_TLALimits[mass][widthVal]
  
  limit = min(exclSigma_dijet,exclSigma_TLA)

  limitDict[mass]["observed_limit_pb"] = limit

# Theoretical cross section
for mass in sorted(limitDict.keys()) :

  sigma_at1 = xSecDict[mass]
  
  # Now we need a prediction for each value of the coupling
  for lambdaVal in lambdaVals :
    #lambdaString = "{0}".format(lambdaVal).replace(".","p")
    sigma_thisVal = lambdaVal*lambdaVal*sigma_at1

    # Prediction is cross section times acceptance times branching ratio
    # Basic version: branching ratio is 1
    #BR = 1.0
    # With new code from Javier, can get actual BR
    logVal = logVals[lambdaVals.index(lambdaVal)]
    BR = rpcrpv_udd_mapping.log10udd_to_BR(logVal,eval(mass))
    if BR < 0 :
      print "Evaluating for",logVal,eval(mass)
      print "   has BR",BR
      continue
    
    acc = dict_acceptances[mass]
    prediction = sigma_thisVal * acc * BR
    name = "predicted_limit_log10lambda{0}_pb".format(logVal)
    limitDict[mass][name] = prediction

    # A point is excluded if the observed limit is lower
    # than the predicted limit
    if prediction > limitDict[mass]["observed_limit_pb"] :
      limitDict[mass]["isExcluded_log10lambda{0}".format(logVal)] = True
    else :
      limitDict[mass]["isExcluded_log10lambda{0}".format(logVal)] = False

# Write results to a text file for Max and Javier.
# As I proceed, check results by printing them out.
with open("dijet_exclusions.txt", 'w') as outfile :
  outfile.write("limit_dict = {\n")
  for mass in sorted([eval(i) for i in limitDict.keys()]) :
    massread = "{0}".format(mass)
    print mass, ":"
    outfile.write("{0} : {1}\n".format(mass,"{"))
    for item in sorted(limitDict[massread].keys()) :
      print "\t",item,":",limitDict[massread][item]
      outfile.write("\t{0} : {1},\n".format(item, limitDict[massread][item]))
    outfile.write("\t},\n")
  outfile.write("}")

# Reformat the results the way they currently want
with open("dijet_exclusions.json", 'w') as outfile :

  outfile.write('[\n')
  for mass in limitDict.keys() :
    thisDict = limitDict[mass]
    obsLimit = thisDict["observed_limit_pb"]
    for logVal in logVals :
      predLimit = thisDict["predicted_limit_log10lambda{0}_pb".format(logVal)]
      ratio = predLimit/obsLimit
      outfile.write('  {\n')
      outfile.write('    "lambdapp": {0},\n'.format(logVal))
      outfile.write('    "mstop": {0},\n'.format(mass))
      outfile.write('    "CLsexp": {0},\n'.format(ratio))
      outfile.write('    "failedstatus": 0,\n')
      outfile.write('    "CLs": {0},\n'.format(ratio))
      outfile.write('  },\n')

  outfile.write(']')




