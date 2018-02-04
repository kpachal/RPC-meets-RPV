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
      highVal = val
  return lowVal, highVal

# Linear extrapolation between value at ylow and value at yhigh
def getExtrapolation(low, ylow, high, yhigh, val) :



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

# Check if di-b can improve on any of these results.
# Here, take the masses from Gaussian limits
# and extrapolate linearly between acceptance and
# efficiency points.
diB_limitDict = {}

# For extrapolating later
knownVals_accTrigEff = sorted([eval(i) for i in dict_acceptances_times_btriggereff.keys()])
print knownVals_accTrigEff
knownVals_taggingEff = sorted([eval(i) for i in dict_taggingEff.keys()])
print knownVals_taggingEff

for mass in sorted(dict_DiBLimits.keys()) :

  # Need a cross section prediction to be useful.
  if not mass in xSecDict.keys() :
    continue
  diB_limitDict[mass] = {}

  # Get excluded values from Gaussian limits
  exclSigma = dict_DiBLimits[mass][widthVal]
  diB_limitDict[mass]["observed_limit_pb"] = exclSigma

  # Now we need a prediction for each value of the coupling
  sigma_at1 = xSecDict[mass]
  for lambdaVal in lambdaVals :

    # Now find numbers to compare to:
    # in this case, xsec * (acceptance * trigger efficiency) * (b-tagging efficiency) * BR
    
    # xsec
    sigma_thisVal = lambdaVal*lambdaVal*sigma_at1

    # (acceptance * trigger efficiency)
    massDown,massUp = getLowHighVals(knownVals_accTrigEff,eval(mass))

    logVal = logVals[lambdaVals.index(lambdaVal)]
    BR = rpcrpv_udd_mapping.log10udd_to_BR(logVal,eval(mass))
    if BR < 0 :
      print "Evaluating for",logVal,eval(mass)
      print "   has BR",BR
      continue

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

# Cross check plot 1: do all my cross sections make sense?
#xsecs = []
#for mass in sorted([eval(i) for i in xSecDict.keys()]) :
#  xsecs.append(xSecDict["{0}".format(mass)])
#plt.plot(xsecs)
#plt.ylabel("cross section [pb]")
#plt.xlabel("stop mass")
#plt.savefig("crossSections_reconstructed.eps")
#plt.yscale('log')
#plt.savefig("crossSections_reconstructed_log.eps")
#plt.clf()

# Cross check plot 2: does my excluded region make sense
# given results shown in Angelo's paper?

# Make axes
#delta = 100
#x = sorted(limitDict.keys())
#y = lambdaVals
#X = []
#Y = []
#Z = []
#for mass in x:
#  for coupling in y :
#    X.append(eval(mass))
#    Y.append(coupling)
#    lambdaString = "{0}".format(coupling).replace(".","p")
#    Z.append(int(limitDict[mass]["isExcluded_lambda{0}".format(lambdaString)]))

# Linearly interpolate to get smooth data
#xi = np.linspace(500,2800,30)
#yi = np.linspace(0.01,1.1,10)

#points = np.random.rand(1000, 2)
#zi = plt.mlab.griddata(X, Y, Z, xi, yi, interp='nn')

#print xi
#print yi
#print zi
#plt.yscale('log')
#plt.ylim([0.01,1.1])
#CS = plt.contour(xi,yi,zi, 1)#, linewidths=0.5, colors='k')
#CS = plt.contourf(xi,yi,zi, 1)#,
                  #vmax=abs(Z).max(), vmin=-abs(Z).max())
#plt.title('2D exclusion contour')
#plt.show()
