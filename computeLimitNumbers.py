from dict_CrossSections import *
from dict_acceptances import *
from dict_DijetLimits_MoriondPaper import *
from dict_TLALimits_ATLAS_CONF_2016_030 import *

import matplotlib.pyplot as plt
import numpy as np

# These are the lambda'' values to scan.
lambdaVals = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]

# This is the width to use. 0.07 is conservative
# according to the estimate in the paper.
widthVal = "0.07"

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
for mass in limitDict.keys() :

  sigma_at1 = xSecDict[mass]
  
  # Now we need a prediction for each value of the coupling
  for lambdaVal in lambdaVals :
    lambdaString = "{0}".format(lambdaVal).replace(".","p")
    sigma_thisVal = lambdaVal*lambdaVal*sigma_at1

    # Prediction is cross section times acceptance times branching ratio
    # Here, branching ratio is 1
    BR = 1.0
    acc = dict_acceptances[mass]
    prediction = sigma_thisVal * acc * BR
    name = "predicted_limit_lambda{0}_pb".format(lambdaString)
    limitDict[mass][name] = prediction

    # A point is excluded if the observed limit is lower
    # than the predicted limit
    if prediction > limitDict[mass]["observed_limit_pb"] :
      limitDict[mass]["isExcluded_lambda{0}".format(lambdaString)] = True
    else :
      limitDict[mass]["isExcluded_lambda{0}".format(lambdaString)] = False

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
xsecs = []
for mass in sorted([eval(i) for i in xSecDict.keys()]) :
  xsecs.append(xSecDict["{0}".format(mass)])
plt.plot(xsecs)
plt.ylabel("cross section [pb]")
plt.xlabel("stop mass")
plt.savefig("crossSections_reconstructed.eps")
plt.yscale('log')
plt.savefig("crossSections_reconstructed_log.eps")
plt.clf()

# Cross check plot 2: does my excluded region make sense
# given results shown in Angelo's paper?

# Make axes
delta = 100
x = sorted(limitDict.keys())
y = lambdaVals
X = []
Y = []
Z = []
for mass in x:
  for coupling in y :
    X.append(eval(mass))
    Y.append(coupling)
    lambdaString = "{0}".format(coupling).replace(".","p")
    Z.append(int(limitDict[mass]["isExcluded_lambda{0}".format(lambdaString)]))

# Linearly interpolate to get smooth data
xi = np.linspace(500,2800,30)
yi = np.linspace(0.01,1.1,10)

#points = np.random.rand(1000, 2)
zi = plt.mlab.griddata(X, Y, Z, xi, yi, interp='nn')

print xi
print yi
print zi
plt.yscale('log')
plt.ylim([0.01,1.1])
CS = plt.contour(xi,yi,zi, 1)#, linewidths=0.5, colors='k')
CS = plt.contourf(xi,yi,zi, 1)#,
                  #vmax=abs(Z).max(), vmin=-abs(Z).max())
plt.title('2D exclusion contour')
#plt.show()
