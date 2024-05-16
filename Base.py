import FileSelector as fls
import Exceptions as exn
import AutoCleaner as acl
import pyabf
import os
import matplotlib.pyplot as plt

#TODO Fully automate the process of correcting data files produced by pClamp, according to the corrections laid out by Dr. Jose Eltit:
#      DONE Find average value of baseline traces
#      DONE Subtract baselines from every trace in a real file
#      DONE Gaussian filter the 405 channel
#      Take the ratio of 470/405
#      DONE Gaussian filter the ratio signal
#      Automate this process
#      Save the final corrections as .txt files that can be easily opened in clampfit 10 for later viewing (if necessary)
#      Do this using a GUI for ease of access for future students and trainees
#   
#      Eventually: Au 
#         
#TODO Data Analysis:
#      Export mean of traces to Excel File, along with corrections from pre-inj average, and bleaching (if available)
#      Define what an "event" is, then anaylze width, height, decay, frequency, and duration of event trains

#Chooses the file to be modified, and the file that baselines come from
userFile = fls.FileSelect("raw", "24508005.abf")
abf = pyabf.ABF(userFile)
userTrace = 25 # int(input("Which trace would you like to see?"))
userChannel = 4 # int(input("Which channel?"))
baselineFile = fls.FileSelect("raw", "24508006.abf")

# Finds baselines and subtracts them from channels 1, 2, 5, and 6 based on the chosen file
baselineSubL = acl.LBaselineGet(baselineFile)
baselineSubR = acl.RBaselineGet(baselineFile)
channelsLeft = [0,1]
channelsRight = [4,5]
subtractLeft = acl.baselineSubtractor(userFile, baselineSubL, channelsLeft)
subtractRight = acl.baselineSubtractor(userFile, baselineSubR, channelsRight)

# Gaussian filters the 405 channels
filteredLeft = acl.wholeTraceGauss(subtractLeft[1])
filteredRight = acl.wholeTraceGauss(subtractRight[1])

#Find ratio of 470/405 channels
ratioSignalLeft = acl.ratio470405(subtractLeft[0], filteredLeft)
ratioSignalRight = acl.ratio470405(subtractRight[0], filteredRight)

# Gaussian filters the ratio signal
finalSignalLeft = acl.wholeTraceGauss(ratioSignalLeft)
finalSignalRight = acl.wholeTraceGauss(ratioSignalRight)

abf.setSweep(sweepNumber= userTrace, channel= userChannel)
plt.plot(abf.sweepX[1000:-2000], abf.sweepY[1000:-2000], color="b", label="original")
plt.plot(abf.sweepX[1000:-2000], subtractRight[0][userTrace][1000:-2000], color="r", label= "Subtracted")
plt.plot(abf.sweepX[1000:-2000], ratioSignalRight[userTrace][1000:-2000], color="y", label= "Ratio-ed")
plt.plot(abf.sweepX[1000:-2000], finalSignalRight[userTrace][1000:-2000], color="g", label= "Processed")
# plt.plot(abf.sweepX[1000:-2000], filteredPlot[userTrace][1000:-2000], color="r", label= "Filtered")
plt.axis([0,15,0,6])
plt.legend()
plt.show()