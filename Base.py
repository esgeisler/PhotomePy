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
#      DONE Take the ratio of 470/405
#      DONE Gaussian filter the ratio signal
#      DONE Automate this process
#      Save the final corrections as .txt files that can be easily opened in clampfit 10 for later viewing (if necessary)
#      Do this using a GUI for ease of access for future students and trainees
#   
#      Eventually: Au 
#         
#TODO Data Analysis:
#      Export mean of traces to Excel File, along with corrections from pre-inj average, and bleaching (if available)
#      Define what an "event" is, then anaylze width, height, decay, frequency, and duration of event trains

#Chooses the file to be modified, and the file that baselines come from
dataType = input("Is the data file raw or processed?")
userFile = fls.FileSelect(dataType, input("What is the name of the file?"))
userTrace = int(input("Which trace would you like to see?"))
userChannel = int(input("Which channel?"))
abf = pyabf.ABF(userFile)
if dataType == "raw":
    baselineFile = fls.FileSelect("raw", input("What is the name of the file to be used for baseline?"))
decision = input("What shall we do today?")

if decision == "baseline":
    baselineSubL = acl.LBaselineGet(baselineFile)
    baselineSubR = acl.RBaselineGet(baselineFile)
    outputString = ("Left - 470: %.2f 405: %.2f\nRight - 470: %.2f 405: %.2f"%(baselineSubL[0], baselineSubL[1], baselineSubR[0], baselineSubR[1]))
    print(outputString)
elif decision == "process":
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
    # plt.plot(abf.sweepX[1000:-2000], abf.sweepY[1000:-2000], color="b", label="original")
    # plt.plot(abf.sweepX[1000:-2000], subtractRight[0][userTrace][1000:-2000], color="r", label= "Subtracted")
    # plt.plot(abf.sweepX[1000:-2000], ratioSignalRight[userTrace][1000:-2000], color="y", label= "Ratio-ed")
    # plt.plot(abf.sweepX[1000:-2000], finalSignalRight[userTrace][1000:-2000], color="g", label= "Processed")
    # plt.plot(abf.sweepX[1000:-2000], filteredPlot[userTrace][1000:-2000], color="r", label= "Filtered")
    # plt.axis([0,15,0,6])
    # plt.legend()
    # plt.show()