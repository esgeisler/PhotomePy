import pyabf
from pyabf import abfWriter
import statistics as stat
import scipy.stats as sciStat
import scipy.ndimage
import numpy as np
import os
import matplotlib.pyplot as plt

# Gets baseline information from 1 min-long recording data taken after trial from the "left" side of the room - channels 1 and 2
def BaselineGet(FileName):
    abf = pyabf.ABF(FileName)
    Channel470Left, Channel405Left, Channel470Right, Channel405Right = 0, 1, 4, 5
    channels = [Channel470Left, Channel405Left, Channel470Right, Channel405Right]
    sweepArrays = {Channel470Left:0.0, Channel405Left:0.0, 
                   Channel470Right:0.0, Channel405Right:0.0}

    for c in channels:
        abf.setSweep(sweepNumber= 0, channel= c)
        sweepArrays[c] = stat.mean(abf.sweepY)

    mean470Left, mean405Left, mean470Right, mean405Right = sweepArrays[Channel470Left], sweepArrays[Channel405Left], sweepArrays[Channel470Right], sweepArrays[Channel405Right]
    return mean470Left, mean405Left, mean470Right, mean405Right

#Function that subtracts baseline from baseline getters from a chosen channel, returns those as full channel trace dictionaries
def baselineSubtractor(fileName, baseline470, baseline405, channelsToSubtract):
    abf = pyabf.ABF(fileName)
    abf.setSweep(0)
    sweepArray470, sweepArray405 = np.zeros((len(abf.sweepList), len(abf.sweepX))), np.zeros((len(abf.sweepList), len(abf.sweepX)))
    #470
    for i in abf.sweepList:
        abf.setSweep(i, channel= channelsToSubtract[0])
        sweepArray470[i] = [x - baseline470 for x in abf.sweepY]
        abf.setSweep(i, channel= channelsToSubtract[1])
        sweepArray405[i] = [x - baseline405 for x in abf.sweepY]
    return sweepArray470, sweepArray405

# Gaussian filters an entire channel with a 40 Hz cutoff freq., as above.
def wholeTraceGauss(signalToFilter):
    sweepArray = np.zeros((len(signalToFilter), len(signalToFilter[0])))
    for i, x in enumerate(signalToFilter):
        filteredSweep = scipy.ndimage.gaussian_filter1d(x, sigma = 16)
        sweepArray[i] = filteredSweep
    return sweepArray

# Divides two channels (470nm/405nm) in one file and returns a complete channel dictionary.
def ratio470405(signal470, signal405):
    ratioSignal = np.zeros((len(signal470), len(signal470[0])))
    for i, x in enumerate(signal470):
        ratioSignal[i] = x/signal405[i]
    return ratioSignal

# PRELIMINARY: Creates a linear regression of the 405 signal to create a predicted signal for motion (to subtract from the 470 signal)
def isoLinReg(fileName, isosbesticChannel, chosenTrace, ratSide):
    abf = pyabf.ABF(fileName)
    abf.setSweep(chosenTrace)
    yArray, xArray = np.zeros((len(abf.sweepList), len(abf.sweepX))), np.zeros((len(abf.sweepList), len(abf.sweepX)))
    for i in abf.sweepList:
        abf.setSweep(i, channel=isosbesticChannel)
        yArray[i] = abf.sweepY
        xArray[i] = abf.sweepX
    samplingFreqMSec = abf.dataPointsPerMs + (1/3)
    samplingFreqSec = samplingFreqMSec * 1000
    seconds = np.arange(0, 47750, samplingFreqSec)
    line = sciStat.linregress(x=xArray[chosenTrace][1000:-1250], y=yArray[chosenTrace][1000:-1250])
    fig = plt.figure()
    motionFig = fig.add_subplot()
    motionFig.plot(yArray[chosenTrace][1000:-1250], 'c.', label="Raw Data")
    motionFig.plot(line.intercept + line.slope*xArray[chosenTrace][1000:-1250], 'k', label=f"Line of Best Fit (R²): {line.rvalue**2:.3f}")
    motionFig.set_title(ratSide)
    plt.xticks(ticks=seconds, labels=range(len(seconds)))
    plt.xlabel("Time (s)")
    plt.ylabel("Fluorescence (AU)")
    plt.legend()
    plt.minorticks_on()
    plt.show()

# Saves cleaned trace file as an ABF file for later viewing in ClampFit 10 "Processed Data", "%s Rat %s Processed Data.abf"%(self.abfDate.strftime("%Y-%m-%d")
def tExport(processedTrace, ratName, experimentDate):
    abfWriter.writeABF1(sweepData= processedTrace, filename= os.path.join(os.getcwd(), "Processed Data", "%s Rat %s Processed Data.abf"%(experimentDate.strftime("%Y-%m-%d"), ratName)), units="V", sampleRateHz= 3333.33)

# Compilation of the other functions in this file
# Gets a baseline, subtracts it from a main signal, gaussian filters the 405 channel, then takes the ratio of the 470/405 channels. Finally, filters the ratio signal
def completeProcessor(experimentFileName, baselineFileName):
    if not os.path.exists(experimentFileName):
        raise FileNotFoundError("main")
    elif not os.path.exists(baselineFileName):
        raise FileNotFoundError("baseline")
    abf = pyabf.ABF(experimentFileName)
    baseline470Left, baseline405Left, baseline470Right, baseline405Right= BaselineGet(baselineFileName)
    channelsLeft, channelsRight = [0,1], [4,5]
    subtract470Left, subtract405Left = baselineSubtractor(experimentFileName, baseline470Left, baseline405Left, channelsLeft)
    subtract470Right, subtract405Right = baselineSubtractor(experimentFileName, baseline470Right, baseline405Right, channelsRight)
    filteredLeft, filteredRight = wholeTraceGauss(subtract405Left), wholeTraceGauss(subtract405Right)
    ratioSignalLeft, ratioSignalRight = ratio470405(subtract470Left, filteredLeft), ratio470405(subtract470Right, filteredRight)
    signalLeft, signalRight = wholeTraceGauss(ratioSignalLeft), wholeTraceGauss(ratioSignalRight)

    finalLeft, finalRight = np.zeros((len(signalLeft), len(signalLeft[0]) - 2250)), np.zeros((len(signalRight), len(signalRight[0]) - 2250))
    for i, x in enumerate(signalLeft):
        finalLeft[i] = x[1000:-1250]
    for i, x in enumerate(signalRight):
        finalRight[i] = x[1000:-1250]
    return finalLeft, finalRight