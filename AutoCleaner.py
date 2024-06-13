import pyabf
from pyabf import abfWriter
import statistics as stat
import scipy.ndimage
import numpy as np
from datetime import datetime

# Gets baseline information from 1 min-long recording data taken after trial from the "left" side of the room - channels 1 and 2
def BaselineGet(FileName):
    abf = pyabf.ABF(FileName)
    sweepDict470Left = np.zeros(len(abf.sweepList))
    sweepDict405Left = np.zeros(len(abf.sweepList))
    sweepArray470Right = np.zeros(len(abf.sweepList))
    sweepArray405Right = np.zeros(len(abf.sweepList))
    Channel470Left = 0
    Channel405Left = 1
    Channel470Right = 4
    Channel405Right = 5

    abf.setSweep(sweepNumber= 0, channel= Channel470Left)
    for sweep470 in abf.sweepList:
        fluor470 = stat.mean(abf.sweepY)
        sweepDict470Left[sweep470] = fluor470
    mean470Left = stat.mean(sweepDict470Left)
    
    abf.setSweep(sweepNumber= 0, channel= Channel405Left)
    for sweep405 in abf.sweepList:
        fluor405 = stat.mean(abf.sweepY)
        sweepDict405Left[sweep405] = fluor405
    mean405Left = stat.mean(sweepDict405Left)
    
    abf.setSweep(sweepNumber= 0, channel= Channel470Right)
    for sweep470 in abf.sweepList:
        fluor470 = stat.mean(abf.sweepY)
        sweepArray470Right[sweep470] = fluor470
    mean470Right = stat.mean(sweepArray470Right)
    
    abf.setSweep(sweepNumber= 0, channel= Channel405Right)
    for sweep405 in abf.sweepList:
        fluor405 = stat.mean(abf.sweepY)
        sweepArray405Right[sweep405] = fluor405
    mean405Right = stat.mean(sweepArray405Right)
    return mean470Left, mean405Left, mean470Right, mean405Right

#Function that subtracts baseline from baseline getters from a chosen channel, returns those as full channel trace dictionaries
def baselineSubtractor(fileName, baseline470, baseline405, channelsToSubtract):
    abf = pyabf.ABF(fileName)
    abf.setSweep(0)
    sweepArray470 = np.zeros((len(abf.sweepList), len(abf.sweepX)))
    sweepArray405 = np.zeros((len(abf.sweepList), len(abf.sweepX)))
    #470
    for sweeps in abf.sweepList:
        abf.setSweep(sweeps, channel= channelsToSubtract[0])
        sweepArray470[sweeps] = [x - baseline470 for x in abf.sweepY]
    #405
    for sweeps in abf.sweepList:
        abf.setSweep(sweeps, channel= channelsToSubtract[1])
        sweepArray405[sweeps] = [x - baseline405 for x in abf.sweepY]
    return sweepArray470, sweepArray405

# Gaussian filters an entire channel with a 40 Hz cutoff freq., as above.
def wholeTraceGauss(signalToFilter):
    sweepArray = np.zeros((len(signalToFilter), len(signalToFilter[0])))
    for x in range(len(signalToFilter)):
            filteredSweep = scipy.ndimage.gaussian_filter1d(signalToFilter[x], sigma = 16)
            sweepArray[x] = filteredSweep
    return sweepArray

# Divides two channels (470nm/405nm) in one file and returns a complete channel dictionary.
def ratio470405(signal470, signal405):
    ratioSignal = np.zeros((len(signal470), len(signal470[0])))
    for x in range(len(signal470)):
        ratioSignal[x] = signal470[x]/signal405[x]
    return ratioSignal

# Saves cleaned trace file as an ABF file for later viewing in ClampFit 10
def tExport(processedTrace, ratName):
    abfWriter.writeABF1(sweepData= processedTrace, filename= "Rat %s Processed Data %s.abf"%(ratName, datetime.today().strftime('%Y-%m-%m')), units="V", sampleRateHz= 3333.33)

# Compilation of the other functions in this file
# Gets a baseline, subtracts it from a main signal, gaussian filters the 405 channel, then takes the ratio of the 470/405 channels. Finally, filters the ratio signal
def completeProcessor(experimentFileName, baselineFileName):
    abf = pyabf.ABF(experimentFileName)
    baseline470Left, baseline405Left, baseline470Right, baseline405Right= BaselineGet(baselineFileName)
    channelsLeft = [0,1]
    channelsRight = [4,5]
    subtract470Left, subtract405Left = baselineSubtractor(experimentFileName, baseline470Left, baseline405Left, channelsLeft)
    subtract470Right, subtract405Right = baselineSubtractor(experimentFileName, baseline470Right, baseline405Right, channelsRight)
    filteredLeft = wholeTraceGauss(subtract405Left)
    filteredRight = wholeTraceGauss(subtract405Right)
    ratioSignalLeft = ratio470405(subtract470Left, filteredLeft)
    ratioSignalRight = ratio470405(subtract470Right, filteredRight)
    signalLeft = wholeTraceGauss(ratioSignalLeft)
    signalRight = wholeTraceGauss(ratioSignalRight)

    finalLeft = np.zeros((len(signalLeft), len(signalLeft[0]) - 2100))
    finalRight = np.zeros((len(signalRight), len(signalRight[0]) - 2100))
    for x in range(len(signalLeft)):
        finalLeft[x] = signalLeft[x][850:-1250]
    for x in range(len(signalLeft)):
        finalRight[x] = signalRight[x][850:-1250]
    return finalLeft, finalRight