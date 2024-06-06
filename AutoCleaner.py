import pyabf
from pyabf import abfWriter
import statistics as stat
import scipy.ndimage
import numpy as np
from datetime import datetime

# Gets baseline information from 1 min-long recording data taken after trial from the "left" side of the room - channels 1 and 2
def LBaselineGet(FileName):
    abf = pyabf.ABF(FileName)
    sweepDict470 = {}
    Channel470 = 0
    sweepDict405 = {}
    Channel405 = 1

    abf.setSweep(sweepNumber = 0, channel = Channel470)
    for sweep470 in abf.sweepList:
        fluor470 = stat.mean(abf.sweepY)
        sweepDict470[sweep470] = fluor470
    mean470 = stat.mean(sweepDict470.values())
    
    abf.setSweep(sweepNumber = 0, channel = Channel405)
    for sweep405 in abf.sweepList:
        fluor405 = stat.mean(abf.sweepY)
        sweepDict405[sweep405] = fluor405
    mean405 = stat.mean(sweepDict405.values())

    finalSubtraction = [mean470, mean405]
    return finalSubtraction

# Gets baseline information from 1 min-long recording data taken after trial from the "right" side of the room - channels 5 and 6
def RBaselineGet(FileName):
    abf = pyabf.ABF(FileName)
    sweepDict470 = {}
    Channel470 = 4
    sweepDict405 = {}
    Channel405 = 5

    abf.setSweep(sweepNumber = 0, channel = Channel470)
    for sweep470 in abf.sweepList:
        fluor470 = stat.mean(abf.sweepY)
        sweepDict470[sweep470] = fluor470
    mean470 = stat.mean(sweepDict470.values())
    
    abf.setSweep(sweepNumber = 0, channel = Channel405)
    for sweep405 in abf.sweepList:
        fluor405 = stat.mean(abf.sweepY)
        sweepDict405[sweep405] = fluor405
    mean405 = stat.mean(sweepDict405.values())

    finalSubtraction = [mean470, mean405]
    return finalSubtraction

#Function that subtracts baseline from baseline getters from a chosen channel, returns those as full channel trace dictionaries
def baselineSubtractor(fileName, baselines, channelsToSubtract):
    abf = pyabf.ABF(fileName)
    sweepDict470 = {}
    sweepDict405 = {}
    if type(baselines) is not list:
        raise TypeError("Values for baseline not a list.")
    #470
    for sweeps in abf.sweepList:
        abf.setSweep(sweeps, channel= channelsToSubtract[0])
        sweepDict470[sweeps] = [x - baselines[0] for x in abf.sweepY]
    #405
    for sweeps in abf.sweepList:
        abf.setSweep(sweeps, channel= channelsToSubtract[1])
        sweepDict405[sweeps] = [x - baselines[1] for x in abf.sweepY]
    subtractedSweeps = [sweepDict470, sweepDict405]
    return subtractedSweeps

# Gaussian filters a single trace with a 40 Hz cutoff frequency, based on pClamp documentation and Calquhon & Sigworth (1995)
def gaussianFilter(FileName, filterChannel, filterSweep):
    abf = pyabf.ABF(FileName)
    abf.setSweep(sweepNumber= filterSweep, channel= filterChannel)
    sweepList = scipy.ndimage.gaussian_filter1d(abf.sweepY, sigma= 16)
    return sweepList

# Gaussian filters an entire channel with a 40 Hz cutoff freq., as above.
def wholeTraceGauss(signalToFilter):
    sweepDict = {}
    for sweeps in signalToFilter:
            filteredSweep = scipy.ndimage.gaussian_filter1d(signalToFilter[sweeps], sigma = 16)
            sweepDict[sweeps] = filteredSweep
    return sweepDict

# Divides two channels (470nm/405nm) in one file and returns a complete channel dictionary.
def ratio470405(signal470, signal405):
    indexRange = len(signal470)
    ratioSignal = {}
    for i in range(indexRange):
        ratioSignal[i] = signal470[i]/signal405[i]
    return ratioSignal

# Saves cleaned trace file as an ABF file for later viewing in ClampFit 10
def tExport(processedTrace, ratName):
    trace = processedTrace.values()
    arrayList = list(trace)
    array = np.array(arrayList)
    abfWriter.writeABF1(sweepData= array, filename= "Rat %s Processed Data %s.abf"%(ratName, datetime.today().strftime('%Y-%m-%m')), units="V", sampleRateHz= 3333.33)

# Compilation of the other functions in this file
# Gets a baseline, subtracts it from a main signal, gaussian filters the 405 channel, then takes the ratio of the 470/405 channels. Finally, filters the ratio signal
def completeProcessor(experimentFileName, baselineFileName):
    abf = pyabf.ABF(experimentFileName)
    baselineSubL = LBaselineGet(baselineFileName)
    baselineSubR = RBaselineGet(baselineFileName)
    channelsLeft = [0,1]
    channelsRight = [4,5]
    subtractLeft = baselineSubtractor(experimentFileName, baselineSubL, channelsLeft)
    subtractRight = baselineSubtractor(experimentFileName, baselineSubR, channelsRight)
    filteredLeft = wholeTraceGauss(subtractLeft[1])
    filteredRight = wholeTraceGauss(subtractRight[1])
    ratioSignalLeft = ratio470405(subtractLeft[0], filteredLeft)
    ratioSignalRight = ratio470405(subtractRight[0], filteredRight)
    finalSignalLeft = wholeTraceGauss(ratioSignalLeft)
    finalSignalRight = wholeTraceGauss(ratioSignalRight)
    return finalSignalLeft, finalSignalRight