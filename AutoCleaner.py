import pyabf
from pyabf import abfWriter
import statistics as stat
import scipy.stats as sciStat
import scipy.ndimage
import scipy.signal as sig
import scipy.optimize as opt
import numpy as np
import os
import matplotlib.pyplot as plt
import csv

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

def wholeTraceMedFilt(signalToFilter):
    sweepArray = np.zeros((len(signalToFilter), len(signalToFilter[0])))
    for i, x in enumerate(signalToFilter):
        medFilteredSweep = sig.medfilt(x, kernel_size=5)
        sweepArray[i] = medFilteredSweep
    return sweepArray

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

# PRELIM: Linear Regression with no plot, across an entire trace
def isoLinReg(signal405, signal470):
    yArray = np.zeros((len(signal405), len(signal405[0][1000:-1250])))
    bumper470 = np.zeros((len(signal470), len(signal470[0][1000:-1250])))
    motionCorrectedSignal = np.zeros((len(signal470), len(signal470[0][1000:-1250])))
    for i, x in enumerate(signal405):
        yArray[i] = x[1000:-1250]
    flatY = yArray.reshape(-1)
    xArray = np.arange(len(flatY))
    line = sciStat.linregress(x=xArray, y=flatY)
    for i, x in enumerate(signal470):
        bumper470[i] = x[1000:-1250]
    convert1d = 0
    for i, x in np.ndenumerate(bumper470):
        motionCorrectedSignal[i[0]][i[1]] = x - (line.intercept + line.slope*xArray[convert1d])
        convert1d +=1
    return motionCorrectedSignal

# PRELIMINARY: Creates a linear regression of the 405 signal to create a predicted signal for motion (to subtract from the 470 signal)
def isoLinRegPlot(fileName, isosbesticChannel, chosenTrace, ratSide):
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

# PRELIMINARY Takes the numbers from deltaF or zCalc and fits a bi-exponential decay function to them
def doubleExpDecayFit(filteredSignal):
    flatSignal = filteredSignal.reshape(-1)
    xDecay = np.arange(0, len(flatSignal))
    yDecay = [i for i in flatSignal]
    p0 = (max(yDecay), -1, max(yDecay), -1, min(yDecay))
    popt, _ = opt.curve_fit(lambda t, a, b, c, d, e: (a * np.exp(b * t)) + (c * np.exp(d * t)) + e, xDecay, yDecay, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, -np.inf, 0, -np.inf, min(yDecay)], 
                                                                ub=[np.inf, 0, np.inf, 0, np.inf]),
                                                                maxfev=1000)
    a, b, c, d, e = popt[0], popt[1], popt[2], popt[3], popt[4]
    squaredDiffs = np.square(yDecay - (a * np.exp(b * xDecay)) + (c * np.exp(d * xDecay)) + e)
    squaredDiffsFromMean = np.square(yDecay - np.mean(yDecay))
    rSquared = 1 - (np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean))
    xDecFinal = np.linspace(np.min(xDecay), np.max(xDecay), len(filteredSignal))
    yDecFinal = (a * np.exp(b * xDecFinal)) + (c * np.exp(d * xDecFinal)) + e
    return yDecFinal, rSquared

def unbleachSignal(filteredSignal, decayFactor):
    unbleachedArray = np.zeros((len(filteredSignal), len(filteredSignal[0])))
    convert1d = 0
    for i, y in enumerate(filteredSignal):
        unbleachedArray[i] = y - decayFactor[convert1d]
        convert1d += 1
    return unbleachedArray

# Reads a CSV file containing all previous data for decay values and averages them, creating an array of average decay
def averageCSV():
    meanDecay405Left, meanDecay470Left, meanDecay405Right, meanDecay470Right = np.arange(0, 47750), np.arange(0, 47750), np.arange(0, 47750), np.arange(0, 47750)
    with open(os.path.join(os.getcwd(), 'Left405Decay.csv'), 'r') as csvfile:
        decayReader = csv.reader(csvfile)
        leftArray = np.zeros((len(decayReader), len(decayReader[0])))
        for i, row in enumerate(decayReader):
            leftArray[i] = row
        meanDecay405Left = leftArray.mean(axis=0)     
    with open(os.path.join(os.getcwd(), 'Right405Decay.csv'), 'r') as csvfile:
        decayReader = csv.reader(csvfile)
        rightArray = np.zeros((len(decayReader), len(decayReader[0])))
        for i, row in enumerate(decayReader):
            rightArray[i] = row
        meanDecay405Right = rightArray.mean(axis=0)
    with open(os.path.join(os.getcwd(), 'Left470Decay.csv'), 'r') as csvfile:
        decayReader = csv.reader(csvfile)
        leftArray = np.zeros((len(decayReader), len(decayReader[0])))
        for i, row in enumerate(decayReader):
            leftArray[i] = row
        meanDecay470Left = leftArray.mean(axis=0)     
    with open(os.path.join(os.getcwd(), 'Right470Decay.csv'), 'r') as csvfile:
        decayReader = csv.reader(csvfile)
        rightArray = np.zeros((len(decayReader), len(decayReader[0])))
        for i, row in enumerate(decayReader):
            rightArray[i] = row
        meanDecay470Right = rightArray.mean(axis=0)
    return meanDecay405Left, meanDecay470Left, meanDecay405Right, meanDecay470Right

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

def newCompleteProcessor(experimentFileName, baselineFileName, salineStatus):
    baseline470Left, baseline405Left, baseline470Right, baseline405Right= BaselineGet(baselineFileName)
    channelsLeft, channelsRight = [0,1], [4,5]
    subtract470Left, subtract405Left = baselineSubtractor(experimentFileName, baseline470Left, baseline405Left, channelsLeft)
    subtract470Right, subtract405Right = baselineSubtractor(experimentFileName, baseline470Right, baseline405Right, channelsRight)
    filtered405Left, filtered405Right = wholeTraceGauss(wholeTraceMedFilt(subtract405Left)), wholeTraceGauss(wholeTraceMedFilt(subtract405Right))
    filtered470Left, filtered470Right = wholeTraceGauss(wholeTraceMedFilt(subtract470Left)), wholeTraceGauss(wholeTraceMedFilt(subtract470Right))
    if salineStatus:
        decayFit405Left, _ = doubleExpDecayFit(filtered405Left)
        decayFit470Left, _ = doubleExpDecayFit(filtered470Left)
        decayFit405Right, _ = doubleExpDecayFit(filtered405Right)
        decayFit470Right, _ = doubleExpDecayFit(filtered470Right)
        unbleached405Left, unbleached405Right = unbleachSignal(filtered405Left, decayFit405Left), unbleachSignal(filtered405Right, decayFit405Right)
        unbleached470Left, unbleached470Right = unbleachSignal(filtered470Left, decayFit470Left), unbleachSignal(filtered470Right, decayFit470Right)
        with open('Left405Decay.csv', 'a') as csvfile:
            decayWriter = csv.writer(csvfile)
            decayWriter.writerow([i for i in decayFit405Left])
        with open('Right405Decay.csv', 'a') as csvfile:
            decayWriter = csv.writer(csvfile)
            decayWriter.writerow([i for i in decayFit405Right])
        with open('Left470Decay.csv', 'a') as csvfile:
            decayWriter = csv.writer(csvfile)
            decayWriter.writerow([i for i in decayFit470Left])
        with open('Right470Decay.csv', 'a') as csvfile:
            decayWriter = csv.writer(csvfile)
            decayWriter.writerow([i for i in decayFit470Right])
    elif not salineStatus:
        decayFit405Left, decayFit470Left, decayFit405Right, decayFit470Right = averageCSV()
        unbleached405Left, unbleached405Right = unbleachSignal(filtered405Left, decayFit405Left), unbleachSignal(filtered405Right, decayFit405Right)
        unbleached470Left, unbleached470Right = unbleachSignal(filtered470Left, decayFit470Left), unbleachSignal(filtered470Right, decayFit470Right)
        signalLeft, signalRight = isoLinReg(unbleached405Left, unbleached470Left), isoLinReg(unbleached405Right, unbleached470Right)
    unbleachedSignals = np.array[unbleached405Left.mean(axis=1), unbleached470Left.mean(axis=1), unbleached405Right.mean(axis=1), unbleached470Right.mean(axis=1)]
    signalLeft, signalRight = isoLinReg(unbleached405Left, unbleached470Left), isoLinReg(unbleached405Right, unbleached470Right)
    

    finalLeft, finalRight = np.zeros((len(signalLeft), len(signalLeft[0]))), np.zeros((len(signalRight), len(signalRight[0])))
    for i, x in enumerate(signalLeft):
        finalLeft[i] = x
    for i, x in enumerate(signalRight):
        finalRight[i] = x
    return finalLeft, finalRight, unbleachedSignals