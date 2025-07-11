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
import pandas as pd
from pathlib import Path
import yaml

pd.options.display.float_format = "{:,.10f}".format

with open("config.yaml") as c:
    userConfig = yaml.safe_load(c)
    leftChannels = userConfig["GENERAL"]["left_rat_channels"]
    rightChannels = userConfig["GENERAL"]["right_rat_channels"]
    start =  userConfig["EVENT_HANDLING"]["trace_start_offset"]
    end =  userConfig["EVENT_HANDLING"]["trace_end_offset"]

# Gets baseline information from 1 min-long recording data taken after trial from the "left" side of the room - channels 1 and 2
def baselineGet(fileName):
    abf = pyabf.ABF(fileName)
    Channel470Left, Channel405Left = leftChannels
    Channel470Right, Channel405Right = rightChannels
    channels = [Channel470Left, Channel405Left, Channel470Right, Channel405Right]
    sweepArrays = {Channel470Left:0.0, Channel405Left:0.0, 
                   Channel470Right:0.0, Channel405Right:0.0}

    for c in channels:
        abf.setSweep(sweepNumber= 0, channel= c)
        sweepArrays[c] = stat.mean(abf.sweepY[start:end])

    mean470Left, mean405Left, mean470Right, mean405Right = sweepArrays[Channel470Left], sweepArrays[Channel405Left], sweepArrays[Channel470Right], sweepArrays[Channel405Right]
    return mean470Left, mean405Left, mean470Right, mean405Right

#Function that subtracts baseline from baseline getters from a chosen channel, returns those as full channel trace dictionaries
def baselineSubtractor(fileName, baseline470, baseline405, channelsToSubtract):
    abf = pyabf.ABF(fileName)
    abf.setSweep(0)
    sweepArray470, sweepArray405 = np.zeros((len(abf.sweepList), len(abf.sweepX) - 2250)), np.zeros((len(abf.sweepList), len(abf.sweepX) - 2250))
    #470
    for i in abf.sweepList:
        abf.setSweep(i, channel= channelsToSubtract[0])
        sweepArray470[i] = abf.sweepY[start:end] - baseline470 
        abf.setSweep(i, channel= channelsToSubtract[1])
        sweepArray405[i] = abf.sweepY[start:end] - baseline405
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
    ratioSignal = signal470/signal405
    return ratioSignal

# PRELIM: Linear Regression with no plot, across an entire trace
def isoLinReg(signal405, signal470):
    motionCorrectedSignal = np.zeros((len(signal470), len(signal470[0])))
    flatY = signal405.reshape(-1)
    xArray = np.arange(len(flatY))
    line = sciStat.linregress(x=xArray, y=flatY)
    convert1d = 0
    for i, x in np.ndenumerate(signal470):
        motionCorrectedSignal[i[0]][i[1]] = x - (line.intercept + line.slope*xArray[convert1d])
        convert1d +=1
    return motionCorrectedSignal

#PRELIM: Estimates the correlation between the 405/470 channels to determine the amount of signal due to movement.
def isoLinRegPlot(fileName, isosbesticChannel, activeChannel, chosenTrace, ratSide):
    abf = pyabf.ABF(fileName)
    abf.setSweep(chosenTrace, channel=isosbesticChannel)
    isosbesticY = abf.sweepY[start:end]
    abf.setSweep(chosenTrace, activeChannel)
    activeY = abf.sweepY[start:end]
    line = sciStat.linregress(x=isosbesticY, y=activeY)
    fig = plt.figure()
    motionFig = fig.add_subplot()
    motionFig.scatter(isosbesticY, activeY)
    motionFig.plot(isosbesticY, line.intercept + line.slope*isosbesticY, 'k', label=f"Line of Best Fit (R²): {line.rvalue**2:.3f}")
    motionFig.set_title(ratSide)
    #plt.xticks(ticks=seconds, labels=range(len(seconds)))
    plt.xlabel("Isosbestic Signal (V)")
    plt.ylabel("Active Signal (V)")
    plt.legend()
    plt.minorticks_on()
    plt.show()

# PRELIMINARY Takes the numbers from deltaF or zCalc and fits a bi-exponential decay function to them
def doubleExpDecayFit(filteredSignal):
    yDecay = filteredSignal.reshape(-1)
    xDecay = np.arange(0, len(yDecay))
    p0 = (max(yDecay), -1, max(yDecay)/2, -1, min(yDecay))
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
    # decayFig = plt.figure()
    # fig = decayFig.add_subplot()
    # fig.plot(xDecFinal, yDecFinal)
    # fig.plot(xDecay, yDecay)
    # plt.show()
    return yDecFinal, rSquared

def doubleExpDecaySingleTrace(filteredSignalArray, samplingFreq):
    yDecFinalArray = np.zeros((len(filteredSignalArray), len(filteredSignalArray[0])))
    rSquaredArray = np.zeros((len(filteredSignalArray), len(filteredSignalArray[0])))
    for i, _ in enumerate(filteredSignalArray):
        yDecay = filteredSignalArray[i]
        xDecay = np.arange(0, len(yDecay))
        decayToY = np.median(yDecay[-500:])
        startY = np.median(yDecay[:500])
        heightDiff = abs(startY - decayToY)
        p0 = (heightDiff/2, -1, heightDiff/2, -1, decayToY)
        popt, _ = opt.curve_fit(lambda t, a, b, c, d, e: (a * np.exp(b * t)) + (c * np.exp(d * t)) + e, xDecay/samplingFreq, yDecay, p0=p0, 
                                bounds=opt.Bounds(lb=[0         ,-np.inf, 0         , -np.inf, decayToY], 
                                                  ub=[heightDiff, 0     , heightDiff,  0,      np.inf]),
                                maxfev=2000)
        a, b, c, d, e = popt[0], popt[1], popt[2], popt[3], popt[4]

        squaredDiffs = np.square(yDecay - ((a * np.exp(b * (xDecay/samplingFreq))) + (c * np.exp(d * (xDecay/samplingFreq))) + e))
        squaredDiffsFromMean = np.square(yDecay - np.mean(yDecay, dtype=np.float64))
        rSquared = 1 - (np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean))
        print(np.sum(squaredDiffs), np.sum(squaredDiffsFromMean))
        print("Trace %i: %f"%(i+1, rSquared))
        xDecFinal = np.linspace(np.min(xDecay), np.max(xDecay), len(filteredSignalArray[i]))
        yDecFinal = (a * np.exp(b * (xDecFinal/samplingFreq))) + (c * np.exp(d * (xDecFinal/samplingFreq))) + e
        yDecFinalArray[i] = yDecFinal
        rSquaredArray[i] = rSquared
    return yDecFinalArray, rSquaredArray

def singleTraceDecayPlot(filteredSignal, fittedDecay, chosenTrace, ratSide):
    isosbesticY = filteredSignal[chosenTrace]
    isosbesticX = np.arange(0, len(filteredSignal[chosenTrace]))
    decayFitFig, ax1 = plt.subplots()
    # ax1.plot(isosbesticX, isosbesticY)
    ax1.plot(isosbesticX, fittedDecay[chosenTrace])    
    decayFitFig.suptitle(ratSide)
    #plt.xticks(ticks=seconds, labels=range(len(seconds)))
    plt.minorticks_on()
    plt.show()

def unbleachSignal(filteredSignal, decayFactor):
    unbleachedArray = np.zeros((len(filteredSignal), len(filteredSignal[0])))
    convert1d = 0
    for i, y in enumerate(filteredSignal):
        try:
            unbleachedArray[i] = y - decayFactor[convert1d]
            convert1d += 1
        except IndexError:
            pass
    return unbleachedArray

# Reads a CSV file containing all previous data for decay values and averages them, creating an array of average decay
def averageCSV(ratNameLeft, ratNameRight):
    meanDecay405Left, meanDecay470Left, meanDecay405Right, meanDecay470Right = np.arange(0, 47750), np.arange(0, 47750), np.arange(0, 47750), np.arange(0, 47750)
    for x in ['Left405Decay.csv', 'Right405Decay.csv', 'Left470Decay.csv', 'Right470Decay.csv']:
        meanReader = pd.read_csv(x)
        meanReader = meanReader.set_index("Experiment", drop=True).rename_axis(None)
        meanIndex = pd.Index(meanReader.index)
        decayArray = pd.DataFrame()
        for row in meanIndex:
            if "Rat %s"%ratNameLeft in row:
                decayArray = pd.concat([decayArray, meanReader.loc[row]], axis=1)
            elif "Rat %s"%ratNameRight in row:
                decayArray = pd.concat([decayArray, meanReader.loc[row]], axis=1)
            else:
                pass
        # deleteHere = np.argwhere(decayArray==np.NaN)
        # decayArray = np.delete(decayArray, deleteHere)
        # decayArray[np.all(np.isnan(decayArray), axis=0)] = 0
        decayArray = decayArray.transpose()
        match x:
            case 'Left405Decay.csv':
                meanDecay405Left = decayArray.mean(axis=0)
            case 'Right405Decay.csv':
                meanDecay405Right = decayArray.mean(axis=0)
            case 'Left470Decay.csv':
                meanDecay470Left = decayArray.mean(axis=0)
            case 'Right470Decay.csv':
                meanDecay470Right = decayArray.mean(axis=0)
            case None:
                raise ValueError("Please Generate a CSV file from a vehicle control first")
            case _:
                raise ValueError("Invalid CSV")
    return meanDecay405Left, meanDecay470Left, meanDecay405Right, meanDecay470Right

# Saves cleaned trace file as an ABF file for later viewing in ClampFit 10 "Processed Data", "%s Rat %s Processed Data.abf"%(self.abfDate.strftime("%Y-%m-%d")
def tExport(processedTrace, ratName, experimentDate):
    abfWriter.writeABF1(sweepData= processedTrace, filename= os.path.join(os.getcwd(), "Processed Data", "%s Rat %s Processed Data.abf"%(experimentDate, ratName)), units="V", sampleRateHz= 3333.33)

def tExportNew(processedTrace, ratName, experimentDate):
    abfWriter.writeABF1(sweepData= processedTrace, filename= os.path.join(os.getcwd(), "Processed Data", "%s Rat %s WIP Processed Data.abf"%(experimentDate, ratName)), units="V", sampleRateHz= 3333.33)