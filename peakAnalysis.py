import scipy.signal as sci
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyabf
import sys
import tracePeaks as trace
import yaml
np.set_printoptions(threshold=sys.maxsize)

with open("config.yaml") as c:
    userConfig = yaml.safe_load(c)
    peakMethod = userConfig["EVENT_HANDLING"]["peak_id_method"]
    peakThreshold = userConfig["EVENT_HANDLING"]["peak_detection_threshold"]
    peakWindow = userConfig["EVENT_HANDLING"]["peak_window"]
    peakTop = userConfig["EVENT_HANDLING"]["peak_top"]
    peakBase = userConfig["EVENT_HANDLING"]["peak_bottom"]

def secondsCalculator(filename):
    abf = pyabf.ABF(filename)
    samplingFreqMSec = abf.dataPointsPerMs + (1/3)
    samplingFreqSec = samplingFreqMSec * 1000
    seconds = np.arange(0, 47750, samplingFreqSec)
    return seconds, samplingFreqMSec, samplingFreqSec

# Returns the number of peaks in the trace with the most peaks
def peakMax(processedSignalArray):
    peakLenArray = np.zeros(np.size(processedSignalArray, 0))
    for i, t in enumerate(processedSignalArray):
        peaks, _ = sci.find_peaks(t, prominence= peakThreshold, wlen=peakWindow)
        peakLenArray[i] = np.size(peaks)
    longestPeak = np.max(peakLenArray)
    return int(longestPeak), peakLenArray

# Finds peaks in a signal and provides their widths, amplitudes, avg. frequencies, and areas across an entire .abf file
def wholeTracePeaks(processedSignalArray, mainFile):
    finalDict = {}
    for z, _ in enumerate(processedSignalArray):
        peaks = trace.TracePeaks(processedSignalArray, mainFile, z)
        peaks.peakFinder(peaks.fullTraceArray)
        peaks.overlapCheck(peaks.peaks)
        peaks.overlapAmplitude()
        peaks.freqSet()
        peaks.areaSet()
        peaks.riseSet()
        peaks.decaySet()
        peakTable = pd.DataFrame(columns= ['Event_Num', 'Peak_Index', 'Peak_Time_Sec', 
                                           'Event_Window_Start', 'Event_Window_End', 
                                           'Amplitude', 'Abs_Amplitude',
                                           'Off_Time_ms', 'Width_at50_ms',
                                           'Frequency', 
                                           'Avg_Area', 'Total_Area',
                                           'Rise_Tau', 'Decay_Tau'])
        peakTable.Event_Num = peaks.peakNum
        peakTable.Peak_Index = peaks.peakIndex
        peakTable.Peak_Time_Sec = peaks.peakLocSec
        peakTable.Event_Window_Start = peaks.leftBounds
        peakTable.Event_Window_End = peaks.rightBounds
        peakTable.Amplitude = peaks.amplitude
        peakTable.Abs_Amplitude = peaks.absoluteAmp
        peakTable.Off_Time_ms = peaks.rightTail
        peakTable.Width_at50_ms = peaks.width
        peakTable.Avg_Area = peaks.meanArea
        peakTable.Rise_Tau = pd.Series(peaks.riseRate)
        peakTable.Decay_Tau = pd.Series(peaks.decayRate)
        if np.size(peaks.peaks) == 0:
            peakTable.Frequency = 0
            peakTable.Total_Area = 0
        else:
            peakTable.Frequency.iat[0] = peaks.frequency
            peakTable.Total_Area.iat[0] = peaks.totArea
        # peakTable.drop(peakTable[peakTable.Peak_Index == 0].index, inplace= True)
        print("Trace #%i done!"%(z+1))
        finalDict[z] = peakTable
    return finalDict

# Appends all of the trace data together into one pandas DataFrame.
def traceProcessor(processedSignal):
    injectionDF = {}
    traceList = []
    for x, traces in enumerate(processedSignal.values()):
        injectionDF[x] = traces
        traceList.append(traces)
    overview = pd.concat(traceList)
    return injectionDF, overview

#Retrieves the peaks of a signal and their properties, then plots them on a graph of the chosen trace
def peakDisplay(processedSignalArray, mainFile, ratSide, currentTrace):
    peaks = trace.TracePeaks(processedSignalArray, mainFile, currentTrace)
    peaks.peakFinder(peaks.fullTraceArray)
    peaks.overlapCheck(peaks.peaks)
    peaks.overlapAmplitude()
    peaks.riseSet()
    peaks.decaySet()
    peakFig, ax = plt.subplots()
    ax.plot(peaks.fullTraceArray)
    ax.plot(peaks.peaks, peaks.fullTraceArray[peaks.peaks], "r.")
    for i, x in enumerate(peaks.peaks):
        ax.plot(int(peaks.trace90Widths[1][i]), peaks.fullTraceArray[int(peaks.trace90Widths[1][i])], 'g.')
        ax.annotate("Trough (%i)"%(i+1), xycoords= 'data', size= 8, horizontalalignment= 'center',
                     xytext = (int(peaks.trace90Widths[1][i]), peaks.fullTraceArray[int(peaks.trace90Widths[1][i])] - 0.3), 
                     xy = (int(peaks.trace90Widths[1][i]), peaks.fullTraceArray[int(peaks.trace90Widths[1][i])] - 0.01),
                     arrowprops=dict(facecolor= 'black', width= 1, headwidth= 5, headlength= 5))
        ax.annotate("p%i"%(i+1), xycoords= 'data', size= 8, horizontalalignment= 'center',
                         xytext= (x, peaks.fullTraceArray[x] + 0.01),
                         xy = (x, peaks.fullTraceArray[x]))
        ax.fill_between(x=np.arange(int(peaks.trace90Widths[0][i]), int(peaks.trace90Widths[1][i])), 
                        y1=peaks.fullTraceArray[int(peaks.trace90Widths[0][i]):int(peaks.trace90Widths[1][i])], 
                        y2=peaks.fullTraceArray[peaks.peaks][i] - peaks.amplitude[i],
                        color="C1", alpha=0.3)
        if np.all(peaks.risePlot[i, 0]):
            ax.plot(peaks.risePlot[i, 0, :], peaks.risePlot[i, 1, :], color="C8")
        else:
            continue
        if np.all(peaks.risePlot[i, 0]):
            ax.plot(peaks.decayPlot[i, 0, :], peaks.decayPlot[i, 1, :], color="C9")
        else:
            continue

    ax.hlines(y=peaks.traceDict["width_heights"], xmin=peaks.traceDict["left_ips"], xmax=peaks.traceDict["right_ips"], color="C6")
    ax.vlines(x=peaks.peaks, ymin=peaks.fullTraceArray[peaks.peaks] - peaks.amplitude, ymax=peaks.fullTraceArray[peaks.peaks], 
              color="C5")
    ax.vlines(x=peaks.peaks, ymin=peaks.fullTraceArray[peaks.peaks] - peaks.absoluteAmp, ymax=peaks.fullTraceArray[peaks.peaks], 
              linestyles="dotted", color="C4")
    peakFig.suptitle("Individual Trace")
    peakFig.set_size_inches(9, 4.5)
    ax.axis([0, 47750, 0, 2])
    ax.tick_params(axis="both", labelsize="large")
    ax.set_xticks(ticks=peaks.seconds, labels=range(len(peaks.seconds)))
    ax.set_xlabel("Time (s)", fontsize="x-large")
    ax.set_ylabel("Fluorescence (AU)", fontsize="x-large")
    peakFig.tight_layout()
    plt.minorticks_on()
    plt.show()
