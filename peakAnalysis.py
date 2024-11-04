import scipy.signal as sci
import scipy.integrate as inte
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyabf
import sys
np.set_printoptions(threshold=sys.maxsize)

# Returns the number of peaks in the trace with the most peaks
def peakMax(processedSignalArray):
    peakArray = np.zeros(np.size(processedSignalArray, 0))
    for i, t in enumerate(processedSignalArray):
        peaks, _ = sci.find_peaks(t, prominence= 0.05, width=0, wlen=20000, rel_height= 0.5)
        peakArray[i] = np.size(peaks)
    longestPeak = np.max(peakArray)
    return int(longestPeak)

# Finds peaks in a signal and provides their widths, amplitudes, avg. frequencies, and areas across an entire .abf file
def wholeTracePeaks(processedSignalArray, mainFile):
    longPeak = peakMax(processedSignalArray)
    abf = pyabf.ABF(mainFile)
    samplingFreqMSec = abf.dataPointsPerMs + (1/3)
    samplingFreqSec = samplingFreqMSec * 1000
    finalDict = {}
    peaksArray = np.zeros((len(abf.sweepList), longPeak))
    peaksDict = {}
    widthArray = np.zeros((len(abf.sweepList), 4, longPeak))
    overlapPeaks = [[0]*longPeak for i in abf.sweepList]
    for index, traces in enumerate(processedSignalArray):
        peaks, peaksDict[index] = sci.find_peaks(traces, prominence= 0.05, width=0, wlen=20000, rel_height= 0.5)
        bottomWidth = sci.peak_widths(traces, peaks, rel_height=1, 
                                      prominence_data=(peaksDict[index]['prominences'], 
                                        peaksDict[index]["left_bases"], peaksDict[index]["right_bases"]), 
                                        wlen=20000)
        peaks = np.pad(peaks, pad_width= (0, longPeak - len(peaks)), mode= 'constant', constant_values= 0)
        bottomWidth = np.array(bottomWidth)
        for i, param in enumerate(bottomWidth):
            param = np.pad(param, pad_width= (0, longPeak - len(param)), mode= 'constant', constant_values= 0)
            widthArray[index][i] = param
        peaksArray[index] = peaks
        for i in peaksDict[index]:
           paddedEntry = np.pad(peaksDict[index][i], pad_width= (0, longPeak - len(peaksDict[index][i])), mode= 'constant', constant_values= 0)
           peaksDict[index][i] = paddedEntry
    for k, checkPeak in np.ndenumerate(peaksArray):
            overlapPeaks[k[0]][k[1]] = np.array([y for y in peaksArray[k[0]] if ((widthArray[k[0]][2][k[1]] < y < widthArray[k[0]][3][k[1]]) and y != checkPeak)])
    for z, x in enumerate(peaksArray):
        peakTable = pd.DataFrame(columns= ['Event_Num', 'Peak_Index', 
                                        'Peak_Time_Sec', 'Event_Window_Start', 
                                        'Event_Window_End', 'Amplitude', 'Off_Time_ms',
                                        'Width_at50_ms','Frequency', 'Area'])
        peakTable.Event_Num = [x + 1 for x, _ in enumerate(x)]
        peakTable.Peak_Index = x
        peakTable.Peak_Time_Sec = ((x/samplingFreqSec) + (z * 30)).round(2)
        peakTable.Event_Window_Start = peaksDict[z]['left_ips'].round(2)
        peakTable.Event_Window_End = peaksDict[z]['right_ips'].round(2)
        peakTable.Amplitude = peaksDict[z]["prominences"].round(3)
        peakTable.Off_Time_ms = ((peaksDict[z]['right_bases'] - x)/(samplingFreqMSec)).round(2)
        peakTable.Width_at50_ms = (peaksDict[z]['widths']/(samplingFreqMSec)).round(2)
        if np.size(x) == 0:
            peakTable.Frequency = 0
        else:
            peakTable.Frequency.iat[0] = round(np.count_nonzero(x)/((len(processedSignalArray[z]) + 2250)/samplingFreqSec), 2) #Peaks/second (15 second trace)
        areaList = []
        for i, peakNum in enumerate(x):
            if len(processedSignalArray[z][int(widthArray[z][2][i]):int(widthArray[z][3][i])]) == 0:
                continue
            peakArea = inte.simpson(y=processedSignalArray[z][int(widthArray[z][2][i]):int(widthArray[z][3][i])], 
                                    x=np.arange(int(widthArray[z][2][i]), int(widthArray[z][3][i])))
            peaksInArea = overlapPeaks[z][i]
            if len(peaksInArea) > 0:
                for areaPeakSub in peaksInArea:
                    j = np.where(x == areaPeakSub)
                    peakArea -= inte.simpson(y=processedSignalArray[z][int(widthArray[z][2][j]):int(widthArray[z][3][j])], 
                                    x=np.arange(int(widthArray[z][2][j]), int(widthArray[z][3][j])))
                areaList.append(peakArea.round(2))
            else:     
                areaList.append(peakArea.round(2))
            
        peakTable.Area = pd.Series(areaList)

        peakTable.drop(peakTable[peakTable.Peak_Index == 0].index, inplace= True)
        finalDict[z] = peakTable
    return finalDict

#TODO fix FutureWarning caused by pre- and postOverview being empty by default.
def traceProcessor(processedSignal):
    injectionDF = {}
    traceList = []
    for x, traces in enumerate(processedSignal.values()):
        injectionDF[x] = traces
        traceList.append(traces)
    overview = pd.concat(traceList)
    return injectionDF, overview

#Retrieves the peaks of a signal and their properties, then plots them on a graph of the chosen trace
def peakDisplay(processedSignalArray, mainFile, ratSide):
    abf = pyabf.ABF(mainFile)
    samplingFreqMSec = abf.dataPointsPerMs + (1/3)
    samplingFreqSec = samplingFreqMSec * 1000
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, width=0, wlen= 20000, rel_height= 0.5)
    peakTable = pd.DataFrame(columns= ['event', 'Peak_Index', 
                                       'PeakTimeSec', 'Event_Window_Start', 
                                       'Event_Window_End','Amplitude'])
    peakTable.event = [x + 1 for x, _ in enumerate(peaks)]
    peakTable.Peak_Index = peaks
    peakTable.PeakTimeSec = peaks/samplingFreqSec
    peakTable.Event_Window_Start = peaksDict['left_ips']
    peakTable.Event_Window_End = peaksDict['right_ips']
    

    widthBottom = sci.peak_widths(processedSignalArray, peaks, rel_height=1, prominence_data=(peaksDict['prominences'], 
                                                                                             peaksDict["left_bases"], peaksDict["right_bases"]), wlen=20000)
    widthHalf = sci.peak_widths(processedSignalArray, peaks, rel_height=0.5, prominence_data=(peaksDict['prominences'], 
                                                                                             peaksDict["left_bases"], peaksDict["right_bases"]), wlen=20000)
    
    fig = plt.figure()
    peakFig = fig.add_subplot()
    peakFig.plot(processedSignalArray)
    peakFig.plot(peaks, processedSignalArray[peaks], "r.")
    peakFig.plot(peaksDict['right_bases'], processedSignalArray[peaksDict['right_bases']], 'g.')
    for i, x in enumerate(processedSignalArray[peaks]):
        peakFig.annotate("Trough", xycoords= 'data', size= 10, horizontalalignment= 'center',
                     xytext = (peaksDict['right_bases'][i], x - 0.3), 
                     xy = (peaksDict['right_bases'][i], processedSignalArray[peaksDict['right_bases'][i]] - 0.01),
                     arrowprops=dict(facecolor= 'black', width= 1, headwidth= 5, headlength= 5))
        peakFig.fill_between(np.arange(int(widthBottom[2][i]), int(widthBottom[3][i])), processedSignalArray[int(widthBottom[2][i]):int(widthBottom[3][i])], 
                             widthBottom[1][i], color="C1", alpha=0.3)
    peakFig.hlines(*widthHalf[1:], color="C6")
    peakFig.hlines(*widthBottom[1:], color="C7")
    peakFig.vlines(x=peaks, ymin=processedSignalArray[peaks] - peaksDict["prominences"], ymax=processedSignalArray[peaks], color="C5")
    peakFig.set_title(ratSide)
    plt.axis([0, 47750, 0, 2])
    seconds = np.arange(0, 47750, samplingFreqSec)
    plt.xticks(ticks=seconds, labels=range(len(seconds)))
    plt.xlabel("Time (s)")
    plt.ylabel("Fluorescence (AU)")
    plt.minorticks_on()
    plt.show()