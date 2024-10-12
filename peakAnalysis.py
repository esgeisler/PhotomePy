import scipy.signal as sci
import scipy.integrate as inte
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyabf

# Returns the number of peaks in the trace with the most peaks
def peakMax(processedSignalArray):
    peakList = []
    for traces in range(len(processedSignalArray)):
        peaks, peaksDict = sci.find_peaks(processedSignalArray[traces], prominence= 0.05, width=0, wlen=20000, rel_height= 0.5)
        peakList.append(peaks)
    longestPeak = max((len(x)) for x in peakList)
    return longestPeak

# Finds peaks in a signal and provides their widths, amplitudes, avg. frequencies, and areas across an entire .abf file
def wholeTracePeaks(processedSignalArray, mainFile):
    longPeak = peakMax(processedSignalArray)
    abf = pyabf.ABF(mainFile)
    samplingFreqMSec = abf.dataPointsPerMs + (1/3)
    samplingFreqSec = samplingFreqMSec * 1000
    finalDict = {}
    peaksArray = np.zeros((len(abf.sweepList), longPeak))
    peaksDict = {}
    for traces in range(len(processedSignalArray)):
        peaks, peaksDict[traces] = sci.find_peaks(processedSignalArray[traces], prominence= 0.05, width=0, wlen=20000, rel_height= 0.5)
        peaks = np.pad(peaks, pad_width= (0, longPeak - len(peaks)), mode= 'constant', constant_values= 0)
        peaksArray[traces] = peaks
        for i in peaksDict[traces]:
           paddedEntry = np.pad(peaksDict[traces][i], pad_width= (0, longPeak - len(peaksDict[traces][i])), mode= 'constant', constant_values= 0)
           peaksDict[traces][i] = paddedEntry
    z = 0
    for x in peaksArray:
        peakTable = pd.DataFrame(columns= ['Event_Num', 'Peak_Index', 
                                        'Peak_Time_Sec', 'Event_Window_Start', 
                                        'Event_Window_End', 'Amplitude', 'Off_Time_ms',
                                        'Width_at50_ms','Frequency', 'Area'])
        peakTable.Event_Num = [x + 1 for x in range(len(x))]
        peakTable.Peak_Index = x
        peakTable.Peak_Time_Sec = ((x/samplingFreqSec) + (z * 30)).round(2)
        peakTable.Event_Window_Start = peaksDict[z]['left_ips'].round(2)
        peakTable.Event_Window_End = peaksDict[z]['right_ips'].round(2)
        peakTable.Amplitude = peaksDict[z]["prominences"].round(3)
        peakTable.Off_Time_ms = ((peaksDict[z]['right_bases'] - x)/(samplingFreqMSec)).round(2)
        peakTable.Width_at50_ms = (peaksDict[z]['widths']/(samplingFreqMSec)).round(2)
        if x.size == 0:
            peakTable.Frequency = 0
        else:
            peakTable.Frequency.iat[0] = round(np.count_nonzero(x)/(len(processedSignalArray)/samplingFreqMSec), 2) #Peaks/second (15 second trace)
        areaList = []
        for i in range(len(x)):
            if len(processedSignalArray[z][int(peaksDict[z]['left_bases'][i]):int(peaksDict[z]['right_bases'][i])]) == 0:
                continue
            peakArea = inte.simpson(y=processedSignalArray[z][int(peaksDict[z]['left_bases'][i]):int(peaksDict[z]['right_bases'][i])], 
                                    x=range(int(peaksDict[z]['left_bases'][i]), int(peaksDict[z]['right_bases'][i])))
            areaList.append(peakArea.round(2))
        areaList = pd.Series(areaList)
        peakTable.Area = areaList

        peakTable.drop(peakTable[peakTable.Peak_Index == 0].index, inplace= True)
        finalDict[z] = peakTable
        z += 1
    return finalDict

#TODO fix FutureWarning caused by pre- and postOverview being empty by default.
def traceProcessor(processedSignal):
    injectionDF = {}
    x = 0
    traceList = []
    for traces in processedSignal.values():
        injectionDF[x] = traces
        traceList.append(traces)
        x += 1
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
                                       'Event_Window_End','Amplitude',
                                        'Frequency', 'Area'])
    peakTable.event = [x for x in range(len(peaks))]
    peakTable.Peak_Index = peaks
    peakTable.PeakTimeSec = peaks/samplingFreqSec
    peakTable.Event_Window_Start = peaksDict['left_ips']
    peakTable.Event_Window_End = peaksDict['right_ips']
    peakTable.Frequency = np.count_nonzero(peaks)/(len(processedSignalArray)/samplingFreqSec)
    areaList = []
    for i in range(len(peaks)):
        if len(processedSignalArray[int(peaksDict['left_bases'][i]):int(peaksDict['right_bases'][i])]) == 0:
            continue
        peakArea = inte.simpson(y=processedSignalArray[int(peaksDict['left_bases'][i]):int(peaksDict['right_bases'][i])], 
                                x=range(int(peaksDict['left_bases'][i]), int(peaksDict['right_bases'][i])))
        areaList.append(peakArea.round(2))
    areaList = pd.Series(areaList)
    peakTable.Area = areaList

    widthBottom = sci.peak_widths(processedSignalArray, peaks, rel_height=1)
    widthHalf = sci.peak_widths(processedSignalArray, peaks, rel_height=0.5)
    
    fig = plt.figure()
    peakFig = fig.add_subplot()
    peakFig.plot(processedSignalArray)
    peakFig.plot(peaks, processedSignalArray[peaks], "r.")
    peakFig.plot(peaksDict['right_bases'], processedSignalArray[peaksDict['right_bases']], 'g.')
    for x in range(len(processedSignalArray[peaks])):
        peakFig.annotate("Trough", xycoords= 'data', size= 10, horizontalalignment= 'center',
                     xytext = (peaksDict['right_bases'][x], processedSignalArray[peaks][x] - 0.3), 
                     xy = (peaksDict['right_bases'][x], processedSignalArray[peaksDict['right_bases'][x]] - 0.01),
                     arrowprops=dict(facecolor= 'black', width= 1, headwidth= 5, headlength= 5))
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