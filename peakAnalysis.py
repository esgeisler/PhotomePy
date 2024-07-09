import scipy.signal as sci
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyabf

#Retrieves the peaks of a single trace and returns a list containing the peaks in a ndarray and their properties in a dictionary
def peakGetter(processedSignalArray):
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, height=0, width=0, wlen= 20000, rel_height= 0.5)
    return [peaks, peaksDict]

# Returns the number of peaks in the trace with the most peaks
def peakMaxxer(processedSignalArray):
    peakList = []
    for traces in range(len(processedSignalArray)):
        peaks, peaksDict = sci.find_peaks(processedSignalArray[traces], prominence= 0.05, height=0, width=0, wlen=20000, rel_height= 0.5)
        peakList.append(peaks)
    longestPeak = max((len(x)) for x in peakList)
    return longestPeak

# Finds peaks, widths, amplitude, frequency in every trace in a signal
def wholeTracePeaks(processedSignalArray, mainFile):
    longPeak = peakMaxxer(processedSignalArray)
    abf = pyabf.ABF(mainFile)
    samplingFreq = int(abf.dataPointsPerMs * 1000)
    finalDict = {}
    peaksArray = np.zeros((len(abf.sweepList), longPeak))
    peaksDict = {}
    for traces in range(len(processedSignalArray)):
        peaks, peaksDict[traces] = sci.find_peaks(processedSignalArray[traces], prominence= 0.05, height=0, width=0, wlen=20000, rel_height= 0.5)
        peaks = np.pad(peaks, pad_width= (0, longPeak - len(peaks)), mode= 'constant', constant_values= 0)
        peaksArray[traces] = peaks
        for i in peaksDict[traces]:
           paddedEntry = np.pad(peaksDict[traces][i], pad_width= (0, longPeak - len(peaksDict[traces][i])), mode= 'constant', constant_values= 0)
           peaksDict[traces][i] = paddedEntry
    z = 0
    for x in peaksArray:
        peakTable = pd.DataFrame(columns= ['Event_Num', 'Peak_Index', 
                                        'Peak_Time_Sec', 'Event_Window_Start', 
                                        'Event_Window_End', 'Amplitude', 'Peak_Decay_ms',
                                        'Width_at50_ms','Frequency'])
        peakTable.Event_Num = [x + 1 for x in range(len(x))]
        peakTable.Peak_Index = x
        peakTable.Peak_Time_Sec = ((x/samplingFreq) + (z * 30)).round(2)
        peakTable.Event_Window_Start = peaksDict[z]['left_ips'].round(2)
        peakTable.Event_Window_End = peaksDict[z]['right_ips'].round(2)
        peakTable.Amplitude = (peaksDict[z]['peak_heights'] - processedSignalArray[z][peaksDict[z]['right_bases']]).round(2)
        peakTable.Peak_Decay_ms = ((peaksDict[z]['right_bases'] - x)/(samplingFreq/1000)).round(2)
        peakTable.Width_at50_ms = (peaksDict[z]['widths']/(samplingFreq/1000)).round(2)
        peakTable.Frequency = round(np.count_nonzero(x)/15, 2) #Peaks/second (15 second trace)

        peakTable.drop(peakTable[peakTable.Peak_Index == 0].index, inplace= True)
        finalDict[z] = peakTable
        z += 1
    return finalDict

#TODO amplitude should be peak to trough height
def traceProcessor(processedSignal, injectionTrace):
    preInjectionDF = {}
    preOverview = pd.DataFrame(columns= ['Event_Num', 'Peak_Index', 
                                        'Peak_Time_Sec', 'Event_Window_Start', 
                                        'Event_Window_End', 'Amplitude', 'Peak_Decay_ms',
                                        'Width_at50_ms','Frequency'])
    postInjectionDF = {}
    postOverview = pd.DataFrame(columns= ['Event_Num', 'Peak_Index', 
                                        'Peak_Time_Sec', 'Event_Window_Start', 
                                        'Event_Window_End', 'Amplitude', 'Peak_Decay_ms',
                                        'Width_at50_ms','Frequency'])
    x = 0
    for traces in processedSignal.values():
        if traces.empty:
            continue
        elif x <= injectionTrace:
            preInjectionDF[x] = traces
            preOverview = pd.concat([preOverview, traces], ignore_index= True)
        elif x > injectionTrace:
            postInjectionDF[x] = traces
            postOverview = pd.concat([postOverview, traces], ignore_index= True)
        x += 1
    
    return preInjectionDF, postInjectionDF, preOverview, postOverview

#Retrieves the peaks of a signal and their properties, then plots them on a graph of the chosen trace
def peakDisplay(processedSignalArray, mainFile, ratSide):
    abf = pyabf.ABF(mainFile)
    samplingFreq = int(abf.dataPointsPerMs * 1000)
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, height=0, width=0, wlen= 20000, rel_height= 0.5)
    peakTable = pd.DataFrame(columns= ['event', 'Peak_Index', 
                                       'PeakTimeSec', 'Event_Window_Start', 
                                       'Event_Window_End', 
                                       'WidthMS'])
    peakTable.event = [x for x in range(len(peaks))]
    peakTable.Peak_Index = peaks
    peakTable.PeakTimeSec = peaks/samplingFreq
    peakTable.Event_Window_Start = peaksDict['left_ips']
    peakTable.Event_Window_End = peaksDict['right_ips']
    peakTable.WidthMS = peaksDict['widths']/(samplingFreq/1000)
    
    fig = plt.figure()
    lad = fig.add_subplot()
    lad.plot(processedSignalArray)
    lad.plot(peaks, processedSignalArray[peaks], "r.")
    lad.plot(peaksDict['right_bases'], processedSignalArray[peaksDict['right_bases']], 'g.')
    for x in range(len(processedSignalArray[peaks])):
        lad.annotate("Trough", xycoords= 'data', size= 10, horizontalalignment= 'center',
                     xytext = (peaksDict['right_bases'][x], processedSignalArray[peaks][x] - 0.3), 
                     xy = (peaksDict['right_bases'][x], processedSignalArray[peaksDict['right_bases'][x]] - 0.01),
                     arrowprops=dict(facecolor= 'black', width= 1, headwidth= 5, headlength= 5)) #, horizontalalignment= 'center', verticalalignment= 'bottom')
    lad.set_title(ratSide)
    plt.axis([0,50000,0.5, 1.5])
    plt.xticks(peaks)
    plt.show()