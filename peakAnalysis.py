import AutoCleaner as acl
import AverageTraces as avg
import scipy.signal as sci
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyabf

#Retrieves the peaks of a single trace and returns a list containing the peaks in a ndarray and their properties in a dictionary
def peakGetter(processedSignalArray):
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, height=0, width=0)
    return [peaks, peaksDict]
    
def wholeTracePeaks(processedSignalArray):
    peaksList = []
    peaksDict = {}
    x = 0
    for traces in processedSignalArray:
        peaks, peaksDict[x] = sci.find_peaks(processedSignalArray[x][850:-1250], prominence= 0.05, height=0)
        x += 1
        peaksList.append(peaks)
    return [peaksList, peaksDict]

#Retrieves the peaks of a signal and plots them on a graph of the chosen trace
#TODO pre and post-trigger window to remove big goofy start and end
def peakDisplay(processedSignalArray, mainFile):
    abf = pyabf.ABF(mainFile)
    samplingFreq = int(abf.dataPointsPerMs * 1000)
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, height=0, width=0)
    peakTable = pd.DataFrame(columns= ['event', 'Peak_Index', 
                                       'PeakTimeSec', 'Event_Window_Start', 
                                       'Event_Window_End', 'Amplitude', 
                                       'WidthMS','Frequency'])
    peakTable.event = [x for x in range(len(peaks))]
    peakTable.Peak_Index = peaks
    peakTable.PeakTimeSec = peaks/samplingFreq
    peakTable.Event_Window_Start = peaksDict['left_ips']
    peakTable.Event_Window_End = peaksDict['right_ips']
    peakTable.Amplitude = peaksDict['peak_heights']
    peakTable.WidthMS = peaksDict['widths']/(samplingFreq/1000)
    peakTable.Frequency = len(peaks)/15 #Peaks/second (15 second trace)
    
    fig = plt.figure()
    lad = fig.add_subplot()
    lad.plot(processedSignalArray)
    lad.plot(peaks, processedSignalArray[peaks], "r.")
    finalTable = plt.table(cellText= peakTable.values, colLabels= peakTable.keys())
    peakTable.round(3)
    finalTable.set_fontsize(40)
    finalTable.scale(1.5, 1.5)
    lad.add_table(finalTable)
    plt.axis([0,50000, -2, 2])
    plt.show()
    

def peakDecay(processedSignalArray):
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, width=0)
    peakdecayList = []
    x = 0
    for peak in peaks:
        peakdecayList.append(peak - peaksDict["right_ips"][x])
        x += 1
    return peakdecayList

def peakAmplitude(processedSignalArray):
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, height=0)
    return peaksDict["peak_heights"]

def peakFreq(processedSignalArray):
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, height=0)
    return len(peaks)

def peakWidth(processedSignalArray):
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05)
    width, widthHeight, leftRight = sci.peak_widths(processedSignalArray, peaks)
    return width

# Plots the original, subtracted, ratio-ed, and processed trace of choice
    # abf.setSweep(sweepNumber= userTrace, channel= userChannel)
    # plt.plot(abf.sweepX[1000:-2000], abf.sweepY[1000:-2000], color="b", label="Original")
    # plt.plot(abf.sweepX[1000:-2000], subtractRight[0][userTrace][1000:-2000], color="r", label= "Subtracted")
    # plt.plot(abf.sweepX[1000:-2000], ratioSignalRight[userTrace][1000:-2000], color="y", label= "Ratio-ed")
    # plt.plot(abf.sweepX[1000:-2000], finalSignalRight[userTrace][1000:-2000], color="g", label= "Processed")
    # plt.axis([0,15,-1,6])
    # plt.legend()
    # plt.show()