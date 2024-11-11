import scipy.signal as sci
import scipy.integrate as inte
import scipy.optimize as opt
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
        peaks, _ = sci.find_peaks(t, prominence= 0.05, wlen=20000)
        peakArray[i] = np.size(peaks)
    longestPeak = np.max(peakArray)
    return int(longestPeak)

def tauFit(tau, a, b, c):
    return a*np.exp(-b/tau) + c



# def riseTauFit(tau, a, b):
#     return a*np.exp(-b * tau)

# Finds peaks in a signal and provides their widths, amplitudes, avg. frequencies, and areas across an entire .abf file
def wholeTracePeaks(processedSignalArray, mainFile):
    longPeak = peakMax(processedSignalArray)
    abf = pyabf.ABF(mainFile)
    traceLen = len(abf.sweepList)
    samplingFreqMSec = abf.dataPointsPerMs + (1/3)
    samplingFreqSec = samplingFreqMSec * 1000
    peaksDict, finalDict = {}, {}
    peaksArray, adjustedArea = np.zeros((traceLen, longPeak)), np.zeros((traceLen, longPeak))
    widthArray = np.zeros((traceLen, 4, longPeak))
    overlapPeaks = [[0]*longPeak for _ in abf.sweepList]
    for index, traces in enumerate(processedSignalArray): # Finds peaks in a signal
        peaks, peaksDict[index] = sci.find_peaks(traces, prominence= 0.05, width=0, wlen=20000, rel_height= 0.5)
        bottomWidth = sci.peak_widths(traces, peaks, rel_height=1, 
                                      prominence_data=(peaksDict[index]['prominences'], peaksDict[index]["left_bases"], 
                                                       peaksDict[index]["right_bases"]), wlen=20000)
        peaks = np.pad(peaks, pad_width= (0, longPeak - len(peaks)), mode= 'constant', constant_values= 0)
        bottomWidth = np.array(bottomWidth)
        for i, param in enumerate(bottomWidth):
            param = np.pad(param, pad_width= (0, longPeak - len(param)), mode= 'constant', constant_values= 0)
            widthArray[index][i] = param
        peaksArray[index] = peaks
        for i in peaksDict[index]:
           paddedEntry = np.pad(peaksDict[index][i], pad_width= (0, longPeak - len(peaksDict[index][i])), mode= 'constant', constant_values= 0)
           peaksDict[index][i] = paddedEntry
    for k, checkPeak in np.ndenumerate(peaksArray): # Finds peaks that have overlap with one another
            overlapPeaks[k[0]][k[1]] = np.array([y for y in peaksArray[k[0]] if ((widthArray[k[0]][2][k[1]] < y < widthArray[k[0]][3][k[1]]) and y != checkPeak)])
    for z, x in enumerate(peaksArray):
        degreeNPeaks = {}
        peakTable = pd.DataFrame(columns= ['Event_Num', 'Peak_Index', 'Peak_Time_Sec', 
                                           'Event_Window_Start', 'Event_Window_End', 
                                           'Amplitude', 
                                           'Off_Time_ms', 'Width_at50_ms',
                                           'Frequency', 
                                           'Avg_Area', 'Total_Area',
                                           'Decay_Tau_exp', 'Rise_Tau_exp'])
        peakTable.Event_Num = [x + 1 for x, _ in enumerate(x)]
        peakTable.Peak_Index = x
        peakTable.Peak_Time_Sec = ((x/samplingFreqSec) + (z * 30)).round(2)
        peakTable.Event_Window_Start = peaksDict[z]['left_ips'].round(2)
        peakTable.Event_Window_End = peaksDict[z]['right_ips'].round(2)
        peakTable.Amplitude = peaksDict[z]["prominences"].round(3)
        peakTable.Off_Time_ms = ((peaksDict[z]['right_bases'] - x)/(samplingFreqMSec)).round(2)
        peakTable.Width_at50_ms = (peaksDict[z]['widths']/(samplingFreqMSec)).round(2)
        
        for i, peakOfDegree in enumerate(x): # Determines the "degree" of a peak; how many peaks overlap with it
            if len(processedSignalArray[z][int(widthArray[z][2][i]):int(widthArray[z][3][i])]) == 0:
                 continue
            degreeNPeaks[peakOfDegree] = len(overlapPeaks[z][i])
        sortedDegreeNPeaks = dict(sorted(degreeNPeaks.items(), key=lambda item: item[1]))
        for _, u in enumerate(sortedDegreeNPeaks): # Generates event area
            i = np.where(x == u)
            if len(processedSignalArray[z][int(widthArray[z][2][i[0][0]]):int(widthArray[z][3][i[0][0]])]) == 0:
                 continue
            peaksInArea = overlapPeaks[z][i[0][0]]
            peakArea = inte.simpson(y=processedSignalArray[z][int(widthArray[z][2][i[0][0]]):int(widthArray[z][3][i[0][0]])], 
                                            x=np.arange(int(widthArray[z][2][i[0][0]]), int(widthArray[z][3][i[0][0]])))
            match degreeNPeaks[u]:
                case 0:
                    adjustedArea[z][i[0][0]] = peakArea.round(2)
                case _:
                    for l in peaksInArea:
                        j = np.where(x == l)
                        peakArea -= adjustedArea[z][j[0][0]]
                    adjustedArea[z][i[0][0]] = peakArea.round(2)
                    
        peakTable.Avg_Area = pd.Series(adjustedArea[z])
        match np.size(x): # Determines frequency and total area of a sweep
            case 0:
                peakTable.Frequency = 0
                peakTable.Total_Area = 0
            case _:
                peakTable.Frequency.iat[0] = round(np.count_nonzero(x)/((len(processedSignalArray[z]) + 2250)/samplingFreqSec), 2) #Peaks/second (15 second trace)
                peakTable.Total_Area.iat[0] = sum(adjustedArea[z])
        peakTable.drop(peakTable[peakTable.Peak_Index == 0].index, inplace= True)
        
        for i, p in enumerate(x):
            if len(processedSignalArray[z][int(widthArray[z][2][i]):int(widthArray[z][3][i])]) == 0:
                 continue
            # print("Hi, I'm trying!")
            # Event Decay
            decayTau = processedSignalArray[z][int(p):peaksDict[z]['right_bases'][i]]
            decayWidth = int(len(decayTau))
            decArray = list(range(0, decayWidth))
            aInitial = 200 #Temp value
            bInitial = 0.5 #Temp value

            popt, pcov = opt.curve_fit(tauFit, decArray, decayTau, p0=(aInitial, bInitial), bounds=opt.Bounds(lb=[0.0, 0.0], ub=[np.inf, np.inf]), maxfev=2000)
            a, b = popt[0], popt[1]
            peakTable.Decay_Tau_exp = abs((1/b)/(samplingFreqMSec))

            # Event Rise
            riseTau = processedSignalArray[z][peaksDict[z]['left_bases'][i]:int(p)]
            riseWidth = int(len(riseTau))
            riseArray = list(range(0, riseWidth))
            aInitial = 200 #Temp value
            bInitial = 0.5 #Temp value

            popt, pcov = opt.curve_fit(tauFit, riseArray, riseTau, p0=(aInitial, bInitial), bounds=opt.Bounds(lb=[0.0, 0.0], ub=[np.inf, np.inf]), maxfev=2000)
            a, b = popt[0], popt[1]
            peakTable.Rise_Tau_exp = abs((1/b)/(samplingFreqMSec))
            # print("Hi, I tried!")


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
    seconds = np.arange(0, 47750, samplingFreqSec)
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, width=0, wlen= 20000, rel_height= 0.5)
    widthBottom = sci.peak_widths(processedSignalArray, peaks, rel_height=1, prominence_data=(peaksDict['prominences'], 
                                                                                             peaksDict["left_bases"], peaksDict["right_bases"]), wlen=20000)
    widthHalf = sci.peak_widths(processedSignalArray, peaks, rel_height=0.5, prominence_data=(peaksDict['prominences'], 
                                                                                             peaksDict["left_bases"], peaksDict["right_bases"]), wlen=20000)
    width10 = sci.peak_widths(processedSignalArray, peaks, rel_height=0.1, prominence_data=(peaksDict['prominences'], 
                                                                                             peaksDict["left_bases"], peaksDict["right_bases"]), wlen=20000)
    width90 = sci.peak_widths(processedSignalArray, peaks, rel_height=0.9, prominence_data=(peaksDict['prominences'], 
                                                                                             peaksDict["left_bases"], peaksDict["right_bases"]), wlen=20000)
    fig = plt.figure()
    peakFig = fig.add_subplot()
    peakFig.plot(processedSignalArray)
    peakFig.plot(peaks, processedSignalArray[peaks], "r.")
    for i, x in enumerate(peaks):
        peakFig.plot(int(width90[3][i]), processedSignalArray[int(width90[3][i])], 'g.')
        peakFig.annotate("Trough for Peak %i"%(i+1), xycoords= 'data', size= 8, horizontalalignment= 'center',
                     xytext = (int(width90[3][i]), processedSignalArray[int(width90[3][i])] - 0.3), 
                     xy = (int(width90[3][i]), processedSignalArray[int(width90[3][i])] - 0.01),
                     arrowprops=dict(facecolor= 'black', width= 1, headwidth= 5, headlength= 5))
        peakFig.annotate("Peak %i"%(i+1), xycoords= 'data', size= 8, horizontalalignment= 'center',
                         xytext= (x, processedSignalArray[x] + 0.01),
                         xy = (x, processedSignalArray[x]))
        peakFig.fill_between(np.arange(int(widthBottom[2][i]), int(widthBottom[3][i])), processedSignalArray[int(widthBottom[2][i]):int(widthBottom[3][i])], 
                             widthBottom[1][i], color="C1", alpha=0.3)
        p0 = (1, -10, processedSignalArray[x])
        # Event Decay
        decayTau = processedSignalArray[int(x):int(width90[3][i])]
        if len(decayTau) == 0:
            continue
        else:
            decArray = list(range(int(x), int(width90[3][i])))
            popt, pcov = opt.curve_fit(tauFit, decArray, decayTau, p0=p0, maxfev=35000)
            print(popt)
            x_dec = np.linspace(np.min(decArray), np.max(decArray), 1000)
            y_dec = tauFit(x_dec, *popt)
            peakFig.plot(x_dec, y_dec, color="C9")
        p0 = (1, 0.3, int(width90[2][i]))
        # Event Rise
        riseTau = processedSignalArray[int(width90[2][i]):x]
        if len(riseTau) == 0:
            continue
        else:
            riseArray = list(range(int(width90[2][i]), x))
            popt, pcov = opt.curve_fit(tauFit, riseArray, riseTau, p0=p0, maxfev=35000)
            x_rise = np.linspace(np.min(riseArray), np.max(riseArray), 1000)
            y_rise = tauFit(x_rise, *popt)
            peakFig.plot(x_rise, y_rise, color="C8")

    
    peakFig.hlines(*widthHalf[1:], color="C6")
    peakFig.hlines(*widthBottom[1:], color="C7")
    peakFig.vlines(x=peaks, ymin=processedSignalArray[peaks] - peaksDict["prominences"], ymax=processedSignalArray[peaks], color="C5")
    peakFig.set_title(ratSide)
    plt.axis([0, 47750, 0, 2])
    plt.xticks(ticks=seconds, labels=range(len(seconds)))
    plt.xlabel("Time (s)")
    plt.ylabel("Fluorescence (AU)")
    plt.minorticks_on()
    plt.show()
