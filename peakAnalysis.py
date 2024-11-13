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

# Exponential decay equation for plotting the changes in fluorescence after a dopamine event
def tauFit(tau, a, b):
    return a*np.exp(-b * tau)

# Finds peaks in a signal and provides their widths, amplitudes, avg. frequencies, and areas across an entire .abf file
def wholeTracePeaks(processedSignalArray, mainFile):
    longPeak = peakMax(processedSignalArray)
    abf = pyabf.ABF(mainFile)
    traceLen = len(abf.sweepList)
    samplingFreqMSec = abf.dataPointsPerMs + (1/3)
    samplingFreqSec = samplingFreqMSec * 1000
    peaksDict, finalDict = {}, {}
    peaksArray, adjustedArea, decayTauList, riseTauList = np.zeros((traceLen, longPeak)), np.zeros((traceLen, longPeak)), np.zeros((traceLen, longPeak)), np.zeros((traceLen, longPeak))
    widthBottomArray = np.zeros((traceLen, 4, longPeak))
    width10Array = np.zeros((traceLen, 4, longPeak))
    width90Array = np.zeros((traceLen, 4, longPeak))
    overlapPeaks = [[0]*longPeak for _ in abf.sweepList]
    for index, traces in enumerate(processedSignalArray): # Finds peaks in a signal
        peaks, peaksDict[index] = sci.find_peaks(traces, prominence= 0.05, width=0, wlen=20000, rel_height= 0.5)
        bottomWidth = sci.peak_widths(traces, peaks, rel_height=1, prominence_data=(peaksDict[index]['prominences'], peaksDict[index]["left_bases"], 
                                                                                     peaksDict[index]["right_bases"]), wlen=20000)
        width10 = sci.peak_widths(traces, peaks, rel_height=0.1, prominence_data=(peaksDict[index]['prominences'], 
                                                                                  peaksDict[index]["left_bases"], peaksDict[index]["right_bases"]), wlen=20000)
        width90 = sci.peak_widths(traces, peaks, rel_height=0.9, prominence_data=(peaksDict[index]['prominences'], 
                                                                                  peaksDict[index]["left_bases"], peaksDict[index]["right_bases"]), wlen=20000)
        peaks = np.pad(peaks, pad_width= (0, longPeak - len(peaks)), mode= 'constant', constant_values= 0)
        bottomWidth, width10, width90 = np.array(bottomWidth), np.array(width10), np.array(width90)
        peaksArray[index] = peaks
        for i, param in enumerate(bottomWidth):
            param = np.pad(param, pad_width= (0, longPeak - len(param)), mode= 'constant', constant_values= 0)
            widthBottomArray[index][i] = param
        for i, param in enumerate(width10):
            param = np.pad(param, pad_width= (0, longPeak - len(param)), mode= 'constant', constant_values= 0)
            width10Array[index][i] = param
        for i, param in enumerate(width90):
            param = np.pad(param, pad_width= (0, longPeak - len(param)), mode= 'constant', constant_values= 0)
            width90Array[index][i] = param
        for i in peaksDict[index]:
           paddedEntry = np.pad(peaksDict[index][i], pad_width= (0, longPeak - len(peaksDict[index][i])), mode= 'constant', constant_values= 0)
           peaksDict[index][i] = paddedEntry
    for k, checkPeak in np.ndenumerate(peaksArray): # Finds peaks that have overlap with one another
            overlapPeaks[k[0]][k[1]] = np.array([y for y in peaksArray[k[0]] if ((widthBottomArray[k[0]][2][k[1]] < y < widthBottomArray[k[0]][3][k[1]]) and y != checkPeak)])
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
            if len(processedSignalArray[z][int(widthBottomArray[z][2][i]):int(widthBottomArray[z][3][i])]) == 0:
                 continue
            degreeNPeaks[peakOfDegree] = len(overlapPeaks[z][i])
        sortedDegreeNPeaks = dict(sorted(degreeNPeaks.items(), key=lambda item: item[1]))
        for _, u in enumerate(sortedDegreeNPeaks): # Generates event area
            i = np.where(x == u)
            if len(processedSignalArray[z][int(widthBottomArray[z][2][i[0][0]]):int(widthBottomArray[z][3][i[0][0]])]) == 0:
                 continue
            peaksInArea = overlapPeaks[z][i[0][0]]
            peakArea = inte.simpson(y=processedSignalArray[z][int(widthBottomArray[z][2][i[0][0]]):int(widthBottomArray[z][3][i[0][0]])], 
                                    x=np.arange(int(widthBottomArray[z][2][i[0][0]]), int(widthBottomArray[z][3][i[0][0]]))/samplingFreqMSec)
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
        for i, p in enumerate(x): # Determines rise and decay time constants (tau) for each peak in a sweep
            if len(processedSignalArray[z][int(widthBottomArray[z][2][i]):int(widthBottomArray[z][3][i])]) == 0:
                continue
            decayTau = np.array(processedSignalArray[z][int(p):int(width90Array[z][3][i])])
            riseTau = np.array(processedSignalArray[z][int(width90Array[z][2][i]):int(p)])
            decWidth = int(len(decayTau))
            riseWidth = int(len(riseTau))
    #TODO implement time analysis across the a data point.
        # Event Decay
            p0 = (processedSignalArray[z][int(p)], -1)
            if len(decayTau) == 0:
                continue
            else:
                try:
                    decArray = np.array(list(range(0, decWidth)), dtype=np.float64)
                    popt, _ = opt.curve_fit(lambda t, a, b: a * np.exp((b/samplingFreqSec) * (t/samplingFreqSec)), decArray, decayTau, p0=p0, 
                                            bounds=opt.Bounds(lb=[-np.inf, -np.inf], 
                                                              ub=[np.inf, np.inf]),
                                                              maxfev=35000, xtol=1e-6, ftol=1e-6)
                    if popt[1] - p == 0:
                        continue
                    decayTauList[z][i] = abs(1/(popt[1]/1000))
                except RuntimeError as e:
                    if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                        decayTauList[z][i] = 0
            # Event Rise
            p0 = (processedSignalArray[z][int(width90Array[z][2][i])], 1)
            if len(riseTau) == 0:
                continue
            else:
                try:
                    riseArray = np.array(list(range(0, riseWidth)), dtype=np.float64)
                    popt, _ = opt.curve_fit(lambda t, a, b: a * np.exp((b/samplingFreqSec) * (t/samplingFreqSec)), riseArray, riseTau, p0=p0, 
                                            bounds=opt.Bounds(lb=[-np.inf, -np.inf], 
                                                              ub=[np.inf, np.inf]),
                                                              maxfev=35000, xtol=1e-6, ftol=1e-6)
                    if popt[1] - p == 0:
                        continue
                    riseTauList[z][i] = abs(1/(popt[1]/1000))
                except RuntimeError as e:
                    if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                        decayTauList[z][i] = 0
            peakTable.Decay_Tau_exp = pd.Series(decayTauList[z])
            peakTable.Rise_Tau_exp = pd.Series(riseTauList[z])
        print("Trace #%i done!"%z)
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
#TODO Remove overlapping peaks from consideration of plotting the line
def peakDisplay(processedSignalArray, mainFile, ratSide):
    abf = pyabf.ABF(mainFile)
    samplingFreqMSec = abf.dataPointsPerMs + (1/3)
    samplingFreqSec = samplingFreqMSec * 1000
    seconds = np.arange(0, 47750, samplingFreqSec)
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, wlen= 20000)
    widthBottom = sci.peak_widths(processedSignalArray, peaks, rel_height=1, prominence_data=(peaksDict['prominences'], 
                                                                                             peaksDict["left_bases"], peaksDict["right_bases"]), wlen=20000)
    widthHalf = sci.peak_widths(processedSignalArray, peaks, rel_height=0.5, prominence_data=(peaksDict['prominences'], 
                                                                                             peaksDict["left_bases"], peaksDict["right_bases"]), wlen=20000)
    width10 = sci.peak_widths(processedSignalArray, peaks, rel_height=0.1, prominence_data=(peaksDict['prominences'], 
                                                                                             peaksDict["left_bases"], peaksDict["right_bases"]), wlen=20000)
    width90 = sci.peak_widths(processedSignalArray, peaks, rel_height=0.9, prominence_data=(peaksDict['prominences'], 
                                                                                             peaksDict["left_bases"], peaksDict["right_bases"]), wlen=20000)
    # for k, checkPeak in enumerate(peaks): # Finds peaks that have overlap with one another
    #         overlapPeaks[k[0]][k[1]] = np.array([y for y in peaksArray[k[0]] if ((widthArray[k[0]][2][k[1]] < y < widthArray[k[0]][3][k[1]]) and y != checkPeak)])                                                                                       
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
        p0 = (processedSignalArray[int(x)], -1)
        # Event Decay
        decayTau = np.array(processedSignalArray[int(x):int(width90[3][i])])
        decWidth = int(len(decayTau))
        if len(decayTau) == 0:
            continue
        else:
            decArray = np.array(list(range(0, decWidth)))
            
            popt, _ = opt.curve_fit(lambda t, a, b: a * np.exp((b/samplingFreqSec) * (t/samplingFreqSec)), decArray, decayTau, p0=p0, 
                                    bounds=opt.Bounds(lb=[-np.inf, -np.inf], 
                                                      ub=[np.inf, np.inf]), 
                                                      maxfev=35000, xtol=1e-6, ftol=1e-6)
            x_dec = np.linspace(np.min(decArray), np.max(decArray), 1000)
            a = popt[0]
            b = popt[1]
            y_dec = a * np.exp((b/samplingFreqSec) * (x_dec/samplingFreqSec))
            peakFig.plot(x_dec+x, y_dec, color="C9")
            print("Tau (Decay):", abs(1/(popt[1]/1000)))
        p0 = (processedSignalArray[int(width90[2][i])], 1)
        # Event Rise
        riseTau = np.array(processedSignalArray[int(width90[2][i]):x])
        riseWidth = int(len(riseTau))
        if len(riseTau) == 0:
            continue
        else:
            riseArray = np.array(list(range(0, riseWidth)))
            popt, _ = opt.curve_fit(lambda t, a, b: a * np.exp((b/samplingFreqSec) * (t/samplingFreqSec)), riseArray, riseTau, p0=p0, bounds=opt.Bounds(lb=[-np.inf, -np.inf], 
                                                                                            ub=[np.inf, np.inf]), maxfev=35000, xtol=1e-6, ftol=1e-6)
            x_rise = np.linspace(np.min(riseArray), np.max(riseArray), 1000)
            a = popt[0]
            b = popt[1]
            y_rise = a * np.exp((b/samplingFreqSec) * (x_rise/samplingFreqSec))
            peakFig.plot(x_rise+width90[2][i], y_rise, color="C8")
            print("Tau (Rise):", abs(1/(popt[1]/1000)))

    
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
