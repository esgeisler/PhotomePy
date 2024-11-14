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
    peaksArray, adjustedArea, riseTauList, decayTauList = (np.zeros((traceLen, longPeak)) for i in range(4))
    widthBottomArray = np.zeros((traceLen, 4, longPeak))
    width10Array = np.zeros((traceLen, 4, longPeak))
    width90Array = np.zeros((traceLen, 4, longPeak))
    overlapPeaks, overlapRise, overlapDecay = [[0]*longPeak for _ in abf.sweepList], [[0]*longPeak for _ in abf.sweepList], [[0]*longPeak for _ in abf.sweepList]
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
            overlapRise[k[0]][k[1]] = np.array([y for y in peaksArray[k[0]] if ((widthBottomArray[k[0]][2][k[1]] < y < checkPeak) and y != checkPeak)])
            overlapDecay[k[0]][k[1]] = np.array([y for y in peaksArray[k[0]] if ((checkPeak < y < widthBottomArray[k[0]][3][k[1]]) and y != checkPeak)])
    for z, x in enumerate(peaksArray):
        degreeNPeaks, decayNPeaks, riseNPeaks = {}, {}, {}
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
            riseNPeaks[peakOfDegree] = len(overlapRise[z][i])
            decayNPeaks[peakOfDegree] = len(overlapDecay[z][i])
        sortedDegreeNPeaks = dict(sorted(degreeNPeaks.items(), key=lambda item: item[1]))
        sortedDegreeRise = dict(sorted(riseNPeaks.items(), key=lambda item: item[1]))
        sortedDegreeDecay = dict(sorted(decayNPeaks.items(), key=lambda item: item[1]))
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
        for _, u in enumerate(sortedDegreeRise): #Generates adjusted rise time constants (tau)
            i = np.where(x == u)
            if len(processedSignalArray[z][int(widthBottomArray[z][2][i[0][0]]):int(widthBottomArray[z][3][i[0][0]])]) == 0:
                continue
            peaksInRise = overlapRise[z][i[0][0]]
            adjustedRiseTau = np.array(processedSignalArray[z][int(width90Array[z][2][i[0][0]]):int(u)])
            riseWidth = int(len(adjustedRiseTau))
            riseArray = np.array(list(range(0, riseWidth)))
            p0 = (processedSignalArray[z][int(width90Array[z][2][i[0][0]])], 1)
            if len(adjustedRiseTau) == 0:
                continue
            else:
                try:
                    match riseNPeaks[u]:
                        case 0:
                            popt, _ = opt.curve_fit(lambda t, a, b: a * np.exp((b/samplingFreqSec) * (t/samplingFreqSec)), riseArray, adjustedRiseTau, p0=p0, 
                                                    bounds=opt.Bounds(lb=[-np.inf, -np.inf], 
                                                                    ub=[np.inf, np.inf]),
                                                                    maxfev=35000, xtol=1e-6, ftol=1e-6)
                            if popt[1] - u == 0:
                                continue
                            riseTauList[z][i[0][0]] = abs(1/(popt[1]/1000))
                        case _:
                            for l in peaksInRise:
                                j = np.where(x == l)
                                adjustedRiseTau = np.delete(adjustedRiseTau, np.s_[int(widthBottomArray[z][2][j[0][0]]):int(widthBottomArray[z][3][j[0][0]])])
                                riseWidth = int(len(adjustedRiseTau))
                                riseArray = np.array(list(range(0, riseWidth)))
                            popt, _ = opt.curve_fit(lambda t, a, b: a * np.exp((b/samplingFreqSec) * (t/samplingFreqSec)), riseArray, adjustedRiseTau, p0=p0, 
                                                    bounds=opt.Bounds(lb=[-np.inf, -np.inf], 
                                                                    ub=[np.inf, np.inf]),
                                                                    maxfev=35000, xtol=1e-6, ftol=1e-6)
                            if popt[1] - u == 0:
                                continue
                            riseTauList[z][i[0][0]] = abs(1/(popt[1]/1000))
                except RuntimeError as e:
                    if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                        riseTauList[z][i[0][0]] = 0 
        peakTable.Rise_Tau_exp = pd.Series(riseTauList[z])
        for _, u in enumerate(sortedDegreeDecay):
            i = np.where(x == u)
            peaksInDecay = overlapDecay[z][i[0][0]]
            adjustedDecayTau = np.array(processedSignalArray[z][int(u):int(width90Array[z][3][i[0][0]])])
            decWidth = int(len(adjustedDecayTau))
            decArray = np.array(list(range(0, decWidth)))
            p0 = (processedSignalArray[z][int(u)], -1)
            if len(adjustedDecayTau) == 0:
                continue
            else:
                try:
                    match decayNPeaks[u]:
                        case 0:
                            popt, _ = opt.curve_fit(lambda t, a, b: a * np.exp((b/samplingFreqSec) * (t/samplingFreqSec)), decArray, adjustedDecayTau, p0=p0, 
                                                    bounds=opt.Bounds(lb=[-np.inf, -np.inf], 
                                                                    ub=[np.inf, np.inf]),
                                                                    maxfev=35000, xtol=1e-6, ftol=1e-6)
                            if popt[1] - u == 0:
                                continue
                            decayTauList[z][i[0][0]] = abs(1/(popt[1]/1000))
                        case _:
                            for l in peaksInDecay:
                                j = np.where(x == l)
                                adjustedDecayTau = np.delete(adjustedDecayTau, np.s_[int(widthBottomArray[z][2][j[0][0]]):int(widthBottomArray[z][3][j[0][0]])])
                                decWidth = int(len(adjustedDecayTau))
                                decArray = np.array(list(range(0, decWidth)))
                            popt, _ = opt.curve_fit(lambda t, a, b: a * np.exp((b/samplingFreqSec) * (t/samplingFreqSec)), decArray, adjustedDecayTau, p0=p0, 
                                                    bounds=opt.Bounds(lb=[-np.inf, -np.inf], 
                                                                    ub=[np.inf, np.inf]),
                                                                    maxfev=35000, xtol=1e-6, ftol=1e-6)
                            if popt[1] - u == 0:
                                continue
                            decayTauList[z][i[0][0]] = abs(1/(popt[1]/1000))
                except RuntimeError as e:
                    if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                        decayTauList[z][i[0][0]] = 0
        peakTable.Decay_Tau_exp = pd.Series(decayTauList[z])
        match np.size(x): # Determines frequency and total area of a sweep
            case 0:
                peakTable.Frequency = 0
                peakTable.Total_Area = 0
            case _:
                peakTable.Frequency.iat[0] = round(np.count_nonzero(x)/((len(processedSignalArray[z]) + 2250)/samplingFreqSec), 2) #Peaks/second (15 second trace)
                peakTable.Total_Area.iat[0] = sum(adjustedArea[z])
        peakTable.drop(peakTable[peakTable.Peak_Index == 0].index, inplace= True)
        print("Trace #%i done!"%(z+1))
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
    decayNPeaks, riseNPeaks = {}, {}
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, wlen= 20000)
    overlapRise, overlapDecay = [0 for _ in range(len(peaks))], [0 for _ in range(len(peaks))]
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
    
    for k, checkPeak in enumerate(peaks): # Finds peaks that have overlap with one another
        overlapRise[k] = [y for y in peaks if ((widthBottom[2][k] < y < checkPeak) and y != checkPeak)]
        overlapDecay[k] = [y for y in peaks if ((checkPeak < y < widthBottom[3][k]) and y != checkPeak)]
    for i, peakOfDegree in enumerate(peaks): # Determines the "degree" of a peak; how many peaks overlap with it
        if len(processedSignalArray[int(widthBottom[2][i]):int(widthBottom[3][i])]) == 0:
                continue
        riseNPeaks[peakOfDegree] = len(overlapRise[i])
        decayNPeaks[peakOfDegree] = len(overlapDecay[i])
    sortedDegreeRise = dict(sorted(riseNPeaks.items(), key=lambda item: item[1]))
    sortedDegreeDecay = dict(sorted(decayNPeaks.items(), key=lambda item: item[1]))
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
    for _, u in enumerate(sortedDegreeRise): #Generates adjusted rise time constants (tau)
        i = np.where(peaks == u)
        if len(processedSignalArray[int(widthBottom[2][i[0][0]]):int(widthBottom[3][i[0][0]])]) == 0:
            continue
        peaksInRise = overlapRise[i[0][0]]
        adjustedRiseTau = np.array(processedSignalArray[int(width90[2][i[0][0]]):int(u)])
        riseWidth = int(len(adjustedRiseTau))
        riseArray = np.array(list(range(0, riseWidth)))
        p0 = (processedSignalArray[int(width90[2][i[0][0]])], 1, processedSignalArray[int(width10[2][i[0][0]])])
        if len(adjustedRiseTau) == 0:
            continue
        else:
            try:
                match riseNPeaks[u]:
                    case 0:
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, riseArray/samplingFreqSec, adjustedRiseTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, 0, 0], 
                                                                ub=[np.inf, np.inf, np.inf]),
                                                                maxfev=1000)
                        # if popt[1] - u == 0:
                        #     continue
                        x_rise = np.linspace(np.min(riseArray), np.max(riseArray), 1000)
                        a = popt[0]
                        b = popt[1]
                        c = popt[2]
                        print(a,b, x_rise[-1]/samplingFreqSec, "rise")
                        y_rise = a * np.exp(b * (x_rise/samplingFreqSec)) + c
                        peakFig.plot(x_rise+width90[2][i], y_rise, color="C8")
                        # print("Tau (Rise):", abs(1/(popt[1]/1000)))
                    case _:
                        for l in peaksInRise:
                            j = np.where(peaks == l)
                            r = np.where(adjustedRiseTau == processedSignalArray[int(widthBottom[2][j[0][0]])])
                            o = np.where(adjustedRiseTau == processedSignalArray[int(widthBottom[3][j[0][0]])])
                            if len(r[0]) == 0 or len(o[0]) == 0:
                                continue
                            averageFloor = np.linspace(int(widthBottom[2][j[0][0]]), int(widthBottom[3][j[0][0]]), len(processedSignalArray[int(widthBottom[2][j[0][0]]):int(widthBottom[3][j[0][0]])]))
                            for y, _ in enumerate(adjustedRiseTau[int(r[0][0]):int(o[0][0])]):
                                adjustedRiseTau[y+r[0][0]] = widthBottom[1][j[0][0]]
                            peakFig.plot(averageFloor, adjustedRiseTau[int(r[0][0]):int(o[0][0])], color="C11")
                        riseWidth = int(len(adjustedRiseTau))
                        riseArray = np.array(list(range(0, riseWidth)))
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp((b * t)) + c, riseArray/samplingFreqSec, adjustedRiseTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, 0, 0], 
                                                                ub=[np.inf, np.inf, np.inf]),
                                                                maxfev=1000)
                        # if popt[1] - u == 0:
                        #     continue
                        x_rise = np.linspace(np.min(riseArray), np.max(riseArray), 1000)
                        a = popt[0]
                        b = popt[1]
                        c = popt[2]
                        print(a,b, x_rise[-1]/samplingFreqSec, "rise")
                        y_rise = a * np.exp(b * (x_rise/samplingFreqSec)) + c
                        peakFig.plot(x_rise+width90[2][i], y_rise, color="C8")
                        # print("Tau (Rise):", abs(1/(popt[1]/1000)))
            except RuntimeError as e:
                if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                    continue 
    for _, u in enumerate(sortedDegreeDecay):
        i = np.where(peaks == u)
        peaksInDecay = overlapDecay[i[0][0]]
        adjustedDecayTau = np.array(processedSignalArray[int(u):int(width90[3][i[0][0]])])
        decWidth = int(len(adjustedDecayTau))
        decArray = np.array(list(range(0, decWidth)))
        p0 = (processedSignalArray[int(u)], -1, processedSignalArray[int(width90[3][i[0][0]])])
        if len(adjustedDecayTau) == 0:
            continue
        else:
            try:
                match decayNPeaks[u]:
                    case 0:
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, decArray/samplingFreqSec, adjustedDecayTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, -np.inf, 0], 
                                                                ub=[np.inf, 0, np.inf]),
                                                                maxfev=1000)
                        # if popt[1] - u == 0:
                        #     continue
                        x_dec = np.linspace(np.min(decArray), np.max(decArray), 1000)
                        a = popt[0]
                        b = popt[1]
                        c = popt[2]
                        print(a,b, x_dec[-1]/samplingFreqSec, "decay")
                        y_dec = a * np.exp(b * ((x_dec/samplingFreqSec))) + c
                        peakFig.plot(x_dec+u, y_dec, color="C9")
                        # print("Tau (Decay):", abs(1/(popt[1]/1000)))
                    case _:
                        for l in peaksInDecay:
                            j = np.where(peaks == l)
                            r = np.where(adjustedDecayTau == processedSignalArray[int(widthBottom[2][j[0][0]])])
                            o = np.where(adjustedDecayTau == processedSignalArray[int(widthBottom[3][j[0][0]])])
                            if len(r[0]) == 0 or len(o[0]) == 0:
                                continue
                            averageFloor = np.linspace(int(widthBottom[2][j[0][0]]), int(widthBottom[3][j[0][0]]), len(processedSignalArray[int(widthBottom[2][j[0][0]]):int(widthBottom[3][j[0][0]])]))
                            for y, _ in enumerate(adjustedDecayTau[int(r[0][0]):int(o[0][0])]):
                                adjustedDecayTau[y+r[0][0]] = widthBottom[1][j[0][0]]
                            peakFig.plot(averageFloor, adjustedDecayTau[int(r[0][0]):int(o[0][0])], color="C11")
                        decWidth = int(len(adjustedDecayTau))
                        decArray = np.array(list(range(0, decWidth)))
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp((b * t)) + c, decArray/samplingFreqSec, adjustedDecayTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, -np.inf, 0], 
                                                                ub=[np.inf, 0, np.inf]),
                                                                maxfev=1000)
                        # if popt[1] - u == 0:
                        #     continue
                        x_dec = np.linspace(np.min(decArray), np.max(decArray), 1000)
                        a = popt[0]
                        b = popt[1]
                        c = popt[2]
                        print(a,b, x_dec[-1]/samplingFreqSec, "decay")
                        y_dec = a * np.exp(b * ((x_dec/samplingFreqSec))) + c
                        peakFig.plot(x_dec+u, y_dec, color="C9")
                        # print("Tau (Decay):", abs(1/(popt[1]/1000)))
            except RuntimeError as e:
                if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                    continue
        # p0 = (processedSignalArray[int(x)], -1)
        # # Event Decay
        # decayTau = np.array(processedSignalArray[int(x):int(width90[3][i])])
        # decWidth = int(len(decayTau))
        # if len(decayTau) == 0:
        #     continue
        # else:
        #     decArray = np.array(list(range(0, decWidth)))
            
        #     popt, _ = opt.curve_fit(lambda t, a, b: a * np.exp((b/samplingFreqSec) * (t/samplingFreqSec)), decArray, decayTau, p0=p0, 
        #                             bounds=opt.Bounds(lb=[-np.inf, -np.inf], 
        #                                               ub=[np.inf, np.inf]), 
        #                                               maxfev=35000, xtol=1e-6, ftol=1e-6)
        #     x_dec = np.linspace(np.min(decArray), np.max(decArray), 1000)
        #     a = popt[0]
        #     b = popt[1]
        #     y_dec = a * np.exp((b/samplingFreqSec) * (x_dec/samplingFreqSec))
        #     peakFig.plot(x_dec+x, y_dec, color="C9")
        #     print("Tau (Decay):", abs(1/(popt[1]/1000)))
        # p0 = (processedSignalArray[int(width90[2][i])], 1)
        # # Event Rise
        # riseTau = np.array(processedSignalArray[int(width90[2][i]):x])
        # riseWidth = int(len(riseTau))
        # if len(riseTau) == 0:
        #     continue
        # else:
        #     riseArray = np.array(list(range(0, riseWidth)))
        #     popt, _ = opt.curve_fit(lambda t, a, b: a * np.exp((b/samplingFreqSec) * (t/samplingFreqSec)), riseArray, riseTau, p0=p0, bounds=opt.Bounds(lb=[-np.inf, -np.inf], 
        #                                                                                     ub=[np.inf, np.inf]), maxfev=35000, xtol=1e-6, ftol=1e-6)
        #     x_rise = np.linspace(np.min(riseArray), np.max(riseArray), 1000)
        #     a = popt[0]
        #     b = popt[1]
        #     y_rise = a * np.exp((b/samplingFreqSec) * (x_rise/samplingFreqSec))
        #     peakFig.plot(x_rise+width90[2][i], y_rise, color="C8")
        #     print("Tau (Rise):", abs(1/(popt[1]/1000)))

    
    peakFig.hlines(*widthHalf[1:], color="C6")
    # peakFig.hlines(*widthBottom[1:], color="C7")
    peakFig.vlines(x=peaks, ymin=processedSignalArray[peaks] - peaksDict["prominences"], ymax=processedSignalArray[peaks], color="C5")
    peakFig.set_title(ratSide)
    plt.axis([0, 47750, 0, 2])
    plt.xticks(ticks=seconds, labels=range(len(seconds)))
    plt.xlabel("Time (s)")
    plt.ylabel("Fluorescence (AU)")
    plt.minorticks_on()
    plt.show()
