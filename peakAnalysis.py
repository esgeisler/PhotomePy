import scipy.signal as sci
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyabf
import sys
np.set_printoptions(threshold=sys.maxsize)

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
        peaks, _ = sci.find_peaks(t, prominence= 0.05, wlen=20000)
        peakLenArray[i] = np.size(peaks)
    longestPeak = np.max(peakLenArray)
    return int(longestPeak), peakLenArray

# Finds peaks in a signal and provides their widths, amplitudes, avg. frequencies, and areas across an entire .abf file
def wholeTracePeaks(processedSignalArray, mainFile):
    finalDict = {}
    for z, _ in enumerate(processedSignalArray):
        peaks = trace.TracePeaks(processedSignalArray, mainFile, z)
        peaks.peakFinder(processedSignalArray[z])
        peaks.overlapCheck(peaks.peaks)
        peakTable = pd.DataFrame(columns= ['Event_Num', 'Peak_Index', 'Peak_Time_Sec', 
                                           'Event_Window_Start', 'Event_Window_End', 
                                           'Amplitude', 
                                           'Off_Time_ms', 'Width_at50_ms',
                                           'Frequency', 
                                           'Avg_Area', 'Total_Area',
                                           'Decay_Tau_exp', 'Rise_Tau_exp'])
        peakTable.Event_Num = peaks.peakNum
        peakTable.Peak_Index = peaks.peakIndex
        peakTable.Peak_Time_Sec = peaks.peakLocSec
        peakTable.Event_Window_Start = peaks.leftBounds
        peakTable.Event_Window_End = peaks.rightBounds
        peakTable.Amplitude = peaks.amplitude
        peakTable.Off_Time_ms = peaks.rightTail
        peakTable.Width_at50_ms = peaks.width

        peaks.freqSet()
        peaks.areaSet()
        peaks.riseSet()
        peaks.decaySet()
        peakTable.Avg_Area = peaks.meanArea
        peakTable.Rise_Tau_exp = pd.Series(peaks.riseRate)
        peakTable.Decay_Tau_exp = pd.Series(peaks.decayRate)
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
def peakDisplay(processedSignalArray, mainFile, ratSide):
    seconds, _, samplingFreqSec = secondsCalculator(mainFile)
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
        overlapRise[k] = [y for y in peaks if ((width90[2][k] < y < checkPeak) and y != checkPeak)]
        overlapDecay[k] = [y for y in peaks if ((checkPeak < y < width90[3][k]) and y != checkPeak)]
    for i, peakOfDegree in enumerate(peaks): # Determines the "degree" of a peak; how many peaks overlap with it
        if len(processedSignalArray[int(widthBottom[2][i]):int(widthBottom[3][i])]) == 0:
                continue
        riseNPeaks[peakOfDegree] = len(overlapRise[i])
        decayNPeaks[peakOfDegree] = len(overlapDecay[i])
    sortedDegreeRise = dict(sorted(riseNPeaks.items(), key=lambda item: item[1]))
    sortedDegreeDecay = dict(sorted(decayNPeaks.items(), key=lambda item: item[1]))
    for i, x in enumerate(peaks):
        peakFig.plot(int(width90[3][i]), processedSignalArray[int(width90[3][i])], 'g.')
        peakFig.annotate("Trough (%i)"%(i+1), xycoords= 'data', size= 8, horizontalalignment= 'center',
                     xytext = (int(width90[3][i]), processedSignalArray[int(width90[3][i])] - 0.3), 
                     xy = (int(width90[3][i]), processedSignalArray[int(width90[3][i])] - 0.01),
                     arrowprops=dict(facecolor= 'black', width= 1, headwidth= 5, headlength= 5))
        peakFig.annotate("p%i"%(i+1), xycoords= 'data', size= 8, horizontalalignment= 'center',
                         xytext= (x, processedSignalArray[x] + 0.01),
                         xy = (x, processedSignalArray[x]))
        peakFig.fill_between(np.arange(int(width90[2][i]), int(width90[3][i])), processedSignalArray[int(width90[2][i]):int(width90[3][i])], 
                             width90[1][i], color="C1", alpha=0.3)
    for _, u in enumerate(sortedDegreeRise): #Generates and plots exponential rise functions for peaks
        i = np.where(peaks == u)
        if len(processedSignalArray[int(widthBottom[2][i[0][0]]):int(widthBottom[3][i[0][0]])]) == 0:
            continue
        peaksInRise = overlapRise[i[0][0]]
        adjustedRiseTau = np.array(processedSignalArray[int(width90[2][i[0][0]]):int(width10[2][i[0][0]])])
        riseWidth = int(len(adjustedRiseTau))
        riseArray = np.array(list(range(0, riseWidth)))
        p0 = (processedSignalArray[int(width10[2][i[0][0]])] - processedSignalArray[int(width90[2][i[0][0]])], 1, processedSignalArray[int(width90[2][i[0][0]])])
        if len(adjustedRiseTau) == 0:
            continue
        else:
            try:
                match riseNPeaks[u]:
                    case 0: # Handles peaks with no overlap
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, riseArray/samplingFreqSec, adjustedRiseTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, 0, processedSignalArray[int(width90[2][i[0][0]])]], 
                                                                ub=[processedSignalArray[int(width10[2][i[0][0]])] - processedSignalArray[int(width90[2][i[0][0]])], np.inf, np.inf]),
                                                                maxfev=1000)
                        x_rise = np.linspace(np.min(riseArray), np.max(riseArray), 1000)
                        a, b, c = popt[0], popt[1], popt[2]
                        squaredDiffs = np.square(adjustedRiseTau - (a * np.exp(b * ((riseArray/samplingFreqSec))) + c))
                        squaredDiffsFromMean = np.square(adjustedRiseTau - np.mean(adjustedRiseTau))
                        rSquared = 1 - (np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean))
                        if rSquared < 0.8:
                            # print("Fit quality poor; not shown")
                            continue
                        print(f"R² = {rSquared}", "rise")
                        y_rise = a * np.exp(b * (x_rise/samplingFreqSec)) + c
                        peakFig.plot(x_rise+width90[2][i], y_rise, color="C8")
                        print("Tau (Rise):", abs(1/popt[1]))
                    case _: #Handles the rest of the peaks
                        # ji = np.where(peaks == peaksInRise[-1])
                        # adjustedRiseTau = np.array(processedSignalArray[int(widthBottom[3][ji[0][0]]):int(width10[2][i[0][0]])])
                        # p0 = (int(processedSignalArray[int(width10[2][i[0][0]])] - processedSignalArray[int(widthBottom[3][ji[0][0]])]), 1, processedSignalArray[int(widthBottom[3][ji[0][0]])])
                        # for l in peaksInRise:
                        #     j = np.where(peaks == l)
                        #     r = np.where(adjustedRiseTau == processedSignalArray[int(width90[2][j[0][0]])])
                        #     o = np.where(adjustedRiseTau == processedSignalArray[int(width90[3][j[0][0]])])
                        #     if len(r[0]) == 0 or len(o[0]) == 0:
                        #         continue
                        #     for y, _ in enumerate(adjustedRiseTau[int(r[0][0]):int(o[0][0])]):
                        #         adjustedRiseTau[y+r[0][0]] = width90[1][j[0][0]]
                        # riseWidth = int(len(adjustedRiseTau))
                        # riseArray = np.array(list(range(0, riseWidth)))
                        # popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp((b * t)) + c, riseArray/samplingFreqSec, adjustedRiseTau, p0=p0, 
                        #                         bounds=opt.Bounds(lb=[0, 0, processedSignalArray[int(widthBottom[3][ji[0][0]])]], 
                        #                                         ub=[np.inf, np.inf, processedSignalArray[int(width10[2][i[0][0]])]]),
                        #                                         maxfev=1000)
                        # x_rise = np.linspace(np.min(riseArray), np.max(riseArray), 1000)
                        # a, b, c = popt[0], popt[1], popt[2]
                        # squaredDiffs = np.square(adjustedRiseTau - (a * np.exp(b * ((riseArray/samplingFreqSec))) + c))
                        # squaredDiffsFromMean = np.square(adjustedRiseTau - np.mean(adjustedRiseTau))
                        # rSquared = 1 - (np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean))
                        # if rSquared < 0.8:
                        #     # print("Fit quality poor; not shown")
                        #     continue
                        # print(f"R² = {rSquared}", "rise")
                        # y_rise = a * np.exp(b * (x_rise/samplingFreqSec)) + c
                        # peakFig.plot(x_rise+widthBottom[3][ji], y_rise, color="C8")
                        # print("Tau (Rise):", abs(1/popt[1]))
                        continue
            except RuntimeError as e:
                if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                    print("Too Many Guesses")
                    continue 
                else:
                    print("something else")
                    continue
            except ValueError as e:
                if str(e) == "'x0' is infeasible":
                    print("Bad Initial Value")
                    continue 
                elif str(e) == "Each lower bound must be strictly less than each upper bound.":
                    print("Lower bound is larger than upper bound.")
                    continue
                
    for _, u in enumerate(sortedDegreeDecay): #Generates and plots expoential decay functions for peaks
        i = np.where(peaks == u)
        peaksInDecay = overlapDecay[i[0][0]]
        adjustedDecayTau = np.array(processedSignalArray[int(width10[3][i[0][0]]):int(width90[3][i[0][0]])])
        decWidth = int(len(adjustedDecayTau))
        decArray = np.array(list(range(0, decWidth)))
        p0 = (processedSignalArray[int(width10[3][i[0][0]])], -1, processedSignalArray[int(width90[3][i[0][0]])])
        if len(adjustedDecayTau) == 0:
            continue
        else:
            try:
                match decayNPeaks[u]:
                    case 0: #Handles peaks with no overlap in their decay slope
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, decArray/samplingFreqSec, adjustedDecayTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, -np.inf, processedSignalArray[int(width90[3][i[0][0]])]], 
                                                                ub=[processedSignalArray[int(width10[3][i[0][0]])], 0, np.inf]),
                                                                maxfev=1000)
                        x_dec = np.linspace(np.min(decArray), np.max(decArray), 1000)
                        a, b, c = popt[0], popt[1], popt[2]
                        y_dec = a * np.exp(b * ((x_dec/samplingFreqSec))) + c
                        squaredDiffs = np.square(adjustedDecayTau - (a * np.exp(b * ((decArray/samplingFreqSec))) + c))
                        squaredDiffsFromMean = np.square(adjustedDecayTau - np.mean(adjustedDecayTau))
                        rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
                        if rSquared < 0.8:
                            # print("Fit quality poor; not shown")
                            continue
                        print(f"R² = {rSquared}", "decay")
                        peakFig.plot(x_dec+int(width10[3][i[0][0]]), y_dec, color="C9")
                        print("Tau (Decay):", abs(1/popt[1]))
                    case _: #Handles the rest
                        # ji = np.where(peaks == peaksInDecay[0])
                        # adjustedDecayTau = np.array(processedSignalArray[int(width10[3][i[0][0]]):int(width90[2][ji[0][0]])])
                        # p0 = (processedSignalArray[int(width10[3][i[0][0]])], -1, processedSignalArray[int(width90[2][ji[0][0]])])
                        # for l in peaksInDecay:
                        #     j = np.where(peaks == l)
                        #     r = np.where(adjustedDecayTau == processedSignalArray[int(width90[2][j[0][0]])])
                        #     o = np.where(adjustedDecayTau == processedSignalArray[int(width90[3][j[0][0]])])
                        #     if len(r[0]) == 0 or len(o[0]) == 0:
                        #         continue
                        #     for y, _ in enumerate(adjustedDecayTau[int(r[0][0]):int(o[0][0])]):
                        #         adjustedDecayTau[y+r[0][0]] = width90[1][j[0][0]]
                        # decWidth = int(len(adjustedDecayTau))
                        # decArray = np.array(list(range(0, decWidth)))
                        # popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp((b * t)) + c, decArray/samplingFreqSec, adjustedDecayTau, p0=p0, 
                        #                         bounds=opt.Bounds(lb=[0, -np.inf, processedSignalArray[int(widthBottom[2][ji[0][0]])]], 
                        #                                         ub=[processedSignalArray[int(width10[3][i[0][0]])], 0, np.inf]),
                        #                                         maxfev=1000)
                        # x_dec = np.linspace(np.min(decArray), np.max(decArray), 1000)
                        # a, b, c = popt[0], popt[1], popt[2]
                        # squaredDiffs = np.square(adjustedDecayTau - (a * np.exp(b * ((decArray/samplingFreqSec))) + c))
                        # squaredDiffsFromMean = np.square(adjustedDecayTau - np.mean(adjustedDecayTau))
                        # rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
                        # if rSquared < 0.8:
                        #     # print("Fit quality poor; not shown")
                        #     continue
                        # print(f"R² = {rSquared}", "decay")
                        # y_dec = a * np.exp(b * ((x_dec/samplingFreqSec))) + c
                        # peakFig.plot(x_dec+int(width10[3][i[0][0]]), y_dec, color="C9")
                        # print("Tau (Decay):", abs(1/popt[1]))
                        continue
            except RuntimeError as e:
                if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                    continue
            except ValueError as e:
                if str(e) == "'x0' is infeasible":
                    continue 

    
    peakFig.hlines(*widthHalf[1:], color="C6")
    # peakFig.hlines(*widthBottom[1:], color="C7")
    peakFig.vlines(x=peaks, ymin=processedSignalArray[peaks] - peaksDict["prominences"] + (0.1 * peaksDict["prominences"]), ymax=processedSignalArray[peaks], color="C5")
    # peakFig.set_title(ratSide)
    peakFig.set_title("Individual Trace")
    fig.set_size_inches(9, 4.5)
    plt.axis([0, 47750, 0, 2])
    plt.yticks(fontsize="large")
    plt.xticks(ticks=seconds, labels=range(len(seconds)), fontsize="large")
    plt.xlabel("Time (s)", fontsize="x-large")
    plt.ylabel("Fluorescence (AU)", fontsize="x-large")
    fig.tight_layout()
    plt.minorticks_on()
    plt.show()
