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
# def tauFit(tau, a, b):
#     return a*np.exp(-b * tau)

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
    widthHalfArray = np.zeros((traceLen, 4, longPeak))
    width90Array = np.zeros((traceLen, 4, longPeak))
    overlapPeaks, overlapRise, overlapDecay = [[0]*longPeak for _ in abf.sweepList], [[0]*longPeak for _ in abf.sweepList], [[0]*longPeak for _ in abf.sweepList]
    for index, traces in enumerate(processedSignalArray): # Finds peaks in a signal
        peaks, peaksDict[index] = sci.find_peaks(traces, prominence= 0.05, width=0, wlen=20000, rel_height= 0.5)
        bottomWidth = sci.peak_widths(traces, peaks, rel_height=1, prominence_data=(peaksDict[index]['prominences'], peaksDict[index]["left_bases"], 
                                                                                     peaksDict[index]["right_bases"]), wlen=20000)
        width10 = sci.peak_widths(traces, peaks, rel_height=0.1, prominence_data=(peaksDict[index]['prominences'], 
                                                                                  peaksDict[index]["left_bases"], peaksDict[index]["right_bases"]), wlen=20000)
        widthHalf = sci.peak_widths(traces, peaks, rel_height=0.5, prominence_data=(peaksDict[index]['prominences'], 
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
        for i, param in enumerate(widthHalf):
            param = np.pad(param, pad_width= (0, longPeak - len(param)), mode= 'constant', constant_values= 0)
            widthHalfArray[index][i] = param
        for i, param in enumerate(width90):
            param = np.pad(param, pad_width= (0, longPeak - len(param)), mode= 'constant', constant_values= 0)
            width90Array[index][i] = param
        for i in peaksDict[index]:
           paddedEntry = np.pad(peaksDict[index][i], pad_width= (0, longPeak - len(peaksDict[index][i])), mode= 'constant', constant_values= 0)
           peaksDict[index][i] = paddedEntry
    # Finds overlapping events by checking if the maxima of a peak is contained within the left and right slopes of another peak
    for k, checkPeak in np.ndenumerate(peaksArray):
            overlapPeaks[k[0]][k[1]] = np.array([y for y in peaksArray[k[0]] if ((widthBottomArray[k[0]][2][k[1]] < y < widthBottomArray[k[0]][3][k[1]]) and y != checkPeak)])
            overlapRise[k[0]][k[1]] = np.array([y for y in peaksArray[k[0]] if ((widthBottomArray[k[0]][2][k[1]] < y < checkPeak) and y != checkPeak)])
            overlapDecay[k[0]][k[1]] = np.array([y for y in peaksArray[k[0]] if ((checkPeak < y < widthBottomArray[k[0]][3][k[1]]) and y != checkPeak)])
    # Declares the dataframe where peak information will be stored
    # Following that, retrieves the simplest-to-calculate metadata: 
    # Event number, Event index and its loc. in seconds, the leftmost and rightmost end of the peak, the amplitude, and the width @ 50% of a peak's prominence.
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
        
        # Determines the "degree" of a peak by counting the amount of peaks found above in the checkpeak for loop for either the rise or decay.
        # Following this, sorts the peaks in order.
        for i, peakOfDegree in enumerate(x): 
            if len(processedSignalArray[z][int(widthBottomArray[z][2][i]):int(widthBottomArray[z][3][i])]) == 0:
                 continue
            degreeNPeaks[peakOfDegree] = len(overlapPeaks[z][i])
            riseNPeaks[peakOfDegree] = len(overlapRise[z][i])
            decayNPeaks[peakOfDegree] = len(overlapDecay[z][i])
        sortedDegreeNPeaks = dict(sorted(degreeNPeaks.items(), key=lambda item: item[1]))
        sortedDegreeRise = dict(sorted(riseNPeaks.items(), key=lambda item: item[1]))
        sortedDegreeDecay = dict(sorted(decayNPeaks.items(), key=lambda item: item[1]))
        # Generates the event area for each peak via the Simpson's approximation of the definite integral.
        # After this, corrects these event areas by subtracting the area of all of the smaller peaks that overlap with that event.
        for _, u in enumerate(sortedDegreeNPeaks):
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

        # Generates adjusted rise time constants (tau)
        # For non-overlapping peaks, they are fit with opt.curve_fit according to the exponential rise equation, with a range from 10-90% of a peak's relative height
        # By taking the reciprocal of the X value, the time constant (tau) is reported.
        # For peaks with overlap in their left slope, these are fit from the right-most base of an overlapping peak to 90% of a peak's height to reduce variability.
        # Any peaks with an R^2 value < 0.8 (80% correlation) are excluded from final reports.
        # If the initial guesses are invalid, or the fit takes >1000 samples, that peak is not calculated and excluded.
        # TODO: Implement changes from PeakDisplay; Update calculations to be  90-10 instead of 100-10
        for _, u in enumerate(sortedDegreeRise): 
            i = np.where(x == u)
            if len(processedSignalArray[z][int(widthBottomArray[z][2][i[0][0]]):int(widthBottomArray[z][3][i[0][0]])]) == 0:
                continue
            peaksInRise = overlapRise[z][i[0][0]]
            adjustedRiseTau = np.array(processedSignalArray[z][int(width90Array[z][2][i[0][0]]):int(u)])
            riseWidth = int(len(adjustedRiseTau))
            riseArray = np.array(list(range(0, riseWidth)))
            p0 = (processedSignalArray[z][int(width10Array[z][2][i[0][0]])], 1, processedSignalArray[z][int(width90Array[z][2][i[0][0]])])
            if len(adjustedRiseTau) == 0:
                continue
            else:
                try:
                    match riseNPeaks[u]:
                        case 0: # Handles peaks with no overlap
                            popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, riseArray/samplingFreqSec, adjustedRiseTau, p0=p0, 
                                                    bounds=opt.Bounds(lb=[0, 0, processedSignalArray[z][int(width90Array[z][2][i[0][0]])]], 
                                                                    ub=[np.inf, np.inf, np.inf]),
                                                                    maxfev=1000)
                            squaredDiffs = np.square(adjustedRiseTau - (popt[0] * np.exp(popt[1] * ((riseArray/samplingFreqSec))) + popt[2]))
                            squaredDiffsFromMean = np.square(adjustedRiseTau - np.mean(adjustedRiseTau))
                            rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
                            if rSquared < 0.8:
                                riseTauList[z][i[0][0]] = np.NaN
                            else:
                                riseTauList[z][i[0][0]] = abs(1/popt[1])
                        case _:
                            adjustedRiseTau = np.array(processedSignalArray[z][int(widthHalfArray[z][2][i[0][0]]):int(u)])
                            p0 = (processedSignalArray[z][int(width10Array[z][2][i[0][0]])], 1, processedSignalArray[z][int(widthHalfArray[z][2][i[0][0]])])
                            for l in peaksInRise:
                                j = np.where(x == l)
                                r = np.where(adjustedRiseTau == processedSignalArray[z][int(widthBottomArray[z][2][j[0][0]])])
                                o = np.where(adjustedRiseTau == processedSignalArray[z][int(widthBottomArray[z][3][j[0][0]])])
                                if len(r[0]) == 0 or len(o[0]) == 0:
                                    continue
                                for y, _ in enumerate(adjustedRiseTau[int(r[0][0]):int(o[0][0])]):
                                    adjustedRiseTau[y+r[0][0]] = widthBottomArray[z][1][j[0][0]]
                            riseWidth = int(len(adjustedRiseTau))
                            riseArray = np.array(list(range(0, riseWidth)))
                            popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp((b * t)) + c, riseArray/samplingFreqSec, adjustedRiseTau, p0=p0, 
                                                    bounds=opt.Bounds(lb=[0, 0, processedSignalArray[z][int(widthHalfArray[z][2][i[0][0]])]], 
                                                                    ub=[np.inf, np.inf, np.inf]),
                                                                    maxfev=1000)
                            squaredDiffs = np.square(adjustedRiseTau - (popt[0] * np.exp(popt[1] * ((riseArray/samplingFreqSec))) + popt[2]))
                            squaredDiffsFromMean = np.square(adjustedRiseTau - np.mean(adjustedRiseTau))
                            rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
                            if rSquared < 0.8:
                                riseTauList[z][i[0][0]] = np.NaN
                            else:
                                riseTauList[z][i[0][0]] = abs(1/popt[1])
                except RuntimeError as e:
                    if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                        riseTauList[z][i[0][0]] = np.NaN 
                except ValueError as e:
                    if str(e) == "'x0' is infeasible":
                        riseTauList[z][i[0][0]] = np.NaN 
        peakTable.Rise_Tau_exp = pd.Series(riseTauList[z])
        # Generated the adjusted decay time constants (tau)
        # This is calculated the same as the rise time, except that the equation is inverted and the curve fit to the right slope of an event.
        # Exclusions are the same as the rise calculations
        # TODO: Implement changes from PeakDisplay; Update calculations to be  90-10 instead of 100-10 
        for _, u in enumerate(sortedDegreeDecay):
            i = np.where(x == u)
            peaksInDecay = overlapDecay[z][i[0][0]]
            adjustedDecayTau = np.array(processedSignalArray[z][int(u):int(width90Array[z][3][i[0][0]])])
            decWidth = int(len(adjustedDecayTau))
            decArray = np.array(list(range(0, decWidth)))
            p0 = (processedSignalArray[z][int(width10Array[z][3][i[0][0]])], -1, processedSignalArray[z][int(width90Array[z][3][i[0][0]])])
            if len(adjustedDecayTau) == 0:
                continue
            else:
                try:
                    match decayNPeaks[u]:
                        case 0:
                            popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, decArray/samplingFreqSec, adjustedDecayTau, p0=p0, 
                                                    bounds=opt.Bounds(lb=[0, -np.inf, processedSignalArray[z][int(width90Array[z][3][i[0][0]])]], 
                                                                    ub=[processedSignalArray[z][int(width10Array[z][3][i[0][0]])], 0, np.inf]),
                                                                    maxfev=1000)
                            squaredDiffs = np.square(adjustedDecayTau - (popt[0] * np.exp(popt[1] * ((decArray/samplingFreqSec))) + popt[2]))
                            squaredDiffsFromMean = np.square(adjustedDecayTau - np.mean(adjustedDecayTau))
                            rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
                            if rSquared < 0.8:
                                decayTauList[z][i[0][0]] = np.NaN
                            else:
                                decayTauList[z][i[0][0]] = abs(1/popt[1])
                        case _:
                            adjustedDecayTau = np.array(processedSignalArray[z][int(u):int(widthHalfArray[z][3][i[0][0]])])
                            p0 = (processedSignalArray[z][int(width10Array[z][3][i[0][0]])], -1, processedSignalArray[z][int(widthHalfArray[z][3][i[0][0]])])
                            for l in peaksInDecay:
                                j = np.where(x == l)
                                r = np.where(adjustedDecayTau == processedSignalArray[z][int(widthBottomArray[z][2][j[0][0]])])
                                o = np.where(adjustedDecayTau == processedSignalArray[z][int(widthBottomArray[z][3][j[0][0]])])
                                if len(r[0]) == 0 or len(o[0]) == 0:
                                    continue
                                for y, _ in enumerate(adjustedDecayTau[int(r[0][0]):int(o[0][0])]):
                                    adjustedDecayTau[y+r[0][0]] = widthBottomArray[z][1][j[0][0]]
                            decWidth = int(len(adjustedDecayTau))
                            decArray = np.array(list(range(0, decWidth)))
                            popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp((b * t)) + c, decArray/samplingFreqSec, adjustedDecayTau, p0=p0, 
                                                    bounds=opt.Bounds(lb=[0, -np.inf, processedSignalArray[z][int(widthHalfArray[z][3][i[0][0]])]], 
                                                                    ub=[processedSignalArray[z][int(width10Array[z][3][i[0][0]])], 0, np.inf]), 
                                                                    maxfev=1000)
                            squaredDiffs = np.square(adjustedDecayTau - (popt[0] * np.exp(popt[1] * ((decArray/samplingFreqSec))) + popt[2]))
                            squaredDiffsFromMean = np.square(adjustedDecayTau - np.mean(adjustedDecayTau))
                            rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
                            if rSquared < 0.8:
                                decayTauList[z][i[0][0]] = np.NaN
                            else:
                                decayTauList[z][i[0][0]] = abs(1/popt[1])
                except RuntimeError as e:
                    if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                        decayTauList[z][i[0][0]] = np.NaN
                except ValueError as e:
                    if str(e) == "'x0' is infeasible":
                        riseTauList[z][i[0][0]] = np.NaN 
        peakTable.Decay_Tau_exp = pd.Series(decayTauList[z])
        # Determines frequency and total area of a sweep
        # Frequency is calculated by taking the # of sweeps in a trace and dividing that by the length of the trace in seconds.
        # Total Area is a sum of all of the (already-adjusted for overlapping smaller peaks) events in a given trace.
        match np.size(x): #
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
        peakFig.annotate("Trough for Peak %i"%(i+1), xycoords= 'data', size= 8, horizontalalignment= 'center',
                     xytext = (int(width90[3][i]), processedSignalArray[int(width90[3][i])] - 0.3), 
                     xy = (int(width90[3][i]), processedSignalArray[int(width90[3][i])] - 0.01),
                     arrowprops=dict(facecolor= 'black', width= 1, headwidth= 5, headlength= 5))
        peakFig.annotate("Peak %i"%(i+1), xycoords= 'data', size= 8, horizontalalignment= 'center',
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
        p0 = (processedSignalArray[int(width90[2][i[0][0]])], 1, processedSignalArray[int(width10[2][i[0][0]])])
        if len(adjustedRiseTau) == 0:
            continue
        else:
            try:
                match riseNPeaks[u]:
                    case 0: # Handles peaks with no overlap
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, riseArray/samplingFreqSec, adjustedRiseTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, 0, processedSignalArray[int(width90[2][i[0][0]])]], 
                                                                ub=[np.inf, np.inf, processedSignalArray[int(width10[2][i[0][0]])]]),
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
                        ji = np.where(peaks == peaksInRise[0])
                        adjustedRiseTau = np.array(processedSignalArray[int(widthBottom[3][ji[0][0]]):int(width10[2][i[0][0]])])
                        p0 = (processedSignalArray[int(widthBottom[3][ji[0][0]])], 1, processedSignalArray[int(width10[2][i[0][0]])])
                        for l in peaksInRise:
                            j = np.where(peaks == l)
                            r = np.where(adjustedRiseTau == processedSignalArray[int(width90[2][j[0][0]])])
                            o = np.where(adjustedRiseTau == processedSignalArray[int(width90[3][j[0][0]])])
                            if len(r[0]) == 0 or len(o[0]) == 0:
                                continue
                            for y, _ in enumerate(adjustedRiseTau[int(r[0][0]):int(o[0][0])]):
                                adjustedRiseTau[y+r[0][0]] = width90[1][j[0][0]]
                        riseWidth = int(len(adjustedRiseTau))
                        riseArray = np.array(list(range(0, riseWidth)))
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp((b * t)) + c, riseArray/samplingFreqSec, adjustedRiseTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, 0, processedSignalArray[int(widthBottom[3][ji[0][0]])]], 
                                                                ub=[np.inf, np.inf, processedSignalArray[int(width10[2][i[0][0]])]]),
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
                        peakFig.plot(x_rise+widthBottom[3][ji], y_rise, color="C8")
                        print("Tau (Rise):", abs(1/popt[1]))
            except RuntimeError as e:
                if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                    continue 
            except ValueError as e:
                if str(e) == "'x0' is infeasible":
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
                        ji = np.where(peaks == peaksInDecay[0])
                        adjustedDecayTau = np.array(processedSignalArray[int(width10[3][i[0][0]]):int(width90[2][ji[0][0]])])
                        p0 = (processedSignalArray[int(width10[3][i[0][0]])], -1, processedSignalArray[int(width90[2][ji[0][0]])])
                        for l in peaksInDecay:
                            j = np.where(peaks == l)
                            r = np.where(adjustedDecayTau == processedSignalArray[int(width90[2][j[0][0]])])
                            o = np.where(adjustedDecayTau == processedSignalArray[int(width90[3][j[0][0]])])
                            if len(r[0]) == 0 or len(o[0]) == 0:
                                continue
                            for y, _ in enumerate(adjustedDecayTau[int(r[0][0]):int(o[0][0])]):
                                adjustedDecayTau[y+r[0][0]] = width90[1][j[0][0]]
                        decWidth = int(len(adjustedDecayTau))
                        decArray = np.array(list(range(0, decWidth)))
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp((b * t)) + c, decArray/samplingFreqSec, adjustedDecayTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, -np.inf, processedSignalArray[int(widthBottom[2][ji[0][0]])]], 
                                                                ub=[processedSignalArray[int(width10[3][i[0][0]])], 0, np.inf]),
                                                                maxfev=1000)
                        x_dec = np.linspace(np.min(decArray), np.max(decArray), 1000)
                        a, b, c = popt[0], popt[1], popt[2]
                        squaredDiffs = np.square(adjustedDecayTau - (a * np.exp(b * ((decArray/samplingFreqSec))) + c))
                        squaredDiffsFromMean = np.square(adjustedDecayTau - np.mean(adjustedDecayTau))
                        rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
                        if rSquared < 0.8:
                            # print("Fit quality poor; not shown")
                            continue
                        print(f"R² = {rSquared}", "decay")
                        y_dec = a * np.exp(b * ((x_dec/samplingFreqSec))) + c
                        peakFig.plot(x_dec+int(width10[3][i[0][0]]), y_dec, color="C9")
                        print("Tau (Decay):", abs(1/popt[1]))
            except RuntimeError as e:
                if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                    continue
            except ValueError as e:
                if str(e) == "'x0' is infeasible":
                    continue 

    
    peakFig.hlines(*widthHalf[1:], color="C6")
    # peakFig.hlines(*widthBottom[1:], color="C7")
    peakFig.vlines(x=peaks, ymin=processedSignalArray[peaks] - peaksDict["prominences"] + (0.1 * peaksDict["prominences"]), ymax=processedSignalArray[peaks], color="C5")
    peakFig.set_title(ratSide)
    plt.axis([0, 47750, 0, 2])
    plt.xticks(ticks=seconds, labels=range(len(seconds)))
    plt.xlabel("Time (s)")
    plt.ylabel("Fluorescence (AU)")
    plt.minorticks_on()
    plt.show()
