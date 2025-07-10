import totalPeaks as top
import numpy as np
import scipy.integrate as inte
import scipy.optimize as opt
import scipy.signal as sci
import yaml

class TracePeaks(top.TotalPeaks):
    def __init__(self, processedSignalArray, mainFile, traceIndex):
        super().__init__(processedSignalArray, mainFile)

        self.processedSignalArray = processedSignalArray
        self.filename = mainFile
        self.traceIndex = traceIndex

        self.numTracePeaks = int(self.numPeaksArray[traceIndex])
        self.fullTraceArray = processedSignalArray[self.traceIndex]
        self.peaks = self.numPeaksArray[traceIndex]
        self.traceDict = {}
        self.traceBottomWidths = np.zeros((2, self.numTracePeaks))
        self.trace90Widths = np.zeros((2, self.numTracePeaks)) # Height at 90% of the relative height (bottom of the peak)
        self.traceHalfWidths = np.zeros((2, self.numTracePeaks))
        self.trace10Widths = np.zeros((2, self.numTracePeaks)) # Height at 10% of the relative height (top of the peak)
        self.degreeNPeaks, self.decayNPeaks, self.riseNPeaks = {}, {}, {} # Key-value pairs of peaks and the number of peaks contained within them, their decay slope, or their rise slope.
        self.overlapPeaks, self.overlapRise, self.overlapDecay = {}, {}, {} # Key-value pairs of peaks and all of the peaks contained within them, their decay, or their rise.
        self.firstTrainPeak, self.lastTrainPeak = {}, {} #Key-value pairs of the left-most or right-most peaks that overlap the rise or decay of a parent peak.

        self.peakNum = np.arange(1, self.numTracePeaks + 1)
        self.peakIndex = self.peaks
        self.peakLocSec = np.zeros(self.numTracePeaks)
        self.leftBounds, self.rightBounds = np.zeros(self.numTracePeaks), np.zeros(self.numTracePeaks)
        self.amplitude = np.zeros(self.numTracePeaks)
        self.absoluteAmp = np.zeros(self.numTracePeaks)
        self.width = np.zeros(self.numTracePeaks)
        self.rightTail = np.zeros(self.numTracePeaks)
        self.leftTail = np.zeros(self.numTracePeaks)
        self.frequency = 0
        self.meanArea = np.zeros(self.numTracePeaks)
        self.totArea = 0
        self.risePlot = np.zeros((self.numTracePeaks, 2, 1000), dtype=np.float64)
        self.riseRate = np.zeros(self.numTracePeaks, dtype=np.float64)
        self.decayPlot = np.zeros((self.numTracePeaks, 2, 1000), dtype=np.float64)
        self.decayRate = np.zeros(self.numTracePeaks, dtype=np.float64)

        with open("config.yaml") as c:
            userConfig = yaml.safe_load(c)
        self.norm_method = userConfig["GENERAL"]["normalization"]
        self.peakMethod = userConfig["EVENT_HANDLING"]["peak_id_method"]
        self.peakThreshold = userConfig["EVENT_HANDLING"]["peak_detection_threshold"]
        self.peakWindow = userConfig["EVENT_HANDLING"]["peak_window"]
        self.peakTop = userConfig["EVENT_HANDLING"]["peak_top"]
        self.peakBase = userConfig["EVENT_HANDLING"]["peak_bottom"]

# Provides the R­­² value to determine the fit of a curve. Currently used in the rise/decay τ calculation to exclude poorly fitted exponential curves.
    def rSquaredGet(self, adjustedTau, slopeArray, a, b, c):
        squaredDiffs = np.square(adjustedTau - (a * np.exp(b * ((slopeArray/self.samplingFreqSec))) + c))
        squaredDiffsFromMean = np.square(adjustedTau - np.mean(adjustedTau))
        rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
        return rSquared

# Locates peaks in a dataset based off of their prominence and the configuration in config.yaml.
# Input: A trace of any length containing continuous photometry signals.
# Output: Updates the tracePeaks peakNum, peakIndex, peakLocSec, leftBounds, rightBounds, amplitude, width, right tail and left tail to match the values from imported dataset.
    def peakFinder(self, singleTrace):
        self.peaks, traceDict = sci.find_peaks(singleTrace, prominence=self.peakThreshold, width=0, wlen=self.peakWindow, rel_height= 0.5)
        width10 = np.array(sci.peak_widths(singleTrace, self.peaks, rel_height=self.peakTop, wlen=self.peakWindow,
                                  prominence_data=(traceDict['prominences'], 
                                                   traceDict["left_bases"], traceDict["right_bases"])))
        widthHalf = np.array(sci.peak_widths(singleTrace, self.peaks, rel_height=0.5, wlen=self.peakWindow,
                                  prominence_data=(traceDict['prominences'], 
                                                   traceDict["left_bases"], traceDict["right_bases"])))
        width90 = np.array(sci.peak_widths(singleTrace, self.peaks, rel_height=self.peakBase, wlen=self.peakWindow,
                                  prominence_data=(traceDict['prominences'], 
                                                   traceDict["left_bases"], traceDict["right_bases"])))
        for i, param in enumerate(width10[2:]):
            self.trace10Widths[i] = param
        for i, param in enumerate(widthHalf[2:]):
            self.traceHalfWidths[i] = param
        for i, param in enumerate(width90[2:]):
            self.trace90Widths[i] = param
        self.traceDict = traceDict

        self.peakNum = np.arange(1, len(self.peaks) + 1)
        self.peakIndex = self.peaks
        self.peakLocSec = ((self.peaks/self.samplingFreqSec) + (self.traceIndex * 30)).round(2)
        self.leftBounds = self.traceDict['left_ips'].round(2)
        self.rightBounds = self.traceDict['right_ips'].round(2)
        self.amplitude = (self.traceDict['prominences'] - (self.peakTop * self.traceDict["prominences"])).round(3)
        self.width = (self.traceDict['widths']/(self.samplingFreqMSec)).round(2)
        self.rightTail = ((self.traceDict['right_bases'] - self.peaks)/(self.samplingFreqMSec)).round(2)
        self.leftTail = ((self.peaks - self.traceDict['left_bases'])/(self.samplingFreqMSec)).round(2)

    def peaksFromZero(self, singleTrace):
        peakCutoff = self.peakThreshold
        for points in singleTrace:
            aboveThresholdSignals = []
            if points >  peakCutoff:
                aboveThresholdSignals.append(points)
        
        return

# Finds overlapping events by checking if the maxima of a peak is contained within the left and right slopes of another peak
# Input: An array of traces containing fiber photometry daya
# Output: Sets values for overlapPeaks, overlapRise, overlapDecay, degreeNPeaks, riseNPeaks, decayNPeaks, degreeRise, and degreeDecay
    def overlapCheck(self, peaksArray):
        for i, peak in enumerate(peaksArray):
            if len(self.fullTraceArray[int(self.trace90Widths[0][i]):int(self.trace90Widths[1][i])]) == 0:
                    continue
            self.overlapPeaks[peak] = [y for y in peaksArray if ((self.trace90Widths[0][i] < y < self.trace90Widths[1][i]) and y != peak)]
            self.overlapRise[peak] = [y for y in peaksArray if ((self.trace90Widths[0][i] < y < peak) and y != peak)]
            self.overlapDecay[peak] = [y for y in peaksArray if ((peak < y < self.trace90Widths[1][i]) and y != peak)]
            self.degreeNPeaks[peak] = len(self.overlapPeaks[peak])
            self.riseNPeaks[peak] = len(self.overlapRise[peak])
            self.decayNPeaks[peak] = len(self.overlapDecay[peak])
        self.degreeNPeaks = dict(sorted(self.degreeNPeaks.items(), key=lambda item: item[1]))
        self.degreeRise = dict(sorted(self.riseNPeaks.items(), key=lambda item: item[1]))
        self.degreeDecay = dict(sorted(self.decayNPeaks.items(), key=lambda item: item[1]))

    def overlapAmplitude(self):
        adjustedAmp = np.zeros(self.numTracePeaks)
        sortedOverlapPeaks = dict(sorted(self.overlapPeaks.items(), key=lambda item: item[1]))
        for i in sortedOverlapPeaks:
            j = np.where(self.peaks == i)
            if not self.overlapPeaks[i]:
                adjustedAmp[j[0][0]] = self.amplitude[j[0][0]].round(3)
            else:
                adjustedAmp[j[0][0]] = self.amplitude[j[0][0]].round(3)
                for x in self.overlapPeaks[i]:
                    h = np.where(self.peaks == x)
                    parentPeak, parentProminence  = self.fullTraceArray[self.peaks][j[0][0]], self.traceDict["prominences"][j[0][0]]
                    parentBase = parentPeak - parentProminence + (self.peakTop * parentProminence)

                    childPeak, childProminence = self.fullTraceArray[self.peaks][h[0][0]], self.traceDict["prominences"][h[0][0]]
                    childBase = childPeak - childProminence + (self.peakTop * childProminence)

                    ampAdjustmentFactor = childBase - parentBase
                    adjustedAmp[h[0][0]] = (childProminence + ampAdjustmentFactor - (self.peakTop * self.traceDict["prominences"][h[0][0]])).round(3)
        self.absoluteAmp = adjustedAmp
            
    def freqSet(self):
        if np.size(self.peaks) == 0:
            self.frequency = 0
        else:
            self.frequency = round(np.count_nonzero(self.peaks)/((len(self.fullTraceArray) + 2250)/self.samplingFreqSec), 2) #Peaks/second (15 second trace)

    def areaSet(self): 
        adjustedArea = np.zeros(self.numTracePeaks)
        for _, u in enumerate(self.degreeNPeaks):
            i = np.where(self.peaks == u)
            if len(self.fullTraceArray[int(self.trace90Widths[0][i[0][0]]):int(self.trace90Widths[1][i[0][0]])]) == 0:
                continue
            peaksInArea = self.overlapPeaks[u] #TODO: Check if updated format of self.overlapPeaks alters the values here
            match self.norm_method:
                case "zSBR":
                    peakArea = inte.simpson(y=self.fullTraceArray[int(self.trace90Widths[0][i[0][0]]):int(self.trace90Widths[1][i[0][0]])]+100, 
                                                x=np.arange(int(self.trace90Widths[0][i[0][0]]), int(self.trace90Widths[1][i[0][0]]))/self.samplingFreqMSec)
                    w = len(np.arange(int(self.trace90Widths[0][i[0][0]]), int(self.trace90Widths[1][i[0][0]])))/self.samplingFreqMSec
                    L = self.fullTraceArray[self.peaks][i[0][0]] + 100 - self.amplitude[i[0][0]]
                    peakArea -= L * w
                case "SBR":
                    peakArea = inte.simpson(y=self.fullTraceArray[int(self.trace90Widths[0][i[0][0]]):int(self.trace90Widths[1][i[0][0]])], 
                                            x=np.arange(int(self.trace90Widths[0][i[0][0]]), int(self.trace90Widths[1][i[0][0]]))/self.samplingFreqMSec)
                    w = (self.trace90Widths[1][i[0][0]] - self.trace90Widths[0][i[0][0]]) /self.samplingFreqMSec
                    L = self.fullTraceArray[self.peaks][i[0][0]] - self.amplitude[i[0][0]]
                    peakArea -= L * w
            if self.degreeNPeaks[u] == 0:
                adjustedArea[i[0][0]] = peakArea.round(2)
            else:
                for l in peaksInArea:
                    j = np.where(self.peaks == l)
                    peakArea -= adjustedArea[j[0][0]]
                adjustedArea[i[0][0]] = peakArea.round(2)
        self.meanArea = adjustedArea
        if np.size(self.peaks) == 0:
            self.totArea = 0
        else:
            self.totArea = sum(adjustedArea)    

    def riseSet(self):
        riseTauList = np.zeros(self.numTracePeaks)
        for _, u in enumerate(self.degreeRise): 
            i = np.where(self.peaks == u)[0][0]
            if len(self.fullTraceArray[int(self.trace90Widths[0][i]):int(self.trace90Widths[1][i])]) == 0:
                continue
            adjustedRiseTau = np.array(self.fullTraceArray[int(self.trace90Widths[0][i]):int(self.trace10Widths[0][i])])
            if len(adjustedRiseTau) == 0:
                continue
            riseWidth = int(len(adjustedRiseTau))
            riseArray = np.arange(0, riseWidth)
            x_rise = np.linspace(np.min(riseArray), np.max(riseArray), 1000)
            p0 = (self.fullTraceArray[int(self.trace10Widths[0][i])], 1, self.fullTraceArray[int(self.trace90Widths[0][i])])
            if len(adjustedRiseTau) == 0:
                continue
            else:
                try:
                    if self.riseNPeaks[u] == 0: # Handles peaks with no overlap
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, riseArray/self.samplingFreqSec, adjustedRiseTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, 0, self.fullTraceArray[int(self.trace90Widths[0][i])]], 
                                                                ub=[np.inf, np.inf, np.inf]),
                                                                maxfev=1000)
                        rSquared = self.rSquaredGet(adjustedRiseTau, riseArray, popt[0], popt[1], popt[2])   
                        y_rise = popt[0] * np.exp(popt[1] * (x_rise/self.samplingFreqSec)) + popt[2]
                        
                        if rSquared < 0.8:
                            riseTauList[i] = np.NaN
                        else:
                            riseTauList[i] = abs(1/popt[1])
                            self.risePlot[i, 0] = x_rise+int(self.trace90Widths[0][i])
                            self.risePlot[i, 1] = y_rise
                    elif self.riseNPeaks[u] > 0: # Handles peaks with overlap, which are the first overlapping peaks
                        h = np.where(self.peaks == self.overlapRise[u][0])[0][0]
                        adjustedRiseTau = np.array(self.fullTraceArray[int(self.trace90Widths[0][i]):int(self.trace10Widths[0][h])])
                        p0 = (self.fullTraceArray[int(self.trace10Widths[0][h])], 1, self.fullTraceArray[int(self.trace90Widths[0][i])])
                        riseWidth = int(len(adjustedRiseTau))
                        riseArray = np.array(list(range(0, riseWidth)))
                        x_rise = np.linspace(np.min(riseArray), np.max(riseArray), 1000)
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp((b * t)) + c, riseArray/self.samplingFreqSec, adjustedRiseTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, 0, self.fullTraceArray[int(self.trace90Widths[0][i])]], 
                                                                ub=[np.inf, np.inf, np.inf]),
                                                                maxfev=1000)
                        rSquared = self.rSquaredGet(adjustedRiseTau, riseArray, popt[0], popt[1], popt[2])   
                        y_rise = popt[0] * np.exp(popt[1] * (x_rise/self.samplingFreqSec)) + popt[2]
                        riseTauList[i] = np.NaN
                        if rSquared < 0.8:
                            riseTauList[h] = np.NaN
                        else:
                            riseTauList[h] = abs(1/popt[1])
                            self.risePlot[h, 0] = x_rise+int(self.trace90Widths[0][i])
                            self.risePlot[h, 1] = y_rise
                        
                    else:
                        riseTauList[i] = np.NaN
                except RuntimeError as e:
                    if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                        riseTauList[i] = np.NaN 
                    else:
                        raise
                except ValueError as e:
                    if str(e) == "'x0' is infeasible.":
                        riseTauList[i] = np.NaN
        self.riseRate = riseTauList

    def decaySet(self):
        decayTauList = np.zeros(self.numTracePeaks)
        for _, u in enumerate(self.degreeDecay):
            i = np.where(self.peaks == u)[0][0]
            adjustedDecayTau = np.array(self.fullTraceArray[int(self.trace10Widths[1][i]):int(self.trace90Widths[1][i])])
            if len(adjustedDecayTau) == 0:
                continue
            decWidth = int(len(adjustedDecayTau))
            decArray = np.arange(0, decWidth)
            x_dec = np.linspace(np.min(decArray), np.max(decArray), 1000)
            p0 = (self.fullTraceArray[int(self.trace10Widths[1][i])], -1, self.fullTraceArray[int(self.trace90Widths[1][i])])
            if len(adjustedDecayTau) == 0:
                continue
            else:
                try:
                    if self.decayNPeaks[u] == 0:
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, decArray/self.samplingFreqSec, adjustedDecayTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, -np.inf, self.fullTraceArray[int(self.trace90Widths[1][i])]], 
                                                                ub=[self.fullTraceArray[int(self.trace10Widths[1][i])], 0, np.inf]),
                                                                maxfev=1000)
                        rSquared = self.rSquaredGet(adjustedDecayTau, decArray, popt[0], popt[1], popt[2])
                        y_dec = popt[0] * np.exp(popt[1] * (x_dec/self.samplingFreqSec)) + popt[2]
                        if rSquared < 0.8:
                            decayTauList[i] = np.NaN
                        else:
                            decayTauList[i] = abs(1/popt[1])
                            self.decayPlot[i, 0] = x_dec+int(self.trace10Widths[1][i])
                            self.decayPlot[i, 1] = y_dec
                    elif self.decayNPeaks[u] > 0:
                        h = np.where(self.peaks == self.overlapDecay[u][-1])[0][0]
                        adjustedDecayTau = np.array(self.fullTraceArray[int(self.trace10Widths[1][h]):int(self.trace90Widths[1][i])])
                        p0 = (self.fullTraceArray[int(self.trace10Widths[1][h])], -1, self.fullTraceArray[int(self.trace90Widths[1][i])])
                        decWidth = int(len(adjustedDecayTau))
                        decArray = np.array(list(range(0, decWidth)))
                        x_dec = np.linspace(np.min(decArray), np.max(decArray), 1000)
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp((b * t)) + c, decArray/self.samplingFreqSec, adjustedDecayTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, -np.inf, self.fullTraceArray[int(self.trace90Widths[1][i])]], 
                                                                ub=[self.fullTraceArray[int(self.trace10Widths[1][h])], 0, np.inf]), 
                                                                maxfev=1000)
                        rSquared = self.rSquaredGet(adjustedDecayTau, decArray, popt[0], popt[1], popt[2])
                        y_dec = popt[0] * np.exp(popt[1] * (x_dec/self.samplingFreqSec)) + popt[2]
                        decayTauList[i] = np.NaN
                        if rSquared < 0.8:
                            decayTauList[h] = np.NaN
                        else:
                            decayTauList[h] = abs(1/popt[1])
                            self.decayPlot[h, 0] = x_dec+int(self.trace10Widths[1][h])
                            self.decayPlot[h, 1] = y_dec
                    else:
                        decayTauList[i] = np.NaN
                except RuntimeError as e:
                    if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                        decayTauList[i] = np.NaN
                    else:
                        raise
                except ValueError as e:
                    if str(e) == "'x0' is infeasible.":
                        decayTauList[i] = np.NaN 
        self.decayRate = decayTauList