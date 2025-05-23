import totalPeaks as top
import numpy as np
import scipy.integrate as inte
import scipy.optimize as opt
import scipy.signal as sci

class TracePeaks(top.TotalPeaks):
    def __init__(self, processedSignalArray, mainFile, traceIndex):
        super().__init__(processedSignalArray, mainFile)

        self.processedSignalArray = processedSignalArray
        self.filename = mainFile
        self.traceIndex = traceIndex

        self.fullTraceArray = processedSignalArray[self.traceIndex]
        self.peaks = np.zeros(self.mostPeaksInTrace)
        self.traceDict = {}
        self.traceBottomWidths = np.zeros((2, self.mostPeaksInTrace))
        self.trace90Widths = np.zeros((2, self.mostPeaksInTrace))
        self.trace10Widths = np.zeros((2, self.mostPeaksInTrace))
        self.degreeNPeaks, self.decayNPeaks, self.riseNPeaks = {}, {}, {}

        self.peakNum = np.arange(1, len(self.peaks) + 1)
        self.peakIndex = self.peaks
        self.peakLocSec = np.zeros(self.mostPeaksInTrace)
        self.leftBounds = np.zeros(self.mostPeaksInTrace)
        self.rightBounds = np.zeros(self.mostPeaksInTrace)
        self.amplitude = np.zeros(self.mostPeaksInTrace)
        self.absoluteAmp = 0 #TODO: Create this
        self.width = np.zeros(self.mostPeaksInTrace)
        self.rightTail = np.zeros(self.mostPeaksInTrace)
        self.leftTail = np.zeros(self.mostPeaksInTrace)
        self.frequency = 0
        self.meanArea = np.zeros(self.mostPeaksInTrace)
        self.totArea = 0
        self.riseRate = np.zeros(self.mostPeaksInTrace, dtype=np.float64)
        self.decayRate = np.zeros(self.mostPeaksInTrace, dtype=np.float64)


    def peakFinder(self, singleTrace):
        peaks, traceDict = sci.find_peaks(singleTrace, prominence= 0.05, width=0, wlen=20000, rel_height= 0.5)
        bottomWidth = np.array(sci.peak_widths(singleTrace, peaks, rel_height=1, wlen=20000,
                                      prominence_data=(traceDict['prominences'], traceDict["left_bases"], 
                                                       traceDict["right_bases"])))
        width10 = np.array(sci.peak_widths(singleTrace, peaks, rel_height=0.1, wlen=20000,
                                  prominence_data=(traceDict['prominences'], 
                                                   traceDict["left_bases"], traceDict["right_bases"])))
        width90 = np.array(sci.peak_widths(singleTrace, peaks, rel_height=0.9, wlen=20000,
                                  prominence_data=(traceDict['prominences'], 
                                                   traceDict["left_bases"], traceDict["right_bases"])))
        self.peaks = np.pad(peaks, pad_width= (0, self.mostPeaksInTrace - len(peaks)), mode= 'constant', constant_values= 0)
        for i, param in enumerate(bottomWidth[2:]):
            param = np.pad(param, pad_width= (0, self.mostPeaksInTrace - len(param)), mode= 'constant', constant_values= 0)
            self.traceBottomWidths[i] = param
        for i, param in enumerate(width10[2:]):
            param = np.pad(param, pad_width= (0, self.mostPeaksInTrace - len(param)), mode= 'constant', constant_values= 0)
            self.trace10Widths[i] = param
        for i, param in enumerate(width90[2:]):
            param = np.pad(param, pad_width= (0, self.mostPeaksInTrace - len(param)), mode= 'constant', constant_values= 0)
            self.trace90Widths[i] = param
        for i in traceDict:
            paddedEntry = np.pad(traceDict[i], pad_width= (0, self.mostPeaksInTrace - len(traceDict[i])), mode= 'constant', constant_values= 0)
            self.traceDict[i] = paddedEntry

        self.peakNum = np.arange(1, len(self.peaks) + 1)
        self.peakIndex = self.peaks
        self.peakLocSec = ((self.peaks/self.samplingFreqSec) + (self.traceIndex * 30)).round(2)
        self.leftBounds = self.traceDict['left_ips'].round(2)
        self.rightBounds = self.traceDict['right_ips'].round(2)
        self.amplitude = self.traceDict['prominences'].round(3)
        self.absoluteAmp = 0 #TODO: Create this
        self.width = (self.traceDict['widths']/(self.samplingFreqMSec)).round(2)
        self.rightTail = ((self.traceDict['right_bases'] - self.peaks)/(self.samplingFreqMSec)).round(2)
        self.leftTail = ((self.peaks - self.traceDict['left_bases'])/(self.samplingFreqMSec)).round(2)
        
    def overlapCheck(self, peaksArray):
        # Finds overlapping events by checking if the maxima of a peak is contained within the left and right slopes of another peak
        for k, checkPeak in enumerate(peaksArray):
            self.overlapPeaks[k] = np.array([y for y in peaksArray if ((self.traceBottomWidths[0][k] < y < self.traceBottomWidths[1][k]) and y != checkPeak)])
            self.overlapRise[k] = np.array([y for y in peaksArray if ((self.traceBottomWidths[0][k] < y < checkPeak) and y != checkPeak)])
            self.overlapDecay[k] = np.array([y for y in peaksArray if ((checkPeak < y < self.traceBottomWidths[1][k]) and y != checkPeak)])
        for i, peakOfDegree in enumerate(peaksArray):
            if len(self.fullTraceArray[int(self.traceBottomWidths[0][i]):int(self.traceBottomWidths[1][i])]) == 0:
                    continue
            self.degreeNPeaks[peakOfDegree] = len(self.overlapPeaks[i])
            self.riseNPeaks[peakOfDegree] = len(self.overlapRise[i])
            self.decayNPeaks[peakOfDegree] = len(self.overlapDecay[i])
        self.degreeNPeaks = dict(sorted(self.degreeNPeaks.items(), key=lambda item: item[1]))
        self.degreeRise = dict(sorted(self.riseNPeaks.items(), key=lambda item: item[1]))
        self.degreeDecay = dict(sorted(self.decayNPeaks.items(), key=lambda item: item[1]))
    
    def freqSet(self):
        if np.size(self.peaks) == 0:
            self.frequency = 0
        else:
            self.frequency = round(np.count_nonzero(self.peaks)/((len(self.fullTraceArray) + 2250)/self.samplingFreqSec), 2) #Peaks/second (15 second trace)

    def areaSet(self):
        adjustedArea = np.zeros_like(self.peaks)
        for _, u in enumerate(self.degreeNPeaks):
            i = np.where(self.peaks == u)
            if len(self.fullTraceArray[int(self.traceBottomWidths[0][i[0][0]]):int(self.traceBottomWidths[1][i[0][0]])]) == 0:
                 continue
            peaksInArea = self.overlapPeaks[i[0][0]]
            peakArea = inte.simpson(y=self.fullTraceArray[int(self.traceBottomWidths[0][i[0][0]]):int(self.traceBottomWidths[1][i[0][0]])], 
                                    x=np.arange(int(self.traceBottomWidths[0][i[0][0]]), int(self.traceBottomWidths[1][i[0][0]]))/self.samplingFreqMSec)
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
            i = np.where(self.peaks == u)
            if len(self.fullTraceArray[int(self.traceBottomWidths[0][i[0][0]]):int(self.traceBottomWidths[1][i[0][0]])]) == 0:
                continue
            adjustedRiseTau = np.array(self.fullTraceArray[int(self.trace90Widths[0][i[0][0]]):int(u)])
            riseWidth = int(len(adjustedRiseTau))
            riseArray = np.arange(0, riseWidth)
            x_rise = np.linspace(np.min(riseArray), np.max(riseArray), 1000)
            p0 = (self.fullTraceArray[int(self.trace10Widths[0][i[0][0]])], 1, self.fullTraceArray[int(self.trace90Widths[0][i[0][0]])])
            if len(adjustedRiseTau) == 0:
                continue
            else:
                try:
                    if self.riseNPeaks[u] == 0: # Handles peaks with no overlap
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, riseArray/self.samplingFreqSec, adjustedRiseTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, 0, self.fullTraceArray[int(self.trace90Widths[0][i[0][0]])]], 
                                                                ub=[np.inf, np.inf, np.inf]),
                                                                maxfev=1000)
                        squaredDiffs = np.square(adjustedRiseTau - (popt[0] * np.exp(popt[1] * ((riseArray/self.samplingFreqSec))) + popt[2]))
                        squaredDiffsFromMean = np.square(adjustedRiseTau - np.mean(adjustedRiseTau))
                        rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
                        if rSquared < 0.8:
                            riseTauList[i[0][0]] = np.NaN
                        else:
                            riseTauList[i[0][0]] = abs(1/popt[1])
                            self.risePlot[i, 0] = x_rise+int(self.trace90Widths[0][i])
                            self.risePlot[i, 1] = y_rise
                            
                    else:
                        # adjustedRiseTau = np.array(self.fullTraceArray[int(self.traceHalfWidths[0][i[0][0]]):int(u)])
                        # p0 = (self.fullTraceArray[int(self.trace10Widths[0][i[0][0]])], 1, self.fullTraceArray[int(self.traceHalfWidths[0][i[0][0]])])
                        # for l in peaksInRise:
                        #     j = np.where(self.peaks == l)
                        #     r = np.where(adjustedRiseTau == self.fullTraceArray[int(self.traceBottomWidths[0][j[0][0]])])
                        #     o = np.where(adjustedRiseTau == self.fullTraceArray[int(self.traceBottomWidths[1][j[0][0]])])
                        #     if len(r[0]) == 0 or len(o[0]) == 0:
                        #         continue
                        #     for y, _ in enumerate(adjustedRiseTau[int(r[0][0]):int(o[0][0])]):
                        #         adjustedRiseTau[y+r[0][0]] = self.traceBottomWidths[1][j[0][0]] #TODO: Check why I need the width_heights array here
                        # riseWidth = int(len(adjustedRiseTau))
                        # riseArray = np.array(list(range(0, riseWidth)))
                        # popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp((b * t)) + c, riseArray/self.samplingFreqSec, adjustedRiseTau, p0=p0, 
                        #                         bounds=opt.Bounds(lb=[0, 0, self.fullTraceArray[int(self.traceBottomWidths[0][i[0][0]])]], 
                        #                                         ub=[np.inf, np.inf, np.inf]),
                        #                                         maxfev=1000)
                        # squaredDiffs = np.square(adjustedRiseTau - (popt[0] * np.exp(popt[1] * ((riseArray/self.samplingFreqSec))) + popt[2]))
                        # squaredDiffsFromMean = np.square(adjustedRiseTau - np.mean(adjustedRiseTau))
                        # rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
                        # if rSquared < 0.8:
                        #     riseTauList[i[0][0]] = np.NaN
                        # else:
                        #     riseTauList[i[0][0]] = abs(1/popt[1])
                        print("what's wrong?")
                        riseTauList[i[0][0]] = np.NaN

                except RuntimeError as e:
                    if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                        riseTauList[i[0][0]] = np.NaN 
                    else:
                        raise
                except ValueError as e:
                    if str(e) == "'x0' is infeasible":
                        riseTauList[i[0][0]] = np.NaN
                    else:
                        raise
        self.riseRate = riseTauList

    def decaySet(self):
        decayTauList = np.zeros(self.numTracePeaks)
        for _, u in enumerate(self.degreeDecay):
            i = np.where(self.peaks == u)
            adjustedDecayTau = np.array(self.fullTraceArray[int(u):int(self.trace90Widths[1][i[0][0]])])
            decWidth = int(len(adjustedDecayTau))
            decArray = np.arange(0, decWidth)
            x_dec = np.linspace(np.min(decArray), np.max(decArray), 1000)
            p0 = (self.fullTraceArray[int(self.trace10Widths[1][i[0][0]])], -1, self.fullTraceArray[int(self.trace90Widths[1][i[0][0]])])
            if len(adjustedDecayTau) == 0:
                continue
            else:
                try:
                    if self.decayNPeaks[u] == 0:
                        popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, decArray/self.samplingFreqSec, adjustedDecayTau, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, -np.inf, self.fullTraceArray[int(self.trace90Widths[1][i[0][0]])]], 
                                                                ub=[self.fullTraceArray[int(self.trace10Widths[1][i[0][0]])], 0, np.inf]),
                                                                maxfev=1000)
                        squaredDiffs = np.square(adjustedDecayTau - (popt[0] * np.exp(popt[1] * ((decArray/self.samplingFreqSec))) + popt[2]))
                        squaredDiffsFromMean = np.square(adjustedDecayTau - np.mean(adjustedDecayTau))
                        rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
                        if rSquared < 0.8:
                            decayTauList[i[0][0]] = np.NaN
                        else:
                            decayTauList[i[0][0]] = abs(1/popt[1])
                            self.decayPlot[i, 0] = x_dec+int(self.trace10Widths[1][i])
                            self.decayPlot[i, 1] = y_dec
                    else:
                        # adjustedDecayTau = np.array(self.fullTraceArray[int(u):int(self.traceHalfWidths[1][i[0][0]])])
                        # p0 = (self.fullTraceArray[int(self.trace10Widths[1][i[0][0]])], -1, self.fullTraceArray[int(self.traceHalfWidths[1][i[0][0]])])
                        # for l in peaksInDecay:
                        #     j = np.where(self.peaks == l)
                        #     r = np.where(adjustedDecayTau == self.fullTraceArray[int(self.traceBottomWidths[0][j[0][0]])])
                        #     o = np.where(adjustedDecayTau == self.fullTraceArray[int(self.traceBottomWidths[1][j[0][0]])])
                        #     if len(r[0]) == 0 or len(o[0]) == 0:
                        #         continue
                        #     for y, _ in enumerate(adjustedDecayTau[int(r[0][0]):int(o[0][0])]):
                        #         adjustedDecayTau[y+r[0][0]] = self.traceBottomWidths[1][j[0][0]] #TODO: Check why I need the width_heights array here
                        # decWidth = int(len(adjustedDecayTau))
                        # decArray = np.array(list(range(0, decWidth)))
                        # popt, _ = opt.curve_fit(lambda t, a, b, c: a * np.exp((b * t)) + c, decArray/self.samplingFreqSec, adjustedDecayTau, p0=p0, 
                        #                         bounds=opt.Bounds(lb=[0, -np.inf, self.fullTraceArray[int(self.traceHalfWidths[1][i[0][0]])]], 
                        #                                         ub=[self.fullTraceArray[int(self.trace10Widths[1][i[0][0]])], 0, np.inf]), 
                        #                                         maxfev=1000)
                        # squaredDiffs = np.square(adjustedDecayTau - (popt[0] * np.exp(popt[1] * ((decArray/self.samplingFreqSec))) + popt[2]))
                        # squaredDiffsFromMean = np.square(adjustedDecayTau - np.mean(adjustedDecayTau))
                        # rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
                        # if rSquared < 0.8:
                        #     decayTauList[i[0][0]] = np.NaN
                        # else:
                        #     decayTauList[i[0][0]] = abs(1/popt[1])
                        decayTauList[i[0][0]] = np.NaN
                except RuntimeError as e:
                    if str(e) == "Optimal parameters not found: The maximum number of function evaluations is exceeded.":
                        decayTauList[i[0][0]] = np.NaN
                except ValueError as e:
                    if str(e) == "'x0' is infeasible":
                        decayTauList[i[0][0]] = np.NaN 
        self.decayRate = decayTauList