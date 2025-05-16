import peakAnalysis as pas
import numpy as np
import scipy.signal as sci
import pyabf

class TotalPeaks():
    def __init__(self, processedSignalArray, mainFile):
        self.processedSignalArray = processedSignalArray
        self.filename = mainFile
        self.abf = pyabf.ABF(mainFile)
        self.numTraces = len(self.abf.sweepList)

        self.mostPeaksInTrace = pas.peakMax(self.processedSignalArray)
        self.peaksArray, self.peaksDict = np.zeros((self.numTraces, self.mostPeaksInTrace)), {}
        self.widthBottomArray, self.width10Array, self.widthHalfArray, self.width90Array = (np.zeros((self.numTraces, 2, self.mostPeaksInTrace)) for _ in range(4))

        self.degreeNPeaks, self.decayNPeaks, self.riseNPeaks = {}, {}, {}
        self.overlapPeaks, self.overlapRise, self.overlapDecay = [[0]*self.mostPeaksInTrace for _ in self.numTraces], [[0]*self.mostPeaksInTrace for _ in self.numTraces], [[0]*self.mostPeaksInTrace for _ in self.numTraces]

    
    def peakWidths(self, processedsignalArray):
        for index, trace in enumerate(processedsignalArray):
            bottomWidth = sci.peak_widths(trace, self.peaksArray[index], rel_height=1, wlen=20000,
                                               prominence_data=(self.peaksDict[index]['prominences'], self.peaksDict[index]["left_bases"], 
                                                                self.peaksDict[index]["right_bases"]))
            width10 = sci.peak_widths(trace, self.peaksArray[index], rel_height=0.1, wlen=20000,
                                           prominence_data=(self.peaksDict[index]['prominences'], 
                                                            self.peaksDict[index]["left_bases"], self.peaksDict[index]["right_bases"]))
            widthHalf = sci.peak_widths(trace, self.peaksArray[index], rel_height=0.5, wlen=20000,
                                             prominence_data=(self.peaksDict[index]['prominences'], 
                                                              self.peaksDict[index]["left_bases"], self.peaksDict[index]["right_bases"]))
            width90 = sci.peak_widths(trace, self.peaksArray[index], rel_height=0.9, wlen=20000,
                                           prominence_data=(self.peaksDict[index]['prominences'], 
                                                            self.peaksDict[index]["left_bases"], self.peaksDict[index]["right_bases"]))
            for i, param in enumerate(bottomWidth[2:]):
                param = np.pad(param, pad_width= (0, self.mostPeaksInTrace - len(param)), mode= 'constant', constant_values= 0)
                self.widthBottomArray[index][i] = param
            for i, param in enumerate(width10[2:]):
                param = np.pad(param, pad_width= (0, self.mostPeaksInTrace - len(param)), mode= 'constant', constant_values= 0)
                self.width10Array[index][i] = param
            for i, param in enumerate(widthHalf[2:]):
                param = np.pad(param, pad_width= (0, self.mostPeaksInTrace - len(param)), mode= 'constant', constant_values= 0)
                self.widthHalfArray[index][i] = param
            for i, param in enumerate(width90[2:]):
                param = np.pad(param, pad_width= (0, self.mostPeaksInTrace - len(param)), mode= 'constant', constant_values= 0)
                self.width90Array[index][i] = param
    
    def peakFinder(self, processedSignalArray):
        for index, traces in enumerate(processedSignalArray): # Finds peaks in a signal
            peaks, self.peaksDict[index] = sci.find_peaks(traces, prominence= 0.05, width=0, wlen=20000, rel_height= 0.5)
            peaks = np.pad(peaks, pad_width= (0, self.mostPeaksInTrace - len(peaks)), mode= 'constant', constant_values= 0)
            bottomWidth, width10, width90 = np.array(bottomWidth), np.array(width10), np.array(width90)
            self.peaksArray[index] = peaks
            for i in self.peaksDict[index]:
                paddedEntry = np.pad(self.peaksDict[index][i], pad_width= (0, self.mostPeaksInTrace - len(self.peaksDict[index][i])), mode= 'constant', constant_values= 0)
                self.peaksDict[index][i] = paddedEntry

    def overlapCheck(self, peaksArray):
        degreeNPeaks, riseNPeaks, decayNPeaks = {}, {}, {}
        # Finds overlapping events by checking if the maxima of a peak is contained within the left and right slopes of another peak
        for k, checkPeak in np.ndenumerate(peaksArray):
            self.overlapPeaks[k[0]][k[1]] = np.array([y for y in peaksArray[k[0]] if ((self.widthBottomArray[k[0]][2][k[1]] < y < self.widthBottomArray[k[0]][3][k[1]]) and y != checkPeak)])
            self.overlapRise[k[0]][k[1]] = np.array([y for y in peaksArray[k[0]] if ((self.widthBottomArray[k[0]][2][k[1]] < y < checkPeak) and y != checkPeak)])
            self.overlapDecay[k[0]][k[1]] = np.array([y for y in peaksArray[k[0]] if ((checkPeak < y < self.widthBottomArray[k[0]][3][k[1]]) and y != checkPeak)])
        for z, x in enumerate(peaksArray):
            for i, peakOfDegree in enumerate(x):
                if len(self.processedSignalArray[z][int(self.widthBottomArray[z][2][i]):int(self.widthBottomArray[z][3][i])]) == 0:
                        continue
                degreeNPeaks[peakOfDegree] = len(self.overlapPeaks[z][i])
                riseNPeaks[peakOfDegree] = len(self.overlapRise[z][i])
                decayNPeaks[peakOfDegree] = len(self.overlapDecay[z][i])
            self.degreeNPeaks = dict(sorted(degreeNPeaks.items(), key=lambda item: item[1]))
            self.degreeRise = dict(sorted(riseNPeaks.items(), key=lambda item: item[1]))
            self.degreeDecay = dict(sorted(decayNPeaks.items(), key=lambda item: item[1]))