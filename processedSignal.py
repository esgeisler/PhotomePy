import pyabf
import numpy as np
import scipy.signal as sci
import peakAnalysis as pas
import AverageTraces as avg

#TODO: Add Z-Score functionality for the individual traces.
class ProcessedTotalSignal():
    def __init__(self, filename, processedSignalArray, ratID, injTrace):
        self.processedSignalArray = processedSignalArray
        self.experimentFile = filename
        self.ratID = ratID
        self.injTrace = injTrace

        self.abf = pyabf.ABF(filename)
        self.numTraces = len(self.abf.sweepList)
        self.traceLength = len(self.processedSignalArray[0])
        self.date = self.abf.abfDateTime.strftime("%Y-%m-%d")

        self.samplingFreqMSec = self.abf.dataPointsPerMs + (1/3)
        self.samplingFreqSec = self.samplingFreqMSec * 1000
        self.seconds = np.arange(0, 47750, self.samplingFreqSec)

        self.wholeTraceAverages, self.wholeTraceMedians = avg.traceAverage(self.processedSignalArray)
        self.preInjAvg, self.preInjStdev = avg.preInjectionAverage(self.processedSignalArray, self.injTrace)
        self.normFluorescence = avg.deltaF(self.wholeTraceAverages, self.preInjAvg)
        self.overallZScore = avg.zCalc(self.wholeTraceAverages, self.processedSignalArray, self.injTrace)

    #     self.mostPeaksInTrace = pas.peakMax(self.processedSignalArray)
    #     self.peaksArray, self.peaksDict, self.widthBottomArray, self.width10Array, self.widthHalfArray, self.width90Array = pas.peakFinder(self.processedSignalArray, self.experimentFile)

    