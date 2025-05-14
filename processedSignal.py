import pyabf
import numpy as np
import scipy.signal as sci
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
        self.traceLength = len(processedSignalArray[0])
        self.date = self.abf.abfDateTime.strftime("%Y-%m-%d")

        self.samplingFreqMSec = self.abf.dataPointsPerMs + (1/3)
        self.samplingFreqSec = self.samplingFreqMSec * 1000
        self.seconds = np.arange(0, 47750, self.samplingFreqSec)

        self.wholeTraceAverages, self.wholeTraceMedians = avg.traceAverage(self.processedSignalArray)
        self.preInjAvg, self.preInjStdev = avg.preInjectionAverage(self.processedSignalArray, self.injTrace)
        self.normFluorescence = avg.deltaF(self.wholeTraceAverages, self.preInjAvg)
        self.overallZScore = avg.zCalc(self.wholeTraceAverages, self.processedSignalArray, self.injTrace)

    def peakMax(self):
        peakArray = np.zeros(np.size(self.processedSignalArray, 0))
        for i, t in enumerate(self.processedSignalArray):
            peaks, _ = sci.find_peaks(t, prominence= 0.05, wlen=20000)
            peakArray[i] = np.size(peaks)
        mostPeaksInTrace = int(np.max(peakArray))
        return mostPeaksInTrace