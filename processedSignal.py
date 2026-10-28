import pyabf
import numpy as np

class ProcessedTotalSignal():
    def __init__(self, filename, processedSignalArray, ratID, injTrace):
        self.processedSignalArray = processedSignalArray
        self.experimentFile = filename
        self.ratID = ratID
        self.injTrace = injTrace

        if self.injTrace <= 1:
            raise ValueError("Injection traces must be larger than 1")
        elif self.injTrace > len(processedSignalArray):
            raise IndexError("Injection traces are larger than total sweeps in signals")

        self.abf = pyabf.ABF(filename)
        self.numTraces = self.abf.sweepCount
        self.traceLength = self.abf.sweepPointCount
        self.date = self.abf.abfDateTime.strftime("%Y-%m-%d")

        

        self.wholeTraceAverages, self.wholeTraceMedians = np.mean(self.processedSignalArray, axis=1, dtype=np.float64), np.median(self.processedSignalArray, axis=1)
        self.preInjAvg = np.mean(processedSignalArray[:self.injTrace + 1], dtype=np.float64)
        self.preInjStdev = np.std(processedSignalArray[:self.injTrace + 1], dtype=np.float64, mean=self.preInjAvg, ddof=1)
        self.normFluorescence = (self.wholeTraceAverages - self.preInjAvg)/self.preInjAvg
        self.overallZScore = (self.wholeTraceAverages - self.preInjAvg)/self.preInjStdev