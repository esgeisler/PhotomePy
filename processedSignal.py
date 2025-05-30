import pyabf
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

        self.wholeTraceAverages, self.wholeTraceMedians = avg.traceAverage(self.processedSignalArray)
        self.preInjAvg, self.preInjStdev = avg.preInjectionAverage(self.processedSignalArray, self.injTrace)
        self.normFluorescence = avg.deltaF(self.wholeTraceAverages, self.preInjAvg)
        self.overallZScore = avg.zCalc(self.wholeTraceAverages, self.processedSignalArray, self.injTrace)