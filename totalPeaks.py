import numpy as np
import pyabf
import peakAnalysis as pas
import yaml

class TotalPeaks():
    def __init__(self, processedSignalArray, mainFile, baselineStDev):
        self.processedSignalArray = processedSignalArray
        self.filename = mainFile
        self.abf = pyabf.ABF(mainFile)
        self.numTraces = self.abf.sweepCount
        
        with open("config.yaml") as c:
            userConfig = yaml.safe_load(c)
            start = userConfig["EVENT_HANDLING"]["trace_start_offset"]
            end = userConfig["EVENT_HANDLING"]["trace_end_offset"]
        removedPointCount = abs(start) + abs(end)

        self.samplingFreqMSec = self.abf.dataPointsPerMs + (1/3)
        self.samplingFreqSec = self.samplingFreqMSec * 1000
        self.seconds = np.arange(0, self.abf.sweepPointCount - removedPointCount, self.samplingFreqSec)

        self.mostPeaksInTrace, self.numPeaksArray = pas.peakMax(self.processedSignalArray, baselineStDev)
        self.peaksArray, self.peaksDict = np.zeros((self.numTraces, self.mostPeaksInTrace)), {}
        self.widthBottomArray, self.width10Array, self.width90Array = (np.zeros((self.numTraces, 2, self.mostPeaksInTrace)) for _ in range(3))