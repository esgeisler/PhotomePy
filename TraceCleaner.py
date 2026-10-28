import yaml
import pandas as pd
import os
import pyabf
import AutoCleaner as acl
import numpy as np
import csv
from pathlib import Path

class TraceCleaner():
    def __init__(self, mainFile, baselineFile):
        if not os.path.exists(mainFile):
            raise FileNotFoundError("main")
        elif not os.path.exists(baselineFile):
            raise FileNotFoundError("baseline")
        self.baselineFile = baselineFile
        self.mainFile = mainFile
        
        with open("config.yaml") as c:
            userConfig = yaml.safe_load(c)
        self.leftChannels = userConfig["GENERAL"]["left_rat_channels"]
        self.rightChannels = userConfig["GENERAL"]["right_rat_channels"]
        self.start = userConfig["EVENT_HANDLING"]["trace_start_offset"]
        self.end = userConfig["EVENT_HANDLING"]["trace_end_offset"]
        self.removedPointCount = abs(self.start) + abs(self.end)
        
        self.abf = pyabf.ABF(self.mainFile)
        self.samplingFreqMSec = self.abf.dataPointsPerMs + (1/3)
        self.samplingFreqSec = self.samplingFreqMSec * 1000

        self.baseline470Left, self.baseline405Left, self.baseline470Right, self.baseline405Right= acl.baselineGet(self.baselineFile)
        self.subtract470Left, self.subtract405Left = acl.baselineSubtractor(self.mainFile, self.baseline470Left, self.baseline405Left, self.leftChannels)
        self.subtract470Right, self.subtract405Right = acl.baselineSubtractor(self.mainFile, self.baseline470Right, self.baseline405Right, self.rightChannels)
        self.filtered405Left, self.filtered405Right = acl.wholeTraceGauss(acl.wholeTraceMedFilt(self.subtract405Left)), acl.wholeTraceGauss(acl.wholeTraceMedFilt(self.subtract405Right))
        self.filtered470Left, self.filtered470Right = acl.wholeTraceMedFilt(self.subtract470Left), acl.wholeTraceMedFilt(self.subtract470Right)

        self.ratioSignalLeft, self.ratioSignalRight = self.filtered470Left/self.filtered405Left, self.filtered470Right/self.filtered405Right
        self.signalLeft, self.signalRight = acl.wholeTraceGauss(self.ratioSignalLeft), acl.wholeTraceGauss(self.ratioSignalRight)

    def completeProcessor(self):
        return self.signalLeft, self.signalRight

# Exact same function as completeProcessor, except that it returns signals at each step of the process.
    def stepwiseProcessor(self):
        return (self.subtract470Left, self.subtract405Left), (self.subtract470Right, self.subtract405Right), (self.filtered470Left, self.filtered405Left), (self.filtered470Right, self.filtered405Right), self.ratioSignalLeft, self.ratioSignalRight, self.signalLeft, self.signalRight
        
# Processor for refReg. 
    def referenceSignalDecayProcessor(self):
        leftFittedArray = np.zeros((self.abf.sweepCount, self.abf.sweepPointCount - self.removedPointCount))
        rightFittedArray = np.zeros((self.abf.sweepCount, self.abf.sweepPointCount - self.removedPointCount))

        # ratioLeft = self.ratioSignalLeft.mean(axis=0)
        # ratioRight = self.ratioSignalRight.mean(axis=0)
        for i in self.abf.sweepList:
            stackedSignalsLeft = np.vstack((self.filtered405Left[i], np.ones(len(self.filtered405Left[i])))).T
            stackedSignalsRight = np.vstack((self.filtered405Right[i], np.ones(len(self.filtered405Left[i])))).T

            gainL, offsetL = np.linalg.lstsq(stackedSignalsLeft, self.filtered470Left[i], rcond=None)[0] 
            gainR, offsetR = np.linalg.lstsq(stackedSignalsRight, self.filtered470Right[i], rcond=None)[0] 

            leftFittedSignal = np.array(list(map(lambda x: x*gainL+offsetL, self.filtered405Left[i])))
            rightFittedSignal = np.array(list(map(lambda x: x*gainR+offsetR, self.filtered405Right[i])))

            self.filtered405Left[i] = leftFittedSignal
            self.filtered405Right[i] = rightFittedSignal

        # self.filtered405Left = acl.wholeTraceGauss(leftFittedArray)
        # self.filtered405Right = acl.wholeTraceGauss(rightFittedArray)

        # self.ratioSignalLeft, self.ratioSignalRight = acl.ratio470405(self.filtered470Left, self.filtered405Left), acl.ratio470405(self.filtered470Right, self.filtered405Right)

        self.ratioSignalLeft = (self.filtered470Left/self.filtered405Left)
        self.ratioSignalRight = (self.filtered470Right/self.filtered405Right)

        self.signalLeft, self.signalRight = acl.wholeTraceGauss(self.ratioSignalLeft), acl.wholeTraceGauss(self.ratioSignalRight)
        
        return self.signalLeft, self.signalRight, self.filtered405Left, self.filtered405Right

# Processor for biexpActive
    def activeSignalDecayProcessor(self):
        decayFitLeft,  rSquaredLeft = acl.doubleExpDecaySingleTrace(self.filtered405Left, self.samplingFreqSec)
        decayFitRight, rSquaredRight = acl.doubleExpDecaySingleTrace(self.filtered405Right, self.samplingFreqSec)

        decayFitLeftScaled = np.zeros_like(decayFitLeft)
        decayFitRightScaled = np.zeros_like(decayFitRight)
        
        for i in self.abf.sweepList:
            stackedSignalsLeft = np.vstack((decayFitLeft[i], np.ones(len(decayFitLeft[i])))).T
            stackedSignalsRight = np.vstack((decayFitRight[i], np.ones(len(decayFitRight[i])))).T

            gainL, offsetL = np.linalg.lstsq(stackedSignalsLeft, self.filtered470Left[i], rcond=None)[0] 
            gainR, offsetR = np.linalg.lstsq(stackedSignalsRight, self.filtered470Right[i], rcond=None)[0] 

            leftFittedSignal = np.array(list(map(lambda x: x*gainL+offsetL, decayFitLeft[i])))
            rightFittedSignal = np.array(list(map(lambda x: x*gainR+offsetR, decayFitRight[i])))

            decayFitLeftScaled[i] = leftFittedSignal - leftFittedSignal[-1]
            decayFitRightScaled[i] = rightFittedSignal - rightFittedSignal[-1]

            decayFitLeft[i] -= decayFitLeft[i][-1]
            decayFitRight[i] -= decayFitRight[i][-1]

        self.filtered470Left -= decayFitLeftScaled
        self.filtered405Left -= decayFitLeft
                
        self.filtered470Right -= decayFitRightScaled
        self.filtered405Right -= decayFitRight

        self.ratioSignalLeft = (self.filtered470Left/self.filtered405Left)
        self.ratioSignalRight = (self.filtered470Right/self.filtered405Right)

        self.signalLeft, self.signalRight = acl.wholeTraceGauss(self.ratioSignalLeft), acl.wholeTraceGauss(self.ratioSignalRight)

        return self.signalLeft, self.signalRight, decayFitLeft, decayFitRight

# Original (Outdated) biexp function
    def newCompleteProcessor(self, salineStatus, ratNameLeft, ratNameRight, experimentDate):
        if salineStatus == 1:
            decayFit405Left, r405Left = acl.doubleExpDecayFit(self.filtered405Left)
            decayFit470Left, r470Left = acl.doubleExpDecayFit(self.filtered470Left)
            decayFit405Right, r405Right = acl.doubleExpDecayFit(self.filtered405Right)
            decayFit470Right, r470Right = acl.doubleExpDecayFit(self.filtered470Right)
            unbleached405Left, unbleached405Right = acl.unbleachSignal(self.filtered405Left, decayFit405Left), acl.unbleachSignal(self.filtered405Right, decayFit405Right)
            unbleached470Left, unbleached470Right = acl.unbleachSignal(self.filtered470Left, decayFit470Left), acl.unbleachSignal(self.filtered470Right, decayFit470Right)
            unbleachedSignals = [decayFit405Left, decayFit470Left, decayFit405Right, decayFit470Right]
            rSquaredDecay = [r405Left, r470Left, r405Right, r470Right]
            csvNames = ['Left405Decay.csv', 'Left470Decay.csv', 'Right405Decay.csv', 'Right470Decay.csv']
            for x in csvNames:
                csvPath = Path(x)
                if csvPath.exists():
                    pass
                elif not csvPath.exists():
                    with open(x, 'w') as csvfile:
                        blankwriter = csv.writer(csvfile)
                        headers = ["Experiment"]
                        match x:
                            case 'Left405Decay.csv':
                                headers.extend(["Trace %i"%(y+1) for y in range(0, len(decayFit405Left))])
                            case 'Left470Decay.csv':
                                headers.extend(["Trace %i"%(y+1) for y in range(0, len(decayFit470Left))])
                            case 'Right405Decay.csv':
                                headers.extend(["Trace %i"%(y+1) for y in range(0, len(decayFit405Right))])
                            case 'Right470Decay.csv':
                                headers.extend(["Trace %i"%(y+1) for y in range(0, len(decayFit470Right))])
                        blankwriter.writerow(headers)
            j = 0
            for x in [decayFit405Left, decayFit470Left, decayFit405Right, decayFit470Right]:
                headers = ["Experiment"]
                headers.extend(["Trace %i"%(y+1) for y in range(0, len(x))])
                containsCheck = True
                inDecayCheck = pd.read_csv(csvNames[j])
                date = experimentDate.strftime("%Y-%m-%d")
                if (inDecayCheck.empty) or (not inDecayCheck.empty and ("%s Rat %s"%(date, ratNameLeft) not in inDecayCheck.index)):
                    labelledFrame = pd.DataFrame({'Experiment':["%s Rat %s"%(date, ratNameLeft)]})
                    unlabelledFrame, storageFrame = pd.DataFrame(), pd.DataFrame()
                    dictOfDecay = {}
                    containsCheck = False
                    for index, entry in enumerate(x):
                        dictOfDecay["Trace %i"%(index+1)] = [entry]
                    storageFrame = unlabelledFrame.assign(**dictOfDecay)
                    concatFrame = pd.concat([labelledFrame, storageFrame], axis=1)
                if not containsCheck:
                    if len(concatFrame.columns) < len(inDecayCheck.columns):
                        commaNum = len(inDecayCheck.columns) - len(concatFrame.columns)
                        paddedComma = pd.DataFrame([[np.nan * commaNum]])
                        headers.extend(["Trace %i"%(y+1) for y in range(len(concatFrame.columns), len(concatFrame.columns)+commaNum)])
                        appendedFrame = pd.concat([concatFrame, paddedComma], ignore_index=True, axis=1)
                    appendedFrame = pd.concat([inDecayCheck, concatFrame], ignore_index=True)
                    appendedFrame = appendedFrame.set_index("Experiment", drop=True).rename_axis(None)                    
                    appendedFrame.to_csv(csvNames[j], header=headers[1:], index_label=headers[0])
                elif containsCheck:
                    continue
                j += 1  
        elif salineStatus == 0:
            decayFit405Left, decayFit470Left, decayFit405Right, decayFit470Right = acl.averageCSV(ratNameLeft, ratNameRight)
            unbleached405Left, unbleached405Right = acl.unbleachSignal(self.filtered405Left, decayFit405Left), acl.unbleachSignal(self.filtered405Right, decayFit405Right)
            unbleached470Left, unbleached470Right = acl.unbleachSignal(self.filtered470Left, decayFit470Left), acl.unbleachSignal(self.filtered470Right, decayFit470Right)
            unbleachedSignals = [decayFit405Left, decayFit470Left, decayFit405Right, decayFit470Right]
            rSquaredDecay = 0
        else:
            raise ValueError(":(")
        
        signalLeft, signalRight = acl.isoLinReg(unbleached405Left, unbleached470Left), acl.isoLinReg(unbleached405Right, unbleached470Right)
        return signalLeft, signalRight, unbleachedSignals, rSquaredDecay