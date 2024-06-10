import statistics as stat
import pandas as pd
import numpy as np

# Averages every sweep's voltage (fluorescence) in a trace, producing a dictionary of 70 averages
def traceAverage(processedSignal):
    meanSignal = np.zeros(len(processedSignal))
    for x in range(len(processedSignal)):
        meanSignal[x] = stat.mean(processedSignal[x][500:-1500])
    return meanSignal

# Averages the sweeps before a certain time point to create a "Pre-Injection Average" of fluorescence
def preInjectionAverage(averagedSignal, injectionTrace):
    preInjList = []
    for x in range(0, injectionTrace):
        preInjList.append(stat.mean(averagedSignal[x][500:-1500]))
    preInjAvg = stat.mean(preInjList)
    return preInjAvg

# Calculates ΔF/F ((trace avg - pre inj avg)/pre inj avg)
def deltaF(averagedSignal, preInjAvg):
    deltaFDivided = np.zeros(len(averagedSignal))
    for sweep in range(len(averagedSignal)):
        deltaFDivided[sweep] = (averagedSignal[sweep] - preInjAvg)/preInjAvg
    return deltaFDivided

# Bleaching correction that returns difference between ΔF/F (drug) and ΔF/F (vehicle)
def bleachingSub(signalDrug, signalVehicle):
    finalSignal = {}
    x = 0
    dictLen = len(signalDrug)
    for x in dictLen:
        finalSignal[x] = signalDrug[x] - signalVehicle[x]
        x += 1
    return finalSignal

# Creates a pandas dataframe that can be exported into excel
def excelExporter(signalAverage, preInjectionAverage, deltaF):
    exportableData = pd.DataFrame({"Trace Number:": range(1, len(signalAverage)+1), "Average Fluorescence": signalAverage, 
                                   "Pre-Injection Average":preInjectionAverage, "ΔF/F": deltaF, "Bleaching Correction": None, })
    return exportableData