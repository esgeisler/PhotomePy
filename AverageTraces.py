import statistics as stat
import numpy as np

# Averages every sweep's voltage (fluorescence) in a trace, producing a dictionary of 70 averages
def traceAverage(processedSignal):
    meanSignal = np.zeros(len(processedSignal))
    for x in range(len(processedSignal)):
        meanSignal[x] = stat.mean(processedSignal[x][1000:-1250])
    return meanSignal

# Averages the sweeps before a certain time point to create a "Pre-Injection Average" of fluorescence
def preInjectionAverage(processedSignal, injectionTrace):
    preInjArray = np.zeros(injectionTrace)
    for x in range(0, injectionTrace):
        preInjArray[x] = stat.mean(processedSignal[x][1000:-1250])
    preInjAvg = stat.mean(preInjArray)
    return preInjAvg

# Calculates ΔF/F ((trace avg - pre inj avg)/pre inj avg)
def deltaF(averagedSignal, preInjAvg):
    deltaFDivided = np.zeros(len(averagedSignal))
    for sweep in range(len(averagedSignal)):
        deltaFDivided[sweep] = (averagedSignal[sweep] - preInjAvg)/preInjAvg
    return deltaFDivided