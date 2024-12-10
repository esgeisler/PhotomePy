import statistics as stat
import numpy as np

# Averages every sweep's voltage (fluorescence) in a trace, producing a dictionary of 70 averages
def traceAverage(processedSignal):
    meanSignal = np.zeros(len(processedSignal))
    medianSignal = np.zeros(len(processedSignal))
    for i, x in enumerate(processedSignal):
        meanSignal[i] = stat.mean(x[1000:-1250])
        medianSignal[i] = stat.median(x[100:-1250])
    return meanSignal, medianSignal

# Averages the sweeps before a certain time point to create a "Pre-Injection Average" of fluorescence
def preInjectionAverage(processedSignal, injectionTrace):
    if injectionTrace <= 1:
        raise ValueError("Injection traces must be larger than 1")
    elif injectionTrace > len(processedSignal):
        raise IndexError("Injection traces are larger than total sweeps in signals")
    else:
        preInjArray = np.zeros(injectionTrace)
        for x in range(0, injectionTrace):
            preInjArray[x] = stat.mean(processedSignal[x][1000:-1250])
        preInjAvg = stat.mean(preInjArray)
        preInjStdDev = stat.stdev(preInjArray)
        return preInjAvg, preInjStdDev

# Calculates ΔF/F ((trace avg - pre inj avg)/pre inj avg)
def deltaF(averagedSignal, preInjAvg):
    deltaFDivided = np.zeros(len(averagedSignal))
    for i, sweep in enumerate(averagedSignal):
        deltaFDivided[i] = (sweep - preInjAvg)/preInjAvg
    return deltaFDivided

# Calculates Z-score ("Standard Score") from the each trace
def zCalc(averagedSignal, processedSignal, injectionTrace):
    zScoreArray = np.zeros(len(processedSignal))
    sampleMean, sampleSigma = preInjectionAverage(processedSignal, injectionTrace)
    for i, sweep in enumerate(averagedSignal):
        zScoreArray[i] = (sweep - sampleMean)/sampleSigma
    return zScoreArray