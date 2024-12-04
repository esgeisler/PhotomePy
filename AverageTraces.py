import statistics as stat
import numpy as np
import scipy.optimize as opt

# Averages every sweep's voltage (fluorescence) in a trace, producing a dictionary of 70 averages
def traceAverage(processedSignal):
    meanSignal = np.zeros(len(processedSignal))
    for i, x in enumerate(processedSignal):
        meanSignal[i] = stat.mean(x[1000:-1250])
    return meanSignal

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

# Takes the numbers from deltaF or zCalc and fits a bi-exponential decay function to them
def doubleExpDecayFit(normalizedSignal):
    p0 = (1,-1,1,-1,1)
    xDecay = np.arange(0, len(normalizedSignal))
    yDecay = [i for i in normalizedSignal]
    popt, _ = opt.curve_fit(lambda t, a, b, c, d, e: (a * np.exp(b * t)) + (c * np.exp(d *t)) + e, xDecay, yDecay, p0=p0, 
                                                bounds=opt.Bounds(lb=[0, -np.inf, 0, -np.inf, min(yDecay)], 
                                                                ub=[np.inf, 0, np.inf, 0, np.inf]),
                                                                maxfev=1000)
    a, b, c, d, e = popt[0], popt[1], popt[2], popt[3], popt[4]
    squaredDiffs = np.square(yDecay - (a * np.exp(b * xDecay)) + (c * np.exp(d * xDecay)) + e)
    squaredDiffsFromMean = np.square(yDecay - np.mean(yDecay))
    rSquared = 1 - (np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean))
    xDecFinal = np.linspace(np.min(xDecay), np.max(xDecay), len(normalizedSignal))
    yDecFinal = (a * np.exp(b * xDecFinal)) + (c * np.exp(d * xDecFinal)) + e
    return xDecFinal, yDecFinal, rSquared