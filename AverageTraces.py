import statistics as stat

# Averages every sweep's voltage (fluorescence) in a trace, producing a dictionary of 70 averages
def traceAverage(processedSignal):
    meanSignal = {}
    x = 0
    for sweep in processedSignal:
        meanSignal[x] = stat.mean(processedSignal[sweep])
        x += 1
    return meanSignal

# Averages the sweeps before a certain time point to create a "Pre-Injection Average" of fluorescence
def preInjectionAverage(averagedSignal, injectionTrace):
    preInjList = []
    x = 0
    for x in range(0, injectionTrace):
        preInjList.append(stat.mean(averagedSignal[x]))
    preInjAvg = stat.mean(preInjList)
    return preInjAvg

# Calculates ΔF/F ((trace avg - pre inj avg)/pre inj avg)
def deltaF(averagedSignal, preInjAvg):
    deltaFDivided = {}
    for sweep in averagedSignal:
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