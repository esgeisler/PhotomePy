import statistics as stat

def traceAverage(processedSignal):
    meanSignal = {}
    x = 0
    for sweep in processedSignal:
        meanSignal[x] = stat.mean(sweep)
        x += 1
    return meanSignal

def preInjectionAverage(averagedSignal, injectionTrace):
    preInjList = []
    for sweep in averagedSignal[:injectionTrace]:
        preInjList.append(sweep)
    preInjAvg = stat.mean(preInjList)
    return preInjAvg