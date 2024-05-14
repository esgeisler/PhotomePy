import pyabf
import statistics as stat
import pyabf.filter
import scipy
import scipy.ndimage

# Gets baseline information from 1 min-long recording data taken after trial from the "left" side of the room - channels 1 and 2
def LBaselineGet(FileName):
    abf = pyabf.ABF(FileName)
    sweepDict470 = {}
    Channel470 = 0
    sweepDict405 = {}
    Channel405 = 1

    abf.setSweep(sweepNumber = 0, channel = Channel470)
    for sweep470 in abf.sweepList:
        fluor470 = stat.mean(abf.sweepY)
        sweepDict470[sweep470] = fluor470
    mean470 = stat.mean(sweepDict470.values())
    
    abf.setSweep(sweepNumber = 0, channel = Channel405)
    for sweep405 in abf.sweepList:
        fluor405 = stat.mean(abf.sweepY)
        sweepDict405[sweep405] = fluor405
    mean405 = stat.mean(sweepDict405.values())

    return mean470, mean405

# Gets baseline information from 1 min-long recording data taken after trial from the "right" side of the room - channels 5 and 6
def RBaselineGet(FileName):
    abf = pyabf.ABF(FileName)
    sweepDict470 = {}
    Channel470 = 4
    sweepDict405 = {}
    Channel405 = 5

    abf.setSweep(sweepNumber = 0, channel = Channel470)
    for sweep470 in abf.sweepList:
        fluor470 = stat.mean(abf.sweepY)
        sweepDict470[sweep470] = fluor470
    mean470 = stat.mean(sweepDict470.values())
    
    abf.setSweep(sweepNumber = 0, channel = Channel405)
    for sweep405 in abf.sweepList:
        fluor405 = stat.mean(abf.sweepY)
        sweepDict405[sweep405] = fluor405
    mean405 = stat.mean(sweepDict405.values())

    return mean470, mean405

# Gaussian filters a single trace with a 40 Hz cutoff frequency, based on pClamp documentation and Calquhon & Sigworth (1995)
def gaussianFilter(FileName, filterChannel, filterSweep):
    abf = pyabf.ABF(FileName)
    abf.setSweep(sweepNumber= filterSweep, channel= filterChannel)
    sweepList = scipy.ndimage.gaussian_filter1d(abf.sweepY, sigma= 22)
    return sweepList

def wholeTraceGauss(fileName, filterChannel):
    abf = pyabf.ABF(fileName)
    sweepDict = {}
    abf.setSweep(sweepNumber= 0, channel= filterChannel)
    for sweeps in abf.sweepList:
        abf.setSweep(sweeps, channel= filterChannel)
        filteredSweep = scipy.ndimage.gaussian_filter1d(abf.sweepY, sigma = 16)
        sweepDict[sweeps] = filteredSweep
    return sweepDict