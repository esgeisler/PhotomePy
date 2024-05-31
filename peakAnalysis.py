import AutoCleaner as acl
import AverageTraces as avg
import scipy.signal as sci
import pyabf
import matplotlib.pyplot as plt
import matplotlib.pyplot

#Retrieves the peaks of a signal and returns a list containing the peaks in a ndarray and their properties in a dictionary
def peakGetter(processedSignalArray):
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, height=0)
    print([peaks, peaksDict])
    

#Retrieves the peaks of a signal and plots them on a graph of the chosen trace
def peakDisplay(processedSignalArray):
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05)
    fig = plt.figure()
    lad = fig.add_subplot()
    lad.plot(processedSignalArray)
    lad.plot(peaks, processedSignalArray[peaks], "r.")
    plt.axis([0,50000, -2, 2])
    plt.show()

def peakDecay(processedSignalArray):
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, width=0)
    peakdecayList = []
    x = 0
    for peak in peaks:
        peakdecayList.append(peak - peaksDict[3][x])
        x += 1
    return peakdecayList

def peakAmplitude(processedSignalArray):
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, height=0)
    return peaksDict[0]

def peakFreq(processedSignalArray):
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05, height=0)
    return len(peaks)

def peakWidth(processedSignalArray):
    peaks, peaksDict = sci.find_peaks(processedSignalArray, prominence= 0.05)
    width, widthHeight, leftRight = sci.peak_widths(processedSignalArray, peaks)
    return width

# Plots the original, subtracted, ratio-ed, and processed trace of choice
    # abf.setSweep(sweepNumber= userTrace, channel= userChannel)
    # plt.plot(abf.sweepX[1000:-2000], abf.sweepY[1000:-2000], color="b", label="Original")
    # plt.plot(abf.sweepX[1000:-2000], subtractRight[0][userTrace][1000:-2000], color="r", label= "Subtracted")
    # plt.plot(abf.sweepX[1000:-2000], ratioSignalRight[userTrace][1000:-2000], color="y", label= "Ratio-ed")
    # plt.plot(abf.sweepX[1000:-2000], finalSignalRight[userTrace][1000:-2000], color="g", label= "Processed")
    # plt.axis([0,15,-1,6])
    # plt.legend()
    # plt.show()