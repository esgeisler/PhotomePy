import AutoCleaner as acl
import AverageTraces as avg
import scipy.signal as sci
import pyabf
import matplotlib.pyplot as plt
import matplotlib.pyplot

#TODO: Write a function that retrieves the number of peaks
def peakGetter(processedSignal):
    peaks, dictOfPeaks = sci.find_peaks(processedSignal, prominence= 0.05)
    fig = plt.figure()
    lad = fig.add_subplot()
    lad.plot(processedSignal)
    lad.plot(peaks, processedSignal[peaks], "r.")
    plt.axis([0,50000, -2, 2])
    plt.show()
    

#TODO: Write a house-keeping function that displays those peaks on a graph
def peakDisplay():
    return


# Plots the original, subtracted, ratio-ed, and processed trace of choice
    # abf.setSweep(sweepNumber= userTrace, channel= userChannel)
    # plt.plot(abf.sweepX[1000:-2000], abf.sweepY[1000:-2000], color="b", label="Original")
    # plt.plot(abf.sweepX[1000:-2000], subtractRight[0][userTrace][1000:-2000], color="r", label= "Subtracted")
    # plt.plot(abf.sweepX[1000:-2000], ratioSignalRight[userTrace][1000:-2000], color="y", label= "Ratio-ed")
    # plt.plot(abf.sweepX[1000:-2000], finalSignalRight[userTrace][1000:-2000], color="g", label= "Processed")
    # plt.axis([0,15,-1,6])
    # plt.legend()
    # plt.show()