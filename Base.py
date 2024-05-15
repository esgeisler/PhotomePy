import FileSelector as fls
import Exceptions as exn
import AutoCleaner as acl
import pyabf
import os
import matplotlib.pyplot as plt

# ToDo = input("What should we do?")

userFile = fls.FileSelect(input("Is the data raw or processed?"), input("What is the name of the file? (XX.abf)"))

abf = pyabf.ABF(userFile)

userTrace = int(input("Which trace would you like to see?"))
userChannel = int(input("Which channel?"))
rawFile = fls.FileSelect("raw", "24508006.abf")
baselineSubL = acl.LBaselineGet(rawFile)
baselineSubR = acl.RBaselineGet(rawFile)
channelsLeft = [0,1]
channelsRight = [4, 5]
abf.setSweep(userTrace, userChannel)

subtractLeft = acl.baselineSubtractor(userFile, baselineSubL, channelsLeft)
subtractRight = acl.baselineSubtractor(userFile, baselineSubR, channelsRight)

plt.plot(abf.sweepX[1000:-2000], abf.sweepY[1000:-2000], color="b", label="original")
plt.plot(abf.sweepX[1000:-2000], subtractRight[1][userTrace][1000:-2000], color="y", label= "Subtracted")
# filteredPlot = acl.wholeTraceGauss(userFile, userChannel)
# plt.plot(abf.sweepX[1000:-2000], filteredPlot[userTrace][1000:-2000], color="r", label= "Filtered")
# plt.axis([0,15,4.6,5.5])
plt.legend()
plt.show()

# print(acl.LBaselineGet(userFile))
# print(acl.RBaselineGet(userFile))