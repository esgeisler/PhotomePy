import AutoCleaner as acl
import pyabf
import matplotlib.pyplot as plt
import AverageTraces as avg

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog



class Main(tk.Frame):
    def __init__(self, master= None, **kwargs):
        super().__init__(master, **kwargs)

        self.experimentFileName = ""
        self.baselinefileName = ""
        self.leftRatName = tk.StringVar()
        self.rightRatName = tk.StringVar()
        self.leftRatInjection = tk.StringVar()
        self.rightRatInjection = tk.StringVar()
        self.ratNameLeft = 0
        self.ratNameRight = 0
        self.ratInjectionLeft = 0
        self.ratInjectionRight = 0

        def fileBrowserExperiment():
            self.experimentFileName = filedialog.askopenfilename(initialdir= "/", title= "Select a Main File", filetypes=(("Axon Binary Fles", "*.abf*"), ("All Files," "*.*")))
            print(self.experimentFileName)
            chosenFileDisplay.insert(tk.END, self.experimentFileName)
            return self.experimentFileName
        
        def fileBrowserBaseline():
            self.baselinefileName = filedialog.askopenfilename(initialdir= "/", title= "Select a Main File", filetypes=(("Axon Binary Fles", "*.abf*"), ("All Files," "*.*")))
            print(self.baselinefileName)
            baselineFileDisplay.insert(tk.END, self.baselinefileName)
            return self.baselinefileName
        
        decisionLabel = tk.Label(self, text="File Process")
        decisionLabel.pack()  

        explorerButton = tk.Button(self, text="Choose a Main File", command= fileBrowserExperiment)
        chosenFileDisplay = tk.Text(self, height= 1, width= 75)
        chosenFileDisplay.pack()
        baselineExplorerButton = tk.Button(self, text="Choose a Baseline File", command= fileBrowserBaseline)
        explorerButton.pack()
        baselineFileDisplay = tk.Text(self, height= 1, width= 75)
        baselineFileDisplay.pack()
        baselineExplorerButton.pack()

        chosenFileDisplay.insert(tk.END, self.experimentFileName)
        baselineFileDisplay.insert(tk.END, self.baselinefileName)

        def onPopSubmit():
            self.ratNameLeft = int(self.leftRatName.get())
            self.ratNameRight = int(self.rightRatName.get())
            self.ratInjectionLeft = int(self.leftRatInjection.get())
            self.ratInjectionRight = int(self.rightRatInjection.get())

        def dataProcessorPop():
            infoPop = tk.Toplevel()
            infoPop.title("Notice")
            leftRatNameFill = tk.Entry(infoPop, textvariable= self.leftRatName)
            rightRatNameFill = tk.Entry(infoPop, textvariable= self.rightRatName)
            leftRatNameFill.pack()
            rightRatNameFill.pack()
            leftRatInjTimeFill = tk.Entry(infoPop, textvariable= self.leftRatInjection)
            rightRatInjTimeFill = tk.Entry(infoPop, textvariable= self.rightRatInjection)
            leftRatInjTimeFill.pack()
            rightRatInjTimeFill.pack()
            submitButton = tk.Button(infoPop, text="Submit", command=lambda:[onPopSubmit(), infoPop.destroy(), dataProcessorReal()])
            submitButton.pack()

        def dataProcessorReal():
            abf = pyabf.ABF(self.experimentFileName)
            baselineSubL = acl.LBaselineGet(self.baselinefileName)
            baselineSubR = acl.RBaselineGet(self.baselinefileName)
            channelsLeft = [0,1]
            channelsRight = [4,5]
            subtractLeft = acl.baselineSubtractor(self.experimentFileName, baselineSubL, channelsLeft)
            subtractRight = acl.baselineSubtractor(self.experimentFileName, baselineSubR, channelsRight)
        # Gaussian filters the 405 channels
            filteredLeft = acl.wholeTraceGauss(subtractLeft[1])
            filteredRight = acl.wholeTraceGauss(subtractRight[1])
        #Find ratio of 470/405 channels
            ratioSignalLeft = acl.ratio470405(subtractLeft[0], filteredLeft)
            ratioSignalRight = acl.ratio470405(subtractRight[0], filteredRight)
        # Gaussian filters the ratio signal
            finalSignalLeft = acl.wholeTraceGauss(ratioSignalLeft)
            finalSignalRight = acl.wholeTraceGauss(ratioSignalRight)
        # Averages the left and right signals
            averageSignalLeft = avg.traceAverage(finalSignalLeft)
            preInjectionAverageLeft = avg.preInjectionAverage(finalSignalLeft, self.ratInjectionLeft)
            fluorescenceLeft = avg.deltaF(averageSignalLeft, preInjectionAverageLeft)

            averageSignalRight = avg.traceAverage(finalSignalRight)
            preInjectionAverageRight = avg.preInjectionAverage(finalSignalRight, self.ratInjectionRight)
            fluorescenceRight = avg.deltaF(averageSignalRight, preInjectionAverageRight)
        # Saves the averaged data to an excel file with the rat's "name"
            ratDataLeft = avg.excelExporter(averageSignalLeft, preInjectionAverageLeft, fluorescenceLeft)
            ratDataRight = avg.excelExporter(averageSignalRight, preInjectionAverageRight, fluorescenceRight)
            filenameLeft = "Rat %i Temp File.xlsx"%(self.ratNameLeft)
            filenameRight = "Rat %i Temp File.xlsx"%(self.ratNameRight)
            ratDataLeft.to_excel(filenameLeft)
            ratDataRight.to_excel(filenameRight)

            puEnd = tk.Toplevel()
            puEnd.wm_title("Results")
            textyBox = tk.Text(puEnd, width= 25, height=1)
            textyBox.pack()
            textyBox.insert(tk.END, "Results Exported to Excel!")
    
        def baselineFinder():
            pyabf.ABF(self.baselinefileName)
            # Opens a pop-up window containing the 4 baselines for doing calculations with
            bpu = tk.Toplevel()
            bpu.wm_title("Baselines")
            baselineSubL = acl.LBaselineGet(self.baselinefileName)
            baselineSubR = acl.RBaselineGet(self.baselinefileName)
            baselineTextBox = tk.Text(bpu, "Left - 470: %.2f 405: %.2f\nRight - 470: %.2f 405: %.2f"%(baselineSubL[0], baselineSubL[1], baselineSubR[0], baselineSubR[1]))
            baselineTextBox.pack()

        runFileButton = tk.Button(self, text="Process a File", command= dataProcessorPop)
        runFileButton.pack()
        baselineGetterButton = tk.Button(self, text="Get the Baselines", command= baselineFinder)
        baselineGetterButton.pack()

def main():
    fp = tk.Tk()
    fp.title("Liz's Data Corrector")
    window = Main(fp)
    window.pack()
    fp.mainloop()

if __name__ == "__main__":
    main()