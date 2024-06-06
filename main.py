import AutoCleaner as acl
import pyabf
import os
import matplotlib.pyplot as plt
import AverageTraces as avg
import peakAnalysis as pas
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
from datetime import datetime

#TODO Change "Process File" button to do averages and traces, export both to excel, and save the modified file as a .abf
class Main(tk.Frame):
    def __init__(self, master= None, **kwargs):
        super().__init__(master, **kwargs)

    # Initialize all of the necessary variables for the GUI
        self.experimentFileName = ""
        self.baselinefileName = ""
        self.leftRatName = tk.StringVar()
        self.rightRatName = tk.StringVar()
        self.leftRatInjection = tk.StringVar()
        self.rightRatInjection = tk.StringVar()
        self.dropValue = tk.StringVar(self, 'Select Trace')
        self.ratNameLeft = 0
        self.ratNameRight = 0
        self.ratInjectionLeft = 0
        self.ratInjectionRight = 0
        self.peaksLeft = []
        self.peaksRight = []
        self.options = []
        self.trace = 0

    # Dropdown menu for selecting the trace used in SingleTracePeaks
        traceDrop = ttk.Combobox(self,  state= 'readonly', width= 9, textvariable= self.dropValue)
        traceDrop['values'] = []
        traceDrop.grid(row= 4, column=2)

    # Sets the values inside the dropdown menu, and sets it to update when a main file is selected
        def traceSelector(event):
            self.trace = (int(traceDrop.get()))
        traceDrop.bind("<<ComboboxSelected>>", traceSelector)
        def dropdownUpdater():
            traceDrop['values'] = self.options

    # Opens the main file, containing data from the session with 2 rats
        def fileBrowserExperiment():
            self.experimentFileName = filedialog.askopenfilename(initialdir= os.getcwd(), title= "Select a Main File", 
                                                                 filetypes=(("Axon Binary Fles", "*.abf*"), ("All Files," "*.*")))
            chosenFileDisplay.insert(tk.END, self.experimentFileName)
            abf = pyabf.ABF(self.experimentFileName)
            self.options = [str(x + 1) for x in abf.sweepList]
            return self.experimentFileName
    
    # Opens the baseline file containing the baseline autofluorescence
        def fileBrowserBaseline():
            self.baselinefileName = filedialog.askopenfilename(initialdir= os.getcwd(), title= "Select a Main File", filetypes=(("Axon Binary Fles", "*.abf*"), ("All Files," "*.*")))
            baselineFileDisplay.insert(tk.END, self.baselinefileName)
            return self.baselinefileName  

    # Initialize all of the remaining buttons for the main GUI window
        explorerButton = tk.Button(self, text="Choose a Main File", command= lambda:[fileBrowserExperiment(), dropdownUpdater()])
        chosenFileDisplay = tk.Text(self, height= 1, width= 50)
        chosenFileDisplay.grid(row= 1, column= 1)
        baselineExplorerButton = tk.Button(self, text="Choose a Baseline File", command= fileBrowserBaseline)
        explorerButton.grid(row= 1, column= 2)
        baselineFileDisplay = tk.Text(self, height= 1, width= 50)
        baselineFileDisplay.grid(row= 2, column= 1)
        baselineExplorerButton.grid(row= 2, column= 2)
    # Fills text boxes with filepath of the main and baseline file chosen by the user
        chosenFileDisplay.insert(tk.END, self.experimentFileName)
        baselineFileDisplay.insert(tk.END, self.baselinefileName)

    # Closes the average processing popup window, saving the values entered.
        def onPopSubmit():
            self.ratNameLeft = int(self.leftRatName.get())
            self.ratNameRight = int(self.rightRatName.get())
            self.ratInjectionLeft = int(self.leftRatInjection.get())
            self.ratInjectionRight = int(self.rightRatInjection.get())

    # Creates a popup window for the user to input the rat's name/number and what trace they were injected during.
        def dataProcessorPop():
            infoPop = tk.Toplevel()
            infoPop.title("Rat Metadata Entry")
            leftRatNameFill = tk.Entry(infoPop, textvariable= self.leftRatName, width= 10)
            leftRatNameLabel = tk.Label(infoPop, text="Enter Left Rat #:")
            rightRatNameFill = tk.Entry(infoPop, textvariable= self.rightRatName, width= 10)
            rightRatNameLabel = tk.Label(infoPop, text= "Enter Right Rat #:")
            leftRatNameLabel.grid(row= 1, column= 1)
            leftRatNameFill.grid(row= 1, column= 2)
            rightRatNameLabel.grid(row= 1, column= 4)
            rightRatNameFill.grid(row= 1, column= 5)
            leftRatInjTimeFill = tk.Entry(infoPop, textvariable= self.leftRatInjection, width = 10)
            leftRatInjTimeLabel = tk.Label(infoPop, text="Enter Left Rat Injection Trace #:")
            rightRatInjTimeFill = tk.Entry(infoPop, textvariable= self.rightRatInjection, width= 10)
            rightRatInjTimeLabel = tk.Label(infoPop, text= "Enter Right Rat Injection Trace #:")
            leftRatInjTimeLabel.grid(row= 2, column= 1)
            leftRatInjTimeFill.grid(row= 2, column= 2)
            rightRatInjTimeLabel.grid(row= 2, column= 4)
            rightRatInjTimeFill.grid(row= 2, column= 5)


            submitButton = tk.Button(infoPop, text="Submit", command=lambda:[onPopSubmit(), infoPop.destroy(), dataProcessorReal()])
            submitButton.grid(row= 3, column=3)

    # Averages the fluorescence of all of the traces, compiles them into an excel sheet with their trace numbers, and calculates the ΔF/F
        def dataProcessorReal():
            finalSignalLeft, finalSignalRight = acl.completeProcessor(self.experimentFileName, self.baselinefileName)
        #     baselineSubL = acl.LBaselineGet(self.baselinefileName)
        #     baselineSubR = acl.RBaselineGet(self.baselinefileName)
        #     channelsLeft = [0,1]
        #     channelsRight = [4,5]
        #     subtractLeft = acl.baselineSubtractor(self.experimentFileName, baselineSubL, channelsLeft)
        #     subtractRight = acl.baselineSubtractor(self.experimentFileName, baselineSubR, channelsRight)
        # # Gaussian filters the 405 channels
        #     filteredLeft = acl.wholeTraceGauss(subtractLeft[1])
        #     filteredRight = acl.wholeTraceGauss(subtractRight[1])
        # #Find ratio of 470/405 channels
        #     ratioSignalLeft = acl.ratio470405(subtractLeft[0], filteredLeft)
        #     ratioSignalRight = acl.ratio470405(subtractRight[0], filteredRight)
        # # Gaussian filters the ratio signal
        #     finalSignalLeft = acl.wholeTraceGauss(ratioSignalLeft)
        #     finalSignalRight = acl.wholeTraceGauss(ratioSignalRight)
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
            filenameLeft = "%s Rat %i Temp File.xlsx"%(datetime.today().strftime('%Y-%m-%m'), self.ratNameLeft)
            filenameRight = "%s Rat %i Temp File.xlsx"%(datetime.today().strftime('%Y-%m-%m'), self.ratNameRight)
            ratDataLeft.to_excel(filenameLeft)
            ratDataRight.to_excel(filenameRight)
            messagebox.showinfo(title= "Data Exporter", message= "Data Exported to Excel!")

    # Retrieves the baseline autofluorescence for the 4 channels analyzed and prints to a message box.
        def baselineFinder():
            pyabf.ABF(self.baselinefileName)
            baselineSubL = acl.LBaselineGet(self.baselinefileName)
            baselineSubR = acl.RBaselineGet(self.baselinefileName)
            messagebox.showinfo(title= "Baselines", message= "Left - 470: %.2f 405: %.2f\nRight - 470: %.2f 405: %.2f"%(baselineSubL[0], baselineSubL[1], baselineSubR[0], baselineSubR[1]))

    # Analyzes the peak decay, amplitude, and frequency of a single trace chosen by the user.
        def singleTracePeaks():
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
        # Left Rat Peak Analysis
            signalValuesLeft = np.array(list(finalSignalLeft.values()))
            self.peaksLeft = pas.peakGetter(signalValuesLeft[self.trace][1000:-1250])
            pas.peakDisplay(signalValuesLeft[self.trace][1000:-1250], self.experimentFileName, "Left Rat")
            plt.close()
        # Right Rat Peak Analysis
            signalValuesRight = np.array(list(finalSignalRight.values()))
            self.peaksRight = pas.peakGetter(signalValuesRight[self.trace][1000:-1250])
            pas.peakDisplay(signalValuesRight[self.trace][1000:-1250], self.experimentFileName, "Right Rat")

    # Analyzes peak decay, amplitude, and frequency across an entire signal containing X traces
        def peakAnalyzer():
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

            signalValuesLeft = np.array(list(finalSignalLeft.values()))
            self.peaksLeft = pas.wholeTracePeaks(signalValuesLeft)

            signalValuesRight = np.array(list(finalSignalRight.values()))
            self.peaksRight = pas.wholeTracePeaks(signalValuesRight)

            preInjectionLeft, postInjectionLeft = pas.traceProcessor(self.peaksLeft)
            preInjectionRight, postInjectionRight = pas.traceProcessor(self.peaksRight)

            preLeft = "%s Rat %i Pre-Injection Peaks.xlsx"%(datetime.today().strftime('%Y-%m-%m'), self.ratNameLeft)
            postLeft = "%s Rat %i Post-Injection Peaks.xlsx"%(datetime.today().strftime('%Y-%m-%m'), self.ratNameLeft)

            preRight = "%s Rat %i Pre-Injection Peaks.xlsx"%(datetime.today().strftime('%Y-%m-%m'), self.ratNameRight)
            postRight = "%s Rat %i Post-Injection Peaks.xlsx"%(datetime.today().strftime('%Y-%m-%m'), self.ratNameRight)
            
            preInjectionLeft.to_excel(preLeft)
            postInjectionLeft.to_excel(postLeft)
            preInjectionRight.to_excel(preRight)
            postInjectionRight.to_excel(postRight)

            messagebox.showinfo(title= "Trace Exporter", message= "Data Exported to Excel!")

            
        runFileButton = tk.Button(self, text="Process a File", command= dataProcessorPop)
        runFileButton.grid(row= 3, column= 1)
        baselineGetterButton = tk.Button(self, text="Get the Baselines", command= baselineFinder)
        baselineGetterButton.grid(row= 3, column= 2)

        testerButton = tk.Button(self, text="Event analysis on a single trace", command= singleTracePeaks)
        testerButton.grid(row=4, column=1)

        peakButton = tk.Button(self, text= "Perform event analysis on an entire signal [WIP]", command= peakAnalyzer)
        peakButton.grid(row=5, column=1)     

def main():
    fp = tk.Tk()
    fp.title("Liz's Data Corrector")
    window = Main(fp)
    window.pack()
    fp.mainloop()

if __name__ == "__main__":
    main()