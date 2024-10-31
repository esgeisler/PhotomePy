import AutoCleaner as acl
import pyabf
import os
import sys
import AverageTraces as avg
import peakAnalysis as pas
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd

class Main(tk.Frame):
    def __init__(self, master= None, **kwargs):
        super().__init__(master, **kwargs)

    # Initialize all of the necessary variables for the GUI
        self.experimentFileName = ""
        self.baselinefileName = ""
        self.leftRatName, self.rightRatName = tk.StringVar(), tk.StringVar()
        self.leftRatInjectionStr, self.rightRatInjectionStr = tk.StringVar(), tk.StringVar()
        self.dropValue = tk.StringVar(self, 'Select Trace')
        self.ratNameLeft, self.ratNameRight = tk.StringVar(), tk.StringVar()
        self.leftRatInjectionInt, self.rightRatInjectionInt = 0, 0
        self.peaksLeft, self.peaksRight = [], []
        self.options = []
        self.trace = 0
        self.abfDate = tk.StringVar()
        self.statusCheck = True
        self.mainStatus = False
        self.traceStatus = False
        self.errorStatus = False

    # Dropdown menu for selecting the trace used in SingleTracePeaks
        traceDrop = ttk.Combobox(self,  state= 'readonly', width= 9, textvariable= self.dropValue)
        traceDrop['values'] = []
        traceDrop.grid(row= 4, column=2)

    # Sets the values inside the dropdown menu, and sets it to update when a main file is selected
        def traceSelector():
            self.trace = (int(traceDrop.get()) - 1)
        traceDrop.bind("<<ComboboxSelected>>", traceSelector)
        def dropdownUpdater():
            traceDrop['values'] = self.options

    # Opens the main file, containing data from the session with 2 rats
        def fileBrowserExperiment():
            x = 0
            while x == 0:
                try:
                    self.experimentFileName = filedialog.askopenfilename(initialdir= os.path.join(os.getcwd(), "Raw Data"), title= "Select a Main File", 
                                                                    filetypes=(("Axon Binary Fles", "*.abf*"), ("All Files," "*.*")))
                    abf = pyabf.ABF(self.experimentFileName)
                    if len(abf.sweepList) <= 1:
                        raise ValueError("Main files must have multiple sweeps")
                    self.options = [str(x + 1) for x in abf.sweepList]
                    self.abfDate = abf.abfDateTime
                except ValueError:
                    messagebox.showerror(title= "Python Error", message= "Please select a main file with multiple sweeps")
                    continue
                except PermissionError as e:
                    if e.errno == 13:
                        x = 1
                    else:
                        raise IOError("Something is wrong with this file")
                else:
                    chosenFileDisplay.insert(tk.END, self.experimentFileName)
                    return self.experimentFileName
    
    # Opens the baseline file containing the baseline autofluorescence
        def fileBrowserBaseline():
            x = 0
            while x == 0:
                try:
                    self.baselinefileName = filedialog.askopenfilename(initialdir= os.path.join(os.getcwd(), "Raw Data"), title= "Select a Baseline File", 
                                                                       filetypes=(("Axon Binary Fles", "*.abf*"), ("All Files," "*.*")))
                    abf = pyabf.ABF(self.baselinefileName)
                    if len(abf.sweepList) > 1:
                        raise ValueError("Baseline files must be single-trace files")
                except ValueError:
                    messagebox.showerror(title= "Python Error", message= "Please select a baseline file with a single sweep")
                    continue
                except PermissionError as e:
                    if e.errno == 13:
                        x = 1
                    else:
                        raise IOError("Something is wrong with this file")
                else:
                    baselineFileDisplay.insert(tk.END, self.baselinefileName)
                    return self.baselinefileName

    # Initialize all of the remaining buttons for the main GUI window
        explorerButton = ttk.Button(self, text="Choose a Main File", command= lambda:[fileBrowserExperiment(), dropdownUpdater(), chosenFileTextUpdate()])
        chosenFileDisplay = tk.Text(self, height= 1, width= 75)
        chosenFileDisplay.grid(row= 1, column= 1)
        baselineExplorerButton = ttk.Button(self, text="Choose a Baseline File", command= lambda:[fileBrowserBaseline(), baselineFileTextUpdate()])
        explorerButton.grid(row= 1, column= 2)
        baselineFileDisplay = tk.Text(self, height= 1, width= 75)
        baselineFileDisplay.grid(row= 2, column= 1)
        baselineExplorerButton.grid(row= 2, column= 2)
    
    # Fills text boxes with filepath of the main and baseline file chosen by the user
        def chosenFileTextUpdate():
            chosenFileDisplay.delete("1.0", tk.END)
            chosenFileDisplay.insert(tk.INSERT, self.experimentFileName)
        def baselineFileTextUpdate():
            baselineFileDisplay.delete("1.0", tk.END)
            baselineFileDisplay.insert(tk.INSERT, self.baselinefileName)

    # Closes the average processing popup window, saving the values entered.
        def onPopSubmit():
            self.ratNameLeft, self.ratNameRight = self.leftRatName.get(), self.rightRatName.get()
            try:
                self.leftRatInjectionInt, self.rightRatInjectionInt = (int(float(self.leftRatInjectionStr.get())) - 1), (int(float(self.rightRatInjectionStr.get())) - 1)
            except ValueError:
                answer = messagebox.askretrycancel(title="Python Error", message="Injection traces must be a number. Please re-enter your trace number", icon="error")
                if answer:
                    self.errorStatus = True
                    self.leftRatInjectionStr, self.rightRatInjectionStr = tk.StringVar(), tk.StringVar()
                    dataProcessorPop()
                    return
                elif not answer:
                    self.destroy()
                    sys.exit(0)
            dataProcessorReal()
            statusChecker()

    # Checks the status of the analysis to determine whether to display a messagebox message.
        def statusChecker():
            while self.statusCheck:
                if self.errorStatus:
                    self.errorStatus = False
                    break
                if not self.traceStatus and not self.mainStatus:
                    continue
                elif self.traceStatus and self.mainStatus:
                    messagebox.showinfo(title= "Trace Exporter", message= "Data Exported to Excel!")
                    self.traceStatus = False
                    self.mainStatus = False
                    break
                
            
    # Creates a popup window for the user to input the rat's name/number and what trace they were injected during.
        def dataProcessorPop():
            infoPop = tk.Toplevel()
            infoPop.title("Rat Metadata Entry")
            leftRatNameFill, rightRatNameFill = ttk.Entry(infoPop, textvariable= self.leftRatName, width= 10), ttk.Entry(infoPop, textvariable= self.rightRatName, width= 10)
            leftRatNameLabel, rightRatNameLabel = ttk.Label(infoPop, text="Enter Left Rat #:"), ttk.Label(infoPop, text= "Enter Right Rat #:")
            leftRatNameLabel.grid(row= 1, column= 1)
            leftRatNameFill.grid(row= 1, column= 2)
            rightRatNameLabel.grid(row= 1, column= 4)
            rightRatNameFill.grid(row= 1, column= 5)
            leftRatInjTimeFill, rightRatInjTimeFill = ttk.Entry(infoPop, textvariable= self.leftRatInjectionStr, width = 10), ttk.Entry(infoPop, textvariable= self.rightRatInjectionStr, width= 10)
            leftRatInjTimeLabel, rightRatInjTimeLabel = ttk.Label(infoPop, text="Enter Left Rat Injection Trace #:"), ttk.Label(infoPop, text= "Enter Right Rat Injection Trace #:")
            leftRatInjTimeLabel.grid(row= 2, column= 1)
            leftRatInjTimeFill.grid(row= 2, column= 2)
            rightRatInjTimeLabel.grid(row= 2, column= 4)
            rightRatInjTimeFill.grid(row= 2, column= 5)

            submitButton = ttk.Button(infoPop, text="Submit", command=lambda:[onPopSubmit(), infoPop.destroy()])
            submitButton.grid(row= 3, column= 3)

    # Retrieves the baseline autofluorescence for the 4 channels analyzed and prints to a message box.
        def baselineFinder():
            try:
                pyabf.ABF(self.baselinefileName)
                subtracted470Left, subtracted405Left, subtracted470Right, subtracted405Right = acl.BaselineGet(self.baselinefileName)
                messagebox.showinfo(title= "Baselines", message= "Left - 470: %.2f 405: %.2f\nRight - 470: %.2f 405: %.2f"%(subtracted470Left, subtracted405Left, subtracted470Right, subtracted405Right))
            except FileNotFoundError:
                answer = messagebox.askretrycancel(title="Python Error", message="No baseline file found. Would you like to select a new file?", icon="error")
                if answer:
                    self.baselinefileName = ""
                    baselineFileTextUpdate()
                    self.errorStatus = True
                    return
                elif not answer:
                    self.destroy()
                    sys.exit(0)
            except PermissionError as e:
                if e.errno == 13:
                    pass
                else:
                    raise IOError("Something is wrong with this file")

    # Averages the fluorescence of all of the traces, compiles them into an excel sheet with their trace numbers, and calculates the ΔF/F
        def dataProcessorReal():
            try:
                finalSignalLeft, finalSignalRight = acl.completeProcessor(self.experimentFileName, self.baselinefileName)
            except FileNotFoundError as e:
                if str(e) == "main":
                    answer = messagebox.askretrycancel(title="Python Error", message="No main file found. Would you like to select a new file?", icon="error")
                    if answer:
                        self.experimentFileName = ""
                        chosenFileTextUpdate()
                        self.errorStatus = True
                        return
                    elif not answer:
                        self.destroy()
                        sys.exit(0)
                elif str(e) == "baseline":
                    answer = messagebox.askretrycancel(title="Python Error", message="No baseline file found. Would you like to select a new file?", icon="error")
                    if answer:
                        self.baselinefileName = ""
                        baselineFileTextUpdate()
                        self.errorStatus = True
                        return
                    elif not answer:
                        self.destroy()
                        sys.exit(0) 
        # Averages the left and right signals
            try:
                averageSignalLeft, averageSignalRight = avg.traceAverage(finalSignalLeft), avg.traceAverage(finalSignalRight)
                preInjectionAverageLeft, preInjectionAverageRight = avg.preInjectionAverage(finalSignalLeft, self.leftRatInjectionInt), avg.preInjectionAverage(finalSignalRight, self.rightRatInjectionInt)
                fluorescenceLeft, fluorescenceRight = avg.deltaF(averageSignalLeft, preInjectionAverageLeft), avg.deltaF(averageSignalRight, preInjectionAverageRight)
            except ValueError as e:
                if str(e) == "Injection traces must be larger than 1":
                    answer = messagebox.askretrycancel(title="Python Error", message="Injection trace values must be greater than 1. Please re-enter your trace number", icon="error")
                    if answer:
                        self.errorStatus = True
                        dataProcessorPop()
                        return
                    elif not answer:
                        self.destroy()
                        sys.exit(0)
                else:
                    raise
            except IndexError as e:
                if str(e) == "Injection traces are larger than total sweeps in signals":
                    answer = messagebox.askretrycancel(title="Python Error", message="Injection trace values cannot be larger than the total traces in a dataset. Please re-enter your trace number", icon="error")
                    if answer:
                        self.errorStatus = True
                        dataProcessorPop()
                        return
                    elif not answer:
                        self.destroy()
                        sys.exit(0)
                else:
                    raise
            else:
            # Analyzes peak decay, amplitude, and frequency across an entire signal containing X traces            
                self.peaksLeft, self.peaksRight = pas.wholeTracePeaks(finalSignalLeft, self.experimentFileName), pas.wholeTracePeaks(finalSignalRight, self.experimentFileName)
                injectionLeft, overviewLeft = pas.traceProcessor(self.peaksLeft)
                injectionRight, overviewRight = pas.traceProcessor(self.peaksRight)

                leftPath = os.path.join(os.getcwd(), "Processed Data", "%s Rat %s Peaks.xlsx"%(self.abfDate.strftime("%Y-%m-%d"), self.ratNameLeft))
                rightPath = os.path.join(os.getcwd(), "Processed Data", "%s Rat %s Peaks.xlsx"%(self.abfDate.strftime("%Y-%m-%d"), self.ratNameRight))
            # Saves the averaged data to an excel file with the rat's "name"
                ratDataLeft = pd.DataFrame({"Trace Number:": range(1, len(averageSignalLeft)+1), "Average Fluorescence": averageSignalLeft, 
                                            "Pre-Injection Average":preInjectionAverageLeft, "ΔF/F": fluorescenceLeft, "Bleaching Correction": None, })
                ratDataRight = pd.DataFrame({"Trace Number:": range(1, len(averageSignalRight)+1), "Average Fluorescence": averageSignalRight, 
                                            "Pre-Injection Average":preInjectionAverageRight, "ΔF/F": fluorescenceRight, "Bleaching Correction": None, })
                filenameLeft = os.path.join(os.getcwd(), "Processed Data", "%s Rat %s Temp File.xlsx"%(self.abfDate.strftime("%Y-%m-%d"), self.ratNameLeft))
                filenameRight = os.path.join(os.getcwd(), "Processed Data", "%s Rat %s Temp File.xlsx"%(self.abfDate.strftime("%Y-%m-%d"), self.ratNameRight))
                ratDataLeft.to_excel(filenameLeft, index= False)
                ratDataRight.to_excel(filenameRight, index= False)
                self.mainStatus = True
            #TODO fix FutureWarning caused by concat being empty by default.
            # Writes trace data to 2 excel files: left and right
                leftWriter, rightWriter = pd.ExcelWriter(leftPath), pd.ExcelWriter(rightPath)
                x = 1
                with leftWriter as writer:
                    # Overview Sheet
                    overviewLeft.to_excel(writer, sheet_name= "All Traces")
                    # Bins of Three
                    leftGroupThree = [list(injectionLeft.values())[i:i+3] for i in range(0, len(injectionLeft), 3)]
                    ampColumn, offTimeColumn, widthColumn, freqColumn, areaColumn = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(),
                    ampList, offTimeList, widthList, freqList, areaList = [], [], [], [], []
                    z = 1
                    for groups in leftGroupThree:
                        concat = pd.concat(groups, ignore_index=True)
                        ampList.append(pd.DataFrame(concat["Amplitude"].rename("Amplitude %i-%i"%(z, z+2)))) 
                        offTimeList.append(pd.DataFrame(concat["Off_Time_ms"].rename("Off Time %i-%i"%(z, z+2))))
                        widthList.append(pd.DataFrame(concat["Width_at50_ms"].rename("Width %i-%i"%(z, z+2))))
                        droppedFreq = pd.DataFrame(concat["Frequency"].rename("Frequency %i-%i"%(z, z+2))).dropna()
                        freqList.append(droppedFreq.reset_index(drop=True))
                        areaList.append(pd.DataFrame(concat["Area"].rename("Area %i-%i"%(z, z+2))))
                        concat.to_excel(writer, sheet_name= "Traces %i-%i"%(z, z+2), index=False)
                        z += 3
                    
                    ampColumn = pd.concat(ampList, axis="columns")
                    offTimeColumn = pd.concat(offTimeList, axis="columns")
                    widthColumn = pd.concat(widthList, axis="columns")
                    freqColumn = pd.concat(freqList, axis="columns", ignore_index=True)
                    areaColumn = pd.concat(areaList, axis="columns")

                    ampColumn.to_excel(writer, sheet_name="Amplitude")
                    offTimeColumn.to_excel(writer, sheet_name="Off Time")
                    widthColumn.to_excel(writer, sheet_name="Width")
                    freqColumn.to_excel(writer, sheet_name="Frequency")
                    areaColumn.to_excel(writer, sheet_name="Peak AUC")
                    # Individual Traces
                    for frames in injectionLeft:
                        injectionLeft[frames].to_excel(writer, sheet_name= "Trace %i"%x, index= False)
                        x += 1
                with rightWriter as writer:
                    # Overview Sheet
                    overviewRight.to_excel(writer, sheet_name= "All Traces")
                    # Bins of Three
                    rightGroupThree = [list(injectionRight.values())[i:i+3] for i in range(0, len(injectionRight), 3)]
                    ampColumn, offTimeColumn, widthColumn, freqColumn, areaColumn = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(),
                    ampList, offTimeList, widthList, freqList, areaList = [], [], [], [], []
                    z = 1
                    for groups in rightGroupThree:
                        concat = pd.concat(groups, ignore_index=True)
                        ampList.append(pd.DataFrame(concat["Amplitude"].rename("Amplitude %i-%i"%(z, z+2))))
                        offTimeList.append(pd.DataFrame(concat["Off_Time_ms"].rename("Off Time %i-%i"%(z, z+2))))
                        widthList.append(pd.DataFrame(concat["Width_at50_ms"].rename("Width %i-%i"%(z, z+2))))
                        droppedFreq = pd.DataFrame(concat["Frequency"].rename("Frequency %i-%i"%(z, z+2))).dropna()
                        freqList.append(droppedFreq.reset_index(drop=True))
                        areaList.append(pd.DataFrame(concat["Area"].rename("Area %i-%i"%(z, z+2))))
                        concat.to_excel(writer, sheet_name= "Traces %i-%i"%(z, z+2), index=False)
                        z += 3

                    ampColumn = pd.concat(ampList, axis="columns")
                    offTimeColumn = pd.concat(offTimeList, axis="columns")
                    widthColumn = pd.concat(widthList, axis="columns")
                    freqColumn = pd.concat(freqList, axis="columns", ignore_index=True)
                    areaColumn = pd.concat(areaList, axis="columns")

                    ampColumn.to_excel(writer, sheet_name="Amplitude")
                    offTimeColumn.to_excel(writer, sheet_name="Off Time")
                    widthColumn.to_excel(writer, sheet_name="Width")
                    freqColumn.to_excel(writer, sheet_name="Frequency")
                    areaColumn.to_excel(writer, sheet_name="Peak AUC")
                    
                    # Individual Traces
                    for frames in injectionRight:
                        injectionRight[frames].to_excel(writer, sheet_name= "Trace %i"%x, index= False)
                        x += 1
                self.traceStatus = True
                acl.tExport(finalSignalLeft, self.ratNameLeft, self.abfDate) #Left
                acl.tExport(finalSignalRight, self.ratNameRight, self.abfDate) #Right

        # Analyzes the peak decay, amplitude, and frequency of a single trace chosen by the user.
        def singleTracePeaks():
            finalSignalLeft, finalSignalRight = acl.completeProcessor(self.experimentFileName, self.baselinefileName)
            pas.peakDisplay(finalSignalLeft[self.trace], self.experimentFileName, "Left Rat")
            pas.peakDisplay(finalSignalRight[self.trace], self.experimentFileName, "Right Rat")

        runFileButton = ttk.Button(self, text="Process a File", command= dataProcessorPop)
        runFileButton.grid(row= 3, column= 1)
        baselineGetterButton = ttk.Button(self, text="Get the Baselines", command= baselineFinder)
        baselineGetterButton.grid(row= 3, column= 2)
        testerButton = ttk.Button(self, text="Event analysis on a single trace", command= singleTracePeaks)
        testerButton.grid(row=4, column=1)

        

def main():
    fp = tk.Tk()
    fp.title("PhotomePy 1.0")
    window = Main(fp)
    window.pack()
    fp.mainloop()

if __name__ == "__main__":
    main()