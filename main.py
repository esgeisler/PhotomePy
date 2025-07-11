import pyabf
import os
import sys
import AutoCleaner as acl
import peakAnalysis as pas
import processedSignal as pro
import TraceCleaner as trace
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import matplotlib.pyplot as plt
import yaml

class Main(tk.Frame):
    def __init__(self, master= None, **kwargs):
        super().__init__(master, **kwargs)

    # Initializes the user configuration from config.yaml
        with open("config.yaml") as c:
            userConfig = yaml.safe_load(c)
        self.norm_method = userConfig["GENERAL"]["normalization"]
        self.bleach_correct = userConfig["GENERAL"]["bleach_correction"]
        self.exportStatus = userConfig["GENERAL"]["export_abf"]
        self.movement_corr = userConfig["GENERAL"]["movement_corr_calc"]
        self.peakMethod = userConfig["EVENT_HANDLING"]["peak_id_method"]
        self.peakThreshold = userConfig["EVENT_HANDLING"]["peak_detection_threshold"]
        self.peakWindow = userConfig["EVENT_HANDLING"]["peak_window"]
        self.peakMax = userConfig["EVENT_HANDLING"]["peak_top"]
        self.peakMin = userConfig["EVENT_HANDLING"]["peak_bottom"]

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
        self.controlStatus = tk.IntVar()

    # Dropdown menu for selecting the trace used in SingleTracePeaks
        traceDrop = ttk.Combobox(self,  state= 'readonly', width= 9, textvariable= self.dropValue)
        traceDrop['values'] = []
        traceDrop.grid(row= 4, column=2)

    # Sets the values inside the dropdown menu, and sets it to update when a main file is selected
        def traceSelector(trace):
            self.trace = (int(traceDrop.get()) - 1)
        traceDrop.bind("<<ComboboxSelected>>", traceSelector)
        def dropdownUpdater():
            traceDrop['values'] = self.options

    # Opens the main file, containing data from the session with 2 rats
        def fileBrowserExperiment():
            t = False
            while not t:
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
                        t = True
                    else:
                        raise IOError("Something is wrong with this file")
                else:
                    chosenFileDisplay.insert(tk.END, self.experimentFileName)
                    return self.experimentFileName
    
    # Opens the baseline file containing the baseline autofluorescence
        def fileBrowserBaseline():
            t = False
            while not t:
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
                        t = True
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
    
        def singleTracePopSubmit():
            self.ratNameLeft, self.ratNameRight = self.leftRatName.get(), self.rightRatName.get()
            try:
                self.leftRatInjectionInt, self.rightRatInjectionInt = (int(float(self.leftRatInjectionStr.get())) - 1), (int(float(self.rightRatInjectionStr.get())) - 1)
            except ValueError:
                answer = messagebox.askretrycancel(title="Python Error", message="Injection traces must be a number. Please re-enter your trace number", icon="error")
                if answer:
                    self.errorStatus = True
                    self.leftRatInjectionStr, self.rightRatInjectionStr = tk.StringVar(), tk.StringVar()
                    singleTracePop()
                    return
                elif not answer:
                    self.destroy()
                    sys.exit(0)
            singleTracePeaks()

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

            submitButton = ttk.Button(infoPop, text="Submit", command=lambda:[infoPop.destroy(), onPopSubmit()])
            submitButton.grid(row= 3, column= 3)
            curveFitCheckBox = ttk.Checkbutton(infoPop, text="Would you like to fit a curve for bleaching correction? \nWARNING: Only check if this is a vehicle control.", variable=self.controlStatus)
            curveFitCheckBox.grid(row=4, column=3)

        def singleTracePop():
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

            submitButton = ttk.Button(infoPop, text="Submit", command=lambda:[infoPop.destroy(), singleTracePopSubmit()])
            submitButton.grid(row= 3, column= 3)

    # Retrieves the baseline autofluorescence for the 4 channels analyzed and prints to a message box.
        def baselineFinder():
            try:
                pyabf.ABF(self.baselinefileName)
                subtracted470Left, subtracted405Left, subtracted470Right, subtracted405Right = acl.baselineGet(self.baselinefileName)
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
            while True:
                try:
                    rawTrace = trace.TraceCleaner(self.experimentFileName, self.baselinefileName)
                    match self.bleach_correct:
                        case "isosbestic":
                            finalSignalLeft, finalSignalRight = rawTrace.completeProcessor()
                        case "biexp": 
                            finalSignalLeft, finalSignalRight, traceDecayFunction, decayRSquared = rawTrace.newCompleteProcessor(self.controlStatus.get(), self.ratNameLeft, self.ratNameRight, self.abfDate)
                        case "biexpActive":
                            finalSignalLeft, finalSignalRight = rawTrace.activeSignalDecayProcessor()
                        case "refReg":
                            finalSignalLeft, finalSignalRight, _, _, = rawTrace.referenceSignalDecayProcessor()
                        case _:
                            raise ValueError("bleach_correct error")
                    leftSignal = pro.ProcessedTotalSignal(self.experimentFileName, finalSignalLeft, self.ratNameLeft, self.leftRatInjectionInt)
                    rightSignal = pro.ProcessedTotalSignal(self.experimentFileName, finalSignalRight, self.ratNameRight, self.rightRatInjectionInt)
                    match self.norm_method:
                        case "SBR":
                            break
                        case "zSBR":
                            leftSignal.processedSignalArray = (leftSignal.processedSignalArray - leftSignal.preInjAvg) / leftSignal.preInjStdev
                            rightSignal.processedSignalArray = (rightSignal.processedSignalArray - rightSignal.preInjAvg) / rightSignal.preInjStdev
                            break
                        case _:
                            raise ValueError("norm error")
                except FileNotFoundError as e:
                    match str(e):
                        case "main":
                            answer = messagebox.askretrycancel(title="Python Error", message="No main file found. Would you like to select a new file?", icon="error")
                            if answer:
                                self.experimentFileName = ""
                                chosenFileTextUpdate()
                                self.errorStatus = True
                                return
                            elif not answer:
                                self.destroy()
                                sys.exit(0)
                        case "baseline":
                            answer = messagebox.askretrycancel(title="Python Error", message="No baseline file found. Would you like to select a new file?", icon="error")
                            if answer:
                                self.baselinefileName = ""
                                baselineFileTextUpdate()
                                self.errorStatus = True
                                return
                            elif not answer:
                                self.destroy()
                                sys.exit(0)
                        case _:
                            raise
                except ValueError as e:
                    match str(e):
                        case "bleach_correct error":
                            print("Values for bleach_correct must be either 'isosbestic' or 'biexp'.\nUsing default value (isosbestic).")
                            self.bleach_correct = "isosbestic"
                        case "norm error":
                            print("Values for normalization must be either 'SBR' or 'zSBR'.\nUsing default value (SBR).")
                            self.norm_method = "SBR"
                        case "Injection traces must be larger than 1":
                            answer = messagebox.askretrycancel(title="Python Error", message="Injection trace values must be greater than 1. Please re-enter your trace number", icon="error")
                            if answer:
                                self.errorStatus = True
                                dataProcessorPop()
                                return
                            elif not answer:
                                self.destroy()
                                sys.exit(0)
                        case _:
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
        # Analyzes peak decay, amplitude, and frequency across an entire signal containing X traces           
            self.peaksLeft, self.peaksRight = pas.wholeTracePeaks(leftSignal.processedSignalArray, self.experimentFileName), pas.wholeTracePeaks(rightSignal.processedSignalArray, self.experimentFileName)
            injectionLeft, overviewLeft = pas.traceProcessor(self.peaksLeft)
            injectionRight, overviewRight = pas.traceProcessor(self.peaksRight)
            match self.bleach_correct:
                case "isosbestic":
                    leftPath = os.path.join(os.getcwd(), "Processed Data", "Default %s Rat %s Peaks.xlsx"%(leftSignal.date, leftSignal.ratID))
                    rightPath = os.path.join(os.getcwd(), "Processed Data", "Default %s Rat %s Peaks.xlsx"%(rightSignal.date, rightSignal.ratID))
                case "refReg":
                    leftPath = os.path.join(os.getcwd(), "Processed Data", "refReg %s Rat %s Peaks.xlsx"%(leftSignal.date, leftSignal.ratID))
                    rightPath = os.path.join(os.getcwd(), "Processed Data", "refReg %s Rat %s Peaks.xlsx"%(rightSignal.date, rightSignal.ratID))
                case "biexpActive":
                    leftPath = os.path.join(os.getcwd(), "Processed Data", "biexpActive %s Rat %s Peaks.xlsx"%(leftSignal.date, leftSignal.ratID))
                    rightPath = os.path.join(os.getcwd(), "Processed Data", "biexpActive %s Rat %s Peaks.xlsx"%(rightSignal.date, rightSignal.ratID))
                case _:
                    raise ValueError("bleach error")
        # Saves the averaged data to an excel file with the rat's "name"
            match self.bleach_correct:
                case "isosbestic":
                    filenameLeft = os.path.join(os.getcwd(), "Processed Data", "%s Rat %s Overall Fluorescence.xlsx"%(leftSignal.date, leftSignal.ratID))
                    filenameRight = os.path.join(os.getcwd(), "Processed Data", "%s Rat %s Overall Fluorescence.xlsx"%(rightSignal.date, rightSignal.ratID))
                case "refReg":
                    filenameLeft = os.path.join(os.getcwd(), "Processed Data", "refReg %s Rat %s Overall Fluorescence.xlsx"%(leftSignal.date, leftSignal.ratID))
                    filenameRight = os.path.join(os.getcwd(), "Processed Data", "refReg %s Rat %s Overall Fluorescence.xlsx"%(rightSignal.date, rightSignal.ratID))
                case "biexpActive":
                    filenameLeft = os.path.join(os.getcwd(), "Processed Data", "biexpActive %s Rat %s Overall Fluorescence.xlsx"%(leftSignal.date, leftSignal.ratID))
                    filenameRight = os.path.join(os.getcwd(), "Processed Data", "biexpActive %s Rat %s Overall Fluorescence.xlsx"%(rightSignal.date, rightSignal.ratID))
                case _:
                    raise ValueError("bleach error")
            leftOverviewWriter, rightOverviewWriter = pd.ExcelWriter(filenameLeft), pd.ExcelWriter(filenameRight)
            leftOrRight = [leftOverviewWriter, rightOverviewWriter]
            for rats in leftOrRight:
                with rats as writer:
                    if rats == leftOverviewWriter:
                        ratData = pd.DataFrame({"Trace Number:": range(1, len(leftSignal.wholeTraceAverages)+1), 
                                                "Mean Fluorescence": leftSignal.wholeTraceAverages, "Median Fluorescence": leftSignal.wholeTraceMedians,
                                                "Pre-Injection Average": leftSignal.preInjAvg, "ΔF/F": leftSignal.normFluorescence, 
                                                "Bleaching Correction": None, "Z-Score": leftSignal.overallZScore
                                                })
                    elif rats == rightOverviewWriter:
                        ratData = pd.DataFrame({"Trace Number:": range(1, len(rightSignal.wholeTraceAverages)+1), 
                                                "Mean Fluorescence": rightSignal.wholeTraceAverages, "Median Fluorescence": rightSignal.wholeTraceMedians,
                                                "Pre-Injection Average": rightSignal.preInjAvg, "ΔF/F": rightSignal.normFluorescence, 
                                                "Bleaching Correction": None, "Z-Score": rightSignal.overallZScore
                                                })
                    ratData.to_excel(writer, index= False, sheet_name="Overview")
                    if self.controlStatus.get() == 1:
                        if rats == leftOverviewWriter:
                            bleach405 = traceDecayFunction[0]
                            bleach470 = traceDecayFunction[1]
                            rsquared405 = decayRSquared[0]
                            rsquared470 = decayRSquared[1]
                        elif rats == rightOverviewWriter:
                            bleach405 = traceDecayFunction[2]
                            bleach470 = traceDecayFunction[3]
                            rsquared405 = decayRSquared[2]
                            rsquared470 = decayRSquared[3]
                        bleachFit = pd.DataFrame({"Trace Number:": range(1, len(bleach405) + 1), 
                                                    "Bleaching Correction (405):": bleach405, "Goodness of fit(R^2):": rsquared405, 
                                                    "Bleaching Correction (470):": bleach470, "Goodness of fit(R^2):": rsquared470
                                                    })
                        bleachFit.to_excel(writer, index=False, sheet_name="Bleaching Correction")
            self.mainStatus = True
        #TODO fix FutureWarning caused by concat being empty by default.
        # Writes trace data to 2 excel files: left and right
            leftWriter, rightWriter = pd.ExcelWriter(leftPath), pd.ExcelWriter(rightPath)
            leftGroupThree, rightGroupThree = [list(injectionLeft.values())[i:i+3] for i in range(0, len(injectionLeft), 3)], [list(injectionRight.values())[i:i+3] for i in range(0, len(injectionRight), 3)]
            peakWriters = [leftWriter, rightWriter]
            for rats in peakWriters:
                with rats as writer:
                    ampColumn, absAmpColumn, offTimeColumn, widthColumn, freqColumn, areaColumn, totAreaColumn, riseColumn, decayColumn = (pd.DataFrame() for _ in range(9))
                    ampList, absAmpList, offTimeList, widthList, freqList, areaList, totAreaList, riseList, decayList = ([] for _ in range(9))
                    z = 1
                    if rats == leftWriter:
                        # Overview Sheet
                        overviewLeft.to_excel(writer, sheet_name= "All Traces")
                        # Bins of Three
                        for groups in leftGroupThree:
                            concat = pd.concat(groups, ignore_index=True)
                            ampList.append(pd.DataFrame(concat["Amplitude"].rename("Amplitude %i-%i"%(z, z+2))))
                            absAmpList.append(pd.DataFrame(concat["Abs_Amplitude"].rename("Absolute Amplitude %i-%i"%(z, z+2))))
                            offTimeList.append(pd.DataFrame(concat["Off_Time_ms"].rename("Off Time %i-%i"%(z, z+2))))
                            widthList.append(pd.DataFrame(concat["Width_at50_ms"].rename("Width %i-%i"%(z, z+2))))
                            droppedFreq = pd.DataFrame(concat["Frequency"].rename("Frequency %i-%i"%(z, z+2))).dropna()
                            freqList.append(droppedFreq.reset_index(drop=True))
                            areaList.append(pd.DataFrame(concat["Avg_Area"].rename("Area %i-%i"%(z, z+2))))
                            droppedArea = pd.DataFrame(concat["Total_Area"].rename("Total Area %i-%i"%(z, z+2))).dropna()
                            totAreaList.append(droppedArea.reset_index(drop=True))
                            riseList.append(pd.DataFrame(concat["Rise_Tau"].rename("Rise Tau %i-%i"%(z, z+2))))
                            decayList.append(pd.DataFrame(concat["Decay_Tau"].rename("Decay Tau %i-%i"%(z, z+2))))
                            concat.to_excel(writer, sheet_name= "Traces %i-%i"%(z, z+2), index=False)
                            z += 3
                    elif rats == rightWriter:
                        # Overview Sheet
                        overviewRight.to_excel(writer, sheet_name= "All Traces")
                        # Bins of Three
                        for groups in rightGroupThree:
                            concat = pd.concat(groups, ignore_index=True)
                            ampList.append(pd.DataFrame(concat["Amplitude"].rename("Amplitude %i-%i"%(z, z+2))))
                            absAmpList.append(pd.DataFrame(concat["Abs_Amplitude"].rename("Absolute Amplitude %i-%i"%(z, z+2))))
                            offTimeList.append(pd.DataFrame(concat["Off_Time_ms"].rename("Off Time %i-%i"%(z, z+2))))
                            widthList.append(pd.DataFrame(concat["Width_at50_ms"].rename("Width %i-%i"%(z, z+2))))
                            droppedFreq = pd.DataFrame(concat["Frequency"].rename("Frequency %i-%i"%(z, z+2))).dropna()
                            freqList.append(droppedFreq.reset_index(drop=True))
                            areaList.append(pd.DataFrame(concat["Avg_Area"].rename("Mean Area %i-%i"%(z, z+2))))
                            droppedArea = pd.DataFrame(concat["Total_Area"].rename("Total Area %i-%i"%(z, z+2))).dropna()
                            totAreaList.append(droppedArea.reset_index(drop=True))
                            riseList.append(pd.DataFrame(concat["Rise_Tau"].rename("Rise Tau %i-%i"%(z, z+2))))
                            decayList.append(pd.DataFrame(concat["Decay_Tau"].rename("Decay Tau %i-%i"%(z, z+2))))
                            concat.to_excel(writer, sheet_name= "Traces %i-%i"%(z, z+2), index=False)
                            z += 3
                    
                    ampColumn = pd.concat(ampList, axis="columns")
                    ampColumn.to_excel(writer, sheet_name="Amplitude")

                    absAmpColumn = pd.concat(absAmpList, axis="columns")
                    absAmpColumn.to_excel(writer, sheet_name="Abs. Amplitude")

                    offTimeColumn = pd.concat(offTimeList, axis="columns")
                    offTimeColumn.to_excel(writer, sheet_name="Off Time")

                    widthColumn = pd.concat(widthList, axis="columns")
                    widthColumn.to_excel(writer, sheet_name="Width")

                    freqColumn = pd.concat(freqList, axis="columns", ignore_index=True)
                    freqColumn.to_excel(writer, sheet_name="Frequency")

                    areaColumn = pd.concat(areaList, axis="columns")
                    areaColumn.to_excel(writer, sheet_name="Mean Area")

                    totAreaColumn = pd.concat(totAreaList, axis="columns", ignore_index=True)
                    totAreaColumn.to_excel(writer, sheet_name="Total Area")

                    riseColumn = pd.concat(riseList, axis="columns")
                    riseColumn.to_excel(writer, sheet_name="Rise Tau")

                    decayColumn = pd.concat(decayList, axis="columns")
                    decayColumn.to_excel(writer, sheet_name="Decay Tau")
        
            self.traceStatus = True
            if self.exportStatus:
                acl.tExport(leftSignal.processedSignalArray, leftSignal.ratID, leftSignal.date) #Left
                acl.tExport(rightSignal.processedSignalArray, rightSignal.ratID, rightSignal.date) #Right 

        # Analyzes the peak decay, amplitude, and frequency of a single trace chosen by the user.
        def singleTracePeaks():
            try:
                rawTrace = trace.TraceCleaner(self.experimentFileName, self.baselinefileName)
                match self.bleach_correct:
                    case "isosbestic":
                        finalSignalLeft, finalSignalRight = rawTrace.completeProcessor()
                    case "biexpActive":
                        finalSignalLeft, finalSignalRight = rawTrace.activeSignalDecayProcessor()
                    case "refReg":
                        finalSignalLeft, finalSignalRight, decayFitLeft, decayFitRight = rawTrace.referenceSignalDecayProcessor()
                leftSignal = pro.ProcessedTotalSignal(self.experimentFileName, finalSignalLeft, self.ratNameLeft, self.leftRatInjectionInt) #TEMP INJECTION VALUE FOR TESTING ONLY
                rightSignal = pro.ProcessedTotalSignal(self.experimentFileName, finalSignalRight, self.ratNameRight, self.rightRatInjectionInt) #TEMP INJECTION VALUE FOR TESTING ONLY
                if self.norm_method == "zSBR":
                    leftSignal.processedSignalArray = (leftSignal.processedSignalArray - leftSignal.preInjAvg) / leftSignal.preInjStdev
                    rightSignal.processedSignalArray = (rightSignal.processedSignalArray - rightSignal.preInjAvg) / rightSignal.preInjStdev
                else:
                    pass
                #finalSignalLeft, finalSignalRight, _, _ = acl.newCompleteProcessor(self.experimentFileName, self.baselinefileName, self.controlStatus.get(), "21", "22", self.abfDate)
                if self.movement_corr: 
                    acl.isoLinRegPlot(self.experimentFileName, 1, 0, self.trace, "Left Rat Motion Correlation")
                    acl.isoLinRegPlot(self.experimentFileName, 5, 4, self.trace, "Right Rat Motion Correlation")
                pas.peakDisplay(leftSignal.processedSignalArray, self.experimentFileName, "Left Rat", self.trace)
                pas.peakDisplay(rightSignal.processedSignalArray, self.experimentFileName, "Right Rat", self.trace)
            except FileNotFoundError as e:
                match str(e):
                    case "main":
                        answer = messagebox.askretrycancel(title="Python Error", message="No main file found. Would you like to select a new file?", icon="error")
                        if answer:
                            self.experimentFileName = ""
                            chosenFileTextUpdate()
                            self.errorStatus = True
                            return
                        elif not answer:
                            self.destroy()
                            sys.exit(0)
                    case "baseline":
                        answer = messagebox.askretrycancel(title="Python Error", message="No baseline file found. Would you like to select a new file?", icon="error")
                        if answer:
                            self.baselinefileName = ""
                            baselineFileTextUpdate()
                            self.errorStatus = True
                            return
                        elif not answer:
                            self.destroy()
                            sys.exit(0)
                    case _:
                        raise  

        #TODO: Adjust the xlim and ylim parameters to be dynamic in relation to the shape of the trace, rather than hard-coded.
        def stepPrinterSingleTrace():
            try:
                leftNoBaseline, rightNoBaseline, leftFiltered, rightFiltered, combinedLeft, combinedRight, finalLeft, finalRight = acl.stepwiseProcessor(self.experimentFileName, self.baselinefileName)
                seconds, _, _ = pas.secondsCalculator(self.experimentFileName)
                # Shows the left rat's raw traces, with no corrections and after subtracting baseline autofluorescence.
                peakFig, ax1 = plt.subplots()

                ax1.set_ylim([leftNoBaseline[0][self.trace][0] - 2, leftNoBaseline[0][self.trace][0] + 1])
                ax1.set_xlim([0, 47750])
                ax1.set_xlabel("Time (s)", fontsize="large")
                ax1.set_xticks(ticks=seconds, labels=range(len(seconds)))
                ax1.set_ylabel("Active Signal (V)", color="g", fontsize="large")
                ax1.autoscale(False)
                ax1.plot(leftNoBaseline[0][self.trace], "g")
                
                ax2 = ax1.twinx()
                ax2.set_ylabel("Isosobestic Signal (V)", color = "b", fontsize="large")
                ax2.autoscale(False)
                ax2.set_ylim([leftNoBaseline[1][self.trace][0] - 1, leftNoBaseline[1][self.trace][0] + 2])
                ax2.plot(leftNoBaseline[1][self.trace], "b")
                
                peakFig.suptitle("Raw Left Rat", fontsize="large")
                peakFig.tight_layout()
                peakFig.set_size_inches(9,3)
                plt.show()

                # Shows the right rat's raw traces, with no corrections and after subtracting baseline autofluorescence.
                peakFig, ax1 = plt.subplots()
                ax1.autoscale(False)
                ax1.set_ylim([rightNoBaseline[0][self.trace][0] - 2, rightNoBaseline[0][self.trace][0] + 1])
                ax1.set_xlim([0, 47750])
                ax1.set_xlabel("Time (s)", fontsize="large")
                ax1.set_xticks(ticks=seconds, labels=range(len(seconds)))
                ax1.set_ylabel("Active Signal (V)", color="g", fontsize="large")
                ax1.plot(rightNoBaseline[0][self.trace], "g")

                ax2 = ax1.twinx()
                ax2.set_ylabel("Isosobestic Signal (V)", color = "b", fontsize="large")
                ax2.autoscale(False)
                ax2.set_ylim([rightNoBaseline[1][self.trace][0] - 0.5, rightNoBaseline[1][self.trace][0] + 0.5])
                ax2.plot(rightNoBaseline[1][self.trace], "b")
                peakFig.suptitle("Raw Right Rat", fontsize="large")
                peakFig.tight_layout()
                peakFig.set_size_inches(9,3)
                plt.show()

                # Shows the left rat's traces after filtering
                peakFig, ax1 = plt.subplots()
                ax1.autoscale(False)
                ax1.set_ylim([leftFiltered[0][self.trace][0] - 2, leftFiltered[0][self.trace][0] + 1])
                ax1.set_xlim([0, 47750])
                ax1.set_xlabel("Time (s)", fontsize="large")
                ax1.set_xticks(ticks=seconds, labels=range(len(seconds)))
                ax1.set_ylabel("Active Signal (V)", color="g", fontsize="large")
                ax1.plot(leftFiltered[0][self.trace], "g")

                ax2 = ax1.twinx()
                ax2.set_ylabel("Isosobestic Signal (V)", color = "b")
                ax2.autoscale(False)
                ax2.set_ylim(leftFiltered[1][self.trace][0] - 1, leftFiltered[1][self.trace][0] + 2)
                ax2.plot(leftFiltered[1][self.trace], "b")
                peakFig.suptitle("Left Rat (Filtered)", fontsize="large")
                peakFig.tight_layout()
                peakFig.set_size_inches(9,3)
                plt.show()

                # Shows the right rat's traces after filtering
                peakFig, ax1 = plt.subplots()
                ax1.autoscale(False)
                ax1.set_ylim([rightFiltered[0][self.trace][0] - 2, rightFiltered[0][self.trace][0] + 1])
                ax1.set_xlim([0, 47750])
                ax1.set_xlabel("Time (s)")
                ax1.set_xticks(ticks=seconds, labels=range(len(seconds)))
                ax1.set_ylabel("Active Signal (V)", color="g")
                ax1.plot(rightFiltered[0][self.trace], "g")

                ax2 = ax1.twinx()
                ax2.set_ylabel("Isosobestic Signal (V)", color = "b", fontsize="large")
                ax2.autoscale(False)
                ax2.set_ylim(rightFiltered[1][self.trace][0] - 0.5, rightFiltered[1][self.trace][0] + 0.5)
                ax2.plot(rightFiltered[1][self.trace], "b")
                peakFig.suptitle("Right Rat (Filtered)", fontsize="large")
                peakFig.tight_layout()
                peakFig.set_size_inches(9,3)
                plt.show()

                # Shows the combined traces of the left rat
                peakFig, ax1 = plt.subplots()
                ax1.autoscale(False)
                ax1.set_ylim([combinedLeft[self.trace][0] - 0.5, combinedLeft[self.trace][0] + 0.5])
                ax1.set_xlim([0, 47750])
                ax1.set_ylabel("Normalized Fluorescence (dF/F)", fontsize="large")
                ax1.set_xlabel("Time (s)", fontsize="large")
                ax1.set_xticks(ticks=seconds, labels=range(len(seconds)))
                ax1.plot(combinedLeft[self.trace], "b")
                peakFig.suptitle("Left Rat (dF/F)", fontsize="large")
                peakFig.tight_layout()
                peakFig.set_size_inches(9,3)
                plt.show()

                # Shows the combined traces of the right rat
                peakFig, ax1 = plt.subplots()
                ax1.autoscale(False)
                ax1.set_ylim([combinedRight[self.trace][0] - 0.25, combinedRight[self.trace][0] + 0.4])
                ax1.set_xlim([0, 47750])
                ax1.set_ylabel("Normalized Fluorescence (dF/F)", fontsize="large")
                ax1.set_xlabel("Time (s)", fontsize="large")
                ax1.set_xticks(ticks=seconds, labels=range(len(seconds)))
                ax1.plot(combinedRight[self.trace], "b")
                peakFig.suptitle("Right Rat (dF/F)", fontsize="large")
                peakFig.tight_layout()
                peakFig.set_size_inches(9,3)
                plt.show()      

                # Shows the final signal of the left rat
                peakFig, ax1 = plt.subplots()
                ax1.autoscale(False)
                ax1.set_ylim([finalLeft[self.trace][0] - 0.5, finalLeft[self.trace][0] + 0.5])
                ax1.set_xlim([0, 47750])
                ax1.set_ylabel("Normalized Fluorescence (dF/F)", fontsize="large")
                ax1.set_xlabel("Time (s)", fontsize="large")
                ax1.set_xticks(ticks=seconds, labels=range(len(seconds)))
                ax1.plot(finalLeft[self.trace], "b")
                peakFig.suptitle("Left Rat Final Signal", fontsize="large")
                peakFig.tight_layout()
                peakFig.set_size_inches(9,3)
                plt.show()

                # Shows the final signal of the right rat
                peakFig, ax1 = plt.subplots()
                ax1.autoscale(False)
                ax1.set_ylim([finalRight[self.trace][0] - 0.25, finalRight[self.trace][0] + 0.4])
                ax1.set_xlim([0, 47750])
                ax1.set_ylabel("Normalized Fluorescence (dF/F)", fontsize="large")
                ax1.set_xlabel("Time (s)", fontsize="large")
                ax1.set_xticks(ticks=seconds, labels=range(len(seconds)))
                ax1.plot(finalRight[self.trace], "b")
                peakFig.suptitle("Right Rat Final Signal", fontsize="large")
                peakFig.tight_layout()
                peakFig.set_size_inches(9,3)
                plt.show()   
            except FileNotFoundError as e:
                match str(e):
                    case "main":
                        answer = messagebox.askretrycancel(title="Python Error", message="No main file found. Would you like to select a new file?", icon="error")
                        if answer:
                            self.experimentFileName = ""
                            chosenFileTextUpdate()
                            self.errorStatus = True
                            return
                        elif not answer:
                            self.destroy()
                            sys.exit(0)
                    case "baseline":
                        answer = messagebox.askretrycancel(title="Python Error", message="No baseline file found. Would you like to select a new file?", icon="error")
                        if answer:
                            self.baselinefileName = ""
                            baselineFileTextUpdate()
                            self.errorStatus = True
                            return
                        elif not answer:
                            self.destroy()
                            sys.exit(0)
                    case _:
                        raise  
            

        runFileButton = ttk.Button(self, text="Process all Traces", command= dataProcessorPop)
        runFileButton.grid(row= 3, column= 1)
        baselineGetterButton = ttk.Button(self, text="Get Baseline Autofluorescence", command= baselineFinder)
        baselineGetterButton.grid(row= 3, column= 2)
        testerButton = ttk.Button(self, text="Process a Single Trace", command= singleTracePop)
        testerButton.grid(row=4, column=1)
        correctionTroubleshootButton = ttk.Button(self, text="Troubleshoot Corrections", command= stepPrinterSingleTrace)
        correctionTroubleshootButton.grid(row=5, column=1)

        

def main():
    fp = tk.Tk()
    fp.title("PhotomePy 1.0")
    window = Main(fp)
    window.pack()
    fp.mainloop()
    

if __name__ == "__main__":
    main()