import AutoCleaner as acl
import pyabf
import os
import AverageTraces as avg
import peakAnalysis as pas
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import scipy.signal as sci

#TODO Peak indices should be summarized over time (peak 2 should be time in trace + 30 sec, peak 3 should be time in trace + 1 minute, etc)
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
        self.ratNameLeft = tk.StringVar()
        self.ratNameRight = tk.StringVar()
        self.ratInjectionLeft = 0
        self.ratInjectionRight = 0
        self.peaksLeft = []
        self.peaksRight = []
        self.options = []
        self.trace = 0
        self.abfDate = tk.StringVar()
        self.statusCheck = True
        self.mainStatus = False
        self.traceStatus = False

    # Dropdown menu for selecting the trace used in SingleTracePeaks
        traceDrop = ttk.Combobox(self,  state= 'readonly', width= 9, textvariable= self.dropValue)
        traceDrop['values'] = []
        traceDrop.grid(row= 4, column=2)

    # Sets the values inside the dropdown menu, and sets it to update when a main file is selected
        def traceSelector(event):
            self.trace = (int(traceDrop.get()) - 1)
        traceDrop.bind("<<ComboboxSelected>>", traceSelector)
        def dropdownUpdater():
            traceDrop['values'] = self.options

    # Opens the main file, containing data from the session with 2 rats
        def fileBrowserExperiment():
            self.experimentFileName = filedialog.askopenfilename(initialdir= os.path.join(os.getcwd(), "Raw Data"), title= "Select a Main File", 
                                                                 filetypes=(("Axon Binary Fles", "*.abf*"), ("All Files," "*.*")))
            chosenFileDisplay.insert(tk.END, self.experimentFileName)
            abf = pyabf.ABF(self.experimentFileName)
            self.options = [str(x + 1) for x in abf.sweepList]
            self.abfDate = abf.abfDateTime
            return self.experimentFileName
    
    # Opens the baseline file containing the baseline autofluorescence
        def fileBrowserBaseline():
            self.baselinefileName = filedialog.askopenfilename(initialdir= os.path.join(os.getcwd(), "Raw Data"), title= "Select a Main File", filetypes=(("Axon Binary Fles", "*.abf*"), ("All Files," "*.*")))
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
            self.ratNameLeft = self.leftRatName.get()
            self.ratNameRight = self.rightRatName.get()
            self.ratInjectionLeft = int(self.leftRatInjection.get()) - 1
            self.ratInjectionRight = int(self.rightRatInjection.get()) - 1

    # Checks the status of the analysis to determine whether to display a messagebox message.
        def statusChecker():
            while self.statusCheck:
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
            leftRatNameFill = ttk.Entry(infoPop, textvariable= self.leftRatName, width= 10)
            leftRatNameLabel = ttk.Label(infoPop, text="Enter Left Rat #:")
            rightRatNameFill = ttk.Entry(infoPop, textvariable= self.rightRatName, width= 10)
            rightRatNameLabel = ttk.Label(infoPop, text= "Enter Right Rat #:")
            leftRatNameLabel.grid(row= 1, column= 1)
            leftRatNameFill.grid(row= 1, column= 2)
            rightRatNameLabel.grid(row= 1, column= 4)
            rightRatNameFill.grid(row= 1, column= 5)
            leftRatInjTimeFill = ttk.Entry(infoPop, textvariable= self.leftRatInjection, width = 10)
            leftRatInjTimeLabel = ttk.Label(infoPop, text="Enter Left Rat Injection Trace #:")
            rightRatInjTimeFill = ttk.Entry(infoPop, textvariable= self.rightRatInjection, width= 10)
            rightRatInjTimeLabel = ttk.Label(infoPop, text= "Enter Right Rat Injection Trace #:")
            leftRatInjTimeLabel.grid(row= 2, column= 1)
            leftRatInjTimeFill.grid(row= 2, column= 2)
            rightRatInjTimeLabel.grid(row= 2, column= 4)
            rightRatInjTimeFill.grid(row= 2, column= 5)


            submitButton = ttk.Button(infoPop, text="Submit", command=lambda:[onPopSubmit(), infoPop.destroy(), dataProcessorReal(), statusChecker()])
            submitButton.grid(row= 3, column= 3)

    # Retrieves the baseline autofluorescence for the 4 channels analyzed and prints to a message box.
        def baselineFinder():
            pyabf.ABF(self.baselinefileName)
            subtracted470Left, subtracted405Left, subtracted470Right, subtracted405Right = acl.BaselineGet(self.baselinefileName)
            messagebox.showinfo(title= "Baselines", message= "Left - 470: %.2f 405: %.2f\nRight - 470: %.2f 405: %.2f"%(subtracted470Left, subtracted405Left, subtracted470Right, subtracted405Right))

    # Averages the fluorescence of all of the traces, compiles them into an excel sheet with their trace numbers, and calculates the ΔF/F
        def dataProcessorReal():
            finalSignalLeft, finalSignalRight = acl.completeProcessor(self.experimentFileName, self.baselinefileName)
        # Averages the left and right signals
            averageSignalLeft = avg.traceAverage(finalSignalLeft)
            preInjectionAverageLeft = avg.preInjectionAverage(finalSignalLeft, self.ratInjectionLeft)
            fluorescenceLeft = avg.deltaF(averageSignalLeft, preInjectionAverageLeft)

            averageSignalRight = avg.traceAverage(finalSignalRight)
            preInjectionAverageRight = avg.preInjectionAverage(finalSignalRight, self.ratInjectionRight)
            fluorescenceRight = avg.deltaF(averageSignalRight, preInjectionAverageRight)
        # Analyzes peak decay, amplitude, and frequency across an entire signal containing X traces            
            self.peaksLeft = pas.wholeTracePeaks(finalSignalLeft, self.experimentFileName)
            self.peaksRight = pas.wholeTracePeaks(finalSignalRight, self.experimentFileName)

            preInjectionLeft, postInjectionLeft, preOverviewLeft, postOverviewLeft = pas.traceProcessor(self.peaksLeft, self.ratInjectionLeft)
            preInjectionRight, postInjectionRight, preOverviewRight, postOverviewRight = pas.traceProcessor(self.peaksRight, self.ratInjectionRight)

            preLeft = os.path.join(os.getcwd(), "Processed Data", "%s Rat %s Pre-Injection Peaks.xlsx"%(self.abfDate.strftime("%Y-%m-%d"), self.ratNameLeft))
            postLeft = os.path.join(os.getcwd(), "Processed Data", "%s Rat %s Post-Injection Peaks.xlsx"%(self.abfDate.strftime("%Y-%m-%d"), self.ratNameLeft))

            preRight = os.path.join(os.getcwd(), "Processed Data", "%s Rat %s Pre-Injection Peaks.xlsx"%(self.abfDate.strftime("%Y-%m-%d"), self.ratNameRight))
            postRight = os.path.join(os.getcwd(), "Processed Data", "%s Rat %s Post-Injection Peaks.xlsx"%(self.abfDate.strftime("%Y-%m-%d"), self.ratNameRight))
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
        # Writes trace data to 2 excel files: Pre-injection and post-injection
            preLeftWriter = pd.ExcelWriter(preLeft)
            postLeftWriter = pd.ExcelWriter(postLeft)
            preRightWriter = pd.ExcelWriter(preRight)
            postRightWriter = pd.ExcelWriter(postRight)
            x = 1
            with preLeftWriter as writer:
                # Overview Sheet
                preOverviewLeft.to_excel(writer, sheet_name= "All Traces")
                # Bins of Three
                preLeftGroupThree = [list(preInjectionLeft.values())[i:i+3] for i in range(0, len(preInjectionLeft), 3)]
                ampColumn = pd.DataFrame()
                offTimeColumn = pd.DataFrame()
                widthColumn = pd.DataFrame()
                freqColumn = pd.DataFrame()
                ampList = []
                offTimeList = []
                widthList = []
                freqList = []
                z = 1
                for groups in preLeftGroupThree:
                    concat = pd.DataFrame(columns= ['Event_Num', 'Peak_Index', 
                                        'Peak_Time_Sec', 'Event_Window_Start', 
                                        'Event_Window_End', 'Amplitude', 'Off_Time_ms',
                                        'Width_at50_ms','Frequency'])
                    for y in groups:
                        concat = pd.concat([concat, y if not y.empty else None], ignore_index=True)
                    ampList.append(concat["Amplitude"].rename("Amplitude %i-%i"%(z, z+2))) 
                    offTimeList.append(concat["Off_Time_ms"].rename("Off Time %i-%i"%(z, z+2)))
                    widthList.append(concat["Width_at50_ms"].rename("Width %i-%i"%(z, z+2)))
                    freqList.append(concat["Frequency"].rename("Frequency %i-%i"%(z, z+2)))
                    concat.to_excel(writer, sheet_name= "Traces %i-%i"%(z, z+2), index=False)
                    print(concat)
                    z += 3
                print(ampList)
                ampColumn = pd.concat(ampList, axis="columns")
                offTimeColumn = pd.concat(offTimeList, axis="columns")
                widthColumn = pd.concat(widthList, axis="columns")
                freqColumn = pd.concat(freqList, axis="columns")

                ampColumn.to_excel(writer, sheet_name="Amplitude")
                offTimeColumn.to_excel(writer, sheet_name="Off Time")
                widthColumn.to_excel(writer, sheet_name="Width")
                freqColumn.to_excel(writer, sheet_name="Frequency")
                # Individual Traces
                for frames in preInjectionLeft:
                    preInjectionLeft[frames].to_excel(writer, sheet_name= "Trace %i"%x, index= False)
                    x += 1
            with postLeftWriter as writer:
                # Overview Sheet
                postOverviewLeft.to_excel(writer, sheet_name= "All Traces")
                # Bins of Three
                postLeftGroupThree = [list(postInjectionLeft.values())[i:i+3] for i in range(0, len(postInjectionLeft), 3)]
                ampColumn = pd.DataFrame()
                offTimeColumn = pd.DataFrame()
                widthColumn = pd.DataFrame()
                freqColumn = pd.DataFrame()
                ampList = []
                offTimeList = []
                widthList = []
                freqList = []
                z = 1
                for groups in postLeftGroupThree:
                    concat = pd.DataFrame(columns= ['Event_Num', 'Peak_Index', 
                                        'Peak_Time_Sec', 'Event_Window_Start', 
                                        'Event_Window_End', 'Amplitude', 'Off_Time_ms',
                                        'Width_at50_ms','Frequency'])
                    for y in groups:
                        concat = pd.concat([concat, y if not y.empty else None], ignore_index=True)
                    ampList.append(pd.DataFrame(concat["Amplitude"].rename("Amplitude %i-%i"%(z+x, z+x+2))))
                    offTimeList.append(pd.DataFrame(concat["Off_Time_ms"].rename("Off Time %i-%i"%(z+x, z+x+2))))
                    widthList.append(pd.DataFrame(concat["Width_at50_ms"].rename("Width %i-%i"%(z+x, z+x+2))))
                    freqList.append(pd.DataFrame(concat["Frequency"].rename("Frequency %i-%i"%(z+x, z+x+2))))
                    concat.to_excel(writer, sheet_name= "Traces %i-%i"%(z+x, z+x+2), index=False)
                    z += 3

                ampColumn = pd.concat(ampList, axis="columns", join="inner")
                offTimeColumn = pd.concat(offTimeList, axis="columns", join="inner")
                widthColumn = pd.concat(widthList, axis="columns", join="inner")
                freqColumn = pd.concat(freqList, axis="columns", join="inner")

                ampColumn.to_excel(writer, sheet_name="Amplitude")
                offTimeColumn.to_excel(writer, sheet_name="Off Time")
                widthColumn.to_excel(writer, sheet_name="Width")
                freqColumn.to_excel(writer, sheet_name="Frequency")
                # Individual Traces
                for frames in postInjectionLeft:
                    postInjectionLeft[frames].to_excel(writer, sheet_name= "Trace %i"%x, index= False)
                    x += 1
            x = 1
            with preRightWriter as writer:
                # Overview Sheet
                preOverviewRight.to_excel(writer, sheet_name= "All Traces")
                # Bins of Three
                preRightGroupThree = [list(preInjectionRight.values())[i:i+3] for i in range(0, len(preInjectionRight), 3)]
                ampColumn = pd.DataFrame()
                offTimeColumn = pd.DataFrame()
                widthColumn = pd.DataFrame()
                freqColumn = pd.DataFrame()
                ampList = []
                offTimeList = []
                widthList = []
                freqList = []
                z = 1
                for groups in preRightGroupThree:
                    concat = pd.DataFrame(columns= ['Event_Num', 'Peak_Index', 
                                        'Peak_Time_Sec', 'Event_Window_Start', 
                                        'Event_Window_End', 'Amplitude', 'Off_Time_ms',
                                        'Width_at50_ms','Frequency'])
                    for y in groups:
                        concat = pd.concat([concat, y if not y.empty else None], ignore_index=True)
                    ampList.append(pd.DataFrame(concat["Amplitude"].rename("Amplitude %i-%i"%(z, z+2))))
                    offTimeList.append(pd.DataFrame(concat["Off_Time_ms"].rename("Off Time %i-%i"%(z, z+2))))
                    widthList.append(pd.DataFrame(concat["Width_at50_ms"].rename("Width %i-%i"%(z, z+2))))
                    freqList.append(pd.DataFrame(concat["Frequency"].rename("Frequency %i-%i"%(z, z+2))))
                    concat.to_excel(writer, sheet_name= "Traces %i-%i"%(z, z+2), index=False)
                    z += 3

                ampColumn = pd.concat(ampList, axis="columns", join="inner")
                offTimeColumn = pd.concat(offTimeList, axis="columns", join="inner")
                widthColumn = pd.concat(widthList, axis="columns", join="inner")
                freqColumn = pd.concat(freqList, axis="columns", join="inner")

                ampColumn.to_excel(writer, sheet_name="Amplitude")
                offTimeColumn.to_excel(writer, sheet_name="Off Time")
                widthColumn.to_excel(writer, sheet_name="Width")
                freqColumn.to_excel(writer, sheet_name="Frequency")
                
                # Individual Traces
                for frames in preInjectionRight:
                    preInjectionRight[frames].to_excel(writer, sheet_name= "Trace %i"%x, index= False)
                    x += 1
            with postRightWriter as writer:
                # Overview Sheet
                postOverviewRight.to_excel(writer, sheet_name= "All Traces")
                # Bins of Three
                postRightGroupThree = [list(postInjectionRight.values())[i:i+3] for i in range(0, len(postInjectionRight), 3)]
                ampColumn = pd.DataFrame()
                offTimeColumn = pd.DataFrame()
                widthColumn = pd.DataFrame()
                freqColumn = pd.DataFrame()
                ampList = []
                offTimeList = []
                widthList = []
                freqList = []
                z = 1
                for groups in postRightGroupThree:
                    concat = pd.DataFrame(columns= ['Event_Num', 'Peak_Index', 
                                        'Peak_Time_Sec', 'Event_Window_Start', 
                                        'Event_Window_End', 'Amplitude', 'Off_Time_ms',
                                        'Width_at50_ms','Frequency'])
                    for y in groups:
                        concat = pd.concat([concat, y if not y.empty else None], ignore_index=True)
                    ampList.append(pd.DataFrame(concat["Amplitude"].rename("Amplitude %i-%i"%(z+x, z+x+2))))
                    offTimeList.append(pd.DataFrame(concat["Off_Time_ms"].rename("Off Time %i-%i"%(z+x, z+x+2))))
                    widthList.append(pd.DataFrame(concat["Width_at50_ms"].rename("Width %i-%i"%(z+x, z+x+2))))
                    freqList.append(pd.DataFrame(concat["Frequency"].rename("Frequency %i-%i"%(z+x, z+x+2))))
                    concat.to_excel(writer, sheet_name= "Traces %i-%i"%(z+x, z+x+2), index=False)
                    z += 3
                
                ampColumn = pd.concat(ampList, axis="columns", join="inner")
                offTimeColumn = pd.concat(offTimeList, axis="columns", join="inner")
                widthColumn = pd.concat(widthList, axis="columns", join="inner")
                freqColumn = pd.concat(freqList, axis="columns", join="inner")

                ampColumn.to_excel(writer, sheet_name="Amplitude")
                offTimeColumn.to_excel(writer, sheet_name="Off Time")
                widthColumn.to_excel(writer, sheet_name="Width")
                freqColumn.to_excel(writer, sheet_name="Frequency")
                # Individual Traces
                for frames in postInjectionRight:
                    postInjectionRight[frames].to_excel(writer, sheet_name= "Trace %i"%x, index= False)
                    x += 1
            self.traceStatus = True
            acl.tExport(finalSignalLeft, self.ratNameLeft, self.abfDate) #Left
            acl.tExport(finalSignalRight, self.ratNameRight, self.abfDate) #Right

        # Analyzes the peak decay, amplitude, and frequency of a single trace chosen by the user.
        def singleTracePeaks():
            finalSignalLeft, finalSignalRight = acl.completeProcessor(self.experimentFileName, self.baselinefileName)
        # Left Rat Peak Analysis
            self.peaksLeft = pas.peakGetter(finalSignalLeft[self.trace])
            pas.peakDisplay(finalSignalLeft[self.trace], self.experimentFileName, "Left Rat")
        # Right Rat Peak Analysis
            self.peaksRight = sci.find_peaks(finalSignalRight[self.trace])
            pas.peakDisplay(finalSignalRight[self.trace], self.experimentFileName, "Right Rat")

        runFileButton = ttk.Button(self, text="Process a File", command= dataProcessorPop)
        runFileButton.grid(row= 3, column= 1)
        baselineGetterButton = ttk.Button(self, text="Get the Baselines", command= baselineFinder)
        baselineGetterButton.grid(row= 3, column= 2)
        testerButton = ttk.Button(self, text="Event analysis on a single trace", command= singleTracePeaks)
        testerButton.grid(row=4, column=1)

        

def main():
    fp = tk.Tk()
    fp.title("Liz's Data Corrector")
    window = Main(fp)
    window.pack()
    fp.mainloop()

if __name__ == "__main__":
    main()