import AutoCleaner as acl
import pyabf
import matplotlib.pyplot as plt
import AverageTraces as avg

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

def fileBrowser():
    fileName = filedialog.askopenfilename(initialdir= "/", title= "Select a Main File", filetypes=(("Axon Binary Fles", "*.abf*"), ("All Files," "*.*")))

class Main(tk.Frame):
    def __init__(self, master= None, **kwargs):
        super().__init__(master, **kwargs)

        decisionLabel = tk.Label(self, text="File Process")
        userChosenProcess = ttk.Combobox(self, state= "readonly", values=["Get Baseline", "Get Averages"])
        userChosenProcess.set("Get Baseline")
        decisionLabel.pack()
        userChosenProcess.pack()     

        explorerButton = tk.Button(self, text="Choose a Main File", command= fileBrowser)
        baselineExplorerButton = tk.Button(self, text="Choose a Baseline File", command= fileBrowser)
        explorerButton.pack()
        baselineExplorerButton.pack()
    
    

# runFileButton = tk.Button(fp, text="Run Program")
# explorerButton.grid(column= 1, row= 2)
# baselineExplorerButton.grid(column= 2, row= 2)
# runFileButton.grid(column=2, row= 3)


# abf = pyabf.ABF(userFile)
# if dataType == "raw":
#     baselineFile = fls.FileSelect("raw", input("What is the name of the file to be used for baseline?"))
# decision =  input("What shall we do today?")

# if decision == "baseline":
#     baselineSubL = acl.LBaselineGet(baselineFile)
#     baselineSubR = acl.RBaselineGet(baselineFile)
#     outputString = ("Left - 470: %.2f 405: %.2f\nRight - 470: %.2f 405: %.2f"%(baselineSubL[0], baselineSubL[1], baselineSubR[0], baselineSubR[1]))
#     print(outputString)
# elif decision == "process":
    #userTrace = int(input("Which trace would you like to see?"))
    #userChannel = int(input("Which channel?"))

def main():
    fp = tk.Tk()
    fp.title("Liz's Data Corrector")
    window = Main(fp)
    window.pack()
    fp.mainloop()

if __name__ == "__main__":
    main()