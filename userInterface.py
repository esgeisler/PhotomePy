import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

def fileBrowser():
    fileName = filedialog.askopenfilename(initialdir= "/", title= "Select a Main File", filetypes=(("Axon Binary Fles", "*.abf*"), ("All Files," "*.*")))


fp = tk.Tk()
fp.title("Liz's Data Corrector")

decisionLabel = tk.Label(fp, text="File Process")
userChosenProcess = ttk.Combobox(fp, state= "readonly", values=["Get Baseline", "Get Averages"])
userChosenProcess.set("Get Baseline")

explorerButton = tk.Button(fp, text="Choose a Main File", command=fileBrowser)
baselineExplorerButton = tk.Button(fp, text="Choose a Baseline File", command=fileBrowser)
runFileButton = tk.Button(fp, text="Run Program")

decisionLabel.grid(column= 1, row= 1)
userChosenProcess.grid(column= 2, row= 1)
explorerButton.grid(column= 1, row= 2)
baselineExplorerButton.grid(column= 2, row= 2)
runFileButton.grid(column=2, row= 3)

fp.mainloop()