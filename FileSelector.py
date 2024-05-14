import os
import Exceptions

#TODO: prevent FileName from being incorrect during while loop
def FileSelect(ProOrRaw, FileName):
    x = 0
    while x == 0:
        if ProOrRaw.lower() == "raw":
            ChosenFile = os.path.join(os.getcwd(), "Raw Data", FileName)
            x = 1
        elif ProOrRaw.lower() == "pro" or "processed":
            ChosenFile = os.path.join(os.getcwd(), "Processed Data", FileName)
            x = 1
        else:
            print("Please type either 'raw' or 'processed'.")
            ProOrRaw = input("Is the data raw or processed?")
    return ChosenFile