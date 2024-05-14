import os
import Exceptions

def FileSelect(ProOrRaw, FileName):
    if ProOrRaw.lower() == "raw":
        ChosenFile = os.getcwd() + "/Raw Data/" + FileName
        return ChosenFile
    elif ProOrRaw.lower() == "pro" or "processed":
        ChosenFile = os.getcwd() + "/Processed Data/" + FileName
        return ChosenFile
    else:
        return FileNotFoundError