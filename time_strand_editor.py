# ----------------------------------------------------
# Flow Physics and Computation Lab
# PI : Prof. Rajat Mittal
#
# Author : Sushrut Kumar
# Date : 8th July 2021
#
# -----------------------------------------------------
# Writes macro file for setting time strands in Tecplot
# First load .plt files in tecplot and then load marker files. After this play the script generated from the code.



import numpy as np
import matplotlib.pyplot as plt

nt = int(20)
ft = 20
time = np.linspace(1,ft,nt) #Initialtime, finaltime, ntsteps
motion = 1 # 0 - stationanry body, 1 - moving body
zone = 1 # 1 - only marker, 2 - marker and flow data
no_body = 1 # If marker files has more than 1 body
no_datasets = len(time)
print(no_datasets)

f = open("edit_time_strand.mcr", "w")
f.write("#!MC 1410 \n")


if zone==1:
    b1 = 1
    b2 = 2
    for ig in range(0,no_datasets):
        i = ig+1
        f.write("$!ExtendedCommand \n")
        f.write("  CommandProcessorID = 'Strand Editor' \n")
        if motion==0:
            f.write("  Command = 'ZoneSet=" + str(i) + ";AssignStrands=TRUE;StrandValue=" + str(time[ig]) + ";AssignSolutionTime=TRUE;TimeValue=" + str(time[ig]) + ";TimeOption=SingleValue;' \n")

        else:    
            if no_body==1:
                f.write("  Command = 'ZoneSet=" + str(i) +  ";AssignStrands=TRUE;StrandValue=" + str(time[ig]) + ";AssignSolutionTime=TRUE;TimeValue=" + str(time[ig]) + ";TimeOption=SingleValue;' \n")
            else:
                f.write("  Command = 'ZoneSet=" + str(b1) + "-" + str(b2) + ";AssignStrands=TRUE;StrandValue=" + str(time[ig]) + ";AssignSolutionTime=TRUE;TimeValue=" + str(time[ig]) + ";TimeOption=SingleValue;' \n")
                b1 = b2+1
                b2 = b1+no_body-1

elif zone==2:
    b1 = 1+no_datasets
    b2 = 1+no_datasets+no_body-1
    for ig in range(0,no_datasets):
        i = ig+1
        f.write("$!ExtendedCommand \n")
        f.write("  CommandProcessorID = 'Strand Editor' \n")
        if motion==0:
            f.write("  Command = 'ZoneSet=" + str(i) + ";AssignStrands=TRUE;StrandValue=" + str(time[ig]) + ";AssignSolutionTime=TRUE;TimeValue=" + str(time[ig]) + ";TimeOption=SingleValue;' \n")

        else:    
            if no_body==1:
                f.write("  Command = 'ZoneSet=" + str(i) + "," + str(i+no_datasets) + ";AssignStrands=TRUE;StrandValue=" + str(time[ig]) + ";AssignSolutionTime=TRUE;TimeValue=" + str(time[ig]) + ";TimeOption=SingleValue;' \n")
            else:
                f.write("  Command = 'ZoneSet=" + str(i) + "," + str(b1) + "-" + str(b2) + ";AssignStrands=TRUE;StrandValue=" + str(time[ig]) + ";AssignSolutionTime=TRUE;TimeValue=" + str(time[ig]) + ";TimeOption=SingleValue;' \n")
                b1 = b2+1
                b2 = b2+no_body
f.close()
