from email import header
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 

matplotlib.rcParams.update({'font.size': 12}) #Default Font Size
matplotlib.rcParams['text.usetex'] = True

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

drag = np.genfromtxt('drag_history.dat',skip_header=4)

plt.figure(1)
plt.plot(drag[2500:,0],drag[2500:,2],'r-',label='$C_D$')
plt.xlabel('$t$')
plt.ylabel('$C_D$')
plt.legend()
plt.show()