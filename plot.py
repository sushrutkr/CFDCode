from email import header
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 

matplotlib.rcParams.update({'font.size': 12}) #Default Font Size
matplotlib.rcParams['text.usetex'] = True

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

Re = 150
tcut = 1000
force = np.genfromtxt('force_history.dat',skip_header=4)

plt.figure(1)
plt.plot(force[tcut:,0],force[tcut:,2],'r-',label='$C_D$',color='dodgerblue')
plt.xlabel(r'$tU_\infty/D$')
plt.ylabel('$C_D$')
plt.title(r'$C_D$ vs $tU_\infty/D$ at $Re$ = '+str(Re))
plt.show()

plt.figure(2)
plt.plot(force[tcut:,0],force[tcut:,4],'b-',label='$C_L$',color='dodgerblue')
plt.xlabel(r'$tU_\infty/D$')
plt.ylabel('$C_L$')
plt.title(r'$C_L$ vs $tU_\infty/D$ at $Re$ = '+str(Re))
plt.show()

probe = np.genfromtxt('probe_output.dat',skip_header=3)
tcut = 1000
plt.figure(3)
plt.plot(probe[tcut:,0],probe[tcut:,1],'dodgerblue',label=r'$u/U_\infty$')
plt.plot(probe[tcut:,0],probe[tcut:,2],'salmon',label=r'$v/U_\infty$')
plt.xlabel(r'$tU_\infty/D$')
plt.ylabel('$u/U_\infty ,v/U_\infty$')
plt.title(r'Probe Output at $Re$ = '+str(Re))
plt.legend()
plt.show()
