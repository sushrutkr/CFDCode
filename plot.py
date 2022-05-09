from email import header
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 

matplotlib.rcParams.update({'font.size': 12}) #Default Font Size
matplotlib.rcParams['text.usetex'] = True

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

tcut = 1000
force = np.genfromtxt('force_history.dat',skip_header=4)

plt.figure(1)
plt.plot(force[tcut:,0],force[tcut:,2],'r-',label='$C_D$')
plt.plot(force[tcut:,0],force[tcut:,4],'b-',label='$C_L$')
plt.xlabel('$t$')
plt.ylabel('$C_D, C_L$')
plt.legend()
plt.show()

# probe = np.genfromtxt('probe_output.dat',skip_header=3)
# tcut = 1000
# plt.figure(1)
# plt.plot(probe[tcut:,0],probe[tcut:,1],'r-',label='$u$')
# plt.plot(probe[tcut:,0],probe[tcut:,2],'b-',label='$v$')
# plt.plot(probe[tcut:,0],probe[tcut:,3],'g-',label='$p$')
# plt.xlabel('$t$')
# plt.ylabel('$u,v,p$')
# plt.legend()
# plt.show()
