from email import header
from lib2to3.pgen2.token import AMPER
from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.fftpack import fft, ifft, fftfreq, fftshift
# from numpy import fft 

matplotlib.rcParams.update({'font.size': 12}) #Default Font Size
matplotlib.rcParams['text.usetex'] = True
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def fourier_transform(n,u,dt):
    xf = fftfreq(n-1, dt)[:n//2]
    yf = fft(u,n-1)[:n//2] # for accuracte periodicity
    yf[1:] = 2*yf[1:]
    yf = (1/(n-1))*yf # for normalization
    return xf, yf


data = np.genfromtxt('probe_output.dat',skip_header=3)

dt = 0.005
t = data[:,0]
v = data[:,2] 
u = data[:,1]

freq_2, amp = fourier_transform(len(t[:]),v,dt)

# amp = fft.fft(v)
# freq_2 = np.fft.fftfreq(len(v))
# threshold = 0.5*max(abs(amp))
# mask = abs(amp)> threshold
# peaks = freq_2[mask]
# print(peaks/0.005)

Re = 150
print(freq_2[np.argmax(amp)])
plt.figure(1)
plt.plot(freq_2,np.abs(amp),color='dodgerblue',linewidth=2)
plt.text(freq_2[np.argmax(amp)]+0.1,np.abs(amp)[np.argmax(amp)],'$St = {:.4f}$'.format(freq_2[np.argmax(amp)]))
# plt.plot([freq_2[np.argmax(amp)],freq_2[np.argmax(amp)]],[1e-5,np.abs(amp)[np.argmax(amp)]], '--', color='lightcoral',linewidth=0.5)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$fD/U_\infty$')
plt.ylabel(r'|A|')
plt.title(r'Fourier Transform of $v$ at $Re$ = '+str(Re))
plt.show()
# freq, amp = fourier_transform(len(t),u)