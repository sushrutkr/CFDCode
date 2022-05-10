import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.ndimage
from matplotlib.gridspec import GridSpec
from scipy.fftpack import fft, ifft, fftfreq, fftshift

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

def get_data(fname,col,nx,ny):
    q_data = np.genfromtxt(fname, skip_header=3)

    x1 = q_data[:,0]
    y1 = q_data[:,1]
    q1 = q_data[:,col-1]

    xq = np.reshape(x1,(nx,ny),order='F') # order F for fortran is import for proper reshape
    yq = np.reshape(y1,(nx,ny),order='F')
    data  = np.reshape(q1,(nx,ny),order='F')
    return xq,yq,data

nx = 246
ny = 160
nx = nx + 2
ny = ny + 2

lx = 10
ly = 8
Re = 150
x_cen = 3
y_cen = 4
dt = 0.005
time = np.array([[5,15],[25,50],[100,150]])

fname_u = 'Data/'+str(Re)+'_'+str(lx)+'_'+str(ly)+'_'+'u'+'.png'
fname_v = 'Data/'+str(Re)+'_'+str(lx)+'_'+str(ly)+'_'+'v'+'.png'
fname_p = 'Data/'+str(Re)+'_'+str(lx)+'_'+str(ly)+'_'+'p'+'.png'
fname_vor = 'Data/'+str(Re)+'_'+str(lx)+'_'+str(ly)+'_'+'vor'+'.png'
fname_st = 'Data/'+str(Re)+'_'+str(lx)+'_'+str(ly)+'_'+'st'+'.png'
fname_probe = 'Data/'+str(Re)+'_'+str(lx)+'_'+str(ly)+'_'+'probe'+'.png'
fname_cl = 'Data/'+str(Re)+'_'+str(lx)+'_'+str(ly)+'_'+'cl'+'.png'
fname_cd = 'Data/'+str(Re)+'_'+str(lx)+'_'+str(ly)+'_'+'cd'+'.png'

#============= Strouhal Number
data = np.genfromtxt('probe_output.dat',skip_header=3)

dt = 0.005
t_probe = data[:,0]
v_probe = data[:,2] 
u_probe = data[:,1]

freq_2, amp = fourier_transform(len(t_probe[:]),v_probe,dt)
plt.figure(1)
plt.plot(freq_2,np.abs(amp),color='dodgerblue',linewidth=2)
plt.text(freq_2[np.argmax(amp)]+0.1,np.abs(amp)[np.argmax(amp)],'$St = {:.4f}$'.format(freq_2[np.argmax(amp)]))
# plt.plot([freq_2[np.argmax(amp)],freq_2[np.argmax(amp)]],[1e-5,np.abs(amp)[np.argmax(amp)]], '--', color='lightcoral',linewidth=0.5)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$fD/U_\infty$')
plt.ylabel(r'A')
plt.title(r'Fourier Transform of $v?U_\infty$ at $Re$ = '+str(Re))
plt.savefig(fname_st,dpi=600,bbox_inches='tight')
plt.close()

#============= Drag Coefficient
tcut = 1000
force = np.genfromtxt('force_history.dat',skip_header=4)

mean_cd = np.mean(force[:-5000,2])
mean_cl = np.mean(force[:-5000,4])
print('mean CD = ',mean_cd)
print('mean CL = ',mean_cl)

plt.figure(1)
plt.plot(force[tcut:,0],force[tcut:,2],label='$C_D$',color='dodgerblue')
plt.text(min(force[-15000:,0]),1.15,r'Mean $C_D = {:.4f}$'.format(mean_cd))
plt.xlabel(r'$tU_\infty/D$')
plt.ylabel('$C_D$')
plt.title(r'$C_D$ vs $tU_\infty/D$ at $Re$ = '+str(Re))
plt.savefig(fname_cd,dpi=600,bbox_inches='tight')
plt.close()

plt.figure(2)
plt.plot(force[tcut:,0],force[tcut:,4],label='$C_L$',color='dodgerblue')
# plt.text(min(force[-5000:,0]),1.15,'$C_L = {:.4f}$'.format(mean_cl))
plt.xlabel(r'$tU_\infty/D$')
plt.ylabel('$C_L$')
plt.title(r'$C_L$ vs $tU_\infty/D$ at $Re$ = '+str(Re))
plt.savefig(fname_cl,dpi=600,bbox_inches='tight')
plt.close()

probe = np.genfromtxt('probe_output.dat',skip_header=3)
tcut = 1000
plt.figure(3)
plt.plot(probe[tcut:,0],probe[tcut:,1],'dodgerblue',label=r'$u/U_\infty$')
plt.plot(probe[tcut:,0],probe[tcut:,2],'salmon',label=r'$v/U_\infty$')
plt.xlabel(r'$tU_\infty/D$')
plt.ylabel('$u/U_\infty ,v/U_\infty$')
plt.title(r'Probe Output at $Re$ = '+str(Re))
plt.legend()
plt.savefig(fname_probe,dpi=600,bbox_inches='tight')
plt.close()
#============== Contours


theta = np.linspace( 0 , 2 * np.pi , 150 )
radius = 0.5
a = radius * np.cos( theta )
b = radius * np.sin( theta )
a = a + x_cen
b = b + y_cen

nrow = len(time[:,0])
ncol = len(time[0,:])
fig = plt.figure(constrained_layout=True, figsize=(10,10))
widths = [5,5,0.1]
gs = GridSpec(nrow, ncol + 1, figure=fig, width_ratios=widths)
axes = []
for i in range(nrow):
    for j in range(ncol):
        fname_read = 'Data/data.{:>07d}.dat'.format(time[i,j])
        x , y , u  = get_data(fname_read,3,nx,ny)
        axes.append(fig.add_subplot(gs[i, j]))
        im = axes[-1].contourf(x,y,u,cmap='Blues')
        axes[-1].plot(a,b,'k-',linewidth=0.5)
        axes[-1].set_title(r'$tU_\infty/D = {:.0f}$'.format(time[i,j]))
        axes[-1].set_aspect('equal', adjustable='box')
        axes[-1].set_xlim(0,lx)
        axes[-1].set_ylim(1,ly-1)
        if j != 0:
            axes[-1].set_yticklabels([])
        if i != nrow - 1:
            axes[-1].set_xticklabels([])
        if j == 0:
            axes[-1].set_ylabel(r'$y/D$')
        if i == nrow - 1:
            axes[-1].set_xlabel(r'$x/D$')
  
axes.append(fig.add_subplot(gs[:, ncol]))
fig.colorbar(im, cax=axes[-1], orientation='vertical')
fig.suptitle(r'$u/U_\infty$ for $Re$ = ' + str(Re))
plt.savefig(fname_u, dpi=600)
plt.close()

fig = plt.figure(constrained_layout=True, figsize=(10,10))
widths = [5,5,0.1]
gs = GridSpec(nrow, ncol + 1, figure=fig, width_ratios=widths)
axes = []
for i in range(nrow):
    for j in range(ncol):
        fname_read = 'Data/data.{:>07d}.dat'.format(time[i,j])
        x , y , v  = get_data(fname_read,4,nx,ny)
        axes.append(fig.add_subplot(gs[i, j]))
        im = axes[-1].contourf(x,y,v,cmap='Blues')
        axes[-1].plot(a,b,'k-',linewidth=0.5)
        axes[-1].set_title(r'$tU_\infty/D = {:.0f}$'.format(time[i,j]))
        axes[-1].set_aspect('equal', adjustable='box')
        axes[-1].set_xlim(0,lx)
        axes[-1].set_ylim(1,ly-1)
        if j != 0:
            axes[-1].set_yticklabels([])
        if i != nrow - 1:
            axes[-1].set_xticklabels([])
        if j == 0:
            axes[-1].set_ylabel(r'$y/D$')
        if i == nrow - 1:
            axes[-1].set_xlabel(r'$x/D$')
  
axes.append(fig.add_subplot(gs[:, ncol]))
fig.colorbar(im, cax=axes[-1], orientation='vertical')
fig.suptitle(r'$v/U_\infty$ for $Re$ = ' + str(Re))
plt.savefig(fname_v, dpi=600)
plt.close()

fig = plt.figure(constrained_layout=True, figsize=(10,10))
widths = [5,5,0.1]
gs = GridSpec(nrow, ncol + 1, figure=fig, width_ratios=widths)
axes = []
for i in range(nrow):
    for j in range(ncol):
        fname_read = 'Data/data.{:>07d}.dat'.format(time[i,j])
        x , y , p  = get_data(fname_read,5,nx,ny)
        p = p - p[2,2]
        axes.append(fig.add_subplot(gs[i, j]))
        im = axes[-1].contourf(x,y,p,cmap='Blues')
        axes[-1].plot(a,b,'k-',linewidth=0.5)
        axes[-1].set_title(r'$tU_\infty/D = {:.0f}$'.format(time[i,j]))
        axes[-1].set_aspect('equal', adjustable='box')
        axes[-1].set_xlim(0,lx)
        axes[-1].set_ylim(1,ly-1)
        if j != 0:
            axes[-1].set_yticklabels([])
        if i != nrow - 1:
            axes[-1].set_xticklabels([])
        if j == 0:
            axes[-1].set_ylabel(r'$y/D$')
        if i == nrow - 1:
            axes[-1].set_xlabel(r'$x/D$')
  
axes.append(fig.add_subplot(gs[:, ncol]))
fig.colorbar(im, cax=axes[-1], orientation='vertical')
fig.suptitle(r'$(p - p_\infty)/\rho U_\infty^2$ for $Re$ = ' + str(Re))
plt.savefig(fname_p, dpi=600)
plt.close()

fig = plt.figure(constrained_layout=True, figsize=(10,10))
widths = [5,5,0.1]
gs = GridSpec(nrow, ncol + 1, figure=fig, width_ratios=widths)
axes = []
for i in range(nrow):
    for j in range(ncol):
        fname_read = 'Data/data.{:>07d}.dat'.format(time[i,j])
        x , y , vor  = get_data(fname_read,6,nx,ny)
        axes.append(fig.add_subplot(gs[i, j]))
        vorticity_levels = [np.min(vor),np.min(vor)/2,-10,-5,-4.5,-4,-3.5,-3,-2.5,-2,-1,
                            1,2,2.5,3,3.5,4,4.5,5,10,np.max(vor)/2,np.max(vor)]
        im = axes[-1].contourf(x,y,vor,cmap='seismic', levels=vorticity_levels)
        axes[-1].plot(a,b,'k-',linewidth=0.5)
        axes[-1].set_title(r'$tU_\infty/D = {:.0f}$'.format(time[i,j]))
        axes[-1].set_aspect('equal', adjustable='box')
        axes[-1].set_xlim(0,lx)
        axes[-1].set_ylim(1,ly-1)
        if j != 0:
            axes[-1].set_yticklabels([])
        if i != nrow - 1:
            axes[-1].set_xticklabels([])
        if j == 0:
            axes[-1].set_ylabel(r'$y/D$')
        if i == nrow - 1:
            axes[-1].set_xlabel(r'$x/D$')
  
axes.append(fig.add_subplot(gs[:, ncol]))
fig.colorbar(im, cax=axes[-1], orientation='vertical')
fig.suptitle(r'$\omega D/U_\infty$ for $Re$ = ' + str(Re))
plt.savefig(fname_vor, dpi=600)
plt.close()

# fig, ax = plt.subplots(3,2,sharex=True,sharey=True,constrained_layout=True)
# fig.suptitle(r'$v/U_\infty$ for $Re$ = ' + str(Re))
# for j in range(0,len(time[0,:])):
#     for i in range(0,len(time[:,0])):
#         fname_read = 'Data/data.{:>07d}.dat'.format(time[i,j])
#         x , y , u  = get_data(fname_read,3,nx,ny)
#         # fname_save = '_{:>07d}.png'.format(i)
#         # fname_u = fname_u + fname_save
#         # fname_v = fname_v + fname_save
#         # fname_p = fname_p + fname_save
#         # fname_vor = fname_vor + fname_save
#         ax[i,j].set_title(r'$tU_\infty/D = {:.0f}$'.format(time[i,j]))
#         ab = ax[i,j].contourf(x,y,u,cmap='Blues')
#         ax[i,j].set_aspect('equal', adjustable='box')
#         ax[i,j].set_xlim(0,lx)
#         ax[i,j].set_ylim(1,ly-1)

# for ax in ax.flat:
#     ax.set(xlabel=r'$x/D$', ylabel=r'$y/D$')
#     print(ax)

# plt.show()  

# a =np.linspace(0,3,10)
# c = np.flip(a)
# a = np.concatenate((a,c))
# b = np.tanh(a)

# plt.scatter(a,1*np.ones(len(a)))
# plt.scatter(b,2*np.ones(len(b)))
# plt.show()
# for i in time:
#     fname_read = 'Data/data.{:>07d}.dat'.format(i)
#     fname_save = '_{:>07d}.png'.format(i)
#     fname_u = fname_u + fname_save
#     fname_v = fname_v + fname_save
#     fname_p = fname_p + fname_save
#     fname_vor = fname_vor + fname_save

#     x , y , u  = get_data(fname_read,3,nx,ny)
#     x , y , v  = get_data(fname_read,4,nx,ny)
#     x , y , p  = get_data(fname_read,5,nx,ny)
#     p = p - p[2,2]
#     x , y , vor = get_data(fname_read,6,nx,ny)
#     x, y, ib = get_data(fname_read,7,nx,ny)

#     u = u*ib
#     v = v*ib
#     vor = vor*ib
#     p = p*ib

#     # y = scipy.ndimage.zoom(y, 3)
#     # x = scipy.ndimage.zoom(x, 3)
#     # vor = scipy.ndimage.zoom(vor, 3)
#     # p = scipy.ndimage.zoom(p, 3)
#     # u = scipy.ndimage.zoom(u, 3)
#     # v = scipy.ndimage.zoom(v, 3)
    
#     plt.figure(1,figsize=(lx+2,ly-2))
#     plt.contourf(x,y,u,cmap='Blues')
#     plt.plot(a,b,'k',linewidth=0.5)
#     plt.title(r'$u/U_\infty$ at $tU_\infty/D = {:.0f}$ for $Re$ = '.format(i) + str(Re))
#     plt.colorbar()
#     plt.xlabel(r'$x/D$')
#     plt.ylabel(r'$y/D$')
#     plt.xlim(0,lx)
#     plt.ylim(1,ly-1)
#     plt.gca().set_aspect('equal', adjustable='box')
#     plt.savefig(fname_u,dpi=600,bbox_inches='tight')

#     plt.figure(2,figsize=(lx+2,ly-2))
#     plt.contourf(x,y,v,cmap='Blues')
#     plt.colorbar()
#     plt.plot(a,b,'k',linewidth=0.5)
#     plt.title(r'$v/U_\infty$ at $tU_\infty/D = {:.0f}$ for $Re$ = '.format(i) + str(Re))
#     plt.xlabel(r'$x/D$')
#     plt.ylabel(r'$y/D$')
#     plt.xlim(0,lx)
#     plt.ylim(1,ly-1)
#     plt.gca().set_aspect('equal', adjustable='box')
#     plt.savefig(fname_v,dpi=600,bbox_inches='tight')

#     plt.figure(3,figsize=(lx+2,ly-2))
#     plt.contourf(x,y,p,cmap='Blues')
#     plt.plot(a,b,'k',linewidth=0.5)
#     plt.title(r'$(p - p_\infty)/\rho U_\infty^2$ at $tU_\infty/D = {:.0f}$ for $Re$ = '.format(i) + str(Re))
#     plt.colorbar()
#     plt.xlabel(r'$x/D$')
#     plt.ylabel(r'$y/D$')
#     plt.xlim(0,lx)
#     plt.ylim(1,ly-1)
#     plt.gca().set_aspect('equal', adjustable='box')
#     plt.savefig(fname_p,dpi=600,bbox_inches='tight')


#     vorticity_levels = [np.min(vor),np.min(vor)/2,-10,-5,-4.5,-4,-3.5,-3,-2.5,-2,
#                         2,2.5,3,3.5,4,4.5,5,10,np.max(vor)/2,np.max(vor)]
#     # vorticity_levels = np.tanh(np.linspace(min(vor),max(vor),31))
#     plt.figure(4,figsize=(lx+2,ly-2))
#     plt.contourf(x,y,vor,cmap='seismic', levels=vorticity_levels)
#     plt.plot(a,b,'k',linewidth=0.5)
#     plt.title(r'$\omega D/U_\infty$ at $tU_\infty/D = {:.0f}$ for $Re$ = '.format(i) + str(Re))
#     # plt.imshow(vor,extent=[x.min(),x.max(),y.min(),y.max()])
#     plt.colorbar()
#     plt.xlabel(r'$x/D$')
#     plt.ylabel(r'$y/D$')
#     plt.xlim(0,lx)
#     plt.ylim(1,ly-1)
#     plt.gca().set_aspect('equal', adjustable='box')
#     plt.savefig(fname_vor,dpi=600,box_inches='tight')
#     ticks = np.linspace(vor.min(),vor.max(),10)
#     plt.close('all')
