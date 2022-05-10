import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import RegularGridInterpolator


def get_data(fname,col,nx,ny):
    q_data = np.genfromtxt(fname, skip_header=3)

    x1 = q_data[:,0]
    y1 = q_data[:,1]
    q1 = q_data[:,col-1]

    xq = np.reshape(x1,(nx,ny),order='F') # order F for fortran is import for proper reshape
    yq = np.reshape(y1,(nx,ny),order='F')
    data  = np.reshape(q1,(nx,ny),order='F')
    return xq,yq,data

tpick = 90
nx_init = 200
ny_init = 128
nx_init = nx_init + 2
ny_init = ny_init + 2 # add ghost cells

nx_final = 200
ny_final = 128
nx_final = nx_final + 2
ny_final = ny_final + 2 # add ghost cells

fname_f = 'Data/data.{:>07d}.dat'.format(tpick)

x_init, y_init, u_init = get_data(fname_f,3,nx_init,ny_init)
x_init, y_init, v_init = get_data(fname_f,4,nx_init,ny_init)
x_init, y_init, p_init = get_data(fname_f,5,nx_init,ny_init)

print(np.max(x_init), np.max(y_init))
x_init_1D = x_init[:,0]
y_init_1D = y_init[0,:]

interpolate_u = RegularGridInterpolator((x_init_1D, y_init_1D), u_init)
interpolate_v = RegularGridInterpolator((x_init_1D, y_init_1D), v_init)
interpolate_p = RegularGridInterpolator((x_init_1D, y_init_1D), p_init)

# Reading face data
xf_final_1D = np.genfromtxt('xgrid.dat')
xf_final_1D = xf_final_1D[:,1]
yf_final_1D = np.genfromtxt('ygrid.dat')
yf_final_1D = yf_final_1D[:,1]

    # ! setting up the cell centers
    # do i = 2,nx-1
    #     x(i) = (xf(i-1)+xf(i))/2.0 
    # end do
    
    # do j = 2,ny-1
    #     y(j) = (yf(j-1)+yf(j))/2.0
    # end do


    # x(1) = -x(2)
    # y(1) = -y(2)
    # x(nx) = xf(nxf) + (xf(nxf) - x(nx-1)) 
    # y(ny) = yf(nyf) + (yf(nyf) - y(ny-1))

# computing the cell centers
x_final_1D = np.zeros(nx_final)
y_final_1D = np.zeros(ny_final)

for i in range(1,nx_final-1):
    x_final_1D[i] = (xf_final_1D[i-1] + xf_final_1D[i])/2.0

for j in range(1,ny_final-1):
    y_final_1D[j] = (yf_final_1D[j-1] + yf_final_1D[j])/2.0

x_final_1D[0] = -x_final_1D[1]
y_final_1D[0] = -y_final_1D[1]
x_final_1D[nx_final-1] = x_final_1D[nx_final-2] + (xf_final_1D[-1] - xf_final_1D[-2])
y_final_1D[ny_final-1] = y_final_1D[ny_final-2] + (yf_final_1D[-1] - yf_final_1D[-2])


x_final, y_final = np.meshgrid(x_final_1D, y_final_1D, indexing='ij')

print(np.max(x_final), np.max(y_final))

ptsq = np.zeros((nx_final,ny_final,2))
for j in range(0,ny_final):
    for i in range(0,nx_final):
        ptsq[i,j,0] = x_final[i,j]
        ptsq[i,j,1] = y_final[i,j]
locq = ptsq.reshape(nx_final*ny_final,2)

print(locq.shape, nx_final*ny_final)
u_final = interpolate_u(locq)
v_final = interpolate_v(locq)
p_final = interpolate_p(locq)

u_final = u_final.reshape(nx_final,ny_final)
v_final = v_final.reshape(nx_final,ny_final)
p_final = p_final.reshape(nx_final,ny_final)


# plt.contourf(x_final,y_final,u_final)
# plt.show()

f = open('data_init_refined.dat','w')
f.write("Title = Data Refined" + "\n")
f.write("Variables = X,Y,U,V,P" + "\n")
f.write("ZONE I = " + str(nx_final) + " J = " + str(ny_final) + " F=POINT" + "\n")
for j in range(0,ny_final):
    for i in range(0,nx_final):
        f.write(str(x_final[i,j]) + "  " + str(y_final[i,j]) + "  " + str(u_final[i,j]) + "  " + str(v_final[i,j]) + "  " + str(p_final[i,j]) + "\n")

f.close()