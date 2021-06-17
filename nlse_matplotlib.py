# coding:utf8

import numpy as np
#from scipy import integrate
import matplotlib.pyplot as plt
from Functions import split_step



#------------------- Grid: --------------------#
#Input pulse
#Tmax = 5 ~ 10ps
T0 = 5E-12 # for pulse Initial width T0 --> Dispersive effects at T0 ~ 1ps  S.64

#T0->10E-12 Beta2 -> 2 P0-> 10E-3 Gamma ->2


# half-width (at 1/e-intensity point)
T_FWHM_g = 2*np.sqrt(np.log(2))*T0
#For Sech pulses:  
C = 0 #Chirp parameter
T_FWHM_sech = 2*np.log(1+np.sqrt(2))*T0
N = 8196 #ammount of points 
T_selected = T0
dt = 100*T_selected/N #25*T_selected/N #the 16 is to get a grid between -8 and 8 for T/T_selected 
T = np.arange(-N/2, N/2)*dt #np.arange(-N/2, N/2) mit len = N

beta2 = 25*1E-24 #ps^2/km  #S.65

#The use of these values in Eq. (3.1.5) shows that the dispersive and nonlinear effects are negligible for L < 50 km if T_selected > 100 ps and P0 ~ 1 mW
#gamma = n2*wo/(speed*Aeff)
gamma = 1 #1/(W*km)
P0 = 10000E-3 #W
#P0 = 1000E-3 #W
print('P0 = {0} W, \u03B3 = {1} 1/(W*km)  , T0 = {2} s , \u03B22 = {3} ps^2/km'.format(P0, gamma, T0, beta2))
#LNL= 1/(gamma*P0)
#Leff = LNL
zmax =2#LD#0.1 # L in km


if beta2 == 0:
    beta3 = 10*1E-36#ps^3/km
    print('N´ = sqrt(L´D/LNL) = : ',np.sqrt((gamma*P0*T_selected**3/np.absolute(beta3))))
    U, UI, UW, UIW, W, z = split_step(beta2, T_selected, T, zmax, gamma, P0, beta3,pulse = 'Gaussian', m = 1, C=0)
else:
    LD = (T_selected**2)/np.absolute(beta2)
    U, UI, UW, UIW, W, z = split_step(beta2, T_selected, T, zmax, gamma, P0 ,pulse = 'Gaussian', m = 1, C=0)

plt.figure(1)
plt.title('NLSE time')
plt.plot(T/T_selected, UI[0], label = 'step = 0')
plt.plot(T/T_selected, UI[49], label = 'step at {0} km'.format( '%.3f' %(z[49])))
plt.plot(T/T_selected, UI[-1], label = 'last step at {0} km'.format('%.3f' %(z[-1])))
plt.xlabel('T/T0')
plt.xlim((-10, 10))
plt.ylabel('Intensity UI')
plt.legend()
plt.grid()


# plt.figure(2)
# plt.title('Phase time')
# plt.plot(T/T_selected, np.angle(U[0]), label = 'step = 0')
# plt.plot(T/T_selected, np.angle(U[49]), label = 'step at {0} km'.format( '%.3f' %(z[49])))
# plt.plot(T/T_selected, np.angle(U[-1]), label = 'last step at {0} km'.format('%.3f' %(z[-1])))
# plt.xlabel('T/T0')
# #plt.xlim((-5, 5))
# plt.ylabel('Phase \u03C6NL')
# plt.legend()
# plt.grid()

plt.figure(3)
plt.title('NLSE freq')
plt.plot(W, UIW[0], label = 'step = 0')
plt.plot(W, UIW[49], label = 'step at {0} km'.format( '%.3f' %(z[49])))
plt.plot(W, UIW[-1], label = 'last step at {0} km'.format('%.3f' %(z[-1])))
plt.xlabel('\u03C9')
plt.xlim((-2E12, 2E12))
plt.ylabel('Intensity UI(\u03C9)')
plt.grid()
plt.legend()

yv, zv = np.meshgrid(T/T_selected, z)
plt.figure(4)
ax = plt.axes(projection='3d')
ax.plot_surface(yv, zv, UI,cmap='jet')
ax.set_title('NLSE Colormap time')
ax.set_xlabel('T/T0')
ax.set_ylabel('z [km]')
ax.set_zlabel('Intensity UI')
plt.grid()
#plt.show()
print((UI).shape)


plt.figure(5)
plt.imshow(UI,extent=[T[0]/T_selected,T[-1]/T_selected,z[0],z[-1]] ,cmap='jet', origin='lower', interpolation='none')
plt.xlim((-5, 5))
plt.title('NLSE Colormap time')
plt.xlabel('T/T0')
plt.ylabel('z [km]')

# plt.imshow(UI.T,extent=[z[0],z[-1],T[0]/T_selected,T[-1]/T_selected] ,cmap='jet')
# plt.ylim((-5, 5))
#plt.contour(zv, yv, UI,cmap='jet')
plt.grid()
plt.show()