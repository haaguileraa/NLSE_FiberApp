# coding:utf8

import numpy as np
#from scipy import integrate
import matplotlib.pyplot as plt
from Functions import Gaussian_pulse_GVD, incident_field



#------------------- Grid: --------------------#
#Input pulse
#Tmax = 5 ~ 10ps
T0 = 1E-12 # for pulse Initial width T0 --> Dispersive effects at T0 ~ 1ps  S.64
# half-width (at 1/e-intensity point)

T_FWHM = 2*np.sqrt(np.log(2))*T0

N = 8196 #ammount of points 

T_selected = T0

dt = 25*T_selected/N #the 16 is to get a grid between -8 and 8 for T/T_selected 
#dt = 8*T_FWHM/N


#T = np.linspace(-8*T_selected,8*T_selected,N)
T = np.arange(-N/2, N/2)*dt #np.arange(-N/2, N/2) mit len = N
# print(T/T_selected)
# print(T[1]-T[0])
# print(1/(T[-1]-T[0]))
# #-----------------------------------------------#
# print(len(T/T_selected))

beta2 = 20 #ps^2/km  #S.65
#The use of these values in Eq. (3.1.5) shows that the dispersive and nonlinear effects are negligible for L < 50 km if T_selected > 100 ps and P0 ~ 1 mW

LD = (T_selected**2)/np.absolute(beta2)
z = 0
z1 = 2*LD
z2 = 4*LD


#---------Using Eq 3.2.7 and 3.2.9--------#

# _, UI1 = Gaussian_pulse_GVD(z,T,T_selected, beta2)
# _, UI2 = Gaussian_pulse_GVD(z1,T,T_selected, beta2)
# _, UI3 = Gaussian_pulse_GVD(z2,T,T_selected, beta2)

#print(T)
# plt.figure(1)
# plt.title('Gaussian using Eq. 3.2.7 and 3.2.9')
# plt.plot(T/T_selected, UI1)
# plt.plot(T/T_selected, UI2)
# plt.plot(T/T_selected, UI3)
# plt.xlabel('T/T_selected')
# plt.xlim((-8, 8))
# plt.ylabel('UI')
# plt.grid()

#For Sech pulses:  
#Chirp parameter
C = 0 
T_FWHM_sech = 2*np.log(1+np.sqrt(2))*T_selected 



#---------Using Eq 3.2.5 and 3.2.6--------#

#Definition Gaussian Pulse's U(0,T) 
UT0_Gaussian = np.exp(-(T**2)/(2*T_selected**2)) #T0 replaced by T_FWHM

_, UI4, UW4, W4 = incident_field(beta2, z, T_selected, T, pulse = 'Gaussian')
_, UI5, UW5, W5 = incident_field(beta2, z1, T_selected, T, pulse = 'Gaussian')
_, UI6, UW6, W6 = incident_field(beta2, z2, T_selected, T, pulse = 'Gaussian')


plt.figure(2)
plt.title('Gaussian using Eq. 3.2.5, 3.2.6 and 3.2.7')
plt.plot(T/T_selected, UI4)
plt.plot(T/T_selected, UI5)
plt.plot(T/T_selected, UI6)
plt.xlabel('T/T0')
plt.xlim((-8, 8))
plt.ylabel('UI_Gaussian')
plt.grid()


UW_4 = np.absolute(UW4)**2
UW_5 = np.absolute(UW5)**2
UW_6 = np.absolute(UW6)**2
plt.figure(4)
plt.title('Frec Gaussian using Eq. 3.2.5, 3.2.6 and 3.2.7')
plt.plot(W4, UW_4)
# plt.plot(W5, UW_5)
# plt.plot(W6, UW_6)
#plt.plot(W4, UW4.real, 'o')
# plt.plot(W5, UW5.real, 'o')
# plt.plot(W6, UW6.real, 'o')
# plt.xlabel('T/T0')
#plt.xlim((-0.5, 0.5))
plt.ylabel('UI_Gaussian')
plt.grid()
plt.show()

# #Definition Sech Pulse's U(0,T) 
# UT0_sech = 1/(np.cosh(T/T_selected))*np.exp(-(1j*C*T**2)/(2*T_selected**2))  #T_selected replaced by T_FWHM

# _, UI7, UW7, W7 = incident_field(beta2, z, T_selected, T, pulse = 'Sech')
# _, UI8, UW8, W9 = incident_field(beta2, z1, T_selected, T, pulse = 'Sech')
# _, UI9, UW9, W8 = incident_field(beta2, z2, T_selected, T, pulse = 'Sech')


# plt.figure(3)
# plt.title('Sech pulse GVD using Eq. 3.2.5, 3.2.6 and 3.2.21')
# plt.plot(T/T_selected, UI7)
# plt.plot(T/T_selected, UI8)
# plt.plot(T/T_selected, UI9)
# plt.xlabel('T/T0')
# plt.xlim((-8, 8))
# plt.ylabel('UI_Sech')
# plt.grid()
# plt.show()
