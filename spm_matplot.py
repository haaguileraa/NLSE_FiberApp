# coding:utf8

import numpy as np
#from scipy import integrate
import matplotlib.pyplot as plt
from Functions import Gaussian_pulse_GVD, incident_field, incident_field_spm



#------------------- Grid: --------------------#
#Input pulse
#Tmax = 5 ~ 10ps
T0 = 200E-12 # for pulse Initial width T0 --> Dispersive effects at T0 ~ 1ps  S.64
# half-width (at 1/e-intensity point)

T_FWHM = 2*np.sqrt(np.log(2))*T0

N = 8196 #ammount of points 

T_selected = T0

dt = 25*T_selected/N #the 16 is to get a grid between -8 and 8 for T/T_selected 
#dt = 8*T_FWHM/N


#T = np.linspace(-8*T_selected,8*T_selected,N)
T = np.arange(-N/2, N/2)*dt #np.arange(-N/2, N/2) mit len = N
#-----------------------------------------------#


#The use of these values in Eq. (3.1.5) shows that the dispersive and nonlinear effects are negligible for L < 50 km if T_selected > 100 ps and P0 ~ 1 mW
C = 0 
T_FWHM_sech = 2*np.log(1+np.sqrt(2))*T0



#gamma = n2*wo/(speed*Aeff)
gamma = 2 #1/(W*km)
P0 = 10E-3
LNL= 1/(gamma*P0)
Leff = LNL

Phi, dw = incident_field_spm(Leff, gamma, P0, T_selected, T, pulse = 'Gaussian', m = 1, C=0)
Phi2, dw2 = incident_field_spm(Leff, gamma, P0, T_selected, T, pulse = 'Gaussian', m = 3, C=0)



plt.figure(1)
plt.title('Temporal variation of SPM-induced phase shift \u03A6_NL')
plt.plot(T/T_selected, Phi, label = 'm = 1')
plt.plot(T/T_selected, Phi2, label = 'm = 3')
plt.xlabel('T/T0')
plt.xlim((-2.5, 2.5))
plt.ylabel('\u03A6_NL')
plt.legend()
plt.grid()


plt.figure(2)
plt.title('frequency chirp \u03B4 \u03C9')
plt.plot(T/T_selected, dw*T_selected, label = 'm = 1' )
plt.plot(T/T_selected, dw2*T_selected, label = 'm = 3')
plt.xlabel('T/T0')
plt.xlim((-2.5, 2.5))
plt.ylabel('\u03B4 \u03C9* T0')
plt.legend()
plt.grid()
plt.show()



Phi_s, dw_s = incident_field_spm(Leff, gamma, P0, T_selected, T, pulse = 'Sech', C=0)


plt.figure(3)
plt.title('Temporal variation of SPM-induced phase shift \u03A6_NL')
plt.plot(T/T_selected, Phi_s)
plt.xlabel('T/T0')
plt.xlim((-5, 5))
plt.ylabel('\u03A6_NL')
plt.grid()

plt.figure(4)
plt.title('frequency chirp \u03B4 \u03C9')
plt.plot(T/T_selected, dw_s*T_selected)
plt.xlabel('T/T0')
plt.xlim((-5, 5))
plt.ylabel('\u03B4 \u03C9* T0')
plt.grid()
plt.show()