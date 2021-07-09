# coding:utf8
import numpy as np
import matplotlib.pyplot as plt
from init_variables import *


T0 = 6E-12 #  duration of input
#for pulse width  --> Dispersive effects at T0 ~ 1ps 
#replaced by usr
T_FWHM_g = 2*np.sqrt(np.log(2))*T0
#For Sech pulses:  
T_FWHM_sech = 2*np.log(1+np.sqrt(2))*T0
N = 8196 #ammount of points 
dt = 100*T0/N #the 100 is to get a grid between -8 and 8 for T/T0   #The number before T0/N sets the time (and freq) frame!
T = np.arange(-N/2, N/2)*dt
T_selected = T0##### Update with a button(?)
T_selected_s = T0
z_initial = 0 #km
zmax = 3.5/3 # km    Values for n/3 -> φNL max = n*π
#gamma = n2*wo/(speed*Aeff)
beta2_initial = 8.3#5.66099
beta3_initial = 10#ps^3/km
gamma_initial = 1 #1/(W*km)
P0 = (3*np.pi)#10E-3
alpha_initial = 0 #dB/km
m0 = 1
#Chirp parameter
C = 0
beta2_initial *= 1E-24
beta3_initial *= 1E-36


###--------SPLIT-STEP---------###
def matp_split_step():
    pulse= Propagation( T_selected, T, m = m0, 
                        C=C, pulsetype = 'Gaussian',
                        solve_type='split_step', 
                        L=zmax, 
                        beta2=beta2_initial,
                        gamma=gamma_initial, 
                        P0=P0)

    plt.figure(1)
    plt.title('NLSE time for a {0} pulse.'.format(pulse.pulsetype))
    plt.plot(T/T_selected, pulse.UI[0], label = 'step = 0')
    plt.plot(T/T_selected, pulse.UI[49], label = 'step at {0} km'.format( '%.3f' %(pulse.z[49])))
    plt.plot(T/T_selected, pulse.UI[-1], label = 'last step at {0} km'.format('%.3f' %(pulse.z[-1])))
    plt.xlabel('T/T0')
    plt.xlim((-10, 10))
    plt.ylabel('Intensity UI')
    plt.legend()
    plt.grid()

    # plt.figure(2)
    # plt.title('Phase time for a {0} pulse.'.format(pulse.pulsetype))
    # plt.plot(T/T_selected, np.angle(pulse.U[0]), label = 'step = 0')
    # plt.plot(T/T_selected, np.angle(pulse.U[49]), label = 'step at {0} km'.format( '%.3f' %(pulse.z[49])))
    # plt.plot(T/T_selected, np.angle(pulse.U[-1]), label = 'last step at {0} km'.format('%.3f' %(pulse.z[-1])))
    # plt.xlabel('T/T0')
    # #plt.xlim((-5, 5))
    # plt.ylabel('Phase \u03C6NL')
    # plt.legend()
    # plt.grid()

    plt.figure(3)
    plt.title('NLSE freq for a {0} pulse.'.format(pulse.pulsetype))
    plt.plot(pulse.W, pulse.UIW[0], label = 'step = 0')
    plt.plot(pulse.W, pulse.UIW[49], label = 'step at {0} km'.format( '%.3f' %(pulse.z[49])))
    plt.plot(pulse.W, pulse.UIW[-1], label = 'last step at {0} km'.format('%.3f' %(pulse.z[-1])))
    plt.xlabel('\u03C9')
    plt.xlim((-2E12, 2E12))
    plt.ylabel('Intensity UI(\u03C9)')
    plt.grid()
    plt.legend()

    yv, zv = np.meshgrid(T/T_selected, pulse.z)
    plt.figure(4)
    ax = plt.axes(projection='3d')
    ax.plot_surface(yv, zv, pulse.UI,cmap='jet')
    ax.set_title('NLSE Colormap time for a {0} pulse.'.format(pulse.pulsetype))
    ax.set_xlabel('T/T0')
    ax.set_ylabel('z [km]')
    ax.set_zlabel('Intensity UI')
    plt.grid()

    wv, zv = np.meshgrid(pulse.W, pulse.z)
    plt.figure(5)
    ax = plt.axes(projection='3d')
    ax.plot_surface(wv, zv, pulse.UIW,cmap='jet')
    ax.set_title('NLSE Colormap freq for a {0} pulse.'.format(pulse.pulsetype))
    ax.set_xlabel('\u03C9-\u03C90')
    ax.set_ylabel('z [km]')
    ax.set_zlabel('Intensity UI')
    plt.grid()
    #print((pulse.UI).shape)

    plt.figure(6)
    plt.imshow(pulse.UI,extent=[T[0]/T_selected,T[-1]/T_selected,pulse.z[0],pulse.z[-1]] ,cmap='jet', origin='lower', interpolation='none', aspect='auto')
    plt.xlim((-5, 5))
    plt.title('NLSE Colormap time for a {0} pulse.'.format(pulse.pulsetype))
    plt.xlabel('T/T0')
    plt.ylabel('z [km]')
    plt.grid()

    plt.figure(7)
    plt.imshow(pulse.UIW,extent=[pulse.W[0],pulse.W[-1],pulse.z[0],pulse.z[-1]] ,cmap='jet', origin='lower', interpolation='none', aspect='auto')
    plt.xlim((-2E12, 2E12))
    plt.title('NLSE Colormap freq for a {0} pulse.'.format(pulse.pulsetype))
    plt.xlabel('\u03C9-\u03C90')
    plt.ylabel('z [km]')
    plt.grid()
    plt.show()
    ###-----------------------###


###---------GVD-----------###
def matp_gvd():
    gvd = Propagation( T_selected, T, m = m0, 
                        C=C, pulsetype = 'Gaussian',
                        solve_type='only_gvd',  
                        beta2=beta2_initial,
                        z0=0,
                        )
    z1 = 2*gvd.compute_LD()
    z2 = 4*gvd.compute_LD()
    gvd1 = Propagation( T_selected, T, m = m0, 
                        C=C, pulsetype = 'Gaussian',
                        solve_type='only_gvd',  
                        beta2=beta2_initial,
                        z0=z1,
                        )
    gvd2 = Propagation( T_selected, T, m = m0, 
                        C=C, pulsetype = 'Gaussian',
                        solve_type='only_gvd',  
                        beta2=beta2_initial,
                        z0=z2,
                        )
    plt.figure(2)
    plt.title('{0} pulse using Eq. 3.2.5, 3.2.6 and 3.2.7.'.format(gvd.pulsetype))
    plt.plot(T/T_selected, gvd.UI)
    plt.plot(T/T_selected, gvd1.UI)
    plt.plot(T/T_selected, gvd2.UI)
    plt.xlabel('T/T0')
    plt.xlim((-8, 8))
    plt.ylabel('UI_Gaussian')
    plt.grid()
    #plt.show()
    plt.figure(3)
    plt.title('Frec for a {0} pulse using Eq. 3.2.5, 3.2.6 and 3.2.7.'.format(gvd.pulsetype))
    plt.plot(gvd.W, gvd.UIW)
    plt.plot(gvd1.W, gvd1.UIW)
    plt.plot(gvd2.W, gvd2.UIW)
    plt.ylabel('|U(z,\u03C9)|^2')
    plt.xlabel('\u03C9-\u03C90')
    plt.grid()
    #plt.show()
    gvd.Gaussian_pulse_GVD()
    gvd1.Gaussian_pulse_GVD()
    gvd2.Gaussian_pulse_GVD()
    plt.figure(4)
    plt.title('Gaussian using Eq. 3.2.7 and 3.2.9')
    plt.plot(T/T_selected, gvd.UI)
    plt.plot(T/T_selected, gvd1.UI)
    plt.plot(T/T_selected, gvd2.UI)
    plt.xlabel('T/T0')
    plt.xlim((-8, 8))
    plt.ylabel('UI_Gaussian')
    plt.grid()
    #plt.show()
    plt.figure(5)
    plt.title('Frec Gaussian using Eq. 3.2.7 and 3.2.9')
    plt.plot(gvd.W, gvd.UIW)
    plt.plot(gvd1.W, gvd1.UIW)
    plt.plot(gvd2.W, gvd2.UIW)
    plt.ylabel('|U(z,\u03C9)|^2')
    plt.xlabel('\u03C9-\u03C90')
    plt.grid()
    plt.show()
###-----------------------###


###---------SPM-----------###
def matp_spm():
    m1 = 1
    m2 = 3
    spm = Propagation( T_selected, T, m = m1, 
                        C=C, pulsetype = 'Gaussian',
                        solve_type='only_spm',  
                        beta2=beta2_initial,
                        gamma=gamma_initial, 
                        P0=P0)
    spm2 = Propagation( T_selected, T, m = m2, 
                        C=C, pulsetype = 'Gaussian',
                        solve_type='only_spm',  
                        beta2=beta2_initial,
                        gamma=gamma_initial, 
                        P0=P0)
    plt.figure(6)
    plt.title('Temporal variation of SPM-induced phase shift \u03A6_NL ') #for a {0} pulse.'.format(spm.pulsetype)
    plt.plot(T/T_selected, spm.Phi_NL, label = 'm = {0}'.format(spm.m))
    plt.plot(T/T_selected, spm2.Phi_NL, label = 'm = {0}'.format(spm2.m))
    plt.xlabel('T/T0')
    plt.xlim((-2.5, 2.5))
    plt.ylabel('\u03A6_NL')
    plt.legend()
    plt.grid()

    plt.figure(7)
    plt.title('frequency chirp \u03B4 \u03C9')
    plt.plot(T/T_selected, spm.delta_w*T_selected, label = 'm = {0}'.format(spm.m))
    plt.plot(T/T_selected, spm2.delta_w*T_selected, label = 'm = {0}'.format(spm2.m))
    plt.xlabel('T/T0')
    plt.xlim((-2.5, 2.5))
    plt.ylabel('\u03B4 \u03C9* T0')
    plt.legend()
    plt.grid()
    plt.show()
###-----------------------###


#Display info:
matp_split_step()
#matp_gvd()
#matp_spm()
