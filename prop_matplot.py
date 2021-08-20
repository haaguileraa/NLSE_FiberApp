# coding:utf8
import numpy as np
import matplotlib.pyplot as plt
from init_variables import *


#T0 = 50E-15 #  duration of input
T0 = 1000E-15 #  duration of input
#for pulse width  --> Dispersive effects at T0 ~ 1ps 
N = 8192 #ammount of points 
dt = 750*T0/N  #The number before T0/N sets the time (and freq) frame!
T = np.arange(-N/2, N/2)*dt
z_initial = 0 #km
zmax = 500#3.5/3 # m    Values for n/3 -> φNL max = n*π
#gamma = n2*wo/(speed*Aeff)
beta2_initial = -10#23.5#8.3#5.66099
beta3_initial = 10#ps^3/km
gamma_initial = 1 #1/(W*km)
P0 = 10#(3*np.pi)#10E-3
alpha_initial = 0 #dB/km
m0 = 1
#Chirp parameter
C = 0
beta2_initial *= 1E-27
beta3_initial *= 1E-39
gamma_initial *= 1E-3
h_step = 0.004

#N=1
'''
T0 = 5E-12 #  duration of input
#for pulse width  --> Dispersive effects at T0 ~ 1ps 
N = 8196 #ammount of points 
dt = 100*T0/N #the 100 is to get a grid between -8 and 8 for T/T0   #The number before T0/N sets the time (and freq) frame!
T = np.arange(-N/2, N/2)*dt
z_initial = 0 #km
zmax = 3.5/3 # km    Values for n/3 -> φNL max = n*π
#gamma = n2*wo/(speed*Aeff)
beta2_initial = 23.5#5.66099
beta3_initial = 10#ps^3/km
gamma_initial = 0.1 #1/(W*km)
P0 = (3*np.pi)#10E-3
alpha_initial = 0 #dB/km
m0 = 1
#Chirp parameter
C = 0
beta2_initial *= 1E-24
beta3_initial *= 1E-36
#'''

###--------SPLIT-STEP---------###
def matp_split_step():
    pulse= Propagation( T0, T, m = m0, 
                        C=C, pulsetype = 'Gaussian',
                        solve_type='split_step', 
                        L=zmax, 
                        beta2=beta2_initial,
                        gamma=gamma_initial, 
                        P0=P0,
                        h_step = h_step)

    plt.figure(1)
    plt.title('NLSE time for a {0} pulse.'.format(pulse.pulsetype))
    plt.plot(T/T0, pulse.UI[0], label = 'step = 0')
    plt.plot(T/T0, pulse.UI[25], label = 'step at {0} m'.format( '%.1f' %(pulse.z[25])))
    plt.plot(T/T0, pulse.UI[-1], label = 'last step at {0} m'.format('%.1f' %(pulse.z[-1])))
    plt.xlabel('T/T0')
    plt.xlim((-15, 15))
    plt.ylabel('Intensity UI')
    plt.legend()
    plt.grid()
    #plt.savefig('G:/Meine Ablage/Bachelorarbeit/Bilder/SSFM/eN1nbt.eps', format='eps')


    # plt.figure(2)
    # plt.title('Phase time for a {0} pulse.'.format(pulse.pulsetype))
    # plt.plot(T/T0, np.angle(pulse.U[0]), label = 'step = 0')
    # plt.plot(T/T0, np.angle(pulse.U[49]), label = 'step at {0} km'.format( '%.3f' %(pulse.z[49])))
    # plt.plot(T/T0, np.angle(pulse.U[-1]), label = 'last step at {0} km'.format('%.3f' %(pulse.z[-1])))
    # plt.xlabel('T/T0')
    # #plt.xlim((-5, 5))
    # plt.ylabel('Phase \u03C6NL')
    # plt.legend()
    # plt.grid()

    plt.figure(3)
    plt.title('NLSE freq for a {0} pulse.'.format(pulse.pulsetype))
    plt.plot(pulse.W, pulse.UIW[0], label = 'step = 0')
    plt.plot(pulse.W, pulse.UIW[25], label = 'step at {0} m'.format( '%.1f' %(pulse.z[25])))
    plt.plot(pulse.W, pulse.UIW[-1], label = 'last step at {0} m'.format('%.1f' %(pulse.z[-1])))
    plt.xlabel('\u03C9-\u03C90')
    plt.xlim((-5E12, 5E12)) #1ps
    #plt.xlim((-0.2E14, 0.2E14)) #1ps
    #plt.xlim((-0.25E14, 0.25E14)) # 500fs
    #plt.xlim((-0.1E15, 0.1E15)) # 50fs
    plt.ylabel('Intensity UI(\u03C9)')
    plt.grid()
    plt.legend()
    #plt.savefig('G:/Meine Ablage/Bachelorarbeit/Bilder/SSFM/eN1nbs.eps', format='eps')
    #plt.savefig('G:/Meine Ablage/Bachelorarbeit/Bilder/SSFM/0pi.eps', format='eps')
    

    yv, zv = np.meshgrid(T/T0, pulse.z)
    

    # plt.figure(4)
    # ax = plt.axes(projection='3d')
    # ax.plot_surface(yv, zv, pulse.UI,cmap='jet')
    # ax.set_title('NLSE Colormap time for a {0} pulse.'.format(pulse.pulsetype))
    # ax.set_xlabel('T/T0')
    # ax.set_ylabel('z [km]')
    # ax.set_zlabel('Intensity UI')
    # plt.grid()

    # wv, zv = np.meshgrid(pulse.W, pulse.z)
    # plt.figure(5)
    # ax = plt.axes(projection='3d')
    # ax.plot_surface(wv, zv, pulse.UIW,cmap='jet')
    # ax.set_title('NLSE Colormap freq for a {0} pulse.'.format(pulse.pulsetype))
    # ax.set_xlabel('\u03C9-\u03C90')
    # ax.set_ylabel('z [km]')
    # ax.set_zlabel('Intensity UI')
    # plt.grid()
    ##print((pulse.UI).shape)

    plt.figure(6)
    plt.imshow(pulse.UI,extent=[T[0]/T0,T[-1]/T0,pulse.z[0],pulse.z[-1]] ,cmap='jet', origin='lower', interpolation='none', aspect='auto')
    plt.xlim((-5, 5))
    plt.title('NLSE Colormap time for a {0} pulse.'.format(pulse.pulsetype))
    plt.xlabel('T/T0')
    plt.ylabel('z [km]')
    plt.grid()
    #plt.savefig('G:/Meine Ablage/Bachelorarbeit/Bilder/SSFM/N1nbt.eps', format='eps')

    plt.figure(7)
    plt.imshow(pulse.UIW,extent=[pulse.W[0],pulse.W[-1],pulse.z[0],pulse.z[-1]] ,cmap='jet', origin='lower', interpolation='none', aspect='auto')
    #plt.xlim((-0.2E14, 0.2E14)) #1ps
    plt.xlim((-5E12, 5E12)) #1ps
    #plt.xlim((-0.25E14, 0.25E14)) # 500fs
    #plt.xlim((-0.1E15, 0.1E15)) # 50fs
    plt.title('NLSE Colormap freq for a {0} pulse.'.format(pulse.pulsetype))
    plt.xlabel('\u03C9-\u03C90')
    plt.ylabel('z [km]')
    plt.grid()
    #plt.savefig('G:/Meine Ablage/Bachelorarbeit/Bilder/SSFM/N1nbs.eps', format='eps')
    plt.show()
###-----------------------###


###---------GVD-----------###
def matp_gvd():
    gvd = Propagation( T0, T, m = m0, 
                        C=C, pulsetype = 'Gaussian',
                        solve_type='only_gvd',  
                        beta2=beta2_initial,
                        z0=0,
                        )
    z1 = 2*gvd.compute_LD()
    z2 = 10*gvd.compute_LD()
    gvd1 = Propagation( T0, T, m = m0, 
                        C=C, pulsetype = 'Gaussian',
                        solve_type='only_gvd',  
                        beta2=beta2_initial,
                        z0=z1,
                        )
    gvd2 = Propagation( T0, T, m = m0, 
                        C=C, pulsetype = 'Gaussian',
                        solve_type='only_gvd',  
                        beta2=beta2_initial,
                        z0=z2,
                        )
    plt.figure(8)
    plt.title('{0} pulse using Eq. 3.2.5, 3.2.6 and 3.2.7.'.format(gvd.pulsetype))
    plt.plot(T/T0, gvd.UI)
    plt.plot(T/T0, gvd1.UI)
    plt.plot(T/T0, gvd2.UI)
    plt.xlabel('T/T0')
    #plt.xlim((-8, 8))
    plt.ylabel('UI_Gaussian')
    plt.grid()
    #plt.show()
    plt.figure(9)
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
    plt.figure(10)
    plt.title('Gaussian using Eq. 3.2.7 and 3.2.9')
    plt.plot(T/T0, gvd.UI)
    plt.plot(T/T0, gvd1.UI)
    plt.plot(T/T0, gvd2.UI)
    plt.xlabel('T/T0')
    #plt.xlim((-8, 8))
    plt.ylabel('UI_Gaussian')
    plt.grid()
    #plt.show()
    plt.figure(11)
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
    spm = Propagation( T0, T, m = m1, 
                        C=C, pulsetype = 'Gaussian',
                        solve_type='only_spm',  
                        beta2=beta2_initial,
                        gamma=gamma_initial, 
                        P0=P0)
    spm2 = Propagation( T0, T, m = m2, 
                        C=C, pulsetype = 'Gaussian',
                        solve_type='only_spm',  
                        beta2=beta2_initial,
                        gamma=gamma_initial, 
                        P0=P0)
    plt.figure(12)
    plt.title('Temporal variation of SPM-induced phase shift \u03A6_NL ') #for a {0} pulse.'.format(spm.pulsetype)
    plt.plot(T/T0, spm.Phi_NL, label = 'm = {0}'.format(spm.m))
    plt.plot(T/T0, spm2.Phi_NL, label = 'm = {0}'.format(spm2.m))
    plt.xlabel('T/T0')
    plt.xlim((-2.5, 2.5))
    plt.ylabel('\u03A6_NL')
    plt.legend()
    plt.grid()
    #plt.savefig('G:/Meine Ablage/Bachelorarbeit/Bilder/spm/shift.eps', format='eps')

    plt.figure(13)
    plt.title('frequency chirp \u03B4 \u03C9')
    plt.plot(T/T0, spm.delta_w*T0, label = 'm = {0}'.format(spm.m))
    plt.plot(T/T0, spm2.delta_w*T0, label = 'm = {0}'.format(spm2.m))
    plt.xlabel('T/T0')
    plt.xlim((-2.5, 2.5))
    plt.ylabel('\u03B4 \u03C9* T0')
    plt.legend()
    plt.grid()
    #plt.savefig('G:/Meine Ablage/Bachelorarbeit/Bilder/spm/chirp.eps', format='eps')
    plt.show()

def matp_spm_sech():
    spm = Propagation( T0, T, 
                        C=C, pulsetype = 'Sech',
                        solve_type='only_spm',  
                        beta2=beta2_initial,
                        gamma=gamma_initial, 
                        P0=P0)
    plt.figure(14)
    plt.title('Temporal variation of SPM-induced phase shift \u03A6_NL (Sech Pulse)') #for a {0} pulse.'.format(spm.pulsetype)
    plt.plot(T/T0, spm.Phi_NL)
    plt.xlabel('T/T0')
    plt.xlim((-5, 5))
    plt.ylabel('\u03A6_NL')
    plt.grid()
    #plt.savefig('G:/Meine Ablage/Bachelorarbeit/Bilder/spm/shift_sech.eps', format='eps')

    plt.figure(15)
    plt.title('frequency chirp \u03B4 \u03C9 (Sech Pulse)')
    plt.plot(T/T0, spm.delta_w*T0)
    plt.xlabel('T/T0')
    plt.xlim((-5, 5))
    plt.ylabel('\u03B4 \u03C9* T0')
    plt.grid()
    #plt.savefig('G:/Meine Ablage/Bachelorarbeit/Bilder/spm/chirp_sech.eps', format='eps')
    plt.show()
###-----------------------###




#Display info:
matp_split_step()
#matp_gvd()
#matp_spm()
#matp_spm_sech()
