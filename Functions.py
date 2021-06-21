# coding:utf8

import numpy as np
import matplotlib.pyplot as plt
import dash_core_components as dcc
from scipy.integrate import cumtrapz, solve_ivp
import plotly.graph_objects as go
from numpy.core.numeric import Inf   

#ifft = +
#fft= -

##------------------------For GVD------------------------##
def incident_field(beta2, z, T_, T, pulse = 'Gaussian', C = 0):
    dT = T[1]-T[0]
    N = len(T)
    if pulse == 'Gaussian':
        UT0 = np.exp(-(T**2)/(2*T_**2))
    elif pulse == 'Sech':
        UT0 = 1/(np.cosh(T/T_))*np.exp(-(1j*C*T**2)/(2*T_**2))
    else:
        raise ValueError("Pulse must be 'Sech' or 'Gaussian'.")
    UW0 =np.fft.fftshift(np.fft.ifft(UT0))
    W = 2*np.pi* np.fft.fftfreq(N, dT)
    W = np.fft.fftshift(W)
    UW = UW0*np.exp((1j*beta2*W**2*z/2))
    UT = np.fft.fft(np.fft.fftshift(UW))
    UI = np.absolute(UT)**2
    return UT, UI, UW, W

def Gaussian_pulse_GVD(z,T,T_, beta2):#T_ replaced by either T0 or TFWHM
    UT = T_/(np.sqrt(T_**2 - 1j*beta2*z))*np.exp(-T**2/(2*(T_**2 - 1j*beta2*z)))
    UI = np.absolute(UT)**2
    return UT, UI


def Sech_pulse_GVD(z,T,T_, C, beta2):
    U0 = 1/(np.cosh(T/T_))*np.exp(-(1j*C*T**2)/(2*T_**2))
    UT = T_/(np.sqrt(T_**2 - 1j*beta2*z))*np.exp(-T**2/(2*(T_**2 - 1j*beta2*z)))
    UI = np.absolute(UT)**2
    return UT, UI
##-------------------------------------------------------##

##-----------NUMERICAL METHODS-----------##
#--------derivatives--------#
def derivative(f,x0,method='central',dx=0.01):
    if method == 'central':
        return (f(x0 + dx) - f(x0 - dx))/(2*dx)
    elif method == 'forward':
        return (f(x0 + dx) - f(x0))/dx
    elif method == 'backward':
        return (f(x0) - f(x0 - dx))/dx
    else:
        raise ValueError("Method must be 'central', 'forward' or 'backward'.") #taken and edited from https://www.math.ubc.ca/~pwalls/math-python/differentiation/differentiation/
#--------------------------#


#Mid-point Method--------#
def mid_step(x0, f, T, *args):
    dt = T[1]-T[0]
    k = len(T)
    X = np.zeros(k)
    X[0]=x0
    for i in range(1,k,1):
        X[i] = X[i-1] + dt*f(T[i-1]+dt/2, *args)
    return X
#------------------------#

##------------------------------------##




##------------------------For SPM------------------------##

#--frequency chirp δω for Gaussian Pulse--#
def delta_g(T, T_, Leff, LNL, m = 1):
    delta_w = ((2*m*Leff)/(T_*LNL))*((T/T_)**(2*m-1))*np.exp(-(T/T_)**(2*m)) #eq (4.1.9)
    return delta_w
#-----------------------------------------#

def incident_field_spm(Leff, gamma, P0, T_, T, pulse = 'Gaussian', m = 1, C=0):
    LNL= 1/(gamma*P0)
    dT = T[1]-T[0]
    if pulse == 'Gaussian':
        UT0 = np.exp(-(T**2)/(2*T_**2))
        if m <= 1:
            Phi_NL = (Leff/LNL)*np.absolute(UT0)**2
            delta_w = delta_g(T, T_, Leff, LNL)
        else: 
            delta_w = delta_g(T, T_, Leff, LNL, m)
            #for array-like data:
            #Phi_NL = -cumtrapz(delta_w, T, initial=0) #Default 'initial' is None, which means no value at x[0] 
            #for function-like input:
            Phi_NL = -mid_step(0, delta_g, T, T_, Leff, LNL, m)

    elif pulse == 'Sech':
        phinl = lambda T: (Leff/LNL)*np.absolute(1/(np.cosh(T/T_))*np.exp(-(1j*C*T**2)/(2*T_**2)))**2
        Phi_NL = phinl(T)
        delta_w = -derivative(phinl, T, dx=dT)

    else:
        raise ValueError("Pulse must be 'Sech' or 'Gaussian'.")

    return Phi_NL, delta_w
    ##-------------------------------------------------------##


    ##Split-Step-Fourier-Method:

def U_noGVD(z,val,gamma,P0,T0,alpha,beta2): 
    '''
    dU(T,z)/dz = N*U(T,z) = i*exp(-alpha*z)/LNL*(|U(T,z)|²)*U(T,z)
    '''
    U = val
    # definition of f = N*U
    #f = 1j*(gamma*P0*T0**2/np.absolute(beta2))*(np.absolute(U)**2)*U
    #f = 1j*(np.exp(-alpha*z)*gamma*P0*T0**2)*(np.absolute(U)**2)*U
    f = 1j*(np.exp(-alpha*z)*gamma*P0*T0**2/np.absolute(beta2))*(np.absolute(U)**2)*U
    return f



def split_step(beta2, T_, T, L, gamma, P0,  beta3=0, alpha = 0, pulse = 'Gaussian', m = 1, C=0):
    dT = T[1]-T[0]  #step size
    points = len(T)      #number of points
    #Length for nonlinearities:
    if gamma != 0 and P0 != 0 :
        LNL= 1/(gamma*P0)
    else:
        LNL = Inf #####BE CAREFUL!!!!!!!!!!!

    #Dispersion Length:
    if beta2 != 0:
        LD = ((T_)**2)/np.absolute(beta2)
    else:
        LD = Inf #####BE CAREFUL!!!!!!!!!!!
    
    #Effective Length:
    if alpha !=0:
        Leff = (1 - np.exp(-alpha*L))/alpha
    else: 
        Leff = max(LNL, LD) #####BE CAREFUL!!!!!!!!!!!
  
    
    #phi_maxnl = gamma*P0*Leff #####BE CAREFUL!!!!!!!!!!! 0*inf...

    if pulse == 'Gaussian':
        UT0 = (np.exp(-(T**2)/(2*T_**2))).astype(complex) #dtype = 'complex' in order to have complex values on solve_ivp
    elif pulse == 'Sech':
        UT0 = (1/(np.cosh(T/T_))*np.exp(-(1j*C*T**2)/(2*T_**2))).astype(complex)
    else:
        raise ValueError("Pulse must be 'Sech' or 'Gaussian'.")
    UW0 =np.fft.fftshift(np.fft.ifft(UT0))
    W = 2*np.pi* np.fft.fftfreq(points, dT)
    W = np.fft.fftshift(W)

    # for normalized NLSE:
    #z = z/LD #normalized distance 
    #tau = T/T_  #normalized time variables  This is following the normalization made on chapter 4 (Not used here)
    #δ²U/z² = [D + N]*U
    #U(z+h,T) = exp(h*D)*exp(h*N)*U(z,T) #first step with z = 0 
    zmax= L  #km  #Max distance that will be evaluated
    h = zmax*0.001 #step size km -> depends on how long is the fiber we are working with
    factork = 0.25 #100 for 1/h a points for each zmax km #just a factor to calculate how many points will be used during SSFM
    k = int(zmax*factork/h)  #how will be z splited (z in km) 
    print('Z max.: ',zmax, 'km')
    print('h size: ',h, 'km')
    print('LD: ',LD, 'km')
    print('LNL: ',LNL, 'km')
    print('N = sqrt(LD/LNL) = : ',np.sqrt(LD/LNL))
    print(k,'number of points evaluated')
    size_array = 100 #Size of the array where we are going to save the values of U(T) and U(w)
    U = np.zeros((size_array,points), dtype=complex) # U(z,T) is divided by certain ammount of steps and it's time length is N for each one of the steps
    UW = np.zeros((size_array,points), dtype=complex)
    U[0] = UT0  #first value for U(0,T)  
    UW[0] = UW0 #first value for U(0,W)  
    #for general def: 
    Dw = -(1j*beta2/2)*(1j*W)**2 + (beta3/6)*(1j*W)**3 - alpha/2 #D = -1j*beta2/2*(δ²/T²) + beta3/6*(δ³/T³) - alpha/2
    #create a way to extract the betas from an array and take them to the frequency domain
    
    #---This is smth in development---#
    #for D on SPM chapter:
    #Dw = -(1j*np.sign(beta2)/(2))*(1j*W)**2 
    #D = np.fft.fft(np.fft.fftshift(Dw))
    #expdw = np.fft.fftshift(np.fft.ifft(np.exp(h*Dw)))
    #expdt = np.fft.fft(np.fft.fftshift(expdw)) #exp(hD)B(z,T) = FT^-1(exp(hD(iw)))FT*B(z,T) 
    #---------------------------------#

    #for i in np.linspace(1,zmax*1000-h,k):  if we pre-defined the resolution (or ammount of steps k) independent of the step size h
    #initial values for split-step
    z=0 # initial value for saving the field
    z_array = np.zeros((size_array))#for saving the z-point where U(T) and U(w) were saved
    m = k//size_array 
    j=1
    #Skip until "Continue":
    '''

    for i in range(0,k-1,1): 
        z = h*i/factork #km
        if i > k-10:
            print('z+h evaluated: ',z+h, 'km')
        #N = 1j*(np.exp(-alpha*(z+h))/LNL)*(np.absolute(Ut)**2)
        #Ut = np.fft.fft(np.fft.fftshift(np.exp((z+h)*Dw)*Uw))
        Ut = np.fft.fft(np.fft.fftshift(np.exp((h)*Dw)*Uw))
        #for N on SPM chapter
        N = 1j*(np.exp(-alpha*(z+h))/LNL)*(np.absolute(Ut)**2)   #NL = i*e^(-alpha*z)*|U|² 
        #N = 1j*(gamma)*(np.absolute(Ut)**2)   #NL = i*e^(-alpha*z)*|U|² 
        #Ut *= np.exp((z+h)*N)
        Ut *= np.exp((h)*N)
        Uw = np.fft.fftshift(np.fft.ifft(Ut))
        #Symmetric        
        # Uw *= np.exp((h/2)*Dw)
        # Ut  = np.fft.fft(np.fft.fftshift(Uw))
        # N = 1j*(np.exp(-alpha*(z+h))/ LNL)*(np.absolute(Ut)**2)
        # Ut *= np.exp(h*N)
        # Uw = np.fft.fftshift(np.fft.ifft(Ut))*np.exp((h/2)*Dw)
        # Ut  = np.fft.fft(np.fft.fftshift(Uw))
        if i == m:
            z_array[j] = z
            U[j] = Ut
            UW[j] = Uw#np.fft.fftshift(np.fft.ifft(U[j]))
            #for the next step
            j +=1
            m = j*(k//size_array)
    
    z_array[-1] = z
    U[-1] = Ut
    UW[-1] = Uw
    UI = np.absolute(U)**2
    
    print('last z evaluated: ',z, 'km')
    print('last z saved: ',z_array[-1], 'km')
    UIW = np.absolute(UW)**2
    return U, UI, UW, UIW, W, z_array

    '''
    #with Runge-Kutta
    '''
    #first from 0 to z_max 
    z_array = np.linspace(0,zmax, size_array)
    sol = solve_ivp ( U_noGVD, # function to integrate
                    (0 , zmax) , # duration [ z_init - z_end ]
                    UT0 , # initial conditions
                    method ='RK23', # method of integration
                    t_eval =z_array, # points where the sol. is stored
                    rtol =1e-8 , # relative accuracy
                    args =[gamma]) # arguments for the function
    U = sol.y.T
    for utemp in U:
        UW[j] = np.fft.fftshift(np.fft.ifft(utemp))
        j+=1
    '''
    #--Continue:

    #Now from z to z+h n times till z_max
    Utemp = UT0[:]
    for i in range(0,k-1,1): 
        z = h*i/factork #km #actual z to evaluate
        sol = solve_ivp ( U_noGVD, # function to integrate
                    (z, z+h) , # duration [ z_init - z_end ]
                    Utemp, # initial conditions
                    method ='RK23', # method of integration
                    #t_eval =np.linspace(z, z+h, size_array), # points where the sol. is stored
                    rtol =1e-8 , # relative accuracy
                    args =[gamma, P0, T_, alpha, beta2]) # arguments for the function
        Utemp = sol.y.T[-1] #value at z+h
        Uwtemp = np.fft.fftshift(np.fft.ifft(Utemp))
        Utemp = np.fft.fft(np.fft.fftshift(np.exp((h)*Dw)*Uwtemp))
        Uwtemp = np.fft.fftshift(np.fft.ifft(Utemp))
        if i == m and j< size_array:
            z_array[j] = z
            U[j] = Utemp
            UW[j] = Uwtemp
            j +=1 #for the next step
            m = (j)*(k//size_array)
            

    z_array[-1] = z
    U[-1] = Utemp
    UW[-1] = Uwtemp
    UI = np.absolute(U)**2
    print('last z evaluated: ',z, 'km')
    print('last z saved: ',z_array[-1], 'km')
    UIW = np.absolute(UW)**2
    return U, UI, UW, UIW, W, z_array


def plot_prop(UI, t,z):
    propagation = go.Figure(data=[go.Heatmap(
                x = t,#np.sort(t),
                y = z,#np.sort(z),
                z = UI,
                #type = 'heatmap',
                zsmooth = "best", #to avoid line at the end of the plot
                colorscale = 'Jet')])
    propagation.update_layout(
                            width=600, height=600,
                            yaxis=dict(range=[z[0], z[-1]],title='z [km]', ), 
                            xaxis=dict(range=[-8, 8],title='T/T0',), 
                            )
    #             shapes = [dict(type="rect",
    #                                             xref="x",
    #                                             yref="y",
    #                                             x0=t[0],
    #                                             y0=z[0],
    #                                             x1=t[-1],
    #                                             y1=z[-1],
    #                                             line_color="LightSeaGreen",
    #                                             line_dash='dash',),],)

    return propagation


###ON PROGRESS it will be finished after all methods are verified
class propagation():
    def __init__(self, T_, T, pulsetype = 'Gaussian' ,soltype = 'split_step'):
        self.pulse = pulsetype
        self.solution_method = soltype
        self.T_selected = T_
        self.T = T
        self.initParam()

    def initParam(self):
        self.dt = self.T[1]-self.T[0]
        self.N = len(self.T)

    def plot_gvd(self, beta2, z, colors, C=0):
        _, UI, UW, W = incident_field(beta2, z, self.T_selected, self.T, pulse = self.pulse, C = C)
        UW_ = np.absolute(UW)**2
        Uz = go.Scatter(x=self.T/self.T_selected ,y=UI, name = self.pulse,
                         line=dict(color=colors['even']))
        Uf = go.Scatter(x=W ,y=UW_, name = self.pulse,
                         line=dict(color=colors['even']))
        return Uz, Uf
    def plot_spm(self, Leff, gamma, P0, colors, m = 1, C = 0):
        Phi, dw = incident_field_spm(Leff, gamma, P0, self.T_selected, self.T, self.pulse, m = 1, C=0)

        phase = go.Scatter(x=self.T/self.T_selected ,y=Phi, name = self.pulse,
                         line=dict(color=colors['even']))
        chirp = go.Scatter(x=self.T/self.T_selected ,y=dw*self.T_selected, name = self.pulse,
                         line=dict(color=colors['even']))
        return phase, chirp


        
# class gvd(propagation):
# class spm(propagation):
# class ssfm(propagation):
