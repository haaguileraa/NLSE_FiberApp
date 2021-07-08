# coding:utf8
import numpy as np
import matplotlib.pyplot as plt
import dash_core_components as dcc
from scipy.integrate import cumtrapz, solve_ivp
import plotly.graph_objects as go
from numpy.core.numeric import Inf  
import pandas as pd


##-----------NUMERICAL METHODS-----------##
#--------Derivatives--------#
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
#----Mid-point Method------#
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



##---------------------SPM---------------------##
#--frequency chirp δω for Gaussian Pulse--#
def delta_g(T, T_, Leff, LNL, m = 1):
    delta_w = ((2*m*Leff)/(T_*LNL))*((T/T_)**(2*m-1))*np.exp(-(T/T_)**(2*m)) #eq (4.1.9) Agrawal
    return delta_w
#-----------------------------------------#
##---------------------------------------------##


#--SPLIT STEP FUNCTIONS--#
def U_noGVD(z,val,gamma,P0,T0,alpha,beta2, beta3): 
    '''
    dU(T,z)/dz = N*U(T,z) = i*exp(-alpha*z)/LNL*(|U(T,z)|²)*U(T,z)
    '''
    tempU = val
    # definition of f = N*U
    #f = 1j*(gamma*P0*T0**2/np.absolute(beta2))*(np.absolute(U)**2)*U
    #f = 1j*(np.exp(-alpha*z)*gamma*P0*T0**2)*(np.absolute(U)**2)*U
    if beta2 != 0:
        f = 1j*(np.exp(-alpha*z)*gamma*P0*T0**2/np.absolute(beta2))*(np.absolute(tempU)**2)*tempU
    elif beta3 != 0 and beta2 == 0:
        f = 1j*(np.exp(-alpha*z)*gamma*P0*T0**3/np.absolute(beta3))*(np.absolute(tempU)**2)*tempU
    else:
        f  = 1j*(np.exp(-alpha*z)*gamma*P0)*(np.absolute(tempU)**2)*tempU
    return f
#-------------------------#


class Pulse:
    def __init__(self,  T_, T,  m = 1, C=0, pulsetype = 'Gaussian'):
        self.pulsetype = pulsetype
        #self.solution_method = soltype
        self.T_selected = T_
        self.T = T
        self.m = m
        self.C = C
        self.initParam()

#L, beta2, gamma, P0,  beta3=0, alpha = 0,

    def initParam(self):
        self.dT = self.T[1]-self.T[0] # step_size
        self.points = len(self.T)
        #Normalized Slow-Varying Envelope
        if self.pulsetype == 'Gaussian':
            self.UT0 = (np.exp(-((1+1j*self.C)/2)*(self.T**2)/(self.T_selected**2))).astype(complex) #dtype = 'complex' in order to have complex values on solve_ivp 
            #No chirp, no Super Gaussian (np.exp(-(T**2)/(2*T_**2))).astype(complex)
        elif self.pulsetype == 'SGaussian':
            self.UT0 = (np.exp(-((1+1j*self.C)/2)*((self.T**2)/(self.T_selected**2))**(2*self.m))).astype(complex) #dtype = 'complex' in order to have complex values on solve_ivp 
            #No chirp, no Super Gaussian (np.exp(-(T**2)/(2*T_**2))).astype(complex)
        elif self.pulsetype == 'Sech':
            self.UT0 = (1/(np.cosh(self.T/self.T_selected))*np.exp(-(1j*self.C*self.T**2)/(2*self.T_selected**2))).astype(complex)
        else:
            raise ValueError("Pulse must be 'Sech', 'Gaussian' or 'SGaussian'.")
        self.UW0 =np.fft.fftshift(np.fft.ifft(self.UT0))
        self.W = 2*np.pi* np.fft.fftfreq(self.points, self.dT)
        self.W = np.fft.fftshift(self.W)








class Propagation(Pulse):
    def __init__(self,T_, T, solve_type= 'incident_field', L=0.1, beta2=0, gamma=0, P0=0,  beta3=0, loss = 0, pulsetype = 'Gaussian', m = 1, C=0, h_step = 0.001, factork = 0.25, size_array = 100):
        Pulse.__init__(self, T_, T, m=m, C=C, pulsetype=pulsetype) 
        self.zmax = L #[km]  #Max distance that will be evaluated
        self.beta2 = beta2 #[ps²/km]
        self.gamma = gamma #[m/W]
        self.P0 = P0 #[W]
        self.beta3= beta3 #[ps³/km]
        self.loss = loss  #[dB/km]
        self.h_step = h_step
        self.factork = factork #100 for 1/h a points for each zmax km #just a factor to calculate how many points will be used during SSFM
        self.size_array = size_array #Size of the array where we are going to save the values of U(T) and U(w)
        self.U = np.zeros((self.size_array,self.points), dtype=complex) # U(z,T) is divided by certain ammount of steps and it's time length is N for each one of the steps
        self.UW = np.zeros((self.size_array,self.points), dtype=complex)
        self.solve_type = solve_type

        if self.solve_type == 'only_gvd':
            self.incident_field()
        elif self.solve_type == 'only_spm':
            self.incident_field_spm()
        elif self.solve_type == 'split_step':
            if self.beta2  == 0 or self.gamma  == 0 or self.P0 == 0:
                print('Either L, beta2, gamma or P0 are zero, this could lead to problems with the calculation of the current method selected!')
            if self.zmax  == 0:
                raise ValueError('Fiber Length cannot be 0 for de desired solution')       
            self.split_step()
        else:
            raise ValueError("solution method must be 'only_gvd', 'only_spm' or 'split_step'.")

    #Lengths
    def compute_LD(self):#Dispersion Length:
        if self.beta2 != 0:
            return ((self.T_selected)**2)/np.absolute(self.beta2)
        else:
            return Inf 

    def compute_LNL(self):#Length for nonlinearities:
        if self.gamma != 0 and self.P0 != 0 :
            return 1/(self.gamma*self.P0)
        else:
            return Inf 

    def compute_Leff(self): #Effective Length:
        if self.alpha !=0:
            return (1 - np.exp(-self.alpha*self.L))/self.alpha
        else: 
            return self.zmax 

    def compute_phimax_nogvd(self):
        Leff = self.compute_Leff()
        return self.gamma*self.P0*Leff #####BE CAREFUL!!!!!!!!!!! 0*inf...


    ##------------------------For GVD------------------------##
    def incident_field(self,z0):  #Function using Eq. 3.2.5 and 3.2.6
        self.UW = self.UW0*np.exp((1j*self.beta2*self.W**2*z0/2))
        self.UT = np.fft.fft(np.fft.fftshift(self.UW))
        self.UI = np.absolute(self.UT)**2
        return self.UT, self.UI, self.UW, self.W

    def Gaussian_pulse_GVD(self,z0): #Function using Eq. 3.2.7 and 3.2.9 Agrawal
        self.UT = self.T_selected/(np.sqrt(self.T_selected**2 - 1j*self.beta2*z0))*np.exp(-self.T**2/(2*(self.T_selected**2 - 1j*self.beta2*z0)))
        self.UI = np.absolute(self.UT)**2
        return self.UT, self.UI


    def Sech_pulse_GVD(self,z0):#Not used 
        self.U0 = 1/(np.cosh(self.T/self.T_selected))*np.exp(-(1j*self.C*self.T**2)/(2*self.T_selected**2))
        self.UT = self.T_selected/(np.sqrt(self.T_selected**2 - 1j*self.beta2*z0))*np.exp(-self.T**2/(2*(self.T_selected**2 - 1j*self.beta2*self.z)))
        self.UI = np.absolute(self.UT)**2
        return self.UT, self.UI
    ##-------------------------------------------------------##
    ##-------------------------------------------------------##


    ##------------------------For SPM------------------------##

    def incident_field_spm(self):
        #Leff = self.compute_Leff()
        LNL = self.compute_LNL()
        Leff = LNL # Normalized to get the plots from the book
        if self.pulsetype == 'Gaussian':
            if self.m <= 1:
                Phi_NL = (Leff/LNL)*np.absolute(self.UT0)**2
                delta_w = delta_g(self.T, self.T_selected, Leff, LNL)
            else: 
                delta_w = delta_g(self.T, self.T_selected, Leff, LNL, self.m)
                #for array-like data
                #Phi_NL = -cumtrapz(delta_w, T, initial=0) #Default 'initial' is None, which means no value at x[0] 
                #for function-like input:
                Phi_NL = -mid_step(0, delta_g, self.T, self.T_selected, Leff, LNL, self.m)

        elif self.pulsetype == 'Sech':
            phinl = lambda t: (Leff/LNL)*np.absolute(1/(np.cosh(t/self.T_selected))*np.exp(-(1j*self.C*t**2)/(2*self.T_selected**2)))**2
            Phi_NL = phinl(self.T)
            delta_w = -derivative(phinl, self.T, dx=self.dT)
        else:
            raise ValueError("Pulse must be 'Sech' or 'Gaussian'.")

        return Phi_NL, delta_w
    ##-------------------------------------------------------##
    ##-------------------------------------------------------##

    ##---------------------SPLIT-STEP------------------------##
    def split_step(self):
        self.alpha = np.log(10**(self.loss/10)) #log y = ln y / ln 10 
        h = self.zmax*self.h_step #step size km -> depends on how long is the fiber we are working with
        k = int(self.zmax*self.factork/h)  #how will be z splited (z in km) 
        print('\nZ max.: ',self.zmax, 'km')
        print('h size: ',h, 'km')
        print('LD: ',self.compute_LD(), 'km')
        print('LNL: ',self.compute_LNL(), 'km')
        print('N = sqrt(LD/LNL) = : ',np.sqrt(self.compute_LD()/self.compute_LNL()))
        print(k,'number of points evaluated')
        #initial values for split-step
        self.U[0] = self.UT0  #first value for U(0,T)  
        self.UW[0] = self.UW0 #first value for U(0,W)  
        #for general def: 
        Dw = -(1j*self.beta2/2)*(1j*self.W)**2 + (self.beta3/6)*(1j*self.W)**3 - self.alpha/2 #D = -1j*beta2/2*(δ²/T²) + beta3/6*(δ³/T³) - alpha/2
        #TODO: create a way to extract the betas from an array and take them to the frequency domain
        
        z=0 # initial value for saving the field
        self.z = np.zeros((self.size_array))#array for saving the nth-point U(T,zn) and U(w,zn)
        mi = k//self.size_array #just a counter
        j=1 #just a counter
        Utemp = self.UT0[:] #variable to update

        #for i in np.linspace(1,zmax*1000-h,k):  if we pre-defined the resolution (or ammount of steps k) independent of the step size h
        for i in range(0,k-1,1): # from z to z+h k times until z_max is reached
            ztemp = h*i/self.factork #km #actual z to evaluate
            sol = solve_ivp (U_noGVD, # function to integrate
                        (z, z+h) , # duration [ z_init - z_end ]
                        Utemp, # initial conditions
                        method ='RK23', # method of integration
                        #t_eval =np.linspace(z, z+h, size_array), # points where the sol. is stored
                        rtol =1e-8 , # relative accuracy
                        args =[self.gamma, self.P0, self.T_selected, self.alpha, self.beta2, self.beta3]) # arguments for the function
            Utemp = sol.y.T[-1] #value at z+h
            Uwtemp = np.fft.fftshift(np.fft.ifft(Utemp))
            Utemp = np.fft.fft(np.fft.fftshift(np.exp((h)*Dw)*Uwtemp))
            Uwtemp = np.fft.fftshift(np.fft.ifft(Utemp))
            if i == mi and j< self.size_array:
                self.z[j] = ztemp
                self.U[j] = Utemp
                self.UW[j] = Uwtemp
                j +=1 #for the next step
                mi = (j)*(k//self.size_array)

        self.z[-1] = ztemp
        self.U[-1] = Utemp
        self.UW[-1] = Uwtemp
        self.UI = np.absolute(self.U)**2
        self.UIW = np.absolute(self.UW)**2
        print('last z evaluated: ',z, 'km')
        print('last z saved: ',self.z[-1], 'km')
        print('\u03B1 = {0}'.format('%.3f' % (self.alpha)))
        print('Leff = {0}'.format('%.3f' % (self.compute_Leff())))
        print('\u03C6 max = {0}*\u03C0'.format('%.3f' % (self.compute_phimax_nogvd()/np.pi)))
        return self.U, self.UI, self.UW, self.UIW, self.W, self.z
    ##-------------------------------------------------------##
    ##-------------------------------------------------------##
    #----- PLOTS -----#
    
    def plot_prop(self):
        propagation = go.Figure(data=[go.Heatmap(
                    x = self.T/self.T_selected,#np.sort(t),
                    y = self.z,#np.sort(z),
                    z = self.UI,
                    #type = 'heatmap',
                    zsmooth = "best",#to avoid line at the end of the plot
                    colorscale = 'Jet')])
        propagation.update_layout(
                                width=600, height=600,
                                yaxis=dict(range=[self.z[0], self.z[-1]],title='z [km]', ), 
                                xaxis=dict(range=[-8, 8],title='T/T0',), 
                                )
        return propagation
