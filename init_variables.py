# coding:utf8
import numpy as np
import matplotlib.pyplot as plt
import dash_core_components as dcc
from scipy.integrate import cumtrapz, solve_ivp
import plotly.graph_objects as go
from numpy.core.numeric import Inf  

#Last edition date:
from datetime import date

today = date.today()
date = today.strftime("%d.%m.%Y")###Edit before deployment with, e.g., '28.06.2021'
#-----------------//

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
    
    tempU = val
    # definition of f = N*U
    #f = 1j*(gamma*P0*T0**2/np.absolute(beta2))*(np.absolute(U)**2)*U
    #f = 1j*(np.exp(-alpha*z)*gamma*P0*T0**2)*(np.absolute(U)**2)*U
    '''
    dU(T,z)/dz = N*U(T,z) = i*exp(-alpha*z)/LNL*(|U(T,z)|²)*U(T,z)
    '''
    f  = 1j*(np.exp(-alpha*z)*gamma*P0)*(np.absolute(tempU)**2)*tempU
    return f
#--------------------CLASS--------------------------------#
class Pulse:
    def __init__(self,  T0, T,  m = 1, C=0, pulsetype = 'Gaussian', wavelength  = 1550E-9):
        self.pulsetype = pulsetype
        #self.solution_method = soltype
        self.wavelength = wavelength
        self.speed = 299792458
        self.w0 = (2*np.pi*self.speed)/(self.wavelength)
        self.T0 = T0
        self.T = T
        self.m = m
        if self.m == 0:
            self.m = 1
            print('m = 0 not allowed, working with 1 instead')
        self.C = C
        self.initParam()
        
        self.colors = { #Colors for the plots and web interface:
            'text': '#111111',
            'background': '#FFFFFF',
            'circle': 'whitesmoke',
            'even': 'darkgreen',
            'odd': 'darkblue',
            'other': 'darkred',
            }


    def initParam(self):
        self.dT = self.T[1]-self.T[0] # step_size
        self.points = len(self.T)
        #Normalized Slow-Varying Envelope
        if self.pulsetype == 'Gaussian':
            self.UT0 = (np.exp(-((1+1j*self.C)/2)*((self.T/self.T0)**(2*self.m)))).astype(complex) #dtype = 'complex' in order to have complex values on solve_ivp 
            # self.T_FWHM = 2*np.sqrt(np.log(2))*self.T0
            # self.T0 = self.T_FWHM
        elif self.pulsetype == 'Sech':
            self.UT0 = (1/(np.cosh(self.T/self.T0))*np.exp(-(1j*self.C*self.T**2)/(2*self.T0**2))).astype(complex)
            # self.T_FWHM = 2*np.log(1+np.sqrt(2))*self.T0
            # self.T0 = self.T_FWHM
        else:
            raise ValueError("Pulse '{0}' not found, it must be 'Sech' or 'Gaussian'.".format(self.pulsetype))
        self.UW0 =np.fft.fftshift(np.fft.ifft(self.UT0))
        self.W = 2*np.pi* np.fft.fftfreq(self.points, self.dT) #+ self.w0
        #2*pi*(-N/2:N/2-1)/(N*dT)
        self.W = np.fft.fftshift(self.W) 

#------------------------------------------------------------#


#------------------------CLASS-------------------------------#
class Propagation(Pulse):
    def __init__(self,T0, T, solve_type= 'incident_field', L=0.1, beta2=0, gamma=0, P0=0,  beta3=0, loss = 0, pulsetype = 'Gaussian', m = 1, C=0, z0=0,h_step = 0.004, size_array = 51, wavelength=1550E-9):
        Pulse.__init__(self, T0, T, m=m, C=C, pulsetype=pulsetype)
        self.wavelength = wavelength
        self.speed = 299792458
        self.w0 = (2*np.pi*self.speed)/(self.wavelength)
        self.zmax = L #[m]  #Max distance that will be evaluated
        self.beta2 = beta2 #[ps²/m]
        self.gamma = gamma #[1/(W m)] !!!!!! 
        self.P0 = P0 #[W]
        self.beta3= beta3 #[ps³/m]
        self.loss = loss*1E-3 #[dB/m]
        self.z0 = z0 #At this point will be evaluated the pulse while using the 'only_gvd' mode
        self.h_step = h_step #step-size of the SSFM
        self.size_array = size_array #Size of the array where we are going to save the values of U(T) and U(w)
        self.U = np.zeros((self.size_array,self.points), dtype=complex) # U(z,T) is divided by certain ammount of steps and it's time length is N for each one of the steps
        self.UW = np.zeros((self.size_array,self.points), dtype=complex)
        self.solve_type = solve_type

        if self.solve_type == 'only_gvd':
            self.incident_field()
        elif self.solve_type == 'gauss_gvd':
            self.Gaussian_pulse_GVD()
        elif self.solve_type == 'only_spm':
            self.incident_field_spm()
        elif self.solve_type == 'split_step':
            if self.beta2  == 0 or self.gamma  == 0 or self.P0 == 0:
                print('Either beta2, gamma or P0 are zero, this could lead to problems with the calculation of the current method selected!')
            if self.zmax  == 0:
                raise ValueError('Fiber Length cannot be 0 for de desired solution')       
            self.split_step()
        else:
            raise ValueError("Solution method '{0}' does not exist, it must be 'only_gvd', 'gauss_gvd','only_spm' or 'split_step'.".format(self.solve_type))

    #Lengths
    def compute_LD(self):#Dispersion Length:
        if self.beta2 != 0:
            return ((self.T0)**2)/np.absolute(self.beta2)
        else:
            return Inf 

    def compute_LNL(self):#Length for nonlinearities:
        if self.gamma != 0 and self.P0 != 0 :
            return 1/(self.gamma*self.P0)
        else:
            return Inf 

    def compute_Leff(self): #Effective Length:
        if self.alpha !=0:
            return (1 - np.exp(-self.alpha*self.zmax))/self.alpha
        else: 
            return self.zmax 

    def compute_phimax_nogvd(self):
        Leff = self.compute_Leff()
        return self.gamma*self.P0*Leff #####BE CAREFUL!!!!!!!!!!! 0*inf...

    ##------------------------For GVD------------------------##
    def incident_field(self):  #Function using Eq. 3.2.5 and 3.2.6
        self.UW = self.UW0*np.exp((1j*self.beta2*self.W**2*self.z0/2))
        self.UT = np.fft.fft(np.fft.fftshift(self.UW))
        self.UI = np.absolute(self.UT)**2
        self.UIW = np.absolute(self.UW)**2
        self.UIW = self.UIW/np.amax(self.UIW) 
        return self.UT, self.UI, self.UW, self.UIW

    def Gaussian_pulse_GVD(self): #Function using Eq. 3.2.7 and 3.2.9 Agrawal
        self.UT = self.T0/(np.sqrt(self.T0**2 - 1j*self.beta2*self.z0))*np.exp(-self.T**2/(2*(self.T0**2 - 1j*self.beta2*self.z0)))
        self.UI = np.absolute(self.UT)**2
        self.UW = np.fft.fftshift(np.fft.ifft(self.UT))
        self.UIW = np.absolute(self.UW)**2
        self.UIW = self.UIW/np.amax(self.UIW) 
        return self.UT, self.UI, self.UW, self.UIW


    def Sech_pulse_GVD(self):#Not used 
        self.U0 = 1/(np.cosh(self.T/self.T0))*np.exp(-(1j*self.C*self.T**2)/(2*self.T0**2))
        self.UT = self.T0/(np.sqrt(self.T0**2 - 1j*self.beta2*self.z0))*np.exp(-self.T**2/(2*(self.T0**2 - 1j*self.beta2*self.z0)))
        self.UI = np.absolute(self.UT)**2
        self.UW = np.fft.fftshift(np.fft.ifft(self.UT))
        self.UIW = np.absolute(self.UW)**2
        self.UIW = self.UIW/np.amax(self.UIW) 
        return self.UT, self.UI, self.UW, self.UIW
    ##-------------------------------------------------------##
    ##-------------------------------------------------------##


    ##------------------------For SPM------------------------##
    def incident_field_spm(self):
        #Leff = self.compute_Leff()
        LNL = self.compute_LNL()
        Leff = LNL # Normalized to get the plots from the book
        if self.pulsetype == 'Gaussian':
            if self.m <= 1: #This conditionals were written when the UT0 was not already generalized for chirped-super-gaussian pulsetypes. Following also the book Arawal
                self.Phi_NL = (Leff/LNL)*np.absolute(self.UT0)**2
                self.delta_w = delta_g(self.T, self.T0, Leff, LNL)
            else: 
                self.delta_w = delta_g(self.T, self.T0, Leff, LNL, self.m)
                #for array-like data
                #Phi_NL = -cumtrapz(delta_w, T, initial=0) #Default 'initial' is None, which means no value at x[0] 
                #for function-like input:
                self.Phi_NL = -mid_step(0, delta_g, self.T, self.T0, Leff, LNL, self.m)

        elif self.pulsetype == 'Sech':
            phinl = lambda t: (Leff/LNL)*np.absolute(1/(np.cosh(t/self.T0))*np.exp(-(1j*self.C*t**2)/(2*self.T0**2)))**2
            self.Phi_NL = phinl(self.T)
            self.delta_w = -derivative(phinl, self.T, dx=self.dT)
        else:
            raise ValueError("Pulse '{0}' not found, it must be 'Sech' or 'Gaussian'.".format(self.pulsetype))

        return self.Phi_NL, self.delta_w
    ##-------------------------------------------------------##
    ##-------------------------------------------------------##

    ##---------------------SPLIT-STEP------------------------##
    def split_step(self):
        self.alpha = np.log(10**(self.loss/10)) #log y = ln y / ln 10 
        h = self.zmax*self.h_step #step size m -> depends on how long is the fiber we are working with
        print('\nZ max.: ',self.zmax, 'm')
        print('h size: ',h, 'm')
        print('LD: ',self.compute_LD(), 'm')
        print('LNL: ',self.compute_LNL(), 'm')
        print('N = sqrt(LD/LNL) = : ',np.sqrt(self.compute_LD()/self.compute_LNL()))
        print(1/self.h_step,'number of points evaluated')
        #initial values for split-step
        self.U[0] = self.UT0  #first value for U(0,T)  
        self.UW[0] = self.UW0 #first value for U(0,W)  
        #for general def (Eq. (2.4.3)): 
        Dw = -(1j*self.beta2/2)*(1j*self.W)**2 + (self.beta3/6)*(1j*self.W)**3 - self.alpha/2 #D = -1j*beta2/2*(δ²/T²) + beta3/6*(δ³/T³) - alpha/2
        #TODO: create a way to extract the betas from an array and take them to the frequency domain
        self.z = np.zeros((self.size_array))#array for saving the nth-point U(T,zn) and U(w,zn)
        mi = int(1/(self.size_array*self.h_step)) #just a counter 
        j = 1
        Utemp = self.UT0[:] #variable to update

        #for i in np.linspace(1,zmax*1000-h,k):  if we pre-defined the resolution (or ammount of steps k) independent of the step size h
        for k in range(0,int(1/self.h_step),1): # from z to z+h k times until z_max is reached
            ztemp = k*h #  [m] actual z to evaluate
            Uwtemp = np.fft.fftshift(np.fft.ifft(Utemp))
            Utemp = np.fft.fft(np.fft.fftshift(np.exp((h/2)*Dw)*Uwtemp))
            sol = solve_ivp (U_noGVD, # function to integrate
                        (ztemp, ztemp+h) , # duration [ z_init - z_end ]
                        Utemp, # initial conditions
                        method ='RK23', # method of integration 23
                        #t_eval =np.linspace(z, z+h, size_array), # points where the sol. is stored
                        rtol =1e-8 , # relative accuracy
                        args =[self.gamma, self.P0, self.T0, self.alpha, self.beta2, self.beta3]) # arguments for the function
            Utemp = sol.y.T[-1] #value at z+h
            Uwtemp = np.fft.fftshift(np.fft.ifft(Utemp))
            Utemp = np.fft.fft(np.fft.fftshift(np.exp((h/2)*Dw)*Uwtemp))
            Uwtemp = np.fft.fftshift(np.fft.ifft(Utemp))
            if k == mi and j< self.size_array:
                self.z[j] = ztemp
                self.U[j] = Utemp
                self.UW[j] = Uwtemp
                j +=1 #for the next step
                mi += int(1/(self.size_array*self.h_step))
        self.z[-1] = ztemp
        self.U[-1] = Utemp
        self.UW[-1] = Uwtemp
        self.UI = np.absolute(self.U)**2
        self.UIW = np.absolute(self.UW)**2
        #self.UIW = self.UIW/np.amax(self.UIW) 
        self.UIW = self.UIW/np.amax(self.UIW[0])
        print('last z evaluated: ',ztemp, 'm')
        print('last z saved: ',self.z[-1], 'm')
        print('\u03B1 = {0}'.format('%.3f' % (self.alpha)))
        print('Leff = {0}'.format('%.3f' % (self.compute_Leff())))
        print('\u03C6 max = {0}*\u03C0'.format('%.3f' % (self.compute_phimax_nogvd()/np.pi)))
        #self.W += self.w0  #Change axis of the functions plot_propagation() and plot_envelope()
        return self.U, self.UI, self.UW, self.UIW, self.W, self.z

    ##-------------------------------------------------------##
    #----- PLOTS -----#
    def plot_propagation(self, mode = 'time', xrange = [-8, 8], size_img = 600):
        if mode == 'time':
            intensity = self.UI
            x = self.T/self.T0#np.sort(t),
            x_title = 'T/T0'
            #updatemenus = []
        elif mode == 'spectrum':
            intensity = self.UIW
            x = self.W
            x_title = '\u03C9-\u03C90'
            xrange = [-0.1E14, 0.1E14]
            # xrange = [-0.1E14+self.w0, 0.1E14+self.w0]
            # updatemenus = [
            #                 dict(
            #                     type="buttons",
            #                     direction="left",
            #                     buttons=list([
            #                         dict(
            #                             args=[{'xaxis.type': 'linear'}],
            #                             label="Linear Scale",
            #                             method="relayout"
            #                         ),
            #                         dict(
            #                             args=[{'xaxis.type': 'log'}],
            #                             label="Log Scale",
            #                             method="relayout"
            #                         )
            #                     ])
            #                 ),
            #             ]
        else:
            raise ValueError("Mode '{0}' not found. Modes available 'time' or 'spectrum'.".format(mode) )
        propagation = go.Figure(data=[go.Heatmap(
                    x = x,
                    y = self.z,#np.sort(z),
                    z = intensity,
                    #type = 'heatmap',
                    zsmooth = "best",#to avoid line at the end of the plot
                    colorscale = 'Jet')])
        propagation.update_layout(
                                width=size_img, height=size_img,
                                yaxis=dict(range=[self.z[0], self.z[-1]],title='z [m]', ), 
                                #xaxis=dict(title=x_title,), 
                                xaxis=dict(range=xrange,title=x_title,), 
                                )
        #propagation.update_layout(updatemenus=updatemenus)
        return propagation
#------------------------------------------------------------#
    def plot_envelope(self, mode = 'time', z0 = 0, xrange = [-8, 8],  yrange = [0, 1.1], size_img = 600):
        if mode == 'time':
            intensity = self.UI
            x = self.T/self.T0
            x_title = 'T/T0'
            y_title = '|U(z,T)|^2'
            #updatemenus = []
        elif mode == 'spectrum':
            intensity = self.UIW
            x = self.W
            x_title = '\u03C9-\u03C90'
            y_title = '|U(z,\u03C9)|^2'
            xrange = [-0.1E14, 0.1E14]
            # xrange = [-0.1E14+self.w0, 0.1E14+self.w0]
            # updatemenus = [
            #                 dict(
            #                     type="buttons",
            #                     direction="left",
            #                     buttons=list([
            #                         dict(
            #                             args=[{'xaxis.type': 'linear'}],
            #                             label="Linear Scale",
            #                             method="relayout"
            #                         ),
            #                         dict(
            #                             args=[{'xaxis.type': 'log'}],
            #                             label="Log Scale",
            #                             method="relayout"
            #                         )
            #                     ])
            #                 ),
            #             ]
        else:
            raise ValueError("Mode '{0}' not found. Modes available 'time' or 'spectrum'.".format(mode) )
        

        scatter = go.Scatter(x=x,y=intensity[z0], name = self.pulsetype,
                                line=dict(color=self.colors['even']))
        
        figure = go.Figure(data=[scatter]).update_layout(
                                width=size_img, height=size_img,
                                plot_bgcolor  = self.colors['background'],
                                paper_bgcolor = self.colors['background'],
                                font= {'color': self.colors['text']},
                                yaxis=dict(range=yrange,title=y_title, ), 
                                #axis=dict(title=x_title,), 
                                xaxis=dict(range=xrange,title=x_title, ), 
                                title = '{0} pulse.'.format(self.pulsetype),
                                )
        #figure.update_layout(updatemenus=updatemenus)
        return figure
        # To plot and prepare the go.Figure for use inside the call back, the object Pulse should be called:

        # env_fig = pulse.plot_envelope()
        # env_graph = dcc.Graph(id='desired_id', #id for callback purposes
        #                         animate=True,
        #                         figure=env_fig.update_layout(
        # ))        

    def plot_envelope_GVD(self, mode = 'time', xrange = [-8, 8],  yrange = [0, 1.1], size_img = 600, color = 'darkgreen' ):
        if mode == 'time':
            intensity = self.UI
            x = self.T/self.T0
            x_title = 'T/T0'
            y_title = '|U(z,T)|^2'
        elif mode == 'spectrum':
            intensity = self.UIW
            x = self.W
            x_title = '\u03C9-\u03C90'
            y_title = '|U(z,\u03C9)|^2'
            xrange = [-0.1E14, 0.1E14]
        else:
            raise ValueError("Mode '{0}' not found. Modes available 'time' or 'spectrum'.".format(mode) )
        

        scatter = go.Scatter(x=x,y=intensity, name = self.pulsetype,
                                line=dict(color=color))
        
        return go.Figure(data=[scatter]).update_layout(
                                width=size_img, height=size_img,
                                plot_bgcolor  = self.colors['background'],
                                paper_bgcolor = self.colors['background'],
                                font= {'color': self.colors['text']},
                                yaxis=dict(range=yrange,title=y_title, ), 
                                #xaxis=dict(title=x_title,), 
                                xaxis=dict(range=xrange,title=x_title, ), 
                                #title = '{0} pulse.'.format(self.pulsetype),
                                )
def plot_shift(pulse1, pulse2, pulse3, xrange = [-2.5, 2.5], yrange = [0, 1.1], size_img = 600):
    p1 = go.Scatter(x=pulse1.T/pulse1.T0 ,y=pulse1.Phi_NL, name = 'Gaussian',
                         line=dict(color='darkgreen'))

    p2 = go.Scatter(x=pulse2.T/pulse2.T0 ,y=pulse2.Phi_NL, name = 'Super-Gaussian',
                        line=dict(color='darkblue'))

    p3 = go.Scatter(x=pulse3.T/pulse3.T0 ,y=pulse3.Phi_NL, name = 'Sech',
                         line=dict(color='darkred'))

    spm_phase = go.Figure(data=[p1,p2,p3]).update_layout( 
        updatemenus = list([
            dict(
                type="buttons",
                active=0,
                buttons=list([   
                    dict(label = 'Gaussian',
                        method = 'update',
                        args = [{'visible': [True, True, False]},
                                {'title': 'Phase: Gaussian Pulse.'}]), 

                    dict(label = 'Sech',
                        method = 'update',
                        args = [{'visible': [False, False, True]},
                                {'title': '''
                                Phase: Sech Pulse.'''}])  
                ]),
            )
        ])
    )


    spm_phase.update_layout( 
                            width=size_img, height=size_img,
                            plot_bgcolor  = '#FFFFFF',
                            paper_bgcolor = '#FFFFFF',
                            font= {
                                    'color': '#111111'},
                            yaxis=dict(range=yrange,title='Phase \u03C6NL', 
                                        ), 
                            xaxis=dict(range=xrange,title='T/T0', 
                                        ), 
                            )
    return spm_phase


def plot_chirp(pulse1, pulse2, pulse3, xrange = [-2.5, 2.5], yrange = [-3, 3], size_img = 600):

    chirp1 = go.Scatter(x=pulse1.T/pulse1.T0 ,y=pulse1.delta_w*pulse1.T0, name = 'Gaussian',
                            line=dict(color='darkgreen'))
    
    chirp2 = go.Scatter(x=pulse2.T/pulse2.T0 ,y=pulse2.delta_w*pulse2.T0, name = 'Super-Gaussian',
                            line=dict(color='darkblue'))

    chirp3 = go.Scatter(x=pulse3.T/pulse3.T0 ,y=pulse3.delta_w*pulse3.T0, name = 'Sech',
                            line=dict(color='darkred'))


    spm_chirp = go.Figure(data=[chirp1,chirp2, chirp3]).update_layout( 
        updatemenus = list([
            dict(
                type="buttons",
                active=0,
                buttons=list([   
                    dict(label = 'Gaussian',
                        method = 'update',
                        args = [{'visible': [True, True, False]},
                                {'title': 'Chirp: Gaussian Pulse.'}]), 

                    dict(label = 'Sech',
                        method = 'update',
                        args = [{'visible': [False, False, True]},
                                {'title': '''
                                Chirp: Sech Pulse.'''}])  
                ]),
            )
        ])
    )


    spm_chirp.update_layout( 
                            width=size_img, height=size_img,
                            plot_bgcolor  = '#FFFFFF',
                            paper_bgcolor = '#FFFFFF',
                            font= {
                                    'color': '#111111'},
                            yaxis=dict(range=yrange,title='frequency chirp \u03B4\u03C9T0', 
                                        ), 
                            xaxis=dict(range=xrange,title='T/T0', 
                                        ), 
                            )
    return spm_chirp

def plot_envelope_x2(pulse1, pulse2, mode = 'time', z0 = 0, xrange = [-8, 8], yrange = [0, 1.1], size_img = 600):
    if mode == 'time':
        intensity1 = pulse1.UI
        intensity2 = pulse2.UI
        x1 = pulse1.T/pulse1.T0
        x2 = pulse2.T/pulse2.T0
        x_title = 'T/T0'
        y_title = '|U(z,T)|^2'
    elif mode == 'spectrum':
        intensity1 = pulse1.UIW
        intensity2 = pulse2.UIW
        x1 = pulse1.W
        x2 = pulse2.W
        x_title = '\u03C9-\u03C90'
        y_title = '|U(z,\u03C9)|^2'
    else:
        raise ValueError("Mode '{0}' not found. Modes available 'time' or 'spectrum'.".format(mode) )
    

    scatter1 = go.Scatter(x=x1,y=intensity1[z0], name = pulse1.pulsetype,
                            line=dict(color=pulse1.colors['even']))
    scatter2 = go.Scatter(x=x2,y=intensity2[z0], name = pulse1.pulsetype,
                            line=dict(color=pulse2.colors['odd']))
    
    env_fig = go.Figure(data=[scatter1, scatter2]).update_layout(    
        updatemenus = list([dict(type="buttons",
                                active=0,
                                buttons=list([dict(label = pulse1.pulsetype,
                                                    method = 'update',
                                                    args = [{'visible': [True, False]},
                                                            {'title': '{0} pulse.'.format(pulse1.pulsetype)}]), 
                                            dict(label = pulse2.pulsetype,
                                                    method = 'update',
                                                    args = [{'visible': [False, True]},
                                                            {'title': '{0} pulse.'.format(pulse2.pulsetype)}])
                            ]),
                        )
                    ])
                )
    env_fig.update_layout( 
                            width=size_img, height=size_img,
                            plot_bgcolor  = pulse1.colors['background'],
                            paper_bgcolor = pulse1.colors['background'],
                            font= {'color': pulse1.colors['text']},
                            yaxis=dict(range=yrange,title=y_title, ), 
                            xaxis=dict(range=xrange,title=x_title, ), 
                            )
    # To plot and prepare the go.Figure for use inside the call back:
    # env_fig = plot_envelope_x2(pulse1, pulse2)
    # env_graph = dcc.Graph(id='desired_id', #id for callback purposes
    #                         animate=True,
    #                         figure=env_fig.update_layout(

    # ))        
    return env_fig

def multiplots_1D( *pulses,  mode = 'time', xrange = [-8, 8],  yrange = [0, 1.1], size_img = 600 ):
    scatters = []
    
    if mode == 'time':
        x_title = 'T/T0'
        y_title = '|U(z,T)|^2'
        for pulse in pulses:
            intensity = pulse.UI
            x = pulse.T/pulse.T0
            scatter = go.Scatter(x=x,y=intensity, name = pulse.pulsetype,)
                                #line=dict(color=pulse.colors['even']))
            scatters.append(scatter)
        return go.Figure(data=scatters).update_layout(width=size_img, height=size_img,
                                                        plot_bgcolor  = pulses[0].colors['background'],
                                                        paper_bgcolor = pulses[0].colors['background'],
                                                        font= {'color': pulses[0].colors['text']},
                                                        yaxis=dict(range=yrange,title=y_title, ), 
                                                        xaxis=dict(range=xrange,title=x_title, ), 
                                                        title = '{0} pulse.'.format(pulses[0].pulsetype),)
    elif mode == 'spectrum':
        x_title = '\u03C9-\u03C90'
        y_title = '|U(z,\u03C9)|^2'
        xrange = [-2E12, 2E12]
        for pulse in pulses:
            intensity = pulse.UIW
            x = pulse.W
            scatter = go.Scatter(x=x,y=intensity, name = pulse.pulsetype,)
                                #line=dict(color=pulse.colors['even']))
            scatters.append(scatter)
        return go.Figure(data=scatters).update_layout(width=size_img, height=size_img,
                                                        plot_bgcolor  = pulses[0].colors['background'],
                                                        paper_bgcolor = pulses[0].colors['background'],
                                                        font= {'color': pulses[0].colors['text']},
                                                        yaxis=dict(range=yrange,title=y_title, ), 
                                                        xaxis=dict(range=xrange,title=x_title, ), 
                                                        title = '{0} pulse.'.format(pulses[0].pulsetype),)
    else:
        raise ValueError("Mode '{0}' not found. Modes available 'time' or 'spectrum'.".format(mode) )


def multiplots_2D(*pulses, mode = 'time', z0 = 0, xrange = [-8, 8],  yrange = [0, 1.1], size_img = 600 ):
    scatters = []
    
    if mode == 'time':
        x_title = 'T/T0'
        y_title = '|U(z,T)|^2'
        for pulse in pulses:
            intensity = pulse.UI
            x = pulse.T/pulse.T0
            scatter = go.Scatter(x=x,y=intensity[z0], name = pulse.pulsetype,
                                line=dict(color=pulse.colors['even']))
            scatters.append(scatter)
        return go.Figure(data=scatters).update_layout(width=size_img, height=size_img,
                                                        plot_bgcolor  = pulses[0].colors['background'],
                                                        paper_bgcolor = pulses[0].colors['background'],
                                                        font= {'color': pulses[0].colors['text']},
                                                        yaxis=dict(range=yrange,title=y_title, ), 
                                                        xaxis=dict(range=xrange,title=x_title, ), 
                                                        title = '{0} pulse.'.format(pulses[0].pulsetype),)
    elif mode == 'spectrum':
        x_title = '\u03C9-\u03C90'
        y_title = '|U(z,\u03C9)|^2'
        xrange = [-2E12, 2E12]
        for pulse in pulses:
            intensity = pulse.UIW
            x = pulse.W
            scatter = go.Scatter(x=x,y=intensity[z0], name = pulse.pulsetype,
                                line=dict(color=pulse.colors['even']))
            scatters.append(scatter)
        return go.Figure(data=scatters).update_layout(width=size_img, height=size_img,
                                                        plot_bgcolor  = pulses[0].colors['background'],
                                                        paper_bgcolor = pulses[0].colors['background'],
                                                        font= {'color': pulses[0].colors['text']},
                                                        yaxis=dict(range=yrange,title=y_title, ), 
                                                        xaxis=dict(range=xrange,title=x_title, ), 
                                                        title = '{0} pulse.'.format(pulses[0].pulsetype),)
    else:
        raise ValueError("Mode '{0}' not found. Modes available 'time' or 'spectrum'.".format(mode) )


def plot_envelope(x,y, pulsetype,colors,mode = 'time', z0 = 0, xrange = [-8, 8],  yrange = [0, 1.1], size_img = 600):
    if mode == 'time':
        x_title = 'T/T0'
        y_title = '|U(z,T)|^2'
        #updatemenus = []
    elif mode == 'spectrum':
        x_title = '\u03C9-\u03C90'
        y_title = '|U(z,\u03C9)|^2'
        xrange = [-0.1E14, 0.1E14]
        # xrange = [-0.1E14+self.w0, 0.1E14+self.w0]
        # updatemenus = [
        #                 dict(
        #                     type="buttons",
        #                     direction="left",
        #                     buttons=list([
        #                         dict(
        #                             args=[{'xaxis.type': 'linear'}],
        #                             label="Linear Scale",
        #                             method="relayout"
        #                         ),
        #                         dict(
        #                             args=[{'xaxis.type': 'log'}],
        #                             label="Log Scale",
        #                             method="relayout"
        #                         )
        #                     ])
        #                 ),
        #             ]
    else:
        raise ValueError("Mode '{0}' not found. Modes available 'time' or 'spectrum'.".format(mode) )
    

    scatter = go.Scatter(x=x,y=y[z0], name = pulsetype,
                            line=dict(color=colors['even']))
    
    figure = go.Figure(data=[scatter]).update_layout(
                            width=size_img, height=size_img,
                            plot_bgcolor  = colors['background'],
                            paper_bgcolor = colors['background'],
                            font= {'color': colors['text']},
                            yaxis=dict(range=yrange,title=y_title, ), 
                            #axis=dict(title=x_title,), 
                            xaxis=dict(range=xrange,title=x_title, ), 
                            title = '{0} pulse.'.format(pulsetype),
                            )
    #figure.update_layout(updatemenus=updatemenus)
    return figure
