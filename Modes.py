#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 10:56:02 2020

@author: njoly

"""

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import cm
from scipy.special import kvp, kv, jv, jn_zeros
from sys import exit
from scipy.constants import mu_0, epsilon_0, speed_of_light
from scipy.interpolate import UnivariateSpline
# from scipy.interpolate import interp1d
import dash_core_components as dcc
import plotly.graph_objects as go   
from plotly.tools import mpl_to_plotly
#import plotly.express as px

class Material:
    """ Class constructing the material:
        - refractive index """
    #order must be 
    def __init__(self, wavelength=0.800, mat="air"):
        self.material=mat
        self.wavelength=wavelength
        self.A = []
        self.B = []
        
        if self.material=="silica": self.initSilica()
        if self.material=="GeO2_63": self.initGeO2_6_3()
        if self.material=="GeO2_193": self.initGeO2_19_3()
        if self.material=="B2O3": self.initB2O3()
        if self.material=="P2O5": self.initP2O5()
        if self.material=="air": self.initAir()
        

    def initGeO2_19_3(self):
        self.A = [0.77347008, 0.4461191, 0.8081698]
        self.B = [0.07646793, 0.12460807, 9.89620331]
        self.order= 0


    def initGeO2_6_3(self):
        self.A = [0.7083952, 0.4203993, 0.8663412]
        self.B = [0.08538421, 0.10248385, 9.89617502]
        self.order= 1

    def initP2O5(self):
        self.A = [0.7058489, 0.4176021, 0.8952753]
        self.B = [0.07212788, 0.11347819, 9.89616138]
        self.order= 2
        
    def initSilica(self):
        self.A = [0.6961663, 0.4079426, 0.8974794]
        self.B = [0.068404, 0.116241, 9.896161] 
        self.order= 3

    def initB2O3(self):
        self.A = [0.6910021, 0.4022430, 0.09439644]
        self.B = [0.22320031, 0.1172887 , 9.89613713]
        self.order= 4
        
    def initAir(self):
        self.A = [0, 0, 0]
        self.B = [0, 0, 0]
        self.order= 5


    def refIndex(self, wavelength):
        """Returns the refractive index for a given wavelength
        [wavelength]=microns"""
        n2=1
        for i in range(3):
            n2 += self.A[i]/(1-(self.B[i]/wavelength)**2)
        return np.sqrt(n2)
    
    def show(self):
        """Returns the information about the material in use
        - current wavelength
        - refractive index
        - name of the material in use """
        n=self.refIndex(self.wavelength)
        print("material is {0}".format(self.material))
        print("Its refractive index at {0} is {1}".format(self.wavelength, n))

class Mode:
    """Class defining the necessary informations requires for preparing
    the spatial mode:
        - core and cladd refracive indices
        - operating wavelength
        - index of the mode (nu in Bessel function) 
        Units: ???"""
        
    def __init__(self, diameter, wavelength, matCore="silica", matCladd="air", nu=1):        
        self.matCore = Material(mat = matCore)#Material(mat="silica")
        self.matCladd = Material(mat = matCladd)
        self._radius=diameter/2.
        self._wavelength=wavelength
        self.nu=nu
        self.initParam()

    def initParam(self):
        self.nCore = self.matCore.refIndex(self._wavelength)
        self.nCladd= self.matCladd.refIndex(self._wavelength)
        self.k0 = 2*np.pi/self._wavelength
        self.V = self.k0*self.radius*np.sqrt(self.nCore**2 - self.nCladd**2)
        self.delta=0.5*(1-(self.nCladd/self.nCore)**2)



    def ploteigen(self, Umin, Umax, N, colors):
        eps=1E-5
        if Umax>self.V: Umax=self.V-eps
        U=np.linspace(Umin,Umax,N)
        #print("nu=",self.nu)
        #print('This U: ',U)
        #print(self.eigenValuePb(U))
        Plot = go.Scatter(x=U ,y=self.eigenValuePb(U), name = 'EigenValuePb',
                            line=dict(color=colors['even']))
        #plot_title = 'Eigenvalue Pb for {0} mode'.format(self.typeMode)
        return Plot#, plot_title

    def plotEigenValuePb(self,Umin,Umax,N, color='b'):
        eps=1E-5
        if Umax>self.V: Umax=self.V-eps
        U=np.linspace(Umin,Umax,N)
        print("nu=",self.nu)
        plt.figure()
        plt.title("Eigenvalue Pb for {0} mode".format(self.typeMode))
        plt.plot(U,self.eigenValuePb(U), color=color)
        plt.axhline(y=0, color='k', linestyle=':')
        plt.ylim(-5,5)
        plt.show()       

    def solveEigenValuePb(self, m):
        self.m = m
        Uinit=self.guessUnit(m)
        self.U = self.solve(Uinit)
        self.W = np.sqrt(self.V * self.V - self.U * self.U)
        self.beta = np.sqrt((self.k0*self.nCore)**2 - (self.U/self.radius)**2) #sqrt((ko*n)^2-alpha^2)
        return self.U
    
    def computeBeta(self,m):
        tmp=self.solveEigenValuePb(m)
        return self.beta
    
    def solve(self, Uinit, nIterMax=10, prec=1E-10):
        error=2.*prec
        U1=Uinit
        U2=U1+error
        y2=self.eigenValuePb(U2)
        k=0
        while ((k<=nIterMax) and (error>prec)):
            y1=self.eigenValuePb(U1)          
            old_y1 = y1
            df = (y2-y1)/(U2-U1)
            y1 = U1-y1/df
            U2=U1
            y2=old_y1
            error=np.abs(U1-y1)
            U1=y1
            k+=1
        return U1
        
    def _getWavelength(self):
        return self._wavelength
    def _setWavelength(self, wavelength):
        self._wavelength=wavelength
        self.initParam() 
    def _getRadius(self):
        return self._radius
    def _setRadius(self, radius):
        self._radius= radius
        self.initParam()

    # def _getpam(self, wavelength, radius):
    #     self._wavelength=wavelength
    #     self._radius = radius
    #     self.initParam() 
    #     return self._radius, self._wavelength

    def show(self):
        print("radius = {0} um".format(self.radius))
        print("lambda = {0} um".format(self.wavelength))
        print("nCl = {0} ; nCo = {1}".format(self.nCladd, self.nCore))
        print("-> delta = {0}".format(self.delta))
        print("V = {0}".format(self.V))
        print("Nbre of available modes: {0}".format(self.nbreOfModes()))
        
    def drawCut(self,xMin=-3,xMax=3, nPts=301):
        X=np.linspace(xMin,xMax,nPts)
        plt.figure()
        plt.plot(X,self.rFieldCut(X))
        plt.show()

    def drawCutH(self,xMin=-3,xMax=3, nPts=301):
        X=np.linspace(xMin,xMax,nPts)
        plt.figure()
        plt.plot(X,self.hFieldCut(X))
        plt.show()        

    def drawCut_dash(self, xMin=-3,xMax=3, nPts=301, colors = {'even':'darkgreen'} ):
        X=np.linspace(xMin,xMax,nPts)
        drawfig = go.Scatter(x=X,y=self.rFieldCut(X), name = 'Cut',
                            line=dict(color=colors['even']))
        return drawfig

    def drawCutH_dash(self, xMin=-3,xMax=3, nPts=301, colors = {'odd':'darkblue'} ):
        X=np.linspace(xMin,xMax,nPts)
        drawfig_H = go.Scatter(x=X,y=self.hFieldCut(X), name = 'Cut_H',
                            line=dict(color=colors['odd']))
        return drawfig_H

    def rFieldCut(self,rho):
        rho/=self.radius
        return np.where(np.abs(rho)<=1, 
                        self.e_rCore(np.abs(rho),0), 
                        self.e_rCladding(np.abs(rho),0))
    def hFieldCut(self,rho):
        rho/=self.radius
        return np.where(np.abs(rho)<=1, 
                        self.h_rCore(np.abs(rho),np.pi/2.), 
                        self.h_rCladding(np.abs(rho),np.pi/2.))
    
    def plotMode(self,xmin=-5,xmax=5, ymin=-5, ymax=5, npts=256):
        # en principe carre => Seult X max est necessaire et x à initialiser 
        # pour faire le mesh
        x=np.linspace(xmin,xmax,npts)
        y=np.linspace(ymin,ymax,npts)
        X,Y=np.meshgrid(x,y)
        Z=np.abs(self.Sz(X,Y))
        plt.figure()
        if (self.typeMode=='EH')or (self.typeMode=="HE"):
            plt.title(r"$|S_z|^2$ for {0}$_{1}$$_{2}$ mode".format(self.typeMode,self.nu, self.m))
        else:
            plt.title(r"$|S_z|^2$ for {0}$_0$$_{1}$ mode".format(self.typeMode,self.m))
        plt.xlabel(r"r/radius")
        plt.ylabel(r"r/radius")
        x = plt.imshow(Z, extent=[xmin,xmax,ymin,ymax], cmap='thermal', origin='lower', interpolation='none')
        plt.colorbar()
        plt.show()


    def plotMode_dash(self,xmin=-5,xmax=5, ymin=-5, ymax=5, npts=256, ax = 's'):
        # en principe carre => Seult X max est necessaire et x à initialiser 
        # pour faire le mesh
        x=np.linspace(xmin,xmax,npts)
        y=np.linspace(ymin,ymax,npts)
        X,Y=np.meshgrid(x,y)
        if ax == 'ey':
            Z=np.abs(self.yFieldComp(X,Y))
        elif ax == 'ez':
            Z=np.abs(self.zFieldComp(X,Y))
        elif ax == 'ex':
            Z=np.abs(self.xFieldComp(X,Y))
        elif ax == 'hy':
            Z=np.abs(self.yHComp(X,Y))
        elif ax == 'hz':
            Z=np.abs(self.zHComp(X,Y))
        elif ax == 'hx':
            Z=np.abs(self.xHComp(X,Y))
        else:
            Z=np.abs(self.Sz(X,Y))
            
        


        if (self.typeMode=='EH')or (self.typeMode=="HE"):
            title = r"$|S_z|^2$ for {0}$_{1}$$_{2}$ mode".format(self.typeMode,self.nu, self.m)
        else:
            title = r"$|S_z|^2$ for {0}$_0$$_{1}$ mode".format(self.typeMode,self.m)
        xlabel = r"r/radius"
        ylabel = r"r/radius"
        
        figure_mode = go.Figure(data=[go.Heatmap(
                x = np.sort(x),
                y = np.sort(y),
                z = Z,
                type = 'heatmap',
                colorscale = 'thermal')])
        figure_mode.update_layout(        
                shapes = [dict(type="circle",
                                                xref="x",
                                                yref="y",
                                                x0=-self._getRadius(),
                                                y0=-self._getRadius(),
                                                x1=self._getRadius(),
                                                y1=self._getRadius(),
                                                line_color="LightSeaGreen",
                                                line_dash='dash',),],)
        return figure_mode, title, [xlabel, xmin, xmax], [ylabel, ymin, ymax]

    def yFieldComp(self,x,y):
        r = lambda x,y : np.sqrt(x**2 + y**2)
        phi = lambda x,y : np.arctan2(y,x)
        return np.where(np.abs(r(x,y)/self.radius)<=1, 
                        np.abs(self.e_rCore(r(x,y)/self.radius,phi(x,y))*np.sin(phi(x,y))
                        +self.e_thetaCore(r(x,y)/self.radius,phi(x,y))*np.cos(phi(x,y))), 
                        np.abs(self.e_rCladding(r(x,y)/self.radius, phi(x,y))*np.sin(phi(x,y))
                        +self.e_thetaCladding(r(x,y)/self.radius, phi(x,y))*np.cos(phi(x,y)))
                        )
                    
    def xFieldComp(self,x,y):
        r = lambda x,y : np.sqrt(x**2 + y**2)
        phi = lambda x,y : np.arctan2(y,x)
        return np.where(np.abs(r(x,y)/self.radius)<=1, 
                        np.abs(self.e_rCore(r(x,y)/self.radius,phi(x,y))*np.cos(phi(x,y))
                        -self.e_thetaCore(r(x,y)/self.radius,phi(x,y))*np.sin(phi(x,y))), 
                        np.abs(self.e_rCladding(r(x,y)/self.radius, phi(x,y))*np.cos(phi(x,y))
                        -self.e_thetaCladding(r(x,y)/self.radius, phi(x,y))*np.sin(phi(x,y)))
                        )
    

    def zFieldComp(self,x,y):
        r = lambda x,y : np.sqrt(x**2 + y**2)
        phi = lambda x,y : np.arctan2(y,x)
        return np.where(np.abs(r(x,y)/self.radius)<=1, 
                        self.e_zCore(r(x,y)/self.radius,phi(x,y)), 
                        self.e_zCladding(r(x,y)/self.radius, phi(x,y))
                        )   

    def xHComp(self,x,y):
        r = lambda x,y : np.sqrt(x**2 + y**2)
        phi = lambda x,y : np.arctan2(y,x)
        return np.where(np.abs(r(x,y)/self.radius)<=1, 
                        self.h_rCore(r(x,y)/self.radius,phi(x,y))*np.cos(phi(x,y))
                        -self.h_thetaCore(r(x,y)/self.radius,phi(x,y))*np.sin(phi(x,y)), 
                       self.h_rCladding(r(x,y)/self.radius, phi(x,y))*np.cos(phi(x,y))
                       -self.h_thetaCladding(r(x,y)/self.radius, phi(x,y))*np.sin(phi(x,y))
                        )    
    
    def yHComp(self,x,y):
        r = lambda x,y : np.sqrt(x**2 + y**2)
        phi = lambda x,y : np.arctan2(y,x)
        return np.where(np.abs(r(x,y)/self.radius)<=1, 
                        self.h_rCore(r(x,y)/self.radius,phi(x,y))*np.sin(phi(x,y))
                        +self.h_thetaCore(r(x,y)/self.radius,phi(x,y))*np.cos(phi(x,y)), 
                       self.h_rCladding(r(x,y)/self.radius, phi(x,y))*np.sin(phi(x,y))
                       +self.h_thetaCladding(r(x,y)/self.radius, phi(x,y))*np.cos(phi(x,y))
                        )
    
    def zHComp(self,x,y):
        r = lambda x,y : np.sqrt(x**2 + y**2)
        phi = lambda x,y : np.arctan2(y,x)
        return np.where(np.abs(r(x,y)/self.radius)<=1, 
                        self.h_zCore(r(x,y)/self.radius,phi(x,y)), 
                        self.h_zCladding(r(x,y)/self.radius, phi(x,y))
                        ) 
    def Sz(self,x,y):
        r = lambda x,y : np.sqrt(x**2 + y**2)
        phi = lambda x,y : np.arctan2(y,x)
        return np.where(np.abs(r(x,y)/self.radius)<=1, 
                        self.SzCore(r(x,y)/self.radius,phi(x,y)), 
                        self.SzCladding(r(x,y)/self.radius, phi(x,y))
                        ) 

    wavelength=property(_getWavelength, _setWavelength)
    radius=property(_getRadius, _setRadius)

class TEMode(Mode):
    """Class for TE modes"""
      
    #radius=0
    def __init__(self, wavelength, diameter, matCore="silica", matCladd="air", nu=0):
        Mode.__init__(self, diameter, wavelength, matCore, matCladd, nu)
        self.typeMode="TE"
                 
    def eigenValuePb(self, U):
        w = np.sqrt(self.V**2 - np.power(U,2))
        y = jv(1,U)/(U*jv(0,U)) + kv(1,w)/(w*kv(0,w))
        #y[y > 7] = None
        #y[y < -7] = None
        return y
    
    def nbreOfModes(self):
        n=1
        while jn_zeros(self.nu, n)[-1]<self.V:
            n+=1
        return n-1
    
    def guessUnit(self,m):
        mMax = self.nbreOfModes()
        if m>mMax:
            print("Mode TE0{0} at lambda = {1} does not exist".format(m, self.wavelength))
            print("lambda is above the cutoff")
            print("the number of possible modes is {0}".format(mMax))
            exit(0)
        lowBoundary = jn_zeros(0,m)[-1]
        highBoundary = jn_zeros(0,m+1)[-1]
        if highBoundary>self.V: highBoundary=self.V
        return 1./2 * (lowBoundary+highBoundary)
        
    def show(self):
        print("This is a TE mode")
        Mode.show(self)
        
    def e_rCore(self, rho, theta): return 0
    def e_thetaCore(self, rho, theta):
        return -jv(1, self.U*rho)/jv(1,self.U)
    def e_zCore(self,rho, theta): return 0
    def e_rCladding(self, rho, theta): return 0
    def e_thetaCladding(self, rho, theta):
        return -kv(1, self.W*rho)/ kv(1, self.W)
    def e_zCaldding(self, rho, theta): return 0
    
    def h_rCore(self, rho, theta):
        prefac = np.sqrt(epsilon_0/mu_0) * self.beta/ self.k0
        return prefac * jv(1,self.U*rho)/jv(1,self.U)
    def h_thetaCore(self, rho, theta): return 0
    def h_zCore(self, rho, theta):
        """Return h_zCore. Careful, in principle e_z is purely imaginary"""
        prefac= np.sqrt(epsilon_0/mu_0) * self.U/(self.k0*self.radius)
        return prefac*jv(0,self.U*rho)/jv(1,self.U)
    def h_rCladding(self, rho, theta): return 0
    def h_thetaCladding(self, rho, theta):
        prefac=np.sqrt(epsilon_0/mu_0) * self.beta/self.k0
        return prefac*kv(1,self.W*rho)/kv(1,self.W)
    def h_zCladding(self, rho, theta):
        prefac= np.sqrt(epsilon_0/mu_0) * self.W/(self.k0 * self.radius)
        return -prefac * kv(0,self.W * rho)/kv(1,self.W)
    
    def SzCore(self, rho, theta):
        prefac = 1./2*np.sqrt(epsilon_0/mu_0) * self.beta/self.k0
        return prefac * np.power(jv(1,self.U*rho)/jv(1,self.U),2)
    def SzCladding(self, rho, theta):
        prefac = 1./2*np.sqrt(epsilon_0/mu_0) * self.beta/self.k0
        return prefac*np.power(kv(1,self.W*rho)/kv(1,self.W),2)

class TMMode(Mode):
    """Class for TM modes"""

    def __init__(self, wavelength, diameter, matCore="silica", matCladd="air", nu=0):
        Mode.__init__(self, diameter, wavelength, matCore, matCladd, nu)
        self.typeMode="TM"
                 
    def eigenValuePb(self, U):
        w = np.sqrt(self.V**2 - np.power(U,2))
        nco2 = self.nCore**2
        ncl2 = self.nCladd**2
        y = nco2*jv(1,U)/(U*jv(0,U)) + ncl2*kv(1,w)/(w*kv(0,w))
        # y[y > 7] = None
        # y[y < -7] = None
        return y
    
    def nbreOfModes(self):
        n=1
        while jn_zeros(self.nu, n)[-1]<self.V:
            n+=1
        return n-1
    
    # def guessUnit(self,m):
    #     lowBoundary=0
    #     n=m
    #     mMax = self.nbreOfModes()
    #     if m>mMax:
    #         print("Mode TH0{0} at lambda = {1} does not exist".format(m, self.wavelength))
    #         print("lambda is above the cutoff")
    #         print("the number of possible {0} modes is {1}".format(self.typeMode, mMax))
    #         exit(0)

    #     # if self.typeMode=="EH":    
    #     #     lowBoundary = jn_zeros(self.nu,n)[-1]
    #     #     highBoundary = jn_zeros(self.nu,n+1)[-1]

    #     # if self.typeMode=="HE":
    #     #     if n!=1:
    #     #         lowBoundary = jn_zeros(self.nu,n-1)[-1]
    #     #         highBoundary = jn_zeros(self.nu,n)[-1]
    #     #     else: highBoundary = jn_zeros(1,1)[-1]
                                         
    #     if highBoundary>self.V: 
    #         highBoundary=self.V
            
    #     return 1./2 * (lowBoundary+highBoundary)

    
    def guessUnit(self,m):
            mMax = self.nbreOfModes()
            if m>mMax:
                print("Mode TH0{0} at lambda = {1} does not exist".format(m, self.wavelength))
                print("lambda is above the cutoff")
                print("the number of possible modes is {0}".format(mMax))
                exit(0)
            lowBoundary = jn_zeros(0,m)[-1]
            highBoundary = jn_zeros(0,m+1)[-1]
            return 1./2 * (lowBoundary+highBoundary)

    def show(self):
        print("This is a TM mode")
        Mode.show(self)
        
    def e_rCore(self, rho, theta):
        return jv(1,self.U*rho)/jv(1,self.U)
    def e_thetaCore(self, rho, theta): return 0
    def e_zCore(self, rho, theta):
        """Careful! in principle e_z is purely imaginary"""
        prefac=self.U /(self.radius * self.beta)
        return prefac*jv(0,self.U*rho)/jv(1,self.U)
    def e_rCladding(self,rho,theta):
        nCo2=self.nCore**2
        nCl2=self.nCladd**2
        return nCo2/nCl2 * kv(1, self.W*rho)/kv(1,self.W)
    def e_thetaCladding(self,rho,theta): return 0           
    def e_zCladding(self,rho,theta):
        nCo2=self.nCore**2
        nCl2=self.nCladd**2
        return -nCo2/nCl2 * self.W/(self.radius*self.beta) * kv(0,self.W*rho)/kv(1,self.W)
    def h_rCore(self, rho, theta): return 0
    def h_thetaCore(self, rho, theta):
        prefac=np.sqrt(epsilon_0/mu_0) * self.k0 * self.nCore**2 / self.beta
        return prefac*jv(1,self.U*rho)/jv(1,self.U)
    def h_zCore(self, rho, theta): return 0
    def h_rCladding(self, rho, theta): return 0
    def h_thetaCladding(self, rho, theta): 
        prefac= np.sqrt(epsilon_0/mu_0) * self.k0 * self.nCore**2 / self.beta
        return prefac*kv(1, self.W*rho)/kv(1,self.W)
    def h_zCladding(self, rho, theta): return 0
    
    def SzCore(self, rho, theta):
        prefac = 1./2*np.sqrt(epsilon_0/mu_0) * self.k0*self.nCore**2/self.beta
        return prefac * np.power(jv(1,self.U*rho)/jv(1,self.U),2)
    def SzCladding(self, rho, theta):
        prefac = 1./2*np.sqrt(epsilon_0/mu_0) * self.k0*self.nCore**2/(self.beta*(1-2.*self.delta))
        return prefac*np.power(kv(1,self.W*rho)/kv(1,self.W),2)


class HybridMode(Mode):
    """Class for HE/EH modes"""
      
    #radius=0
    def __init__(self, wavelength, diameter, matCore="silica", matCladd="air", nu=1, typeMode="HE"):
        Mode.__init__(self, diameter, wavelength, matCore, matCladd, nu)
        self.typeMode=typeMode
        
    def eigenValuePb(self, U):
        #print('U eigen:',U)
        if self.typeMode=='HE': res=self.HE_mode(U)
        else: res=self.EH_mode(U)
        return res
    
    def solveEigenValuePb(self, m):
        self.m = m
        Uinit=self.guessUnit(m)
        self.U = self.solve(Uinit)
        #print('U init: ', Uinit )
        self.W = np.sqrt(self.V * self.V - self.U * self.U)
        self.beta = np.sqrt((self.k0*self.nCore)**2 - (self.U/self.radius)**2)   
        self.compute_bi()
        self.computeFi()
        self.compute_ai()
        #print("U solve=",self.U)
        return self.U
    
    def nbreOfModes(self):
        n=1
        while jn_zeros(self.nu, n)[-1]<self.V:
            n+=1
        if (self.typeMode=="EH"):
            return n-1
        else:    
            return n
    
    def HE_mode(self, U):
        w = np.sqrt(self.V**2-U*U)
        # w_nan = np.isnan(w)
        # if w_nan.all() != False: #cleaning nan values
        #     w[w_nan] = 0.0 #with zero value
        nco2 = self.nCore**2
        ncl2 = self.nCladd**2
        R2 = ((nco2 - ncl2)/(2.*nco2))**2 *  np.power(kvp(self.nu,w,1)/(w*kv(self.nu,w)),2)+ self.nu**2 *(1./(U*U) + (ncl2/nco2)/(w*w)) * (1./(U*U) + 1./(w*w))                                                          
        KK = (nco2 + ncl2)/(2.*nco2) * kvp(self.nu,w,1) / (w*kv(self.nu,w))       
        y = jv(self.nu-1,U)/(U*jv(self.nu,U)) + KK - self.nu/(U*U) + np.sqrt(R2)
        #y[y.all() > 5] = np.NaN
        #y[y < -7] = None
        return y
    def EH_mode(self, U):
        w = np.sqrt(self.V**2-U*U)
        # w_nan = np.isnan(w)
        # if w_nan.all() != False: #cleaning nan values
        #     w[w_nan] = 0.0 #with zero value
        nco2 = self.nCore**2
        ncl2 = self.nCladd**2
        R2 = ((nco2 - ncl2)/(2.*nco2))**2 *  np.power(kvp(self.nu,w,1)/(w*kv(self.nu,w)),2)+ self.nu**2 *(1./(U*U) + (ncl2/nco2)/(w*w)) * (1./(U*U) + 1./(w*w))                                                          
        KK = (nco2 + ncl2)/(2.*nco2) * kvp(self.nu,w,1) / (w*kv(self.nu,w)) 
        y = jv(self.nu+1,U)/(U*jv(self.nu,U)) - KK - self.nu/(U*U) + np.sqrt(R2)
        #y[y.all() > 5] = None
        #y[y < -7] = None
        return y
    def guessUnit(self,m):
        lowBoundary=0
        n=m
        if self.typeMode=='EH':n+=1
        mMax = self.nbreOfModes()
        if m>mMax:
            print("Mode TH0{0} at lambda = {1} does not exist".format(m, self.wavelength))
            print("lambda is above the cutoff")
            print("the number of possible {0} modes is {1}".format(self.typeMode, mMax))
            exit(0)
            
        if (n>1):
            lowBoundary = jn_zeros(1,n-1)[-1]
            highBoundary = jn_zeros(1,n)[-1]
            
        if (n==1):
            highBoundary = jn_zeros(1,1)[-1]
                       
        if highBoundary>self.V: 
            highBoundary=self.V
        return 1./2 * (lowBoundary+highBoundary)
 
    def compute_ai(self):
        self.ai=np.zeros(6)
        self.ai[0]=(self.F2-1)/2.
        self.ai[1]=(self.F2+1)/2.
        self.ai[2]=(self.F1-1)/2.
        self.ai[3]=(self.F1+1)/2.
        self.ai[4]=(self.F1-1+2.*self.delta)/2.
        self.ai[5]=(self.F1+1-2.*self.delta)/2.
    
    def compute_bi(self):
        self.b1 = 1./(2.*self.U)*(jv(self.nu-1,self.U) - jv(self.nu+1, self.U))/jv(self.nu, self.U)
        self.b2 = -1./(2.*self.W)*((kv(self.nu-1,self.W)) + kv(self.nu+1, self.W))/kv(self.nu, self.W)
 
    def computeFi(self):
        p=np.power(self.U*self.W/self.V,2)
        self.F1=p*(self.b1+(1-2.*self.delta)*self.b2)/self.nu
        self.F2=self.nu/(self.b1 + self.b2)/p
                  
    def show(self):
        print("This is an hybrid mode {0}".format(self.typeMode))
        Mode.show(self)
        
    def e_rCore(self, rho, theta):
        """We have already rho=1 at boundary"""
        num=-(self.ai[0]*jv(self.nu-1, self.U*rho) + self.ai[1]*jv(self.nu+1, self.U*rho))
        denom=jv(self.nu, self.U)
        angular = np.cos(self.nu*theta)
        return (num/denom)*angular
    def e_thetaCore(self, rho, theta):
        """We have already rho=1 at boundary"""
        num=-(self.ai[0]*jv(self.nu-1, self.U*rho) - self.ai[1]*jv(self.nu+1, self.U*rho))
        denom=jv(self.nu, self.U)
        angular=-np.sin(self.nu*theta)
        return (num/denom)*angular
    def e_zCore(self, rho, theta):
        """In principle, e_z is purely imaginary"""
        num= - self.U*jv(self.nu, self.U*rho) # * I 
        denom=jv(self.nu, self.U)*self.radius*self.beta 
        angular= np.cos(self.nu*theta) # even mode
        return (num/denom)*angular
    
    def e_rCladding(self, rho, theta):
        num=-(self.U/self.W)*(self.ai[0]*kv(self.nu-1, self.W*rho)-self.ai[1]*kv(self.nu+1, self.W*rho))
        denom = kv(self.nu, self.W)
        angular = np.cos(self.nu*theta)
        return num/denom*angular
    def e_thetaCladding(self, rho, theta):
        num=-(self.U/self.W)*(self.ai[0]*kv(self.nu-1, self.W*rho)+self.ai[1]*kv(self.nu+1, self.W*rho))
        denom = kv(self.nu, self.W)
        angular = -np.sin(self.nu*theta)  # even mode
        return num/denom*angular
    
    def e_zCladding(self, rho, theta):
        num= -self.U*kv(self.nu, self.W*rho)
        denom=kv(self.nu, self.W)*self.radius*self.beta
        angular=np.cos(self.nu*theta)
        return num/denom*angular
    
    def h_rCore(self,rho,theta):
        num=self.ai[2]*jv(self.nu-1, self.U*rho)-self.ai[3]*jv(self.nu+1, self.U*rho)
        denom=jv(self.nu, self.U)
        prefac=np.sqrt(epsilon_0/mu_0) * self.k0 * self.nCore**2 / self.beta
        angular=-np.sin(self.nu*theta)
        return prefac*num/denom*angular
    def h_thetaCore(self, rho, theta):
        num=self.ai[2]*jv(self.nu-1, self.U*rho)+self.ai[3]*jv(self.nu+1, self.U*rho)
        denom=jv(self.nu, self.U)
        prefac=np.sqrt(epsilon_0/mu_0) * self.k0 * self.nCore**2 / self.beta
        angular = np.cos(self.nu*theta)
        return -prefac*num/denom*angular
    def h_zCore(self,rho,theta):
        num=self.U*self.F2*jv(self.nu, self.U*rho)
        denom=self.k0*self.radius*jv(self.nu, self.U)
        angular=-np.sin(self.nu*theta)
        prefac = np.sqrt(epsilon_0/mu_0)
        return -prefac*num/denom*angular
    def h_rCladding(self, rho, theta):
        num=self.ai[4]*kv(self.nu-1, self.W*rho)+self.ai[5]*kv(self.nu+1, self.W*rho)
        denom=kv(self.nu, self.W)
        prefac=np.sqrt(epsilon_0/mu_0) * self.k0 * self.nCore**2 / self.beta*self.U/self.W
        angular=-np.sin(self.nu*theta)
        return prefac*num/denom*angular
        
    def h_thetaCladding(self,rho,theta):
        num=self.ai[4]*kv(self.nu-1, self.W*rho)-self.ai[5]*kv(self.nu+1, self.W*rho)
        denom=kv(self.nu, self.W)
        prefac=np.sqrt(epsilon_0/mu_0) * self.k0 * self.nCore**2 / self.beta*self.U/self.W
        angular=np.cos(self.nu*theta)
        return -prefac*num/denom*angular  

    def h_zCladding(self, rho, theta):
        num=self.U*self.F2*kv(self.nu, self.W*rho)
        denom=self.k0*self.radius*kv(self.nu, self.W)
        angular= -np.sin(self.nu*theta)
        prefac=np.sqrt(epsilon_0/mu_0)
        return -prefac*num/denom*angular
  
    def SzCore(self,rho,theta):
        factor=np.sqrt(epsilon_0/mu_0)*self.k0*self.nCore**2
        evenOrOdd=+1
        a=self.ai[0]*self.ai[2]*np.power(jv(self.nu-1, self.U*rho),2)
        b=self.ai[1]*self.ai[3]*np.power(jv(self.nu+1,self.U*rho),2)
        c=(1-self.F1*self.F2)/2. * jv(self.nu-1, self.U*rho)*jv(self.nu+1, self.U*rho)*np.cos(2.*self.nu*theta)
        return factor * (a+b-evenOrOdd*c)/np.power(jv(self.nu,self.U),2)/self.beta
    
    def SzCladding(self, rho, theta):
        factor=np.sqrt(epsilon_0/mu_0)*self.k0*self.nCore**2
        evenOrOdd=+1
        a=self.ai[0]*self.ai[4]*np.power(kv(self.nu-1,self.W*rho),2)
        b=self.ai[1]*self.ai[5]*np.power(kv(self.nu+1,self.W*rho),2)
        c=(1-2*self.delta-self.F1*self.F2)/2.*kv(self.nu-1,self.W*rho)*kv(self.nu+1, self.W*rho)*np.cos(2*self.nu*theta)
        return factor*(a+b+evenOrOdd*c)/np.power(kv(self.nu,self.W),2)*np.power(self.U/self.W,2)/self.beta


def calculate_CB_mode(new_wavelength, new_diameter, new_mode, 
                new_plot_mode, new_mcore, new_mcladd, new_nu):
    if new_mode == None:  # For error while deselect of mode type, so new_mode = None
        new_mode = 'HE'
        my_newMode = HybridMode(new_wavelength, new_diameter, matCore=new_mcore, matCladd=new_mcladd, typeMode=new_mode, nu=new_nu)

    elif new_mode == 'TE':
        my_newMode = TEMode(new_wavelength, new_diameter, matCore=new_mcore, matCladd=new_mcladd)

    elif new_mode == 'TM':
        my_newMode = TMMode(new_wavelength, new_diameter, matCore=new_mcore, matCladd=new_mcladd)

    else:
        my_newMode = HybridMode(new_wavelength, new_diameter, matCore=new_mcore, matCladd=new_mcladd, typeMode=new_mode, nu=new_nu)

    return my_newMode

def matDropdowns(newMode):
    order_core = newMode.matCore.order
    order_cladd = newMode.matCladd.order

    availableMaterials=[
        {'label': 'GeO2_193', 'value': 'GeO2_193'},
        {'label': 'GeO2_63', 'value': 'GeO2_63'},
        {'label': 'P2O5', 'value': 'P2O5'},
        {'label': 'silica', 'value': 'silica'}, 
        {'label': 'B2O3', 'value': 'B2O3'}, #commented due to its behaviour
        {'label': 'air', 'value': 'air'},
        #add more
        ]
    return [availableMaterials[:order_cladd], availableMaterials[order_core+1:]] # == materials core, materials clad
    #return [availableMaterials[:order_cladd], availableMaterials[order_core+1:]] # == materials core, materials clad
    #return [availableMaterials[:order_cladd], [availableMaterials[-1]]] # == materials core, materials clad



def plot_beta_2(diameter_d, n_min, n_max, N, dependence ='omega', m=1, mode = HybridMode, mode_d ='HE'):
    beta_i = []
    if dependence == 'omega':
        w = np.linspace(n_min,n_max,N)
        wavelengths = 2*np.pi*speed_of_light/(w*10**9)
    elif dependence == 'lambda':
        wavelengths = np.linspace(n_max, n_min,N)
        w = 2*np.pi*speed_of_light/(wavelengths*10**9)
    else:
        print('Error: Select a valid option: either omega or lambda')
        pass
    for wavelength in wavelengths:
        def_mode=mode(wavelength, diameter_d, typeMode=mode_d, nu=1)
        beta_i.append(def_mode.computeBeta(m))
    
    #spl = UnivariateSpline(x, y, k=Degree of the smoothing spline, s= Positive smoothing factor used to choose the number of knots)
    #inter = UnivariateSpline(w, beta_i)#, k=2, s=0)
    inter = UnivariateSpline(np.flip(wavelengths), np.flip(beta_i))#, k=2, s=0)
    
    #Manually change the amount of smoothing:
    #inter.set_smoothing_factor(0)
    
    new_wl = np.linspace(1.7, 0.4,200)#micron
    #new_omega = 2*np.pi*speed_of_light/(new_wl*10**9)

    #beta_new = inter(new_omega)
    
        
    gruppenlaufzeit = inter.derivative()
    dispersion = gruppenlaufzeit.derivative()
    gruppenlaufzeit = dispersion.derivative() #Gruppenlaufzeit (linearer Phasenanteil)

    # beta_2= go.Scatter(x=np.flip(new_wl) ,y=dispersion(np.flip(new_wl)), name = 'Beta 2',
    #                         line=dict(color='red'))
    beta_2= go.Scatter(x=np.flip(new_wl) ,y=dispersion(np.flip(new_wl)), name = 'Beta 2',
                            line=dict(color='red'))
    Plot = go.Figure(data=[ beta_2])
    return Plot

#only beta 2

# myMode = HybridMode(1.596, 3, typeMode='HE', nu=1)        
# myMode.show()
# myMode.computeBeta(2)
# myMode.plotMode()

# myMode.plotEigenValuePb(0.5,30,250)      

# # print("=========")
# m=1
# print("beta for mode {0} = {1}".format(m,myMode.computeBeta(m)))
# # myMode.plotMode(-2,2,-2,2,512)  

# # myTEMode = TEMode(0.532,5)
# # myTEMode.show()
# # #myTEMode.plotEigenValuePb(0.5,30,250)
# # myTEMode.computeBeta(2)
# # myTEMode.plotMode()

# # # print("=========")
# # myTMMode = TMMode(0.532,5)
# # myTMMode.show()
# # #myTEMode.plotEigenValuePb(0.5,30,250)
# # myTMMode.computeBeta(3)
# # myTMMode.plotMode()


# print("=============")
# beta=myMode.computeBeta(1)
# print("beta = {0}".format(beta))
# myMode.drawCut()

