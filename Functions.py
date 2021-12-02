# coding:utf8

import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

#--------To send and plot--------#
#space within we are going to work (size of the vectors we will plot)
space=np.arange(0.01,((1/2)*np.pi)-0.01,0.01)


### flags and method added
#IMPORTANT: a = d/2 
def create_constants(a,lambda_1,n1,n2):
    ko=2*np.pi/lambda_1 # = 2*pi*frec/c ##   ----->> with d = 2a
    #ko=np.pi/lambda_1 # = *pi*frec/c     ##   ----->> with d = a
    V = np.sqrt(((ko*a)**2)*((n1**2)-(n2**2)))
    N = int(2*V/np.pi+ 1) #number of modes
    #N = int(V//(np.pi)+ 1) #number of modes
    flags = [] 

    #Array U 2D that contains N*2 vectors, each one with size of vector "space" to plot and find the solutions:
    #U=np.zeros([N*2, space.shape[0]])
    U=np.zeros([N, space.shape[0]])

    ##for each mode are created the vectors U1, U2, etc., within the amount of modes.
    #for i in range(0,int(N)*2,1):#for i in range(0,int(N)*2,1):
    #Possible by indexing the array ->future Optimization fix
    for i in range(0,N,1):
        if i % 2 == 0:
            flags.append(True)
            if (((i/2)+1/2)*np.pi) < V:
                U[i] = np.arange(((i/2)*np.pi)+0.01,(((i/2)+1/2)*np.pi)-0.01,0.01, dtype=np.float64 )
            else:
                U[i] = np.linspace(((i/2)*np.pi)+0.01,V, space.shape[0], dtype=np.float64 )
            
        else:
            flags.append(False)
            if (((i+1)/2)*np.pi) < V:
                U[i] = np.arange(((((i+1)/2)-1/2)*np.pi)+0.01,(((i+1)/2)*np.pi)-0.01,0.01, dtype=np.float64 )
            else:
                U[i] = np.linspace(((((i+1)/2)-1/2)*np.pi)+0.01,V, space.shape[0], dtype=np.float64 )
            
    alpha = U/a        
    return alpha, flags, ko, N, U, V


#^^^^^^^^To send and plot^^^^^^^^#


#if G = sqrt(V^2-U^2) 
#and W = U*tan(U) for even modes 
#and W = -U*cotan(U) for odd modes, then:
#f = G - W

####--------------------------FUNCTIONS--------------------------####

#Function 'G' = sqrt(V^2-U^2)
def G(a,V,alpha):
    VG = np.sqrt((V**2)-((alpha*a)**2))
    VG_nan = np.isnan(VG)
    if VG_nan.all() != False: #cleaning nan values
        VG[VG_nan] = 0.0 #with zero value
    return VG

#---functions W for even and odd modes---#
def W_even(a,alpha):
    return alpha*a*np.tan(alpha*a)
    
def W_odd(a,alpha):
    return -alpha*a*(1/np.tan(alpha*a))   
#--------------------------------------#
 
#----function for calculate G-W in either even or odd mode ----#
def f_TE(alpha_val,a,flag,V):
    if flag:
        return  G(a,V,alpha_val)-W_even(a,alpha_val)  # sqrt(V^2-U^2 )-W for find roots (even)
    if not flag:
        return G(a,V,alpha_val)-W_odd(a,alpha_val) #for odd mode
#-------------------------------------------#       


#----function for calculate G-W in either even or odd mode ----#
def f_TM(alpha_val,a,flag,n1,n2,V):
    if flag:
        return  G(a,V,alpha_val)-((n2**2)/(n1**2))*W_even(a,alpha_val)
    if not flag:
        return G(a,V,alpha_val)-((n2**2)/(n1**2))*W_odd(a,alpha_val)
#-------------------------------------------#       


#------Functions to obtain beta from either alpha or gamma---------#
def beta(ko,n1,alpha): #obtain beta from alpha
    return np.sqrt((ko**2)*(n1**2)-(alpha**2))/ko  #Working with normalized propagation constant
    
def beta2(a,ko,n2,W): #obtain beta from gamma and gamma from W   #there is a problem with this function
    gamma = W/a
    return np.sqrt((gamma**2)+(ko**2)*(n2**2))/ko
#-----------------------------------------------------------------#


#------------Bisection def for obtaining root alpha ---------#
def bisection(alpha,a,U1,U2,tolerance,f,flag,ko,n1,n2,V):
    val1 = U1/a
    val2 = U2/a
    while (np.abs(val1-val2) >=tolerance):
        c = (val1+val2)/2.0
        if f == f_TE:
            prod = f(val1,a,flag,V)*f(c,a,flag,V)
        if f == f_TM:
            prod = f(val1,a,flag,n1,n2,V)*f(c,a,flag,n1,n2,V)
        if prod > tolerance:
            val1 = c
        else:
            if prod < tolerance:
                val2 = c
    #print('alpha = c:', c )                
    print('beta_alpha:', beta(ko,n1,c))
    print('beta_gamma:', beta2(a,ko,n2,G(a,V,c)))
    return beta(ko,n1,c)
#------------------End Bisection---------------------------#


#-----------Add all Beta answers to an array:----------------
#--For TE Mode---# 
def beta_TE(a, alpha, beta, flags, ko, n1, V):
    beta_w = np.zeros(alpha.shape[0])
    for i in range(0,alpha.shape[0],1):
        x_lim = alpha[i,-20]-0.001
        solv_ev = fsolve(f_TE,x_lim,args=(a,flags[i],V))  #fsolve(func, x0, args=())
        beta_w[i]=beta(ko,n1,solv_ev[0])
    return beta_w
#----------------#
#--For TM Mode---#
def beta_TM(a, alpha, beta, flags, ko, n1, n2, V,):
    beta_w = np.zeros(alpha.shape[0])
    for i in range(0,alpha.shape[0],1):
        x_lim = alpha[i,-20]-0.001
        solv_ev = fsolve(f_TM,x_lim,args=(a,flags[i], n1, n2, V))  #fsolve(func, x0, args=())
        beta_w[i]=beta(ko,n1,solv_ev[0])
    return beta_w
#----------------#
#-------------------------------------------------------------#

#Obtaining values from beta and Alpha for plot of the wave functions
#(only TE Mode)
def beta_U(a, alpha, beta, flags, ko, n1, V):
    beta_w = np.zeros(alpha.shape[0])
    alpha_solved = np.zeros(alpha.shape[0])
    for i in range(0,alpha.shape[0],1):
        x_lim = alpha[i,-20]-0.001
        solv_ev = fsolve(f_TE,x_lim,args=(a,flags[i],V))  #fsolve(func, x0, args=())
        alpha_solved[i] = solv_ev[0]
        beta_w[i]=beta(ko,n1,solv_ev[0])
    return beta_w, alpha_solved
#----------------------------------------------------------------#





#--------Derivatives from beta interpolation--------# (Not finished)
#(Taylor Reihenfolge)
def beta_tay(beta_w,lambda_1):
    n_eff= (beta_w*lambda_1)/(2*np.pi)
    #wo = 3*2*np.pi/lambda_1  #*(10**8)
    #beta0 = 2*np.pi/lambda_1

    #invert vetor order:
    beta_w = np.flip(beta_w)#[::-1]
    beta_2w = beta_w + beta_w[-1]
    #beta_w = np.append(beta_w, beta0)  #add beta0 value
    beta_w= np.append(beta_w, beta_2w) #beta with central value Beta0  (-beta_2, -beta_1, beta0, +beta_1, +beta_2 )

    w = 3*beta_w #*(10**8)
    inter=interp1d(w, beta_w)

    new_w=np.arange(w[0],w[-1],0.2)

    
    beta1 = (beta_w[1]-beta_w[0])/(w[1]-w[0])
    beta2 = (beta_w[2]-2*beta_w[1]+beta_w[0])/((w[1]-w[0])**2)

    #backward
    # beta1 = (beta_w[0]-beta_w[1])/(w[0]-w[1])
    # beta2 = (beta_w[0]-2*beta_w[1]+beta_w[2])/((w[0]-w[1])**2)
    # beta1 = (wo-beta_w[0])/(wo-w[0])
    # beta2 = (wo-2*beta_w[0]+beta_w[1])/((wo-w[0])**2)

    lambda_f =2*np.pi/beta_w

    return beta1, beta2, w, lambda_f, new_w, inter(new_w)
#---------------------------------------------------------------------------#


#-----------Used on Plot Beta vs V ----------------------------------------#
def modes(a, beta_mode, n1, n2, normalized_freq):

    lambda_V = 2*np.pi*a*np.sqrt((n1**2)-(n2**2))/normalized_freq #for 2*a = d
    #lambda_V = np.pi*a*np.sqrt((n1**2)-(n2**2))/normalized_freq #for a = d

    V =  np.amax(normalized_freq)

    N_max = int(2*V/np.pi+ 1) #number of modes
    
    T = np.zeros([N_max, normalized_freq.shape[0]]) + n2

    for i in range(0, lambda_V.shape[0], 1):
        alpha, flags, ko, N, U, V = create_constants(a,lambda_V[i],n1,n2)
        if beta_mode == beta_TE:
            beta_values = beta_TE(a, alpha, beta, flags, ko, n1, V)
        if beta_mode == beta_TM:
            beta_values = beta_TM(a, alpha, beta, flags, ko, n1, n2, V)
        for Num in range(0, N, 1):
            T[Num,i]= beta_values[Num]    

    return T, lambda_V
#--------------------------------------------------------------------------#

#####-------------------------------------####
#####------------------------------####
#####----------------------####
#####--------------####
