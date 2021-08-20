# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from numpy.core.numeric import Inf
import plotly.graph_objects as go   
import numpy as np
from init_variables import *
from init_variables import date

# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

from app import app

colors = { #Definition of colors for future use
    #'background': '#111111',
    #'text': '#FFFFFF',
    'text': '#111111',
    'background': '#FFFFFF',
    'circle': 'whitesmoke',
    'even': 'darkgreen',
    'odd': 'darkblue',
    'other': 'darkred', 
}

last_update = 'Last update: ' + date

#------------------Definitons--------------#

#------------------- Grid: --------------------#
#Input pulse
#Tmax = 5 ~ 10ps
T0 = 1E-12 #  duration of input
#for pulse width  --> Dispersive effects at T0 ~ 1ps 
N = 8196 #ammount of points 
dt = 750*T0/N #the 16 is to get a grid between -8 and 8 for T/T0   #The number before T0/N sets the time (and freq) frame!
T = np.arange(-N/2, N/2)*dt
#---------------- Parameters: -----------------#
zmax = 3.5/3 # km    Values for n/3 -> φNL max = n*π
#gamma = n2*wo/(speed*Aeff)
beta2_initial = 8.3#5.66099
beta3_initial = 10#ps^3/km
gamma_initial = 1 #1/(W*km)
P0 = (3*np.pi)#10E-3
alpha_initial = 0 #dB/km
z0 = 0
m0 = 2
C0 = 0 #Chirp parameter
pulsetype = 'Gaussian'# or Sech
#-------------------------------------------#

#-------- SLIDERS---------------

gamma_slider = dcc.Slider(
        id='gamma_slid',
        min=0.001,
        max=1,
        step=0.001,
        value=gamma_initial,
        marks={
        0.001: {'label': '0.001', 'style': {'color': colors['text']}},
        1: {'label': '1', 'style': {'color': colors['text']}}},
    )

PPower_slider = dcc.Slider(
        id='PP_slid',
        min=0.001,
        max=1,
        step=0.0001,
        value=P0,
        marks={
        0.001: {'label': '1 mW', 'style': {'color': colors['text']}},
        1: {'label': '1W', 'style': {'color': colors['text']}}},
    )

ms_slider = dcc.Slider(
        id='ms_slid',
        min=2,
        max=10,
        step=1,
        value=m0,
        marks={
        2: {'label': '2', 'style': {'color': colors['text']}},
        10: {'label': '10', 'style': {'color': colors['text']}}},
    )

# C_slider = dcc.Slider(
#         id='cs_slid',
#         min=0,
#         max=10,
#         step=1,
#         value=C0,
#         marks={
#         0: {'label': '0', 'style': {'color': colors['text']}},
#         10: {'label': '10', 'style': {'color': colors['text']}}},
#     )

#---------------------------------------------------------------#
#///////////////////////////////////////////////////////////////#
#-------------------                        --------------------#
spm1 = Propagation( T0, T, m = 1, 
                        C=C0, pulsetype = 'Gaussian',
                        solve_type='only_spm',  
                        beta2=beta2_initial,
                        gamma=gamma_initial, 
                        P0=P0)
spm2 = Propagation( T0, T, m = m0, 
                        C=C0, pulsetype = 'Gaussian',
                        solve_type='only_spm',  
                        beta2=beta2_initial,
                        gamma=gamma_initial, 
                        P0=P0)
spm3 = Propagation( T0, T, m = m0, 
                        C=C0, pulsetype = 'Sech',
                        solve_type='only_spm',  
                        beta2=beta2_initial,
                        gamma=gamma_initial, 
                        P0=P0)
LNL = spm1.compute_LNL()
spm_phase = plot_shift(spm1, spm2,spm3)
spm_chirp = plot_chirp(spm1, spm2,spm3)

phase_graph = dcc.Graph(id='phase_spm_plot',
                        animate=True,
                        figure=spm_phase.update_layout())

chirp_graph = dcc.Graph(id='chirp_spm_plot',
                        animate=True,
                        figure=spm_chirp.update_layout())
#---------------------------------------------------------------#

#-------Final Layout for plotting and display of sliders and constants-------#
layout = html.Div(style={'backgroundColor': colors['background']},
    children=[
        html.Div(id='dpf-display-value'),
        dcc.Link('Go to GVD effect', href='/apps/gvd'), ##EDIT LINK TO OTHER PAGES
        html.Div(id='modesf-display-value'),
        dcc.Link('Go to Split-Step for NLSE', href='/apps/nlse'),  
        html.Div(className='rows', 
        children=[
            html.Div(className='four columns div-user-controls', style={'backgroundColor': colors['background']}, 
            children = [
                html.P('Front-End for understanding the SPM contribution', style={'color': colors['text']}),
                html.H3('\u03B3: '+str(gamma_initial)+ ' [1/W/m]', id='gamma_val',style={'color': colors['text']}),
                gamma_slider,
                html.H3('P0: '+str(P0)+r'W', id='p0_val', style={'color': colors['text']}),
                PPower_slider,
                html.H3('LNL: '+ str( '%.3f' %  (LNL)) + ' m', id = 'LNL_display'),
                html.H3('Leff = LNL'),
                html.H3('m: '+ str(m0), id = 'm_display'),
                ms_slider,
                #html.H3('C: '+ str(C0), id = 'c_display'),
                #C_slider,
            ]),  # Define the left element
            html.Div(className='eight columns div-for-charts bg-grey', style={'backgroundColor': colors['background']},  
            children = [
                     html.H1(children='Self-Phase Modulation',
                            style={
                            'textAlign': 'center',
                            'color': colors['text']
                    }),
                    html.Div(style={
                        'textAlign': 'center',
                        'color': colors['text']
                    }),
                    phase_graph,
                    chirp_graph,
            ])  
        ]),
        html.Footer('Joly Nicolas, Aguilera Hernan. Max-Planck-Institut', style={'color': colors['text']}),
        html.Footer(last_update, style={'color': colors['text']}),
    ])

@app.callback( 
    [Output('gamma_val', 'children'),
    Output('p0_val', 'children'),
    Output('LNL_display', 'children'),
    Output('m_display', 'children'),
    #Output('c_display', 'children'),
    Output('phase_spm_plot', 'figure'),
    Output('chirp_spm_plot', 'figure'),
    ],
    [Input('gamma_slid', 'value'),
    Input('PP_slid', 'value'),
    Input('ms_slid', 'value'),
    #Input('cs_slid', 'value'),
    ]
    )

# function to update with the callback:
def update_plot(new_gamma, new_pp, mn):#,cn):
    
    spm1 = Propagation( T0, T, m = 1, 
                            #C=cn, pulsetype = 'Gaussian',
                            C=C0, pulsetype = 'Gaussian',
                            solve_type='only_spm',
                            gamma=new_gamma, 
                            P0=new_pp)
    spm2 = Propagation( T0, T, m = mn, 
                            #C=cn, pulsetype = 'Gaussian',
                            C=C0, pulsetype = 'Gaussian',
                            solve_type='only_spm',
                            gamma=new_gamma, 
                            P0=new_pp)
    spm3 = Propagation( T0, T,
                            #C=cn, pulsetype = 'Sech',
                            C=C0, pulsetype = 'Sech',
                            solve_type='only_spm',
                            gamma=new_gamma, 
                            P0=new_pp)
    LNL = spm1.compute_LNL()
    spm_phase = plot_shift(spm1, spm2,spm3)
    spm_chirp = plot_chirp(spm1, spm2,spm3)

    return ('\u03B3: '+str(new_gamma)+ ' [1/W/m]',
            'P0: '+str(new_pp)+r' W',
            'LNL: '+ str( '%.3f' %  (LNL)) + ' m',
            'm: '+ str(mn),
            #'C: '+ str((cn)),
            spm_phase,
            spm_chirp
            )
