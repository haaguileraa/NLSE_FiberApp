# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
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
dt =  750*T0/N#100*T0/N #the 16 is to get a grid between -8 and 8 for T/T0   #The number before T0/N sets the time (and freq) frame!
T = np.arange(-N/2, N/2)*dt
#---------------- Parameters: -----------------#
zmax = 3.5/3 # km    Values for n/3 -> φNL max = n*π
beta2_initial = 8.3#5.66099
alpha_initial = 0 #dB/km
z0 = 0
m0 = 1
C0 = 0 #Chirp parameter
pulsetype = 'Gaussian'# or Sech
#-------------------------------------------#

#-------- SLIDERS---------------
beta2_slider = dcc.Slider(
        id='beta2_slid',
        min=-20.0,
        max=20.0,
        step=0.1,
        value=beta2_initial,
        marks={
        -20: {'label': '-20', 'style': {'color': colors['text']}},
        0: {'label': '0', 'style': {'color': colors['text']}},
        20: {'label': '20', 'style': {'color': colors['text']}}},
    )

m_slider = dcc.Slider(
        id='mg_slid',
        min=1,
        max=10,
        step=1,
        value=m0,
        marks={
        1: {'label': '1', 'style': {'color': colors['text']}},
        10: {'label': '10', 'style': {'color': colors['text']}}},
    )

C_slider = dcc.Slider(
        id='cg_slid',
        min=-5,
        max=5,
        step=1,
        value=C0,
        marks={
        -5: {'label': '-5', 'style': {'color': colors['text']}},
        5: {'label': '5', 'style': {'color': colors['text']}}},
    )

Lange_slider = dcc.Slider(
        id='lange_slid',
        min=0,
        max=2.5,
        step=0.001,
        value=z0,
        marks={
        0: {'label': '0', 'style': {'color': colors['text']}},
        2.5: {'label': '2.5 km', 'style': {'color': colors['text']}}},
    )

beta2_initial *= 1E-24

pulse = Propagation( T0, T, m = m0, 
                        C=C0, pulsetype = 'Gaussian',
                        solve_type='only_gvd',  
                        beta2=beta2_initial,
                        z0=z0,
                        )
LD = pulse.compute_LD()
z1 = 2*LD
z2 = 4*LD

pulse1 = Propagation( T0, T, m = m0, 
                        C=C0, 
                        solve_type='gauss_gvd',  
                        beta2=beta2_initial,
                        z0=z0,
                        )

pulse2 = Propagation( T0, T, m = m0, 
                    C=C0, 
                    solve_type='gauss_gvd',  
                    beta2=beta2_initial,
                    z0=z1,
                    )
pulse3 = Propagation( T0, T, m = m0, 
                    C=C0,
                    solve_type='only_gvd',  
                    beta2=beta2_initial,
                    z0=z2,
                    )
#-------------------Using Eq 3.2.5 and 3.2.6--------------------#
gvd_t = pulse.plot_envelope_GVD(mode = 'time') # envelope plot in time
gvd_w = pulse.plot_envelope_GVD(mode = 'spectrum') # envelope plot in time  
gvd_graph = dcc.Graph(id='gvdt', #id for callback purposes
                        animate=True,
                        figure=gvd_t.update_layout(
))
gvd_graph_frec = dcc.Graph(id='gvdw', #id for callback purposes
                        animate=True,
                        figure=gvd_w.update_layout(
))
#---------------------------------------------------------------#
#///////////////////////////////////////////////////////////////#
#-----------------------Eq. 3.2.7 and 3.2.9---------------------#

gvd_t1 = go.Scatter(x=pulse1.T/pulse1.T0,y=pulse1.UI,name = 'z = 0',
                                line=dict(color='darkgreen'))
gvd_t2 = go.Scatter(x=pulse1.T/pulse1.T0,y=pulse2.UI,name = 'z = 2*LD',
                                line=dict(color='darkblue'))
gvd_t3 = go.Scatter(x=pulse1.T/pulse1.T0,y=pulse3.UI, name = 'z = 4*LD',
                                line=dict(color='darkred'))

env_fig = go.Figure(data=[gvd_t1, gvd_t2, gvd_t3]).update_layout()

plot_gvd = dcc.Graph(id='gvdt1', #id for callback purposes
                        animate=True,
                        figure=env_fig.update_layout( 
                            width=600, height=600,
                            plot_bgcolor  = colors['background'],
                            paper_bgcolor = colors['background'],
                            font= {
                                    'color': colors['text']},
                            yaxis=dict(range=[0, 1.1],title='|U(z,T)|^2', 
                                        ), 
                            xaxis=dict(range=[-8, 8],title='T/T0', 
                                        ), 
                            title= '''
                                    Dispersion-induced broadening of a Gaussian pulse.
                                    ''',
                            )
          )
#---------------------------------------------------------------#

#-------Final Layout for plotting and display of sliders and constants-------#
layout = html.Div(style={'backgroundColor': colors['background']},
    children=[
        html.Div(id='dpf-display-value'),
        dcc.Link('Go to SPM effect', href='/apps/spm'), ##EDIT LINK TO OTHER PAGES
        html.Div(id='modesf-display-value'),
        dcc.Link('Go to Split-Step for NLSE', href='/apps/nlse'),
        
        html.Div(className='rows', 
        children=[
            html.Div(className='four columns div-user-controls', style={'backgroundColor': colors['background']}, 
            children = [
                html.P('Front-End for understanding the GVD contribution', style={'color': colors['text']}),
                html.H3(children=[html.Var('\u03B22: '+str('%.3f' % (beta2_initial))+ ' ps', id='beta2_val'),
                html.Sup(2), html.Var('/km')],
                style={'color': colors['text']}),
                beta2_slider,
                html.H3('LD: '+ str( '%.3f' %  (LD)) + ' km', id = 'LD_display'),
                html.H3('m: '+ str(m0), id = 'mg_val'),
                m_slider,
                html.H3('C: '+ str(C0), id = 'cg_val'),
                C_slider,
                html.H3('z: '+str(z0)+r'km', id='z_val', style={'color': colors['text']}),
                Lange_slider,
                daq.ToggleSwitch(id='switchgvd',value=False,label='Type of pulse: ', labelPosition='bottom'),
            ]),  # Define the left element
            html.Div(className='eight columns div-for-charts bg-grey', style={'backgroundColor': colors['background']},  
            children = [
                     html.H1(children='Group-Velocity Dispersion',
                            style={
                            'textAlign': 'center',
                            'color': colors['text']
                    }),

                    html.Div(style={
                        'textAlign': 'center',
                        'color': colors['text']
                    }),
                    html.H2('Using Eq. 3.2.5 and 3.2.6 Agrawal'),
                    #html.H3('Time'),
                    gvd_graph,
                    #html.H3('Spectrum'),
                    gvd_graph_frec,
                    html.H2('Using Eq. 3.2.7 and 3.2.9 Agrawal'),
                    plot_gvd,                  
            ])  
        ]),
        html.Footer('Joly Nicolas, Aguilera Hernan. Max-Planck-Institut', style={'color': colors['text']}),
        html.Footer(last_update, style={'color': colors['text']}),
    ])

@app.callback( 
    [Output('beta2_val', 'children'),
    Output('z_val', 'children'),
    Output('mg_val', 'children'),
    Output('cg_val', 'children'),
    Output('LD_display', 'children'),
    Output('gvdt', 'figure'),
    Output('gvdw', 'figure'),
    Output('switchgvd', 'label'),
    ],
    [Input('beta2_slid', 'value'),
    Input('mg_slid', 'value'),
    Input('cg_slid', 'value'),
    Input('lange_slid', 'value'),
    Input('switchgvd', 'value')]
    )

# function to update with the callback:
def update_plot(new_beta2, mn, cn, new_z, switch):
    new_beta2 *= 1E-24
    if switch: pulsetype = 'Sech'
    else: pulsetype = 'Gaussian'
    pulse = Propagation( T0, T, m = mn, 
                        C=cn, pulsetype = pulsetype,
                        solve_type='only_gvd',  
                        beta2=new_beta2,
                        z0=new_z,
                        )

    return ('\u03B22: '+ str( '%.3f' %  (new_beta2*1E+24)) + ' ps',
            'z: '+str(new_z)+r'km',
            'm: '+ str(mn),
            'C: '+ str((cn)),
            'LD: '+ str( '%.3f' %  (pulse.compute_LD())) + ' km',
            pulse.plot_envelope_GVD(mode = 'time'), # envelope plot in time,
            pulse.plot_envelope_GVD(mode = 'spectrum'), # envelope plot in spectrum 
            'Type of pulse: {0}'.format(pulsetype)
            )
