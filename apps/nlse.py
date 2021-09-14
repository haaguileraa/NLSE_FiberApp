# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
#import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from numpy.core.numeric import Inf
import plotly.graph_objects as go   
import numpy as np
from init_variables import *
from init_variables import date

#------------------Definitons--------------#



from app import app


last_update = 'Last update: ' + date
colors = { #Colors for the plots and web interface:
            'text': '#111111',
            'background': '#FFFFFF',
            'circle': 'whitesmoke',
            'even': 'darkgreen',
            'odd': 'darkblue',
            'other': 'darkred',
            }

#------------------- Grid: --------------------#
#Input pulse
#Tmax = 5 ~ 10ps -> much less GVD -> LD >> L and LNL
T0 = 1E-12 #  duration of input #slider 50fs -1ps  -grid-8ps-8ps
#for pulse width  --> Dispersive effects at T0 ~ 1ps 
N = 8192 #ammount of points 
dt = 750*T0/N #the 16 is to get a grid between -8 and 8 for T/T0   #The number before T0/N sets the time (and freq) frame!
T = np.arange(-N/2, N/2)*dt
#---------------- Parameters: -----------------#
zmax = 500#m #3.5/3*1E2 # km    Values for n/3 -> φNL max = n*π
#gamma = n2*wo/(speed*Aeff)
beta2_initial = 8.3#5.66099 #ps^2/km
beta3_initial = 0#ps^3/km
gamma_initial = 1 #1/(W*km) 
P0 = (3*np.pi)#10E-3 # W
alpha_initial = 0 #dB/km
z0 = 0 #initial point for pulse envelope
m0 = 1 #for Gaussian/Super-Gaussian pulses
C0 = 0 #Chirp parameter
pulsetype = 'Gaussian'# or Sech
size_array = 13
lastz = int(size_array-1) #last value to create the slider: from 0 to size_array-1 -> indey never out of bounds
#-------------------------------------------#

#-------- SLIDERS---------------
T_slider = dcc.Slider(
        id='T_slid',
        min=50,
        max=1000,
        step=10,
        value=zmax,
        marks={
        50: {'label': '50 fs', 'style': {'color': colors['text']}},
        1000: {'label': '{0} ps'.format('%.2f' % (1)), 'style': {'color': colors['text']}}},
    )



alpha_slider = dcc.Slider(
        id='alpha_slider',
        min=0,
        max=10,
        step=0.001,
        value=alpha_initial,
        marks={
        0: {'label': '0', 'style': {'color': colors['text']}},
        10: {'label': '10', 'style': {'color': colors['text']}}},
    )

beta2_slider = dcc.Slider(
        id='beta2_slider',
        min=-25.0,
        max=25.0,
        step=0.1,
        value=beta2_initial,
        marks={
        -25: {'label': '-25', 'style': {'color': colors['text']}},
        0: {'label': '0', 'style': {'color': colors['text']}},
        25: {'label': '25', 'style': {'color': colors['text']}}},
    )

beta3_slider = dcc.Slider(
        id='beta3_slid',
        min=-20.0,
        max=20.0,
        step=0.1,
        value=beta3_initial,
        marks={
        -20: {'label': '-20', 'style': {'color': colors['text']}},
        0: {'label': '0', 'style': {'color': colors['text']}},
        20: {'label': '20', 'style': {'color': colors['text']}}},
    )

z_slider = dcc.Slider(
        id='z_slid',
        min=0,
        max=lastz,
        step=1,
        value=z0,
        marks={
        0: {'label': '0', 'style': {'color': colors['text']}},
        lastz: {'label': '{0} m'.format('%.2f' % (zmax)), 'style': {'color': colors['text']}}},
    )

L_slider = dcc.Slider(
        id='L_slid',
        min=0.1,
        max=500,
        step=0.1,
        value=zmax,
        marks={
        0.1: {'label': '0.1', 'style': {'color': colors['text']}},
        500: {'label': '{0} m'.format('%.2f' % (500)), 'style': {'color': colors['text']}}},
    )

gamma_slider = dcc.Slider(
        id='gamma_slider',
        min=0.01,
        max=10,
        step=0.01,
        value=gamma_initial,
        marks={
        0.01: {'label': '0.01', 'style': {'color': colors['text']}},
        10: {'label': '10', 'style': {'color': colors['text']}}},
    )

Power_slider = dcc.Slider(
        id='P0_slid',
        min=0.001,
        max=10,
        step=0.0001,
        value=P0,
        marks={
        0.001: {'label': '1 mW', 'style': {'color': colors['text']}},
        10: {'label': '10W', 'style': {'color': colors['text']}}},
    )

m_slider = dcc.Slider(
        id='m_slid',
        min=1,
        max=10,
        step=1,
        value=m0,
        marks={
        1: {'label': '1', 'style': {'color': colors['text']}},
        10: {'label': '10', 'style': {'color': colors['text']}}},
    )

C_slider = dcc.Slider(
        id='c_slid',
        min=-5,
        max=5,
        step=1,
        value=C0,
        marks={
        -5: {'label': '-5', 'style': {'color': colors['text']}},
        5: {'label': '5', 'style': {'color': colors['text']}}},
    )


#-----------------------------------------------------------------------#
beta2_initial *= 1E-27 # 1E-24 * 1E-3 (ps²/m)
beta3_initial *= 1E-39 # 1E-36 * 1E-3 (ps³/m)
gamma_initial *= 1E-3 # (1/(W m))

#if beta2_initial == 0:
pulse = Propagation(T0, T, m = m0, 
                    C=C0, pulsetype = pulsetype,
                    solve_type='split_step', 
                    L=zmax, 
                    beta2=beta2_initial,
                    beta3=beta3_initial,
                    gamma=gamma_initial, 
                    P0=P0,
                    z0=z0,
                    size_array = size_array)
# else:
#     pulse = Propagation(T0, T, m = m0, 
#                         C=C0, pulsetype = pulsetype,
#                         solve_type='split_step', 
#                         L=zmax, 
#                         beta2=beta2_initial,
#                         gamma=gamma_initial, 
#                         P0=P0,
#                         z0=z0,
#                         size_array = size_array)
UI = pulse.UI
UIW = pulse.UIW
Z = pulse.z
W = pulse.W
LNL = pulse.compute_LNL()
LD = pulse.compute_LD()

env_fig_t = pulse.plot_envelope(mode = 'time') # envelope plot in time 

env_graph_t = dcc.Graph(id='envelopet', #id for callback purposes
                        animate=True,
                        figure=env_fig_t.update_layout(
))

env_fig_w = pulse.plot_envelope(mode = 'spectrum') # envelope plot in spectrum 

env_graph_w = dcc.Graph(id='envelopew', #id for callback purposes
                        animate=True,
                        figure=env_fig_w.update_layout(
))  

import sys

print('SIZE: ',sys.getsizeof(UI),sys.getsizeof(UIW), sys.getsizeof(Z))
print('Size: ', UI.shape, UIW.shape, Z.shape)
print('Type: ', UI.dtype, UIW.dtype, Z.dtype)

# #It's generally safe to store up to 2MB in most environments, and 5~10MB in most desktop-only applications.
# # The memory store reverts to the default on every page refresh
# dcc.Store(id='memory'),
# # The local store will take the initial data
# # only the first time the page is loaded
# # and keep it until it is cleared.
# dcc.Store(id='local', storage_type='local'),
# # Same as the local store but will lose the data
# # when the browser/tab closes.
# dcc.Store(id='session', storage_type='session'),


prop_t = pulse.plot_propagation(mode = 'time') # propagation plot in time 
prop_w = pulse.plot_propagation(mode = 'spectrum') # propagation plot in spectrum 
   
#-------Final Layout for plotting and display of sliders and constants-------#
layout = html.Div(style={'backgroundColor': colors['background']},
    children=[
        dcc.Store(id='session', storage_type='local'), #problems with 'session' and 'local'. 'Memory' to slow
        html.Div(id='dpf-display-value'),
        dcc.Link('Go to SPM effect', href='/apps/spm'), ##EDIT LINK TO OTHER PAGES
        html.Div(id='modesf-display-value'),
        dcc.Link('Go to GVD effect', href='/apps/gvd'), ##EDIT LINK TO OTHER PAGES
        
        html.Div(className='rows', 
        children=[
            html.Div(className='four columns div-user-controls', style={'backgroundColor': colors['background']}, 
            children = [
                html.P('Front-End for NLSE', style={'color': colors['text']}),
                html.P('Pulse with T0 = 1ps', style={'color': colors['text']}),
                html.H3('\u03B1: '+str(alpha_initial)+ ' [dB/km]', id='alpha_v',style={'color': colors['text']}),
                alpha_slider,
                html.H3(children=[html.Var('\u03B22: '+str( '%.1f' % (beta2_initial*1E24))+ ' ps', id='beta2_v'),
                html.Sup(2), html.Var('/km')],
                style={'color': colors['text']}),
                beta2_slider,
                html.H3(children=[html.Var('\u03B23: '+str( '%.1f' % (beta3_initial*1E36))+ ' ps', id='beta3_v'),
                html.Sup(3), html.Var('/km')],
                style={'color': colors['text']}),
                beta3_slider,
                html.H3('\u03B3: '+str(gamma_initial)+ ' [ 1/(W km)]', id='gamma_v',style={'color': colors['text']}),
                gamma_slider,
                html.H3('P0: '+str('%.3f' % (P0))+r' W', id='p0_v', style={'color': colors['text']}),
                Power_slider,
                html.H3('m: '+ str(m0), id = 'm_val'),
                m_slider,
                html.H3('C: '+ str(C0), id = 'c_val'),
                C_slider,
                html.H5('LNL: '+ str( '%.3f' %  (LNL)) + ' m', id = 'LNL_v'),
                html.H5('LD: '+ str( '%.3f' %  (LD)) + ' m', id = 'LD_v'),
                html.H5('N = sqrt(LD/LNL)= '+ str( '%.3f' %  np.sqrt(LD/LNL)), id = 'Ng_v'),
                html.H3('L: '+str('%.4f' % (zmax))+r'm', id='L_', style={'color': colors['text']}),
                L_slider,
                daq.ToggleSwitch(id='switch',value=False,label='Type of pulse (Sech or Gaussian): ', labelPosition='top'),
                html.Button('Calculate', id='calculate', n_clicks=0),
            ]),  # Define the left element
            html.Div(className='eight columns div-for-charts bg-grey', style={'backgroundColor': colors['background']},  
            children = [
                     html.H1(children='Solution for the NLSE',
                            style={
                            'textAlign': 'center',
                            'color': colors['text']
                    }),
                    html.Div(style={
                        'textAlign': 'center',
                        'color': colors['text']
                    }),
                    html.H3('Colormap time'),
                    dcc.Loading(id="loading1",children=[html.Div([dcc.Graph(id='propagationt', #id for callback purposes
                                                                            animate=True,
                                                                            figure=prop_t.update_layout()
                                                                            )   ])
                                                                ],type="circle",),                                       
                    html.H3('Colormap spectrum'),
                    dcc.Loading(id="loading2",children=[html.Div([dcc.Graph(id='propagationw', #id for callback purposes
                                                                            animate=True,
                                                                            figure=prop_w.update_layout()
                                                                            )   ])
                                                                ],type="circle",), 
                    html.H3('Envelope time'),   
                    env_graph_t,
                    html.H3('z: '+str(pulse.z[z0])+r'm', id='z_', style={'color': colors['text']}),
                    z_slider,
                    html.H3('Envelope spectrum'),  
                    env_graph_w,          
            ])  
        ]),
        html.Button("Download Parameters as Text", id="download"), dcc.Download(id="download-text"),
        html.Footer('Joly Nicolas, Aguilera Hernan. Max-Planck-Institut', style={'color': colors['text']}),
        html.Footer(last_update, style={'color': colors['text']}),
    ])

@app.callback(
    dash.dependencies.Output('switch', 'label'),

    [dash.dependencies.Input('switch', 'value')])
def update_output(value):
    if value: return 'Type of pulse (Sech or Gaussian): Sech'
    else:  return 'Type of pulse (Sech or Gaussian): Gaussian'
    

@app.callback( 
    [Output('beta2_v', 'children'),
    Output('beta3_v', 'children'),
    Output('gamma_v', 'children'),
    Output('p0_v', 'children'),
    Output('LD_v', 'children'),
    Output('LNL_v', 'children'),
    Output('Ng_v', 'children'),
    ],
    [Input('beta2_slider', 'value'),
    Input('beta3_slid', 'value'),
    Input('gamma_slider', 'value'),
    Input('P0_slid', 'value'),
    ]
    )

# function to update with the callback:
#def update_left(new_alpha, new_beta2, new_beta3, new_zmax, new_gamma, new_p0, new_m):
def update_left(new_beta2, new_beta3, new_gamma, new_p0):
    if new_beta2 != 0:
        LD = ((T0)**2)/np.absolute(new_beta2*1E-27)  # Einheit!!!: km #HERE 1E-24 due to the square of ps^2
        mes = 'LD '+ str( '%.3f' %  (LD)) + ' m'
    elif new_beta3 != 0 and new_beta2 == 0:
        LD = ((T0)**3)/np.absolute(new_beta3*1E-39)  # Einheit!!!: km #HERE 1E-36 due to ps^3
        mes = "L'D: "+ str( '%.3f' %  (LD)) + ' m'
    else:
        LD = Inf
        mes = 'LD: '+ str( '%.3f' %  (LD)) + ' m' ,
    LNL= 1/((new_gamma*1E-3)*new_p0)             
    return ['\u03B22: '+ str( '%.2f' %  (new_beta2)) + ' ps',
            '\u03B23: '+ str( '%.2f' %  (new_beta3)) + ' ps',
            '\u03B3: '+str('%.3f' % (new_gamma))+ ' [1/(W km)]',
            'P0: '+str('%.3f' % (new_p0))+r' W',
            mes,  
            'LNL: '+ str( '%.2f' %  (LNL)) + ' m',
            'N: '+ str( '%.2f' %  (np.sqrt(LD/LNL)))
            ]

@app.callback( 
    [Output('alpha_v', 'children'),
    Output('m_val', 'children'),
    Output('c_val', 'children'),
    Output('L_', 'children'),
    ],
    [Input('alpha_slider', 'value'),
    Input('m_slid', 'value'),
    Input('c_slid', 'value'),
    Input('L_slid', 'value'),
    ]
    )
# function to update with the callback:
def update_left(new_alpha, new_m, new_c, new_zmax):                     
    return ['\u03B1: '+str('%.3f' % (new_alpha))+ ' [dB/km]',
            'm: '+ str((new_m)),
            'C: '+ str((new_c)),
            'L: '+str('%.2f' % (new_zmax))+r'm',
            ]

@app.callback( 
    [
    Output('loading1', 'children'),
    Output('loading2', 'children'),
    Output('z_slid', 'marks'),
    Output('z_slid', 'value'),
    Output('session', 'data'),
    ],
    [Input('calculate', 'n_clicks'),],
    [
    State('alpha_slider', 'value'),
    State('beta2_slider', 'value'),
    State('beta3_slid', 'value'),
    State('gamma_slider', 'value'),
    State('P0_slid', 'value'),
    State('m_slid', 'value'),
    State('c_slid', 'value'),
    State('L_slid', 'value'),
    State('switch', 'value'),
    State('session', 'data'),
    ]
    )
# function to update with the callback:
def update_plots(n_clicks,new_alpha,new_beta2, new_beta3, new_gamma, new_p0, new_m, new_c, new_L, switch, data):
    if n_clicks is None:
        raise dash.exceptions.PreventUpdate
    else:
        pulsetype = 'Sech' if switch else 'Gaussian'
        new_beta2 *= 1E-27
        new_beta3 *= 1E-39
        new_gamma *= 1E-3
#T0, T, solve_type= 'incident_field', L=0.1, beta2=0, gamma=0, P0=0,  beta3=0, loss = 0, pulsetype = 'Gaussian', m = 1, C=0
        #global pulse
        #if new_beta2 == 0:
        pulse = Propagation(T0, T, m = new_m, 
                        C=new_c, pulsetype = pulsetype,
                        solve_type='split_step', 
                        L=new_L, 
                        beta2=new_beta2,
                        beta3=new_beta3,
                        gamma=new_gamma, 
                        P0=new_p0,
                        loss = new_alpha,
                        size_array = size_array)
        # else:
        #     pulse = Propagation(T0, T, m = new_m, 
        #                     C=new_c, pulsetype = pulsetype,
        #                     solve_type='split_step', 
        #                     L=new_L, 
        #                     beta2=new_beta2,
        #                     gamma=new_gamma, 
        #                     P0=new_p0,
        #                     loss = new_alpha,
        #                     size_array = size_array)
        UI = pulse.UI
        UIW = pulse.UIW
        Z = pulse.z
        data = {'ui': UI, 'uiw': UIW, 'zi': Z}
        marks={
        0: {'label': '0', 'style': {'color': colors['text']}},
        lastz: {'label': '{0} m'.format( '%.2f' % (pulse.z[-1])), 'style': {'color': colors['text']}}}
        prop_t = pulse.plot_propagation(mode = 'time')
        prop_w = pulse.plot_propagation(mode = 'spectrum')
        return [dcc.Graph(figure=prop_t.update_layout()), 
                dcc.Graph(figure=prop_w.update_layout()), 
                marks,0,data]

@app.callback( 
    [Output('envelopet', 'figure'),
    Output('envelopew', 'figure'),
    Output('z_', 'children'),
    ],
    [Input('z_slid', 'value'),
    State('session', 'data'),
    State('switch', 'value'),
    ])

def update_envelope(new_z,data,switch):
    pulsetype = 'Sech' if switch else 'Gaussian'
    Zi = data.get('zi', Z)
    return [plot_envelope(T/T0,data.get('ui', UI), pulsetype,colors,mode = 'time', z0 = new_z),
            plot_envelope(W,data.get('uiw', UIW), pulsetype,colors,mode = 'spectrum', z0 = new_z),
            #'z: '+str('%.2f' % (pulse.z[new_z]))+r'm',
            'z: '+str('%.2f' % (Zi[new_z]))+r'm',
            ]

@app.callback(Output("download-text", "data"), Input("download", "n_clicks"), 
    [State('alpha_slider', 'value'),
    State('beta2_slider', 'value'),
    State('beta3_slid', 'value'),
    State('gamma_slider', 'value'),
    State('P0_slid', 'value'),
    State('m_slid', 'value'),
    State('c_slid', 'value'),
    State('L_slid', 'value'),
    State('z_slid', 'value'),
    State('switch', 'label'),
    State('session', 'data'),]) #save local->download the pandas csv

def func(n_clicks, new_alpha, new_beta2, new_beta3, new_gamma, new_p0, new_m, new_C, new_L, new_z0, pulse, data):
    if n_clicks is None:
        raise dash.exceptions.PreventUpdate
    else:
        Zi = data.get('zi', Z)
        return dict(content=
        '''
    Loss = {0} dB/km
    Dispersion = {1} ps²/km
    Third-order Dispersion = {2} ps³/km
    Gamma = {3} 1/(W km)
    Peak power = {4} W
    m = {5}
    Initial chirp = {6}
    Fiber length = {7} m
    Current z for envelope = {8} m
    for a {9}
    with T0 = {10} s
        '''.format(str('%.3f' % (new_alpha)), str('%.3f' % (new_beta2)),str('%.3f' % (new_beta3)),str('%.3f' % (new_gamma)), str('%.3f' % (new_p0)), new_m, new_C, new_L, Zi[new_z0], pulse, T0), 
        filename="parameters.txt")
