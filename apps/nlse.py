# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from numpy.core.numeric import Inf
import plotly.graph_objects as go   
import numpy as np
from Functions import split_step, plot_prop

#------------------Definitons--------------#


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

last_update = 'Last update: 14.06.2021'

#------------------- Grid: --------------------#
#Input pulse
#Tmax = 5 ~ 10ps
T0 = 5E-12 #  duration of input
#for pulse width  --> Dispersive effects at T0 ~ 1ps 
#replaced by usr

T_FWHM_g = 2*np.sqrt(np.log(2))*T0

#For Sech pulses:  
#Chirp parameter
C = 0 
T_FWHM_sech = 2*np.log(1+np.sqrt(2))*T0

N = 8196 #ammount of points 
dt = 100*T0/N #the 16 is to get a grid between -8 and 8 for T/T0 
#dt = 8*T_FWHM/N

T = np.arange(-N/2, N/2)*dt


# T_selected_g = T_FWHM_g##### Update with a button(?)
# T_selected_s = T_FWHM_sech

T_selected_g = T0##### Update with a button(?)
T_selected_s = T0

z_initial = 0 #km
zmax = 2.5 # km
#gamma = n2*wo/(speed*Aeff)
beta2_initial = 5.66099
beta3_initial = 10#ps^3/km
gamma_initial = 2 #1/(W*km)
P0 = 10E-3
m0 = 3
C = 0



# def create_z_slider(z, colors):
#     z_slider = dcc.Slider(
#         min=z[0],
#         max=z[-1],
#         step=z[1]-z[0],
#         value=z[0],
#         marks={
#         z[0]: {'label': '0', 'style': {'color': colors['text']}},
#         z[-1]: {'label': '{0} km'.format(z[-1]), 'style': {'color': colors['text']}}},
#     )
#     return z_slider

def create_z_slider(z, colors):
    last = len(z)-1
    z_slider = dcc.Slider(
        min=0,
        max=last,
        step=1,
        value=0,
        marks={
        0: {'label': '0', 'style': {'color': colors['text']}},
        last: {'label': '{0} km'.format(z[-1]), 'style': {'color': colors['text']}}},
    )
    return z_slider

#-------------------------------------------#

#-------- SLIDERS---------------

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
        max=99,
        step=1,
        value=z_initial,
        marks={
        0: {'label': '0', 'style': {'color': colors['text']}},
        99: {'label': '{0} km'.format(zmax), 'style': {'color': colors['text']}}},
    )

L_slider = dcc.Slider(
        id='L_slid',
        min=0.001,
        max=5,
        step=0.001,
        value=zmax,
        marks={
        0: {'label': '0', 'style': {'color': colors['text']}},
        5: {'label': '{0} km'.format(5), 'style': {'color': colors['text']}}},
    )

gamma_slider = dcc.Slider(
        id='gamma_slider',
        min=0.01,
        max=2.0,
        step=0.01,
        value=gamma_initial,
        marks={
        0.01: {'label': '0.01', 'style': {'color': colors['text']}},
        2: {'label': '2', 'style': {'color': colors['text']}}},
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
        id='mg_slid',
        min=2,
        max=10,
        step=1,
        value=m0,
        marks={
        2: {'label': '2', 'style': {'color': colors['text']}},
        10: {'label': '10', 'style': {'color': colors['text']}}},
    )

#-----------------------------------------------------------------------#
beta2_initial *= 1E-24
beta3_initial *= 1E-36
LD_g = ((T_selected_g)**2)/np.absolute(beta2_initial*1E-24)  # Einheit!!!: km #HERE 1E-24 due to the square of ps^2
LD_s = ((T_selected_s)**2)/np.absolute(beta2_initial*1E-24)  # Einheit!!!: km #HERE 1E-24 due to the square of ps^2
LNL= 1/(gamma_initial*P0)

if beta2_initial == 0:
    print('N´ = sqrt(L´D/LNL) = : ',np.sqrt((gamma_initial*P0*T_selected_g**3/np.absolute(beta3_initial))))
    U, UI, UW, UIW, W, z = split_step(beta2_initial, T_selected_g, T, zmax, gamma_initial, P0, beta3_initial,pulse = 'Gaussian', m = 1, C=0)
    Us, UIs, UWs, UIWs, Ws, zs = split_step(beta2_initial, T_selected_s, T, zmax, gamma_initial, P0, beta3_initial,pulse = 'Sech', m = 1, C=0)
else:
    LD = (T_selected_g**2)/np.absolute(beta2_initial)
    U, UI, UW, UIW, W, z = split_step(beta2_initial, T_selected_g, T, zmax, gamma_initial, P0 ,pulse = 'Gaussian', m = 1, C=0)
    Us, UIs, UWs, UIWs, Ws, zs = split_step(beta2_initial, T_selected_s, T, zmax, gamma_initial, P0 ,pulse = 'Sech', m = 1, C=0)

#PLOTS:
Utg = go.Scatter(x=T/T_selected_g ,y=UI[0], name = 'Gaussian',
                         line=dict(color=colors['even']))

Uts = go.Scatter(x=T/T_selected_s ,y=UIs[0], name = 'Hyperbolic-Secant',
                         line=dict(color=colors['odd']))
nlset_fig = go.Figure(data=[Utg,Uts]).update_layout(
#      
    updatemenus = list([
        dict(
            type="buttons",
            active=0,
            buttons=list([   
                dict(label = 'Gaussian',
                    method = 'update',
                    #args = [{'visible': [True, True, True, False, False, False]},
                    args = [{'visible': [True, False]},
                            {'title': 'NLSE: Gaussian Pulse.'}]), # using Eq. 3.2.5 and 3.2.6

                dict(label = 'Sech',
                    method = 'update',
                    args = [{'visible': [False, True]},
                            {'title': '''
                            NLSE: Sech Pulse.'''}])  #using Eq. 3.2.5 and 3.2.6
            ]),
        )
    ])
)
nlset_fig.update_layout( 
                        width=600, height=600,
                        plot_bgcolor  = colors['background'],
                        paper_bgcolor = colors['background'],
                        font= {
                                'color': colors['text']},
                        yaxis=dict(range=[0, 1.1],title='|U(z,T)|^2', 
                                    ), 
                        xaxis=dict(range=[-8, 8],title='T/T0', 
                                    ), 
                        
                        )


nlset_graph = dcc.Graph(id='nlset_plot',
                        animate=True,
                        figure=nlset_fig.update_layout(

))


Ufg = go.Scatter(x=W ,y=UIW[0], name = 'Gaussian',
                         line=dict(color=colors['even']))

Ufs = go.Scatter(x=Ws ,y=UIWs[0], name = ' Hyperbolic-Secant',
                         line=dict(color=colors['odd']))

nlsef_fig = go.Figure(data=[Ufg,Ufs]).update_layout(
#      
    updatemenus = list([
        dict(
            type="buttons",
            active=0,
            buttons=list([   
                dict(label = 'Gaussian',
                    method = 'update',
                    #args = [{'visible': [True, True, True, False, False, False]},
                    args = [{'visible': [True, False]},
                            {'title': 'NLSE freq: Gaussian Pulse.'}]), # using Eq. 3.2.5 and 3.2.6

                dict(label = 'Sech',
                    method = 'update',
                    args = [{'visible': [False, True]},
                            {'title': '''
                            NLSE freq: Sech Pulse.'''}])  #using Eq. 3.2.5 and 3.2.6
            ]),
        )
    ])
)

nlsef_fig.update_layout( 
                        width=600, height=600,
                        plot_bgcolor  = colors['background'],
                        paper_bgcolor = colors['background'],
                        font= {
                                'color': colors['text']},
                        yaxis=dict(title='|U(z,\u03C9)|^2', 
                                    ), 
                        xaxis=dict(range=[-2E12, 2E12],title='\u03C9-\u03C90', #rangemode='nonzero',
                                    ),
                        )


nlsef_graph = dcc.Graph(id='nlsef_plot',
                        animate=True,
                        figure=nlsef_fig.update_layout(

))


prop_plot_g = plot_prop(UI, T/T_selected_g,z)
prop_plot_s = plot_prop(UIs, T/T_selected_s,zs)

propagation_g = dcc.Graph(id='propagation_gauss',
                            figure=prop_plot_g)

propagation_s = dcc.Graph(id='propagation_sech',
                            figure=prop_plot_s)

#-------Final Layout for plotting and display of sliders and constants-------#
layout = html.Div(style={'backgroundColor': colors['background']},
    children=[
        html.Div(id='dpf-display-value'),
        dcc.Link('Go to SPM effect', href='/apps/spm'), ##EDIT LINK TO OTHER PAGES
        html.Div(id='modesf-display-value'),
        dcc.Link('Go to GVD effect', href='/apps/gvd'), ##EDIT LINK TO OTHER PAGES
        
        html.Div(className='rows', 
        children=[
            html.Div(className='four columns div-user-controls', style={'backgroundColor': colors['background']}, 
            children = [
                html.P('Front-End for NLSE', style={'color': colors['text']}),
                html.H3(children=[html.Var('\u03B22: '+str( '%.3f' % (beta2_initial))+ ' ps', id='beta2_v'),
                html.Sup(2), html.Var('/km')],
                style={'color': colors['text']}),
                beta2_slider,
                html.H3(children=[html.Var('\u03B23: '+str('%.3f' % (beta3_initial))+ ' ', id='beta3_v'),
                html.Sup(2), html.Var(' ')],
                style={'color': colors['text']}),
                beta3_slider,
                html.H3('\u03B3: '+str(gamma_initial)+ ' [m/W]', id='gamma_v',style={'color': colors['text']}),
                gamma_slider,
                html.H3('P0: '+str(P0)+r'W', id='p0_v', style={'color': colors['text']}),
                Power_slider,
                html.H3('m: '+ str(m0), id = 'm_val'),
                m_slider,
                html.H3('LD Gauss: '+ str( '%.3f' %  (LD_g)) + ' km', id = 'LD_g'),
                html.H3('LD Sech: '+ str( '%.3f' %  (LD_s)) + ' km', id = 'LD_s'),
                html.H3('LNL: '+ str( '%.3f' %  (LNL)) + ' km', id = 'LNL_v'),
                html.H3('L: '+str(zmax)+r'km', id='L_', style={'color': colors['text']}),
                L_slider,
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
                    html.H3('NLSE Colormap time (Gauss)'),                 
                    propagation_g,
                    html.H3('NLSE Colormap time (Sech)'),                 
                    propagation_s,
                    html.H3('NLSE time'),
                    nlset_graph,
                    html.H3('z: '+str(z[z_initial])+r'km', id='z_', style={'color': colors['text']}),
                    z_slider,
                    html.H3('NLSE freq'),
                    nlsef_graph,
                    
                              
            ])  
        ]),
        html.Footer('Joly Nicolas, Aguilera Hernan. Max-Planck-Institut', style={'color': colors['text']}),
        html.Footer(last_update, style={'color': colors['text']}),
    ])

    
@app.callback( 
    [Output('beta2_v', 'children'),
    Output('beta3_v', 'children'),
    Output('gamma_v', 'children'),
    Output('p0_v', 'children'),
    Output('m_val', 'children'),
    Output('LD_g', 'children'),
    Output('LD_s', 'children'),
    Output('LNL_v', 'children'),
    #Output('z_', 'children'),
    Output('L_', 'children'),
    ],
    [Input('beta2_slider', 'value'),
    Input('beta3_slid', 'value'),
    #Input('z_slid', 'value'),
    Input('L_slid', 'value'),
    Input('gamma_slider', 'value'),
    Input('P0_slid', 'value'),
    Input('mg_slid', 'value'),
    ]
    )

# function to update with the callback:
def update_left(new_beta2, new_beta3, new_zmax, new_gamma, new_p0, new_m):
    LD_g = ((T_selected_g)**2)/np.absolute(new_beta2*1E-24)  # Einheit!!!: km #HERE 1E-24 due to the square of ps^2
    LD_s = ((T_selected_s)**2)/np.absolute(new_beta2*1E-24)  # Einheit!!!: km #HERE 1E-24 due to the square of ps^2
    
    LNL_n= 1/(new_gamma*new_p0)

                        
    return ['\u03B22: '+ str( '%.3f' %  (new_beta2)) + ' ps',
            '\u03B23: '+ str( '%.3f' %  (new_beta3)) + ' ps',
            '\u03B3: '+str(new_gamma)+ ' [m/W]',
            'P0: '+str(new_p0)+r' W',
            'm: '+ str(new_m),
            'LD Gauss: '+ str( '%.3f' %  (LD_g)) + ' km',
            'LD Sech: '+ str( '%.3f' %  (LD_s)) + ' km',
            'LNL: '+ str( '%.3f' %  (LNL_n)) + ' km',
            #'z: '+str(z[new_z])+r'km',
            'L: '+str(new_zmax)+r'km',
            ]


@app.callback( 
    [
    Output('propagation_gauss', 'figure'),
    Output('propagation_sech', 'figure'),
    Output('z_slid', 'marks'),
    Output('z_slid', 'value')
    ],
    [Input('calculate', 'n_clicks'),
    Input('beta2_slider', 'value'),
    Input('beta3_slid', 'value'),
    Input('L_slid', 'value'),
    Input('gamma_slider', 'value'),
    Input('P0_slid', 'value'),
    Input('mg_slid', 'value'),
    ],
    [
    State(component_id='propagation_gauss', component_property='figure'),
    State(component_id='propagation_sech', component_property='figure'),
    State(component_id='z_slid', component_property='marks'),
    State(component_id='z_slid', component_property='value'),]
    )

# function to update with the callback:
def update_plots(btn1,new_beta2, new_beta3, new_zmax, new_gamma, new_p0, new_m, plot1, plot2, slidd, vaal): 
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    new_beta2 *= 1E-24
    new_beta3 *= 1E-36
    global UI, UIW, W, z, UIs, UIWs, Ws, zs
    if 'calculate' in changed_id:
        if beta2_initial == 0:
            _, UI, _, UIW, W, z = split_step(new_beta2, T_selected_g, T, new_zmax, new_gamma, new_p0, new_beta3, pulse = 'Gaussian', m = new_m, C=0)
            _, UIs, _, UIWs, Ws, zs = split_step(new_beta2, T_selected_g, T, new_zmax, new_gamma, new_p0, new_beta3,pulse = 'Sech', C=0)
        else:
            _, UI, _, UIW, W, z = split_step(new_beta2, T_selected_g, T, new_zmax, new_gamma, new_p0 ,pulse = 'Gaussian', m = 1, C=0)
            _, UIs, _, UIWs, Ws, zs = split_step(new_beta2, T_selected_g, T, new_zmax, new_gamma, new_p0,pulse = 'Sech', m = 1, C=0)
        
        marks={
        0: {'label': '0', 'style': {'color': colors['text']}},
        99: {'label': '{0} km'.format( '%.3f' % (z[-1])), 'style': {'color': colors['text']}}}

        
        prop_plot_g = plot_prop(UI, T/T_selected_g,z)
        prop_plot_s = plot_prop(UIs, T/T_selected_s,zs)
        
        return [prop_plot_g,
                prop_plot_s,
                marks,0]
    else: return [plot1, plot2, slidd,vaal]
    

@app.callback( 
    [Output('nlset_plot', 'figure'),
    Output('nlsef_plot', 'figure'),
    Output('z_', 'children'),
    ],
    [Input('z_slid', 'value'),
    ])

def update_envelope(new_z):
    Utg = go.Scatter(x=T/T_selected_g ,y=UI[new_z], name = 'Gaussian',
                        line=dict(color=colors['even']))

    Uts = go.Scatter(x=T/T_selected_s ,y=UIs[new_z], name = 'Hyperbolic-Secant',
                        line=dict(color=colors['odd']))
    Ufg = go.Scatter(x=W ,y=UIW[new_z], name = 'Gaussian',
                        line=dict(color=colors['even']))

    Ufs = go.Scatter(x=Ws ,y=UIWs[new_z], name = ' Hyperbolic-Secant',
                        line=dict(color=colors['odd']))
    return [
            go.Figure(data=[Utg,Uts]).update_layout(),
            go.Figure(data=[Ufg,Ufs]).update_layout(),
            'z: '+str(z[new_z])+r'km',
            ]