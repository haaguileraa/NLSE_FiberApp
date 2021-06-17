# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from numpy.core.numeric import Inf
import plotly.graph_objects as go   
import numpy as np
from Functions import incident_field_spm


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

#------------------Definitons--------------#

#------------------- Grid: --------------------#
#Input pulse
#Tmax = 5 ~ 10ps
T0 = 200E-12 # for pulse Initial width T0 --> Dispersive effects at T0 ~ 1ps  S.64
# half-width (at 1/e-intensity point)

N = 8196 #ammount of points 
dt = 100*T0/N #the 16 is to get a grid between -8 and 8 for T/T0 
#dt = 8*T_FWHM/N


#T = np.linspace(-8*T0,8*T0,N)
T = np.arange(-N/2, N/2)*dt


T_FWHM_g = 2*np.sqrt(np.log(2))*T0
T_FWHM_sech = 2*np.log(1+np.sqrt(2))*T0


T_selected_g = T0##### Update with a button(?)

T_selected_s = T0##### Update with a button(?)


#gamma = n2*wo/(speed*Aeff)
gamma_initial = 1 #1/(W*km)
P0 = 10E-3
LNL= 1/(gamma_initial*P0)
Leff = LNL
m0 = 3
C = 0
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


m_slider = dcc.Slider(
        id='m_slid',
        min=2,
        max=10,
        step=1,
        value=m0,
        marks={
        2: {'label': '2', 'style': {'color': colors['text']}},
        10: {'label': '10', 'style': {'color': colors['text']}}},
    )









#---------------------------------------------------------------#
#///////////////////////////////////////////////////////////////#
#-------------------                        --------------------#
Phi, dw = incident_field_spm(Leff, gamma_initial, P0, T_selected_g, T, pulse = 'Gaussian', m = 1, C=0)
Phi2, dw2 = incident_field_spm(Leff, gamma_initial, P0, T_selected_g, T, pulse = 'Gaussian', m = m0, C=0)
Phi_s, dw_s = incident_field_spm(Leff, gamma_initial, P0, T_selected_s, T, pulse = 'Sech', C=0)

gauss_p1 = go.Scatter(x=T/T_selected_g ,y=Phi, name = 'Gaussian',
                         line=dict(color=colors['even']))


gauss_p2 = go.Scatter(x=T/T_selected_g ,y=Phi2, name = 'Super-Gaussian',
                         line=dict(color=colors['odd']))

gauss_chirp1 = go.Scatter(x=T/T_selected_g ,y=dw*T_selected_g, name = 'Gaussian',
                         line=dict(color=colors['even']))


gauss_chirp2 = go.Scatter(x=T/T_selected_g ,y=dw2*T_selected_g, name = 'Super-Gaussian',
                         line=dict(color=colors['odd']))


sech_p = go.Scatter(x=T/T_selected_s ,y=Phi_s, name = 'Sech',
                         line=dict(color=colors['other']))

sech_chirp = go.Scatter(x=T/T_selected_s ,y=dw_s*T_selected_s, name = 'Sech',
                         line=dict(color=colors['other']))


spm_phase = go.Figure(data=[gauss_p1,gauss_p2, sech_p ]).update_layout( 
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
                        width=600, height=600,
                        plot_bgcolor  = colors['background'],
                        paper_bgcolor = colors['background'],
                        font= {
                                'color': colors['text']},
                        yaxis=dict(range=[0, 1.1],title='Phase \u03C6NL', 
                                    ), 
                        xaxis=dict(range=[-2.5, 2.5],title='T/T0', 
                                    ), 
                        )


phase_graph = dcc.Graph(id='phase_spm_plot',
                        animate=True,
                        figure=spm_phase.update_layout(

))



spm_chirp = go.Figure(data=[gauss_chirp1,gauss_chirp2, sech_chirp ]).update_layout( 
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
                        width=600, height=600,
                        plot_bgcolor  = colors['background'],
                        paper_bgcolor = colors['background'],
                        font= {
                                'color': colors['text']},
                        yaxis=dict(range=[-3, 3],title='frequency chirp \u03B4\u03C9T0', 
                                    ), 
                        xaxis=dict(range=[-2.5, 2.5],title='T/T0', 
                                    ), 
                        )


chirp_graph = dcc.Graph(id='chirp_spm_plot',
                        animate=True,
                        figure=spm_chirp.update_layout(

))

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
                html.H3('m: '+ str(m0), id = 'm_display'),
                m_slider,
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
    Output('phase_spm_plot', 'figure'),
    Output('chirp_spm_plot', 'figure'),
    ],
    [Input('gamma_slid', 'value'),
    Input('PP_slid', 'value'),
    Input('m_slid', 'value')
    ]
    )

# function to update with the callback:
def update_plot(new_gamma, new_pp, mn):
    
    if new_gamma != 0 and new_pp != 0 :
        LNL_n = 1/(new_gamma*new_pp)
    else:
        LNL_n = Inf
    Leff_n = LNL_n
    Phi_n, dw_n = incident_field_spm(Leff_n, new_gamma, new_pp, T_selected_g, T, pulse = 'Gaussian', m = 1, C=0)
    Phi2_n, dw2_n = incident_field_spm(Leff_n, new_gamma, new_pp, T_selected_g, T, pulse = 'Gaussian', m = mn, C=0)
    Phi_s_n, dw_s_n = incident_field_spm(Leff_n, new_gamma, new_pp, T_selected_s, T, pulse = 'Sech', C=0)

    gauss_p1_n = go.Scatter(x=T/T_selected_g ,y=Phi_n, name = 'Gaussian',
                            line=dict(color=colors['even']))


    gauss_p2_n = go.Scatter(x=T/T_selected_g ,y=Phi2_n, name = 'Super-Gaussian',
                            line=dict(color=colors['odd']))

    gauss_chirp1_n = go.Scatter(x=T/T_selected_g ,y=dw_n*T_selected_g, name = 'Gaussian',
                            line=dict(color=colors['even']))


    gauss_chirp2_n = go.Scatter(x=T/T_selected_g ,y=dw2_n*T_selected_g, name = 'Super-Gaussian',
                            line=dict(color=colors['odd']))


    sech_p_n = go.Scatter(x=T/T_selected_s ,y=Phi_s_n, name = 'Sech',
                            line=dict(color=colors['other']))

    sech_chirp_n = go.Scatter(x=T/T_selected_s ,y=dw_s_n*T_selected_s, name = 'Sech',
                            line=dict(color=colors['other']))

    spm_chirp_n = go.Figure(data=[gauss_chirp1_n,gauss_chirp2_n, sech_chirp_n ]).update_layout()
    spm_phase_n = go.Figure(data=[gauss_p1_n,gauss_p2_n, sech_p_n ]).update_layout()


    return ('\u03B3: '+str(new_gamma)+ ' [1/W/m]',
            'P0: '+str(new_pp)+r' W',
            'LNL: '+ str( '%.3f' %  (LNL_n)) + ' m',
            'm: '+ str(mn),
            spm_phase_n,
            spm_chirp_n
            )