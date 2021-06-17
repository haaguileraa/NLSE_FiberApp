# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from numpy.core.numeric import Inf
import plotly.graph_objects as go   
import numpy as np
from Functions import Gaussian_pulse_GVD, incident_field


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
T0 = 1E-12 #  duration of input
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


T_selected_g = T_FWHM_g##### Update with a button(?)
T_selected_s = T_FWHM_sech


beta2_initial = 5.66099
L_initial = 0

LD = ((T0)**2)/np.absolute(beta2_initial*1E-24)  # Einheit!!!: km #HERE 1E-24 due to the square of ps^2
C = 0
#-------------------------------------------#


z = 0
z1 = 2*LD
z2 = 4*LD



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


Lange_slider = dcc.Slider(
        id='lange_slid',
        min=0,
        max=2.5,
        step=0.001,
        value=L_initial,
        marks={
        0: {'label': '0', 'style': {'color': colors['text']}},
        2.5: {'label': '2.5 km', 'style': {'color': colors['text']}}},
    )


beta2_initial *= 1E-24








#-----------------------Eq. 3.2.7 and 3.2.9---------------------#


_, UI1 = Gaussian_pulse_GVD(z,T,T_selected_g, beta2_initial)
_, UI2 = Gaussian_pulse_GVD(z1,T,T_selected_g, beta2_initial)
_, UI3 = Gaussian_pulse_GVD(z2,T,T_selected_g, beta2_initial)

#print(UI1)
Uz0 = go.Scatter(x=T/T0 ,y=UI1, name = 'z = 0',
                         line=dict(color=colors['even']))
U2LD = go.Scatter(x=T/T0 ,y=UI2, name = 'z = 2*LD',
                         line=dict(color=colors['odd']))
U4LD = go.Scatter(x=T/T0 ,y=UI3, name = 'z = 4*LD',
                         line=dict(color=colors['other']))

gvd_z = go.Figure(data=[Uz0, U2LD, U4LD])

plot_gvd = dcc.Graph(id='gauss',
          animate=True,
          figure=gvd_z.update_layout( 
                            width=600, height=600,
                            plot_bgcolor  = colors['background'],
                            paper_bgcolor = colors['background'],
                            font= {
                                    'color': colors['text']},
                            yaxis=dict(range=[0, 1.1],title='|U(z,T)|^2', 
                                        # gridwidth=0.8, gridcolor='#121212', 
                                        # zeroline=True,zerolinewidth=0.8, 
                                        # zerolinecolor='#121212'
                                        ), 
                            xaxis=dict(range=[-8, 8],title='T/T0', 
                                        # gridwidth=0.8, gridcolor='#121212',
                                        # zeroline=True,zerolinewidth=0.8, 
                                        # zerolinecolor='#121212'
                                        ), 
                            title= '''
                                    Dispersion-induced broadening of a Gaussian pulse.
                                    ''',
                            )
          )
#---------------------------------------------------------------#
#///////////////////////////////////////////////////////////////#
#-------------------Using Eq 3.2.5 and 3.2.6--------------------#
_, UI4, UW4, W4 = incident_field(beta2_initial, z, T_selected_g, T, pulse = 'Gaussian')
# _, UI5, UW5, W5 = incident_field(beta2_initial, z1, T_selected_g, T, pulse = 'Gaussian')
# _, UI6, UW6, W6 = incident_field(beta2_initial, z2, T_selected_g, T, pulse = 'Gaussian')
_, UI7, UW7, W7 = incident_field(beta2_initial, z, T_selected_s, T, pulse = 'Sech', C = C)
# _, UI8, UW8, W8 = incident_field(beta2_initial, z1, T_selected_s, T, pulse = 'Sech', C = C)
# _, UI9, UW9, W9 = incident_field(beta2_initial, z2, T_selected_s, T, pulse = 'Sech', C = C)


Uz4 = go.Scatter(x=T/T_selected_g ,y=UI4, name = 'Gaussian',
                         line=dict(color=colors['even']))
# Uz5 = go.Scatter(x=T/T_selected_g ,y=UI5, name = 'z = 2*LD',
#                          line=dict(color=colors['odd']))
# Uz6 = go.Scatter(x=T/T_selected_g ,y=UI6, name = 'z = 4*LD',
#                          line=dict(color=colors['other']))


Uz7 = go.Scatter(x=T/T_selected_s ,y=UI7, name = 'Hyperbolic-Secant',
                         line=dict(color=colors['odd']))
# Uz8 = go.Scatter(x=T/T_selected_s ,y=UI8, name = 'z = 2*LD',
#                          line=dict(color=colors['odd']))
# Uz9 = go.Scatter(x=T/T_selected_s ,y=UI9, name = 'z = 4*LD',
#                          line=dict(color=colors['other']))




#gvd_fig = go.Figure(data = [Uz4, Uz5, Uz6, Uz7, Uz8, Uz9]).update_layout(     
gvd_fig = go.Figure(data=[Uz4,Uz7]).update_layout(
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
                            {'title': 'Dispersion-induced broadening: Gaussian Pulse.'}]), # using Eq. 3.2.5 and 3.2.6

                dict(label = 'Sech',
                    method = 'update',
                    args = [{'visible': [False, True]},
                            {'title': '''
                            Dispersion-induced broadening: Sech Pulse.'''}])  #using Eq. 3.2.5 and 3.2.6
            ]),
        )
    ])
)

#gvd_fig.update(data = [Uz4, Uz5, Uz6, Uz7, Uz8, Uz9])
#gvd_z = go.Figure(data=[Uz4, Uz5, Uz6, Uz7, Uz8, Uz9])

gvd_fig.update_layout( 
                        width=600, height=600,
                        plot_bgcolor  = colors['background'],
                        paper_bgcolor = colors['background'],
                        font= {
                                'color': colors['text']},
                        yaxis=dict(range=[0, 1.1],title='|U(z,T)|^2', 
                                    # gridwidth=0.8, gridcolor='#121212', 
                                    # zeroline=True,zerolinewidth=0.8, 
                                    # zerolinecolor='#121212'
                                    ), 
                        xaxis=dict(range=[-5, 5],title='T/T0', 
                        #xaxis=dict(range=[T[0]/T_selected_g-1, T[-1]/T_selected_g+1],title='T/T0', 
                                    # gridwidth=0.8, gridcolor='#121212',
                                    # zeroline=True,zerolinewidth=0.8, 
                                    # zerolinecolor='#121212'
                                    ), 
                        
                        )


gvd_graph = dcc.Graph(id='sech_gauss_plot',
                        animate=True,
                        figure=gvd_fig.update_layout(

))


UW_4 = np.absolute(UW4)**2
UW_7 = np.absolute(UW7)**2

Uf4 = go.Scatter(x=W4 ,y=UW_4, name = 'Gaussian',
                         line=dict(color=colors['even']))

Uf7 = go.Scatter(x=W7 ,y=UW_7, name = ' Hyperbolic-Secant',
                         line=dict(color=colors['odd']))

gvd_fig_frec = go.Figure(data=[Uf4,Uf7]).update_layout(
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
                            {'title': 'Dispersion-induced broadening: Gaussian Pulse.'}]), # using Eq. 3.2.5 and 3.2.6

                dict(label = 'Sech',
                    method = 'update',
                    args = [{'visible': [False, True]},
                            {'title': '''
                            Dispersion-induced broadening: Sech Pulse.'''}])  #using Eq. 3.2.5 and 3.2.6
            ]),
        )
    ])
)

gvd_fig_frec.update_layout( 
                        width=600, height=600,
                        plot_bgcolor  = colors['background'],
                        paper_bgcolor = colors['background'],
                        font= {
                                'color': colors['text']},
                        yaxis=dict(title='|U(z,\u03C9)|^2', 
                                    ), 
                        xaxis=dict(title='\u03C9-\u03C90', #rangemode='nonzero',
                                    ),
                        )


gvd_graph_frec = dcc.Graph(id='sech_gauss_plot_f',
                        animate=True,
                        figure=gvd_fig_frec.update_layout(

))

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
                html.H3('z: '+str(L_initial)+r'km', id='z_val', style={'color': colors['text']}),
                Lange_slider,
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
                    html.H3('Using Eq. 3.2.5 and 3.2.6 Agrawal'),
                    gvd_graph,
                    gvd_graph_frec,
                    html.H3('Using Eq. 3.2.7 and 3.2.9 Agrawal'),
                    plot_gvd,
                    
                    
            ])  
        ]),
        html.Footer('Joly Nicolas, Aguilera Hernan. Max-Planck-Institut', style={'color': colors['text']}),
        html.Footer(last_update, style={'color': colors['text']}),
    ])

@app.callback( 
    [Output('beta2_val', 'children'),
    Output('z_val', 'children'),
    Output('LD_display', 'children'),
    Output('sech_gauss_plot', 'figure'),
    Output('sech_gauss_plot_f', 'figure'),
    ],
    [Input('beta2_slid', 'value'),
    Input('lange_slid', 'value'),]
    )

# function to update with the callback:
def update_plot(new_beta2, new_z):
    new_beta2 *= 1E-24
    if new_beta2 != 0:
        LD_n = ((T0)**2)/np.absolute(new_beta2)
    else:
        LD_n = Inf
        #-------------------Using Eq 3.2.5 and 3.2.6--------------------#
    _, UI4_n, UW4_n, W4_n = incident_field(new_beta2, new_z, T_selected_g, T, pulse = 'Gaussian')#*1E-24
    _, UI7_n, UW7_n, W7_n = incident_field(new_beta2, new_z, T_selected_s, T, pulse = 'Sech', C = C)


    Uz4_n = go.Scatter(x=T/T_selected_g ,y=UI4_n, name = 'Gaussian z = '+str(new_z),
                            line=dict(color=colors['even']))

    Uz7_n = go.Scatter(x=T/T_selected_s ,y=UI7_n, name = 'Sech z = '+str(new_z),
                            line=dict(color=colors['odd']))

    new_sec_Gauss = go.Figure(data=[Uz4_n, Uz7_n]).update_layout()

    UW_4_n = np.absolute(UW4_n)**2
    UW_7_n = np.absolute(UW7_n)**2
    Uf4_n = go.Scatter(x=W4_n ,y=UW_4_n, name = 'Gaussian',
                         line=dict(color=colors['even']))

    Uf7_n = go.Scatter(x=W7_n ,y=UW_7_n, name = ' Hyperbolic-Secant',
                         line=dict(color=colors['odd']))

    new_sec_Gauss_freq = go.Figure(data=[Uf4_n, Uf7_n]).update_layout()



    return ('\u03B22: '+ str( '%.3f' %  (new_beta2*1E+24)) + ' ps',
            'z: '+str(new_z)+r'km',
            'LD: '+ str( '%.3f' %  (LD_n)) + ' km',
            new_sec_Gauss,
            new_sec_Gauss_freq
            )
