# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objects as go   

from Functions import create_constants, W_even, W_odd, G, beta, beta_TE, beta_U
import numpy as np

#---------Constants--------#
#lambda_1 = 0.7 #expresed in mic.
lambda_1 = 1.55 #expresed in mic.
a = 2.5 #expresed in mic.
n1 =1.6
n2=1.5

#------Get first results----#
def calculations(a,ll,n1,n2):
    alpha, flags, ko, N, U, V = create_constants(a,ll,n1,n2)

    beta_w, alpha_solv = beta_U(a, alpha, beta, flags, ko, n1, V)
    U_solv = alpha_solv*a
    U_plot = U.flatten()
    V_plot = G(a,V,U.flatten()/a)
    ue = U[::2].flatten()
    weven = W_even(a,alpha[::2])
    for i in weven:
        i[-1] = None

    weven = weven.flatten()

    uo = U[1::2].flatten()
    wodd = W_odd(a,alpha[1::2])
    for i in wodd:
        i[-1] = None
    
    wodd= wodd.flatten()
    #print(beta_w)
    return U_plot, V_plot, beta_w, U_solv, ue, weven, uo, wodd, V

U_plot, V_plot, beta_w, U_solv, ue, weven, uo, wodd, V = calculations(a,lambda_1,n1,n2)

U_11 =  "Root for mode {0} = {1}".format(1, '%.3f' % U_solv[0]) 
beta_11 = "\u03B2 for mode {0} = {1}".format(1, '%.3f' %  beta_w[0])
# wee = np.isnan(weven) # weven[wee] = None # woo = np.isnan(wodd) # wodd[woo] = None # weven[weven==0] = None # wodd[wodd==0] = None

# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
from app import app

colors = { #Definition of colors for future use
    #'background': '#111111',
    #'text': '#FFFFFF',
    'text': '#111111',
    'background': '#FFFFFF',
    'circle': 'darkred',
    'even': 'darkgreen',
    'odd': 'darkblue'    
}
#-------- Sliders:
lambda_slider = dcc.Slider(
        id='lambda_slid_pw',
        min=0.7,
        max=1.8,
        step=0.05,
        value=lambda_1,
        marks={
        0.7: {'label': '0.7', 'style': {'color': colors['text']}},
        #lambda_1: {'label': str(lambda_1), 'style': {'color': colors['text']}},
        1.8: {'label': '1.8', 'style': {'color': colors['text']}}},
    )

d_slider = dcc.Slider(
        id='d_slid_pw',
        min=1,
        max=5,
        step=0.5,
        value=a,
        marks={
        1: {'label': '1', 'style': {'color': colors['text']}},
        5: {'label': str(5), 'style': {'color': colors['text']}}},
        #a: {'label': str(a), 'style': {'color': colors['text']}}},
    )
#--------------------------------------------------------------------------------#

def list_b(beta):
    num_modes = len(beta)
    mod_list = []
    for i in range(1,num_modes+1):
        dictModes = {'label': str(i), 'value': i}
        mod_list.append(dictModes)

    return mod_list

#----------Dropdowns-----------------------------------#
#-----Number of modes
num_mode_selec = dcc.Dropdown(id='n_modes_dropdown',
            options=list_b(beta_w),
            value=1,
            placeholder='Select a number'
    )











#----add lines for future plot-----#
ff= go.Scatter(x=U_plot,y=V_plot, name = 'V',
                         line=dict(color=colors['circle'], width=0.8, dash='dash'))

even_plot= go.Scatter(x=ue,y=weven, name = 'even',
                         line=dict(color=colors['even']))

odd_plot= go.Scatter(x=uo,y=wodd, name = 'odd',
                         line=dict(color=colors['odd']))

custom_F = go.Figure(data=[ff, even_plot, odd_plot])

eigen_plot=dcc.Graph(id='second_option_pw',
          animate=True,
          figure=custom_F.update_layout( 
                            width=600, height=600,
                            plot_bgcolor  = colors['background'],
                            paper_bgcolor = colors['background'],
                            font= {
                                    'color': colors['text']},
                            yaxis=dict(zeroline=True,zerolinewidth=0.5, zerolinecolor='lightslategrey',range=[0, int(V)+1],
                                        title='W',showgrid=True, gridwidth=0.5, gridcolor='lightslategrey'), 
                            xaxis=dict(zeroline=True,zerolinewidth=0.5, zerolinecolor='lightslategrey',range=[0, int(V)+1],
                                        title='U',showgrid=True, gridwidth=0.5, gridcolor='lightslategrey'), 
                            title='Visualization of graphic solution of the eigenvalue problem',
                            )
          )
#--------------------------------------#
#-------Final Layout for plotting and display of sliders and constants-------#
layout = html.Div(style={'backgroundColor': colors['background']},
    children=[    
        html.Div(id='resdp-display-value'),
        dcc.Link('Go to Fiber page', href='/apps/results'),
        html.Div(id='modesdp-display-value'),    
        html.Div(className='rows', 
        children=[
            html.Div(className='four columns div-user-controls', style={'backgroundColor': colors['background']},
            children = [
                html.Div(),
                dcc.Link('Go to Modes for planar waveguide', href='/apps/modes_dsply'),
                html.Div(),
                dcc.Link('Go to GVD effect', href='/apps/gvd'),
                html.Div(),
                dcc.Link('Go to SPM effect', href='/apps/spm'), ##EDIT LINK TO OTHER PAGES
                html.Div(),
                dcc.Link('Go to Split-Step for NLSE', href='/apps/nlse'),
                html.H2('Parameters', style={'color': colors['text']}),
                html.H6('n1: '+str(n1), style={'color': colors['text']}),
                html.H6('n2: '+str(n2), style={'color': colors['text']}),
                html.H3('\u03BB: '+str(lambda_1)+' \u03BCm', id='lambda_val_pw', style={'color': colors['text']}),
                lambda_slider,
                html.H3('a: '+str(a)+' \u03BCm', id='d_val_pw', style={'color': colors['text']}),
                d_slider,
                html.H6('Select desired mode number:', style={'color': colors['text']}),
                num_mode_selec,
                #html.H3('Beta:', style={'color': colors['text']}),
                html.H3(beta_11, id='beta_value_pw', style={'color': colors['text']}),
                html.H3(U_11, id='alpha_value', style={'color': colors['text']}),
                
            ]),  # Define the left element
            html.Div(className='eight columns div-for-charts bg-grey', style={'backgroundColor': colors['background']},  
            children = [
                     html.H1(children='Eigenvalue Problem',
                            style={
                            'textAlign': 'center',
                            'color': colors['text']
                    }),

                    html.Div(style={
                        'textAlign': 'center',
                        'color': colors['text']
                    }),
                    eigen_plot,        #Using Plotly
            ])  
        ]),
    ])
#-------------------------------------------------------------------------#
# callback for the sliders:
@app.callback( 
    [Output('second_option_pw', 'figure'),
    Output('lambda_val_pw', 'children'),
    Output('d_val_pw', 'children'),
    Output('beta_value_pw', 'children'),
    Output('alpha_value', 'children'),
    Output('n_modes_dropdown','options')],
    [Input('lambda_slid_pw', 'value'),
    Input('d_slid_pw', 'value'),
    Input('n_modes_dropdown', 'value')])
# function to update with the callback:
def update_plot(new_lambda, new_d, new_beta):
    U_plot_n, V_plot_n, beta_w_n, U_solv_n, ue_n, weven_n, uo_n, wodd_n, V_n = calculations(new_d,new_lambda,n1,n2)

    if new_beta == None: # For error while deselect of number, so new_number = None
        new_beta = 1
    if new_beta > len(beta_w_n): #For error while coming back from high order modes
        new_beta = len(beta_w_n)
    
    mod_list = list_b(beta_w_n)
    
    U_sn =  "Root for mode {0} = {1}".format(1, '%.3f' % U_solv_n[new_beta-1]) 
    beta_sn = "\u03B2 for mode {0} = {1}".format(1, '%.3f' %  beta_w_n[new_beta-1])

    ff= go.Scatter(x=U_plot_n,y=V_plot_n)

    even_plot= go.Scatter(x=ue_n,y=weven_n)

    odd_plot= go.Scatter(x=uo_n,y=wodd_n)
                            
    custom_F = go.Figure(data=[ff, even_plot, odd_plot]).update_layout( 
                            yaxis=dict(range=[0, int(V_n)+1],title='W'),#,showgrid=True, gridwidth=0.8, gridcolor='#121212'), 
                            xaxis=dict(range=[0, int(V_n)+1],title='U')#,showgrid=True, gridwidth=0.8, gridcolor='#121212')
                            )




    return (custom_F, 
            '\u03BB: '+str(new_lambda) +' \u03BCm', 
            'a: '+str(new_d) +' \u03BCm',
            beta_sn, 
            U_sn,
            mod_list )

if __name__ == '__main__':
    app.run_server(debug=True)