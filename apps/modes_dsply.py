# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objects as go   

from Functions import create_constants, f_TE, bisection, beta, beta_TE, beta_TM, beta_tay, W_even, W_odd, G, modes
import numpy as np

#---------Constants--------#
#lambda_1 = 0.7 #expresed in mic.
lambda_1 = 1.55 #expresed in mic.
a = 2.5 #expresed in mic.
n1 =2
n2=1.5


normalized_freq = np.arange(0.4,6.5,0.4)

TE, lambda_V = modes(a, beta_TE, n1, n2, normalized_freq)
TM, lambda_V = modes(a, beta_TM, n1, n2, normalized_freq)


alpha, flags, ko, N, U, largest_V = create_constants(a,np.amin(lambda_V),n1,n2)

for i in TE:
        i[-1] = None
for i in TM:
        i[-1] = None

normalized_freq =  np.ones((TE.shape[0], normalized_freq.shape[0]))*normalized_freq

normalized_freq, TE, TM = normalized_freq.flatten(), TE.flatten(), TM.flatten()

# wee = np.isnan(weven) # weven[wee] = None # woo = np.isnan(wodd) # wodd[woo] = None # weven[weven==0] = None # wodd[wodd==0] = None

# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

colors = { #Definition of colors for future use
    #'background': '#111111',
    #'text': '#FFFFFF',
    'text': '#111111',
    'background': '#FFFFFF',
    'circle': 'whitesmoke',
    'even': 'darkgreen',
    'odd': 'darkblue'    
}
#-------- Sliders:

TE_modes = go.Scatter(x=normalized_freq ,y=TE, name = 'TE',
                         line=dict(color=colors['even']))

TM_modes = go.Scatter(x=normalized_freq ,y=TM, name = 'TM',
                         line=dict(color=colors['odd']))


custom_modes = go.Figure(data=[TE_modes, TM_modes])

brech_norm_fre = dcc.Graph(id='norm_brech',
          animate=True,
          figure=custom_modes.update_layout( 
                            width=600, height=600,
                            plot_bgcolor  = colors['background'],
                            paper_bgcolor = colors['background'],
                            font= {
                                    'color': colors['text']},
                            yaxis=dict(range=[n2+0.01, n1],title='''
                                        \u03B7<sub>eff</sub>
                                        ''',showgrid=True, gridwidth=0.8, gridcolor='#121212'), 
                            xaxis=dict(range=[normalized_freq[0], normalized_freq[-1]],title='V',showgrid=True, gridwidth=0.8, gridcolor='#121212'), 
                            title= '''
                                    Plot of  \u03B7<sub>eff</sub> vs. V
                                    ''',#'\u03B7 eff vs V',
                            )
          )

#--------------------------------------#
#-------Final Layout for plotting and display of sliders and constants-------#
layout = html.Div(style={'backgroundColor': colors['background']},
    children=[
        html.Div(id='resmodes-display-value'),
        dcc.Link('Go to Fiber page', href='/apps/results'),
        html.Div(id='dpmodes-display-value'),
        dcc.Link('Go to planar waveguide page', href='/apps/dash_plot'),
        html.Div(className='rows', 
        children=[
            html.Div(className='four columns div-user-controls', style={'backgroundColor': colors['background']}, 
            children = [
                html.H2('Parameters', style={'color': colors['text']}),
                html.H6('n1: '+str(n1), style={'color': colors['text']}),
                html.H6('n2: '+str(n2), style={'color': colors['text']}),
                html.P('Effective index versus normalized frequency V of a typical symetric waveguide', style={'color': colors['text']}),

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
                    brech_norm_fre        #Using Plotly
            ])  
        ]),
        html.Div(),
        dcc.Link('Go to GVD effect', href='/apps/gvd'),
        html.Div(),
        dcc.Link('Go to SPM effect', href='/apps/spm'), ##EDIT LINK TO OTHER PAGES
        html.Div(),
        dcc.Link('Go to Split-Step for NLSE', href='/apps/nlse'),
        
    ])
#-------------------------------------------------------------------------#
# callback for the sliders:
# @app.callback( 
#     [Output('second_option', 'figure'),
#     Output('lambda_val', 'children'),
#     Output('d_val', 'children'),
#     Output('beta_value', 'children')],
#     [Input('lambda_slid', 'value'),
#     Input('d_slid', 'value')])
# # function to update with the callback:
# def update_plot(new_lambda, new_d):
#     U_plot_n, V_plot_n, beta_w_n, ue_n, weven_n, uo_n, wodd_n, V_n = calculations(new_d,new_lambda,n1,n2)

#     ff= go.Scatter(x=U_plot_n,y=V_plot_n)

#     even_plot= go.Scatter(x=ue_n,y=weven_n)

#     odd_plot= go.Scatter(x=uo_n,y=wodd_n)
                            
#     custom_F = go.Figure(data=[ff, even_plot, odd_plot]).update_layout( 
#                             yaxis=dict(range=[0, int(V_n)+1],title='W',showgrid=True, gridwidth=0.8, gridcolor='#121212'), 
#                             xaxis=dict(range=[0, int(V_n)+1],title='U',showgrid=True, gridwidth=0.8, gridcolor='#121212')
#                             )
#     return custom_F, '\u03BB: '+str(new_lambda) +' \u03BCm', 'a: '+str(new_d) +' \u03BCm', '%.4f' % beta_w_n[0]

# if __name__ == '__main__':
#     app.run_server(debug=True)