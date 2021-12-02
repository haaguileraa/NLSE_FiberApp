# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objects as go   
import numpy as np
from Modes import*


# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']


from app import app

colors = { #Definition of colors for future use
    #'background': '#111111',
    #'text': '#FFFFFF',
    'text': '#111111',
    'background': '#FFFFFF',
    'circle': 'whitesmoke',
    'even': 'darkgreen',
    'odd': 'darkblue'    
}

diameter_1 = 3.0
radius_i = diameter_1/2 + 1 #for Poynting vector plot
wavelength_1 = 1.596
initial_mode = 'HE'

#Definition of myMode
myMode = HybridMode(wavelength_1, diameter_1, typeMode=initial_mode, nu=1)
U_erst =  "Root for mode {0} = {1}".format(1, '%.5f' % myMode.solveEigenValuePb(1)) 
beta_erst = "\u03B2 for mode {0} = {1}".format(1, '%.5f' %  myMode.computeBeta(1))

last_update = 'Last update: 16.11.2021'

def list_m(mode):
    num_modes = mode.nbreOfModes()
    mod_list = []
    for i in range(1,num_modes+1):
        dictModes = {'label': str(i), 'value': i}
        mod_list.append(dictModes)

    return mod_list






#-------- SLIDERS---------------

wavelength_slider = dcc.Slider(
        id='wavelength_slid',
        min=0.4,
        max=2.0,
        step=0.004,
        value=wavelength_1,
        marks={
        0.4: {'label': '0.4', 'style': {'color': colors['text']}},
        #wavelength_1: {'label': str(wavelength_1), 'style': {'color': colors['text']}},
        2: {'label': '2.0', 'style': {'color': colors['text']}}},
    )


diameter_slider = dcc.Slider(
        id='diameter_slid',
        min=0.75,
        max=4.0,
        step=0.5,
        value=diameter_1,
        marks={
        0.75: {'label': '0.75', 'style': {'color': colors['text']}},
        4: {'label': '4.0', 'style': {'color': colors['text']}}},
    )


#-----------------DROPDOWNS------------------

#---------------Modes Selector:--------------
mode_selec = dcc.Dropdown(id='mode_dropdown',
            options=[
            {'label': 'HE', 'value': 'HE'},
            {'label': 'EH', 'value': 'EH'},
            {'label': 'TE', 'value': 'TE'},
            {'label': 'TM', 'value': 'TM'},
            ],
            value=initial_mode,
            placeholder="Select a mode"
    )
#--------------------//-----------------------
#-------------Material selectors--------------
matcore_selec = dcc.Dropdown(id='matcore_dropdown',
            options=[
            {'label': 'GeO2_193', 'value': 'GeO2_193'},
            {'label': 'GeO2_63', 'value': 'GeO2_63'},
            {'label': 'P2O5', 'value': 'P2O5'},
            {'label': 'silica', 'value': 'silica'}, 
            {'label': 'B2O3', 'value': 'B2O3'}, #commented due to its behaviour
            {'label': 'air', 'value': 'air'},
            #add more
            ],
            value='GeO2_193',
            placeholder="Select a material"
    )

matcladd_selec = dcc.Dropdown(id='matcladd_dropdown',
            options=[
           
            {'label': 'GeO2_63', 'value': 'GeO2_63'},
            {'label': 'P2O5', 'value': 'P2O5'},
            {'label': 'silica', 'value': 'silica'},
            {'label': 'B2O3', 'value': 'B2O3'}, #commented due to its behaviour 
            {'label': 'air', 'value': 'air'},
            #add more
            ],
            value='air',
            placeholder="Select a material"
    )
#--------------------//-----------------------
#---------------Number of modes---------------
number_mode_selec = dcc.Dropdown(id='number_mode_dropdown',
            options=list_m(myMode),
            value=1,
            placeholder='Select a number'
    )
#--------------------//-----------------------
#--------------MODE SELECTION-----------------
selection_mode = dcc.Dropdown(id='selec_ax',
            options=[
            {'label': 'Ex', 'value': 'ex'},
            {'label': 'Ey', 'value': 'ey'},
            {'label': 'Ez', 'value': 'ez'},
            {'label': 'Hx', 'value': 'hx'},
            {'label': 'Hy', 'value': 'hy'},
            {'label': 'Hz', 'value': 'hz'},
            {'label': 'Poynting Vector', 'value': 's'},
            ],
            value='s',
            placeholder="Select type"
    )

#--------------------//-----------------------
#-------------------NU-----------------------
selection_nu = dcc.Dropdown(id='nu_dropdown',
            options=[
            {'label': '1', 'value': 1},
            {'label': '2', 'value': 2},
            {'label': '3', 'value': 3},
            {'label': '4', 'value': 4},
            {'label': '5', 'value': 5},
            ],
            value=1,
            placeholder="Select index"
    )
#--------------------//-----------------------
#-------------------CUT-----------------------
cut_selec = dcc.Dropdown(id='cut_dropdown',
            options=[
            {'label': 'Cut', 'value': 'cut'},
            {'label': 'Cut H', 'value': 'cuth'},
            #{'label': 'TM', 'value': 'TM'},
            ],
            value=None,
            placeholder="Select Cut type"
    )
#--------------------//-----------------------
#--------------------//-----------------------

#------------PLOTS-------------------

Plot = myMode.ploteigen(0.5,30,250,colors)
Fig = go.Figure(data=[Plot])
plot_eigen = dcc.Graph(id='plot_eigen',
            animate= True,
            figure=Fig.update_layout( 
                                width=600, height=600,
                                plot_bgcolor  = colors['background'],
                                paper_bgcolor = colors['background'],
                                font= {
                                        'color': colors['text']},
                                yaxis=dict(title = 'Eigenvalue Pb', zeroline=True,zerolinewidth=0.8, zerolinecolor='#CCCCCC',
                                            mirror=True,
                                            ticks='outside',
                                            showline=True,
                                            linecolor='lightslategrey', 
                                            range=[-5,5],showgrid=True, gridwidth=0.8, gridcolor='#CCCCCC'), #,title='$\u03B7_{eff}$'
                                xaxis=dict(title = 'U',showgrid=True, gridwidth=0.8, gridcolor='#CCCCCC',
                                            autorange= True,
                                            mirror=True,
                                            ticks='outside',
                                            showline=True,
                                            linecolor='lightslategrey', ),# 
                                title='Eigenvalue Pb',
                                ))


modeplot, mode_title, xdata, ydata = myMode.plotMode_dash(xmin=-radius_i,xmax=radius_i, ymin=-radius_i, ymax=radius_i)
plot_mode = dcc.Graph(id='mode_plot',
            #style={"width": "75%", "display": "inline-block"},
            figure=modeplot.update_layout(
                                width=600, height=600,
                                font= {
                                        'color': colors['text']},
                                yaxis=dict(title = ydata[0]),
                                xaxis=dict(title = xdata[0]),
                                ))

cutplot = go.Figure(data=[])
cut_pl = dcc.Graph(id='cut_plot',
                figure=cutplot.update_layout( 
                                width=600, height=600,
                                # paper_bgcolor = colors['background'],
                                # plot_bgcolor  = colors['background'],
                                # yaxis=dict(showgrid=True, gridwidth=0.8, gridcolor='#CCCCCC'), #,title='$\u03B7_{eff}$'
                                # xaxis=dict(showgrid=True, gridwidth=0.8, gridcolor='#CCCCCC'), 
                                ))


beta_plot = plot_beta_2(diameter_1, 0.7, 1.8, 50, dependence = 'lambda')
beta_pl =dcc.Graph(id='second_option_pw',
          animate=True,
          figure=beta_plot.update_layout( 
                            width=600, height=600,
                            plot_bgcolor  = colors['background'],
                            paper_bgcolor = colors['background'],
                            font= {
                                    'color': colors['text']},
                            yaxis=dict(zeroline=True,zerolinewidth=0.5, zerolinecolor='lightslategrey',
                                        mirror=True,
                                        ticks='outside',
                                        showline=True,
                                        linecolor='lightslategrey',
                                        #title='\u03B2 2 (ps^2/km) ',showgrid=True, gridwidth=0.5, gridcolor='lightslategrey'), 
                                        title='D (ps^2/km) ',showgrid=True, gridwidth=0.5, gridcolor='lightslategrey'), 
                            xaxis=dict(zeroline=True,zerolinewidth=0.5, zerolinecolor='lightslategrey',
                                        mirror=True,
                                        ticks='outside',
                                        showline=True,
                                        linecolor='lightslategrey',
                                        title='Wavelength \u03BB (\u03BCm)',showgrid=True, gridwidth=0.5, gridcolor='lightslategrey'), 
                            title='Plot of beta',
                            )
          )
#--------------------------------------#
                            

#-------Final Layout for plotting and display of sliders and constants-------#
layout = html.Div(style={'backgroundColor': colors['background']},
    children=[
        html.Div(id='dpf-display-value'),
        dcc.Link('Go to planar waveguide page', href='/apps/dash_plot'),
        html.Div(id='modesf-display-value'),
        dcc.Link('Go to Modes for planar waveguide page', href='/apps/modes_dsply'),
        html.Div(className='rows', 
        children=[
            html.Div(className='four columns div-user-controls', style={'backgroundColor': colors['background']}, 
            children = [
                html.P('Front-End for understanding the plot inside original modes.py file', style={'color': colors['text']}),
                html.H6('Select material for the core:', style={'color': colors['text']}),
                matcore_selec,
                html.H6('Select material for the cladding:', style={'color': colors['text']}),
                matcladd_selec,
                html.H6('Select Mode:', style={'color': colors['text']}),
                mode_selec,
                html.H3('D: '+str(diameter_1)+ ' \u03BCm', id='diameter_val', style={'color': colors['text']}), #Ro: \u03C1
                diameter_slider,
                html.H3('\u03BB: '+str(wavelength_1)+' \u03BCm', id='wavelength_val', style={'color': colors['text']}),
                wavelength_slider,
                html.H6('Select desired nu:', style={'color': colors['text']}),
                selection_nu,
                html.H6('Select desired mode number:', style={'color': colors['text']}),
                number_mode_selec,
                html.H3(U_erst, id='root_value', style={'color': colors['text']}),
                html.H3(beta_erst, id='beta_value', style={'color': colors['text']}),
                
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
                    plot_eigen,
                    selection_mode,        
                    plot_mode,
                    cut_selec,
                    cut_pl,
                    #beta_pl,
            ])  
        ]),
        html.Div(),
        dcc.Link('Go to GVD effect', href='/apps/gvd'),
        html.Div(),
        dcc.Link('Go to SPM effect', href='/apps/spm'), ##EDIT LINK TO OTHER PAGES
        html.Div(),
        dcc.Link('Go to Split-Step for NLSE', href='/apps/nlse'),
        html.Footer('Joly Nicolas, Aguilera Hernan. Max-Planck-Institut', style={'color': colors['text']}),
        html.Footer(last_update, style={'color': colors['text']}),
    ])


@app.callback( 
    [Output('plot_eigen', 'figure'),
    Output('mode_plot', 'figure'),
    Output('wavelength_val', 'children'),
    Output('diameter_val', 'children'),
    Output('number_mode_dropdown','options'),
    Output('matcore_dropdown','options'),
    Output('matcladd_dropdown','options'),
    Output('selec_ax','options'),
    Output('root_value','children'),
    Output('beta_value','children'),
    Output('cut_plot','figure'),
    ],

    [Input('wavelength_slid', 'value'),
    Input('diameter_slid', 'value'),
    Input('mode_dropdown', 'value'),
    Input('selec_ax', 'value'),
    Input('number_mode_dropdown', 'value'),
    Input('cut_dropdown', 'value'),
    Input('matcore_dropdown', 'value'),
    Input('matcladd_dropdown', 'value'),
    Input('nu_dropdown', 'value'),
    ]
    )

# function to update with the callback:
def update_plot(new_wavelength, new_diameter, new_mode, 
                new_plot_mode, new_number, new_cut, new_mcore, 
                new_mcladd,new_nu):
    
    if new_mcore == None:
        new_mcore = 'GeO2_193'
    if new_mcladd == None:
        new_mcladd = 'air'
    
    if new_mode == 'TE':
        if new_plot_mode == 'ez':
            new_plot_mode = 'ex'
        opt_selec_ax = [{'label': 'Ex', 'value': 'ex'}, {'label': 'Ey', 'value': 'ey'},
                        {'label': 'Hx', 'value': 'hx'}, {'label': 'Hy', 'value': 'hy'},
                        {'label': 'Hz', 'value': 'hz'}, {'label': 'Poynting Vector', 'value': 's'},]
    elif new_mode == 'TM':
        if new_plot_mode == 'hz':
            new_plot_mode = 'hx'
        opt_selec_ax = [{'label': 'Ex', 'value': 'ex'}, {'label': 'Ey', 'value': 'ey'},
                        {'label': 'Ez', 'value': 'ez'}, {'label': 'Hx', 'value': 'hx'},
                        {'label': 'Hy', 'value': 'hy'}, {'label': 'Poynting Vector', 'value': 's'},]
    else: 
        opt_selec_ax = [{'label': 'Ex', 'value': 'ex'}, {'label': 'Ey', 'value': 'ey'},
                        {'label': 'Ez', 'value': 'ez'}, {'label': 'Hx', 'value': 'hx'},
                        {'label': 'Hy', 'value': 'hy'}, {'label': 'Hz', 'value': 'hz'},
                        {'label': 'Poynting Vector', 'value': 's'},]
    my_newMode = calculate_CB_mode(new_wavelength, new_diameter, new_mode, 
                    new_plot_mode, new_mcore, new_mcladd,new_nu)

    newMaterials = matDropdowns(my_newMode)


    new_num_modes = my_newMode.nbreOfModes()

    if new_number == None: # For error while deselect of number, so new_number = None
        new_number = 1
    elif new_number > new_num_modes: #For error while coming back from high order modes
        new_number = new_num_modes

    mod_list = list_m(my_newMode)

    Plot_new = my_newMode.ploteigen(0.5,30,250,colors) #, plot_title_new
    Fig_n = go.Figure(data=[Plot_new]).update_layout(xaxis = dict(autorange= True,)).update()
    new_root =  "Root for mode {0} = {1}".format(new_number, '%.5f' % my_newMode.solveEigenValuePb(new_number)) 
    beta = "\u03B2 for mode {0} = {1}".format(new_number, '%.5f' % my_newMode.computeBeta(new_number))
    
    radius_plot = new_diameter/2+1

    new_modeplot, _, _, _ = my_newMode.plotMode_dash(xmin=-radius_plot,xmax=radius_plot, ymin=-radius_plot, ymax=radius_plot,ax=new_plot_mode)

    if new_cut == 'cut':
        new_cut_plot = go.Figure(data=[my_newMode.drawCut_dash()]).update_layout(title= 'R Field Cut',
                                width=600, height=600, 
                                paper_bgcolor = colors['background'],
                                plot_bgcolor  = colors['background'],
                                yaxis=dict(showgrid=True, gridwidth=0.8, gridcolor='#CCCCCC',
                                            zeroline=True,zerolinewidth=0.8, zerolinecolor='#CCCCCC'), #,title='$\u03B7_{eff}$'
                                xaxis=dict(showgrid=True, gridwidth=0.8, gridcolor='#CCCCCC',
                                            zeroline=True,zerolinewidth=0.8, zerolinecolor='#CCCCCC'),).update()
    elif new_cut == 'cuth':
        new_cut_plot = go.Figure(data=[my_newMode.drawCutH_dash()]).update_layout(title= 'H Field Cut ',
                                width=600, height=600, 
                                paper_bgcolor = colors['background'],
                                plot_bgcolor  = colors['background'],
                                yaxis=dict(showgrid=True, gridwidth=0.8, gridcolor='#CCCCCC',
                                            zeroline=True,zerolinewidth=0.8, zerolinecolor='#CCCCCC'), #,title='$\u03B7_{eff}$'
                                xaxis=dict(showgrid=True, gridwidth=0.8, gridcolor='#CCCCCC',
                                            zeroline=True,zerolinewidth=0.8, zerolinecolor='#CCCCCC'),).update()
    else:
        new_cut_plot = go.Figure(data=[])

    return (Fig_n,
            new_modeplot,
            '\u03BB: '+str(new_wavelength) +' \u03BCm', 
            'd: '+str(new_diameter) +' \u03BCm', 
            mod_list,
            newMaterials[0],
            newMaterials[1],
            opt_selec_ax,
            new_root,
            beta,
            new_cut_plot
            )
